#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Construct and certify the fixed-K, variable-shear active-v2 inner inlet."""

import argparse
import math
from pathlib import Path
import sys
import time

import physical_oracle as po
from active_v2_inlet_profile import MODEL, SCHEMA, TRANSPORT, Profile
from inner_k_columns import planar_response, stress
from inner_k_fixture_oracle import atomic_out
from inner_k_mechanism import read
import numpy as np
from numpy.polynomial.chebyshev import chebfit, chebval, chebder, chebint
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


class Column:
    def __init__(self, ref, deadline):
        self.ref, self.deadline = ref, deadline
        self.c = float(ref.cutoff(ref.Kin))
        self.panels = sorted({0.,self.c,.4*ref.H,ref.h,.75*ref.H,ref.H})
        self.trace = []

    def shear(self, r, p, T):
        a = self.ref
        k = a.Kin*p
        po.require(k>0 and T>0 and math.isfinite(k+T),'nonpositive trial k or inlet shear')
        E = a.A/max(r,self.c)
        u = (r/a.H-.4)/.35
        chi = 0. if u<=0 else 1. if u>=1 else u**3*(10-15*u+6*u*u)
        aux = a.a*(r-a.zv)
        alpha = a.nu+(1-chi)*aux
        if chi==0.:
            G, tangent, Tk = T/alpha, alpha, 0.
        else:
            # With s=q/N, the simple-shear cubic gives
            # N=1.8(1+s^2)/(1-.35s^2) and -R_xr/k=.6s/(1+s^2).
            # s in [0,1] parametrizes the increasing full WJ shear branch.
            def reduced(s):
                z = s*s
                return 1.8*s*(1+z)/(1-.35*z), .6*s/(1+z)
            upper = alpha*E/k*(3.6/.65)+chi*k*.3
            po.require(T<=upper,f'shear branch capacity {upper:.17g} < {T:.17g} at r={r:.17g}')
            def imbalance(s):
                q,f = reduced(s)
                return alpha*E*q/k+chi*k*f-T
            s = brentq(imbalance,0.,1.,xtol=2e-15)
            q,f = reduced(s)
            z = s*s
            q_s = 1.8*((1+z)/(1-.35*z)+2.7*z/(1-.35*z)**2)
            fq = .6*(1-z)/(1+z)**2/q_s
            G = E*q/k
            tangent = alpha+chi*k*k/E*fq
            Tk = chi*a.Kin*(f+q*fq)
        diffusion = a.nu+a.cd*k*k/(a.sk*E)
        return G,E,diffusion,tangent,Tk,aux,chi

    def integrate(self, unknowns, tight=False, initial=None, lower=None):
        a = self.ref
        T = float(unknowns[0])*a.tau
        forward = initial is not None
        value = np.asarray(initial if forward else [1.,unknowns[1],0.,0.,0.,0.,1.,0.,0.])
        rtol,maxstep = (5e-13,.0125) if tight else (2e-11,.025)
        method = 'RK45' if forward else 'DOP853'
        panels = []

        def rhs(y, v):
            if time.monotonic()>=self.deadline:
                raise TimeoutError('inlet IVP exceeded its cooperative bound')
            p,Q = v[:2]
            G,E,D,tangent,Tk,_,_ = self.shear(y*a.H,p,T)
            production = (T-a.nu*G)*G
            result = [a.H*a.A*Q/(a.Kin*D),a.H*(E-production)/a.A,a.H*G]
            if forward:
                return result
            dG = np.array([a.tau,0.])-Tk*v[3:5]
            dG /= tangent
            dP = (T-2*a.nu*G)*dG+np.array([a.tau*G,0.])
            dD = 2*(D-a.nu)*v[3:5]/p
            return np.r_[result,a.H*a.A/(a.Kin*D)*(v[5:7]-Q*dD/D),-a.H*dP/a.A,a.H*dG]

        po.require(lower is None or (not forward and lower in self.panels),
                   'requested stop is not a backward material break')
        y = np.asarray([r for r in self.panels if lower is None or r>=lower])/a.H
        if not forward:
            y = y[::-1]
        for start,end in zip(y[:-1],y[1:]):
            ode = solve_ivp(rhs,(start,end),value,method=method,rtol=rtol,atol=.04*rtol,
                            max_step=maxstep,dense_output=True)
            po.require(ode.success and ode.t[-1]==end,'inlet IVP failed on an admitted panel')
            value = ode.y[:,-1]
            panels.append((min(start,end),max(start,end),ode.sol))

        def sample(r):
            z = r/a.H
            for lo,hi,sol in panels:
                if lo-4e-15<=z<=hi+4e-15:
                    return sol(z)[:3]
            raise ValueError('column witness outside [0,H]')
        return value,sample

    def residual(self, value, target):
        return np.array([value[1],-value[2]-target]),np.array([value[5:7],-value[7:9]])

    def correct(self, x, target, tight):
        value, sample = self.integrate(x,tight)
        for iteration in range(11):
            f,J = self.residual(value,target)
            norm = np.linalg.norm(f,np.inf)
            self.trace.append(dict(target=target,unknowns=x.tolist(),tight=tight,iteration=iteration,
                                   defect=f.tolist(),condition=float(np.linalg.cond(J))))
            if norm<=4e-13:
                return x,value,sample,J
            step = np.linalg.solve(J,-f)
            for divisor in (1,2,4,8,16,32):
                trial = x+step/divisor
                try:
                    v,s = self.integrate(trial,tight)
                except ValueError as error:
                    self.trace.append(dict(target=target,unknowns=trial.tolist(),rejection=str(error)))
                    continue
                test = self.residual(v,target)[0]
                if np.linalg.norm(test,np.inf)<norm:
                    x,value,sample = trial,v,s
                    break
                self.trace.append(dict(target=target,unknowns=trial.tolist(),rejection='no endpoint decrease',
                                       defect=test.tolist()))
            else:
                raise ValueError('no admissible decreasing shooting step')
        raise ValueError('shooting correction exceeded 11 iterations')


def independent(column, x, terminal, sample, target, require_admissible=True, forward_control=True):
    ref = column.ref
    T = ref.tau*float(x[0])
    QH = ref.A*float(x[1])
    r = np.unique(np.r_[np.linspace(0.,ref.H,513),column.panels])

    def fields(rr, s, anchor):
        y = np.array([s(v) for v in rr])
        kin = y[:,0]
        actions = np.array([column.shear(v,p,T) for v,p in zip(rr,kin)])
        G,E,D = actions[:,:3].T
        return y[:,2]+anchor,ref.Kin*kin,G,ref.A*y[:,1]/D,E,actions

    U,k,G,kr,E,actions = fields(r,sample,target)
    endpoints = dict(U0=float(target+terminal[2]),Q0_over_A=float(terminal[1]))
    profile_error = None
    if forward_control:
        forward_end, forward = column.integrate(x,tight=True,initial=[terminal[0],0.,0.])
        Uf,kf,Gf,krf,_,_ = fields(r,forward,0.)
        profile_error = float(max(np.max(abs(U-Uf)),np.max(abs(k-kf))/ref.Kin,np.max(abs(G-Gf)),
                                  np.max(abs(kr-krf))*ref.H/ref.Kin))
        endpoints.update(forward_UH_error=float(forward_end[2]-target),
                         forward_KH_error_over_Kin=float(forward_end[0]-1.),
                         forward_QH_error_over_A=float(forward_end[1]-x[1]))
        po.require(profile_error<=1e-8,'independent forward profile gate failed')
    guard, _ = planar_response(k,E,G,actions[:,5],actions[:,6],ref.nu)
    admissible = (min(k)>0 and min(E)>0 and min(G)>0 and min(actions[:,3])>0 and
               guard['minimum_covariance_over_k']>1e-10 and guard['minimum_principal_over_mu']>1e-8 and
               guard['tangent_refinement_difference']<=1e-8)
    po.require(admissible or not require_admissible,'sampled amended-v2 inlet admissibility failed')
    audits = []
    for order in (32,48,64):
        q,w = leggauss(order)
        rr = np.concatenate([(a+b)/2+(b-a)*q/2 for a,b in zip(column.panels,column.panels[1:])])
        weight = np.concatenate([w*(b-a)/2 for a,b in zip(column.panels,column.panels[1:])])
        _,kk,gg,_,ee,info = fields(rr,sample,target)
        grad = np.zeros((len(rr),3,3))
        grad[:,0,2] = gg
        tensor = stress(kk,ee,grad,info[:,5],info[:,6])
        P = -tensor[:,0,2]*gg
        ip,ie = float(weight@P),float(weight@ee)
        audits.append(dict(order=order,production=ip,destruction=ie,
            budget_error_over_A=(ip-ie+QH-ref.A*terminal[1])/ref.A,
            shear_error_over_tau=float(np.max(abs(ref.nu*gg-tensor[:,0,2]-T)))/ref.tau))
    po.require(max(map(abs,endpoints.values()))<=1e-10 and
               max(abs(b[key]) for b in audits for key in ('budget_error_over_A','shear_error_over_tau'))<=1e-10,
               'independent fixed-endpoint or full-stress gate failed')

    g_out = ref.utau/(ref.kappa*ref.H)
    grad = np.zeros((2,3,3))
    grad[:,0,2] = [G[-1],g_out]
    outer_stress = po.full_wj(np.full(2,ref.Kin),np.full(2,E[-1]),grad)[0]
    traction = (ref.nu*grad[:,0,2]-outer_stress[:,0,2]).tolist()
    interface = dict(r=ref.H,U=float(U[-1]),k=float(k[-1]),E=float(E[-1]),
                     U_r=grad[:,0,2].tolist(),k_r=[float(kr[-1]),0.],E_r=[-float(E[-1])/ref.H]*2,
                     finite_nu_total_shear=traction,shear_inner_minus_outer=traction[0]-traction[1],
                     stress_normal_component=outer_stress[:,2,2].tolist(),Q=[QH,0.])
    return dict(T_in=T,T_in_over_tau=T/ref.tau,Q_H=QH,Q_H_over_A=float(x[1]),
                endpoint_errors=endpoints,quadrature=audits,independent_forward_profile_error=profile_error,
                min_k=float(min(k)),min_E=float(min(E)),min_shear_tangent=float(min(actions[:,3])),
                guard=guard,sampled_admissibility_pass=bool(admissible),interface=interface), fields


def make_profile(column, solution, target, fields, degree):
    ref = column.ref
    # Subdivide the steep near-ground segment and keep every constitutive break.
    edges = sorted(set(column.panels+[.001,.0025,.005,.01,.02,.045,.085,.13]))
    nodes = np.cos(np.pi*np.arange(degree+1)/degree)
    panels = []
    for lo,hi in zip(edges[:-1],edges[1:]):
        t = (lo+hi)/2+(hi-lo)*nodes/2
        U,k,G,kr = fields(t,solution,target)[:4]
        # Integrate the independently sampled ODE derivatives. Differentiating
        # a high-order fit to rounded dense-output values amplifies roundoff,
        # especially on the thin wall panels. The tiny constant adjustment
        # keeps the verified ODE endpoint values exact to rounding.
        def interpolate(value, derivative):
            coeffs = chebint(chebfit(nodes,derivative,degree),k=[value[-1]],lbnd=-1.,scl=(hi-lo)/2)
            correction = (value[0]-chebval(1.,coeffs))/2
            coeffs[:2] += correction
            return coeffs.tolist()
        panels.append(dict(lo=lo,hi=hi,U=interpolate(U,G),k=interpolate(k,kr)))
    result = dict(schema=SCHEMA,model=MODEL,transport=TRANSPORT,reference_table_sha256=po.TABLE_SHA,
                  H=ref.H,Kin=ref.Kin,nu=ref.nu,cutoff=column.c,panels=panels)
    r = np.unique(np.r_[np.linspace(0,ref.H,2049),edges])
    error = {}
    U,k,G,kr = fields(r,solution,target)[:4]
    o = Profile(result,ref).evaluate(r)
    for key,exact,scale in (('U',U,1.),('k',k,ref.Kin),('G',G,1.),('kr',kr,ref.Kin/ref.H)):
        error[key] = float(np.max(abs(o[key]-exact)))/scale
    joins = []
    for a,b in zip(panels[:-1],panels[1:]):
        data = []
        for p,t in ((a,1.),(b,-1.)):
            data.append([chebval(t,p['U']),chebval(t,p['k']),
                         2*chebval(t,chebder(p['U']))/(p['hi']-p['lo']),
                         2*chebval(t,chebder(p['k']))/(p['hi']-p['lo'])])
        joins.append(float(np.max(abs(np.subtract(*data))/[1.,ref.Kin,1.,ref.Kin/ref.H])))
    po.require(max(*error.values(),*joins)<=1e-8,
               f'immutable panel profile exceeds representation budget: {error}, join={max(joins)}')
    return result,dict(derivative_degree=degree,samples=len(r),maximum_scaled_error=error,
                       maximum_scaled_join_jump=max(joins))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--table',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--seconds',type=float,required=True)
    args = p.parse_args()
    po.require(sys.platform.startswith('linux') and args.seconds>0 and
               po.identity(args.table)['sha256']==po.TABLE_SHA,'requires confirmed Linux calibration')
    ref = po.Reference(read(args.table))
    column = Column(ref,time.monotonic()+args.seconds)
    args.out.mkdir()
    try:
        x = np.array([1.,0.])
        first,_ = column.integrate(x)
        canonical = -float(first[2])
        target = float(ref.utau/ref.kappa*math.log(ref.H/ref.z0))
        po.require(abs(canonical-.42150691439960103)<1e-9,'canonical v2 column changed')
        for value in np.linspace(canonical,target,9):
            x,_,coarse,_ = column.correct(x,float(value),False)
        xu,terminal,sample,J = column.correct(x,target,True)
        radius = np.unique(np.r_[np.linspace(0,ref.H,513),column.panels])
        change = float(max(np.max(abs(sample(r)-coarse(r))) for r in radius))
        check,fields = independent(column,xu,terminal,sample,target)
        delta = 5e-6
        fd = np.column_stack([sum(w*column.residual(column.integrate(xu+m*d,True)[0],target)[0]
            for m,w in ((-2,1.),(-1,-8.),(1,8.),(2,-1.)))/(12*delta) for d in delta*np.eye(2)])
        jac_error = float(np.max(abs(fd-J)))
        po.require(jac_error<=1e-8 and change<=1e-8,'shooting sensitivity or IVP refinement failed')
        artifact, representation = make_profile(column,sample,target,fields,18)
        artifact.update(T_in=check['T_in'],Q_H=check['Q_H'],U_H=target)
        atomic_out(args.out/'profile.json',artifact)
        atomic_out(args.out/'assessment.json',dict(schema='WJH27MatchedInletColumnAssessmentV1',
            status='PASS',table=po.identity(args.table),profile=po.identity(args.out/'profile.json'),
            canonical_UH=canonical,canonical_Q0_over_A=float(first[1]),coarse_to_refined_change=change,
            jacobian=J.tolist(),independent_jacobian=fd.tolist(),jacobian_error=jac_error,
            jacobian_condition=float(np.linalg.cond(J)),representation=representation,**check))
    finally:
        atomic_out(args.out/'shooting.json',column.trace)
    print('ACTIVE_V2_MATCHED_INNER_COLUMN_PASS',flush=True)


if __name__=='__main__':
    main()
