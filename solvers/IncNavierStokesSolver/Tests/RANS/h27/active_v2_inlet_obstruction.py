#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check the prescribed branch limit using fixed-shear and two-ended solves."""

import argparse
import math
from pathlib import Path
import sys
import time

import physical_oracle as po
from active_v2_inlet import Column, independent
from inner_k_fixture_oracle import atomic_out
from inner_k_mechanism import read
from reference_oracle import Q_PEAK, shear_law
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, root


def lower_match(column, shear):
    a = column.ref
    r = .4*a.H
    y = a.a*r/a.D0
    # Substitution B=d+a r in integral(T^2 nu_A/(nu+nu_A)^2).
    J = math.log1p(y)-(a.nu/a.D0)*y/(1+y)
    P = shear**2/a.a*J
    E = a.A*(1+math.log(r/column.c))
    q,w = leggauss(64)
    direct_P, direct_E = 0.,0.
    for lo,hi in ((0.,column.c),(column.c,r)):
        rr = (lo+hi)/2+(hi-lo)*q/2
        aux = a.a*(rr-a.zv)
        direct_P += (hi-lo)/2*float(w@(aux*(shear/(a.nu+aux))**2))
        direct_E += (hi-lo)/2*float(w@(a.A/np.maximum(rr,column.c)))
    error = max(abs(direct_P-P),abs(direct_E-E))/a.A
    po.require(error<=1e-10,'lower fixed-shear quadrature disagrees with its rederived integral')
    return (E-P)/a.A,-2*shear*a.tau/a.a*J/a.A,error


def fixed_shear(column, t, guess, tight):
    a = column.ref
    T = a.tau*t
    stop = .4*a.H
    q_target, dtarget, integral_error = lower_match(column,T)
    upper_E = math.log(a.H/stop)
    lo = q_target+upper_E-(a.H-stop)*T*T/(4*a.nu*a.A)
    hi = q_target+upper_E
    cache = {}
    def residual(qh):
        key = float(qh)
        if key not in cache:
            try:
                v,_ = column.integrate([t,key],tight=tight,lower=stop)
                cache[key] = v
            except ValueError as error:
                cache[key] = str(error)
        value = cache[key]
        if isinstance(value,str):
            raise ValueError(value)
        return float(value[1]-q_target)

    probes = np.unique(np.r_[np.linspace(lo,hi,41),
                     [guess+e for e in (0.,-.1,.1,-.04,.04,-.02,.02,-.005,.005)
                      if lo<=guess+e<=hi]])
    for q in probes:
        try:
            residual(q)
        except ValueError:
            continue
    rows = [dict(Q_H_over_A=q,rejection=v) if isinstance(v,str) else
            dict(Q_H_over_A=q,defect=float(v[1]-q_target),slope=float(v[6]))
            for q,v in sorted(cache.items())]
    brackets = [(u['Q_H_over_A'],v['Q_H_over_A']) for u,v in zip(rows[:-1],rows[1:])
                if 'defect' in u and 'defect' in v and u['defect']*v['defect']<=0]
    result = dict(T_over_tau=t,T=T,limit_Q_H_over_A=[lo,hi],brackets=brackets,trials=rows,
                  lower_integral_error_over_A=integral_error)
    controls = []
    for left,right in brackets:
        try:
            qh = brentq(residual,left,right,xtol=2e-14)
            residual(qh)
            v = cache[qh]
            full, sample = column.integrate([t,qh],tight=tight)
        except ValueError as error:
            controls.append(dict(bracket=[left,right],rejection=str(error)))
            continue
        po.require(abs(full[1])<=1e-10,'fixed-shear solution misses the independent ground endpoint')
        dq = -(v[5]-dtarget)/v[6]
        G0 = T/a.D0
        controls.append(dict(Q_H_over_A=qh,U_H=float(-full[2]),k_ground=float(full[0]*a.Kin),
                             G_ground=G0,ground_residual_over_A=float(full[1]),
                             dUH_d_T_over_tau=float(-(full[7]+full[8]*dq)),
                             upper_match_slope=float(v[6])))
    return {**result,'solutions':controls}


def two_ended(column, t, k0, qh, near_anchor, primary_sample):
    """Independent split RK45 integration with the companion-root shear law."""
    a = column.ref
    T = a.tau*t
    cut = column.c/a.H
    break_up = [.4,a.h/a.H,.75,1.]
    break_low = [0.,cut,.4]

    def rhs(y,v):
        po.require(time.monotonic()<column.deadline,'two-ended audit exceeded bound')
        k = a.Kin*v[0]
        po.require(k>0,'independent shot lost positivity')
        r = y*a.H
        aux = float(a.aux(a.Kin,r))
        E = a.A/max(r,column.c)
        sigma = np.clip((y-.4)/.35,0.,1.)
        chi = sigma**3*(10.-15.*sigma+6.*sigma*sigma)
        alpha = a.nu+(1.-chi)*aux
        if chi==0.:
            G = T/alpha
            turb = aux*G
        else:
            def residual(q):
                return (alpha*E/k)*q+chi*k*shear_law(q)[0]-T
            po.require(residual(Q_PEAK)>=0,'companion-root shot left increasing shear')
            q = brentq(residual,0.,Q_PEAK,xtol=1e-13)
            G = (E/k)*q
            turb = (1.-chi)*aux*G+chi*k*shear_law(q)[0]
        return [a.H*a.A*v[1]/(a.Kin*(a.nu+a.cd*k*k/(a.sk*E))),
                a.H*(E-turb*G)/a.A,a.H*G]

    def propagate(points,start):
        polynomials = []
        for lo,hi in zip(points[:-1],points[1:]):
            sol = solve_ivp(rhs,(lo,hi),start,method='RK45',rtol=1e-12,atol=4e-14,max_step=.0125,dense_output=True)
            po.require(sol.success and sol.t[-1]==hi,'independent split IVP did not reach the match')
            polynomials.append((min(lo,hi),max(lo,hi),sol.sol))
            start = sol.y[:,-1]
        return start,polynomials

    def pairs(v):
        return propagate(break_low,[v[0],0.,0.]),propagate(break_up[::-1],[1.,v[1],0.])

    def residual(v):
        lower,upper = pairs(v)
        return (lower[0]-upper[0])[:2]

    initial = np.array([k0/a.Kin,qh])
    first = residual(initial)
    if np.max(abs(first))>2e-11:
        solved = root(residual,initial,tol=1e-9)
        po.require(np.max(abs(solved.fun))<=1e-10,'two-ended state matching failed')
        x = solved.x
    else:
        # The independently integrated endpoint defects already pass.
        x = initial
    (lower,lower_polys),(upper,upper_polys) = pairs(x)
    U_H = float(lower[2]-upper[2])
    po.require(abs(U_H-near_anchor)<=1e-9 and np.max(abs((lower-upper)[:2]))<=1e-10,
               'fixed-shear and two-ended BVPs disagree')
    errors = np.zeros(5)
    for shift, panels in ((0.,lower_polys),(U_H,upper_polys)):
        for lo,hi,coefficients in panels:
            for y in np.linspace(lo,hi,41):
                r = y*a.H
                current = coefficients(y)
                prior = primary_sample(r)
                dg = rhs(y,current)
                G,E,D = column.shear(r,prior[0],T)[:3]
                delta = [current[0]-prior[0],current[1]-prior[1],
                         current[2]+shift-prior[2]-near_anchor,dg[2]/a.H-G,
                         dg[0]-a.H*a.A*prior[1]/(a.Kin*D)]
                errors = np.maximum(errors,abs(np.asarray(delta)))
    po.require(np.max(errors)<=1e-8,'two-ended profile and derivative refinement differ')
    return dict(T=T,U_H=U_H,sampled_ground_k=float(x[0]*a.Kin),Q_H_over_A=float(x[1]),
                match_residual=(lower-upper)[:2].tolist(),difference_in_U_H=U_H-near_anchor,
                fixed_scaled_profile_difference=dict(zip(('k','Q','U','U_r','k_r'),errors.tolist())))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--table',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--seconds',type=float,required=True)
    p.add_argument('--cached-sweep',type=Path)
    p.add_argument('--cached-sweep-sha256')
    args = p.parse_args()
    po.require(sys.platform.startswith('linux') and args.seconds>0 and
               po.identity(args.table)['sha256']==po.TABLE_SHA,'requires confirmed Linux calibration')
    a = po.Reference(read(args.table))
    column = Column(a,time.monotonic()+args.seconds)
    args.out.mkdir()
    EH = a.A/a.H
    T_cap = .3*a.Kin+a.nu*(EH/a.Kin)*Q_PEAK
    t_near = T_cap*(1.-1e-10)/a.tau
    target = a.utau/a.kappa*math.log(a.H/a.z0)
    cached_identity = None
    sweep = []
    if args.cached_sweep:
        cached_identity = po.identity(args.cached_sweep)
        po.require(cached_identity['sha256']==args.cached_sweep_sha256,'stored sweep binding differs')
        sweep = read(args.cached_sweep)
        po.require(len(sweep)==27 and sweep[0]['T_over_tau']==t_near,
                   'stored sweep limit or coverage differs')
    approach = []
    guess = -.2226
    try:
        samples = [t_near,*np.linspace(.985,t_near-.0005,17)[::-1],.97,.94,.90,.85,.75,.6,.4,.2,.05]
        if not sweep:
            for index,t in enumerate(samples):
                result = fixed_shear(column,float(t),guess,tight=index==0)
                sweep.append(result)
                if len(result['solutions'])==1:
                    guess = result['solutions'][0]['Q_H_over_A']
        guess = sweep[0]['solutions'][0]['Q_H_over_A']
        for gap in (1e-5,1e-7,1e-9,1e-12):
            t = T_cap*(1-gap)/a.tau
            approach.append(fixed_shear(column,t,guess,tight=True))
        edge = approach[-1]
        po.require(len(edge['solutions'])==1,'fixed-shear witness bracket at H is unresolved')
        item = edge['solutions'][0]
        qh = item['Q_H_over_A']
        x = np.array([edge['T_over_tau'],qh])
        terminal,sample = column.integrate(x,tight=True)
        # Enforce k(H) directly in the independent upper shot. At a fractional
        # distance of 1e-12 from the cap, a generic forward RK stage can cross
        # the constitutive limit even when the boundary-value solution does not.
        metrics,_ = independent(column,x,terminal,sample,item['U_H'],
                                require_admissible=False,forward_control=False)
        check = two_ended(column,edge['T_over_tau'],item['k_ground'],qh,item['U_H'],sample)
        tracked = [v for row in sweep+approach for v in row['solutions'] if 'U_H' in v]
        po.require(item['U_H']<target and all(v['U_H']<target for v in tracked),
                   'fixed-shear scan found a possible terminal inlet')
        report = dict(schema='WJH27MatchedInletBranchAssessmentV1',status='BLOCKED',
                      table=po.identity(args.table),cutoff=column.c,
                      exact_shear_ceiling_at_H=T_cap,maximum_q_at_H=Q_PEAK,
                      target_U_H=target,last_column=item,anchor_gap=target-item['U_H'],
                      highest_sampled_U_H=max(v['U_H'] for v in tracked),
                      minimum_sampled_UH_sensitivity=min(v['dUH_d_T_over_tau'] for v in tracked),
                      sampled_shears=len(sweep)+len(approach),sampled_roots=len(tracked),checked_branch=metrics,
                      cached_sweep=cached_identity,
                      edge_approach=[r['solutions'] for r in approach],
                      two_ended=check,
                      interpretation='Canonical increasing-shear branch reached the H traction bound with a positive velocity gap. The fixed-shear scan and two-ended solve test this obstruction without a global-uniqueness claim.')
        atomic_out(args.out/'branch-assessment.json',report)
    finally:
        atomic_out(args.out/'fixed-shear-trials.json',dict(scan=sweep if not cached_identity else [],
                                                          cached_sweep=cached_identity,edge_approach=approach))
    print('ACTIVE_V2_INLET_BRANCH_OBSTRUCTION_CHECKED',flush=True)


if __name__=='__main__':
    main()
