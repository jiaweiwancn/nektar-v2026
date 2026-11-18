#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Solve and audit four fixed-choice shear columns with free Q(H) and U(H)."""

import argparse
import json
import math
from pathlib import Path
import sys
import time

from inner_k_mechanism import ROWS, TOL, read
from physical_oracle import Reference, TABLE_SHA, full_wj, identity, require
from reference_oracle import Column, Parameters
import mpmath as mp
import numpy as np
from numpy.polynomial import Polynomial
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


FACTORS = (.90, .925, .95, .975, 1., 1.025, 1.05, 1.075, 1.10)
Q_SCAN = np.linspace(-.4, .4, 17)
DERIVATIVE_LIMIT = 1e-8


def historical(ref, K, r):
    require(math.isfinite(K) and ref.lo <= K <= ref.hi, 'historical column K outside qualified interval')
    root = math.sqrt(K)
    X = ref.h*root/(70*ref.nu)
    logX = math.log(X)
    m1, m2 = .5-.4*X**.7, 1-.79*X**(-.28)
    require(X > 1 and m2 > max(m1, 0.), 'historical negative roughness branch or derivative changed')
    exponent = m2*logX
    multiplier = 10.7*ref.kappa/math.sqrt(ref.qe)
    a = ref.kappa*math.sqrt(K/ref.qe)
    d = ref.nu*(1+multiplier*math.expm1(exponent))
    dK = ref.nu*multiplier*math.exp(exponent)*(m2+.2212*X**(-.28)*logX)/(2*K)
    B = d+a*r
    BK = dK+a*r/(2*K)
    J = math.log1p(a*r/d)+ref.nu/B-ref.nu/d
    JK = BK*(B-ref.nu)/B**2-dK*(d-ref.nu)/d**2
    return dict(a=a, d=d, D=B, D_K=BK, m2_minus_m1=m2-m1, X=X, J=J, J_K=JK)


def cutoff(ref, K):
    row = historical(ref, K, ref.H)
    c = ref.H*math.exp(1-row['J'])
    require(0 < c < .4*ref.H, 'historical cutoff beyond lower layer')
    return c, -c*row['J_K']


class ShootingColumn:
    def __init__(self, reference, K, law, deadline, tight=False):
        self.ref = reference
        self.core = Column(Parameters(K=K, H=reference.H, h=reference.h, nu=reference.nu,
                                      cd=reference.cd, sigma_k=reference.sk, inner_diffusion=law))
        require(self.core.zv < 0 and reference.sk == 1., 'unqualified column branch or diffusion factor')
        self.deadline = deadline
        self.rtol, self.step = (5e-13, .0125) if tight else (2e-11, .025)
        self.cache = {}

    def state(self, y, k_ratio, c):
        if time.monotonic() >= self.deadline:
            raise TimeoutError('column exceeded its cooperative allowance')
        return self.core.state(float(y), float(k_ratio), math.log(c/self.ref.H))

    def integrate(self, qh, c, breaks):
        p, amp, us = self.core.p, self.core.A, self.core.us
        panels = []
        # Four and five carry derivatives of k/K and Q/A with respect to Q(H)/A.
        v = np.array([1., qh, 0., 0., 1.])

        def rhs(y, v):
            s = self.state(y, v[0], c)
            g, q, chi, tangent = (s[key] for key in ('g', 'q', 'chi', 'shear_tangent'))
            require(tangent > 0, 'selected shear Jacobian is not positive')
            dg = -chi*p.K*(s['shear_function']+q*s['shear_derivative'])*v[3]/tangent
            dD = 2*(s['diffusion']-p.nu)*v[3]/v[0] if p.inner_diffusion == 'scalar' else 0.
            denom = s['diffusion']
            return [p.H*amp*v[1]/(p.K*denom), p.H*(s['e']-s['production'])/amp, p.H*g/us,
                    p.H*amp*(v[4]*denom-v[1]*dD)/(p.K*denom**2),
                    -p.H*(self.core.t-2*p.nu*g)*dg/amp]

        for right, left in zip(breaks[:-1], breaks[1:]):
            out = solve_ivp(rhs, (right, left), v, method='DOP853', rtol=self.rtol,
                            atol=.04*self.rtol, max_step=self.step, dense_output=True)
            require(out.success and out.t[-1] == left, 'column IVP exited before panel end')
            panels.append((left, right, out.sol))
            v = out.y[:, -1]

        def sample(y):
            for lo, hi, dense in panels:
                if lo-1e-14 <= y <= hi+1e-14:
                    return dense(y)
            raise ValueError('radius outside solved column panels')
        return v, sample

    def ground(self, q, c):
        if q not in self.cache:
            try:
                # E=A/r for every admitted choice on these common panels.
                result, _ = self.integrate(q, .1*self.ref.H, [1., .75, self.ref.h/self.ref.H, .4])
                self.cache[q] = (float(result[1]), float(result[4]))
            except ValueError as error:
                self.cache[q] = str(error)
        value = self.cache[q]
        if isinstance(value, str):
            raise ValueError(value)
        lower_source = 1+math.log(.4*self.ref.H/c)-historical(self.ref, self.core.p.K, .4*self.ref.H)['J']
        return value[0]-lower_source, value[1]

    def shoot(self, c):
        sweep = []
        for q in Q_SCAN:
            try:
                f, df = self.ground(float(q), c)
                sweep.append(dict(q=float(q), ground=f, sensitivity=df))
            except ValueError as error:
                sweep.append(dict(q=float(q), rejection=str(error)))
        brackets = [(a['q'], b['q']) for a,b in zip(sweep[:-1], sweep[1:])
                    if 'ground' in a and 'ground' in b and a['ground']*b['ground'] <= 0]
        require(len(brackets) == 1, f'unresolved number of adjacent admissible root brackets: {brackets}')
        q = brentq(lambda q: self.ground(q,c)[0], *brackets[0], xtol=2e-14, rtol=2e-14)
        residual, slope = self.ground(q, c)
        supported_slopes = [x['sensitivity'] for x in sweep if 'sensitivity' in x]
        require(min(supported_slopes+[slope]) > 0 and abs(residual) <= TOL,
                'sampled column endpoint-root slope or ground residual failed')
        dq = 2e-5
        centered = sum(w*self.ground(q+j*dq,c)[0] for j,w in ((-2,1),(-1,-8),(1,8),(2,-1)))/(12*dq)
        require(abs(centered-slope) <= DERIVATIVE_LIMIT, 'independent endpoint-root sensitivity differs')
        breaks = sorted({0., c/self.ref.H, .4, self.ref.h/self.ref.H, .75, 1.}, reverse=True)
        v, sample = self.integrate(q, c, breaks)
        require(abs(v[1]) <= TOL and abs(v[4]-slope) <= DERIVATIVE_LIMIT, 'whole-column ground checks failed')
        return dict(QH_over_A=q, bracket=brackets[0], scan=sweep,
                    ground_Q_over_A=float(v[1]), sensitivity=slope,
                    minimum_scanned_sensitivity=min(supported_slopes), differentiated_sensitivity=centered,
                    UH=float(-v[2]*self.core.us)), breaks, sample


def stress(k, e, gradient, auxiliary, chi):
    n = np.array([0., 0., 1.])
    tangential = np.eye(3)-np.outer(n,n)
    M = gradient+gradient.transpose(0,2,1)
    shear = -auxiliary[:,None]*(M@n@tangential)
    R = (2*k/3)[:,None,None]*np.eye(3)+shear[:,:,None]*n+n[None,:,None]*shear[:,None,:]
    active = chi > 0
    if np.any(active):
        R[active] = (1-chi[active,None,None])*R[active]+chi[active,None,None]*full_wj(k[active], e[active], gradient[active])[0]
    return R


def planar_response(k, e, g, auxiliary, chi, nu):
    G = np.zeros((len(k),3,3))
    G[:,0,2] = g
    R = stress(k,e,G,auxiliary,chi)
    mu = nu-R[:,0,2]/g
    production = -R[:,0,2]*g
    derivative_rows = []
    for step in (1e-3, 5e-4):
        h = step*e/k
        derivatives = []
        for i,j in ((0,0), (0,2), (2,0)):
            direction = np.zeros((3,3))
            direction[i,j] = 1.
            if i == j:
                direction[2,2] = -1.
            dR = sum(w*stress(k,e,G+(m*h)[:,None,None]*direction,auxiliary,chi)
                     for m,w in ((-2,1), (-1,-8), (1,8), (2,-1)))/(12*h[:,None,None])
            derivatives.append(nu*(direction+direction.T)-dR)
        derivative_rows.append(np.asarray(derivatives)/mu[None,:,None,None])
    fine = derivative_rows[-1]
    err = float(np.max(np.abs(fine-derivative_rows[0])))
    # With s/c=t, the divergence-free rank-one variation has coordinates
    # (-t,-t^2,1)/(1+t^2). The stress contraction supplies a second quadratic.
    minima = []
    for d in fine.transpose(1,0,2,3):
        factors = [Polynomial([0.,-1.]), Polynomial([0.,0.,-1.]), Polynomial([1.])]
        p = sum(f*Polynomial([S[0,2], S[2,2]-S[0,0], -S[0,2]]) for f,S in zip(factors,d))
        coeff = np.pad(p.coef, (0,5-len(p.coef)))
        station = Polynomial((p.deriv()*Polynomial([1.,0.,1.])-p*Polynomial([0.,4.])).coef[:5])
        values = [coeff[0],coeff[4]]+[p(t.real)/(1+t.real*t.real)**2 for t in station.roots() if abs(t.imag) <= 1e-9]
        minima.append(min(values))
    return dict(minimum_covariance_over_k=float(np.min(np.linalg.eigvalsh(R/k[:,None,None]))),
                minimum_principal_over_mu=float(min(minima)), tangent_refinement_difference=err,
                radii_checked=len(k)), production


def measure(column, c, result, breaks, sample, order, count):
    ref = column.ref
    nodes, weights = leggauss(order)
    lower, upper = np.asarray(breaks[1:]), np.asarray(breaks[:-1])
    y = np.r_[((upper[:,None]+lower[:,None]+(upper-lower)[:,None]*nodes)/2).ravel(),
              np.unique(np.r_[np.linspace(0,1,count), breaks, ref.old_c/ref.H])]
    weight = np.r_[((upper-lower)[:,None]*weights/2).ravel(), np.zeros(len(y)-len(lower)*order)]
    v = np.asarray([sample(x) for x in y])
    k = v[:,0]*column.core.p.K
    U = result['UH']+column.core.us*v[:,2]
    old = ref.inlet(y*ref.H)
    terms = {key:[] for key in ('production','diffusion','shear_tangent','traction_error','e','g','auxiliary','chi')}
    for x, ratio in zip(y, v[:,0]):
        row = column.state(x, ratio, c)
        for key in terms:
            terms[key].append(row[key])
    f = {key:np.asarray(values) for key,values in terms.items()}
    # The independent derivative/guard grid includes both cutoff ends and all
    # junctions and is deliberately finite in wall-normal position.
    pick = np.unique(np.r_[np.searchsorted(y[-(len(y)-len(lower)*order):], np.linspace(0,1,129)) + len(lower)*order,
                           [int(np.argmin(abs(y-x))) for x in breaks]])
    guard, production = planar_response(*(a[pick] for a in (k, f['e'], f['g'], f['auxiliary'], f['chi'])), ref.nu)
    amp_scale = column.core.A/ref.H
    root_error = float(np.max(np.abs(production-f['production'][pick]))/amp_scale)
    error = float(np.max(abs(f['traction_error']))/column.core.t)
    production_integral = float(weight@f['production'])/amp_scale
    energy_integral = float(weight@f['e'])/amp_scale
    budget = energy_integral-production_integral-(result['QH_over_A']-result['ground_Q_over_A'])
    ea_exact = 1+math.log(ref.H/c)
    require(max(error, abs(budget), abs(ea_exact-energy_integral), root_error) <= TOL, 'independent column shear or budget differs')
    require(float(np.min(f['shear_tangent'])) > 0 and float(np.min(k)) > 0 and float(np.min(f['e'])) > 0 and
            guard['minimum_covariance_over_k'] > 1e-10 and guard['minimum_principal_over_mu'] > 1e-8 and
            guard['tangent_refinement_difference'] <= DERIVATIVE_LIMIT, 'sampled amended-v2 admissibility failed')
    norms = {}
    for name, delta in (('k_over_Kin', (k-old['k'])/ref.Kin), ('U', U-old['U'])):
        at = int(np.argmax(np.abs(delta)))
        norms[name] = dict(RMS=math.sqrt(float(weight@(delta*delta))), MAX=float(abs(delta[at])),
                           signed_peak=float(delta[at]), r=float(y[at]*ref.H))
    return dict(cutoff=c, quadrature=order, samples=len(y), minimum_k=float(np.min(k)), minimum_e=float(np.min(f['e'])),
                minimum_shear_tangent=float(np.min(f['shear_tangent'])), traction_error_over_tau=error,
                tensor_production_error_over_A_per_H=root_error, budget_error_over_A=budget,
                P_integral_over_A=production_integral, E_integral_over_A=energy_integral, integrated_E_error=ea_exact-energy_integral,
                admissibility=guard, norms=norms,
                ground_k=float(k[-(len(y)-len(lower)*order)]), **result), y


def derivative_audit(ref, K):
    c, dc = cutoff(ref, K)
    action = historical(ref, K, .31*ref.H)
    with mp.workdps(50):
        nu, H, h, kap, qe = map(mp.mpf, (ref.nu,ref.H,ref.h,ref.kappa,ref.qe))

        def dcoef(q, r):
            X = h*mp.sqrt(q)/(70*nu)
            v = (1-mp.mpf('.79')*mp.power(X,mp.mpf('-.28')))*mp.log(X)
            zv = -mp.mpf('10.7')*nu*mp.expm1(v)/mp.sqrt(q)
            return nu+kap*mp.sqrt(q/qe)*(r-zv)

        def historical_c(q):
            d, B = dcoef(q,0), dcoef(q,H)
            return H*mp.exp(1-mp.log(B/d)-nu/B+nu/d)
        values = [float(historical_c(K)), float(mp.diff(historical_c,K)),
                  float(dcoef(K,.31*H)), float(mp.diff(lambda q:dcoef(q,.31*H),K))]
    errors = [abs(c-values[0])/ref.H, abs(dc-values[1])*ref.Kin/ref.H,
              abs(action['D']-values[2])/action['D'], abs(action['D_K']-values[3])*ref.Kin/action['D']]
    require(max(errors) <= TOL and dc < 0, 'historical coefficient or derivative identity differs')
    return dict(cutoff=c, cutoff_dK=dc, normalized_action_errors=errors, **action)


def run(ref, deadline, report):
    groups, derivatives, fine = {}, [], {}
    report.update(groups=groups, derivatives=derivatives, kin_refined=fine)
    for factor in FACTORS:
        K = ref.lo if factor == FACTORS[0] else ref.hi if factor == FACTORS[-1] else factor*ref.Kin
        dc = derivative_audit(ref,K)
        derivatives.append(dict(K_over_Kin=factor, **dc))
        choices = dict(v2=float(ref.cutoff(K)), historical=dc['cutoff'])
        solvers = {law:ShootingColumn(ref,K,law,deadline) for law in ('scalar','auxiliary')}
        rows = {}
        groups[str(factor)] = rows
        for name,(law,choice) in ROWS.items():
            c = choices[choice]
            solver = solvers[law]
            root, breaks, sample = solver.shoot(c)
            row, _ = measure(solver,c,root,breaks,sample,48,257)
            if name == 'A':
                require(abs(row['QH_over_A']) <= TOL, 'v2 table no longer recovers its zero endpoint flux')
            if factor == 1.:
                refined = ShootingColumn(ref,K,law,deadline,tight=True)
                rr, rb, rs = refined.shoot(c)
                precise, y = measure(refined,c,rr,rb,rs,80,1025)
                a,b = (np.asarray([fun(x) for x in y]) for fun in (sample,rs))
                common_k = float(np.max(abs(a[:,0]-b[:,0])))
                common_u = float(np.max(abs(root['UH']-rr['UH']+solver.core.us*(a[:,2]-b[:,2]))))
                rms, max_change = (max(abs(precise['norms'][name][kind]-row['norms'][name][kind])
                                       for name in row['norms']) for kind in ('RMS','MAX'))
                require(max(common_k, common_u, rms) <= 2e-9 and max_change <= 2e-5 and
                        abs(root['QH_over_A']-rr['QH_over_A']) <= TOL, 'Kin column refinement unresolved')
                fine[name] = dict(result=precise, common_k_MAX=common_k, common_U_MAX=common_u,
                                  norm_RMS_change=rms, sampled_MAX_change=max_change)
            rows[name] = row
        print('COLUMN_K_FACTOR',factor, 'QH/A', [rows[name]['QH_over_A'] for name in ROWS], flush=True)
    invalid = []
    for K in (ref.lo*(1-1e-6), ref.hi*(1+1e-6), math.nan, 0.):
        try:
            cutoff(ref,K)
        except ValueError as error:
            invalid.append(dict(K=repr(K), rejection=str(error)))
        else:
            raise ValueError('historical cutoff accepted an invalid K')
    return groups, derivatives, fine, invalid


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('table','mechanism','out'):
        p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--seconds',type=float,required=True)
    args = p.parse_args()
    require(sys.platform.startswith('linux') and 0 < args.seconds <= 490 and not args.out.exists(),
            'create-only bounded Linux column execution required')
    started = time.monotonic()
    require(identity(args.table)['sha256'] == TABLE_SHA and read(args.mechanism)['MECHANISM_DIAGNOSIS'] == 'PASS',
            'qualified frozen F0 and v2 inputs required')
    report = dict(schema='WJInnerKFixedChoiceColumnsV1', TASK3_READY=False, table=identity(args.table),
                  mechanism=identity(args.mechanism), script=identity(__file__),
                  brackets=[float(Q_SCAN[0]),float(Q_SCAN[-1])], scope='Finite K/radius and branch checks; counterfactual columns')
    try:
        groups, derivatives, fine, invalid = run(Reference(read(args.table)), started+args.seconds, report)
        d = fine['D']['result']['norms']
        supported = d['k_over_Kin']['RMS'] <= .05 and d['k_over_Kin']['MAX'] <= .10 and d['U']['RMS'] <= .01 and d['U']['MAX'] <= .02
        report.update(groups=groups, derivatives=derivatives, kin_refined=fine, invalid_K_controls=invalid,
                      COLUMN_REMEDY='SUPPORTED' if supported else 'REJECTED')
    except (ValueError, TimeoutError) as error:
        report.update(COLUMN_REMEDY='UNRESOLVED', exception=str(error))
    report['elapsed_seconds'] = time.monotonic()-started
    with args.out.open('x') as stream:
        json.dump(report,stream,indent=2,allow_nan=False)
        stream.write('\n')
    print('COLUMN_REMEDY',report['COLUMN_REMEDY'],report.get('exception',''),flush=True)
    if 'exception' in report:
        raise SystemExit(2)


if __name__ == '__main__':
    main()
