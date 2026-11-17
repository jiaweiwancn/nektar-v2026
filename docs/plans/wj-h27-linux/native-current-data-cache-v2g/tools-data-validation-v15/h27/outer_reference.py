#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded finite-nu stationary outer-column feasibility test; Linux only.

This is a local parallel outer BVP with fixed inlet/interface scalar values
and the retained top k/epsilon fluxes. It is not a 2D developed reference.
No states, roots, stresses or scientific parameters are clipped or floored.
"""
import argparse
import json
import math
from dataclasses import asdict
from pathlib import Path
import time

import numpy as np
from scipy.integrate import quad, solve_bvp
from scipy.optimize import brentq

from reference_oracle import Column, Parameters, CP, QE, Q_PEAK, KAPPA


def shear(q):
    """Vector physical N root, Newton from its rigorous upper bracket."""
    q = np.asarray(q, dtype=float)
    if np.any(~np.isfinite(q)) or np.any(q < 0) or np.any(q > Q_PEAK*(1+1e-14)):
        raise ValueError('q outside declared first shear branch')
    n = CP+math.sqrt(1.35)*q
    for _ in range(12):
        f = (n-CP)*(n*n+q*q)-1.35*n*q*q
        fn = 3*n*n-2*CP*n-.35*q*q
        n = n-f/fn
    residual = (n-CP)*(n*n+q*q)-1.35*n*q*q
    scale = CP**3+n**3+q*q*n
    if np.any(n < CP*(1-1e-13)) or np.any(abs(residual)/scale > 2e-13):
        raise ValueError('physical N Newton residual or branch failure')
    return .6*n*q/(n*n+q*q), n


class Outer:
    def __init__(self, parameters, top, deadline):
        self.p = parameters
        self.column = Column(parameters)
        self.K, self.A = parameters.K, self.column.A
        self.t = self.K/QE
        self.rtop = top
        self.B = parameters.cd*self.K*self.K/1.3
        self.L = math.log(top/parameters.H)
        self.deadline = deadline
        if not top > parameters.H:
            raise ValueError('outer top must exceed matching height')

    def check_time(self):
        if time.monotonic() >= self.deadline:
            raise TimeoutError('whole-program elapsed diagnostic allowance exhausted')

    def constitutive(self, x, values, nu):
        self.check_time()
        p = self.p
        if np.any(~np.isfinite(values)) or np.any(values[0] <= 0) or np.any(values[2] <= 0):
            raise ValueError('nonpositive/nonfinite BVP iterate; no state clamp')
        radius = p.H*np.exp(x)
        k = self.K*values[0]
        epsilon = self.A*values[2]/radius
        def traction(q):
            return nu*epsilon*q/k+k*shear(q)[0]
        low = np.zeros_like(radius)
        high = np.full_like(radius, Q_PEAK)
        if np.any(traction(high) < self.t*(1-1e-13)):
            raise ValueError('no total-traction root on q<=q_peak')
        for _ in range(50):
            middle = (low+high)/2
            left = traction(middle) < self.t
            low = np.where(left, middle, low)
            high = np.where(left, high, middle)
        q = (low+high)/2
        shear_value, n = shear(q)
        gradient = epsilon*q/k
        production = k*shear_value*gradient
        dk = nu+p.cd*k*k/(p.sigma_k*epsilon)
        de = nu+p.cd*k*k/(1.3*epsilon)
        traction_error = nu*gradient+k*shear_value-self.t
        if np.any(abs(traction_error)/self.t > 3e-13):
            raise ValueError('bounded traction solve residual')
        return radius, k, epsilon, q, n, gradient, production, dk, de

    def rhs(self, x, values, nu):
        radius, k, epsilon, q, n, gradient, production, dk, de = self.constitutive(x, values, nu)
        pscaled = radius*production/self.A
        result = np.empty_like(values)
        result[0] = radius*self.A*values[1]/(self.K*dk)
        result[1] = values[2]-pscaled
        result[2] = values[2]+self.B*radius*values[3]/(de*self.A)
        result[3] = (values[3]+self.A*self.A/(self.B*self.K)
                     *(1.82*values[2]-1.44*pscaled)*values[2]/values[0])
        return result

    def bc(self, left, right, nu):
        return np.array([left[0]-1, left[2]-1, right[1],
                         right[3]+1+nu*self.A/(self.B*self.rtop)])

    def solve(self, nodes, tolerance):
        grid = np.linspace(0, self.L, nodes)
        initial = np.tile(np.array([[1.], [0.], [1.], [-1.]]), (1, nodes))
        attempts = []
        self.last_attempts = attempts
        current = None
        fraction = 0.0
        step = .125
        while True:
            target = 0.0 if current is None else min(1.0, fraction+step)
            nu = target*self.p.nu
            started = time.monotonic()
            try:
                self.check_time()
                guess = initial if current is None else current.sol(grid)
                solution = solve_bvp(lambda x, y: self.rhs(x, y, nu),
                                     lambda a, b: self.bc(a, b, nu),
                                     grid, guess, tol=tolerance, bc_tol=tolerance*.1,
                                     max_nodes=8001, verbose=0)
                if not solution.success:
                    raise ValueError('solve_bvp: '+solution.message)
                bc_error = float(np.max(abs(self.bc(solution.y[:, 0], solution.y[:, -1], nu))))
                values = self.constitutive(solution.x, solution.y, nu)
                attempts.append(dict(nu_fraction=target, success=True,
                                     seconds=time.monotonic()-started, nodes=len(solution.x),
                                     iterations=solution.niter, bc_max=bc_error,
                                     max_collocation_residual=float(np.max(solution.rms_residuals)),
                                     QH_over_A=float(solution.y[1, 0]),
                                     minimum_k_over_K=float(np.min(solution.y[0])),
                                     maximum_q=float(np.max(values[3]))))
                current, fraction = solution, target
                if fraction == 1.0:
                    break
                step = min(.25, step*1.5)
            except TimeoutError as error:
                attempts.append(dict(nu_fraction=target, success=False,
                                     seconds=time.monotonic()-started, reason=str(error)))
                return dict(status='DIAGNOSTIC_TIME_ALLOWANCE', initial_nodes=nodes,
                            tolerance=tolerance, attempts=attempts)
            except (ValueError, RuntimeError, FloatingPointError) as error:
                attempts.append(dict(nu_fraction=target, success=False,
                                     seconds=time.monotonic()-started, reason=str(error)))
                progress_threshold = max(64*np.finfo(float).eps*max(1., abs(fraction)), tolerance*.1)
                if current is None or step <= progress_threshold or target == fraction:
                    return dict(status='REJECTED', initial_nodes=nodes, tolerance=tolerance,
                                continuation_stop_reason=('initial high-Re BVP failed' if current is None
                                    else 'fractional nu step is at or below tol/10 or floating-point progress scale'),
                                minimum_distinguishable_nu_fraction_step=progress_threshold,
                                attempts=attempts)
                step /= 2
                continue
        probe = np.linspace(0, self.L, 129)
        values = current.sol(probe)
        radius, k, epsilon, q, n, gradient, production, dk, de = self.constitutive(probe, values, self.p.nu)
        derivative = current.sol(probe, 1)
        independent_ode = self.rhs(probe, values, self.p.nu)
        def balance(x):
            y = current.sol(np.array([x]))
            rr, kk, ee, qq, nn, gg, pp, dkk, dee = self.constitutive(np.array([x]), y, self.p.nu)
            return float(rr[0]*(pp[0]-ee[0])/self.A)
        integral, error = quad(balance, 0, self.L, epsabs=2e-10, epsrel=2e-9, limit=100)
        qh = float(current.y[1, 0])
        samples = [dict(x=float(probe[i]), r=float(radius[i]), k=float(k[i]), epsilon=float(epsilon[i]),
                        q=float(q[i]), G=float(gradient[i]), P=float(production[i]),
                        Qk_over_A=float(values[1, i]), scaled_epsilon_flux=float(values[3, i]))
                   for i in range(0, len(probe), 8)]
        output = dict(status='SOLVED', initial_nodes=nodes, tolerance=tolerance, attempts=attempts,
                      final_nodes=len(current.x), QH_over_A=qh, QH=qh*self.A,
                      integrated_P_minus_E_over_A=float(integral),
                      endpoint_budget_error=qh-float(integral), quadrature_error_estimate=float(error),
                      maximum_sampled_ode_defect=float(np.max(abs(derivative-independent_ode))),
                      min_k_over_K=float(np.min(values[0])), max_k_over_K=float(np.max(values[0])),
                      min_q=float(np.min(q)), max_q=float(np.max(q)),
                      epsilon_H=float(epsilon[0]), epsilon_trace_relative_error=float(epsilon[0]*self.p.H/self.A-1),
                      top_k_flux_over_A=float(values[1, -1]),
                      top_epsilon_flux=float(values[3, -1]*self.B/self.rtop),
                      prescribed_top_epsilon_flux=-self.B/self.rtop-self.p.nu*self.A/self.rtop**2,
                      samples=samples)
        return output


def inner_followup(parameters, qh, check_time):
    from dataclasses import replace
    p = replace(parameters, endpoint_flux_over_A=qh)
    column = Column(p)
    # A local search is deliberately not a global uniqueness assertion.
    center = math.log(.19696577127569226)-qh
    lower, upper = max(-4., center-.45), min(-.05, center+.45)
    if lower >= upper:
        return dict(status='LOCAL_CUTOFF_SEARCH_OUT_OF_RANGE', QH_over_A=qh)
    scans, roots = [], []
    for lam in np.linspace(lower, upper, 13):
        try:
            check_time()
        except TimeoutError as error:
            return dict(status='DIAGNOSTIC_TIME_ALLOWANCE', reason=str(error),
                        QH_over_A=qh, scans=scans, roots=roots)
        try:
            trial = column.solve_trial(float(lam))
            scans.append(dict(logc=float(lam), residual=trial['ground_flux']))
        except TimeoutError as error:
            return dict(status='DIAGNOSTIC_TIME_ALLOWANCE', reason=str(error),
                        QH_over_A=qh, scans=scans, roots=roots)
        except (ValueError, RuntimeError, OverflowError) as error:
            scans.append(dict(logc=float(lam), failure=str(error)))
    for left, right in zip(scans[:-1], scans[1:]):
        if 'failure' in left or 'failure' in right or left['residual']*right['residual'] > 0:
            continue
        try:
            def residual(value):
                check_time()
                return column.solve_trial(value)['ground_flux']
            lam = brentq(residual,
                         left['logc'], right['logc'], xtol=2e-12, rtol=1e-13)
            trial = column.solve_trial(lam)
            production = 0.0
            for lo, hi in zip(trial['breaks'][:-1], trial['breaks'][1:]):
                production += quad(lambda y: column.state(y, trial['sample'](y)[0], lam)['production'],
                                   lo, hi, epsabs=column.A*2e-12, epsrel=2e-11, limit=100)[0]*p.H
            ground = column.state(0, trial['sample'](0)[0], lam)
            roots.append(dict(c=trial['c'], c_over_H=trial['c']/p.H,
                              ground_flux_over_A=trial['ground_flux'],
                              budget_error_over_A=1-lam-production/column.A-qh,
                              derivative=trial['ground_flux_derivative'], anchor_U=trial['anchor'],
                              ground=ground))
        except TimeoutError as error:
            return dict(status='DIAGNOSTIC_TIME_ALLOWANCE', reason=str(error),
                        QH_over_A=qh, scans=scans, roots=roots)
        except (ValueError, RuntimeError, OverflowError) as error:
            roots.append(dict(bracket=[left['logc'], right['logc']], failure=str(error)))
    d = p.nu+column.auxiliary(0)
    ap = column.a*.4*p.H
    bound = math.log1p(ap/d)+p.nu/(d+ap)-p.nu/d+math.log(2.5)
    gamma = math.exp(1-bound)
    return dict(status='LOCAL_ROOTS_RECORDED', QH_over_A=qh, parameters=asdict(p),
                scans=scans, roots=roots,
                nonpositive_flux_collar_bound=dict(applicable=qh <= 0,
                    gamma_min=gamma, tail_branch_cutoff_lower_bound_over_H=gamma*math.exp(-qh),
                    whole_family_cutoff_lower_bound_over_H=min(.4, gamma*math.exp(-qh)),
                    caveat='Requires positive k diffusion and first shear branch; tail derivation uses c<.4H.'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--top', type=float, default=20.-.0875)
    parser.add_argument('--skip-inner-root', action='store_true')
    parser.add_argument('--diagnostic-seconds', type=float, default=140.)
    args = parser.parse_args()
    started = time.monotonic()
    if not math.isfinite(args.diagnostic_seconds) or args.diagnostic_seconds <= 0:
        parser.error('positive finite elapsed diagnostic allowance required')
    parameters = Parameters()
    problem = Outer(parameters, args.top, started+args.diagnostic_seconds)
    runs = []
    for nodes, tolerance in [(101, 2e-7), (201, 2e-9)]:
        try:
            runs.append(problem.solve(nodes, tolerance))
        except TimeoutError as error:
            runs.append(dict(status='DIAGNOSTIC_TIME_ALLOWANCE', initial_nodes=nodes,
                             tolerance=tolerance, reason=str(error), attempts=problem.last_attempts))
    payload = dict(schema='WJH27FiniteNuOuterReferenceV1', parameters=asdict(parameters),
                   scope='Stationary parallel outer reference only; not the same-inlet 2D developed flow',
                   boundary_data=dict(H=parameters.H, top_model_radius=args.top, kH=parameters.K,
                       epsilonH=problem.A/parameters.H, total_traction=problem.t, top_k_gradient_flux=0,
                       top_epsilon_gradient_flux=-problem.B/args.top-parameters.nu*problem.A/args.top**2),
                   variables=['k/K', 'Qk/A', 'r*epsilon/A', 'r*De*epsilon_prime/B'],
                   diagnostic_seconds=args.diagnostic_seconds, runs=runs)
    if all(run['status'] == 'SOLVED' for run in runs):
        payload['QH_over_A_refinement_difference'] = runs[1]['QH_over_A']-runs[0]['QH_over_A']
        if not args.skip_inner_root:
            try:
                payload['inner_followup'] = inner_followup(parameters, runs[1]['QH_over_A'], problem.check_time)
            except TimeoutError as error:
                payload['inner_followup'] = dict(status='DIAGNOSTIC_TIME_ALLOWANCE', reason=str(error))
    payload['seconds'] = time.monotonic()-started
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(payload, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(output=str(args.output), statuses=[run['status'] for run in runs],
                          QH_over_A=[run.get('QH_over_A') for run in runs], seconds=payload['seconds'])))


if __name__ == '__main__':
    main()
