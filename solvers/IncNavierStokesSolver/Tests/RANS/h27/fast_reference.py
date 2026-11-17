#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Generate and independently validate bounded H27-v2 reference tables on Linux.

Only counterfactual Q_ref(H)=0 calibration is tabulated. No instantaneous
actual-flow production budget is constrained. Geometry, inlet and tail
amplitude are fixed. A time allowance replaces a fixed number of retries.
"""
import os
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_key] = '1'

import argparse
from dataclasses import asdict, replace
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from numpy.polynomial import chebyshev as cheb
from scipy.integrate import quad, solve_ivp

import reference_oracle as oracle
from reference_oracle import Column, Parameters, QE, KAPPA, CL

NE = 81.0/20.0
K_LOW, K_HIGH = .9, 1.1
BUDGET_TOL = 1e-10
DERIVATIVE_TOL = 1e-8
TRACTION_TOL = 1e-10


def identity(path):
    path = Path(path)
    data = path.read_bytes()
    return dict(path=str(path), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


class Rejection(ValueError):
    pass


class Generator:
    def __init__(self, seconds):
        self.started = time.monotonic()
        self.deadline = self.started+seconds
        self.p = Parameters(rtol=2e-12)
        self.kin = self.p.K
        self.cache = {}
        self.full_cache = {}
        self.trials = []
        self.payload = dict(schema='WJH27V2ReferenceTablesV1', status='IN_PROGRESS',
            model='IsotropicLowerWallFullWJCollarV2',
            scope='Counterfactual zero-endpoint-k-flux reference only; never actual-flow budget forcing',
            certification='Finite held-out validation of a bounded interpolation; no continuum error or full principal claim',
            parameters=asdict(self.p),
            constants=dict(q_e=QE, N_e=NE, kappa=KAPPA, c_length=CL),
            K_interval=[K_LOW*self.kin, K_HIGH*self.kin], K_factor_interval=[K_LOW, K_HIGH],
            guard=dict(K_out_of_range='reject; no clamping', cutoff_required='0<c<.4H',
                       reference_endpoint_flux_over_A=0, sigma_e=1.3,
                       full_principal_qualification='Separate native v2 qualification; not claimed here'),
            acceptance=dict(full_IVP_ground_flux_over_A=BUDGET_TOL,
                            independent_budget_over_A=BUDGET_TOL,
                            K_in_times_logcutoff_derivative_error=DERIVATIVE_TOL,
                            outlet_traction_error_over_tau=TRACTION_TOL),
            coefficient_convention='Ascending Chebyshev T_j coefficients; x=(2*argument-lo-hi)/(hi-lo)',
            diagnostic_seconds=seconds, trials=self.trials,
            sources=dict(generator=identity(__file__), reference_oracle=identity(oracle.__file__)))
        self.bind_inlet()

    def check_time(self):
        if time.monotonic() >= self.deadline:
            raise TimeoutError('whole-generator elapsed diagnostic allowance exhausted')

    def column(self, factor):
        if not math.isfinite(factor) or not K_LOW <= factor <= K_HIGH:
            raise Rejection('reference K is outside declared interval')
        parent = self
        class TimedColumn(Column):
            def state(self, *args):
                parent.check_time()
                return super().state(*args)
        return TimedColumn(replace(self.p, K=self.kin*factor))

    def bind_inlet(self):
        repository = Path(__file__).resolve().parents[5]
        self.payload['sources']['amended_model_contract'] = identity(
            repository/'docs/plans/wj-h27-amended-candidate-v2.md')
        path = repository/'docs/plans/wj-h5-linux/calibration.json'
        candidates = json.loads(path.read_text())['candidates']
        old = next(row for row in candidates if row['H_w_mm'] == 6.)
        current = dict(K=self.p.K, H_w=self.p.H, h_s=self.p.h, nu=self.p.nu,
                       A=self.p.K**1.5/CL)
        errors = {key: abs(value-old[key])/abs(old[key]) for key, value in current.items()}
        epsilon_error = abs(current['A']/self.p.H-old['A']/old['H_w'])/(old['A']/old['H_w'])
        if max([epsilon_error]+list(errors.values())) > 1e-13:
            raise Rejection('fixed inlet/roughness/height/molecular-viscosity binding mismatch')
        self.rtop = old['top']-old['d0']
        self.tau = self.kin/QE
        self.payload['inlet_binding'] = dict(calibration=identity(path), relative_errors=errors,
            epsilon_corner_relative_error=epsilon_error, K_in=self.kin,
            A_in=current['A'], epsilon_in_H=current['A']/self.p.H,
            u_tau=old['u_tau'], d0=old['d0'], top_physical=old['top'], top_model_radius=self.rtop,
            Hw_mm=old['H_w_mm'], physical_interface_mm=old['physical_z_H_mm'],
            roughness=old['h_s'], roughness_mm=old['h_s_mm'],
            inlet_scope='Retained full-height high-Re-derived application inlet; not replaced by reference anchor')

    def fast(self, factor):
        self.check_time()
        key = float(factor).hex()
        if key in self.cache:
            return self.cache[key]
        began = time.monotonic()
        record = dict(kind='fast_collar', K_factor=float(factor))
        try:
            column = self.column(float(factor))
            if not column.zv < 0:
                raise Rejection('reference is outside selected negative-intercept rough branch')
            p, A = column.p, column.A
            value = np.array([1., 0., 0., 0.])
            panels = []
            def rhs(y, state):
                v = column.state(float(y), float(state[0]), math.log(.1))
                return [p.H*A*state[1]/(p.K*v['diffusion']),
                        p.H*(v['e']-v['production'])/A,
                        p.H*v['g']/column.us, p.H*v['production']/A]
            for hi, lo in [(1., .75), (.75, .4)]:
                solution = solve_ivp(rhs, (hi, lo), value, method='DOP853',
                                     rtol=2e-13, atol=8e-15, max_step=.02, dense_output=True)
                if not solution.success:
                    raise Rejection('fast collar integration: '+solution.message)
                value = solution.y[:, -1]
                panels.append((lo, hi, solution))
            def sample(y):
                for lo, hi, solution in panels:
                    if lo-1e-14 <= y <= hi+1e-14:
                        return solution.sol(y)
                raise Rejection('fast collar sample outside [.4,1]')
            d = p.nu+column.auxiliary(0)
            span = column.a*.4*p.H
            production_lower = math.log1p(span/d)+p.nu/(d+span)-p.nu/d
            production_total = production_lower-float(value[3])
            logc = 1-production_total  # Q_ref(H)/A is exactly zero.
            cutoff = p.H*math.exp(logc)
            if not 0 < cutoff < .4*p.H:
                raise Rejection('derived cutoff leaves strict c<.4H branch: c/H='+str(cutoff/p.H))
            rows = [column.state(float(y), float(sample(y)[0]), logc)
                    for y in np.linspace(.4, 1, 65)]
            if (min(row['k'] for row in rows) <= 0 or min(row['g'] for row in rows) <= 0
                    or min(row['shear_tangent'] for row in rows) <= 0):
                raise Rejection('nonpositive collar k or first-branch shear tangent')
            anchor = column.t/column.a*math.log1p(span/d)-float(value[2])*column.us
            result = dict(column=column, sample=sample, logc=logc, cutoff=cutoff,
                          production_over_A=production_total, anchor=anchor)
            self.cache[key] = result
            record.update(status='PASS', c_over_H=cutoff/p.H,
                          min_collar_k_over_K=min(row['k']/p.K for row in rows),
                          minimum_collar_shear=min(row['g'] for row in rows),
                          max_collar_q=max(row['q'] for row in rows),
                          minimum_shear_tangent=min(row['shear_tangent'] for row in rows))
            return result
        except Exception as error:
            record.update(status='REJECTED', reason=str(error))
            raise
        finally:
            record['seconds'] = time.monotonic()-began
            self.trials.append(record)

    def fit_logcutoff(self, degree):
        nodes = np.cos(np.pi*(np.arange(degree+1)+.5)/(degree+1))
        values = [self.fast(float(1+.1*x))['logc'] for x in nodes]
        return cheb.chebfit(nodes, values, degree)

    def full_validation(self, factor, logc):
        self.check_time()
        record = dict(kind='complete_IVP_and_independent_GK', K_factor=float(factor), logc=float(logc))
        began = time.monotonic()
        try:
            column = self.column(float(factor))
            trial = column.solve_trial(float(logc))
            integrals = [0., 0.]
            for lo, hi in zip(trial['breaks'][:-1], trial['breaks'][1:]):
                for j, name in enumerate(('production', 'e')):
                    integrals[j] += column.p.H*quad(
                        lambda y: column.state(float(y), float(trial['sample'](y)[0]), float(logc))[name],
                        lo, hi, epsabs=column.A*2e-12, epsrel=2e-11, limit=100)[0]
            probes = sorted(set(np.linspace(0, 1, 129).tolist()+trial['breaks']))
            states = [column.state(y, float(trial['sample'](y)[0]), float(logc)) for y in probes]
            row = dict(K_factor=float(factor), logc=float(logc), c_over_H=math.exp(logc),
                       ground_flux_over_A=trial['ground_flux'],
                       independent_budget_over_A=(integrals[1]-integrals[0])/column.A,
                       independent_sink_identity_error=integrals[1]/column.A-(1-logc),
                       complete_IVP_logcutoff_slope=trial['ground_flux_derivative'],
                       min_k_over_K=min(v['k']/column.p.K for v in states),
                       minimum_shear=min(v['g'] for v in states),
                       minimum_shear_tangent=min(v['shear_tangent'] for v in states),
                       epsilon_corner_relative_error=abs(column.state(1, 1., logc)['e']*column.p.H/column.A-1))
            if (row['min_k_over_K'] <= 0 or row['minimum_shear'] <= 0
                    or row['minimum_shear_tangent'] <= 0):
                raise Rejection('complete reference column violates positive k/shear/first-branch tangent')
            self.full_cache[float(factor).hex()] = (column, trial)
            record.update(status='RECORDED', **{k: v for k, v in row.items() if k != 'K_factor'})
            return row
        except Exception as error:
            record.update(status='REJECTED', reason=str(error))
            raise
        finally:
            record['seconds'] = time.monotonic()-began
            self.trials.append(record)

    def derivative_check(self, coefficients, factor):
        derivative = cheb.chebder(coefficients)/.1  # derivative w.r.t. K/K_in
        predicted = float(cheb.chebval((factor-1)/.1, derivative))
        values = []
        for step in (.001, .0005):
            if factor-2*step >= K_LOW and factor+2*step <= K_HIGH:
                fd = (self.fast(factor-2*step)['logc']-8*self.fast(factor-step)['logc']
                      +8*self.fast(factor+step)['logc']-self.fast(factor+2*step)['logc'])/(12*step)
            else:
                sign = 1. if factor < 1 else -1.
                samples = [self.fast(float(factor+sign*j*step))['logc'] for j in range(5)]
                fd = sign*np.dot([-25, 48, -36, 16, -3], samples)/(12*step)
            values.append(dict(step_in_K_factor=step, FD_derivative=float(fd),
                               K_in_times_derivative_error=abs(float(fd)-predicted)))
        return dict(K_factor=factor, derivative_with_respect_to_K_factor=predicted, checks=values)

    def cutoff_table(self):
        heldout = np.linspace(K_LOW, K_HIGH, 33)
        heldout[16] = 1.  # Pin the inlet case used by independent outlet checks.
        coarse = self.fit_logcutoff(14)
        degree = 20
        convergence = []
        self.payload['cutoff_table'] = dict(refinement=convergence)
        while True:
            self.check_time()
            coefficients = self.fit_logcutoff(degree)
            direct_errors = [abs(float(cheb.chebval((float(f)-1)/.1, coefficients))-self.fast(float(f))['logc'])
                             for f in heldout]
            coarse_errors = [abs(float(cheb.chebval((float(f)-1)/.1, coarse))-self.fast(float(f))['logc'])
                             for f in heldout]
            convergence.append(dict(degree=degree, maximum_logc_error=max(direct_errors)))
            if max(direct_errors) <= 1e-12:
                break
            degree += 8
        validations = []
        self.payload['cutoff_table']['validations'] = validations
        for factor in heldout:
            validations.append(self.full_validation(
                float(factor), float(cheb.chebval((float(factor)-1)/.1, coefficients))))
        worst = max(max(abs(row['ground_flux_over_A']), abs(row['independent_budget_over_A']))
                    for row in validations)
        derivatives = [self.derivative_check(coefficients, float(heldout[i]))
                       for i in (0, 5, 11, 16, 21, 27, 32)]
        derivative_error = max(check['K_in_times_derivative_error'] for row in derivatives for check in row['checks'])
        self.payload['cutoff_table'] = dict(argument='K', domain=[K_LOW*self.kin, K_HIGH*self.kin],
            value='log(c/H)', degree=degree, coefficients=coefficients.tolist(),
            logcutoff_dK_coefficients=(cheb.chebder(coefficients)/(.1*self.kin)).tolist(),
            coarse_degree=14, coarse_coefficients=coarse.tolist(),
            coarse_maximum_logc_error=max(coarse_errors), refinement=convergence,
            heldout_K_count=len(validations), validations=validations, derivative_validation=derivatives,
            maximum_budget_or_ground_flux_error=worst,
            maximum_K_in_times_derivative_error=derivative_error,
            exact_conditional_logcutoff_residual_slope=1.)
        if worst > BUDGET_TOL or derivative_error > DERIVATIVE_TOL:
            raise Rejection('cutoff table independent validation tolerance failed')

    def traction(self, r, complete=False):
        self.check_time()
        if not 0 <= r <= self.rtop:
            raise Rejection('outlet radius outside fixed physical domain')
        if r <= .4*self.p.H:
            return self.kin/NE
        y = r/self.p.H
        if r <= self.p.H:
            if complete:
                column, trial = self.full_cache[float(1.).hex()]
                ratio = float(trial['sample'](y)[0])
                state = column.state(y, ratio, trial['logc'])
            else:
                ref = self.fast(1.)
                state = ref['column'].state(y, float(ref['sample'](y)[0]), ref['logc'])
        else:
            # Full-WJ outer counterfactual with fixed k=K_in, E=A/r and total
            # shear tau0. This is not the stationary outer scalar BVP.
            state = self.fast(1.)['column'].state(y, 1., math.log(.1))
        diagonal = 3*state['q']**2/(5*(state['n']**2+state['q']**2))
        return self.kin/NE-2*state['chi']*state['k']*diagonal

    def outlet_table(self):
        segments = [dict(name='wall_transition', argument='r', domain=[.4*self.p.H, .75*self.p.H]),
                    dict(name='full_WJ_inner_collar', argument='r', domain=[.75*self.p.H, self.p.H]),
                    dict(name='outer', argument='log(r/H)', domain=[0., math.log(self.rtop/self.p.H)])]
        self.payload['outlet_traction'] = dict(normal=[1, 0, 0], component_order=['x', 'z'],
            pressure_reference='p/rho+R_zz= -K_in/N_e; gravity excluded from this reduced-pressure traction',
            formula='t_x=K_in/N_e-2*chi*k_ref*3*q^2/[5*(N^2+q^2)], t_z=tau0',
            lower_constant=dict(domain=[0., .4*self.p.H], tx=self.kin/NE, tz=self.tau),
            segments=segments,
            outer_scope='Counterfactual k=K_in, epsilon=A_in/r, finite-nu constant total shear; not stationary scalar balance')
        for segment in segments:
            lo, hi = segment['domain']
            def radius(argument):
                if segment['argument'] == 'r':
                    return argument
                # Exact arithmetic endpoints, not clamping a scientific state.
                if argument == lo:
                    return self.p.H
                if argument == hi:
                    return self.rtop
                return self.p.H*math.exp(argument)
            heldout = np.linspace(lo, hi, 129)
            degree, attempts = 24, []
            segment['refinement'] = attempts
            while True:
                self.check_time()
                attempt = dict(kind='outlet_Chebyshev', segment=segment['name'], degree=degree,
                               status='IN_PROGRESS')
                self.trials.append(attempt)
                nodes = np.cos(np.pi*(np.arange(degree+1)+.5)/(degree+1))
                values = [self.traction(radius(float((lo+hi)/2+(hi-lo)*x/2))) for x in nodes]
                coefficients = cheb.chebfit(nodes, values, degree)
                errors = []
                for argument in heldout:
                    estimate = float(cheb.chebval((2*argument-lo-hi)/(hi-lo), coefficients))
                    expected = self.traction(radius(float(argument)), complete=True)
                    errors.append(abs(estimate-expected)/self.tau)
                attempts.append(dict(degree=degree, maximum_tx_error_over_tau=max(errors)))
                attempt.update(status='PASS' if max(errors) <= TRACTION_TOL else 'REFINE',
                               maximum_tx_error_over_tau=max(errors))
                if max(errors) <= TRACTION_TOL:
                    break
                degree += 8
            segment.update(degree=degree, tx_coefficients=coefficients.tolist(),
                           tz_constant=self.tau, refinement=attempts,
                           heldout_count=len(heldout), maximum_tx_error_over_tau=max(errors),
                           endpoint_values=[self.traction(radius(lo), True), self.traction(radius(hi), True)])
        def at(segment, argument):
            lo, hi = segment['domain']
            return float(cheb.chebval((2*argument-lo-hi)/(hi-lo), segment['tx_coefficients']))
        continuity = dict(at_0p4H_over_tau=abs(at(segments[0], .4*self.p.H)-self.kin/NE)/self.tau,
            at_0p75H_over_tau=abs(at(segments[0], .75*self.p.H)-at(segments[1], .75*self.p.H))/self.tau,
            at_H_over_tau=abs(at(segments[1], self.p.H)-at(segments[2], 0.))/self.tau)
        self.payload['outlet_traction']['continuity'] = continuity
        if max(continuity.values()) > TRACTION_TOL:
            raise Rejection('outlet table one-sided endpoint continuity tolerance failed')

    def execute(self):
        self.cutoff_table()
        self.outlet_table()
        self.payload['status'] = 'PASS'
        self.payload['derivation'] = (
            'On 0<c<.4H, collar E=A/r and backward k,Q from H are independent of c; '
            'below .4H, P_aux and G are independent of k,E,c. The new full-WJ-collar '
            'budget gives log(c/H)=1-P_ref_budget/A-Q_ref(H)/A. Q_ref(H)=0 only '
            'for this counterfactual reference; actual storage/transport/flux remain unconstrained.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--diagnostic-seconds', type=float, default=280.)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('immutable output already exists')
    if not math.isfinite(args.diagnostic_seconds) or args.diagnostic_seconds <= 0:
        parser.error('positive finite diagnostic allowance required')
    generator = None
    try:
        generator = Generator(args.diagnostic_seconds)
        generator.execute()
        payload = generator.payload
    except Exception as error:
        payload = generator.payload if generator else dict(schema='WJH27V2ReferenceTablesV1')
        payload.update(status='NOT_MET', failure=dict(type=type(error).__name__, reason=str(error)))
        for trial in payload.get('trials', []):
            if trial.get('status') == 'IN_PROGRESS':
                trial.update(status='REJECTED', reason=str(error))
    if generator:
        payload['elapsed_seconds'] = time.monotonic()-generator.started
    encoded = json.dumps(payload, indent=2, allow_nan=False)+'\n'
    if len(encoded.encode()) >= 1024*1024:
        raise RuntimeError('bounded output exceeds 1 MiB; no truncated success artifact written')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        stream.write(encoded)
    print(json.dumps(dict(output=str(args.output), status=payload['status'],
                          seconds=payload.get('elapsed_seconds'), failure=payload.get('failure'))))
    if payload['status'] != 'PASS':
        raise SystemExit(2)


if __name__ == '__main__':
    main()
