#!/usr/bin/env python3
"""Linux original-CLS transverse failure witness with independent weak oracles.

This is a reduced linearized Fourier-x/native-CG-z experiment with three
physical velocity components. It is not a nonlinear 3D patch trajectory.
Native matrices come from the compiled original stress tangent. Independent
dimensional indexed tensors, Jacobi/GL integration and SciPy evolution build
the reference answer. A passing verification records a FAILED stability gate.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import uuid

sys.dont_write_bytecode = True
THREADS = {name: '1' for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                                  'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')}
os.environ.update(THREADS)
import numpy as np
import scipy
from scipy.linalg import expm, lu_factor, lu_solve
from scipy.optimize import brentq
import validate_cls_epawf_m1 as constitutive
from validate_cls_epawf_m1 import indexed
from validate_cls_epawf_m3 import basis

H, ZA, ZB, US, FINAL, SPATIAL_DT = 2., 4., 12., .4, .05, .00015625
MESHES = ((3, 2), (5, 2), (7, 2), (3, 4), (3, 8), (7, 8))
WAVES = (1., 2., 4.)
AMPLITUDES = (1e-5, 1e-4, 1e-3)
TIME_STEPS = (.0025, .00125, .000625)
NOTE_SHA = 'd35e437c2f2f8b25a46f7a3f285c6f1ad8b48d7b554aa1d161315f59eeb6de0f'


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            value.update(block)
    return value.hexdigest()


def experiment_plan():
    return dict(schema=1, milestone='M7', frozen_before_native_execution=True,
                experiment='reduced linearized Fourier-x / native Modified_A continuous Galerkin-z transverse original CLS',
                physical_components=3, perturbation='delta U=(0,v(x,z,t),0), partial_y=0, delta k=delta epsilon=0',
                support=dict(H=H, z=[ZA, ZB], boundary='homogeneous Dirichlet at both band endpoints, zero extension outside band',
                             qualification='compact H1 native envelope strictly above H; smooth target, not an evolving whole-domain support claim'),
                initial_envelope='exp(4-1/(s*(1-s))) for 0<s=(z-4)/8<1, zero otherwise; constrained physical L2 projection',
                meshes=[dict(degree=p, cells=n) for p, n in MESHES], xi=list(WAVES),
                wavelengths=[2*math.pi/xi for xi in WAVES], amplitudes=list(AMPLITUDES),
                native_spatial_time_step=SPATIAL_DT, final_time=FINAL,
                temporal=dict(mesh=dict(degree=7, cells=8), xi=4., amplitudes=[1e-4],
                              orders=[1, 2], steps=list(TIME_STEPS), reference='independent matrix exponential'),
                spatial_reference=dict(degree=9, cells=16, quadrature=64,
                    interpretation='same fine BDF2 step for actual native field comparisons; separate exact semidiscrete exponential control'),
                partition='all original tensor diffusion and base advection implicit; explicit derivative identically zero',
                acceptance=dict(mass_and_native_rows=2e-12, original_tangent=3e-9, complete_weak_matrix=3e-9,
                    projected_initial_finite_rule=1e-10, native_time_recurrence=5e-11,
                    covariance_sample_identity=2e-10, minimum_covariance='strictly positive at all sampled initial/final points/phases',
                    turbulence_linear_forcing=1e-10, spatial_fine_to_coarse_error_ratio=.5,
                    time_error_reduction_minimum={'IMEX1': 1.6, 'IMEX2': 3.},
                    growing_energy='gain>1 at all three wavelengths, retained as dt decreases',
                    amplitude_scaled_gain=2e-10, restricted_chi_eta_operator_identity=1e-12),
                excluded=['nonlinear 3D trajectory', 'uniform long-time stability', 'constitutive retuning',
                          'molecular outer diffusion', 'filter/SVV', 'positivity floors', 'new inner epsilon PDE'])


class Audit:
    def __init__(self):
        self.checks, self.errors = 0, {}

    def require(self, good, message):
        self.checks += 1
        if not good:
            raise RuntimeError(message)

    def near(self, actual, expected, name, tolerance=1e-10, scale=1.):
        actual, expected = np.asarray(actual), np.asarray(expected)
        self.require(actual.shape == expected.shape and np.all(np.isfinite(actual)) and
                     np.all(np.isfinite(expected)), 'Invalid '+name)
        value = float(np.max(np.abs(actual-expected), initial=0.))/max(scale, float(np.max(np.abs(expected), initial=0.)))
        self.errors[name] = max(self.errors.get(name, 0.), value)
        self.require(value <= tolerance, f'{name}: {value} > {tolerance}')
        return value


class Runner:
    def __init__(self, build, output, audit, receipt):
        self.build, self.output, self.audit, self.receipt = build, output, audit, receipt
        self.environment = dict(os.environ, **THREADS, PYTHONDONTWRITEBYTECODE='1')
        self.environment['LD_LIBRARY_PATH'] = str(build/'geometry-lib')+':'+self.environment.get('LD_LIBRARY_PATH', '')
        receipt['environment'] = dict(THREADS, TMPDIR=self.environment.get('TMPDIR'),
                                      LD_LIBRARY_PATH=self.environment['LD_LIBRARY_PATH'],
                                      cpu_affinity=sorted(os.sched_getaffinity(0)))

    def save(self):
        self.receipt.update(checks=self.audit.checks, errors=self.audit.errors)
        (self.output/'validation.json').write_text(json.dumps(self.receipt, indent=2, allow_nan=False)+'\n')

    def command(self, name, argv, cwd, timeout=600):
        log, start = self.output/(name+'.log'), time.monotonic()
        argv = list(map(str, argv))
        code, failure = 127, None
        try:
            with log.open('w') as stream:
                result = subprocess.run(argv, cwd=cwd, env=self.environment, stdout=stream,
                                        stderr=subprocess.STDOUT, timeout=timeout, check=False)
                code = result.returncode
        except subprocess.TimeoutExpired:
            code, failure = 124, 'Native transverse command timed out'
        except OSError as error:
            failure = str(error)
        if failure:
            with log.open('a') as stream:
                stream.write('\n'+failure+'\n')
        self.receipt['commands'].append(dict(name=name, argv=argv, cwd=str(cwd), ranks=1, threads=THREADS,
                                             exit_code=code, seconds=time.monotonic()-start,
                                             log=str(log), log_sha256=digest(log)))
        self.save()
        self.audit.require(code == 0, f'{name}: exit {code}; see {log}')
        return log.read_text().strip()


def model_reference():
    def c(q):
        return .3/(1+.35*q**1.5)*(-np.expm1(-.36*np.exp(.75*q)))
    Q = brentq(lambda q: c(q)*q*q-1, 1., 5., xtol=1e-14)
    def epsilon_balance(kappa):
        ratio = kappa*Q**1.5/2.5
        return Q*kappa*kappa/1.3+.83*max(ratio-1, 0)*ratio*ratio-.48
    kappa = brentq(epsilon_balance, .1, .8, xtol=1e-14)
    K, nu, rough = Q*US*US, 1.5e-5, .1
    x = rough*math.sqrt(K)/(70*nu)
    power = max(.5-.4*x**.7, 1-.79*x**(-.28))
    zv, slope = nu*10.7*(1-x**power)/math.sqrt(K), .2295*math.sqrt(K)
    split = min(H, max(0., zv))
    J0 = split/nu
    if H > split:
        J0 += math.log1p(slope*(H-split)/(nu+slope*(split-zv)))/slope
    return dict(Q=Q, kappa=kappa, C=c(Q), K=K, matching=US*US*J0)


def envelope(z):
    s = (np.asarray(z)-ZA)/(ZB-ZA)
    result = np.zeros_like(s)
    inside = (s > 0)&(s < 1)
    result[inside] = np.exp(4-1/(s[inside]*(1-s[inside])))
    return result


class Galerkin:
    def __init__(self, p, cells, xi, reference, order=64):
        self.p, self.cells, self.xi, self.reference = p, cells, xi, reference
        self.size = cells*p+1
        self.free = [i for i in range(self.size) if i not in (0, cells)]
        self.index = {value: i for i, value in enumerate(self.free)}
        self.maps = [[e, e+1]+[cells+1+e*(p-1)+j-2 for j in range(2, p+1)] for e in range(cells)]
        self.mass = np.zeros((len(self.free), len(self.free)))
        self.real, self.imaginary, self.pde_imaginary = [np.zeros_like(self.mass) for _ in range(3)]
        rhs = np.zeros(len(self.free))
        z, w = np.polynomial.legendre.leggauss(order)
        P, D = basis(p, z)
        width = (ZB-ZA)/cells
        for e, ids in enumerate(self.maps):
            points, weights = ZA+width*(e+(z+1)/2), w*width/2
            phi, derivative = P, 2*D/width
            nut = reference['kappa']*US*points
            a, b = (1-20*reference['C'])*nut, .36*reference['Q']*nut
            velocity = reference['matching']+US/reference['kappa']*np.log(points/H)
            m = phi.T@(weights[:, None]*phi)
            A = -xi*xi*phi.T@((weights*a)[:, None]*phi)-derivative.T@((weights*nut)[:, None]*derivative)
            B = -xi*(phi.T@((weights*b)[:, None]*derivative)+derivative.T@((weights*b)[:, None]*phi)+
                     phi.T@((weights*velocity)[:, None]*phi))
            C = xi*phi.T@((weights*(.36*reference['Q']*reference['kappa']*US-velocity))[:, None]*phi)
            load = phi.T@(weights*envelope(points))
            for i, gi in enumerate(ids):
                if gi not in self.index:
                    continue
                ii = self.index[gi]
                rhs[ii] += load[i]
                for j, gj in enumerate(ids):
                    if gj in self.index:
                        jj = self.index[gj]
                        self.mass[ii, jj] += m[i, j]
                        self.real[ii, jj] += A[i, j]
                        self.imaginary[ii, jj] += B[i, j]
                        self.pde_imaginary[ii, jj] += C[i, j]
        self.initial = np.linalg.solve(self.mass, rhs).astype(complex)
        self.load = self.real+1j*self.imaginary
        self.generator = np.linalg.solve(self.mass, self.load)

    def energy(self, state):
        return .5*float(np.vdot(state, self.mass@state).real)

    def exact(self):
        return expm(FINAL*self.generator)@self.initial

    def evolve(self, order, dt, amplitude=1.):
        count = round(FINAL/dt)
        y, old = amplitude*self.initial.copy(), None
        factors = {}
        values, histories = [y.copy()], []
        for step in range(count):
            lam = dt if order == 1 or step == 0 else 2*dt/3
            history = y.copy() if order == 1 or step == 0 else (4*y-old)/3
            if lam not in factors:
                factors[lam] = lu_factor(self.mass-lam*self.load)
            result = lu_solve(factors[lam], self.mass@history)
            histories.append(history.copy())
            old, y = y, result
            values.append(y.copy())
        return values, histories

    def values(self, state, points, derivative=False):
        points = np.asarray(points)
        cells = np.minimum(self.cells-1, np.maximum(0, ((points-ZA)/(ZB-ZA)*self.cells).astype(int)))
        width, result = (ZB-ZA)/self.cells, np.zeros(len(points), dtype=complex)
        for e in range(self.cells):
            indices = np.flatnonzero(cells == e)
            phi, grad = basis(self.p, 2*(points[indices]-ZA-e*width)/width-1)
            rows = 2*grad/width if derivative else phi
            coefficients = np.array([state[self.index[i]] if i in self.index else 0 for i in self.maps[e]])
            result[indices] = rows@coefficients
        return result


def paired(values):
    return np.r_[np.asarray(values).real, np.asarray(values).imag]


def unpaired(values):
    values = np.asarray(values)
    return values[:len(values)//2]+1j*values[len(values)//2:]


def field_error(a, ca, b, cb):
    # Union breakpoints handles both independent refinements without sampling
    # across a polynomial derivative discontinuity as one smooth panel.
    breaks = sorted(set(np.linspace(ZA, ZB, a.cells+1))|set(np.linspace(ZA, ZB, b.cells+1)))
    z, w = np.polynomial.legendre.leggauss(32)
    square = norm = 0.
    for lo, hi in zip(breaks[:-1], breaks[1:]):
        points, weights = lo+(z+1)*(hi-lo)/2, w*(hi-lo)/2
        av, bv = a.values(ca, points), b.values(cb, points)
        square += float(weights@np.abs(av-bv)**2)
        norm += float(weights@np.abs(bv)**2)
    return math.sqrt(max(0., square)/norm)


def covariance(galerkin, state):
    reference, minima, velocities = galerkin.reference, [], []
    nodes, _ = np.polynomial.legendre.leggauss(max(24, 2*galerkin.p+6))
    width = (ZB-ZA)/galerkin.cells
    points = np.concatenate([ZA+width*(e+(nodes+1)/2) for e in range(galerkin.cells)])
    value, derivative = galerkin.values(state, points), galerkin.values(state, points, True)
    for z, v, dv in zip(points, value, derivative):
        epsilon, G = US**3/(reference['kappa']*z), US/(reference['kappa']*z)
        for phase in (0., np.pi/4, np.pi/2, 3*np.pi/4):
            gradient = np.zeros((3, 3))
            gradient[0, 2] = G
            gradient[1, 0] = -galerkin.xi*(v.real*np.sin(phase)+v.imag*np.cos(phase))
            gradient[1, 2] = dv.real*np.cos(phase)-dv.imag*np.sin(phase)
            tensor = indexed(gradient, reference['K'], epsilon, 'physical')['stress'].real
            minima.append(float(np.linalg.eigvalsh(tensor)[0]))
            velocities.append(abs(v.real*np.cos(phase)-v.imag*np.sin(phase)))
    return [min(minima), max(velocities)]


def constitutive_checks(raw, reference, audit):
    records = []
    for sample in raw['constitutive_samples']:
        z, step, actual = sample['z'], sample['relative_step'], sample['value']
        epsilon, G = US**3/(reference['kappa']*z), US/(reference['kappa']*z)
        gradient = np.zeros((3, 3)); gradient[0, 2] = G
        base = indexed(gradient, reference['K'], epsilon, 'physical')
        block, forcing = np.zeros((2, 2)), []
        for column, axis in enumerate((0, 2)):
            direction = np.zeros((3, 3)); direction[1, axis] = 1
            tangent = indexed(gradient+1j*1e-30*direction, reference['K'], epsilon, 'strain')
            block[:, column] = -tangent['stress'].imag[1, [0, 2]]/1e-30
            forcing.append(float(tangent['production'].imag/1e-30))
        closed = float(base['nut'].real)*np.array([[1-20*reference['C'], -.36*reference['Q']],
                                                   [.36*reference['Q'], 1.]])
        audit.near(block, closed, 'indexed_closed_transverse_block', 2e-12)
        # The coarser physical differences demonstrate their O(step^2)
        # truncation; the matrices below use only the finest native step.
        bound = 3e-9+30*step*step
        error = audit.near(actual['block'], block, 'native_transverse_tangent', bound)
        audit.near(actual['nut'], float(base['nut'].real), 'neutral_eddy_viscosity', 2e-12)
        audit.near(actual['production'], epsilon, 'neutral_positive_production', 2e-12)
        audit.near(actual['covariance_minimum'], float(np.linalg.eigvalsh(base['stress'].real)[0]),
                   'neutral_covariance', 2e-12)
        audit.near(actual['linear_turbulence_forcing'], np.zeros(4), 'unforced_linear_turbulence', 1e-10)
        audit.near(forcing, np.zeros(2), 'independent_zero_production_derivative', 1e-12)
        audit.require(block[0, 0] < 0 and actual['covariance_minimum'] > 0,
                      'Original adverse diffusion or covariance positivity was lost')
        records.append(dict(z=z, relative_step=step, tangent_error=error,
                            adverse_ratio=float(block[0, 0]/base['nut'].real), covariance_minimum=actual['covariance_minimum']))
    audit.require(len(records) == 12, 'Missing physical constitutive refinement samples')
    return records


def verify(raw, audit):
    audit.require(raw['schema'] == 1 and raw['success'], 'Native reduced experiment failed')
    reference = model_reference()
    expected = [H, ZA, ZB, US, reference['K'], reference['Q'], reference['C'], reference['kappa'], reference['matching'], FINAL, SPATIAL_DT]
    audit.near(raw['parameters'], expected, 'independent_reference_inputs', 3e-12)
    audit.require(ZA > H and raw['time_partition'].startswith('all advection'), 'Wrong support/time partition')
    stress = constitutive_checks(raw, reference, audit)
    operators, native = {}, {record['id']: record for record in raw['operators']}
    audit.require(len(native) == 18, 'Missing p/h/wavelength operator matrix')
    initial_projection_errors = []
    for p, cells in MESHES:
        for xi in WAVES:
            name = f'p{p}-n{cells}-xi{int(xi)}'
            actual = native[name]
            space = Galerkin(p, cells, xi, reference)
            finite = Galerkin(p, cells, xi, reference, max(24, 2*p+6))
            operators[name] = space
            audit.require(actual['degree'] == p and actual['cells'] == cells and actual['size'] == space.size and
                          actual['maps'] == space.maps and actual['free_ids'] == space.free, 'Native CG map/boundary mismatch')
            audit.near(actual['mass'], space.mass, 'native_physical_mass', 2e-12)
            audit.near(actual['real_load'], space.real, 'complete_real_transverse_weak_matrix', 3e-9)
            audit.near(actual['imaginary_load'], space.imaginary, 'complete_imaginary_transverse_weak_matrix', 3e-9)
            audit.near(space.imaginary, space.pde_imaginary, 'mixed_stress_derivatives_cancel', 2e-12)
            audit.near(actual['initial_unit_coefficients'], paired(finite.initial), 'native_initial_finite_rule_projection', 1e-10)
            initial_projection_errors.append(dict(id=name, relative_error=field_error(space, unpaired(actual['initial_unit_coefficients']), space, space.initial)))
            audit.near(actual['linear_turbulence_forcing'], np.zeros(4), 'all_native_linear_turbulence_forcing', 1e-10)
            audit.require(actual['base_covariance_minimum'] > 0 and actual['initial_unit_energy_rate'] > 0,
                          'Native growing mode did not have positive covariance')
            for row in actual['native_rows']:
                element = actual['maps'].index(row['ids'])
                width = (ZB-ZA)/cells
                phi, derivative = basis(p, 2*(row['z']-ZA-element*width)/width-1)
                audit.near(row['phi'], phi[0], 'native_Modified_A_basis', 2e-12)
                audit.near(row['derivative'], 2*derivative[0]/width, 'native_Modified_A_derivative', 2e-12)

    trial_records, actual_final, exact_cache = [], {}, {}
    audit.require(len(raw['trials']) == 26, 'Unexpected frozen native trial count')
    for trial in raw['trials']:
        name, order, dt, amplitude = trial['operator'], trial['order'], trial['dt'], trial['amplitude']
        space = operators[name]
        # The independent finite-rule initial projection checks the native
        # declared projection; no exported native coefficient constructs it.
        finite = Galerkin(space.p, space.cells, space.xi, reference, max(24, 2*space.p+6))
        original_initial = space.initial
        space.initial = finite.initial
        values, histories = space.evolve(order, dt, amplitude)
        space.initial = original_initial
        audit.near(trial['initial_state'], paired(values[0]), 'native_initial_amplitude', 1e-10, amplitude)
        audit.require(len(trial['trajectory']) == len(values)-1, 'Native step count differs')
        for step, record in enumerate(trial['trajectory'], 1):
            audit.near(record['time'], step*dt, 'native_time', 3e-13)
            audit.near(record['stage_time'], step*dt, 'native_stage_epoch', 3e-13)
            audit.near(record['lambda'], dt if order == 1 or step == 1 else 2*dt/3, 'native_IMEX_lambda', 3e-13)
            audit.near(record['stage_history'], paired(histories[step-1]), 'native_complete_history_combination', 5e-11, amplitude)
            audit.near(record['state'], paired(values[step]), 'native_complete_IMEX_recurrence', 5e-11, amplitude)
            audit.near(record['energy'], space.energy(values[step]), 'native_mass_energy', 2e-10, amplitude*amplitude)
            audit.near(record['energy_rate'], float(np.vdot(values[step], space.load@values[step]).real),
                       'native_transverse_energy_identity', 3e-9, amplitude*amplitude)
            audit.require(record['stage_residual'] <= 1e-10*amplitude, 'Unconverged native implicit linear stage')
        final = unpaired(trial['trajectory'][-1]['state'])
        audit.near(trial['initial_covariance'], covariance(space, unpaired(trial['initial_state'])), 'finite_initial_covariance', 2e-10)
        audit.near(trial['final_covariance'], covariance(space, final), 'finite_final_covariance', 2e-10)
        audit.require(trial['initial_covariance'][0] > 0 and trial['final_covariance'][0] > 0,
                      'Finite perturbation lost positive covariance in this sample set')
        expected_history = [paired(values[-1-slot]) for slot in range(order)]+[np.zeros_like(paired(values[-1])) for _ in range(order)]
        audit.near(trial['native_history'], expected_history, 'every_native_history_slot', 5e-11, amplitude)
        audit.near(trial['native_history_times'], [FINAL-slot*dt for slot in range(order)]+[dt]*order, 'native_history_times', 3e-13)
        gain = space.energy(final)/space.energy(unpaired(trial['initial_state']))
        audit.near(trial['gain'], gain, 'native_energy_gain', 2e-10)
        audit.require(gain > 1, 'The adverse transverse energy growth was removed')
        key = (name, order, dt, amplitude)
        audit.require(key not in actual_final, 'Duplicate trial')
        actual_final[key] = final
        if name not in exact_cache:
            exact_cache[name] = expm(FINAL*space.generator)@finite.initial
        exact = amplitude*exact_cache[name]
        difference = final-exact
        relative_error = math.sqrt(max(0., float(np.vdot(difference, space.mass@difference).real))/
                                   float(np.vdot(amplitude*finite.initial, space.mass@(amplitude*finite.initial)).real))
        wavelength = 2*math.pi/space.xi
        center_velocity = reference['matching']+US/reference['kappa']*math.log((ZA+ZB)/(2*H))
        trial_records.append(dict(operator=name, order=order, dt=dt, amplitude=amplitude,
            gain=gain, relative_error_to_semidiscrete_exponential=relative_error,
            initial_covariance_minimum=trial['initial_covariance'][0], final_covariance_minimum=trial['final_covariance'][0],
            wavelength=wavelength, flow_time=wavelength/center_velocity, turnover_time=H/US,
            duration_in_flow_times=FINAL/(wavelength/center_velocity), duration_in_turnover_times=FINAL/(H/US)))

    spatial = []
    for xi in WAVES:
        high = Galerkin(9, 16, xi, reference)
        high_bdf = high.evolve(2, SPATIAL_DT, 1e-4)[0][-1]
        high_exact = high.exact()
        records = []
        for p, cells in MESHES:
            name = f'p{p}-n{cells}-xi{int(xi)}'
            space = operators[name]
            final = actual_final[(name, 2, SPATIAL_DT, 1e-4)]
            error = field_error(space, final, high, high_bdf)
            exact_error = field_error(space, space.exact(), high, high_exact)
            records.append(dict(degree=p, cells=cells, actual_native_relative_field_error=error,
                                independent_semidiscrete_spatial_error=exact_error,
                                actual_native_energy_gain=space.energy(final)/space.energy(unpaired(next(
                                    trial['initial_state'] for trial in raw['trials'] if trial['operator'] == name and trial['dt'] == SPATIAL_DT and trial['amplitude'] == 1e-4))),
                                independent_reference_gain=high.energy(high_bdf)/high.energy(1e-4*high.initial)))
        by_mesh = {(record['degree'], record['cells']): record for record in records}
        coarse = by_mesh[(3, 2)]['actual_native_relative_field_error']
        audit.require(by_mesh[(7, 2)]['actual_native_relative_field_error'] < .5*coarse,
                      'Independent p refinement did not reduce actual native transverse field error')
        audit.require(by_mesh[(3, 8)]['actual_native_relative_field_error'] < .5*coarse,
                      'Independent h refinement did not reduce actual native transverse field error')
        audit.require(by_mesh[(7, 2)]['independent_semidiscrete_spatial_error'] < .5*by_mesh[(3, 2)]['independent_semidiscrete_spatial_error'] and
                      by_mesh[(3, 8)]['independent_semidiscrete_spatial_error'] < .5*by_mesh[(3, 2)]['independent_semidiscrete_spatial_error'],
                      'PDE weak-operator spatial refinement failed separately from time error')
        audit.require(high.energy(high_bdf) > high.energy(1e-4*high.initial), 'Refined PDE weak operator lost adverse growth')
        spatial.append(dict(xi=xi, reference=dict(degree=9, cells=16, dt=SPATIAL_DT), records=records))

    temporal = []
    for order in (1, 2):
        records = [next(record for record in trial_records if record['operator'] == 'p7-n8-xi4' and
                        record['order'] == order and record['dt'] == dt and record['amplitude'] == 1e-4) for dt in TIME_STEPS]
        error = [record['relative_error_to_semidiscrete_exponential'] for record in records]
        ratio = [error[j]/error[j+1] for j in range(2)]
        audit.require(min(ratio) > (1.6 if order == 1 else 3.), 'Native time refinement did not recover its declared order')
        temporal.append(dict(order=order, records=records, actual_error_reduction=ratio))
    gains = [next(record['gain'] for record in trial_records if record['operator'] == 'p7-n8-xi4' and
                  record['order'] == 2 and record['dt'] == SPATIAL_DT and record['amplitude'] == amplitude) for amplitude in AMPLITUDES]
    audit.near(gains, np.full(3, gains[1]), 'linear_amplitude_energy_gain', 2e-10)
    controls = raw['restricted_controls']
    audit.require({(record['chi'], record['eta']) for record in controls} == {(c, e) for c in (0, .5, 1) for e in (0, 1)}, 'Missing restricted chi/eta controls')
    for record in controls:
        audit.near(record['matching_velocity'], reference['matching'], 'original_AWF_momentum_chi_independence', 2e-12)
        audit.near(record['maximum_operator_difference'], 0., 'restricted_outer_chi_eta_identity', 1e-12)
    return dict(physical_stability_gate='failed_as_specified', constitutive=stress, trials=trial_records,
                spatial=spatial, temporal=temporal, initial_projection_control=initial_projection_errors,
                amplitude=dict(values=list(AMPLITUDES), energy_gains=gains), restricted_controls=controls,
                control_scope=raw['control_scope'],
                covariance_scope='actual original nonlinear covariance sampled at native volume GL points and four Fourier phases, initial and final times; not a continuum minimum certificate',
                experiment_scope='reduced linearized outer-band Fourier/CG experiment; no nonlinear 3D trajectory or long-time stability claim')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', type=Path, required=True)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if platform.system() != 'Linux':
        raise RuntimeError('Scientific transverse computations belong on Linux')
    build = args.build.resolve()
    default_root = build.parent if str(build).startswith('/media/') else Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised')
    output = args.output.resolve() if args.output else default_root/('m7-transverse-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S')+'-'+uuid.uuid4().hex[:8])
    output.mkdir(parents=True, exist_ok=False)
    plan = experiment_plan()
    (output/'experiment-plan.json').write_text(json.dumps(plan, indent=2)+'\n')
    audit = Audit()
    receipt = dict(schema=1, milestone='M7', success=False, verification_status='failed', started_utc=datetime.now(timezone.utc).isoformat(),
                   commands=[], failures=[], experiment_plan=str(output/'experiment-plan.json'), experiment_plan_sha256=digest(output/'experiment-plan.json'),
                   physical_stability_gate='failed_as_specified', versions=dict(numpy=np.__version__, scipy=scipy.__version__))
    run = Runner(build, output, audit, receipt)
    source = Path(__file__).resolve().parents[4]
    pins = {}
    try:
        receipt['source_commit'] = run.command('git-start-revision', ['git', '-c', 'gc.auto=0', '-C', source, 'rev-parse', 'HEAD'], source).strip()
        status = run.command('git-start-status', ['git', '-c', 'gc.auto=0', '-C', source, 'status', '--porcelain=v1', '--untracked-files=all'], source)
        audit.require(not status, 'Source checkout is not clean before transverse qualification')
        receipt['source_clean_start'] = True
        audit.require(digest(source/'docs/literature/cls_neutral_abl_equilibrium_awf_note_revised.tex') == NOTE_SHA, 'Frozen revised note changed')
        for name in ('RANSConstitutiveProbe', 'RANSTransverseEvolutionProbe'):
            path = build/name
            audit.require(path.is_file(), 'Missing native target '+name)
            pins[name] = dict(path=str(path), sha256=digest(path))
        receipt['programs'] = pins
        kernel = output/'constitutive.json'
        run.command('compiled-original-constitutive', [build/'RANSConstitutiveProbe', '--output', kernel], output)
        data = json.loads(kernel.read_text())
        receipt['reused_constitutive_verification'] = dict(tangents=constitutive.tangents(data), reference=constitutive.profiles(data),
                                                         native_report_sha256=digest(kernel), checks=constitutive.checks)
        native = output/'native-transverse.json'
        run.command('native-reduced-transverse-evolution', [build/'RANSTransverseEvolutionProbe', '--output', native], output)
        receipt['native_report'] = dict(path=str(native), sha256=digest(native))
        receipt['sections'] = verify(json.loads(native.read_text()), audit)
        receipt['success'] = True
    except Exception as error:
        receipt['failures'].append(str(error))
    try:
        end = run.command('git-end-revision', ['git', '-c', 'gc.auto=0', '-C', source, 'rev-parse', 'HEAD'], source)
        status = run.command('git-end-status', ['git', '-c', 'gc.auto=0', '-C', source, 'status', '--porcelain=v1', '--untracked-files=all'], source)
        receipt.update(source_commit_end=end, source_clean_end=not status)
        audit.require(not status and end == receipt.get('source_commit'), 'Source changed during transverse qualification')
        for name, record in pins.items():
            record['sha256_end'] = digest(record['path'])
            audit.require(record['sha256_end'] == record['sha256'], 'Native program changed during qualification: '+name)
    except Exception as error:
        receipt['success'] = False
        receipt['failures'].append(str(error))
    receipt.update(finished_utc=datetime.now(timezone.utc).isoformat(), verification_status='passed' if receipt['success'] else 'failed')
    run.save()
    print(json.dumps(dict(success=receipt['success'], checks=audit.checks, physical_stability_gate=receipt['physical_stability_gate'],
                          failures=receipt['failures'], output=str(output))), flush=True)
    return 0 if receipt['success'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
