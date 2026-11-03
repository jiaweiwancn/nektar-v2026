#!/usr/bin/env python3
"""Independent native transverse qualification of explicit M8 CLS variants.

This repeats the M7 reduced Fourier-x/native-CG-z experiment with an explicitly
selected common cubic factor and prescribed gamma(z)=GammaScale*kappa*u_star*z.
The independent indexed constitutive law, full mixed weak operator, native
IMEX recurrences and matrix exponential all use the selected variant. Passing
means that this reference-local transverse mode decays under the declared
refinements. It does not qualify a global closure or a nonlinear trajectory.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
import uuid

sys.dont_write_bytecode = True
THREADS = {name: '1' for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                                  'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')}
os.environ.update(THREADS)
import numpy as np
import scipy
from scipy.linalg import expm
import validate_cls_epawf_m7_transverse as legacy
from validate_cls_epawf_m1 import indexed as original_indexed
from validate_cls_epawf_m3 import basis

H, ZA, ZB, US, FINAL, SPATIAL_DT = 2., 4., 12., .4, .05, .00015625
MESHES, WAVES, AMPLITUDES, TIME_STEPS = legacy.MESHES, legacy.WAVES, legacy.AMPLITUDES, legacy.TIME_STEPS
Audit, Runner = legacy.Audit, legacy.Runner
paired, unpaired, field_error = legacy.paired, legacy.unpaired, legacy.field_error
digest = legacy.digest
FROZEN = {
    'docs/literature/cls_neutral_abl_equilibrium_awf_note_revised.tex': legacy.NOTE_SHA,
    'docs/plans/cls-epawf-revised-contract-v1.json': '4e7bd75022409b8dc2ef31cc73425779b646a97535760b31eeea4ec49087de8a',
    'library/SolverUtils/Turbulence/RANSModels.cpp': '86a7154ba3cb27e4d3319aa42d156a228fd1b180f35f9b217159822ec9d66b05',
    'library/SolverUtils/Turbulence/AnalyticalWallFunction.cpp': '5baf29819fed6b81ad1b042bc09258a3e0779d15de13076b688bf32f80ddcf7e'}


def experiment_plan(gamma_scale):
    plan = legacy.experiment_plan()
    plan.update(milestone='M8', experiment='explicit variant, reduced linearized Fourier/native-CG reference-local transverse decay',
                variants=[dict(beta=.4, gamma_scale=0.), dict(beta=.45, gamma_scale=0.),
                          dict(beta=.4, gamma_scale=gamma_scale)],
                stress='R_beta_gamma=R_original+(beta-1)*R_cubic-gamma*(S-Sref); common beta and fixed physical gamma',
                prescribed_gamma='GammaScale*kappa_c*u_star*z; independent of evolving k,epsilon,gradient',
                prescribed_reference_strain='Sref_xz=Sref_zx=u_star/(kappa_c*z), all other entries zero',
                independent_operator='full indexed stress complex-step tangent; both mixed derivatives retained before independent cancellation',
                original_control='source-bound M7 receipt is retained as a separate original-law growth witness; no duplicate original campaign',
                qualification_scope='reference-local reduced transverse mode only; not a global ellipticity or coupled stability certificate')
    plan['acceptance'].pop('growing_energy')
    plan['acceptance']['decaying_energy'] = '0<gain<1 at every native trial and independently refined weak PDE reference; negative instantaneous energy rate'
    plan['acceptance']['positive_reference_transverse_diffusion'] = True
    plan['excluded'] = ['nonlinear 3D trajectory', 'uniform or global constitutive stability', 'certified off-reference envelope',
                        'whole-patch chi/eta equivalence', 'outer molecular diffusion', 'filter/SVV', 'positivity floors', 'new inner epsilon PDE']
    return plan


def indexed(gradient, k, epsilon, branch, beta, gamma, reference_strain):
    """Independent dimensional indexed CLS, with explicit selected increments."""
    g = np.asarray(gradient, dtype=complex)
    values = original_indexed(g, k, epsilon, branch)
    stress = values['stress']+(beta-1)*values['cubic']-gamma*(g+g.T-reference_strain)
    return dict(values, stress=stress, production=-np.sum(stress*g), cubic=beta*values['cubic'])


def reference_state(z, reference):
    epsilon, shear = US**3/(reference['kappa']*z), US/(reference['kappa']*z)
    gradient = np.zeros((3, 3)); gradient[0, 2] = shear
    return epsilon, gradient, gradient+gradient.T


def tangent(z, reference, beta, gamma_scale):
    epsilon, gradient, strain = reference_state(z, reference)
    gamma = gamma_scale*reference['kappa']*US*z
    block, forcing = np.zeros((2, 2)), []
    for column, axis in enumerate((0, 2)):
        direction = np.zeros((3, 3)); direction[1, axis] = 1
        value = indexed(gradient+1j*1e-30*direction, reference['K'], epsilon, 'strain', beta, gamma, strain)
        block[:, column] = -value['stress'].imag[1, [0, 2]]/1e-30
        forcing.append(float(value['production'].imag/1e-30))
    base = indexed(gradient, reference['K'], epsilon, 'physical', beta, gamma, strain)
    return block, forcing, base


class Galerkin(legacy.Galerkin):
    """Independent Jacobi CG assembly; native arrays never construct its rows."""
    def __init__(self, p, cells, xi, reference, beta, gamma_scale, order=64):
        self.p, self.cells, self.xi, self.reference = p, cells, xi, reference
        self.beta, self.gamma_scale = beta, gamma_scale
        self.size = cells*p+1
        self.free = [i for i in range(self.size) if i not in (0, cells)]
        self.index = {value: i for i, value in enumerate(self.free)}
        self.maps = [[e, e+1]+[cells+1+e*(p-1)+j-2 for j in range(2, p+1)] for e in range(cells)]
        self.mass = np.zeros((len(self.free), len(self.free)))
        self.real, self.imaginary, self.pde_imaginary = [np.zeros_like(self.mass) for _ in range(3)]
        rhs = np.zeros(len(self.free)); z, w = np.polynomial.legendre.leggauss(order)
        phi, local_derivative = basis(p, z); width = (ZB-ZA)/cells
        derivative = 2*local_derivative/width
        for e, ids in enumerate(self.maps):
            points, weights = ZA+width*(e+(z+1)/2), w*width/2
            blocks = np.array([tangent(float(point), reference, beta, gamma_scale)[0] for point in points])
            velocity = reference['matching']+US/reference['kappa']*np.log(points/H)
            m = phi.T@(weights[:, None]*phi)
            A = -xi*xi*phi.T@((weights*blocks[:, 0, 0])[:, None]*phi)-derivative.T@((weights*blocks[:, 1, 1])[:, None]*derivative)
            B = xi*(phi.T@((weights*blocks[:, 0, 1])[:, None]*derivative)-
                    derivative.T@((weights*blocks[:, 1, 0])[:, None]*phi)-phi.T@((weights*velocity)[:, None]*phi))
            # Separate integration-by-parts form: beta and prescribed gamma
            # leave the quadratic mixed coefficient and base advection intact.
            C = xi*phi.T@((weights*(.36*reference['Q']*reference['kappa']*US-velocity))[:, None]*phi)
            load = phi.T@(weights*legacy.envelope(points))
            for i, gi in enumerate(ids):
                if gi not in self.index:
                    continue
                ii = self.index[gi]; rhs[ii] += load[i]
                for j, gj in enumerate(ids):
                    if gj in self.index:
                        jj = self.index[gj]
                        self.mass[ii, jj] += m[i, j]; self.real[ii, jj] += A[i, j]
                        self.imaginary[ii, jj] += B[i, j]; self.pde_imaginary[ii, jj] += C[i, j]
        self.initial = np.linalg.solve(self.mass, rhs).astype(complex)
        self.load = self.real+1j*self.imaginary
        self.generator = np.linalg.solve(self.mass, self.load)


def covariance(space, state):
    reference, minima, velocities = space.reference, [], []
    nodes, _ = np.polynomial.legendre.leggauss(max(24, 2*space.p+6))
    width = (ZB-ZA)/space.cells
    points = np.concatenate([ZA+width*(e+(nodes+1)/2) for e in range(space.cells)])
    values, derivatives = space.values(state, points), space.values(state, points, True)
    for z, value, derivative in zip(points, values, derivatives):
        epsilon, gradient, strain = reference_state(z, reference)
        gamma = space.gamma_scale*reference['kappa']*US*z
        for phase in (0., np.pi/4, np.pi/2, 3*np.pi/4):
            G = gradient.copy()
            G[1, 0] = -space.xi*(value.real*np.sin(phase)+value.imag*np.cos(phase))
            G[1, 2] = derivative.real*np.cos(phase)-derivative.imag*np.sin(phase)
            stress = indexed(G, reference['K'], epsilon, 'physical', space.beta, gamma, strain)['stress'].real
            minima.append(float(np.linalg.eigvalsh(stress)[0]))
            velocities.append(abs(value.real*np.cos(phase)-value.imag*np.sin(phase)))
    return [min(minima), max(velocities)]


def verify(raw, audit, beta, gamma_scale):
    audit.require(raw['schema'] == 1 and raw['success'], 'Native selected-variant experiment failed')
    audit.near(raw['variant'], [beta, gamma_scale], 'explicit_selected_variant', 0.)
    reference = legacy.model_reference()
    expected = [H, ZA, ZB, US, reference['K'], reference['Q'], reference['C'], reference['kappa'], reference['matching'], FINAL, SPATIAL_DT]
    audit.near(raw['parameters'], expected, 'independent_unchanged_reference_inputs', 3e-12)
    audit.require(ZA > H and raw['time_partition'].startswith('all advection'), 'Wrong support or implicit partition')
    records = []
    for sample in raw['constitutive_samples']:
        z, step, actual = sample['z'], sample['relative_step'], sample['value']
        block, forcing, base = tangent(z, reference, beta, gamma_scale)
        closed = float(base['nut'].real)*np.array([[1-20*beta*reference['C']+gamma_scale, -.36*reference['Q']],
                                                  [.36*reference['Q'], 1+gamma_scale]])
        audit.near(block, closed, 'full_indexed_selected_transverse_block', 2e-12)
        error = audit.near(actual['block'], block, 'native_selected_transverse_tangent', 1e-10+30*step*step)
        if step == 1e-6:
            audit.require(error <= 1e-10, 'Finest variant tangent exceeds isolated-identity tolerance')
        epsilon, gradient, _ = reference_state(z, reference)
        original = original_indexed(gradient, reference['K'], epsilon, 'physical')
        audit.near(base['stress'].real, original['stress'].real, 'reference_stress_preserved', 2e-12)
        audit.near(actual['production'], epsilon, 'reference_production_preserved', 2e-12)
        audit.near(actual['nut'], original['nut'].real, 'original_diffusion_coefficient_retained', 2e-12)
        audit.near(actual['covariance_minimum'], np.linalg.eigvalsh(base['stress'].real)[0], 'reference_covariance', 2e-12)
        audit.near(actual['linear_turbulence_forcing'], np.zeros(4), 'compiled_zero_linear_turbulence_forcing', 1e-10)
        audit.near(forcing, np.zeros(2), 'indexed_zero_linear_turbulence_forcing', 1e-12)
        audit.require(block[0, 0] > 0 and block[1, 1] > 0 and actual['covariance_minimum'] > 0,
                      'Reference transverse diffusion or covariance is not positive')
        records.append(dict(z=z, relative_step=step, tangent_error=error,
                            diagonal_ratios=[float(block[0, 0]/base['nut'].real), float(block[1, 1]/base['nut'].real)]))
    audit.require(len(records) == 12, 'Missing constitutive difference refinement')

    operators, natives = {}, {r['id']: r for r in raw['operators']}
    audit.require(len(natives) == 18, 'Missing p/h/wavenumber native matrix')
    projections = []
    for p, cells in MESHES:
        for xi in WAVES:
            name = f'p{p}-n{cells}-xi{int(xi)}'; actual = natives[name]
            space = Galerkin(p, cells, xi, reference, beta, gamma_scale)
            finite = Galerkin(p, cells, xi, reference, beta, gamma_scale, max(24, 2*p+6))
            operators[name] = space
            audit.require(actual['degree'] == p and actual['cells'] == cells and actual['size'] == space.size and
                          actual['maps'] == space.maps and actual['free_ids'] == space.free, 'Native CG topology mismatch')
            audit.near(actual['mass'], space.mass, 'native_physical_mass', 2e-12)
            audit.near(actual['real_load'], space.real, 'complete_real_selected_weak_operator', 1e-10)
            audit.near(actual['imaginary_load'], space.imaginary, 'complete_imaginary_selected_weak_operator', 1e-10)
            audit.near(space.imaginary, space.pde_imaginary, 'independent_mixed_derivative_cancellation', 2e-12)
            audit.near(actual['initial_unit_coefficients'], paired(finite.initial), 'native_finite_rule_initial_projection', 1e-10)
            audit.near(actual['linear_turbulence_forcing'], np.zeros(4), 'all_compiled_zero_turbulence_forcing', 1e-10)
            audit.require(actual['base_covariance_minimum'] > 0 and actual['initial_unit_energy_rate'] < 0 and
                          np.linalg.eigvalsh(space.real)[-1] < 0, 'Selected weak mode is not dissipative')
            projections.append(dict(id=name, relative_error=field_error(space, unpaired(actual['initial_unit_coefficients']), space, space.initial)))
            for row in actual['native_rows']:
                element = actual['maps'].index(row['ids']); width = (ZB-ZA)/cells
                phi, derivative = basis(p, 2*(row['z']-ZA-element*width)/width-1)
                audit.near(row['phi'], phi[0], 'native_Modified_A_basis', 2e-12)
                audit.near(row['derivative'], 2*derivative[0]/width, 'native_Modified_A_derivative', 2e-12)

    trials, finals, exact_cache, finite_cache = [], {}, {}, {}
    audit.require(len(raw['trials']) == 26, 'Frozen native trial count changed')
    for trial in raw['trials']:
        name, order, dt, amplitude = trial['operator'], trial['order'], trial['dt'], trial['amplitude']
        space = operators[name]
        if name not in finite_cache:
            finite_cache[name] = Galerkin(space.p, space.cells, space.xi, reference, beta, gamma_scale, max(24, 2*space.p+6))
        finite = finite_cache[name]
        initial = space.initial; space.initial = finite.initial
        values, histories = space.evolve(order, dt, amplitude); space.initial = initial
        audit.near(trial['initial_state'], paired(values[0]), 'native_initial_amplitude', 1e-10, amplitude)
        audit.require(len(trial['trajectory']) == len(values)-1, 'Native step count differs')
        for step, row in enumerate(trial['trajectory'], 1):
            audit.near(row['time'], step*dt, 'native_time', 3e-13)
            audit.near(row['stage_time'], step*dt, 'native_stage_epoch', 3e-13)
            audit.near(row['lambda'], dt if order == 1 or step == 1 else 2*dt/3, 'native_IMEX_lambda', 3e-13)
            audit.near(row['stage_history'], paired(histories[step-1]), 'complete_native_stage_history', 5e-11, amplitude)
            audit.near(row['state'], paired(values[step]), 'complete_native_IMEX_recurrence', 5e-11, amplitude)
            audit.near(row['energy'], space.energy(values[step]), 'native_physical_energy', 1e-10, amplitude*amplitude)
            audit.near(row['energy_rate'], np.vdot(values[step], space.load@values[step]).real,
                       'native_complete_energy_identity', 1e-10, amplitude*amplitude)
            audit.require(row['stage_residual'] <= 1e-10*amplitude and row['energy_rate'] < 0,
                          'Unconverged or nondissipative selected native stage')
        final = unpaired(trial['trajectory'][-1]['state'])
        audit.near(trial['initial_covariance'], covariance(space, unpaired(trial['initial_state'])), 'finite_initial_covariance', 1e-10)
        audit.near(trial['final_covariance'], covariance(space, final), 'finite_final_covariance', 1e-10)
        audit.require(trial['initial_covariance'][0] > 0 and trial['final_covariance'][0] > 0, 'Sampled finite covariance is not positive')
        expected_history = [paired(values[-1-slot]) for slot in range(order)]+[np.zeros_like(paired(values[-1])) for _ in range(order)]
        audit.near(trial['native_history'], expected_history, 'all_native_history_slots', 5e-11, amplitude)
        audit.near(trial['native_history_times'], [FINAL-slot*dt for slot in range(order)]+[dt]*order, 'all_native_history_epochs', 3e-13)
        gain = space.energy(final)/space.energy(unpaired(trial['initial_state']))
        audit.near(trial['gain'], gain, 'actual_native_energy_gain', 1e-10)
        audit.require(0 < gain < 1, 'Actual selected native mode did not decay')
        key = (name, order, dt, amplitude); audit.require(key not in finals, 'Duplicate native trial'); finals[key] = final
        if name not in exact_cache:
            exact_cache[name] = expm(FINAL*space.generator)@finite.initial
        difference = final-amplitude*exact_cache[name]
        error = math.sqrt(max(0., float(np.vdot(difference, space.mass@difference).real))/
                          float(np.vdot(amplitude*finite.initial, space.mass@(amplitude*finite.initial)).real))
        wavelength = 2*math.pi/space.xi
        velocity = reference['matching']+US/reference['kappa']*math.log((ZA+ZB)/(2*H))
        trials.append(dict(operator=name, order=order, dt=dt, amplitude=amplitude, gain=gain,
            relative_error_to_semidiscrete_exponential=error, initial_covariance_minimum=trial['initial_covariance'][0],
            final_covariance_minimum=trial['final_covariance'][0], wavelength=wavelength, flow_time=wavelength/velocity,
            turnover_time=H/US, duration_in_flow_times=FINAL/(wavelength/velocity), duration_in_turnover_times=FINAL/(H/US)))

    spatial = []
    for xi in WAVES:
        high = Galerkin(9, 16, xi, reference, beta, gamma_scale)
        high_bdf, high_exact = high.evolve(2, SPATIAL_DT, 1e-4)[0][-1], high.exact()
        rows = []
        for p, cells in MESHES:
            name = f'p{p}-n{cells}-xi{int(xi)}'; space = operators[name]
            final = finals[(name, 2, SPATIAL_DT, 1e-4)]
            rows.append(dict(degree=p, cells=cells, actual_native_relative_field_error=field_error(space, final, high, high_bdf),
                independent_semidiscrete_spatial_error=field_error(space, space.exact(), high, high_exact),
                actual_native_energy_gain=space.energy(final)/space.energy(1e-4*finite_cache[name].initial)))
        lookup = {(r['degree'], r['cells']): r for r in rows}
        for key in ('actual_native_relative_field_error', 'independent_semidiscrete_spatial_error'):
            audit.require(lookup[(7, 2)][key] < .5*lookup[(3, 2)][key] and lookup[(3, 8)][key] < .5*lookup[(3, 2)][key],
                          'Independent p/h refinement did not reduce selected mode error: '+key)
        reference_gain = high.energy(high_bdf)/high.energy(1e-4*high.initial)
        audit.require(0 < reference_gain < 1, 'Refined independent weak PDE mode did not decay')
        spatial.append(dict(xi=xi, reference=dict(degree=9, cells=16, dt=SPATIAL_DT), reference_gain=reference_gain, records=rows))
    temporal = []
    for order in (1, 2):
        rows = [next(r for r in trials if r['operator'] == 'p7-n8-xi4' and r['order'] == order and
                     r['dt'] == dt and r['amplitude'] == 1e-4) for dt in TIME_STEPS]
        errors = [r['relative_error_to_semidiscrete_exponential'] for r in rows]
        ratios = [errors[i]/errors[i+1] for i in range(2)]
        audit.require(min(ratios) > (1.6 if order == 1 else 3.), 'Selected native temporal refinement did not recover order')
        temporal.append(dict(order=order, records=rows, actual_error_reduction=ratios))
    gains = [next(r['gain'] for r in trials if r['operator'] == 'p7-n8-xi4' and r['order'] == 2 and
                  r['dt'] == SPATIAL_DT and r['amplitude'] == amplitude) for amplitude in AMPLITUDES]
    audit.near(gains, np.full(3, gains[1]), 'linear_amplitude_energy_gain', 1e-10)
    controls = raw['restricted_controls']
    audit.require({(r['chi'], r['eta']) for r in controls} == {(c, e) for c in (0, .5, 1) for e in (0, 1)}, 'Missing restricted controls')
    for row in controls:
        audit.near(row['matching_velocity'], reference['matching'], 'original_AWF_matching_retained', 2e-12)
        audit.near(row['maximum_operator_difference'], 0., 'restricted_outer_chi_eta_identity', 1e-12)
    return dict(beta=beta, gamma_scale=gamma_scale, qualification='reference_local_reduced_transverse_decay',
                global_or_coupled_stability='not_qualified', constitutive=records, trials=trials, spatial=spatial, temporal=temporal,
                initial_projection_control=projections, amplitude=dict(values=list(AMPLITUDES), energy_gains=gains),
                restricted_controls=controls, control_scope=raw['control_scope'],
                covariance_scope='selected nonlinear stress sampled at native volume GL points and four phases at initial/final epochs; no continuum minimum claim')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', type=Path, required=True)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--gamma-scale', type=float, default=.25, help='Declared prescribed-gamma sensitivity; no global envelope is inferred')
    parser.add_argument('--original-control', type=Path, help='Existing qualified original M7 transverse validation.json')
    args = parser.parse_args()
    if platform.system() != 'Linux':
        raise RuntimeError('Scientific transverse validation belongs on Linux')
    if not math.isfinite(args.gamma_scale) or args.gamma_scale <= 0:
        raise ValueError('Gamma sensitivity scale must be finite and strictly positive')
    source, build = Path(__file__).resolve().parents[4], args.build.resolve()
    base = build.parent if str(build).startswith('/media/') else Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised')
    output = args.output.resolve() if args.output else base/('m8-transverse-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S')+'-'+uuid.uuid4().hex[:8])
    output.mkdir(parents=True, exist_ok=False)
    (output/'experiment-plan.json').write_text(json.dumps(experiment_plan(args.gamma_scale), indent=2)+'\n')
    audit = Audit()
    receipt = dict(schema=1, milestone='M8', success=False, verification_status='failed',
        started_utc=datetime.now(timezone.utc).isoformat(), commands=[], failures=[], variants=[],
        experiment_plan=str(output/'experiment-plan.json'), experiment_plan_sha256=digest(output/'experiment-plan.json'),
        physical_stability_gate='reference_local_transverse_not_yet_verified', global_or_coupled_stability='not_qualified',
        versions=dict(numpy=np.__version__, scipy=scipy.__version__), script_sha256=digest(__file__))
    run = Runner(build, output, audit, receipt); pin = None
    try:
        receipt['source_commit'] = run.command('git-start-revision', ['git', '-c', 'gc.auto=0', '-C', source, 'rev-parse', 'HEAD'], source)
        status = run.command('git-start-status', ['git', '-c', 'gc.auto=0', '-C', source, 'status', '--porcelain=v1', '--untracked-files=all'], source)
        audit.require(not status, 'Source checkout is dirty before M8 transverse qualification'); receipt['source_clean_start'] = True
        for name, expected in FROZEN.items():
            audit.require(digest(source/name) == expected, 'Frozen original numerical input changed: '+name)
        receipt['frozen_inputs'] = FROZEN
        original_path = (args.original_control or source/'docs/plans/cls-epawf-revised-evidence/m7/transverse-validation.json').resolve()
        original = json.loads(original_path.read_text())
        audit.require(original['success'] and original['verification_status'] == 'passed' and
                      original['physical_stability_gate'] == 'failed_as_specified', 'Original-law control is not qualified as a growth witness')
        receipt['original_control'] = dict(path=str(original_path), sha256=digest(original_path), source_commit=original['source_commit'],
            checks=original['checks'], physical_stability_gate=original['physical_stability_gate'], scope='reused original M7 witness, not a newly rerun legacy suite')
        program = build/'RANSTransverseEvolutionProbe'; audit.require(program.is_file(), 'Missing native transverse executable')
        pin = dict(path=str(program), sha256=digest(program)); receipt['program'] = pin
        for index, (beta, gamma_scale) in enumerate(((.4, 0.), (.45, 0.), (.4, args.gamma_scale)), 1):
            native = output/f'variant-{index}-native.json'
            run.command(f'variant-{index}-native', [program, '--output', native, '--beta', repr(beta), '--gamma-scale', repr(gamma_scale)], output)
            section = verify(json.loads(native.read_text()), audit, beta, gamma_scale)
            receipt['variants'].append(dict(native_report=str(native), native_report_sha256=digest(native), sections=section)); run.save()
        receipt.update(success=True, physical_stability_gate='reference_local_transverse_decay_verified')
    except Exception as error:
        receipt['failures'].append(str(error))
    try:
        end = run.command('git-end-revision', ['git', '-c', 'gc.auto=0', '-C', source, 'rev-parse', 'HEAD'], source)
        status = run.command('git-end-status', ['git', '-c', 'gc.auto=0', '-C', source, 'status', '--porcelain=v1', '--untracked-files=all'], source)
        receipt.update(source_commit_end=end, source_clean_end=not status)
        audit.require(not status and end == receipt.get('source_commit'), 'Source changed during M8 transverse qualification')
        if pin:
            pin['sha256_end'] = digest(pin['path']); audit.require(pin['sha256'] == pin['sha256_end'], 'Native transverse executable changed')
    except Exception as error:
        receipt['success'] = False; receipt['failures'].append(str(error))
    receipt.update(finished_utc=datetime.now(timezone.utc).isoformat(), verification_status='passed' if receipt['success'] else 'failed')
    run.save()
    print(json.dumps(dict(success=receipt['success'], checks=audit.checks, failures=receipt['failures'], output=str(output))), flush=True)
    return 0 if receipt['success'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
