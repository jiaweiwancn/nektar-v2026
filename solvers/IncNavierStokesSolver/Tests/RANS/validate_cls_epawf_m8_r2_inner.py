#!/usr/bin/env python3
"""Focused compiled R2 kernel and ACTIVE actual EvolutionPatch verification.

The eight fixed native fixtures are operator evaluations, never trajectories or
qualified physical controls.  Both inactive identity and active full weak
loads are reconstructed independently. The legal rK=1 cubic-depletion case
is uncapped. A separately retained, pressure-active cap fixture supplies an
additional all-row weak/mass/trace check without a new native call. Neither
fixture is a qualified physical control. Numerical execution is Linux-only.

With --reuse-run alone this performs read-only reanalysis, leaving missing
fixtures unperformed. --continue-missing separately authorizes execution of
only those fixed missing fixtures with the identical retained probe binaries.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import traceback
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
THREADS = {name: '1' for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
    'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')}
os.environ.update(THREADS)
import numpy as np

from cls_epawf_m3_cases import set_info, set_parameter, write
from cls_epawf_m4_cases import tensor_mesh
from cls_epawf_m8_r2_oracle import indexed_inner, inner_mean
from cls_epawf_m8_r1_kernel import inner_value as r1_value
from cls_epawf_m7_energy_oracle import Audit, _fields, _horizontal, _sample, _velocity
from cls_epawf_m8_native_oracle import ray_source, selected_stress
from cls_epawf_m8_r2_preflight import weak_loads, selected_direct, selected_interface
from cls_epawf_m6_operator_oracle import _segment_projection


def stamp():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def reserve(path):
    for index in range(1000):
        candidate = path if index == 0 else path.with_name(path.name+'-'+str(index))
        try:
            candidate.mkdir(parents=True, exist_ok=False)
            return candidate
        except FileExistsError:
            pass
    raise RuntimeError('Could not reserve a fresh result directory')


class Verification:
    def __init__(self):
        self.checks = 0
        self.errors = {}

    def require(self, condition, name):
        self.checks += 1
        if not condition:
            raise ValueError(name)

    def near(self, actual, expected, name, tolerance=2e-10):
        actual, expected = np.asarray(actual), np.asarray(expected)
        self.require(actual.shape == expected.shape, name+' shape')
        # Component scaling prevents a large source or count from hiding a
        # pressure, covariance or small-field discrepancy.
        error = float(np.max(np.abs(actual-expected)/np.maximum(1., np.abs(expected)), initial=0.))
        self.errors[name] = max(error, self.errors.get(name, 0.))
        self.require(np.all(np.isfinite(actual)) and error <= tolerance,
                     name+': '+str(error)+' > '+str(tolerance))
        return error


def kernel(report, audit):
    cases = []
    for case in report['cases']:
        name, selector = case['name'], case['selector']
        G, H = np.asarray(case['gradient']), np.asarray(case['direction'])
        k, nua, dk, dn = case['k'], case['nuA'], case['dk'], case['dnuA']
        expected = indexed_inner(G, k, nua, selector)
        native = case['value']
        for key in ('stress', 'turbulent', 'strain', 'production', 's', 'a', 'c', 't', 'ta', 'tc', 'radial', 'tangential'):
            audit.near(native[key], np.real(expected[key]), 'compiled_indexed_'+key, 2e-12)
        mean = inner_mean(G, k, nua, case['nu'], case['pi'], selector)
        audit.near(native['mean_stress'], mean.real, 'compiled_mean_stress', 2e-12)
        physical = -(case['pi']-2*k/3)*np.eye(3)+case['nu']*(G+G.T)-expected['stress']
        for n in (np.array([1., 0., 0.]), np.array([0., 0., 1.]), np.array([1., 2., -1.])/np.sqrt(6)):
            audit.near(np.asarray(native['mean_stress'])@n, physical.real@n, 'actual_full_traction', 2e-12)
        if selector == 'SaturatedR2':
            historical = r1_value(G, k, nua)
            audit.near(native['stress'], historical['R'], 'compiled_frozen_R1_stress', 2e-12)
            audit.near(native['turbulent'], historical['F'], 'compiled_frozen_R1_F', 2e-12)
        directional_only = selector == 'SaturatedR2' and name in ('joint_zero_k_strain', 'joint_zero_k_viscosity')
        audit.require(case['action']['directional_only'] == directional_only,
                      'joint-origin derivative classification '+name)
        if not case['one_sided'] and name != 'saturation_join':
            step = 1e-30
            complex_value = indexed_inner(G+1j*step*H, k+1j*step*dk, nua+1j*step*dn, selector)
            for key in ('stress', 'turbulent', 'strain', 'production'):
                audit.near(case['action'][key], np.imag(complex_value[key])/step,
                           'compiled_full_action_'+key, 3e-11)
        ladders = []
        for sample in case['differences']:
            h = sample['step']
            for sign, side in ((1, 'plus'), (-1, 'minus')):
                if side not in sample:
                    continue
                value = indexed_inner(G+sign*h*H, k+sign*h*dk, nua+sign*h*dn, selector)
                for key in ('stress', 'turbulent', 'production'):
                    audit.near(sample[side][key], np.real(value[key]), 'compiled_perturbed_'+key, 3e-12)
            errors = []
            for key in ('stress', 'turbulent', 'production'):
                other = sample.get('minus', native)[key]
                divisor = h if case['one_sided'] else 2*h
                difference = (np.asarray(sample['plus'][key])-np.asarray(other))/divisor
                action = np.asarray(case['action'][key])
                errors.append(float(np.max(np.abs(difference-action)/np.maximum(1., np.abs(action)), initial=0.)))
            ladders.append(dict(step=h, maximum_action_error=max(errors)))
        audit.require(min(row['maximum_action_error'] for row in ladders) < 3e-6,
                      'No converging admissible action ladder: '+selector+'/'+name)
        if selector == 'SaturatedR2' and k == 0:
            audit.near(native['stress'], np.zeros((3, 3)), 'exact_zero_k_R', 0.)
        cases.append(dict(name=name, selector=selector, active=native['active'],
                          directional_only=directional_only, difference_ladder=ladders))
    audit.require({entry['name'] for entry in report['rejections']} == {
        'negative_k', 'negative_nuA', 'unknown_selector', 'outward_zero_k',
        'outward_zero_nuA', 'nonfinite_gradient'}, 'Missing kernel domain rejections')
    return dict(cases=cases, native_checks=report['checks'],
                scope='Compiled complete inner law and admissible directional actions; no PDE/discrete-positivity inference.')


def session(template, destination, selector, fixture, chi, rK):
    root = ET.parse(template).getroot()
    tensor_mesh(root, [0., 1., 2., 4.] if fixture == 4 else [0., 2., 4.],
                3, 2, 1., 1., points=8)
    for key, value in {'RANSWallTransition': 'Fixed', 'RANSWallMap': 'Flat',
            'RANSWallVariant': 'Regularized', 'RANSWallDestruction': 'Capped',
            'RANSWallInnerStress': selector}.items():
        set_info(root, key, value)
    for key, value in {'RANSWallHeight': 2, 'RANSWallSampleHeight': rK,
            'RANSWallTailStart': 1, 'RANSWallChi': chi, 'RANSWallEta': 0,
            'RANSWallCubicBeta': .4, 'RANSWallRegularizationGamma': 0,
            'RANSWallRegularizationScale': 1, 'RANSWallPositivityTheta': .2,
            'RANSWallPositivityDelta': .1, 'RANSWallTangentialOrder': 8,
            'RANSWallNormalOrder': 16, 'RANSWallProbeInnerStressFixture': fixture}.items():
        set_parameter(root, key, value)
    return write(root, destination)


def constrained(raw, values, input_state, constraints):
    result = np.asarray(values).copy()
    for row, value in constraints:
        result[int(row)] = input_state[int(row)]-value
    return result


def lifted_traction(raw, fields, state, audit):
    """Independent order-12 traction on the actual lifted operator state.

    Audit() evaluates its supplied state. The probe supplies raw projected
    state to Audit but lifted state to Interface/EnergyBudget/physical weak
    loads. Preserve and verify both, instead of comparing distinct states.
    """
    d, off, H, nu = raw['dimension'], raw['offsets'], raw['parameters'][1], raw['parameters'][0]
    integral = np.zeros(6)
    maximum_trace_difference = 0.

    def at(point):
        ray = ray_source(raw, fields, state, point)
        profile, metric = ray['profile'], ray['metric']
        _, Gi = _velocity(fields, state, off, d, [point], side=-1)
        _, Go = _velocity(fields, state, off, d, [point], side=1)
        k = _sample(fields, state, off, d, [point], side=-1)[3][0]
        epsilon = _sample(fields, state, off, d+1, [point], side=1)[3][0]
        pi = _sample(fields, state, off, d+2, [point], side=-1)[3][0]
        po = _sample(fields, state, off, d+3, [point], side=1)[3][0]
        nt = float(profile.a)*max(metric['height']-float(profile.zv), 0.)
        selector = raw['inner_stress_model']
        Ti = inner_mean(Gi[0], k, nt, nu, pi, selector).real
        R, _ = selected_stress(raw, Go[0], k, epsilon, point)
        To = -R+(2*k/3-po)*np.eye(3)
        dependent = ray['lam']*float(np.asarray(profile.fields([metric['height']]), float)[2, 0])
        return Ti@metric['matching_area'], To@metric['matching_area'], epsilon, dependent

    for point, weight in _horizontal(fields, d, H, 12):
        ti, to, epsilon, dependent = at(point)
        integral[:3] += weight*ti;integral[3:] += weight*to
        maximum_trace_difference = max(maximum_trace_difference, abs(epsilon-dependent))
    audit.near(np.asarray(raw['interface']['integrals'])[6:12], integral,
               'actual_lifted_interface_all_traction_integrals', 3e-10)
    audit.near(raw['interface']['maximum_epsilon_projection_difference'], maximum_trace_difference,
               'actual_lifted_projected_trace_difference', 3e-10)
    for sample in raw['interface']['samples']:
        ti, to, epsilon, dependent = at(np.asarray(sample['chart']))
        audit.near(sample['traction_inner_area'], ti, 'actual_lifted_inner_traction_point', 3e-10)
        audit.near(sample['traction_outer_area'], to, 'actual_lifted_outer_traction_point', 3e-10)
        audit.near(sample['epsilon'], epsilon, 'actual_lifted_outer_epsilon_point', 3e-10)
        audit.near(sample['dependent_epsilon'], dependent, 'actual_lambda_EpreH_point', 3e-10)
    return dict(state_key='lifted_state', quadrature=12, traction_integrals=integral.tolist(),
                maximum_trace_projection_difference=maximum_trace_difference,
                point_checks=len(raw['interface']['samples']))


def native(raw, audit):
    audit.require(raw['inner_stress_fixture'] in (1, 2, 3, 4) and raw['unforced_alpha0_fixture']
                  and raw['alpha'] == 0 and not raw['upstream_transition'], 'Not the declared verification-only native fixture')
    fields, state = _fields(raw), np.asarray(raw['lifted_state'])
    expected, metadata = weak_loads(raw, fields, state)
    audit.near(raw['raw_residual'], expected, 'actual_complete_active_weak_operator', 2e-10)
    audit.near(np.asarray(raw['region_residuals'][0])+np.asarray(raw['region_residuals'][1]),
               expected, 'actual_regional_weak_operator_sum', 2e-10)
    old = Audit()
    direct = selected_direct(raw, fields, state, old)
    # The existing probe deliberately preserves the projected state and its
    # supplied-state Audit separately from the lifted physical operator. No
    # production/audit result is changed by this postprocessor correction.
    audit_key = 'lifted_state' if raw['upstream_transition'] else 'state'
    audit_state = np.asarray(raw[audit_key])
    supplied_direct = selected_direct(raw, fields, audit_state, old)
    observed = raw['constitutive_audit']
    audit.near(observed['wall_sources'], supplied_direct['wall_sources'], 'supplied_state_selected_wall_P_and_sources', 3e-10)
    audit.near(observed['robust_source'], supplied_direct['robust_source'], 'supplied_state_capped_budget_and_trace', 3e-10)
    audit_interface = selected_interface(raw, fields, audit_state)
    audit.near(observed['interface'], audit_interface, 'supplied_state_full_normal_tangential_traction', 3e-10)
    actual_traction = lifted_traction(raw, fields, state, audit)
    audit.near(observed['source_volume_inner'][3], supplied_direct['inner_covariance_minimum'], 'supplied_state_selected_covariance_including_zero', 3e-10)
    audit.near(np.asarray(raw['energy_budget'])[[2, 3, 4]], [direct['stress_work_minus_k_source'],
        direct['resolved_advection'], direct['roughness_force_work']], 'actual_same_stress_energy_cancellation', 3e-10)
    if raw['inner_stress_model'] == 'SaturatedR2':
        audit.near(observed['inner_stress'], supplied_direct['inner_stress'], 'supplied_state_saturation_volume_and_range', 3e-10)
        if raw['inner_stress_fixture'] == 2:
            audit.require(observed['inner_stress'][0] == 0, 'Inactive fixture unexpectedly saturates')
        else:
            audit.require(observed['inner_stress'][0] > 0, 'ACTIVE native fixture did not exercise saturation')
    if raw['inner_stress_fixture'] == 3:
        audit.require(raw['parameters'][3] == 1 and observed['robust_source'][4] > 0
                      and observed['robust_source'][5] == 0 and observed['robust_source'][6] == 1,
                      'Legal depleted chi=0 fixture must limit locally but remain uncapped')
    if raw['inner_stress_fixture'] == 4:
        zero = observed['inner_zero_covariance']
        audit.require(zero[0] > 0, 'The actual covariance diagnostic skipped exact inner k=0')
        if raw['inner_stress_model'] == 'SaturatedR2':
            audit.require(zero[1] == 0 and zero[2] == 0,
                          'Selected R2 covariance is not exactly zero at native k=0')
        else:
            audit.require(zero[2] < 0, 'Linear zero-k negative covariance control was not observed')
        audit.require(not raw['r2_coupled_action']['samples'] and 'skipped' in raw['r2_coupled_action'],
                      'Zero-k fixture must not take an outward central dk perturbation')
    actions = []
    for sample in raw['r2_coupled_action']['samples']:
        calculated = {}
        for side in ('plus', 'minus'):
            value, _ = weak_loads(raw, fields, np.asarray(sample['lifted_'+side]))
            audit.near(sample['raw_'+side], value, 'actual_mixed_perturbed_weak_operator', 3e-10)
            calculated[side] = constrained(raw, value, sample['input_'+side], sample['constraints_'+side])
        expected_action = (calculated['plus']-calculated['minus'])/(2*sample['step'])
        error = audit.near(sample['native_action'], expected_action,
                           'actual_mixed_local_K_pressure_trace_action', 8e-6)
        actions.append(dict(step=sample['step'], independent_action_error=error))
    return dict(inner_stress_model=raw['inner_stress_model'], fixture=raw['inner_stress_fixture'],
        chi=raw['chi'], weak=metadata, direct=direct, coupled_action=actions,
        supplied_audit_state_key=audit_key, supplied_audit_direct=supplied_direct,
        supplied_audit_interface=audit_interface.tolist(), lifted_traction=actual_traction,
        physical_signature=raw['physical_signature'],
        old_oracle_checks=old.checks, old_oracle_errors=old.errors,
        qualification='Operator verification only; these active states are outside the fixed outer box. No stage/trajectory or qualified physical solution is claimed.')


def retained_pressure_cap(directory, producing_source, audit):
    """Independently check the already executed, same-epoch capped operator.

    This consumes exact existing native files only. The capped case has active
    state-dependent chi; it is explicitly distinct from the uncapped chi=0
    CurvedWallProbe fixture and does not establish a cap-active chi=0 test.
    """
    directory = directory.resolve()
    validation_path, receipt_path = directory/'validation.json', directory.parent/'receipt.json'
    validation, receipt = (json.loads(path.read_text()) for path in (validation_path, receipt_path))
    audit.require(validation['success'] and validation['source_commit'] == producing_source,
                  'Retained pressure-cap verification is not passed at the producing binary epoch')
    audit.require(receipt['suite'] == 'pressure' and receipt['status'] == 'passed'
                  and receipt['runtime_unchanged'] and receipt['source_clean_start'] and receipt['source_clean_end']
                  and receipt['source_commit'] == receipt['source_commit_end'] == producing_source,
                  'Retained pressure-cap wrapper lacks a clean unchanged producing epoch')
    pin_paths = [directory.parent/'runtime-before.json', directory.parent/'runtime-after.json']
    before, after = (json.loads(path.read_text()) for path in pin_paths)
    for pin in (before, after):
        audit.require(pin['success'] and pin['clean'] and pin['source_commit'] == producing_source,
                      'Retained pressure-cap runtime fingerprint failed')
    for key in ('source_inputs', 'runtime_artifacts', 'runtime_bindings', 'runtime_environment',
                'python_packages', 'tool_images'):
        audit.require(before[key] == after[key], 'Retained pressure-cap runtime changed: '+key)
    program = Path(validation['program']['path'])
    audit.require(sha(program) == validation['program']['sha256'], 'Retained pressure-cap program bytes changed')
    raw_path = directory/'pressure-active/probe.json'
    command = next(row for row in validation['commands']
                   if Path(row['command'][-1]).resolve() == raw_path)
    audit.require(command['exit_code'] == 0 and command['report_sha256'] == sha(raw_path),
                  'Retained pressure-cap native report bytes changed')
    raw = json.loads(raw_path.read_text())
    audit.require(raw['fixture'] == 'pressure_sensor_and_capped_trace_active'
                  and raw['inner_stress_model'] == 'SaturatedR2' and raw['operator_verification_only']
                  and not raw['physical_qualification'] and raw['parameters'][1:5] == [2, .1, 1, 1]
                  and raw['robust_parameters'] == [.4, 0, .2, .1, 0, 1, 1],
                  'Retained pressure-cap fixture differs from the fixed legal contract')
    fields, state, history = _fields(raw), np.asarray(raw['state']), np.asarray(raw['history'])
    expected, metadata = weak_loads(raw, fields, state)
    mass = np.zeros(len(state))
    for block in range(raw['dimension']+2):
        ids = raw['offsets'][block]+np.arange(fields[block].size)
        mass[ids] = raw['alpha']*fields[block].mass()@(state[ids]-history[ids])
    expected = constrained(raw, expected+mass, state, raw['constraints'])
    residual_error = audit.near(raw['base_residual'], expected,
        'retained_cap_actual_complete_weak_mass_constraints', 3e-10)

    rays = {}

    def ray(chart, baseline=False):
        key = (float(chart[0]), bool(baseline))
        if key not in rays:
            view = dict(raw, upstream_transition=False, chi=0.) if baseline else raw
            rays[key] = ray_source(view, fields, state, chart)
        return rays[key]

    cap_count, interior_chi_count, ray_records = 0, 0, []
    H, d, offsets = raw['parameters'][1], raw['dimension'], raw['offsets']
    for sample in raw['samples']:
        value = ray(np.asarray(sample['chart']))
        profile, measure = value['profile'], value['measure']
        pre = float(np.asarray(profile.fields([H]), float)[2, 0])
        audit.near([sample['K'], sample['b'], sample['m'], sample['lambda']],
            [float(profile.k), value['b'], value['m'], value['lam']], 'retained_cap_full_ray_normalization', 3e-10)
        audit.near([sample['pre_trace'], sample['actual_trace']],
            [pre, value['lam']*pre], 'retained_cap_actual_lambda_Epre_trace', 3e-10)
        pre_e, pre_d = float(measure@value['pre_e']), float(measure@value['pre_d'])
        e, destruction = float(measure@value['e']), float(measure@value['destruction'])
        audit.near([pre_e, e, float(measure@value['redistribution'])],
            [pre_d, destruction, 0.], 'retained_cap_pre_and_actual_budget_cancellation', 3e-10)
        cap_count += int(value['b'] < raw['robust_parameters'][3] and 0 < value['lam'] < 1)
        interior_chi_count += int(0 < float(profile.chi) < 1)
        ray_records.append(dict(chart=sample['chart'], chi=float(profile.chi), b=value['b'],
            lambda_value=value['lam'], pre_E_budget=pre_e, pre_D_budget=pre_d,
            actual_E_budget=e, actual_D_budget=destruction))
    audit.require(cap_count > 0 and interior_chi_count > 0,
                  'Retained pressure fixture must exercise both attenuation and interior chi')

    # Native projection is baseline(chi=0) plus the current-minus-baseline
    # correction; downstream baseline-face correction modes are constrained.
    tables, baseline_ids = [{}, {}], set()
    for element in fields[d+1].elements:
        if abs(element['lower'][2]-H) > 1e-11:
            continue
        nodes, transform = _segment_projection(element['degree'])
        xs = element['lower'][0]+(nodes+1)*(element['upper'][0]-element['lower'][0])/2
        data = []
        for x in xs:
            chart = np.array([x, 0., H])
            base, current = ray(chart, True), ray(chart)
            eb = base['lam']*float(np.asarray(base['profile'].fields([H]), float)[2, 0])
            ec = current['lam']*float(np.asarray(current['profile'].fields([H]), float)[2, 0])
            data.append([eb, ec-eb])
        data = np.asarray(data).T
        count = element['degree']+1
        for part, table in enumerate(tables):
            for index, sign, coefficient in zip(element['ids'][:count], element['sign'][:count], transform@data[part]):
                index, coefficient = int(index), float(sign*coefficient)
                if index in table:
                    audit.near(table[index], coefficient, 'retained_cap_shared_projected_trace_mode', 3e-10)
                table[index] = coefficient
        if element['lower'][0] >= raw['transition_controls'][4]-1e-11:
            baseline_ids.update(map(int, element['ids'][:count]))
    constraints = {int(row): value for row, value in raw['constraints']}
    for index in baseline_ids:
        tables[1][index] = 0.
    rows = sorted(tables[0])
    projected = [tables[0][i]+tables[1][i] for i in rows]
    audit.near([constraints[offsets[d+1]+i] for i in rows], projected,
               'retained_cap_actual_projected_constraint_trace', 3e-10)
    audit.near(state[offsets[d+1]+np.asarray(rows, int)], projected,
               'retained_cap_actual_lifted_projected_trace', 3e-10)
    files = [validation_path, receipt_path, *pin_paths, raw_path,
             directory/'pressure-active/case.xml', directory/'pressure-active/run.log']
    return dict(producing_source_commit=producing_source, native_calls=0,
        input_files={str(path): dict(bytes=path.stat().st_size, sha256=sha(path)) for path in files},
        pressure_probe_sha256=sha(program), base_residual_error=residual_error, alpha=raw['alpha'],
        weak=metadata, capped_rays=cap_count, interior_chi_rays=interior_chi_count,
        matching_trace_rows=len(rows), rays=ray_records,
        scope='Analysis of retained same-binary pressure-active capped operator, exact time mass and projected trace. Distinct from the legal uncapped chi=0 depletion fixture; no cap-active chi=0 native claim or solver trajectory.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--kernel-probe', required=True, type=Path)
    parser.add_argument('--wall-probe', required=True, type=Path)
    parser.add_argument('--template', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--reuse-run', type=Path,
                        help='Reanalyze exact retained kernel/native files; execute no probes unless --continue-missing is supplied')
    parser.add_argument('--continue-missing', action='store_true',
                        help='Execute only missing fixed fixtures with the identical retained probe binaries')
    parser.add_argument('--producing-source-commit', default='',
                        help='Actual compiled probe source epoch, distinct from the current Python postprocessor epoch')
    parser.add_argument('--capped-pressure-run', required=True, type=Path,
                        help='Passed same-producing-epoch pressure-suite results; retained cap full residual/trace analysis only')
    args = parser.parse_args()
    if platform.system() != 'Linux':
        raise SystemExit('Run numerical verification only on the authoritative Linux host')
    if args.continue_missing and args.reuse_run is None:
        raise SystemExit('--continue-missing requires the retained --reuse-run')
    if args.reuse_run is not None and not args.producing_source_commit:
        raise SystemExit('Retained replay requires explicit --producing-source-commit')
    output = reserve(args.output)
    started = time.monotonic()
    report = dict(schema=1, phase='M8-R2 B compiled inner law and ACTIVE native weak operator',
        started_utc=stamp(), numerical_workers=1, numerical_threads=THREADS,
        helper_sha256=sha(__file__), kernel_probe_sha256=sha(args.kernel_probe),
        wall_probe_sha256=sha(args.wall_probe), template_sha256=sha(args.template),
        commands=[], cases=[], failures=[], reused_inputs={}, unperformed_fixtures=[],
        producing_source_commit=args.producing_source_commit or None,
        trajectory_launched=False, solver_stage_launched=False)
    report['imported_source_files'] = {str(Path(module.__file__).resolve()): sha(module.__file__)
        for module in list(sys.modules.values()) if getattr(module, '__file__', None)
        and '/Tests/RANS/' in str(Path(module.__file__).resolve())}
    audit = Verification()

    def run(command, directory):
        directory.mkdir(parents=True, exist_ok=True)
        entry = dict(command=list(map(str, command)), cwd=str(directory))
        report['commands'].append(entry)
        with (directory/'run.log').open('w') as log:
            completed = subprocess.run(entry['command'], cwd=directory, env=os.environ.copy(),
                                       stdout=log, stderr=subprocess.STDOUT, timeout=900)
        entry.update(exit_code=completed.returncode, log_sha256=sha(directory/'run.log'))
        audit.require(completed.returncode == 0, 'Compiled probe failed; see '+str(directory/'run.log'))

    try:
        if args.reuse_run is None:
            kernel_file = output/'kernel.json'
            run([args.kernel_probe.resolve(), '--output', kernel_file.resolve()], output/'kernel-run')
        else:
            previous_path = args.reuse_run/'validation.json'
            previous = json.loads(previous_path.read_text())
            audit.require({'kernel_probe_sha256', 'wall_probe_sha256'} <= set(previous),
                          '--reuse-run must be the verifier results directory containing kernel.json and the probe-hash validation.json, not the supervisor directory')
            audit.require(sha(args.kernel_probe) == previous['kernel_probe_sha256'] and
                          sha(args.wall_probe) == previous['wall_probe_sha256'],
                          'Retained continuation requires identical compiled probe bytes; use a new explicit epoch otherwise')
            kernel_file = args.reuse_run/'kernel.json'
            report['reused_inputs'][str(previous_path)] = sha(previous_path)
            report['reused_inputs'][str(kernel_file)] = sha(kernel_file)
            report['retained_original_success'] = previous.get('success')
            report['retained_original_success_field_present'] = 'success' in previous
            report['retained_original_failures'] = previous.get('failures', [])
            report['kernel_native_reexecution'] = False
        report['kernel'] = kernel(json.loads(kernel_file.read_text()), audit)
        report['kernel_sha256'] = sha(kernel_file)
        fixed = [('linear-active-chi0', 'Linear', 1, 0, 1),
                 ('saturated-active-chi0', 'SaturatedR2', 1, 0, 1),
                 ('saturated-active-chi1', 'SaturatedR2', 1, 1, 1),
                 ('linear-inactive-chi0', 'Linear', 2, 0, 1),
                 ('saturated-inactive-chi0', 'SaturatedR2', 2, 0, 1),
                 ('saturated-depleted-chi0', 'SaturatedR2', 3, 0, 1),
                 ('linear-zero-local-k', 'Linear', 4, 0, 1),
                 ('saturated-zero-local-k', 'SaturatedR2', 4, 0, 1)]
        raw_cases = {}
        for name, selector, fixture, chi, rK in fixed:
            retained = args.reuse_run/name/'native.json' if args.reuse_run is not None else None
            if retained is not None and retained.exists():
                raw_path, xml = retained, retained.parent/'case.xml'
                report['reused_inputs'][str(raw_path)] = sha(raw_path)
                report['reused_inputs'][str(xml)] = sha(xml)
                reexecuted = False
            elif args.reuse_run is not None and not args.continue_missing:
                report['unperformed_fixtures'].append(name)
                continue
            else:
                directory = output/name;directory.mkdir()
                xml = session(args.template, directory/'case.xml', selector, fixture, chi, rK)
                raw_path = directory/'native.json'
                run([args.wall_probe.resolve(), xml.resolve(), raw_path.resolve()], directory)
                reexecuted = True
            raw = json.loads(raw_path.read_text());raw_cases[name] = raw
            audit.require(raw['inner_stress_model'] == selector and raw['inner_stress_fixture'] == fixture
                          and raw['chi'] == chi and raw['parameters'][3] == rK,
                          'Retained/native fixture identity differs from the fixed suite: '+name)
            value = native(raw, audit)
            report['cases'].append(dict(name=name, input_raw=str(raw_path), input_xml=str(xml),
                xml_sha256=sha(xml), raw_sha256=sha(raw_path), native_reexecution=reexecuted,
                producing_source_commit=args.producing_source_commit or None, result=value))
            print(name+' independently verified', flush=True)
        if {'linear-inactive-chi0', 'saturated-inactive-chi0'} <= set(raw_cases):
            audit.near(raw_cases['linear-inactive-chi0']['lifted_state'],
                       raw_cases['saturated-inactive-chi0']['lifted_state'], 'same_actual_inactive_state', 0.)
            audit.near(raw_cases['linear-inactive-chi0']['raw_residual'],
                       raw_cases['saturated-inactive-chi0']['raw_residual'], 'inactive_native_operator_identity', 3e-13)
        if {'linear-zero-local-k', 'saturated-zero-local-k'} <= set(raw_cases):
            audit.near(raw_cases['linear-zero-local-k']['lifted_state'],
                       raw_cases['saturated-zero-local-k']['lifted_state'], 'same_actual_zero_fixture_state', 0.)
        for name, raw in raw_cases.items():
            audit.require(('M8-R2:InnerStress:SaturatedR2:v1:rho=0.8' in raw['physical_signature']) ==
                          (raw['inner_stress_model'] == 'SaturatedR2'), 'Selector/checkpoint identity '+name)
        report['retained_pressure_cap'] = retained_pressure_cap(
            args.capped_pressure_run, args.producing_source_commit, audit)
        report['available_input_reanalysis_success'] = True
        report['verification_complete'] = not report['unperformed_fixtures']
        report['success'] = report['verification_complete']
        report['status'] = 'passed' if report['verification_complete'] else 'available retained inputs verified; remaining fixtures explicitly unperformed'
    except Exception as error:
        report['success'] = False
        report['available_input_reanalysis_success'] = False
        report['verification_complete'] = False
        report['failures'].append(dict(type=type(error).__name__, message=str(error), traceback=traceback.format_exc()))
        print(traceback.format_exc(), file=sys.stderr)
    report.update(checks=audit.checks, errors=audit.errors, finished_utc=stamp(), seconds=time.monotonic()-started)
    changed = [path for path, value in report['imported_source_files'].items() if sha(path) != value]
    if changed:
        report['success'] = report['available_input_reanalysis_success'] = report['verification_complete'] = False
        report['failures'].append(dict(type='SourceChanged', paths=changed))
    report['scope'] = 'Constitutive and actual native weak/operator-action verification only. Does not satisfy fixed outer-box coverage, coupled-stage convergence, full discrete positivity, curved MMS, long-time or physical validation gates.'
    save(output/'validation.json', report)
    print(json.dumps(dict(success=report['success'], checks=audit.checks, output=str(output/'validation.json'))), flush=True)
    return 0 if not report['failures'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
