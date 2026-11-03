#!/usr/bin/env python3
"""Focused compiled R2 kernel and ACTIVE actual EvolutionPatch verification.

The six fixed native fixtures are operator evaluations, never trajectories or
qualified physical controls.  Both inactive identity and active full weak
loads are reconstructed independently.  The capped fixture has rK=H solely
to exercise the source/trace mechanism, and is excluded from stage/control
qualification.  All numerical execution belongs on Linux, one worker.
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
from cls_epawf_m7_energy_oracle import Audit, _fields
from cls_epawf_m8_r2_preflight import weak_loads, selected_direct, selected_interface


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
    tensor_mesh(root, [0., 2., 4.], 3, 2, 1., 1., points=8)
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


def native(raw, audit):
    audit.require(raw['inner_stress_fixture'] in (1, 2, 3) and raw['unforced_alpha0_fixture']
                  and raw['alpha'] == 0 and not raw['upstream_transition'], 'Not the declared verification-only native fixture')
    fields, state = _fields(raw), np.asarray(raw['lifted_state'])
    expected, metadata = weak_loads(raw, fields, state)
    audit.near(raw['raw_residual'], expected, 'actual_complete_active_weak_operator', 2e-10)
    audit.near(np.asarray(raw['region_residuals'][0])+np.asarray(raw['region_residuals'][1]),
               expected, 'actual_regional_weak_operator_sum', 2e-10)
    old = Audit()
    direct = selected_direct(raw, fields, state, old)
    observed = raw['constitutive_audit']
    audit.near(observed['wall_sources'], direct['wall_sources'], 'actual_selected_wall_P_and_sources', 3e-10)
    audit.near(observed['robust_source'], direct['robust_source'], 'actual_capped_budget_and_trace', 3e-10)
    audit.near(observed['interface'], selected_interface(raw, fields, state), 'actual_full_normal_tangential_traction', 3e-10)
    audit.near(observed['source_volume_inner'][3], direct['inner_covariance_minimum'], 'actual_selected_covariance_including_zero', 3e-10)
    audit.near(np.asarray(raw['energy_budget'])[[2, 3, 4]], [direct['stress_work_minus_k_source'],
        direct['resolved_advection'], direct['roughness_force_work']], 'actual_same_stress_energy_cancellation', 3e-10)
    if raw['inner_stress_model'] == 'SaturatedR2':
        audit.near(observed['inner_stress'], direct['inner_stress'], 'actual_saturation_volume_and_range', 3e-10)
        if raw['inner_stress_fixture'] == 2:
            audit.require(observed['inner_stress'][0] == 0, 'Inactive fixture unexpectedly saturates')
        else:
            audit.require(observed['inner_stress'][0] > 0, 'ACTIVE native fixture did not exercise saturation')
    if raw['inner_stress_fixture'] == 3:
        audit.require(observed['robust_source'][5] > 0 and observed['robust_source'][6] < 1,
                      'Declared capped native fixture did not exercise the normalization cap')
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
        physical_signature=raw['physical_signature'],
        old_oracle_checks=old.checks, old_oracle_errors=old.errors,
        qualification='Operator verification only; these active states are outside the fixed outer box. No stage/trajectory or qualified physical solution is claimed.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--kernel-probe', required=True, type=Path)
    parser.add_argument('--wall-probe', required=True, type=Path)
    parser.add_argument('--template', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if platform.system() != 'Linux':
        raise SystemExit('Run numerical verification only on the authoritative Linux host')
    output = reserve(args.output)
    started = time.monotonic()
    report = dict(schema=1, phase='M8-R2 B compiled inner law and ACTIVE native weak operator',
        started_utc=stamp(), numerical_workers=1, numerical_threads=THREADS,
        helper_sha256=sha(__file__), kernel_probe_sha256=sha(args.kernel_probe),
        wall_probe_sha256=sha(args.wall_probe), template_sha256=sha(args.template),
        commands=[], cases=[], failures=[], trajectory_launched=False, solver_stage_launched=False)
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
        kernel_file = output/'kernel.json'
        run([args.kernel_probe.resolve(), '--output', kernel_file.resolve()], output/'kernel-run')
        report['kernel'] = kernel(json.loads(kernel_file.read_text()), audit)
        report['kernel_sha256'] = sha(kernel_file)
        fixed = [('linear-active-chi0', 'Linear', 1, 0, 1),
                 ('saturated-active-chi0', 'SaturatedR2', 1, 0, 1),
                 ('saturated-active-chi1', 'SaturatedR2', 1, 1, 1),
                 ('linear-inactive-chi0', 'Linear', 2, 0, 1),
                 ('saturated-inactive-chi0', 'SaturatedR2', 2, 0, 1),
                 ('saturated-capped-chi0', 'SaturatedR2', 3, 0, 2)]
        raw_cases = {}
        for name, selector, fixture, chi, rK in fixed:
            directory = output/name;directory.mkdir()
            xml = session(args.template, directory/'case.xml', selector, fixture, chi, rK)
            raw_path = directory/'native.json'
            run([args.wall_probe.resolve(), xml.resolve(), raw_path.resolve()], directory)
            raw = json.loads(raw_path.read_text());raw_cases[name] = raw
            value = native(raw, audit)
            report['cases'].append(dict(name=name, xml_sha256=sha(xml), raw_sha256=sha(raw_path), result=value))
            print(name+' independently verified', flush=True)
        audit.near(raw_cases['linear-inactive-chi0']['lifted_state'],
                   raw_cases['saturated-inactive-chi0']['lifted_state'], 'same_actual_inactive_state', 0.)
        audit.near(raw_cases['linear-inactive-chi0']['raw_residual'],
                   raw_cases['saturated-inactive-chi0']['raw_residual'], 'inactive_native_operator_identity', 3e-13)
        for name, raw in raw_cases.items():
            audit.require(('M8-R2:InnerStress:SaturatedR2:v1:rho=0.8' in raw['physical_signature']) ==
                          (raw['inner_stress_model'] == 'SaturatedR2'), 'Selector/checkpoint identity '+name)
        report['success'] = True
    except Exception as error:
        report['success'] = False
        report['failures'].append(dict(type=type(error).__name__, message=str(error), traceback=traceback.format_exc()))
        print(traceback.format_exc(), file=sys.stderr)
    report.update(checks=audit.checks, errors=audit.errors, finished_utc=stamp(), seconds=time.monotonic()-started)
    report['scope'] = 'Constitutive and actual native weak/operator-action verification only. Does not satisfy fixed outer-box coverage, coupled-stage convergence, full discrete positivity, curved MMS, long-time or physical validation gates.'
    save(output/'validation.json', report)
    print(json.dumps(dict(success=report['success'], checks=audit.checks, output=str(output/'validation.json'))), flush=True)
    return 0 if report['success'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
