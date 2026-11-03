#!/usr/bin/env python3
"""Read-only bounded R1 controls analysis; never launches a native program.

Reuses the source-bound independent t=0 preflight through exact coefficient,
space, physical-input and raw-load identity. Independently recomputes every
recorded nonlinear residual norm and trial decision, native GLL scalar values,
and exported offending-point tensors. Internal Newton acceptance is never
called a converged physical stage or a completed time step.
"""
import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import xml.etree.ElementTree as ET

THREADS = {key: '1' for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
    'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')}
os.environ.update(THREADS)
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
sys.dont_write_bytecode = True


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def stamp():
    return datetime.now(timezone.utc).isoformat()


def canonical(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def vector_hash(value):
    return hashlib.sha256(struct.pack('<'+str(len(value))+'d', *value)).hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def compact_physical(value):
    if not value:
        return value
    result = dict(value)
    for key in ('outer_branch_pattern', 'limiter_branch_pattern', 'cap_branch_pattern'):
        if key in result:
            pattern = result.pop(key)
            result[key+'_summary'] = dict(characters=len(pattern), sha256=hashlib.sha256(pattern.encode()).hexdigest(),
                counts={letter: pattern.count(letter) for letter in sorted(set(pattern))})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root, out = args.root.resolve(), args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    source, build, controls = root/'source', root/'build', root/'controls-001'
    r1 = root.parent
    producer_path = controls/'receipt.json'
    producer = json.loads(producer_path.read_text())
    assert producer['status'] == 'bounded_controls_completed' and producer['runtime_unchanged']
    assert producer['source_commit'] == '360cbb7ac728b575cdf9af21e6d0963a1d88acb9'
    tool_pins = {str(path):sha(path) for path in Path(__file__).resolve().parent.glob('*.py')}
    lock = (root/'validation.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_SH)
    runtime_path = controls/'runtime-after.json'
    runtime = json.loads(runtime_path.read_text())
    for key, value in runtime['runtime_environment'].items():
        if value is None:
            os.environ.pop(key, None)
        else:
            os.environ[key] = value
    os.environ.update(THREADS)
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    os.environ['TMPDIR'] = str(root/'tmp')
    os.sched_setaffinity(0, {49})
    fp = source/'docs/plans/cls-epawf-revised-evidence/m7/fingerprint-runtime.py'
    subprocess.run(['python3', '-B', str(fp), '--root', str(root), '--output', str(out/'runtime-before.json'),
                    '--compare', str(runtime_path)], check=True)
    sys.path.insert(0, str(source/'solvers/IncNavierStokesSolver/Tests/RANS'))
    import numpy as np
    from cls_epawf_m8_r1_nonlinear import Audit, iterations
    from cls_epawf_m8_r1_preflight import inner_values
    from cls_epawf_m8_state_coverage import gll, tensor_points
    from cls_epawf_m6_operator_oracle import _segment_projection
    from cls_epawf_m8_r1_cases import physics, DIAGNOSTIC_INFO

    copied, input_pins = [], {}
    def read(path):
        path = Path(path); data = path.read_bytes(); input_pins[str(path)] = hashlib.sha256(data).hexdigest()
        return json.loads(data)
    def copy(path, name):
        path = Path(path); target = out/name; target.write_bytes(path.read_bytes())
        assert sha(target) == sha(path)
        copied.append(dict(file=name, original=str(path), sha256=sha(path), bytes=path.stat().st_size,
                           transformation='byte-for-byte copy'))
    report = dict(schema=1, milestone='M8-R1', phase='final controlled diagnostics', status='running',
        started_utc=stamp(), source_commit=producer['source_commit'], numerical_threads=THREADS, affinity=[49],
        solver_sha256=producer['executable_sha256'], producer_receipt=dict(path=str(producer_path), sha256=sha(producer_path)),
        producing_runtime=dict(before=str(controls/'runtime-before.json'), before_sha256=sha(controls/'runtime-before.json'),
            after=str(runtime_path), after_sha256=sha(runtime_path), unchanged=True),
        native_evolution_budget=dict(slots=[5, 6], additional_allowed=0, analyzer_solver_calls=0),
        qualification='All continuation unqualified from t=0; frozen local region already exited',
        cases=[], results=[], events=[], failures=[])
    focused = read(root/'focused-checks-001/receipt.json')
    assert focused['status'] == 'passed' and focused['runtime_unchanged']
    preflight = r1/'preflight-001/cases'
    preflight_audit = r1/'preflight-completion-001'
    checks_total = 0

    class ScalarRows:
        def __init__(self, audit):
            self.audit, self.raw, self.fields = audit, audit.raw, audit.fields
            self.d, self.off = audit.d, audit.off
            self.rows = []
            self.trace_ids = set()
            for e in self.fields[self.d+1].elements:
                if abs(e['lower'][2]-audit.H) < 1e-12:
                    self.trace_ids.update(e['ids'][:e['degree']+1])
            self.external = {int(row): value for row, value in self.raw['reference_constraints']
                if not (self.off[self.d+1] <= int(row) < self.off[self.d+2] and int(row)-self.off[self.d+1] in self.trace_ids)}
            for e in self.fields[self.d].elements:
                points = tensor_points(e, self.fields[self.d].axes, [gll(8), gll(8)])
                ids, phi, _ = self.fields[self.d].rows(e, points)
                epsilon_rows = None
                if e['lower'][2] >= audit.H:
                    ei, ep, _ = self.fields[self.d+1].rows(self.fields[self.d+1].by_id[e['geometry']], points)
                    epsilon_rows = ei, ep
                self.rows.append((e, points, ids, phi, epsilon_rows))
            self.source_x = []
            nodes = np.polynomial.legendre.leggauss(int(self.raw['quadrature_orders'][0]))[0]
            self.trace_x = []
            for e in self.fields[self.d+1].elements:
                if abs(e['lower'][2]-audit.H) < 1e-12:
                    a, b = e['lower'][0], e['upper'][0]
                    self.source_x.extend((a+(nodes+1)*(b-a)/2).tolist())
                    self.trace_x.extend((a+(_segment_projection(e['degree'])[0]+1)*(b-a)/2).tolist())
        def lifted(self, state):
            result = np.asarray(state, float).copy()
            for row, value in self.external.items(): result[row] = value
            return result
        def scan(self, state):
            effective = self.lifted(state)
            best = {key: None for key in ('inner_k', 'outer_k', 'outer_epsilon')}
            zero = 0
            for e, points, ids, phi, erows in self.rows:
                c = effective[self.off[self.d]+ids]
                constant = self.fields[self.d].constant()[ids]
                centre = float(c[0]/constant[0])
                delta = c-centre*constant
                values = np.asarray(np.longdouble(centre)+phi.astype(np.longdouble)@delta.astype(np.longdouble), float)
                ordinary = phi@c
                key = 'inner_k' if e['upper'][2] <= self.audit.H else 'outer_k'
                index = int(np.argmin(values))
                if key == 'inner_k': zero += int(np.count_nonzero(values == 0))
                if best[key] is None or values[index] < best[key]['value']:
                    bound = 64*np.finfo(float).eps*(abs(centre)+float(np.abs(phi[index])@np.abs(delta)))
                    best[key] = dict(value=float(values[index]), ordinary_dot=float(ordinary[index]),
                        chart=points[index].tolist(), geometry=e['geometry'], arithmetic_roundoff_scale=bound)
                if erows is not None:
                    ei, ep = erows; values = ep@effective[self.off[self.d+1]+ei]; index = int(np.argmin(values))
                    if best['outer_epsilon'] is None or values[index] < best['outer_epsilon']['value']:
                        best['outer_epsilon'] = dict(value=float(values[index]), chart=points[index].tolist(), geometry=e['geometry'])
            return dict(**best, exact_zero_inner_k_samples=zero,
                scope='Independent centered long-double Jacobi reconstruction at actual native GLL points; roundoff scale is not an interval certificate')
        def ray_kind(self, context):
            if not context or 'ray_chart' not in context: return 'not_exported'
            x = context['ray_chart'][0]
            trace = min(abs(x-value) for value in self.trace_x) < 1e-12
            source = min(abs(x-value) for value in self.source_x) < 1e-12
            return 'matching_trace_GLL_normalization_ray' if trace and not source else 'volume_source_Gauss_normalization_ray' if source and not trace else 'shared_or_unclassified_ray'
        def inner_point(self, state, point, identity=None):
            effective = self.lifted(state); point = np.asarray(point, float)
            identity = self.fields[0].element(point)['geometry'] if identity is None else identity
            G, k, K, nt, R = inner_values(self.raw, self.fields, effective, identity, point[None, :])
            S = G[0]+G[0].T; values, vectors = np.linalg.eigh(R[0]); vector = vectors[:, 0]
            if vector[np.argmax(np.abs(vector))] < 0: vector = -vector
            return dict(chart=point.tolist(), geometry=identity, k=float(k[0]), K=float(K[0]), k_over_K=float(k[0]/K[0]),
                nu_A=float(nt[0]), G=G[0].tolist(), S=S.tolist(), maximum_strain_eigenvalue=float(np.linalg.eigvalsh(S)[-1]),
                covariance_minimum=float(values[0]), covariance_eigenvector=vector.tolist(),
                r_R=float(3*nt[0]*np.linalg.eigvalsh(S)[-1]/(2*k[0])) if k[0] > 0 else None,
                transported_k_nonnegative=bool(k[0] >= 0),
                interpretation='Formal inner algebra outside physical k domain' if k[0] < 0 else 'Independent pointwise original inner AWF tensor')

    try:
        focused_results = []
        for name, pre in (('original-physical-stop', 'original-base'), ('resolved-local-stop', 'resolved-base')):
            p = root/'focused-checks-001'/name
            raw = read(p/'case.evolution.json'); before = read(preflight/pre/'case.evolution.json')
            event = raw['physical_range_events'][0]
            assert len(raw['physical_range_events']) == 1 and event['time'] == 0 and event['alpha'] == 0
            assert event['stopped_before_history_advance'] and raw['steps'] == 0 and raw['time'] == 0
            assert raw['accepted_states'] == [] and raw['audit_times'] == []
            identical = event['state'] == before['accepted_states'][0] and event['history'] == before['accepted_states'][0]
            history_identical = raw['native_history'] == before['native_history'] and raw['native_history_times'] == before['native_history_times']
            assert identical and history_identical and raw['rejected_state'] == event['state']
            focused_results.append(dict(name=name, guard=raw['physical_range_guard'], reasons=event['reasons'],
                XML_sha256=sha(p/'case.xml'), report_sha256=sha(p/'case.evolution.json'),
                matching_preflight_report_sha256=sha(preflight/pre/'case.evolution.json'),
                state_decoded_float64_LE_sha256=vector_hash(event['state']), history_decoded_float64_LE_sha256=vector_hash(event['history']),
                native_history_slot_sha256=[vector_hash(slot) for slot in raw['native_history']],
                initial_native_audit=event['audit'],
                independent_preflight_audit=dict(path=str(preflight_audit/(pre+'-audit.json')),
                    sha256=sha(preflight_audit/(pre+'-audit.json'))),
                time=0, steps=0, accepted_states=0, startup_history_and_event_state_identical_to_preflight=True,
                stopped_before_history_advance=True))
        save(out/'controls-focused-guard-identities.json', dict(schema=1, checks=focused_results, additional_solver_calls=0))

        for name, label in (('resolved-dt005', 'dt005'), ('resolved-dt0025', 'dt0025')):
            p = controls/'cases'/name; raw_path, trace_path = p/'case.evolution.json', p/'case.nonlinear.json'
            raw, trace = read(raw_path), read(trace_path)
            prior = read(preflight/'resolved-base/case.evolution.json')
            initial_audit = read(preflight_audit/'resolved-base-audit.json')
            xml_path = p/'case.xml'; input_pins[str(xml_path)] = sha(xml_path)
            settings = {node.get('PROPERTY'):node.get('VALUE') for node in ET.parse(xml_path).getroot().findall('CONDITIONS/SOLVERINFO/I')}
            assert all(settings.get(key)==value for key,value in DIAGNOSTIC_INFO.items())
            assert settings['RANSWallRestart']=='None' and settings['RANSWallCheckpoint']=='None'
            current_physics, prior_physics = physics(ET.parse(xml_path).getroot(), True), physics(ET.parse(preflight/'resolved-base/case.xml').getroot(), True)
            for key in DIAGNOSTIC_INFO:current_physics['solver_info'].pop(key.upper())
            current_physics['parameters'].pop('TIMESTEP');prior_physics['parameters'].pop('TIMESTEP')
            assert current_physics==prior_physics
            audit = Audit(raw, source); audit.gauge_variable_scale = audit.us/audit.H
            audit.require(raw['time'] == 0 and raw['steps'] == 0 and len(raw['accepted_states']) == 1,
                          'Control advanced despite the first-stage failure')
            audit.require(raw['stage_audits'] == [] and raw['step_audits'] == [], 'A converged stage was recorded unexpectedly')
            audit.require(raw['accepted_states'][0] == prior['accepted_states'][0], 'Initial coefficients changed')
            audit.require(raw['fields'] == prior['fields'] and raw['parameters'] == prior['parameters'], 'Physical space or inputs changed')
            audit.near(raw['raw_reference_loads'], prior['raw_reference_loads'], 'identical_initial_unanchored_raw_loads', 1e-12)
            audit.require(raw['native_history'] == prior['native_history'], 'Failed first stage changed startup history values')
            audit.near(raw['native_history_times'], [0., -raw['dt'], raw['dt'], raw['dt']], 'fresh_IMEX2_history_time_labels', 0.)
            audit.near(trace['time'], raw['dt'], 'first_stage_time', 0.)
            audit.near(trace['alpha'], 1/raw['dt'], 'IMEX1_startup_alpha_inside_native_IMEX2', 0.)
            audit.near(trace['history'], raw['accepted_states'][0], 'fresh_startup_implicit_history', 0.)
            audit.near(trace['dimensionless_variable_scales'], audit.scales(), 'actual_inverse_time_gauge_variable_scale')
            event = raw['physical_range_events'][0]
            audit.require(len(raw['physical_range_events']) == 1 and event['time'] == 0 and not event['stopped_before_history_advance'],
                          'Unexpected converged physical range stop')
            audit.require(event['reasons'] == ['outer_gradient_box_exit'], 'Unexpected initial range event')
            audit.require(event['state'] == raw['accepted_states'][0] and event['history'] == event['state'], 'Initial event payload differs')
            audit.require(trace['fresh_jacobian_retry'] is None and trace['fresh_jacobian'] is None,
                          'Unallocated same-stage retry was executed')
            matrix = np.asarray(trace['cached_jacobian'])
            builds = {row['jacobian_build'] for row in trace['iterations']}
            checked = iterations(trace['iterations'], audit, matrix if len(builds) == 1 else None)
            scalar = ScalarRows(audit); physical_events = []
            for native, row in zip(trace['iterations'], checked['iterations']):
                state = np.asarray(native['state']); scan = scalar.scan(state); physical = native.get('physical', {})
                row.update(state_decoded_float64_LE_sha256=vector_hash(state), residual_decoded_float64_LE_sha256=vector_hash(native['residual']),
                    step_decoded_float64_LE_sha256=vector_hash(native['linear_step']['step']), native_GLL_independent=scan)
                if physical.get('status') == 'evaluated_source_quadrature':
                    for key, field in (('minimum_native_inner_k','inner_k'), ('minimum_native_outer_k','outer_k'), ('minimum_native_outer_epsilon','outer_epsilon')):
                        audit.near(scan[field]['value'], physical[key], 'independent_'+key)
                    point = physical['worst_inner_point']
                    detail = scalar.inner_point(state, point[:3])
                    audit.near([detail['k'],detail['K'],detail['k_over_K'],detail['nu_A'],detail['maximum_strain_eigenvalue'],detail['covariance_minimum']], point[3:9], 'independent_internal_source_tensor')
                    audit.near(detail['G'], physical['worst_inner_gradient'], 'independent_internal_G')
                    audit.near(detail['S'], physical['worst_inner_strain'], 'independent_internal_S')
                    row['independent_worst_source_point'] = detail
                    row['initial_at_same_source_point'] = scalar.inner_point(raw['accepted_states'][0], point[:3])
                row['independent_native_inner_minimum_point'] = scalar.inner_point(state, scan['inner_k']['chart'], scan['inner_k']['geometry'])
                row['physical'] = compact_physical(row['physical'])
                for native_trial, trial in zip(native['trials'], row['trials']):
                    trial_state = state+native_trial['fraction']*np.asarray(native['linear_step']['step'])
                    trial['reconstructed_state_float64_LE_sha256'] = vector_hash(trial_state)
                    trial['state_hash_scope'] = 'Reconstructed from exported base/fraction/step; native trial vector not separately exported'
                    if native_trial['evaluation_succeeded']:
                        trial['physical'] = compact_physical(trial.get('physical'))
                    else:
                        context = native_trial.get('invalid_context')
                        trial['normalization_ray_family'] = scalar.ray_kind(context)
                        if context and 'chart' in context:
                            trial['base_at_same_guard_coordinate'] = scalar.inner_point(state, context['chart'])
                if not physical_events and physical.get('minimum_inner_covariance', 1.) < 0:
                    physical_events.append(dict(kind='first_negative_source_covariance_internal_Newton_iterate', iteration=native['iteration'],
                        stage_time=trace['time'], physical_stage_accepted=False, detail=row.get('independent_worst_source_point')))
            first_native_negative = next((row for row in checked['iterations'] if row['native_GLL_independent']['inner_k']['value'] < 0), None)
            exceptions = [trial for row in checked['iterations'] for trial in row['trials'] if trial['classification'] == 'unevaluable_residual_or_norm']
            invalid_ray_counts = {kind:sum(trial['normalization_ray_family']==kind for trial in exceptions)
                                  for kind in sorted({trial['normalization_ray_family'] for trial in exceptions})}
            audit.require(len(trace['iterations'][-1]['trials']) == 21 and not trace['iterations'][-1]['accepted'], 'Final search is incomplete')
            last = checked['iterations'][-1]
            summary = dict(name=name, dt=raw['dt'], stage_time=trace['time'], alpha=trace['alpha'],
                native_solver_seconds=next(row['seconds'] for row in producer['commands'] if row['name']==name),
                native_exit_code=next(row['exit_code'] for row in producer['commands'] if row['name']==name),
                last_accepted_time=0, accepted_states=1, completed_steps=0,
                classification='Newton residual-domain/globalization failure before any converged audited stage; not a physical guard stop',
                failure=raw['failure'], Newton_iterations=len(trace['iterations']), classification_counts=checked['classification_counts'],
                invalid_normalization_ray_counts=invalid_ray_counts,
                initial_event=dict(time=0, reasons=event['reasons'], stopped=False,
                    state_sha256=vector_hash(event['state']), history_sha256=vector_hash(event['history'])),
                first_internal_negative_covariance=physical_events[0] if physical_events else None,
                first_internal_native_GLL_negative_k_iteration=first_native_negative['iteration'] if first_native_negative else None,
                final_internal_native_GLL_k=last['native_GLL_independent']['inner_k'],
                final_internal_source_observation=last['physical'],
                final_all_fractions_unevaluable=all(not trial['evaluation_succeeded'] for trial in trace['iterations'][-1]['trials']),
                final_smallest_fraction=trace['iterations'][-1]['trials'][-1]['fraction'],
                first_invalid_trial=exceptions[0] if exceptions else None, last_invalid_trial=exceptions[-1] if exceptions else None,
                cached_Jacobian_builds=sorted(builds), cached_Jacobian_decoded_float64_LE_sha256=hashlib.sha256(matrix.astype('<f8').tobytes()).hexdigest(),
                cached_Jacobian_state_sha256=vector_hash(trace['cached_jacobian_state']),
                cached_Jacobian_history_sha256=vector_hash(trace['cached_jacobian_history']),
                cached_Jacobian_time=trace['cached_jacobian_time'],cached_Jacobian_iteration=trace['cached_jacobian_iteration'],
                last_Jacobian_age=last['jacobian_age'], final_scaled_residual=last['residual_norms'],
                initial_state_identical_to_refined_preflight=True, initial_history_values_identical_to_preflight=True,
                fresh_valid_dt_startup=True, checkpoint_exported=bool(raw.get('checkpoint')),
                initial_raw_unanchored_residual=initial_audit['result']['raw_unanchored_residual'],
                source_covariance_note='Internal nonlinear iterates are not accepted stages. Source-ray observations do not replace the stronger finite guard inventory.',
                no_fresh_J_comparison_in_this_control=True, no_refresh_policy_causality_proved=True,
                initial_independent_preflight=dict(path=str(preflight_audit/'resolved-base-audit.json'), sha256=sha(preflight_audit/'resolved-base-audit.json')),
                inputs=dict(XML=dict(path=str(p/'case.xml'), sha256=sha(p/'case.xml')),
                    report=dict(path=str(raw_path), sha256=sha(raw_path)), trace=dict(path=str(trace_path), sha256=sha(trace_path), bytes=trace_path.stat().st_size)),
                histories=dict(native_slot_sha256=[vector_hash(slot) for slot in raw['native_history']], native_times=raw['native_history_times'],
                    failed_stage_history_sha256=vector_hash(trace['history']), failed_iterate_sha256=vector_hash(trace['state'])),
                independent_checks=audit.checks, errors=audit.errors)
            save(out/('controls-'+label+'-iterations.json'), dict(schema=1, case=name,
                trace_sha256=sha(trace_path), original_policy=checked,
                interpretation='Every internal Newton iterate and frozen backtracking trial retained by hash/classification; complete native vectors and matrix remain in the original trace'))
            summary['iteration_ledger'] = 'controls-'+label+'-iterations.json'
            report['cases'].append(summary)
            report['results'].append({key:summary[key] for key in ('name','dt','classification','Newton_iterations','classification_counts',
                'invalid_normalization_ray_counts','first_internal_native_GLL_negative_k_iteration','final_smallest_fraction')})
            report['events'].extend([dict(case=name, time=0, type='initial_local_box_exit', accepted_initial_state=True, physical_stop=False),
                dict(case=name, time=trace['time'], type='first_stage_Newton_stop', completed_step=False,
                    iterations=len(trace['iterations']), physical_guard_stop=False)])
            checks_total += audit.checks
            save(out/'controls-ledger.json', report)
        report['focused_guard_identity_ledger'] = 'controls-focused-guard-identities.json'
        report['backend_exact_zero_k_caveat'] = dict(
            limitation='Native Audit allows inner k=0 but evaluates covariance only for k>0. The linear inner tensor remains defined at k=0 for positive sampled K; the guard is not a complete exact-zero covariance guarantee.',
            observed_exact_zero_native_GLL_samples=sum(row['native_GLL_independent']['exact_zero_inner_k_samples']
                for case in report['cases'] for row in read(out/case['iteration_ledger'])['original_policy']['iterations']),
            source_changed_after_run=False, no_field_floor_or_clipping=True,
            inference='Negative internal native-GLL samples are directly observed; no converged physical candidate reached the guard.')
        report['independent_checks'] = checks_total
        report['native_command_timings'] = producer['commands']
        report['status'] = 'analysis_passed_unsuccessful_controls_preserved'
    except Exception as error:
        report.update(status='analysis_failed', failures=report['failures']+[dict(type=type(error).__name__, message=str(error))])
    finally:
        subprocess.run(['python3', '-B', str(fp), '--root', str(root), '--output', str(out/'runtime-after.json'),
                        '--compare', str(out/'runtime-before.json')], check=True)
        after = json.loads((out/'runtime-after.json').read_text())
        report['analysis_runtime_unchanged'] = after['success'] and after['identical_qualification_inputs_and_runtime']
        report['input_hashes_unchanged'] = all(sha(path)==digest for path,digest in input_pins.items())
        report['frozen_tool_hashes_unchanged'] = all(sha(path)==digest for path,digest in tool_pins.items())
        report['frozen_tool_pins'] = tool_pins
        report['analysis_recipe_sha256'] = sha(__file__)
        report['source_module_pins'] = {str(Path(module.__file__).resolve()):sha(module.__file__) for module in list(sys.modules.values())
            if getattr(module,'__file__',None) and (str(Path(module.__file__).resolve()).startswith(str(source)) or
                Path(module.__file__).resolve().parent == Path(__file__).resolve().parent)}
        report['finished_utc'] = stamp()
        save(out/'controls-ledger.json', report)
    if report['status'] != 'analysis_passed_unsuccessful_controls_preserved':
        print(json.dumps(dict(status=report['status'], failures=report['failures'])), flush=True)
        return 1
    copy(producer_path, 'controls-runner.json')
    copy(controls/'cases/trajectory-plan.json', 'controls-plan.json')
    copy(root/'focused-checks-001/receipt.json', 'controls-focused-runner.json')
    copy(Path(__file__), 'controls-analysis.py')
    for name in ('resolved-dt005','resolved-dt0025'):
        copy(controls/'cases'/name/'case.xml', 'controls-'+('dt005' if name.endswith('dt005') else 'dt0025')+'.xml')
    save(out/'controls-curation-index.json',dict(schema=1,exact_copies=copied,
        raw_evidence_retained='Full original native reports, nonlinear traces, coefficient histories, Jacobians and source/runtime manifests remain on Linux; all paths/hashes are in controls-ledger.json'))
    lines=['| Case | Time | Event | Observation |','|---|---:|---|---|']
    for case in report['cases']:
        lines.append(f"| {case['name']} | 0 | Local region exit | Initial state accepted only as unqualified diagnostic; no physical guard stop |")
        counts=case['classification_counts']
        lines.append(f"| {case['name']} | {case['stage_time']} | Newton stop | {counts['accepted_descent']} internal descents, {counts['unevaluable_residual_or_norm']} exceptions, {counts['valid_non_descent']} valid non-descents; no completed step |")
    (out/'controls-events.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(dict(status=report['status'],checks=checks_total,results=report['results'],output=str(out))),flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
