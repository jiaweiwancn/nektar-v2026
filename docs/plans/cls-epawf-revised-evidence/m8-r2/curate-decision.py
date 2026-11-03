#!/usr/bin/env python3
"""Bind final R2 evidence and copy stage gates; stdlib metadata/hashes only.

No solver, scientific library, operator reconstruction, Git mutation or remote
command is invoked. Missing evidence stays pending. --final requires both
original-target stage records and report-bound independent analyses, but never
infers that a conditional frozen-scalar run occurred. --bounded-final instead
finalizes a source-bound terminal observation; it never supplies complete-stage
readiness, convergence or physical qualification from an incomplete record.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path


REVIEWED = '5186d63bd825aadd198488d93cff02fc9496a232'
PRODUCER = '9f20b3053c303e8805f5183f668a8ec468d9656b'
RUNTIME_SECTIONS = ('source_inputs', 'runtime_artifacts', 'runtime_bindings',
    'runtime_environment', 'python_packages', 'tool_images')
FROZEN = {
    'candidate-contract-v2.md': '77e20f90cf931f3572d9a9eed3bd9658eb3d2f44ca177b7e05e735a6a81c3bdd',
    'candidate-contract-v2.json': '7cbc782674178dfd7366df5e31a9152f0d2fa6c139fd25c4002f66a8f8524e80',
    'operational-policy-v1.md': 'c9a00043ae91812ea3a60702285576b51f3b0f7c0d52409fb7b99fe41410dab0',
    'operational-policy-v1.json': '1d630c518f1dadfcae4c7b6d0663c82180a2928f46d822bf71c465e9ea272f53',
}
SUMMARY_KEYS = ('schema', 'milestone', 'status', 'success', 'analysis_success',
    'checks', 'independent_checks', 'native_kernel_checks', 'native_fixtures',
    'source_commit', 'source_commit_end', 'producing_source_commit', 'producer_source',
    'analysis_source_commit', 'source_clean_start', 'source_clean_end', 'clean',
    'runtime_unchanged', 'full_runtime_source_inputs_unchanged', 'input_hashes_unchanged',
    'mode', 'suite', 'native_calls', 'operator_reconstructions', 'trajectory_calls',
    'weak_residual_calls', 'current_action_gate_passed', 'exact_affine_weak_continuity_gate_passed',
    'initial_derivative_gates_passed', 'analytic_controls_met', 'fixed_absolute_gate_met',
    'localization_classification', 'stored_J_minus_analytic', 'stored_J_roughness_support_error',
    'stored_J_outer_support_error', 'original_geometry_identity_gate',
    'diagnosis_completed', 'original_failed_full_analysis_remains_failed')


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for part in iter(lambda: stream.read(1024*1024), b''):
            h.update(part)
    return h.hexdigest()


def canonical(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
        ensure_ascii=True, allow_nan=False).encode()).hexdigest()


def compact_failure(value):
    if isinstance(value, dict):
        return {key: value[key] for key in ('label', 'where', 'type', 'message') if key in value}
    return value


def pick(value, keys):
    return {key: value[key] for key in keys if key in value}


def tri_all(values):
    if any(value is False for value in values):
        return False
    return True if values and all(value is True for value in values) else None


class Evidence:
    def __init__(self):
        self.files, self.roles, self.data, self.issues = {}, {}, {}, []

    def add(self, role, path, required=False, expected=None):
        path = Path(path).resolve(); name = str(path)
        self.roles[role] = name
        if name not in self.files:
            item = dict(path=name, available=path.is_file())
            self.files[name] = item
            if item['available']:
                before = path.stat()
                item.update(bytes=before.st_size, sha256=digest(path))
                after = path.stat()
                item['stable_during_hash'] = (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns)
                if path.suffix == '.json':
                    try:
                        value = json.loads(path.read_text())
                        self.data[name] = value
                        if isinstance(value, dict):
                            item['observed'] = pick(value, SUMMARY_KEYS)
                            failures = value.get('failures', [])
                            item['failures_count'] = len(failures) if isinstance(failures, list) else None
                            item['first_failures'] = [compact_failure(x) for x in failures[:3]] if isinstance(failures, list) else []
                    except (ValueError, OSError) as error:
                        item['parse_error'] = str(error)
                if not item['stable_during_hash']:
                    self.issues.append(dict(role=role, kind='evidence_changed_during_hash', path=name))
        item = self.files[name]
        if required and not item['available']:
            self.issues.append(dict(role=role, kind='required_evidence_missing', path=name))
        if required and item.get('parse_error') and not role.startswith('retained_attempt_'):
            self.issues.append(dict(role=role, kind='required_JSON_unreadable', path=name))
        if expected is not None:
            item['expected_sha256'] = expected
            item['expected_hash_matches'] = item.get('sha256') == expected
            if item['available'] and not item['expected_hash_matches']:
                self.issues.append(dict(role=role, kind='pinned_evidence_hash_mismatch', path=name))
        return self.data.get(name)

    def reference(self, role, value, required=False):
        if not isinstance(value, dict) or not value.get('path'):
            if required:
                self.issues.append(dict(role=role, kind='missing_reference_pin'))
            return None
        return self.add(role, value['path'], required, value.get('sha256'))


def runtime_binding(evidence, role, value, baseline):
    raw = evidence.reference(role, value, required=True)
    if not isinstance(raw, dict):
        return None
    flags = [raw.get('success') is True, raw.get('clean') is True, raw.get('source_commit') == PRODUCER]
    if not isinstance(baseline, dict) or not all(key in baseline and key in raw for key in RUNTIME_SECTIONS):
        return None
    flags.extend(raw.get(key) == baseline.get(key) for key in RUNTIME_SECTIONS)
    return all(flags)


def analysis_entries(evidence, paths, label, runtime):
    results = []
    for index, path in enumerate(paths):
        role = label+'_analysis_'+str(index)
        data = evidence.add(role, path, required=True)
        if not isinstance(data, dict):
            continue
        receipt = None
        if data.get('input_receipt'):
            receipt = evidence.add(role+'_receipt', data['input_receipt'], required=True)
        if receipt is None and 'inputs' in data and 'driver_sha256' in data:
            receipt = data
        producing_ok = None
        if isinstance(receipt, dict) and isinstance(runtime, dict) and all(key in runtime for key in RUNTIME_SECTIONS):
            producing_ok = (receipt.get('producing_source_commit') == PRODUCER and
                receipt.get('source_clean_start') is True and receipt.get('source_clean_end') is True and
                receipt.get('source_commit_end') == receipt.get('source_commit') and
                all(receipt.get('producing_runtime_sections', {}).get(key) == canonical(runtime[key])
                    for key in RUNTIME_SECTIONS))
        results.append(dict(role=role, data=data, receipt=receipt, producer_binding=producing_ok))
    return results


def bound_analysis(entry, case):
    receipt = entry.get('receipt')
    if not isinstance(receipt, dict):
        return False
    inputs = receipt.get('inputs', {})
    for key in ('report', 'trace', 'xml'):
        pin = case.get(key)
        if not isinstance(pin, dict) or not pin.get('sha256'):
            return False
        direct = inputs.get(pin['path'])
        if not isinstance(direct, dict) or direct.get('sha256') != pin['sha256']:
            return False
    return True


def analysis_target(entry):
    data = entry['data']
    levels = [row for row in data.get('stage_levels', []) if row.get('alpha') == 400 and row.get('attempted')]
    target = levels[0] if len(levels) == 1 else {}
    independent = target.get('independent', {})
    physical = independent.get('post_solve_physical', independent.get('physical', {}))
    common = physical.get('common_inventory_physical_qualification') or {}
    extra = physical.get('physical_ranges') or {}
    return dict(role=entry['role'], producer_binding=entry['producer_binding'],
        analysis_receipt_success=(entry.get('receipt') or {}).get('success'),
        analysis_success=data.get('analysis_success'), checks=data.get('checks'),
        current_action_gate_passed=data.get('current_action_gate_passed'),
        exact_affine_weak_continuity_gate_passed=data.get('exact_affine_weak_continuity_gate_passed'),
        target_converged=data.get('target_converged'),
        target_physical_range_passed=data.get('target_physical_range_passed'),
        scalar_domain_admissible=independent.get('scalar_domain', {}).get('admissible'),
        common_inventory_qualification=pick(common, ('fixed_outer_box_covered', 'strict_covariance_nonnegative',
            'declared_native_covariance_roundoff_guard', 'minimum_lambda', 'qualified_on_common_inventory')),
        additional_declared_finite_inventory=extra.get('finite_joint_preflight_passed'),
        target_analysis_failed=independent.get('analysis_failed'),
        failures_count=len(data.get('failures', [])),
        first_failures=[compact_failure(x) for x in data.get('failures', [])[:3]],
        complete_target_analysis=bool(target) and isinstance(data.get('protocol'), dict))


def stage_family(evidence, summary_path, paths, label, runtime):
    if summary_path is None:
        return dict(supplied=False, cases=[], analyses_supplied=len(paths))
    data = evidence.add(label+'_summary', summary_path, required=True)
    if not isinstance(data, dict):
        return dict(supplied=True, available=False, cases=[])
    evidence.reference(label+'_native_receipt', data.get('run_receipt'), required=True)
    before = runtime_binding(evidence, label+'_runtime_before', data.get('runtime_before'), runtime)
    after = runtime_binding(evidence, label+'_runtime_after', data.get('runtime_after'), runtime)
    entries = analysis_entries(evidence, paths, label, runtime)
    result = dict(supplied=True, mode=data.get('mode'), producer_source=data.get('producer_source'),
        producer_source_end=data.get('producer_source_end'), source_clean=data.get('source_clean'),
        runtime_unchanged=data.get('runtime_unchanged'), runtime_before_matches_producer=before,
        runtime_after_matches_producer=after, target=data.get('target'),
        allowed_alpha_levels=data.get('allowed_alpha_levels'), cases=[])
    for index, case in enumerate(data.get('cases', [])):
        for key in ('xml', 'report', 'trace'):
            # Bind the exact pins already checked by the native curator and
            # independent analyzer; do not reread multi-gigabyte traces here.
            if case.get(key):
                result.setdefault('native_file_pins', []).append(dict(case=case.get('inner_stress'), kind=key, **case[key]))
        targets = [row for row in case.get('levels', []) if row.get('alpha') == 400 and row.get('attempted')]
        target = targets[0] if len(targets) == 1 else {}
        matches = [entry for entry in entries if bound_analysis(entry, case)]
        checked = [analysis_target(entry) for entry in matches]
        domain = target.get('scalar_domain', {})
        row = dict(inner_stress=case.get('inner_stress'), mode=data.get('mode'),
            native_exit_code=case.get('native_exit_code'), native_seconds=case.get('native_seconds'),
            target_attempted=case.get('target_attempted'), target_record_present=bool(target),
            native_target_true_coupled_convergence=case.get('target_true_coupled_convergence'),
            native_target_physical_range_flag=case.get('target_native_physical_range_flag'),
            native_time=case.get('native_time'), native_steps=case.get('native_steps'),
            accepted_physical_states=case.get('accepted_physical_states'),
            initial_state_sha256=case.get('initial_state_float64_LE_sha256'),
            target=pick(target, ('alpha', 'time', 'converged', 'true_coupled_residual_converged',
                'frozen_velocity_pressure_diagnostic', 'iterations', 'failure', 'residual_verification_failure',
                'qualification_audit_failure', 'scaled_block_norms', 'physical_range_passed',
                'native_history_advanced', 'state_sha256', 'history_sha256', 'residual_sha256')),
            scalar_domain=pick(domain, ('admissible', 'phase', 'failure', 'minimum_dimensionless_margin',
                'limiting_inventory', 'limiting_point', 'minimum_inner_k', 'minimum_outer_k',
                'minimum_outer_epsilon', 'minimum_sampled_K', 'minimum_projected_trace_epsilon', 'minimum_lambda')),
            native_target_audit=pick(target.get('audit', {}), ('r2_qualified', 'r2_outer_range', 'r2_covariance',
                'robust_source', 'inner_stress', 'inner_constitutive_modification', 'interface', 'energy')),
            independent_analyses=checked, report_bound_analysis_count=len(matches))
        frozen = data.get('mode') == 'frozen-scalar' or target.get('frozen_velocity_pressure_diagnostic') is True
        row['frozen_scalar_is_not_a_coupled_root'] = frozen
        eligible = []
        for check in checked:
            eligible.append(tri_all([check['producer_binding'], check['analysis_receipt_success'], check['analysis_success'],
                check['current_action_gate_passed'], check['exact_affine_weak_continuity_gate_passed'],
                check['scalar_domain_admissible'], check['common_inventory_qualification'].get('qualified_on_common_inventory'),
                check['additional_declared_finite_inventory']]))
        independent = True if any(x is True for x in eligible) else False if eligible and all(x is False for x in eligible) else None
        row['all_supplied_metadata_qualification_gates'] = (False if frozen else tri_all([
            row['native_target_true_coupled_convergence'], row['native_target_physical_range_flag'],
            domain.get('admissible'), independent, before, after, data.get('runtime_unchanged')]))
        if not target or not checked:
            row['result_from_available_records'] = 'pending original-target record or report-bound final independent analysis'
        elif frozen:
            row['result_from_available_records'] = 'frozen-scalar diagnostic only; inspect its solved scalar convergence separately'
        elif row['all_supplied_metadata_qualification_gates'] is True:
            row['result_from_available_records'] = 'original-target root passes all supplied finite-set metadata gates'
        elif row['native_target_true_coupled_convergence'] is True:
            row['result_from_available_records'] = 'native target converged; independent/physical qualification is failed or incomplete'
        else:
            row['result_from_available_records'] = 'bounded target search did not supply an independently qualified coupled root'
        result['cases'].append(row)
    result['unmatched_analysis_roles'] = [entry['role'] for entry in entries if not any(bound_analysis(entry, case) for case in data.get('cases', []))]
    return result


def stage_assessment(evidence, path, runtime):
    """Import the strict assessor's terminal evidence without rerunning it."""
    if path is None:
        return dict(supplied=False, bounded_evidence_ready=False, cases=[]), []
    start = len(evidence.issues)
    wrapper = evidence.add('strict_stage_assessment', path, required=True)
    if not isinstance(wrapper, dict):
        return dict(supplied=True, bounded_evidence_ready=False, cases=[]), ['Terminal stage assessment unavailable']
    data = wrapper.get('assessment', {})
    pinned = lambda value: (isinstance(value, dict) and bool(value.get('path')) and
        isinstance(value.get('sha256'), str) and len(value['sha256']) == 64 and
        all(char in '0123456789abcdef' for char in value['sha256']))
    if not all(pinned(value) for value in [wrapper.get('assessor')]+[
            data.get(key) for key in ('original_receipt', 'preflight', 'original', 'stage_plan')]):
        evidence.issues.append(dict(kind='strict_stage_required_reference_hash_missing'))
    evidence.reference('strict_stage_assessor', wrapper.get('assessor'), required=True)
    receipt = evidence.reference('strict_stage_original_receipt', data.get('original_receipt'), required=True)
    for key in ('preflight', 'original', 'stage_plan'):
        evidence.reference('strict_stage_'+key, data.get(key), required=True)
    # The strict assessor deliberately uses UTF-8, ensure_ascii=False hashes.
    # Keep the existing analyzer/runtime canonical encoding above unchanged.
    unicode_hash = lambda value: hashlib.sha256(json.dumps(value, sort_keys=True,
        ensure_ascii=False, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
    descriptor = data.get('runtime', {})
    runtime_ok = (isinstance(runtime, dict) and all(key in runtime for key in RUNTIME_SECTIONS) and
        descriptor.get('source_commit') == PRODUCER and
        descriptor.get('source_root') == runtime.get('source_root') and
        descriptor.get('build_root') == runtime.get('build_root') and
        all(descriptor.get('sections', {}).get(key) == unicode_hash(runtime[key]) for key in RUNTIME_SECTIONS))
    solver_key = 'build/IncNavierStokesSolver'
    runtime_ok = bool(runtime_ok and descriptor.get('solver', {}).get('sha256') ==
        runtime.get('runtime_artifacts', {}).get(solver_key))
    if not (wrapper.get('success') is True and wrapper.get('failures') == [] and
            wrapper.get('native_calls') == 0 and data.get('source_commit') == PRODUCER and
            data.get('mode') == 'matched' and data.get('family_status') in ('completed', 'incomplete', 'failed') and
            isinstance(receipt, dict) and receipt.get('status') == data.get('family_status') and
            receipt.get('source_commit') == PRODUCER and receipt.get('finished_utc') and runtime_ok):
        evidence.issues.append(dict(kind='strict_stage_assessment_binding_or_terminal_scope_failed'))
    rows, unknowns = [], []
    for case in data.get('cases', []):
        law = case.get('inner_stress')
        if (case.get('target_attempted') is not True or 'converged_target_evidence' not in case or
                type(case.get('converged_target_evidence')) not in (bool, type(None)) or
                not isinstance(case.get('issues'), list) or case.get('issues')):
            evidence.issues.append(dict(kind='strict_stage_case_attempt_or_status_unresolved',
                inner_stress=law, issues=case.get('issues')))
        for key in ('xml', 'log'):
            if not pinned(case.get(key)):
                evidence.issues.append(dict(kind='strict_stage_case_reference_hash_missing', inner_stress=law, reference=key))
            evidence.reference('strict_stage_'+str(law)+'_'+key, case.get(key), required=True)
        trace = case.get('trace_evidence', {})
        row = pick(case, ('inner_stress', 'preflight_h0_sha256', 'unknowns', 'evidence_status',
            'target_start_observed_in_log', 'complete_target_solve_record', 'target_attempted',
            'converged_target_evidence', 'full_coupled_convergence_evidence', 'scalar_solve_convergence_evidence',
            'wall_cap_stop_observed', 'wall_cap_interpretation', 'native_postqualification_report_available',
            'native_history_report_available', 'native_history_unchanged_verified', 'issues'))
        row.update(native_command=case.get('native_command'), log_evidence=case.get('log_evidence'),
            final_report_evidence=case.get('final_report_evidence'),
            complete_prefix=pick(trace, ('complete_prefix_bytes', 'complete_prefix_sha256',
                'complete_records', 'unpublished_tail')),
            retained_native_pins=pick(case, ('report', 'trace', 'matrices')),
            independent_complete_target_analysis_supplied_by_assessment=False)
        if case.get('complete_target_solve_record') is not True:
            unknowns.append(str(law)+': complete original-target solve record unavailable')
        if case.get('converged_target_evidence') is None:
            unknowns.append(str(law)+': original-target convergence status unavailable')
        if case.get('native_postqualification_report_available') is not True:
            unknowns.append(str(law)+': original-target physical/accounting postqualification report unavailable')
        if case.get('native_history_unchanged_verified') is None:
            unknowns.append(str(law)+': final native-history equality not verified from a report')
        rows.append(row)
    if (len(rows) != 2 or {row.get('inner_stress') for row in rows} != {'Linear', 'SaturatedR2'} or
            len({row.get('preflight_h0_sha256') for row in rows}) != 1 or
            any(not row.get('preflight_h0_sha256') for row in rows)):
        evidence.issues.append(dict(kind='strict_stage_matched_family_or_history_identity_failed'))
    ready = len(evidence.issues) == start
    positive = any(row.get('converged_target_evidence') is True for row in rows)
    incomplete = any(row.get('complete_target_solve_record') is not True for row in rows)
    outcome = ('pending terminal source-bound stage assessment' if not ready else
        'positive target-convergence evidence retained; qualification remains separate' if positive else
        'bounded target attempts observed; completion unavailable' if incomplete else
        'completed bounded target solve records report no convergence'
            if rows and all(row.get('converged_target_evidence') is False for row in rows) else
        'bounded target observations recorded; convergence/qualification remains pending')
    return dict(supplied=True, bounded_evidence_ready=ready, bounded_task_outcome=outcome,
        source_commit=data.get('source_commit'), family_status=data.get('family_status'),
        original_receipt=data.get('original_receipt'), runtime_binding_passed=runtime_ok,
        bound_file_count=len(data.get('bound_files', {})),
        bound_file_inventory_sha256=unicode_hash(data.get('bound_files', {})),
        bound_file_policy='Pins are cryptographically bound by the strict assessment; huge matrices/traces are not rehashed here',
        conditional_frozen_scalar_allowed=data.get('conditional_frozen_scalar_allowed'),
        fallback_execution_inferred=False, cases=rows), unknowns


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True, help='Linux m8-r2 experiment root')
    parser.add_argument('--documents', type=Path, required=True, help='Frozen/current m8-r2 evidence document directory')
    parser.add_argument('--stage-summary', type=Path)
    parser.add_argument('--stage-analysis', type=Path, action='append', default=[])
    parser.add_argument('--stage-assessment', type=Path)
    parser.add_argument('--scalar-summary', type=Path)
    parser.add_argument('--scalar-analysis', type=Path, action='append', default=[])
    parser.add_argument('--retained-attempt', type=Path, action='append', default=[])
    parser.add_argument('--extra', action='append', default=[], metavar='ROLE=PATH')
    closure = parser.add_mutually_exclusive_group()
    closure.add_argument('--final', action='store_true')
    closure.add_argument('--bounded-final', action='store_true',
        help='Finalize a strict bounded observation, without complete-stage or physical qualification')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Preserve existing decision ledgers: '+str(args.output))
    root, docs, e = args.root.resolve(), args.documents.resolve(), Evidence()
    epoch = root/'epoch-002'
    for name, expected in FROZEN.items():
        e.add('frozen_'+name, docs/name, True, expected)
    e.add('preserved_R1_contract', docs.parent/'m8-r1/candidate-contract.md', True,
        'c635154db92a07be5bf27b5b50817340e9d6b5301d778c6a1099ac1d220c8879')
    e.add('preserved_R1_preregistered', docs.parent/'m8-r1/kernel-contract-preregistered.md', True,
        '082a791fb086a08bd83285e68ae5e05923fdb688c7919410263bb894a6f22313')
    for name in ('operational-policy-addendum-v2.md', 'operational-policy-addendum-v2.json',
            'operational-policy-addendum-v3.md', 'inner-stress-caller-audit.md', 'pressure-action-origin.md',
            'inner-verification-result.json', 'inner-zero-ratio-diagnosis.json', 'cap-fixture-screen.json',
            'analyzer-static-result.json', 'preflight-001-stop.json', 'preflight-002-summary.json'):
        e.add('document_'+name, docs/name, required=True)
    canonical_inputs = {
        'build': epoch/'build-receipt.json', 'runtime': epoch/'runtime-fingerprint.json',
        'inner_wrapper': epoch/'inner-reanalysis-004/receipt.json',
        'inner_validation': epoch/'inner-reanalysis-004/results/validation.json',
        'pressure_wrapper': epoch/'pressure-verification-001/receipt.json',
        'pressure_validation': epoch/'pressure-verification-001/results/validation.json',
        'historical_pressure_reanalysis': root/'pressure-reanalysis-002/validation.json',
        'legacy': epoch/'legacy-verification-001/receipt.json',
        'static': epoch/'analyzer-static-002/receipt.json',
        'preflight_wrapper': epoch/'preflight-002/receipt.json',
    }
    raw = {name: e.add(name, path, required=True) for name, path in canonical_inputs.items()}
    for name, value in raw.items():
        if not isinstance(value, dict):
            e.issues.append(dict(role=name, kind='canonical_JSON_unavailable'))
    runtime = raw['runtime']
    if isinstance(runtime, dict) and runtime.get('source_commit') != PRODUCER:
        e.issues.append(dict(kind='wrong_producing_source', actual=runtime.get('source_commit'), expected=PRODUCER))
    if isinstance(runtime, dict) and (runtime.get('success') is not True or runtime.get('clean') is not True or
            not all(key in runtime for key in RUNTIME_SECTIONS)):
        e.issues.append(dict(kind='incomplete_or_failed_producing_runtime_pin'))
    producer = dict(commit=PRODUCER, runtime_sections={name: canonical(runtime[name]) for name in RUNTIME_SECTIONS}
        if isinstance(runtime, dict) and all(name in runtime for name in RUNTIME_SECTIONS) else None)
    if isinstance(runtime, dict):
        keys = [key for key in runtime.get('runtime_artifacts', {}) if Path(key).name == 'IncNavierStokesSolver']
        if len(keys) == 1:
            key = keys[0]
            producer['solver'] = dict(path=runtime.get('runtime_resolved_paths', {}).get(key), sha256=runtime['runtime_artifacts'][key])
        else:
            e.issues.append(dict(kind='native_solver_image_not_unique_in_runtime_pin', count=len(keys)))
    attempts = set(args.retained_attempt)
    for pattern in ('*/*receipt*.json', '*/*/*receipt*.json'):
        attempts.update(path for path in root.glob(pattern) if not {'source', 'build', '.git'}.intersection(path.relative_to(root).parts))
    for index, path in enumerate(sorted(attempts, key=str)):
        e.add('retained_attempt_'+str(index), path, required=True)
    for item in args.extra:
        name, path = item.split('=', 1)
        e.add('extra_'+name, path, required=True)
    matched = stage_family(e, args.stage_summary, args.stage_analysis, 'matched', runtime)
    scalar = stage_family(e, args.scalar_summary, args.scalar_analysis, 'scalar', runtime)
    assessed, remaining_unknowns = stage_assessment(e, args.stage_assessment, runtime)
    contract = e.data.get(str((docs/'candidate-contract-v2.json').resolve()), {})
    if contract.get('reviewed_source') != REVIEWED:
        e.issues.append(dict(kind='reviewed_baseline_disagrees_with_frozen_contract'))
    b_summary = e.data.get(str((docs/'inner-verification-result.json').resolve()), {})
    e_summary = e.data.get(str((docs/'preflight-002-summary.json').resolve()), {})
    gate_overview = dict(
        B_inner=dict(reported_status=b_summary.get('status'), independent_checks=b_summary.get('independent_checks'),
            native_kernel_checks=b_summary.get('native_kernel_checks'), validation_success=(raw['inner_validation'] or {}).get('success'),
            scope='Focused kernel/operator evidence; unresolved ratio diagnostic is retained in its supplied status'),
        C_pressure=pick(raw['pressure_validation'] or {}, ('success', 'checks', 'fixed_action_tolerance')),
        E_initial=pick(e_summary, ('runner_status', 'outer6_coverage_decision', 'outer9_needed',
            'final_provenance_passed', 'source_runtime_inputs_unchanged', 'input_hashes_unchanged')),
        static=pick(raw['static'] or {}, ('success', 'checks', 'scope')),
        focused_legacy=pick(raw['legacy'] or {}, ('status', 'source_commit', 'runtime_unchanged')))
    final_rows = matched.get('cases', [])
    final_ready = (matched.get('supplied') is True and len(final_rows) == 2 and
        {row.get('inner_stress') for row in final_rows} == {'Linear', 'SaturatedR2'} and
        all(row.get('target_record_present') and row.get('target_attempted') is True and
            row.get('report_bound_analysis_count', 0) > 0 and any(check.get('complete_target_analysis')
                for check in row['independent_analyses']) for row in final_rows))
    if args.final and not final_ready:
        e.issues.append(dict(kind='final_original_target_evidence_incomplete',
            reason='Both matched alpha400 records and report-bound final independent analyses are required'))
    if args.bounded_final and not assessed.get('bounded_evidence_ready'):
        e.issues.append(dict(kind='bounded_terminal_observation_incomplete',
            reason='A source-bound terminal matched assessment with both actual target attempts and no issues is required'))
    if not final_ready:
        remaining_unknowns.append('Complete original-target independent residual/action/accounting qualification is not established')
    for path, item in e.files.items():
        if item.get('available') and item.get('sha256') != digest(path):
            e.issues.append(dict(kind='evidence_changed_during_curation', path=path))
    payload = dict(schema=1, milestone='M8-R2 decision ledger', created_utc=datetime.now(timezone.utc).isoformat(),
        curator=dict(path=str(Path(__file__).resolve()), sha256=digest(__file__)),
        method='stdlib metadata/hash extraction only; no scientific/native calls',
        scientific_or_native_calls=0, reviewed_baseline=REVIEWED, producing_epoch=producer,
        fixed_coefficients=contract.get('fixed_coefficients'), fixed_outer_box=contract.get('outer', {}).get('fixed_box'),
        historical_preservation_from_build=pick(raw['build'] or {}, ('preserved_sources', 'preserved_images',
            'preserved_sources_unchanged', 'preserved_images_unchanged')),
        final_requested=args.final, final_stage_evidence_ready=final_ready,
        bounded_final_requested=args.bounded_final,
        bounded_task_evidence_ready=assessed.get('bounded_evidence_ready') is True and not bool(e.issues),
        bounded_task_outcome=assessed.get('bounded_task_outcome', 'terminal stage assessment not supplied'),
        stage_assessment=assessed, remaining_unknowns=remaining_unknowns,
        ledger_integrity_passed=not bool(e.issues), issues=e.issues,
        evidence_roles=e.roles, evidence_files=list(e.files.values()),
        retained_receipt_inventory_count=len(attempts), matched_stage_gate_table=matched,
        prior_gate_overview=gate_overview,
        frozen_scalar_diagnostic=scalar, fallback_inferred=False,
        M8_closed=False, physical_trajectory_qualification=False,
        main_stage_result=('See bounded_task_outcome; complete-stage readiness remains false'
            if assessed.get('bounded_evidence_ready') else 'pending complete final target records/analyses') if not final_ready else
            'See matched stage gate table; native convergence, action accuracy and physical applicability remain separate',
        still_open=['genuine curved MMS', 'active-cap curved/variable-H verification',
            'fully discrete positivity', 'long-time/refinement qualification', 'scalability', 'paired terrain'],
        scope='No failed bounded search proves root absence; no auxiliary state is a physical trajectory sample. Missing evidence is not a passing gate.')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(payload, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps(dict(output=str(args.output), sha256=digest(args.output),
        integrity=payload['ledger_integrity_passed'], final_stage_evidence_ready=final_ready, issues=len(e.issues))))
    return 0 if payload['ledger_integrity_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
