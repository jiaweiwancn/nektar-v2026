#!/usr/bin/env python3
"""Bind completed frozen-scalar artifacts after the R2 stdout supervisor failure.

Read-only stdlib metadata/hashes only. A separate receipt preserves the failed
supervisor receipt and the lost native exit status/duration as unknown. The
fresh post-completion fingerprint must already have been captured externally.
No solver, fingerprint, scientific oracle, process or Git command is launched.
"""
import argparse
from datetime import datetime, timezone
import importlib.util
import json
from pathlib import Path
import sys


PRODUCER = '9f20b3053c303e8805f5183f668a8ec468d9656b'
KIND = 'R2FrozenScalarAdditiveCompletion:v1'


def evidence_module():
    path = Path(__file__).with_name('assess-stage-evidence.py')
    spec = importlib.util.spec_from_file_location('r2_completion_metadata_assessor', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def utc(value):
    result = datetime.fromisoformat(value)
    if result.tzinfo is None:
        raise ValueError('Recovery metadata timestamp lacks a timezone')
    return result


def collect(receipt_path, completion_check_path, expected_source=PRODUCER):
    """Deterministic recovery facts, rechecked independently by the assessor."""
    e = evidence_module()
    require = e.require
    receipt_path, completion_check_path = (Path(path).resolve() for path in
        (receipt_path, completion_check_path))
    run = receipt_path.parent
    bindings = e.Bindings()
    original = bindings.file(receipt_path)
    receipt = e.load(receipt_path)
    check_pin = bindings.file(completion_check_path)
    check = e.load(completion_check_path)
    require(receipt.get('schema') == 2 and receipt.get('status') == 'failed' and
        receipt.get('mode') == 'frozen-scalar' and receipt.get('source_commit') == expected_source and
        receipt.get('source_clean_start') is True and receipt.get('failure', {}).get('type') == 'BrokenPipeError',
        'Not the supported failed frozen-scalar stdout supervisor')
    require(not any(key in receipt for key in ('finished_utc', 'source_commit_end', 'source_clean_end', 'runtime_unchanged')),
        'Original supervisor already has terminal provenance; use the normal assessor')
    require(check.get('schema') == 1 and check.get('success') is True and check.get('native_calls') == 0 and
        check.get('completed_report_present_before_fingerprint') is True and
        check.get('no_producing_solver_process_observed') is True and check.get('artifacts_unchanged') is True and
        check.get('source_clean_start') is True and check.get('source_clean_end') is True and
        check.get('source_commit_start') == check.get('source_commit_end') == expected_source,
        'Missing successful fresh post-completion observation/source check')
    require(check.get('native_process_exit_code', 'missing') is None and
        check.get('native_process_elapsed_seconds', 'missing') is None,
        'Lost native process exit status/duration must remain null')
    require(check.get('numerical_threads') == {key:'1' for key in e.THREAD_KEYS} and
        isinstance(check.get('cpu_affinity'), list) and len(check['cpu_affinity']) == 1,
        'Completion check changed single-worker/thread scope')
    require(len(receipt.get('cases', [])) == 1 and receipt['cases'][0].get('inner_stress') == 'SaturatedR2' and
        receipt['cases'][0].get('status') == 'started', 'Unsupported original native case state')
    row = receipt['cases'][0]
    require(not any(key in row for key in ('native_exit_code', 'native_seconds', 'report_path', 'trace_path')),
        'Original native observations are already complete or contain fabricated process status')
    names = [command.get('name') for command in receipt.get('commands', [])]
    require(names == ['fingerprint-before', 'prepare', 'native-SaturatedR2', 'fingerprint-after'],
        'Unexpected/repeated original command family')
    commands = {command['name']: command for command in receipt['commands']}
    native = commands['native-SaturatedR2']
    require(not any(key in native for key in ('exit_code', 'seconds', 'log_path', 'log_sha256')),
        'Lost native command completion must not be replaced with an exit status/timing')
    for name in ('fingerprint-before', 'prepare', 'fingerprint-after'):
        command = commands[name]
        require(command.get('exit_code') == 0, 'Original metadata command did not complete: '+name)
        bindings.file(command['log_path'], command['log_sha256'])
    xml = Path(row['xml_path']).resolve()
    require(xml.is_relative_to(run/'cases') and e.sha(xml) == row['xml_sha256'], 'Original case XML changed')
    report_path, trace_path = xml.with_name('case.evolution.json'), xml.with_name('case.feasible.jsonl')
    log_path = run/'native-SaturatedR2.log'
    before_path = run/'runtime-before.json'
    required = {str(path) for path in (receipt_path, before_path, xml, report_path, trace_path, log_path)}
    preexisting = check.get('preexisting_completed_artifacts', [])
    require(len(preexisting) == len(required) and {item.get('path') for item in preexisting} == required,
        'Completion check does not bind the exact original/complete artifact set')
    for pin in preexisting:
        require(bindings.file(pin['path']) == pin, 'Completion-check artifact changed: '+pin['path'])
    before = e.load(before_path)
    post_pin = check.get('postcompletion_runtime')
    require(isinstance(post_pin, dict) and bindings.file(post_pin['path']) == post_pin,
        'Fresh post-completion runtime pin changed or is missing')
    post_path = Path(post_pin['path']).resolve()
    require(post_path.is_relative_to(completion_check_path.parent) and post_path != run/'runtime-after.json',
        'Early supervisor fingerprint is not a fresh post-completion fingerprint')
    post = e.load(post_path)
    require(post.get('compared_fingerprint') == str(before_path) and
        post.get('compared_fingerprint_sha256') == e.sha(before_path) and
        post.get('identical_qualification_inputs_and_runtime') is True and post.get('differences') == {},
        'Fresh fingerprint was not an exact comparison against actual before-runtime')
    require(utc(check['started_utc']) <= utc(post['captured_utc']) <= utc(check['finished_utc']),
        'Runtime fingerprint was not captured during the fresh completed-artifact check')
    fingerprint = check['fingerprint_command']
    cpu = check['cpu_affinity'][0]
    expected_fp = Path(before['source_root'])/'docs/plans/cls-epawf-revised-evidence/m7/fingerprint-runtime.py'
    require(fingerprint.get('argv') == ['taskset', '-c', str(cpu), 'python3', '-B', str(expected_fp),
        '--root', str(Path(before['build_root']).parent), '--output', str(post_path), '--compare', str(before_path)] and
        fingerprint.get('exit_code') == 0, 'Fresh fingerprint command/scope changed')
    require(bindings.file(fingerprint['log']['path']) == fingerprint['log'], 'Fresh fingerprint log changed')
    bindings.file(expected_fp)
    runtime = e.runtime_identity(before, post, receipt, bindings)
    early_after = bindings.file(run/'runtime-after.json')
    require(native.get('argv') == ['taskset', '-c', '0', str(Path(before['build_root'])/'IncNavierStokesSolver'), str(xml)] and
        Path(native['cwd']).resolve() == xml.parent and native.get('timeout_seconds') == e.WALL_CAP,
        'Original native launch differs from the frozen scalar case')

    # Reassess the exact original matched evidence using the normal path. The
    # stored authorization assessor remains bound by its original source hash;
    # an additive helper edit does not rewrite that historical authorization.
    gate = receipt['matched_evidence_gate']
    matched_pin, assessment_pin = gate['original_matched_receipt'], gate['assessment']
    require(bindings.file(matched_pin['path']) == matched_pin and
        receipt['matched_receipt'] == {key: matched_pin[key] for key in ('path', 'sha256')},
        'Original matched prerequisite receipt changed')
    require(bindings.file(assessment_pin['path']) == assessment_pin, 'Original matched assessment changed')
    stored = e.load(assessment_pin['path'])
    require(stored.get('success') is True, 'Original matched assessment did not succeed')
    require(bindings.file(stored['assessor']['path']) == stored['assessor'], 'Historical authorization assessor changed')
    matched = e.assess(matched_pin['path'], expected_source)
    require(e.canonical(matched) == e.canonical(stored['assessment']) == gate['independently_rechecked_metadata_sha256'] and
        matched['mode'] == 'matched' and matched['conditional_frozen_scalar_allowed'] is True and
        gate.get('both_original_targets_observed') is True and gate.get('both_targets_remain_unresolved') is True and
        gate.get('coupled_target_retry_permitted') is False,
        'Actual original matched targets no longer satisfy the conditional scalar prerequisite')
    for key, path_key, hash_key in (('preflight', 'preflight_path', 'preflight_sha256'),
                                    ('original', 'original_path', 'original_sha256')):
        require(matched[key] == bindings.file(receipt[path_key], receipt[hash_key]),
            'Frozen diagnostic changed matched '+key)
    require(matched['representation'] == receipt['representation'] and
        matched['runtime']['sections'] == runtime['sections'], 'Frozen diagnostic changed matched representation/runtime')
    preflight = e.load(receipt['preflight_path'])
    prior = [item for item in preflight['cases'] if item.get('inner_stress') == 'SaturatedR2' and item.get('quadrature') == 'base']
    require(len(prior) == 1, 'Ambiguous scalar base preflight')
    prior = prior[0]
    bindings.file(prior['xml'], prior['xml_sha256'])
    prior_key, xml_key = e.session_identity(prior['xml']), e.session_identity(xml)
    require(prior_key['stage_test'] == 'Off' and xml_key['stage_test'] == 'FrozenScalarR2' and
        prior_key['protected_xml_sha256'] == xml_key['protected_xml_sha256'], 'Frozen case changed protected preflight physics')
    bindings.file(prior['report'], prior['report_sha256'])
    initial = e.load(prior['report'])
    h0, fields, n = e.vector_sha(initial['initial']), e.canonical(initial['fields']), len(initial['initial'])
    require(h0 == prior['initial_state_float64_LE_sha256'] and n <= 3000, 'Actual preflight history/size changed')
    trace, referenced = e.trace_prefix(trace_path, 'SaturatedR2', h0, n, 'frozen-scalar', bindings)
    require(trace['unpublished_tail'] is None and len(trace['complete_records']) == 1 and
        trace['complete_records'][0]['alpha'] == e.TARGET_ALPHA, 'No unique complete frozen target record')
    raw = e.load(report_path)
    report = e.final_report(raw, 'SaturatedR2', h0, fields, 'frozen-scalar')
    require(report['target_level_present'] is True and report['target_attempted'] is True and
        report['native_history_unchanged_verified'] is True, 'Final native report/stage/history is incomplete')
    record = trace['complete_records'][0]
    level = next(level for level in raw['stage_test_levels'] if level['alpha'] == e.TARGET_ALPHA)
    require(e.vector_sha(level['state']) == record['state_sha256'] and
        level['converged'] is record['converged'] and level['scaled_block_norms'] == record['scaled_block_norms'],
        'Complete final stage and complete JSONL terminal record differ')
    matrices = []
    for path in sorted(xml.parent.glob('case.feasible.jsonl.solve*.J*.bin')):
        require(e.MATRIX_NAME.fullmatch(path.name) and path.stat().st_size == 8*n*n, 'Malformed completed matrix artifact')
        pin = bindings.file(path)
        matrices.append(dict(pin, shape=[n, n], referenced_by_complete_record=str(path.resolve()) in referenced))
    log_pin = bindings.file(log_path)
    report_pin, trace_pin = bindings.file(report_path), bindings.file(trace_path)
    bindings.finish()
    return dict(kind=KIND, original_receipt=original, completion_check=check_pin,
        producing_source_commit=expected_source, source_clean_postcompletion=True,
        before_runtime=bindings.file(before_path), early_supervisor_after_runtime=early_after,
        post_completion_runtime=post_pin, producing_runtime=runtime,
        matched_prerequisite=dict(original_receipt=matched_pin, original_assessment=assessment_pin,
            reassessed_metadata_sha256=e.canonical(matched), both_original_targets_remain_unresolved=True),
        native_command=dict(native, exit_code=None, seconds=None, log_path=str(log_path), log_sha256=log_pin['sha256']),
        native_case=dict(row, status='recovered_artifacts_observed', native_exit_code=None, native_seconds=None,
            report_path=str(report_path), report_sha256=report_pin['sha256'],
            trace_path=str(trace_path), trace_sha256=trace_pin['sha256'],
            target_attempted=True, target_converged=raw['stage_test_target_converged'],
            target_physical_range=raw['stage_test_target_physical_range']),
        complete_trace=trace, completed_final_report_metadata=report, matrices=matrices,
        native_log=log_pin, native_report=report_pin, native_trace=trace_pin,
        native_process_exit_code=None, native_process_elapsed_seconds=None,
        native_wall_cap_enforcement_after_supervisor_failure_verified=False,
        lost_process_status_policy='No native exit code or duration inferred from completed artifacts or solver timing text',
        original_supervisor_status=receipt['status'], original_supervisor_failure=receipt['failure'])


def receipt_view(original, recovery):
    """In-memory assessment view only; never replaces the original receipt."""
    view = dict(original)
    view.update(cases=[recovery['native_case']],
        commands=[recovery['native_command'] if row['name'] == 'native-SaturatedR2' else row for row in original['commands']])
    return view


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt', type=Path, required=True)
    parser.add_argument('--completion-check', type=Path, required=True)
    parser.add_argument('--expected-source', default=PRODUCER)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    require = evidence_module().require
    require(not args.output.exists(), 'Recovery output exists; preserve prior observations')
    e = evidence_module()
    result = dict(schema=1, kind=KIND, success=False, created_utc=datetime.now(timezone.utc).isoformat(),
        collector=e.identity(__file__), assessor=e.identity(Path(__file__).with_name('assess-stage-evidence.py')),
        native_calls=0, scientific_reconstructions=0, failures=[])
    try:
        result['recovery'] = collect(args.receipt, args.completion_check, args.expected_source)
        result['success'] = True
    except Exception as error:
        result['failures'].append(dict(type=type(error).__name__, message=str(error)))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(success=result['success'], output=str(args.output), failures=result['failures'])))
    return 0 if result['success'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
