#!/usr/bin/env python3
"""Assess published R2 stage evidence using metadata only; never run a solver.

A wall-cap exit is an observation about capacity. It is neither a nonlinear
failure classification nor evidence that a coupled root does not exist. In
particular, a complete converged target JSONL record remains positive evidence
when a final evolution/postqualification report was not published.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import struct
import sys
import xml.etree.ElementTree as ET


RUNTIME_KEYS = ('source_inputs', 'runtime_artifacts', 'runtime_bindings',
                'runtime_environment', 'python_packages', 'tool_images')
THREAD_KEYS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
               'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')
TARGET_ALPHA = 400.
TIME = .0025
TOLERANCE = 2e-12
MAXIMUM_UPDATES = 12
WALL_CAP = 14400
LAWS = ('Linear', 'SaturatedR2')
MATRIX_NAME = re.compile(r'case\.feasible\.jsonl\.solve([1-9][0-9]*)\.J([0-9]+)\.bin\Z')
NUMBER = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
START = re.compile(r'^EP-AWF ScalarFeasibleR2 fresh Jacobian iteration (\d+) alpha ('+NUMBER+r') residual ('+NUMBER+r')\s*$')
SUMMARY = re.compile(r'^EP-AWF R2 stage test alpha ('+NUMBER+r') converged ([01]) physical range ([01]) history unchanged\s*$')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            result.update(block)
    return result.hexdigest()


def canonical(value):
    result = hashlib.sha256()
    encoder = json.JSONEncoder(sort_keys=True, ensure_ascii=False, separators=(',', ':'), allow_nan=False)
    for block in encoder.iterencode(value):
        result.update(block.encode())
    return result.hexdigest()


def vector_sha(values):
    require(isinstance(values, list) and all(type(v) in (int, float) and math.isfinite(v) for v in values),
            'A retained coefficient/history vector is not finite')
    result = hashlib.sha256()
    for index in range(0, len(values), 8192):
        chunk = values[index:index+8192]
        result.update(struct.pack('<'+str(len(chunk))+'d', *chunk))
    return result.hexdigest()


def load(path):
    with Path(path).open() as stream:
        return json.load(stream, parse_constant=lambda value: (_ for _ in ()).throw(ValueError('Nonfinite JSON: '+value)))


def identity(path):
    path = Path(path).resolve()
    return dict(path=str(path), bytes=path.stat().st_size, sha256=sha(path))


class Bindings:
    """Hash every read input, including unpublished matrices; detect later changes."""
    def __init__(self):
        self.files = {}

    def file(self, path, expected=None):
        item = identity(path)
        if expected is not None:
            require(item['sha256'] == expected, 'Changed bound file: '+item['path'])
        old = self.files.setdefault(item['path'], item)
        require(old == item, 'Input changed during metadata assessment: '+item['path'])
        return dict(item)

    def finish(self):
        for item in self.files.values():
            require(identity(item['path']) == item, 'Input changed during metadata assessment: '+item['path'])


def session_identity(path):
    """The preflight and stage XML may differ only in the stage-test selector."""
    root = ET.parse(path).getroot()
    selectors = [item for item in root.findall('.//SOLVERINFO/I')
                 if item.get('PROPERTY', '').upper() == 'RANSWALLSTAGETEST']
    require(len(selectors) == 1, 'Missing or repeated stage-test XML selector')
    selector = selectors[0].get('VALUE')
    selectors[0].set('VALUE', 'Off')
    # ElementTree canonicalization removes pretty-print whitespace only; no
    # numerical expression or physical input is evaluated or normalized away.
    normalized = ET.canonicalize(ET.tostring(root, encoding='unicode'), strip_text=True)
    return dict(stage_test=selector, protected_xml_sha256=hashlib.sha256(normalized.encode()).hexdigest())


def runtime_identity(before, after, receipt, bindings):
    for pin in (before, after):
        require(pin.get('success') is True and pin.get('clean') is True, 'Unsuccessful/unclean runtime fingerprint')
        require(pin.get('source_commit') == receipt['source_commit'], 'Runtime/source revision differs')
    result = {}
    for key in RUNTIME_KEYS:
        require(isinstance(before.get(key), dict) and before[key] == after.get(key), 'Runtime section changed: '+key)
        result[key] = canonical(before[key])
    require(before.get('source_root') == after.get('source_root') and
            before.get('build_root') == after.get('build_root') and
            before.get('runtime_resolved_paths') == after.get('runtime_resolved_paths'), 'Runtime roots/resolved image paths changed')
    for name, expected in before['source_inputs'].items():
        path = Path(before['source_root'])/name
        require(path.resolve().is_relative_to(Path(before['source_root']).resolve()), 'Source input escapes its root')
        bindings.file(path, expected)
    resolved = before['runtime_resolved_paths']
    require(set(resolved) == set(before['runtime_artifacts']), 'Runtime image paths are incomplete')
    for name, expected in before['runtime_artifacts'].items():
        bindings.file(resolved[name], expected)
    for tool in before['tool_images'].values():
        bindings.file(tool['resolved'], tool['sha256'])
    solver_key = 'build/IncNavierStokesSolver'
    require(solver_key in resolved, 'Producing solver image is absent')
    solver = Path(resolved[solver_key]).resolve()
    require(solver == (Path(before['build_root'])/'IncNavierStokesSolver').resolve(), 'Producing solver path differs')
    return dict(sections=result, source_root=before['source_root'], build_root=before['build_root'],
                source_commit=before['source_commit'], solver=bindings.file(solver),
                source_input_count=len(before['source_inputs']), runtime_image_count=len(resolved),
                package_and_environment_scope='Complete before/after fingerprint sections are equal; no packages or native images are executed')


def read_log(path):
    starts, summaries = [], []
    with Path(path).open(errors='strict') as stream:
        for number, line in enumerate(stream, 1):
            match = START.fullmatch(line.rstrip('\r\n'))
            if match and float(match[2]) == TARGET_ALPHA:
                starts.append(dict(line=number, iteration=int(match[1]), alpha=float(match[2]), residual=float(match[3])))
            match = SUMMARY.fullmatch(line.rstrip('\r\n'))
            if match and float(match[1]) == TARGET_ALPHA:
                summaries.append(dict(line=number, alpha=float(match[1]), converged=match[2] == '1',
                                      physical_range=match[3] == '1', history_unchanged_declaration=True))
    return dict(target_start_lines=starts, target_summary_lines=summaries)


def trace_prefix(path, law, h0, n, mode, bindings):
    """Only newline-terminated records are complete published JSONL records."""
    records, referenced = [], set()
    prefix = hashlib.sha256()
    prefix_bytes = 0
    tail = None
    allowed = [TARGET_ALPHA] if mode == 'frozen-scalar' else [3200., 1600., 800., TARGET_ALPHA]
    last_position = -1
    with Path(path).open('rb') as stream:
        index = 0
        while True:
            line = stream.readline()
            if not line:
                break
            if not line.endswith(b'\n'):
                tail = dict(offset=prefix_bytes, bytes=len(line), sha256=hashlib.sha256(line).hexdigest(),
                            interpretation='Unpublished tail; no completed-solve classification is inferred')
                break
            require(line.strip(), 'Empty complete JSONL record')
            record = json.loads(line, parse_constant=lambda value: (_ for _ in ()).throw(ValueError('Nonfinite JSONL: '+value)))
            require(record.get('schema') == 1 and record.get('policy') == 'ScalarFeasibleR2:v1', 'Unknown native solve record schema/policy')
            require(record.get('inner_stress_model') == law and record.get('time') == TIME, 'Solve record law/time differs')
            alpha = record.get('alpha')
            require(alpha in allowed and allowed.index(alpha) > last_position, 'Repeated/out-of-order/unregistered solve level')
            last_position = allowed.index(alpha)
            require(type(record.get('converged')) is bool and record.get('maximum_iterations') == MAXIMUM_UPDATES and
                    record.get('tolerance') == TOLERANCE, 'Changed solve tolerance/update cap or missing convergence flag')
            require(record.get('history_advanced') is False and record.get('residual_unchanged') is True and
                    record.get('frozen_velocity_pressure_diagnostic') is (mode == 'frozen-scalar'), 'Changed solve/history scope')
            require(len(record['state']) == n and vector_sha(record['history']) == h0, 'Native solve shape/history is not actual preflight h0')
            iterations = record['iterations']
            require(len(iterations) <= MAXIMUM_UPDATES, 'Solve exceeds the fixed Newton update ceiling')
            record_matrices = set()
            for iteration in iterations:
                matrix = Path(iteration['jacobian_binary_file'])
                matrix = matrix if matrix.is_absolute() else Path(path).parent/matrix
                require(matrix.resolve().parent == Path(path).resolve().parent and MATRIX_NAME.fullmatch(matrix.name),
                        'Referenced matrix is not a local declared native Jacobian')
                require(iteration['jacobian_binary_shape'] == [n, n], 'Referenced matrix shape differs')
                require(int(MATRIX_NAME.fullmatch(matrix.name)[1]) == index+1, 'Referenced matrix solve number differs')
                require(matrix.stat().st_size == 8*n*n, 'Native Jacobian byte length differs')
                bindings.file(matrix)
                referenced.add(str(matrix.resolve()))
                record_matrices.add(str(matrix.resolve()))
            entry = {key: record.get(key) for key in ('alpha', 'time', 'converged', 'terminal', 'tolerance',
                'maximum_iterations', 'history_advanced', 'residual_unchanged', 'frozen_velocity_pressure_diagnostic', 'scaled_block_norms')}
            entry.update(index=index, offset=prefix_bytes, bytes=len(line), raw_line_sha256=hashlib.sha256(line).hexdigest(),
                         canonical_record_sha256=canonical(record), state_sha256=vector_sha(record['state']),
                         history_sha256=h0, iteration_count=len(iterations), matrix_files=sorted(record_matrices))
            records.append(entry)
            prefix.update(line)
            prefix_bytes += len(line)
            index += 1
    return dict(complete_records=records, complete_prefix_bytes=prefix_bytes, complete_prefix_sha256=prefix.hexdigest(),
                unpublished_tail=tail), referenced


def final_report(raw, law, h0, fields, mode):
    require(raw.get('inner_stress_model') == law, 'Final report law differs')
    require(vector_sha(raw['initial']) == h0 and canonical(raw['fields']) == fields, 'Final report initial state/native field descriptors differ')
    require(raw.get('stage_test') == ('FrozenScalarR2' if mode == 'frozen-scalar' else 'MatchedR2'), 'Final report stage scope differs')
    target = [level for level in raw.get('stage_test_levels', []) if level.get('alpha') == TARGET_ALPHA]
    require(len(target) <= 1, 'Repeated final-report target level')
    levels = []
    for level in raw.get('stage_test_levels', []):
        if level.get('attempted'):
            require(level.get('time') == TIME and vector_sha(level['history']) == h0 and
                    level.get('native_history_advanced') is False, 'Final level changed target time/history')
        levels.append({key: level.get(key) for key in ('alpha', 'attempted', 'converged', 'true_coupled_residual_converged',
            'physical_range_passed', 'residual_verification_failure', 'qualification_audit_failure', 'failure', 'scaled_block_norms')})
    snapshot = raw.get('stage_test_initial_native_history')
    history_verified = None
    if snapshot is not None and all(key in raw for key in ('native_history', 'native_history_times', 'state', 'accepted_states')):
        history_verified = (snapshot['time'] == raw.get('time') == 0 and snapshot['steps'] == raw.get('steps') == 0 and
            vector_sha(snapshot['state']) == vector_sha(raw['state']) == h0 and
            snapshot['history'] == raw['native_history'] and snapshot['times'] == raw['native_history_times'] and
            len(raw['accepted_states']) == 1 and vector_sha(raw['accepted_states'][0]) == h0)
    target_level = target[0] if target else {}
    audit = target_level.get('audit', {})
    audit_keys = ('energy', 'interface', 'strong_continuity', 'r2_qualified', 'r2_outer_range', 'r2_covariance', 'r2_scalar_domain')
    audit_complete = bool(target_level.get('attempted') and not target_level.get('qualification_audit_failure') and
                          isinstance(audit, dict) and all(key in audit for key in audit_keys))
    return dict(levels=levels, target_level_present=bool(target),
                target_attempted=target_level.get('attempted'),
                target_solver_converged=target_level.get('converged'),
                target_reported_converged=raw.get('stage_test_target_converged'),
                target_true_coupled_residual_converged=target_level.get('true_coupled_residual_converged'),
                target_native_physical_range=raw.get('stage_test_target_physical_range'),
                native_postqualification_report_available=audit_complete,
                native_history_report_available=snapshot is not None, native_history_unchanged_verified=history_verified,
                native_time=raw.get('time'), native_steps=raw.get('steps'), failure=raw.get('failure'))


def assess(receipt_path, expected_source=None):
    receipt_path = Path(receipt_path).resolve()
    run = receipt_path.parent
    bindings = Bindings()
    receipt_file = bindings.file(receipt_path)
    receipt = load(receipt_path)
    require(receipt.get('schema') == 2 and receipt.get('status') in ('completed', 'incomplete', 'failed') and receipt.get('finished_utc'),
            'Require a terminal stages-family receipt; do not assess an active supervisor')
    mode = receipt['mode']
    require(mode in ('matched', 'frozen-scalar'), 'Unknown stage-family mode')
    require(receipt.get('runtime_unchanged') is True and receipt.get('source_clean_start') is True and
            receipt.get('source_clean_end') is True and receipt.get('source_commit_end') == receipt['source_commit'],
            'Unbound/changed native source or runtime')
    if expected_source is not None:
        require(receipt['source_commit'] == expected_source, 'Producing source is not the requested pin')
    require(receipt.get('target_dt') == TIME and receipt.get('target_time') == TIME and
            receipt.get('maximum_Newton_updates_per_level') == MAXIMUM_UPDATES and
            receipt.get('maximum_seconds_per_native_case') == WALL_CAP, 'Changed target/update/wall-cap contract')
    require(receipt.get('ranks') == 1 and receipt.get('maximum_concurrent_workers') == 1 and
            receipt.get('numerical_threads') == {key:'1' for key in THREAD_KEYS}, 'Changed native worker/thread scope')
    allowed = [TARGET_ALPHA] if mode == 'frozen-scalar' else [3200., 1600., 800., TARGET_ALPHA]
    require(receipt['alpha_levels'] == allowed, 'Changed declared level sequence')
    require(receipt['matrix_encoding'] == dict(dtype='IEEE754 binary64', byte_order=sys.byteorder, layout='row-major'),
            'Unsupported/changed native matrix encoding')
    before_path, after_path = run/'runtime-before.json', run/'runtime-after.json'
    bindings.file(before_path); bindings.file(after_path)
    before, after = load(before_path), load(after_path)
    runtime = runtime_identity(before, after, receipt, bindings)
    for name in ('fingerprint-before', 'fingerprint-after', 'prepare'):
        commands = [row for row in receipt['commands'] if row.get('name') == name]
        require(len(commands) == 1 and commands[0].get('exit_code') == 0, 'Missing/failed completed metadata command: '+name)
        bindings.file(commands[0]['log_path'], commands[0]['log_sha256'])
    preflight_path = Path(receipt['preflight_path']).resolve()
    preflight_file = bindings.file(preflight_path, receipt['preflight_sha256'])
    preflight = load(preflight_path)
    require(preflight.get('status') == 'passed' and preflight.get('full_runtime_source_inputs_unchanged') is True and
            preflight.get('source_commit') == receipt['source_commit'] and
            preflight.get('matched_initial_state_identity_all_quadratures') is True, 'Actual matched preflight is not pinned/passed')
    initial = {law:[row for row in preflight['cases'] if row.get('inner_stress') == law and row.get('quadrature') == 'base'] for law in LAWS}
    require(all(len(rows) == 1 for rows in initial.values()), 'Missing/ambiguous law-specific base preflight')
    original_file = bindings.file(receipt['original_path'], receipt['original_sha256'])
    plan_path = run/'cases/stage-plan.json'
    plan_file = bindings.file(plan_path, receipt['stage_plan_sha256'])
    plan = load(plan_path)
    require(plan.get('alpha_levels') == allowed and plan.get('target_alpha') == TARGET_ALPHA and
            plan.get('maximum_Newton_updates_per_level') == MAXIMUM_UPDATES and plan.get('time_and_native_history_advanced') is False,
            'Prepared stage plan differs from fixed level/history contract')
    require(plan['original']['sha256'] == original_file['sha256'] and Path(plan['original']['path']).resolve() == Path(original_file['path']),
            'Prepared plan original XML differs')
    laws = LAWS if mode == 'matched' else ('SaturatedR2',)
    require(len(plan['cases']) == len(laws) and set(row['inner_stress'] for row in plan['cases']) == set(laws), 'Unexpected planned law family')
    command_names = [row.get('name') for row in receipt['commands']]
    require(len(set(command_names)) == len(command_names) and set(command_names).issubset(
        {'fingerprint-before', 'fingerprint-after', 'prepare'} | {'native-'+law for law in laws}),
        'Unexpected/repeated family command; no coupled retry is permitted')
    require(len({row['inner_stress'] for row in receipt['cases']}) == len(receipt['cases']) and
            set(row['inner_stress'] for row in receipt['cases']).issubset(laws), 'Repeated/unregistered native case')
    case_results = []
    for case in plan['cases']:
        law = case['inner_stress']
        xml = Path(case['xml']).resolve()
        require(xml.is_relative_to(run/'cases') and case.get('representation') == receipt['representation'], 'Prepared XML/representation escapes family')
        xml_file = bindings.file(xml, case['xml_sha256'])
        xml_key = session_identity(xml)
        require(xml_key['stage_test'] == ('FrozenScalarR2' if mode == 'frozen-scalar' else 'MatchedR2'), 'Changed native stage selector')
        prior = initial[law][0]
        bindings.file(prior['xml'], prior['xml_sha256'])
        prior_xml_key = session_identity(prior['xml'])
        require(prior_xml_key['stage_test'] == 'Off' and prior_xml_key['protected_xml_sha256'] == xml_key['protected_xml_sha256'],
                'Native stage XML differs from actual law-specific preflight beyond the stage selector')
        bindings.file(prior['report'], prior['report_sha256'])
        raw = load(prior['report'])
        h0 = vector_sha(raw['initial'])
        require(h0 == prior['initial_state_float64_LE_sha256'] and raw.get('time') == 0 and raw.get('steps') == 0 and
                len(raw['accepted_states']) == 1 and vector_sha(raw['accepted_states'][0]) == h0, 'Actual preflight h0 is inconsistent')
        n, fields = len(raw['initial']), canonical(raw['fields'])
        require(n == case['cost']['unknowns'] and n <= 3000, 'Changed native field size/dense guard')
        item = dict(inner_stress=law, xml=xml_file, protected_xml_sha256=xml_key['protected_xml_sha256'],
                    preflight_h0_sha256=h0, native_field_descriptors_sha256=fields, unknowns=n,
                    target_start_observed_in_log=False, complete_target_solve_record=False,
                    target_attempted=None, converged_target_evidence=None,
                    scalar_solve_convergence_evidence=None, full_coupled_convergence_evidence=None,
                    native_postqualification_report_available=False, native_history_report_available=False,
                    native_history_unchanged_verified=None, issues=[], matrices=[])
        rows = [row for row in receipt['cases'] if row['inner_stress'] == law]
        commands = [row for row in receipt['commands'] if row.get('name') == 'native-'+law]
        require(len(commands) <= 1, 'A coupled native case was retried')
        if not rows:
            require(not commands, 'Native command has no corresponding case receipt')
            item.update(evidence_status='no_native_case_record', unresolved_evidence=True,
                        conditional_frozen_scalar_case_eligible=False)
            case_results.append(item)
            continue
        row = rows[0]
        require(Path(row['xml_path']).resolve() == xml and row['xml_sha256'] == xml_file['sha256'], 'Receipt/planned XML differs')
        item['receipt_target_declarations'] = {key:row.get(key) for key in
            ('target_attempted', 'target_converged', 'target_physical_range')}
        if not commands:
            item.update(evidence_status='no_completed_native_command', unresolved_evidence=True,
                        conditional_frozen_scalar_case_eligible=False)
            case_results.append(item)
            continue
        command = commands[0]
        require(command['argv'] == ['taskset', '-c', '0', str(Path(before['build_root'])/'IncNavierStokesSolver'), str(xml)] and
                Path(command['cwd']).resolve() == xml.parent and command.get('timeout_seconds') == WALL_CAP and
                type(command.get('exit_code')) is int and type(command.get('seconds')) in (int, float) and math.isfinite(command['seconds']),
                'Unexpected native argv/cwd/timeout or unfinished command')
        require(row.get('native_exit_code') == command['exit_code'] and row.get('native_seconds') == command['seconds'], 'Native row/command timing differs')
        log = bindings.file(command['log_path'], command['log_sha256'])
        log_evidence = read_log(log['path'])
        item.update(native_command=command, log=log, log_evidence=log_evidence,
                    wall_cap_stop_observed=command['exit_code'] == 124 and command['seconds'] >= WALL_CAP,
                    wall_cap_interpretation='Recorded exit/timing only; no nonlinear-failure or root-absence inference')
        report_path, trace_path = xml.with_name('case.evolution.json'), xml.with_name('case.feasible.jsonl')
        report = None
        if row.get('report_path'):
            require(Path(row['report_path']).resolve() == report_path, 'Noncanonical final report path')
            item['report'] = bindings.file(report_path, row['report_sha256'])
            report = final_report(load(report_path), law, h0, fields, mode)
            item['final_report_evidence'] = report
            for key in ('native_postqualification_report_available', 'native_history_report_available', 'native_history_unchanged_verified'):
                item[key] = report[key]
        else:
            require(not report_path.exists(), 'Final report exists but is not bound by completed family receipt')
        trace = dict(complete_records=[], complete_prefix_bytes=0, unpublished_tail=None)
        referenced = set()
        if row.get('trace_path'):
            require(Path(row['trace_path']).resolve() == trace_path, 'Noncanonical native trace path')
            item['trace'] = bindings.file(trace_path, row['trace_sha256'])
            trace, referenced = trace_prefix(trace_path, law, h0, n, mode, bindings)
        else:
            require(not trace_path.exists(), 'Native trace exists but is not bound by completed family receipt')
        item['trace_evidence'] = trace
        for path in sorted(xml.parent.glob('case.feasible.jsonl.solve*.J*.bin')):
            require(MATRIX_NAME.fullmatch(path.name) and path.stat().st_size == 8*n*n, 'Malformed native matrix artifact')
            matrix = bindings.file(path)
            matrix.update(shape=[n, n], referenced_by_complete_record=str(path.resolve()) in referenced,
                          scope='Stored matrix only; no derivative/linear-solve accuracy is assessed')
            item['matrices'].append(matrix)
        targets = [record for record in trace['complete_records'] if record['alpha'] == TARGET_ALPHA]
        starts, summaries = log_evidence['target_start_lines'], log_evidence['target_summary_lines']
        item['target_start_observed_in_log'] = bool(starts or summaries)
        item['complete_target_solve_record'] = bool(targets)
        attempted = bool(starts or summaries or targets or (report and report['target_attempted'] is True))
        item['target_attempted'] = True if attempted else (False if report and report['target_attempted'] is False else None)
        if mode == 'frozen-scalar':
            # The frozen solve tests only scalar rows. The stage-summary log,
            # receipt target flag and postqualification flags still test the
            # full coupled residual; false there does not contradict a true
            # frozen scalar solve. Never merge these two convergence scopes.
            flags = [record['converged'] for record in targets]
            coupled_flags = [record['converged'] for record in summaries]
            if type(row.get('target_converged')) is bool:
                coupled_flags.append(row['target_converged'])
            if report:
                if type(report['target_solver_converged']) is bool:
                    flags.append(report['target_solver_converged'])
                coupled_flags.extend(value for value in (report['target_reported_converged'],
                    report['target_true_coupled_residual_converged']) if type(value) is bool)
            item['scalar_solve_convergence_evidence'] = True if True in flags else (False if flags else None)
            item['full_coupled_convergence_evidence'] = True if True in coupled_flags else (False if coupled_flags else None)
            if True in coupled_flags and False in coupled_flags:
                item['issues'].append('Conflicting full coupled convergence declarations; positive coupled evidence is retained')
        else:
            flags = [record['converged'] for record in targets]+[record['converged'] for record in summaries]
            if row.get('target_converged') is True:
                # Never discard positive receipt evidence even if another
                # native artifact is incomplete or inconsistent. A false
                # receipt flag alone establishes no completed nonconvergence.
                flags.append(True)
            if report:
                flags.extend(value for value in (report['target_solver_converged'], report['target_reported_converged'],
                    report['target_true_coupled_residual_converged']) if type(value) is bool)
            item['full_coupled_convergence_evidence'] = True if True in flags else (False if flags else None)
        item['converged_target_evidence'] = True if True in flags else (False if flags else None)
        if True in flags and False in flags:
            item['issues'].append('Conflicting scalar solve convergence declarations; positive scalar evidence is retained'
                if mode == 'frozen-scalar' else
                'Conflicting convergence declarations; positive evidence is retained and blocks conditional fallback')
        if report and report['target_attempted'] is False and attempted:
            item['issues'].append('Target-attempt declarations conflict with actual log/solve evidence')
        if trace['unpublished_tail'] is not None:
            item['issues'].append('Incomplete unpublished JSONL tail; target completion may be unavailable')
        if item['native_history_unchanged_verified'] is False:
            item['issues'].append('Available final native history differs from actual initial snapshot')
        if item['converged_target_evidence'] is True:
            status = 'frozen_scalar_convergence_evidence_available' if mode == 'frozen-scalar' else 'convergence_evidence_available'
        elif targets:
            status = 'complete_frozen_scalar_target_record_not_converged' if mode == 'frozen-scalar' else 'complete_target_record_not_converged'
        elif item['target_attempted'] is True:
            status = 'target_attempt_observed_completion_unavailable'
        else:
            status = 'target_attempt_not_established'
        item.update(evidence_status=status, unresolved_evidence=item['converged_target_evidence'] is not True,
                    full_coupled_root_unresolved=item['full_coupled_convergence_evidence'] is not True,
                    conditional_frozen_scalar_case_eligible=(mode == 'matched' and item['target_attempted'] is True and
                        item['converged_target_evidence'] is not True and not item['issues']))
        case_results.append(item)
    require(len({item['preflight_h0_sha256'] for item in case_results}) == 1, 'Matched actual h0 differs across laws')
    eligible = mode == 'matched' and len(case_results) == 2 and all(item['conditional_frozen_scalar_case_eligible'] for item in case_results)
    bindings.finish()
    return dict(schema=1, method='Strict metadata/hash assessment; no native calls, operator evaluation or numerical qualification',
                native_calls=0, operator_reconstructions=0, original_receipt=receipt_file,
                family_status=receipt['status'], mode=mode, representation=receipt['representation'],
                source_commit=receipt['source_commit'], runtime=runtime, preflight=preflight_file,
                original=original_file, stage_plan=plan_file, target=dict(alpha=TARGET_ALPHA, time=TIME, dt=TIME),
                target_convergence_scope=('converged_target_evidence describes the frozen scalar solve; full coupled evidence is separate and scalar convergence is no coupled-root claim' if mode == 'frozen-scalar'
                    else 'Native coupled-solve declarations; independent derivative/accounting qualification remains separate'),
                cases=case_results, conditional_frozen_scalar_allowed=eligible,
                conditional_scope='At most the separately authorized frozen-scalar diagnostic; never a coupled target retry',
                unknowns_policy='Missing completion/history/audits remain unavailable; a wall cap proves no nonlinear failure or root absence',
                bound_files=bindings.files)


def check_frozen_scalar(assessment_path, receipt_path, *, expected_source, preflight, original, representation, runtime_pin):
    """Future recipe guard: reassess exact files, do not trust a hand-edited flag."""
    current = assess(receipt_path, expected_source)
    if assessment_path is not None:
        stored = load(assessment_path)
        require(stored.get('success') is True and stored.get('assessor', {}).get('sha256') == sha(__file__),
                'Assessment did not succeed with this exact metadata assessor')
        require(canonical(stored.get('assessment')) == canonical(current), 'Assessment/original receipt or bound metadata changed')
    require(current['mode'] == 'matched' and current['conditional_frozen_scalar_allowed'],
            'Both original coupled targets must actually have been attempted and remain unresolved; convergence/unknown attempt blocks fallback')
    for key, path in (('preflight', preflight), ('original', original)):
        require(current[key] == identity(path), 'Conditional diagnostic changed the matched '+key)
    require(current['representation'] == representation, 'Conditional diagnostic changed the matched representation')
    pin = load(runtime_pin)
    require(pin.get('success') is True and pin.get('clean') is True and pin.get('source_commit') == expected_source and
            {key:canonical(pin[key]) for key in RUNTIME_KEYS} == current['runtime']['sections'],
            'Conditional diagnostic changed the matched native runtime/source')
    return dict(original_matched_receipt=current['original_receipt'],
                assessment=identity(assessment_path) if assessment_path is not None else None,
                independently_rechecked_metadata_sha256=canonical(current),
                both_original_targets_observed=True, both_targets_remain_unresolved=True,
                coupled_target_retry_permitted=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt', type=Path, required=True)
    parser.add_argument('--expected-source')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists(), 'Assessment output already exists; retain original observations')
    result = dict(schema=1, success=False, created_utc=datetime.now(timezone.utc).isoformat(),
                  assessor=identity(__file__), native_calls=0, failures=[])
    try:
        result['assessment'] = assess(args.receipt, args.expected_source)
        result['success'] = True
    except Exception as error:
        result['failures'].append(dict(type=type(error).__name__, message=str(error)))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(success=result['success'], output=str(args.output),
                          conditional_frozen_scalar_allowed=result.get('assessment', {}).get('conditional_frozen_scalar_allowed'),
                          failures=result['failures'])))
    return 0 if result['success'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
