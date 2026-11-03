#!/usr/bin/env python3
"""One gated four-case outer6 R2 preflight; never evolves a native history.

Root supplies the next producing epoch and passed operator/action receipts.
The script performs no source synchronization, build, retry, stage solve or
representation search. Every child inherits CPU49 and six one-thread settings.
"""
import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import platform
import signal
import struct
import subprocess
import sys
import time


THREADS = {key: '1' for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
    'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')}
FP_FIELDS = ('source_inputs', 'runtime_artifacts', 'runtime_bindings',
    'python_packages', 'tool_images', 'runtime_environment')


def sha(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            result.update(chunk)
    return result.hexdigest()


def canonical(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False,
        separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def vector_sha(values):
    return hashlib.sha256(struct.pack('<'+str(len(values))+'d', *values)).hexdigest()


def stamp():
    return datetime.now(timezone.utc).isoformat()


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def transition_comparison(reference, candidate):
    """Native observations on common physical Surface(...,8) anchors.

    New-state sensors are independently checked by its preflight oracle. This
    compact comparison itself quotes source-bound native data, not a new model.
    """
    old, new = reference['transition_samples'], candidate['transition_samples']
    require(len(old) == len(new) and len(old) > 0, 'Transition anchor count changed')
    rows = []
    names = ('velocity', 'increment', 'normal', 'k', 'epsilon', 'pressure', 'convection')
    for before in old:
        matches = [after for after in new if max(abs(a-b) for a, b in zip(before['chart'], after['chart'])) <= 1e-11]
        require(len(matches) == 1, 'Transition physical anchor is absent or ambiguous')
        after = matches[0]
        sensor_change = {name: max(abs(a-b) for a, b in zip(old_group, new_group))
            for name, old_group, new_group in zip(names, before['sensors'], after['sensors'])}
        rows.append(dict(chart=before['chart'],
            maximum_coordinate_difference=max(abs(a-b) for a, b in zip(before['chart'], after['chart'])),
            old_sensors=before['sensors'], new_sensors=after['sensors'],
            maximum_sensor_change_by_group=sensor_change,
            old_gates=before['gates'], new_gates=after['gates'],
            chi_change=after['gates'][0]-before['gates'][0],
            old_input=before['input'], new_input=after['input']))
    return dict(common_physical_anchor_count=len(rows),
        maximum_absolute_chi_change=max(abs(row['chi_change']) for row in rows),
        maximum_sensor_change_by_group={name: max(row['maximum_sensor_change_by_group'][name] for row in rows) for name in names},
        samples=rows,
        interpretation='Observed t0 representation difference at unchanged physical sensors/transition; no causal conclusion inferred from mesh support alone')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True, help='Linux experiment root containing source/ and build/')
    parser.add_argument('--expected-source', required=True)
    parser.add_argument('--analysis-source', type=Path, required=True,
        help='Clean frozen source for corrected independent oracles; never compiled by this script')
    parser.add_argument('--expected-analysis-source', required=True)
    parser.add_argument('--reuse-preflight', type=Path,
        help='Preserved prior preflight directory; reuse only its verified completed Linear/base initialization')
    parser.add_argument('--build-receipt', type=Path, required=True)
    parser.add_argument('--runtime-pin', type=Path, required=True)
    parser.add_argument('--inner-gate', type=Path, required=True)
    parser.add_argument('--pressure-gate', type=Path, required=True)
    parser.add_argument('--original', type=Path, required=True)
    parser.add_argument('--r1-report', type=Path, required=True, help='Actual resolved R1 zero-step native report, read-only')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(platform.system() == 'Linux', 'Native preflight belongs on authoritative Linux')
    root, out = args.root.resolve(), args.output.resolve()
    source, build, analysis = root/'source', root/'build', args.analysis_source.resolve()
    require(out.is_relative_to(root) and not out.is_relative_to(source) and not out.is_relative_to(analysis),
        'Generated output must be outside both source checkouts')
    out.mkdir(parents=True, exist_ok=False)
    lock = (root/'validation.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_SH)
    os.sched_setaffinity(0, {49})
    env = dict(os.environ, **THREADS, PYTHONDONTWRITEBYTECODE='1',
        TMPDIR=str(root/'tmp'), LD_LIBRARY_PATH=str(build/'geometry-lib'))
    Path(env['TMPDIR']).mkdir(exist_ok=True)
    tests = source/'solvers/IncNavierStokesSolver/Tests/RANS'
    fingerprint = source/'docs/plans/cls-epawf-revised-evidence/m7/fingerprint-runtime.py'
    program, generator = build/'IncNavierStokesSolver', tests/'cls_epawf_m8_r2_cases.py'
    oracle = analysis/'solvers/IncNavierStokesSolver/Tests/RANS/cls_epawf_m8_r2_preflight.py'
    git = lambda *words: subprocess.check_output(['git', '-c', 'gc.auto=0', '-C', str(source), *words], text=True).strip()
    analysis_git = lambda *words: subprocess.check_output(['git', '-c', 'gc.auto=0', '-C', str(analysis), *words], text=True).strip()
    input_paths = dict(original=args.original.resolve(), r1_report=args.r1_report.resolve(),
        build_receipt=args.build_receipt.resolve(), runtime_pin=args.runtime_pin.resolve(),
        inner_gate=args.inner_gate.resolve(), pressure_gate=args.pressure_gate.resolve(),
        generator=generator, oracle=oracle, fingerprint=fingerprint, solver=program)
    if args.reuse_preflight:
        previous = args.reuse_preflight.resolve()
        input_paths.update(previous_receipt=previous/'receipt.json',
            previous_runtime_before=previous/'runtime-before.json', previous_runtime_after=previous/'runtime-after.json')
    pins = {name: dict(path=str(path), sha256=sha(path)) for name, path in input_paths.items()}
    receipt = dict(schema=2, milestone='M8-R2', role='E outer6 representation preflight',
        status='running', started_utc=stamp(), source_commit=args.expected_source,
        analysis_source=dict(path=str(analysis), commit=args.expected_analysis_source),
        runner=dict(path=str(Path(__file__).resolve()), sha256=sha(__file__)), inputs=pins,
        maximum_workers=1, ranks=1, affinity=sorted(os.sched_getaffinity(0)), numerical_threads=THREADS,
        environment={key: value for key, value in env.items() if key in (*THREADS,
            'LD_LIBRARY_PATH', 'PYTHONPATH', 'PYTHONHOME', 'PYTHONDONTWRITEBYTECODE', 'TMPDIR',
            'LD_PRELOAD', 'LD_AUDIT', 'LANG', 'LC_ALL') or key.startswith(('OMPI_MCA_', 'PMIX_MCA_', 'UCX_', 'I_MPI_', 'MPICH_'))},
        planned_zero_step_initializations=4, actual_zero_step_initializations=0, reused_zero_step_initializations=0,
        native_evolution_calls=0, native_stage_calls=0, extra_Jacobian_probes=0,
        representation='outer6 only', native_GLL_points=[8, 8], source_quadratures=[[8, 16], [12, 24]],
        commands=[], cases=[])

    def save():
        temporary = out/'receipt.json.tmp'
        temporary.write_text(json.dumps(receipt, indent=2, allow_nan=False)+'\n')
        temporary.replace(out/'receipt.json')

    def command(name, argv, cwd, timeout):
        argv = list(map(str, argv))
        log, started = out/(name+'.log'), time.monotonic()
        row = dict(name=name, argv=argv, cwd=str(cwd), timeout_seconds=timeout,
            argv_cwd_environment_sha256=canonical(dict(argv=argv, cwd=str(cwd), environment=receipt['environment'])))
        receipt['commands'].append(row); save()
        with log.open('w') as stream:
            process = subprocess.Popen(argv, cwd=cwd, env=env, stdout=stream,
                stderr=subprocess.STDOUT, start_new_session=True)
            row['pid'] = process.pid; save()
            last = started
            while process.poll() is None:
                now = time.monotonic()
                if now-started > timeout:
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(10)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL); process.wait()
                    row['exit_code'] = 124
                    break
                if now-last >= 30:
                    print(name+': alive '+str(round(now-started, 1))+' s; '+str(log.stat().st_size)+' log bytes', flush=True)
                    last = now
                time.sleep(1)
            else:
                row['exit_code'] = process.returncode
        row.update(seconds=time.monotonic()-started, log=str(log), log_sha256=sha(log))
        save()
        print(name+': exit '+str(row['exit_code'])+' ('+str(round(row['seconds'], 3))+' s)', flush=True)
        if row['exit_code']:
            print(log.read_text()[-3000:], flush=True)
        return row

    save()
    before_exists = False
    try:
        require(git('rev-parse', 'HEAD') == args.expected_source and not git('status', '--porcelain'), 'Producing source is not the frozen clean revision')
        require(analysis_git('rev-parse', 'HEAD') == args.expected_analysis_source and not analysis_git('status', '--porcelain'),
            'Analysis source is not the frozen clean revision')
        built = json.loads(args.build_receipt.read_text())
        require(built.get('status') == 'passed' and built.get('source_commit') == args.expected_source, 'Producing build receipt is not passed at the pinned source')
        for kind, path in (('inner', args.inner_gate), ('pressure', args.pressure_gate)):
            gate = json.loads(path.read_text())
            require(gate.get('status') == 'passed' and gate.get('suite') == kind and gate.get('source_commit') == args.expected_source,
                'Required focused '+kind+' operator/action gate is not passed at the producing source')
            require(gate.get('runtime_unchanged') is True, 'Focused '+kind+' gate lacks unchanged runtime provenance')
        before = command('fingerprint-before', ['python3', '-B', fingerprint, '--root', root, '--output', out/'runtime-before.json'], out, 240)
        require(before['exit_code'] == 0, 'Preflight runtime/source fingerprint failed')
        before_exists = True
        captured, pinned = json.loads((out/'runtime-before.json').read_text()), json.loads(args.runtime_pin.read_text())
        require(captured.get('success') and pinned.get('success'), 'Missing successful runtime fingerprint')
        require(captured['source_commit'] == pinned['source_commit'] == args.expected_source, 'Runtime pin has another source')
        for name in FP_FIELDS:
            require(captured[name] == pinned[name], 'Producing runtime differs from build pin: '+name)
        receipt['build_pin_full_comparison'] = dict(identical_fields=list(FP_FIELDS), passed=True)
        cpp_inputs = {path: digest for path, digest in captured['source_inputs'].items()
            if Path(path).suffix.lower() in ('.cpp', '.cc', '.c', '.cxx', '.h', '.hpp', '.hxx', '.inl')}
        require(all((analysis/path).is_file() and sha(analysis/path) == digest for path, digest in cpp_inputs.items()),
            'Analysis-source C/C++ inputs differ from the actual producing epoch')
        receipt['analysis_Cxx_equivalence_to_producer'] = dict(passed=True, checked_files=len(cpp_inputs),
            producing_commit=args.expected_source, analysis_commit=args.expected_analysis_source,
            scope='C/C++ only; corrected Python oracles are independently pinned and may differ')
        reusable = None
        if args.reuse_preflight:
            prior = json.loads((previous/'receipt.json').read_text())
            require(prior.get('source_commit') == prior.get('source_commit_end') == args.expected_source,
                'Prior initialization does not identify the producing source')
            for key in ('full_runtime_source_inputs_unchanged', 'source_clean_end', 'input_hashes_unchanged'):
                require(prior.get(key) is True, 'Prior preflight lacks verified unchanged provenance: '+key)
            require(prior['inputs']['solver']['sha256'] == pins['solver']['sha256'] and
                prior['inputs']['generator']['sha256'] == pins['generator']['sha256'] and
                prior['inputs']['original']['sha256'] == pins['original']['sha256'],
                'Prior solver, native generator or original XML differs')
            for filename in ('runtime-before.json', 'runtime-after.json'):
                old_runtime = json.loads((previous/filename).read_text())
                require(old_runtime.get('success') and old_runtime['source_commit'] == args.expected_source,
                    'Prior runtime fingerprint failed or has another source')
                for name in FP_FIELDS:
                    require(old_runtime[name] == captured[name], 'Prior/current native epoch differs: '+name)
            rows = [case for case in prior['cases'] if case['name'] == 'outer6-linear-base']
            require(len(rows) == 1 and rows[0].get('solver_exit_code') == 0,
                'Only a completed successful prior Linear/base initialization may be reused')
            reusable = rows[0]
            require(sha(reusable['report']) == reusable['report_sha256'] and
                sha(reusable['xml']) == reusable['xml_sha256'], 'Retained native report or XML changed')
            input_paths.update(reused_native_report=Path(reusable['report']), reused_native_xml=Path(reusable['xml']))
            for name in ('reused_native_report', 'reused_native_xml'):
                pins[name] = dict(path=str(input_paths[name]), sha256=sha(input_paths[name]))
            receipt['reuse_native_epoch_verification'] = dict(passed=True,
                previous_receipt=pins['previous_receipt'], prior_wrapper_status=prior['status'],
                prior_before_after_runtime_identical_to_current=True,
                solver_image_sha256=pins['solver']['sha256'], report_sha256=reusable['report_sha256'],
                native_report_copied=False, completed_initialization_rerun=False)
            save()
        old = json.loads(args.r1_report.read_text())
        require(old.get('time') == 0 and old.get('steps') == 0 and old.get('audit_times') == [0], 'R1 comparison must be the actual initial zero-step report')
        require(len(old['state']) == 1571 and old['eta'] == 0 and not old['mapped'],
            'R1 comparison is not the improved1571 unanchored flat representation')
        prepared = command('prepare', ['python3', '-B', generator, '--original', args.original,
            '--output', out/'cases', '--representation', 'outer6', '--laws', 'matched', '--family', 'preflight'], out, 60)
        require(prepared['exit_code'] == 0, 'Frozen outer6 session generation failed')
        plan = out/'cases/preflight-plan.json'
        matrix = json.loads(plan.read_text())
        require(len(matrix['cases']) == 4, 'Unexpected preflight case count')
        receipt['case_plan_sha256'] = sha(plan); save()
        reports = {}
        for case in matrix['cases']:
            xml, name = Path(case['xml']), case['name']
            row = dict(name=name, xml=str(xml), xml_sha256=sha(xml), cost=case['cost'],
                inner_stress=case['inner_stress'], quadrature=case['quadrature'], steps=0, status='prepared')
            receipt['cases'].append(row); save()
            require(row['xml_sha256'] == case['xml_sha256'], 'Prepared XML hash changed')
            if reusable is not None and name == reusable['name']:
                require(row['xml_sha256'] == reusable['xml_sha256'],
                    'Freshly generated case XML is not byte-identical to retained successful initialization')
                raw_path = Path(reusable['report'])
                row.update(solver_exit_code=0, initialization_wall_seconds=reusable['initialization_wall_seconds'],
                    native_execution='reused completed native initialization; independent audit rerun', new_solver_calls=0,
                    original_case_xml=str(reusable['xml']), original_native_report_sha256=reusable['report_sha256'],
                    initialization_timing_origin='retained prior actual native execution')
                receipt['reused_zero_step_initializations'] += 1
            else:
                receipt['actual_zero_step_initializations'] += 1; save()
                solved = command('initialize-'+name, [program, xml], xml.parent, 1200)
                raw_path = xml.with_name('case.evolution.json')
                row.update(solver_exit_code=solved['exit_code'], initialization_wall_seconds=solved['seconds'],
                    native_execution='fresh native initialization', new_solver_calls=1,
                    initialization_timing_origin='this wrapper execution')
            if raw_path.exists():
                row.update(report=str(raw_path), report_sha256=sha(raw_path))
            save()
            require(row['solver_exit_code'] == 0, 'Zero-step native initialization failed; no retry: '+name)
            raw = json.loads(raw_path.read_text())
            require(raw['time'] == 0 and raw['steps'] == 0 and raw['audit_times'] == [0] and len(raw['accepted_states']) == 1,
                'Preflight advanced or omitted an accepted history')
            require(not raw.get('failure') and not raw['stage_audits'] and not raw['step_audits'] and not raw.get('stage_test_levels'),
                'Zero-step preflight contains stage/evolution work or a failure')
            require(raw['ranks'] == 1 and raw['dt'] == .0025 and raw['eta'] == 0, 'Changed fixed preflight configuration')
            for key in ('parameters', 'ustar', 'direction', 'chi', 'eta', 'map_parameters', 'transition_controls', 'transition_weights'):
                require(raw[key] == old[key], 'R1/new physical or transition input differs: '+key)
            row.update(native_history_count=len(raw['native_history']), native_history_times=raw['native_history_times'],
                initial_state_float64_LE_sha256=vector_sha(raw['accepted_states'][0]),
                native_history_canonical_sha256=canonical(raw['native_history']),
                native_separate_initial_evaluation_seconds=raw.get('initial_evaluation_seconds'),
                native_separate_timing_available='initial_evaluation_seconds' in raw,
                fresh_Jacobian_cost='not measured in preflight; first actual auxiliary level only')
            target = xml.parent/'preflight.json'
            audited = command('audit-'+name, ['python3', '-B', oracle, '--raw', raw_path, '--session', xml,
                '--plan', plan, '--source', analysis, '--expected-source', args.expected_analysis_source, '--output', target], xml.parent, 1800)
            if target.exists():
                row.update(audit=str(target), audit_sha256=sha(target))
            save()
            require(audited['exit_code'] == 0, 'Independent preflight implementation audit failed; no retry: '+name)
            answer = json.loads(target.read_text())
            require(answer['status'] == 'passed' and answer['input_hashes_unchanged'], 'Independent audit did not finish with unchanged inputs')
            result = answer['result']
            require(result['native_inventory_independently_crosschecked'], 'Actual native scalar inventory was not independently crosschecked')
            row.update(status='implementation_audited', finite_observed_joint_preflight=result['finite_observed_joint_preflight'],
                fixed_outer_box_margins=result['finite_ranges']['fixed_outer_box_margins'],
                sampled_K_over_K0=result['finite_ranges']['sampled_K_over_K0'],
                inner_a_over_c=result['finite_ranges']['inner_a_over_c'],
                stricter_implementation_sampled_K_box_margin=result['stricter_implementation_sampled_K_box_margin'],
                independent_covariance=[result['finite_ranges']['covariance_inner_minimum'], result['finite_ranges']['covariance_outer_minimum']],
                inner_stress_source_summary=result['independent_source_energy']['inner_stress'])
            comparison = transition_comparison(old, raw)
            transition_path = xml.parent/'transition-versus-R1.json'
            transition_path.write_text(json.dumps(dict(schema=1,
                old_report=pins['r1_report'], new_report=dict(path=str(raw_path), sha256=sha(raw_path)),
                comparison=comparison), indent=2, allow_nan=False)+'\n')
            row.update(transition_comparison=str(transition_path), transition_comparison_sha256=sha(transition_path),
                maximum_initial_chi_change_from_R1=comparison['maximum_absolute_chi_change'])
            reports[(case['inner_stress'], case['quadrature'])] = raw
            save()
        identities = []
        for quadrature in ('base', 'high'):
            old_state, new_state = (reports[(law, quadrature)]['accepted_states'][0] for law in ('Linear', 'SaturatedR2'))
            require(len(old_state) == len(new_state) == 2165, 'Unexpected matched initial shape')
            identities.append(dict(quadrature=quadrature, bitwise_equal=old_state == new_state,
                maximum_absolute_coefficient_difference=max(abs(a-b) for a, b in zip(old_state, new_state)),
                linear_sha256=vector_sha(old_state), saturated_sha256=vector_sha(new_state)))
        receipt['cross_law_actual_initial_state_identity'] = identities
        receipt['observed_joint_coverage_all_cases'] = all(row['finite_observed_joint_preflight'] for row in receipt['cases'])
        receipt['matched_initial_state_identity_all_quadratures'] = all(row['bitwise_equal'] for row in identities)
        receipt['status'] = 'passed'; save()
    except Exception as error:
        receipt.update(status='failed', failure=dict(type=type(error).__name__, message=str(error)))
        save()
    finally:
        if before_exists:
            checked = command('fingerprint-after', ['python3', '-B', fingerprint, '--root', root,
                '--output', out/'runtime-after.json', '--compare', out/'runtime-before.json'], out, 240)
            after = json.loads((out/'runtime-after.json').read_text()) if (out/'runtime-after.json').exists() else {}
            receipt['full_runtime_source_inputs_unchanged'] = (checked['exit_code'] == 0 and after.get('success') is True and
                after.get('identical_qualification_inputs_and_runtime') is True)
        else:
            receipt['full_runtime_source_inputs_unchanged'] = False
        receipt['source_clean_end'] = not bool(git('status', '--porcelain'))
        receipt['source_commit_end'] = git('rev-parse', 'HEAD')
        receipt['analysis_source_clean_end'] = not bool(analysis_git('status', '--porcelain'))
        receipt['analysis_source_commit_end'] = analysis_git('rev-parse', 'HEAD')
        receipt['input_hashes_unchanged'] = all(sha(value['path']) == value['sha256'] for value in pins.values())
        if not (receipt['full_runtime_source_inputs_unchanged'] and receipt['source_clean_end'] and
                receipt['source_commit_end'] == args.expected_source and receipt['input_hashes_unchanged'] and
                receipt['analysis_source_clean_end'] and receipt['analysis_source_commit_end'] == args.expected_analysis_source):
            receipt['status'] = 'failed'
        receipt['finished_utc'] = stamp(); save()
        print(json.dumps(dict(status=receipt['status'], failure=receipt.get('failure'),
            zero_step_initializations=receipt['actual_zero_step_initializations'],
            observed_joint_coverage=receipt.get('observed_joint_coverage_all_cases'), output=str(out))), flush=True)
    return 0 if receipt['status'] == 'passed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
