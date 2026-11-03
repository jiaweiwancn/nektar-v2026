"""Transfer the retained initial k-action diagnosis to the actual target J0.

Only stored matrices and previously qualified analytic/mass actions are used.
No native calls, fresh residual endpoints, nonlinear retries or terminal-state
inferences are made. The target initial guess is bound through the pinned
native control flow and the completed preceding failed auxiliary record.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import time
import traceback
import xml.etree.ElementTree as ET

PRODUCER = '9f20b3053c303e8805f5183f668a8ec468d9656b'
INITIAL_SHA = '52491ec3974757cf4b177b0502360171e5715448c68ab1ceb19b8bf564e38977'
DIAGNOSIS_SHA = '2224a4e387a92d0a8a160e817bb9265534a3d7291384faf025d23ddea8b6ecf0'
THREADS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
    'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')
RUNTIME_KEYS = ('source_inputs', 'runtime_artifacts', 'runtime_bindings',
    'runtime_environment', 'python_packages', 'tool_images')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def canonical(value):
    digest = hashlib.sha256()
    for block in json.JSONEncoder(sort_keys=True, separators=(',', ':'), allow_nan=False).iterencode(value):
        digest.update(block.encode())
    return digest.hexdigest()


def identity(path):
    path = Path(path).resolve()
    return dict(bytes=path.stat().st_size, sha256=sha(path))


def git(source, *words):
    return subprocess.check_output(['git', '-c', 'gc.auto=0', '-c', 'pack.threads=1',
        '-c', 'checkout.workers=1', '-C', str(source), *words], text=True).strip()


def save(path, value):
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')
    temporary.replace(path)


def stamp():
    return datetime.now(timezone.utc).isoformat()


def option(argv, name):
    return Path(argv[argv.index('--'+name)+1]).resolve()


def selected_family(receipt):
    """The other law may still run; bind this completed case projection only."""
    rows = [row for row in receipt['cases'] if row['inner_stress'] == 'Linear']
    commands = [row for row in receipt['commands'] if row['name'] == 'native-Linear']
    require(len(rows) == len(commands) == 1, 'Missing/repeated actual Linear native case')
    keys = ('mode', 'source_commit', 'source_clean_start', 'preflight_path', 'preflight_sha256',
        'original_path', 'original_sha256', 'representation', 'target_dt', 'target_time', 'alpha_levels',
        'maximum_Newton_updates_per_level', 'maximum_seconds_per_native_case', 'ranks',
        'numerical_threads', 'maximum_concurrent_workers', 'matrix_encoding', 'stage_plan_sha256')
    return {key:receipt[key] for key in keys} | dict(case=rows[0], command=commands[0])


def control_flow(producing, pin):
    names = ('RANSEquilibriumWallEvolution.cpp', 'RANSWallEvolutionPatch.cpp',
             'RANSWallEvolutionPatch.h', 'RANSWallEvolutionFeasible.cpp')
    files = {}
    for name in names:
        relative = 'solvers/IncNavierStokesSolver/EquationSystems/'+name
        path = producing/relative
        require(sha(path) == pin['source_inputs'][relative], 'Native control-flow source differs from the producing pin: '+name)
        files[str(path)] = identity(path)
    source = (producing/'solvers/IncNavierStokesSolver/EquationSystems/RANSEquilibriumWallEvolution.cpp').read_text()
    begin = source.index('void RANSEquilibriumWallEvolution::RunRestrictedStageTest()')
    end = source.index('void RANSEquilibriumWallEvolution::Checkpoint()', begin)
    body = source[begin:end]
    compact = re.sub(r'\s+', '', body)
    required = ('constVectorhistory=m_state;autoguess=history;boolskipAuxiliary=false;',
        'std::vector<double>{3200.,1600.,800.,400.}', 'if(skipAuxiliary&&alpha!=400.)',
        'autocandidate=guess;boolconverged=false;', 'catch(conststd::exception&error)',
        'skipAuxiliary=true;', 'if(converged)guess=candidate;',
        'm_patch->Solve(guess,history,alpha,targetTime,m_implicitData,false,iterations,m_tolerance,m_maximumIterations)')
    require(all(value in compact for value in required), 'Reviewed failed-auxiliary to target h0 control flow differs')
    require(len(re.findall(r'\bguess\s*=', body)) == 2, 'Unexpected additional guess assignment')
    patch = (producing/'solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.cpp').read_text()
    require(re.search(r'RANSWallEvolutionPatch::Solve\(Vector\s+s,const Vector\s*&history', patch),
            'Coupled solve no longer takes the guess by value')
    return files, dict(function_sha256=hashlib.sha256(body.encode()).hexdigest(),
        initial_guess='guess=history=m_state=h0', failed_auxiliary_guess_assignment=False,
        last_feasible_candidate_is_diagnostic_only_after_failure=True,
        remaining_auxiliary_levels_skipped=True, target_alpha=400., target_time=.0025,
        conclusion='For the observed failed first auxiliary solve, the native call for400 receives initial h0 by value',
        scope='Source-control-flow binding; no directly published target state vector or final history is claimed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('source', 'initial-report', 'initial-receipt', 'diagnosis-report', 'diagnosis-receipt',
                 'family-receipt', 'matrix', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--expected-commit', required=True)
    parser.add_argument('--cpu', type=int, choices=(49, 50), default=49)
    args = parser.parse_args()
    require(platform.system() == 'Linux', 'Scientific matrix-action transfer belongs to authoritative Linux')
    args.source = args.source.resolve(); args.output.mkdir(parents=True, exist_ok=False)
    receipt = dict(schema=1, success=False, native_calls=0, fresh_residual_endpoints=0,
        started_utc=stamp(), command_argv=sys.argv, driver_sha256=sha(__file__),
        source_commit=git(args.source, 'rev-parse', 'HEAD'),
        source_clean_start=not bool(git(args.source, 'status', '--porcelain')),
        cpu_affinity=sorted(os.sched_getaffinity(0)), numerical_threads={key:os.environ.get(key) for key in THREADS},
        inputs={}, failures=[], solved_path_assessed=False, target_stage_status_assessed=False,
        final_native_history_slots_assessed=False, physical_trajectory_qualification=False, M8_closed=False)

    def bind(path, expected=None):
        path = Path(path).resolve(); item = identity(path)
        if expected is not None:
            require(item['sha256'] == expected, 'Bound input changed: '+str(path))
        old = receipt['inputs'].setdefault(str(path), item)
        require(old['bytes'] == item['bytes'] and old['sha256'] == item['sha256'], 'Input changed during transfer check: '+str(path))

    for path in (args.initial_report, args.initial_receipt, args.diagnosis_report, args.diagnosis_receipt, args.matrix, Path(__file__)):
        bind(path)
    save(args.output/'receipt.json', receipt); started = time.monotonic()
    result = None
    try:
        require(receipt['source_commit'] == args.expected_commit and receipt['source_clean_start'], 'Unpinned read-only driver')
        require(receipt['cpu_affinity'] == [args.cpu] and all(value == '1' for value in receipt['numerical_threads'].values()),
                'Require the allocated CPU49/50 and six one-thread settings')
        require(sha(args.initial_report) == INITIAL_SHA and sha(args.diagnosis_report) == DIAGNOSIS_SHA,
                'This bounded transfer requires the exact retained initial and component reports')
        previous = json.loads(args.initial_receipt.read_text()); component_receipt = json.loads(args.diagnosis_receipt.read_text())
        require(previous['success'] and component_receipt['success'] and previous['producing_source_commit'] == PRODUCER,
                'Prior independent analysis/provenance did not complete')
        for old_receipt in (previous, component_receipt):
            for path, item in old_receipt['inputs'].items():
                bind(path, item['sha256'])
                require(Path(path).stat().st_size == item['bytes'], 'Prior input byte length differs')
            for path, expected in old_receipt['imported_source_files_start'].items():
                bind(path, expected)
        initial = json.loads(args.initial_report.read_text()); component = json.loads(args.diagnosis_report.read_text())
        require(initial['analysis_success'] and component['analytic_controls_met'] and
                initial['initial_derivative_gates_passed'] is False and component['fixed_absolute_gate_met'] is False,
                'The original failed gates or analytic controls changed')
        old_argv = previous['command_argv']
        raw_path, preflight_xml, stage_xml, old_matrix, runtime_path = (option(old_argv, name) for name in
            ('preflight', 'preflight-session', 'stage-session', 'matrix', 'runtime-pin'))
        scientific_source = option(old_argv, 'source')
        require(git(scientific_source, 'rev-parse', 'HEAD') == previous['source_commit'] and
                not git(scientific_source, 'status', '--porcelain'), 'Retained scientific-oracle source changed')
        pin = json.loads(runtime_path.read_text()); producing = Path(pin['source_root'])
        require(pin['success'] and pin['clean'] and pin['source_commit'] == PRODUCER and
                git(producing, 'rev-parse', 'HEAD') == PRODUCER and not git(producing, 'status', '--porcelain'),
                'Producing native source/runtime is not the clean original9f20 epoch')
        solver = previous['producing_solver']; bind(solver['path'], solver['sha256'])
        require(solver['sha256'] == pin['runtime_artifacts']['build/IncNavierStokesSolver'], 'Producing solver/pin differs')
        receipt.update(producing_source_commit=PRODUCER, producing_solver=solver,
            producing_runtime_sections={key:canonical(pin[key]) for key in RUNTIME_KEYS},
            scientific_source=dict(path=str(scientific_source), commit=previous['source_commit']))
        producer_files, flow = control_flow(producing, pin)
        for path, item in producer_files.items(): bind(path, item['sha256'])
        receipt['producing_control_flow_sources'] = producer_files
        family_bytes = args.family_receipt.read_bytes(); family = json.loads(family_bytes)
        selected = selected_family(family); native = selected['command']; case = selected['case']
        require(selected['mode'] == 'matched' and selected['source_commit'] == PRODUCER and
                selected['target_dt'] == selected['target_time'] == .0025 and selected['alpha_levels'] == [3200.,1600.,800.,400.] and
                selected['maximum_Newton_updates_per_level'] == 12 and selected['maximum_seconds_per_native_case'] == 14400,
                'Actual matched-stage contract changed')
        require(native['exit_code'] == case['native_exit_code'] == 124 and native['seconds'] >= 14400 and
                native['seconds'] == case['native_seconds'] and native['timeout_seconds'] == 14400,
                'Selected Linear case has not completed its declared capacity stop')
        require(Path(case['xml_path']).resolve() == stage_xml and
                native['argv'] == ['taskset','-c','0',str(Path(pin['build_root'])/'IncNavierStokesSolver'),str(stage_xml)] and
                Path(native['cwd']).resolve() == stage_xml.parent, 'Actual native case argv/XML differs')
        bind(stage_xml, case['xml_sha256']); bind(native['log_path'], native['log_sha256'])
        log = Path(native['log_path']).read_bytes()
        prefix = rb'EP-AWF ScalarFeasibleR2 fresh Jacobian iteration 0 alpha '
        observed = {}
        for alpha in (3200,400):
            matches = list(re.finditer(prefix+str(alpha).encode()+rb'(?:\.0*)? residual ([0-9eE+.-]+)[^\n]*\n',log))
            require(len(matches) == 1, 'No unique actual iteration0 log for alpha'+str(alpha))
            observed[str(alpha)] = dict(byte_offset=matches[0].start(), printed_residual=matches[0][1].decode())
        require(observed['3200']['byte_offset'] < observed['400']['byte_offset'] and
                observed['3200']['printed_residual'] == observed['400']['printed_residual'], 'Target initial residual log does not match h0 control-flow evidence')
        require(args.matrix.resolve() == stage_xml.with_name('case.feasible.jsonl.solve2.J0.bin') and
                old_matrix == stage_xml.with_name('case.feasible.jsonl.solve1.J0.bin'), 'Wrong target or initial matrix')
        trace = Path(case['trace_path']); bind(trace, case['trace_sha256'])
        require(trace.resolve() == stage_xml.with_name('case.feasible.jsonl'), 'Wrong actual native JSONL path')
        with trace.open('rb') as stream:
            line = stream.readline(); tail = stream.read()
        require(line.endswith(b'\n') and (not tail or not tail.endswith(b'\n')), 'Require one complete preceding auxiliary record; a complete target record must be assessed separately')
        auxiliary = json.loads(line)
        require(auxiliary['alpha'] == 3200. and auxiliary['time'] == .0025 and auxiliary['converged'] is False and
                auxiliary['history_advanced'] is False and auxiliary['residual_unchanged'] is True and
                auxiliary['inner_stress_model'] == 'Linear' and auxiliary['frozen_velocity_pressure_diagnostic'] is False and
                auxiliary['maximum_iterations'] == 12 and auxiliary['tolerance'] == 2e-12,
                'Preceding native record does not establish the failed coupled auxiliary level')
        receipt['selected_family_metadata'] = dict(path=str(args.family_receipt.resolve()),
            original_receipt_sha256_at_start=hashlib.sha256(family_bytes).hexdigest(),
            projection_sha256=canonical(selected), projection=selected,
            other_case_append_permitted=True, selected_Linear_metadata_must_remain_identical=True)
        receipt['preceding_complete_record'] = dict(bytes=len(line), raw_line_sha256=hashlib.sha256(line).hexdigest(),
            canonical_record_sha256=canonical(auxiliary), converged=False, alpha=3200.,
            unpublished_tail_bytes=len(tail), unpublished_tail_sha256=hashlib.sha256(tail).hexdigest())
        save(args.output/'receipt.json', receipt)
        sys.dont_write_bytecode = True
        sys.path.insert(0, str(scientific_source/'solvers/IncNavierStokesSolver/Tests/RANS'))
        import numpy as np
        import cls_epawf_m8_r2_stage_analysis as stage
        from cls_epawf_m8_r2_cases import normalized_physics
        imported = {str(Path(module.__file__).resolve()):sha(module.__file__)
            for module in list(sys.modules.values()) if getattr(module,'__file__',None) and
            str(Path(module.__file__).resolve()).startswith(str(scientific_source))}
        for path, expected in imported.items():
            require(previous['imported_source_files_start'].get(path) == expected, 'Scientific import differs from retained initial analysis: '+path)
            bind(path, expected)
        receipt['imported_source_files_start'] = imported
        raw = json.loads(raw_path.read_text()); h0 = np.asarray(raw['initial'])
        require(stage.vector_hash(h0) == initial['initial_h0_float64_LE_sha256'] and
                stage.vector_hash(auxiliary['history']) == stage.vector_hash(h0) and
                stage.vector_hash(auxiliary['iterations'][0]['state']) == stage.vector_hash(h0),
                'Actual completed auxiliary history/initial iterate is not the retained preflight h0')
        require(Path(auxiliary['iterations'][0]['jacobian_binary_file']).name == old_matrix.name,
                'Auxiliary iteration0 does not refer to the retained initial matrix')
        before_xml, target_xml = (ET.parse(path).getroot() for path in (preflight_xml,stage_xml))
        require(normalized_physics(before_xml) == normalized_physics(target_xml) and
                canonical(normalized_physics(target_xml)) == initial['protected_physics_sha256'], 'Protected original physics changed')
        geometry = {}
        for name in ('GEOMETRY','EXPANSIONS'):
            left, right = stage.xml_structure(before_xml.find(name)), stage.xml_structure(target_xml.find(name))
            require(left == right, 'Target geometry/expansions differ from source-bound h0')
            geometry[name] = canonical(right)
        analyzer = stage.Analyzer(raw,stage_xml,scientific_source,args.matrix.resolve().parent)
        require(analyzer.size == 2165 and raw['inner_stress_model'] == 'Linear', 'Wrong bounded native space/law')
        constraints = analyzer.constraints(h0,raw['reference_constraints'])
        direction = np.asarray(initial['outer_k_action']['direction'])
        mass = np.asarray(component['exact_mass_action']); analytic = np.asarray(component['complete_initial_analytic_action'])
        require(len(direction) == len(mass) == len(analytic) == analyzer.size, 'Retained action shape differs')
        auxiliary_J = analyzer.matrix(dict(jacobian_binary_file=str(old_matrix),jacobian_binary_shape=[2165,2165]))
        target_J = analyzer.matrix(dict(jacobian_binary_file=str(args.matrix),jacobian_binary_shape=[2165,2165]))
        old_action, target_action = auxiliary_J@direction, target_J@direction
        analyzer.audit.near(old_action,initial['outer_k_action']['actual_current_J_action'],'retained_auxiliary_matrix_action',2e-12)
        independent_mass = analyzer.time_mass(direction,np.zeros(analyzer.size),3200.,constraints)
        analyzer.audit.near(independent_mass,mass,'retained_exact_alpha3200_mass_action',2e-12)
        target_analytic = analytic-(7./8)*mass
        mass_shift_error = target_action-old_action+(7./8)*mass
        target_error = target_action-target_analytic
        target_norm = analyzer.audit.norm_record(target_error,constraints)
        shift_norm = analyzer.audit.norm_record(mass_shift_error,constraints)
        continuity = analyzer.affine_continuity_columns(target_J,False)
        pressure = []
        for row in initial['pressure_actions']:
            require(row['fixed_1e_minus7_gate_met'], 'Retained independent analytic pressure gate did not pass')
            p = np.asarray(row['direction']); delta = target_J@p-auxiliary_J@p
            analyzer.audit.near(analyzer.time_mass(p,np.zeros(analyzer.size),3200.,constraints),np.zeros(analyzer.size),
                'pressure_direction_has_no_time_mass',0.)
            change = analyzer.audit.norm_record(delta,constraints)
            bound = change['maximum']+row['absolute_physical_error']['maximum']
            pressure.append(dict(kind=row['kind'],target_minus_auxiliary_action=change,
                retained_independent_analytic_auxiliary_error=row['absolute_physical_error'],
                target_independent_analytic_error_upper_bound_by_triangle=float(bound),
                fixed_1e_minus7_gate_met=bool(bound<=1e-7),
                scope='Same h0/time and zero pressure mass; exact analytic reference is the bound retained pressure chain, not a new FD endpoint'))
        mask_errors = {}
        for name in ('roughness_velocity_load_rows','outer_velocity_load_rows'):
            value = np.zeros(analyzer.size); ids = component['masks'][name]; value[ids] = target_error[ids]
            mask_errors[name] = analyzer.audit.norm_record(value,constraints)
        result = dict(schema=1,native_calls=0,fresh_residual_endpoints=0,
            scope='Initial alpha400 solve2.J0 action only, source-bound to h0; no terminal target state/root, history or physical qualification',
            target_initial_guess_binding=flow, actual_log_iteration0=observed,
            initial_h0_float64_LE_sha256=stage.vector_hash(h0), native_field_descriptors_sha256=canonical(raw['fields']),
            protected_physics_sha256=initial['protected_physics_sha256'],geometry_and_expansions=geometry,
            direction_float64_LE_sha256=stage.vector_hash(direction), direction=direction.tolist(),
            original_initial_derivative_gates_passed=initial['initial_derivative_gates_passed'],
            original_outer_k_FD_gate=initial['outer_k_action']['gate'],
            original_component_absolute_gate_met=component['fixed_absolute_gate_met'],
            original_auxiliary_J_minus_analytic=component['stored_J_minus_analytic'],
            original_failed_gates_preserved=True, auxiliary_alpha=3200.,target_alpha=400.,physical_stage_time=.0025,
            exact_mass_change='alpha400/alpha3200=1/8; analytic400=analytic3200−(7/8)*mass3200',
            target_minus_auxiliary_minus_exact_mass_change=shift_norm,
            mass_change_control_met=bool(shift_norm['maximum']<=1e-10),
            target_J_minus_transferred_independent_analytic=target_norm,
            target_support_errors=mask_errors, fixed_absolute_action_tolerance=1e-7,
            target_initial_k_action_gate_met=bool(target_norm['maximum']<=1e-7),
            exact_affine_continuity_columns=continuity, pressure_action_transfer=pressure,
            target_k_action=target_action.tolist(),target_independent_analytic_action=target_analytic.tolist(),
            no_rigorous_floating_point_error_proof=True,rigorous_derivative_defect_proven=False,
            nonlinear_terminal_cause_assessed=False,target_stage_status_assessed=False,final_native_history_slots_assessed=False,
            coupled_retry_requested=False,physical_trajectory_qualification=False,M8_closed=False,
            checks=analyzer.audit.checks,errors=analyzer.audit.errors,dense_matrix_files=analyzer.matrix_files)
        save(args.output/'validation.json',result)
        receipt.update(success=True,target_initial_k_action_gate_met=result['target_initial_k_action_gate_met'],
            mass_change_control_met=result['mass_change_control_met'])
        for path,item in analyzer.matrix_files.items(): bind(path,item['sha256'])
    except Exception as error:
        receipt['success'] = False
        receipt['failures'].append(dict(type=type(error).__name__,message=str(error),traceback=traceback.format_exc()))
    for path,item in receipt['inputs'].items():
        if identity(path) != item:
            receipt['success'] = False; receipt['failures'].append(dict(type='ProvenanceError',message='Bound input changed: '+path))
    if 'selected_family_metadata' in receipt:
        family_end = args.family_receipt.read_bytes()
        receipt['selected_family_metadata']['original_receipt_sha256_at_end'] = hashlib.sha256(family_end).hexdigest()
        if canonical(selected_family(json.loads(family_end))) != receipt['selected_family_metadata']['projection_sha256']:
            receipt['success'] = False; receipt['failures'].append(dict(type='ProvenanceError',message='Selected completed Linear metadata changed'))
    receipt.update(source_commit_end=git(args.source,'rev-parse','HEAD'),
        source_clean_end=not bool(git(args.source,'status','--porcelain')),seconds=time.monotonic()-started,finished_utc=stamp())
    if receipt['source_commit_end'] != receipt['source_commit'] or not receipt['source_clean_end']:
        receipt['success'] = False; receipt['failures'].append(dict(type='ProvenanceError',message='Read-only driver source changed'))
    save(args.output/'receipt.json',receipt)
    print(json.dumps({key:receipt.get(key) for key in ('success','target_initial_k_action_gate_met','mass_change_control_met','seconds','failures')}))
    return 0 if receipt['success'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
