"""Read-only derivative verification of the retained first R2 auxiliary J0.

No nonlinear record, solved path or target stage is inferred from this check.
Independent scientific calculations use the existing stage-analysis oracles.
All numerical work runs on the allocated Linux CPU with one thread.
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


THREADS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
    'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')
PRODUCING_COMMIT = '9f20b3053c303e8805f5183f668a8ec468d9656b'


def stamp():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def canonical(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
        allow_nan=False).encode()).hexdigest()


def git(source, *arguments):
    return subprocess.check_output(['git', '-c', 'gc.auto=0', '-c', 'pack.threads=1',
        '-c', 'checkout.workers=1', '-C', str(source), *arguments], text=True).strip()


def save(path, value):
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('source', 'preflight', 'preflight-session', 'stage-session',
            'log', 'matrix', 'runtime-pin', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--expected-commit', required=True)
    parser.add_argument('--cpu', type=int, default=48, help='Explicitly allocated single Linux CPU (default:48)')
    args = parser.parse_args()
    if platform.system() != 'Linux':
        raise SystemExit('Retained J0 scientific verification belongs to Linux')
    args.source = args.source.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    helper = args.source/'solvers/IncNavierStokesSolver/Tests/RANS/cls_epawf_m8_r2_stage_analysis.py'
    files = [args.preflight, args.preflight_session, args.stage_session,
        args.matrix, args.runtime_pin, helper, Path(__file__)]
    receipt = dict(schema=1, success=False, started_utc=stamp(), native_calls=0,
        scope='Initial state, first auxiliary alpha3200 J0 derivative verification only',
        command_argv=sys.argv, driver_sha256=sha(__file__),
        source_commit=git(args.source, 'rev-parse', 'HEAD'),
        source_clean_start=not bool(git(args.source, 'status', '--porcelain')),
        cpu_affinity=sorted(os.sched_getaffinity(0)),
        selected_cpu=args.cpu,
        numerical_threads={key:os.environ.get(key) for key in THREADS},
        inputs={str(path.resolve()):dict(bytes=path.stat().st_size, sha256=sha(path)) for path in files},
        solved_path_assessed=False, target_stage_status_assessed=False,
        physical_trajectory_qualification=False, M8_closed=False, failures=[])
    save(args.output/'receipt.json', receipt)
    started = time.monotonic()
    result = None
    try:
        if receipt['source_commit'] != args.expected_commit or not receipt['source_clean_start']:
            raise ValueError('Analysis source differs from its requested clean pin')
        if args.cpu < 0 or receipt['cpu_affinity'] != [args.cpu] or any(receipt['numerical_threads'][key] != '1' for key in THREADS):
            raise ValueError('Initial-J0 verification requires the explicitly allocated CPU and six numerical thread settings equal to one')
        log_bytes = args.log.read_bytes()
        expression = rb'EP-AWF ScalarFeasibleR2 fresh Jacobian iteration 0 alpha 3200(?:\.0*)? residual ([0-9eE+.-]+)[^\n]*\n'
        matches = list(re.finditer(expression, log_bytes))
        if len(matches) != 1:
            raise ValueError('Actual stage log does not uniquely identify first iteration0 alpha3200')
        match = matches[0]
        prefix = log_bytes[:match.end()]
        snapshot = args.output/'initial-native-log-prefix.bin'
        snapshot.write_bytes(prefix)
        receipt['active_log_prefix'] = dict(path=str(args.log.resolve()), bytes=len(prefix),
            sha256=sha(snapshot), snapshot=str(snapshot), original_log_bytes=len(log_bytes),
            initial_printed_residual=match.group(1).decode(), later_append_permitted=True)
        receipt['inputs'][str(snapshot.resolve())] = dict(bytes=len(prefix), sha256=sha(snapshot))
        save(args.output/'receipt.json', receipt)
        pin = json.loads(args.runtime_pin.read_text())
        if not pin['success'] or not pin['clean'] or pin['source_commit'] != PRODUCING_COMMIT:
            raise ValueError('Producing runtime is not the pinned clean 9f20 epoch')
        producing = Path(pin['source_root'])
        if git(producing, 'rev-parse', 'HEAD') != PRODUCING_COMMIT or git(producing, 'status', '--porcelain'):
            raise ValueError('Producing source changed after the native run pin')
        solver_keys = [key for key in pin['runtime_artifacts'] if Path(key).name == 'IncNavierStokesSolver']
        if len(solver_keys) != 1:
            raise ValueError('No unique producing solver image')
        solver = Path(pin['runtime_resolved_paths'][solver_keys[0]])
        if sha(solver) != pin['runtime_artifacts'][solver_keys[0]]:
            raise ValueError('Actual producing solver changed')
        receipt['producing_solver'] = dict(path=str(solver), sha256=sha(solver), bytes=solver.stat().st_size)
        receipt['producing_source_commit'] = PRODUCING_COMMIT
        receipt['producing_runtime_sections'] = {key:canonical(pin[key]) for key in
            ('source_inputs', 'runtime_artifacts', 'runtime_bindings', 'runtime_environment', 'python_packages', 'tool_images')}
        producer_files = {}
        for name in ('RANSWallEvolutionFeasible.cpp', 'RANSEquilibriumWallEvolution.cpp'):
            relative = 'solvers/IncNavierStokesSolver/EquationSystems/'+name
            path = producing/relative
            if sha(path) != pin['source_inputs'][relative]:
                raise ValueError('Producing initialization/Jacobian source differs from runtime pin')
            producer_files[relative] = sha(path)
        receipt['producing_state_and_matrix_source_files'] = producer_files
        sys.dont_write_bytecode = True
        sys.path.insert(0, str(helper.parent))
        import numpy as np
        import cls_epawf_m8_r2_stage_analysis as analysis
        from cls_epawf_m8_r2_cases import normalized_physics, parameters, solver_info
        from validate_cls_epawf_m8_r2_pressure import PressureOracle
        imported = {str(Path(module.__file__).resolve()):sha(module.__file__)
            for module in list(sys.modules.values()) if getattr(module, '__file__', None)
            and str(Path(module.__file__).resolve()).startswith(str(args.source))}
        receipt['imported_source_files_start'] = imported
        before, stage = (ET.parse(path).getroot() for path in (args.preflight_session, args.stage_session))
        if normalized_physics(before) != normalized_physics(stage):
            raise ValueError('Stage changed protected preflight physics')
        for tag in ('GEOMETRY', 'EXPANSIONS'):
            if analysis.xml_structure(before.find(tag)) != analysis.xml_structure(stage.find(tag)):
                raise ValueError('Stage changed source-bound native '+tag)
        info, params = solver_info(stage), parameters(stage)
        raw = json.loads(args.preflight.read_text())
        if not (raw['dimension'] == 2 and raw['ranks'] == 1 and not raw['mapped'] and
                raw['eta'] == 0 and raw['steps'] == 0 and raw['time'] == 0 and raw['dt'] == .0025):
            raise ValueError('Preflight is not the frozen serial flat initialization')
        if info.get('RANSWALLSTAGETEST') != 'MatchedR2' or info.get('RANSWALLNONLINEARPOLICY') != 'ScalarFeasibleR2':
            raise ValueError('J0 case did not select the actual coupled matched-stage prototype')
        if info['RANSWALLINNERSTRESS'] != raw['inner_stress_model'] or int(params['NUMSTEPS']) != 0 or float(params['TIMESTEP']) != .0025:
            raise ValueError('Selected law or physical-step controls differ')
        if int(params['RANSWALLMAXITERATIONS']) != 12 or float(params.get('RANSWALLNONLINEARTOLERANCE', '2e-12')) != 2e-12:
            raise ValueError('Original residual tolerance or preregistered update cap changed')
        if raw['quadrature_orders'] != list(map(int, (params['RANSWALLTANGENTIALORDER'], params['RANSWALLNORMALORDER']))):
            raise ValueError('Initial native quadrature differs')
        if args.matrix.resolve().parent != args.stage_session.resolve().parent or not args.matrix.name.endswith('.solve1.J0.bin'):
            raise ValueError('Matrix is not the first retained J0 beside the actual stage XML')
        analyzer = analysis.Analyzer(raw, args.stage_session, args.source, args.matrix.parent)
        if analyzer.size != 2165:
            raise ValueError('Initial-J0 check is restricted to the frozen outer6 2165-unknown space')
        h0 = np.asarray(raw['initial'], float)
        if not np.array_equal(h0, np.asarray(raw['state'])):
            raise ValueError('Preflight stored state differs from its actual initial coefficients')
        domain = analyzer.domain(h0, raw['initial_scalar_domain'], retain=True)
        if not domain['admissible']:
            raise ValueError('Initial state is outside the independently checked common scalar domain')
        constraints = analyzer.constraints(h0, raw['reference_constraints'])
        matrix = analyzer.matrix(dict(jacobian_binary_file=str(args.matrix.resolve()),
            jacobian_binary_shape=[2165, 2165]))
        result = dict(schema=1, scope=receipt['scope'], initial_h0_float64_LE_sha256=analysis.vector_hash(h0),
            native_field_descriptors_sha256=canonical(raw['fields']),
            protected_physics_sha256=canonical(normalized_physics(stage)),
            initial_h0_binding='Source-bound identical native geometry/expansions/physics and producing initialization code; completed trace history comparison is not yet available',
            completed_native_nonlinear_record_used=False, source_initial_iteration=0, auxiliary_alpha=3200.,
            physical_stage_time=.0025, original_history_is_initial_h0=True,
            scalar_domain=domain, independent_static_coverage=analyzer.static_coverage,
            continuity_values=analyzer.affine_continuity_values(raw['fixed_reference'], raw['raw_reference_loads']),
            continuity_columns=analyzer.affine_continuity_columns(matrix, False),
            outer_k_action=analyzer.attempt('initial-J0 independent outer-k action', lambda:
                analyzer.independent_outer_k_action(matrix, h0, h0, 3200., constraints, False, domain)),
            pressure_actions=[], solved_path_assessed=False, target_stage_status_assessed=False,
            final_native_history_slots_assessed=False, physical_trajectory_qualification=False, M8_closed=False)
        view = dict(raw, state=h0.tolist(), constraints=[[int(row), float(value)] for row, value in sorted(constraints.items())])
        oracle = PressureOracle(view, analyzer.external_lift(h0))
        for kind in ('inner_pressure', 'outer_pressure', 'gauge', 'pressure_and_gauge'):
            def pressure_check(kind=kind):
                direction = np.zeros(analyzer.size)
                for block in range(analyzer.d+2, len(analyzer.fields)):
                    if kind not in ('pressure_and_gauge', 'inner_pressure' if block == analyzer.d+2 else 'outer_pressure'):
                        continue
                    for index in range(analyzer.fields[block].size):
                        row = analyzer.offsets[block]+index
                        direction[row] = analyzer.audit.us**2*np.sin(1.317*(row+1))
                if kind in ('gauge', 'pressure_and_gauge'):
                    direction[-1] = analyzer.audit.us/analyzer.H*np.sin(1.317*analyzer.size)
                action = matrix@direction
                affine = oracle.affine(direction)
                nonlinear, records = oracle.source_and_trace(direction)
                error = analyzer.audit.norm_record(action-affine-nonlinear, constraints)
                return dict(kind=kind, direction=direction.tolist(), absolute_physical_error=error,
                    fixed_1e_minus7_gate_met=bool(error['maximum'] <= analysis.ACTION_TOLERANCE),
                    nonlinear_source_and_trace_norm=analyzer.audit.norm_record(nonlinear, constraints),
                    ray_records=records, exact_affine_weak_continuity=analyzer.affine_continuity_action(direction, action),
                    comparison='Independent exact affine pressure/roughness/gauge plus full analytic sensor/profile/capped-source/projected-trace chain')
            result['pressure_actions'].append(analyzer.attempt('initial-J0 analytic '+kind, pressure_check))
        result['analysis_success'] = not bool(analyzer.failures)
        result['failures'] = analyzer.failures
        result['checks'], result['errors'] = analyzer.audit.checks, analyzer.audit.errors
        result['dense_matrix_files'] = analyzer.matrix_files
        result['continuity_assembly'] = analyzer.continuity_metadata
        result['initial_derivative_gates_passed'] = bool(result['continuity_columns']['fixed_1e_minus7_gate_met'] and
            result['outer_k_action'].get('gate', {}).get('gate_met') and
            all(row.get('fixed_1e_minus7_gate_met', False) for row in result['pressure_actions']))
        result['qualification'] = 'Only initial auxiliary J0 derivatives are assessed; no completed solve, target residual, target root, physical path or stability qualification'
        save(args.output/'validation.json', result)
        for path, item in analyzer.matrix_files.items():
            receipt['inputs'][path] = item
        receipt['success'] = result['analysis_success']
        receipt['initial_derivative_gates_passed'] = result['initial_derivative_gates_passed']
        receipt['success_definition'] = 'Read-only analysis completed without provenance/implementation exceptions; derivative gate outcomes are separate and retained'
        for path, expected in imported.items():
            if sha(path) != expected:
                raise ValueError('Imported analysis source changed: '+path)
    except Exception as error:
        receipt['success'] = False
        receipt['failures'].append(dict(type=type(error).__name__, message=str(error), traceback=traceback.format_exc()))
    for path, expected in receipt['inputs'].items():
        if sha(path) != expected['sha256'] or Path(path).stat().st_size != expected['bytes']:
            receipt['success'] = False
            receipt['failures'].append(dict(type='ProvenanceError', message='Retained input changed: '+path))
    if 'active_log_prefix' in receipt:
        prefix = receipt['active_log_prefix']
        with args.log.open('rb') as stream:
            prefix['sha256_end'] = hashlib.sha256(stream.read(prefix['bytes'])).hexdigest()
        prefix['log_bytes_end'] = args.log.stat().st_size
        if prefix['sha256_end'] != prefix['sha256']:
            receipt['success'] = False
            receipt['failures'].append(dict(type='ProvenanceError', message='Initial native log prefix changed'))
    receipt.update(source_commit_end=git(args.source, 'rev-parse', 'HEAD'),
        source_clean_end=not bool(git(args.source, 'status', '--porcelain')), seconds=time.monotonic()-started,
        finished_utc=stamp())
    if receipt['source_commit_end'] != receipt['source_commit'] or not receipt['source_clean_end']:
        receipt['success'] = False
        receipt['failures'].append(dict(type='ProvenanceError', message='Analysis source changed'))
    save(args.output/'receipt.json', receipt)
    print(json.dumps({key:receipt.get(key) for key in ('success', 'initial_derivative_gates_passed', 'seconds', 'failures')}))
    return 0 if receipt['success'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
