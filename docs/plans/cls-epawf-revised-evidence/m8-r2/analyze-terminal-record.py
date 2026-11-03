#!/usr/bin/env python3
"""Read-only terminal-state audit of exactly one complete native R2 record.

This deliberately omits Newton paths, matrices, directional endpoints and
retries. Each phase publishes an atomic checkpoint so an externally capped
attempt retains its completed observations without becoming a full-path or
native-history qualification. Scientific execution belongs to Linux.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import traceback
import xml.etree.ElementTree as ET

PRODUCER = '9f20b3053c303e8805f5183f668a8ec468d9656b'
THREADS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
           'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')
RUNTIME_KEYS = ('source_inputs', 'runtime_artifacts', 'runtime_bindings',
                'runtime_environment', 'python_packages', 'tool_images')
SCOPE = ('One completed native record: terminal common scalar inventory, constraints, '
         'true coupled residual and finite-set accounting only. No iteration path, '
         'Jacobian/action endpoints, target-search protocol, native-history audit, '
         'physical trajectory, root-absence proof or M8 closure.')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def stamp():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def canonical(value):
    digest = hashlib.sha256()
    for block in json.JSONEncoder(sort_keys=True, separators=(',', ':'), allow_nan=False).iterencode(value):
        digest.update(block.encode())
    return digest.hexdigest()


def file_identity(path):
    path = Path(path).resolve()
    return dict(path=str(path), bytes=path.stat().st_size, sha256=sha(path))


def git(source, *words):
    return subprocess.check_output(['git', '-c', 'gc.auto=0', '-c', 'pack.threads=1',
        '-c', 'checkout.workers=1', '-C', str(source), *words], text=True).strip()


def save(path, value):
    """Same-filesystem replacement leaves either the old or new whole JSON."""
    temporary = path.with_suffix(path.suffix+'.tmp')
    with temporary.open('w') as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)


def prefix_sha(path, count):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        while count:
            block = stream.read(min(count, 1024*1024))
            require(bool(block), 'Published native prefix was truncated')
            digest.update(block)
            count -= len(block)
    return digest.hexdigest()


def snapshot_record(trace, index, output):
    require(index >= 0, 'Negative native record index')
    prefix, count, offset = hashlib.sha256(), 0, 0
    with Path(trace).open('rb') as stream:
        for _ in range(index+1):
            line = stream.readline()
            require(bool(line) and line.endswith(b'\n'), 'Selected native record is not completely published')
            offset = count
            count += len(line)
            prefix.update(line)
    path = output/('native-terminal-record-'+str(index)+'.raw.jsonl')
    temporary = path.with_suffix(path.suffix+'.tmp')
    with temporary.open('wb') as stream:
        stream.write(line)
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)
    record = json.loads(line)
    require(isinstance(record, dict), 'Native record is not an object')
    # Preserve the exact complete record on disk, but do not inspect or retain
    # its iterations in the scientific analysis. The count is metadata only.
    iteration_count = len(record.get('iterations', []))
    selected = {key:record[key] for key in ('schema', 'policy', 'inner_stress_model',
        'time', 'alpha', 'tolerance', 'maximum_iterations', 'converged', 'terminal',
        'frozen_velocity_pressure_diagnostic', 'history_advanced', 'residual_unchanged',
        'domain', 'state', 'history', 'constraints', 'residual', 'scaled_block_norms')}
    return selected, dict(path=str(Path(trace).resolve()), prefix_bytes=count,
        prefix_sha256=prefix.hexdigest(), selected_offset=offset, selected_bytes=len(line),
        selected_sha256=sha(path), snapshot=str(path), original_trace_bytes=Path(trace).stat().st_size,
        later_append_permitted=True, native_iteration_entries_present=iteration_count,
        iteration_entries_analyzed=0, canonical_terminal_fields_sha256=canonical(selected))


def validate_record(raw, record, preflight_session, stage_session, analysis):
    """Same input binding as verify_auxiliary, without invoking its path audit."""
    from cls_epawf_m8_r2_cases import normalized_physics, parameters, solver_info
    before, stage = (ET.parse(path).getroot() for path in (preflight_session, stage_session))
    require(normalized_physics(before) == normalized_physics(stage), 'Terminal stage changed protected preflight physics')
    geometry = {}
    for name in ('GEOMETRY', 'EXPANSIONS'):
        left, right = analysis.xml_structure(before.find(name)), analysis.xml_structure(stage.find(name))
        require(left == right, 'Terminal stage changed source-bound native '+name)
        geometry[name] = canonical(right)
    info, params = solver_info(stage), parameters(stage)
    mode = info.get('RANSWALLSTAGETEST')
    require(mode in ('MatchedR2', 'FrozenScalarR2'), 'Not a restricted native R2 stage session')
    require(record['frozen_velocity_pressure_diagnostic'] is (mode == 'FrozenScalarR2'), 'Record/XML frozen-scalar scope differs')
    require(mode != 'FrozenScalarR2' or record['alpha'] == 400., 'Frozen-scalar record changed the original target alpha')
    require(info.get('RANSWALLNONLINEARPOLICY') == raw['nonlinear_policy'] == 'ScalarFeasibleR2', 'Changed nonlinear policy')
    require(float(params['TIMESTEP']) == .0025 and int(params['NUMSTEPS']) == 0 and
        raw['dt'] == .0025 and raw['steps'] == 0 and raw['time'] == 0 and raw['ranks'] == 1 and
        raw['dimension'] == 2 and not raw['mapped'] and raw['eta'] == 0, 'Changed fixed serial flat startup/forcing scope')
    require(info['RANSWALLINNERSTRESS'] == record['inner_stress_model'] == raw['inner_stress_model'], 'Changed selected inner law')
    require(list(map(int, (params['RANSWALLTANGENTIALORDER'], params['RANSWALLNORMALORDER']))) == raw['quadrature_orders'],
            'Native preflight and selected record source orders differ')
    require(analysis.vector_hash(record['history']) == analysis.vector_hash(raw['initial']), 'Selected record history differs from actual preflight h0')
    require(record['schema'] == 1 and record['policy'] == 'ScalarFeasibleR2:v1' and
        record['alpha'] in analysis.LEVELS and record['time'] == .0025 and
        record['maximum_iterations'] == int(params.get('RANSWALLMAXITERATIONS', '20')) == 12 and
        record['tolerance'] == float(params.get('RANSWALLNONLINEARTOLERANCE', '2e-12')) == 2e-12,
        'Selected terminal record changed schema, stage/time/tolerance or the bounded update cap')
    require(type(record['converged']) is bool and record['history_advanced'] is False and record['residual_unchanged'] is True,
            'Missing convergence observation, changed residual or declared history advancement')
    forbidden = {'BODYFORCEINNER', 'BODYFORCEOUTER', 'SCALARSOURCEINNER', 'SCALARSOURCEOUTER'}
    require(not any(item.get('NAME', '').upper() in forbidden for item in stage.findall('.//FUNCTION')),
            'Terminal accounting is restricted to the fixed unforced sessions')
    return dict(mode=mode, protected_physics_sha256=canonical(normalized_physics(stage)),
        geometry_and_expansions=geometry, preflight_h0_float64_LE_sha256=analysis.vector_hash(raw['initial']),
        native_field_descriptors_sha256=canonical(raw['fields']), alpha=record['alpha'], time=record['time'],
        native_converged_observation=record['converged'], native_terminal_observation=record['terminal'],
        native_frozen_scalar_diagnostic=record['frozen_velocity_pressure_diagnostic'],
        full_native_protocol_verified=False, final_native_history_slots_assessed=False)


def optional_native_audit(path, stage_session, raw, record, analysis):
    """Only a bound, matching native postqualification table is crosschecked."""
    if path is None:
        return None, dict(available=False, reason='No native audit report supplied; independent accounting only')
    require(Path(path).resolve() == Path(stage_session).resolve().with_name('case.evolution.json'),
            'Optional native audit is not the actual sibling evolution report')
    report = json.loads(Path(path).read_text())
    require(report['inner_stress_model'] == raw['inner_stress_model'] and
        canonical(report['fields']) == canonical(raw['fields']) and
        analysis.vector_hash(report['initial']) == analysis.vector_hash(raw['initial']),
        'Optional native report changes actual fields/law/h0')
    levels = [level for level in report.get('stage_test_levels', []) if level.get('attempted') and
        level['alpha'] == record['alpha'] and level['time'] == record['time'] and
        analysis.vector_hash(level['state']) == analysis.vector_hash(record['state']) and
        analysis.vector_hash(level['history']) == analysis.vector_hash(record['history'])]
    require(len(levels) == 1, 'No unique native audit of the exact selected terminal state/history/alpha/time')
    level = levels[0]
    require(not level.get('qualification_audit_failure'), 'Optional native postqualification audit did not complete')
    audit = level.get('audit')
    required = {'energy', 'interface', 'strong_continuity', 'robust_source', 'wall_sources',
        'constitutive_modification', 'r2_outer_range', 'r2_covariance', 'r2_qualified'}
    require(isinstance(audit, dict) and required.issubset(audit), 'Optional native accounting table is incomplete')
    return audit, dict(available=True, report=file_identity(path), native_budget_keys=sorted(audit),
        coefficient_state_and_history_match=True, native_postqualification_completed=True,
        physical_range_flag_is_separate=level.get('physical_range_passed'),
        full_native_history_audit_performed=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('source', 'preflight', 'preflight-session', 'stage-session', 'trace', 'runtime-pin', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--expected-commit', required=True)
    parser.add_argument('--record-index', type=int, default=0)
    parser.add_argument('--cpu', type=int, default=49, help='Explicitly allocated single Linux CPU (default:49)')
    parser.add_argument('--native-report', type=Path, help='Optional actual sibling report with a completed matching postqualification table')
    args = parser.parse_args()
    require(platform.system() == 'Linux', 'Terminal-state scientific analysis belongs to authoritative Linux')
    args.source = args.source.resolve(); args.output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    helper = args.source/'solvers/IncNavierStokesSolver/Tests/RANS/cls_epawf_m8_r2_stage_analysis.py'
    receipt = dict(schema=1, success=False, status='preparing', scope=SCOPE, started_utc=stamp(),
        native_calls=0, fresh_residual_endpoints=0, analyzed_iteration_paths=0, matrix_files={},
        derivative_accuracy_assessed=False, full_native_case_protocol_assessed=False,
        final_native_history_slots_assessed=False, physical_trajectory_qualification=False,
        target_search_protocol_assessed=False, root_absence_proven=False, M8_closed=False,
        source_commit=git(args.source, 'rev-parse', 'HEAD'),
        source_clean_start=not bool(git(args.source, 'status', '--porcelain')),
        numerical_threads={key:os.environ.get(key) for key in THREADS},
        cpu_affinity=sorted(os.sched_getaffinity(0)), selected_cpu=args.cpu,
        declared_external_capacity_seconds=900, capacity_enforced_by_external_supervisor=True,
        command_argv=sys.argv, wrapper_sha256=sha(__file__), record_index=args.record_index,
        inputs={}, phases=[], failures=[])
    results, objects = {}, {}

    def bind(path, expected=None):
        item = file_identity(path)
        if expected is not None:
            require(item['sha256'] == expected, 'Bound input changed: '+item['path'])
        old = receipt['inputs'].setdefault(item['path'], item)
        require(old == item, 'Input changed during terminal audit: '+item['path'])
        return item

    def flush():
        receipt['elapsed_seconds'] = time.monotonic()-started
        save(args.output/'receipt.json', receipt)
        save(args.output/'validation.json', dict(schema=1, scope=SCOPE, status=receipt['status'],
            success=receipt['success'], analysis_success=receipt['success'], phases=results,
            failures=receipt['failures'], native_calls=0, iteration_path_assessed=False,
            current_action_or_derivative_gates_assessed=False, full_native_case_protocol_assessed=False,
            final_native_history_slots_assessed=False, physical_trajectory_qualification=False, M8_closed=False))

    def phase(name, function, fatal=False):
        index = len(receipt['phases'])+1
        path = args.output/('phase-'+str(index).zfill(2)+'-'+name+'.json')
        row = dict(name=name, status='running', started_utc=stamp(), started_elapsed_seconds=time.monotonic()-started,
            checkpoint=str(path))
        receipt['phases'].append(row); receipt['status']='running'; receipt['current_phase']=name
        save(path, dict(schema=1, scope=SCOPE, phase=row, result=None))
        flush(); phase_started=time.monotonic()
        try:
            answer = function()
            row.update(status='completed', finished_utc=stamp(), seconds=time.monotonic()-phase_started)
            save(path, dict(schema=1, scope=SCOPE, phase=dict(row), result=answer))
            results[name] = answer
            row['completed_checkpoint'] = file_identity(path)
            flush()
            print(json.dumps(dict(phase=name,status='completed',seconds=row['seconds'])),flush=True)
            return answer
        except Exception as error:
            failure = dict(phase=name,type=type(error).__name__,message=str(error),traceback=traceback.format_exc())
            receipt['failures'].append(failure)
            row.update(status='failed',finished_utc=stamp(),seconds=time.monotonic()-phase_started)
            save(path,dict(schema=1,scope=SCOPE,phase=dict(row),failure=failure))
            row['completed_checkpoint'] = file_identity(path)
            flush()
            print(json.dumps(dict(phase=name,status='failed',failure=failure['message'])),flush=True)
            if fatal:
                raise
            return None

    def provenance():
        require(receipt['source_commit'] == args.expected_commit and receipt['source_clean_start'], 'Driver is not the requested clean source pin')
        require(args.cpu >= 0 and receipt['cpu_affinity'] == [args.cpu] and
            all(value == '1' for value in receipt['numerical_threads'].values()), 'Wrong allocated CPU/thread settings')
        for path in (args.preflight,args.preflight_session,args.stage_session,args.runtime_pin,helper,Path(__file__)):
            bind(path)
        if args.native_report is not None:
            bind(args.native_report)
        selected,prefix = snapshot_record(args.trace,args.record_index,args.output)
        bind(prefix['snapshot'],prefix['selected_sha256'])
        receipt['active_trace_prefix'] = prefix
        pin = json.loads(args.runtime_pin.read_text())
        require(pin.get('success') is True and pin.get('clean') is True and pin.get('source_commit') == PRODUCER,
                'Producing runtime is not the successful clean original9f20 pin')
        producing = Path(pin['source_root'])
        require(git(producing,'rev-parse','HEAD') == PRODUCER and not git(producing,'status','--porcelain'),
                'Producing native source changed')
        key='build/IncNavierStokesSolver'
        require(key in pin['runtime_artifacts'] and key in pin['runtime_resolved_paths'], 'No producing solver image')
        solver = bind(pin['runtime_resolved_paths'][key],pin['runtime_artifacts'][key])
        for name in ('RANSEquilibriumWallEvolution.cpp','RANSWallEvolutionPatch.cpp','RANSWallEvolutionFeasible.cpp'):
            relative='solvers/IncNavierStokesSolver/EquationSystems/'+name
            bind(producing/relative,pin['source_inputs'][relative])
        receipt.update(producing_source_commit=PRODUCER,producing_source_path=str(producing),producing_solver=solver,
            producing_runtime_sections={key:canonical(pin[key]) for key in RUNTIME_KEYS})
        sys.dont_write_bytecode=True;sys.path.insert(0,str(helper.parent))
        import numpy as np
        import cls_epawf_m8_r2_stage_analysis as analysis
        raw=json.loads(args.preflight.read_text())
        binding=validate_record(raw,selected,args.preflight_session,args.stage_session,analysis)
        native,native_binding=optional_native_audit(args.native_report,args.stage_session,raw,selected,analysis)
        imported={str(Path(module.__file__).resolve()):sha(module.__file__) for module in list(sys.modules.values())
            if getattr(module,'__file__',None) and str(Path(module.__file__).resolve()).startswith(str(args.source))}
        for path,expected in imported.items():bind(path,expected)
        receipt['imported_source_files_start']=imported
        objects.update(np=np,analysis=analysis,raw=raw,record=selected,native=native,
            state=np.asarray(selected['state'],float),history=np.asarray(selected['history'],float))
        return dict(input_binding=binding,record_snapshot=prefix,native_audit=native_binding,
            selected_state_float64_LE_sha256=analysis.vector_hash(selected['state']),
            selected_history_float64_LE_sha256=analysis.vector_hash(selected['history']),
            selected_native_residual_float64_LE_sha256=analysis.vector_hash(selected['residual']),
            producing_solver=solver,producing_source_commit=PRODUCER,
            matrix_files_used=0,iteration_paths_analyzed=0)

    def constructor():
        analysis,raw,record=objects['analysis'],objects['raw'],objects['record']
        analyzer=analysis.Analyzer(dict(raw,stage_test_target_time=record['time']),args.stage_session,args.source)
        require(len(record['state']) == len(record['history']) == analyzer.size and
                len(record['residual']) in (0,analyzer.size),
                'Terminal state/history/residual shape differs from the bound native fields')
        objects['analyzer']=analyzer;objects['effective']=analyzer.external_lift(objects['state'])
        return dict(size=analyzer.size,common_static_inventory=analyzer.static_coverage,
            gauge_multiplier_scale=analyzer.audit.gauge_variable_scale,
            fields_bound_to_actual_preflight=True,no_matrix_loaded=True,
            native_terminal_residual_export_available=bool(record['residual']))

    def common_domain():
        return objects['analyzer'].domain(objects['state'],objects['record']['domain'],retain=True)

    def constraints():
        analyzer,record=objects['analyzer'],objects['record']
        current=analyzer.constraints(objects['state'],record['constraints'])
        objects['constraints']=current
        return dict(actual_native_constraints_crosschecked=True,
            entries=[[int(row),float(value)] for row,value in sorted(current.items())],
            independent_projected_epsilon_H_retained_as_target=True,
            effective_state_substitutes_external_conditions_only=True)

    def residual():
        analyzer,record=objects['analyzer'],objects['record']
        require(len(record['residual'])==analyzer.size and bool(record['scaled_block_norms']),
                'Native record did not export a complete terminal residual/norm table; common-domain observations remain separate')
        answer,unmasked,current=analyzer.weak(objects['state'],objects['history'],record['alpha'],
            record['residual'],record['constraints'])
        norms=analyzer.audit.norms(record['residual'],current)
        analyzer.audit.near(norms,record['scaled_block_norms'],'terminal_actual_true_mass_dual_norms',2e-10)
        maximum=float(max(norms));scalar=float(max(norms[analyzer.d],norms[analyzer.d+1],norms[-1]))
        independent_norms=answer['physical_norms']
        independent_maximum=float(max(independent_norms))
        independent_scalar=float(max(independent_norms[analyzer.d],
            independent_norms[analyzer.d+1],independent_norms[-1]))
        solved=scalar if record['frozen_velocity_pressure_diagnostic'] else maximum
        analyzer.audit.require(not record['converged'] or solved<record['tolerance'],
            'Native convergence flag does not meet the actual solved residual tolerance')
        objects.update(unmasked=unmasked,constraints=current)
        return dict(weak_residual=answer,actual_full_coupled_norms=norms.tolist(),
            actual_full_coupled_maximum=maximum,actual_scalar_block_maximum=scalar,tolerance=record['tolerance'],
            independent_full_coupled_maximum=independent_maximum,
            independent_scalar_block_maximum=independent_scalar,
            independent_full_coupled_residual_tolerance_met=bool(independent_maximum<record['tolerance']),
            independent_scalar_block_tolerance_met=bool(independent_scalar<record['tolerance']),
            independent_terminal_residual_reconstructions=1,
            full_coupled_residual_tolerance_met=bool(maximum<record['tolerance']),
            solved_scalar_or_coupled_tolerance_met=bool(solved<record['tolerance']),
            frozen_scalar_scope=record['frozen_velocity_pressure_diagnostic'],
            native_converged_observation=record['converged'],
            convergence_of_a_record_is_not_full_case_or_physical_qualification=True)

    def direct():
        analyzer=objects['analyzer']
        answer=objects['analysis'].selected_direct(objects['raw'],analyzer.fields,objects['effective'],analyzer.audit)
        objects['direct']=answer
        return dict(independent_source_and_direct_accounting=answer,
            native_source_budget_table_available=objects['native'] is not None)

    def interface():
        analyzer=objects['analyzer']
        answer=objects['analysis'].selected_interface(objects['raw'],analyzer.fields,objects['effective'])
        objects['interface']=answer
        return dict(full_total_traction_and_k_flux=answer.tolist(),
            tangential_traction_RMS=float(math.sqrt(max(0.,answer[3]**2-answer[2]**2))),
            native_traction_budget_table_available=objects['native'] is not None,
            pressure_jump_is_not_a_mechanical_pressure_continuity_condition=True)

    def continuity():
        answer=objects['analyzer'].continuity(objects['effective']);objects['continuity']=answer
        return dict(independent_strong_continuity=answer.tolist(),
            weak_continuity_checked_separately_in_residual_phase=True,
            native_strong_continuity_table_available=objects['native'] is not None)

    def energy():
        analyzer,analysis,record=objects['analyzer'],objects['analysis'],objects['record']
        storage,stored=analysis._mass(objects['raw'],analyzer.fields,objects['effective'],objects['history'],record['alpha'])
        source=objects['direct']
        subtotal=storage+source['stress_work_minus_k_source']+source['resolved_advection']+source['roughness_force_work']
        answer=dict(storage=float(storage),stored_energy=float(stored),
            stress_work_minus_k_source=source['stress_work_minus_k_source'],resolved_advection=source['resolved_advection'],
            roughness_force_work=source['roughness_force_work'],prescribed_body_work=0.,direct_sum=float(subtotal),anchor_work=0.,
            native_energy_budget_table_available=objects['native'] is not None,
            independent_unmasked_load_available='unmasked' in objects,
            unmasked_load_origin='Independent weak oracle; native JSONL exports the constrained residual only')
        if 'unmasked' in objects and 'constraints' in objects:
            work=analysis._work(objects['effective'],objects['unmasked'],analyzer.fields,analyzer.offsets,analyzer.d)
            reaction=analysis._work(objects['effective'],objects['unmasked'],analyzer.fields,analyzer.offsets,analyzer.d,set(objects['constraints']))
            answer.update(actual_unmasked_weak_work=float(work),essential_reaction=float(reaction),
                free_weak_work=float(work-reaction),direct_remainder=float(work-subtotal))
        objects['energy']=answer
        return answer

    def ranges():
        return objects['analyzer'].restricted_ranges(objects['state'],objects['record']['domain'],native_audit=None)

    def native_crosschecks():
        native=objects['native']
        if native is None:
            return dict(performed=False,reason='No optional actual native audit table; all budgets above are independent terminal-state accounting')
        analyzer,direct=objects['analyzer'],objects['direct']
        names=[]
        for key in ('robust_source','wall_sources','constitutive_modification','inner_stress','inner_constitutive_modification'):
            if key in native:
                analyzer.audit.near(native[key],direct[key],'terminal_native_'+key,3e-10);names.append(key)
        for key in ('interface','strong_continuity'):
            independent=objects['interface' if key=='interface' else 'continuity']
            analyzer.audit.near(native[key],independent,'terminal_native_'+key,3e-10);names.append(key)
        en=objects['energy']
        require('actual_unmasked_weak_work' in en,'Native energy crosscheck lacks independently assembled unmasked loads')
        expected=[en[key] for key in ('actual_unmasked_weak_work','storage','stress_work_minus_k_source',
            'resolved_advection','roughness_force_work','prescribed_body_work','direct_sum','direct_remainder',
            'stored_energy','essential_reaction','anchor_work','free_weak_work')]
        analyzer.audit.near(native['energy'],expected,'terminal_native_complete_energy_and_reaction',3e-10);names.append('energy')
        common=results['common_inventory_physical_ranges']
        analyzer.audit.near(native['r2_outer_range'],common['outer_range'],'terminal_native_fixed_outer_range',3e-10)
        analyzer.audit.near(native['r2_covariance'],common['covariance_minima'],'terminal_native_common_covariance',3e-10)
        analyzer.audit.require(bool(native['r2_qualified'][0])==common['qualified_on_common_inventory'],
            'Native and independent common-inventory qualification differ')
        return dict(performed=True,native_budget_keys_crosschecked=names+['r2_outer_range','r2_covariance','r2_qualified'],
            source_table='Optional exact-state sibling native postqualification report',
            derivative_path_or_final_history_qualification=False)

    flush()
    try:
        phase('provenance_and_complete_record',provenance,fatal=True)
        phase('constructor_and_static_inventory',constructor,fatal=True)
        domain=phase('common_scalar_domain',common_domain)
        if domain is not None and domain.get('admissible'):
            current_constraints=phase('current_projected_constraints',constraints)
            if current_constraints is not None:
                phase('true_coupled_residual',residual)
            else:
                receipt['residual_unavailable_reason']='Independent current constraint crosscheck did not complete; it is not retried'
            phase('selected_source_and_direct_accounting',direct)
            phase('matching_full_traction_and_k_flux',interface)
            phase('strong_continuity',continuity)
            if 'direct' in objects:
                phase('physical_mass_and_complete_energy',energy)
            phase('common_inventory_physical_ranges',ranges)
            ready=all(key in objects for key in ('direct','interface','continuity','energy')) and 'common_inventory_physical_ranges' in results
            if objects['native'] is None or ready:
                phase('optional_native_budget_crosschecks',native_crosschecks)
        else:
            receipt['accounting_unavailable_reason']='Independent common domain did not complete/admit the selected terminal state; no residual/source call is promoted'
    except Exception:
        # phase() already preserved the exact exception and completed pieces.
        pass
    try:
        for item in receipt['inputs'].values():
            require(file_identity(item['path'])==item,'Bound input changed: '+item['path'])
        if 'active_trace_prefix' in receipt:
            prefix=receipt['active_trace_prefix']
            prefix.update(prefix_sha256_end=prefix_sha(args.trace,prefix['prefix_bytes']),trace_bytes_end=args.trace.stat().st_size)
            require(prefix['prefix_sha256_end']==prefix['prefix_sha256'],'Selected complete JSONL prefix changed')
        receipt.update(source_commit_end=git(args.source,'rev-parse','HEAD'),
            source_clean_end=not bool(git(args.source,'status','--porcelain')))
        require(receipt['source_commit_end']==receipt['source_commit'] and receipt['source_clean_end'],'Analysis driver source changed')
        if 'producing_source_path' in receipt:
            producing=Path(receipt['producing_source_path'])
            require(git(producing,'rev-parse','HEAD')==PRODUCER and not git(producing,'status','--porcelain'),
                'Producing source changed during terminal-only audit')
    except Exception as error:
        receipt['failures'].append(dict(phase='final_provenance',type=type(error).__name__,message=str(error),traceback=traceback.format_exc()))
    analyzer=objects.get('analyzer')
    if analyzer is not None:
        receipt['checks']=analyzer.audit.checks;receipt['errors']=analyzer.audit.errors
        require(not analyzer.matrix_files and not analyzer.gates,'Terminal-only wrapper unexpectedly used matrices/current actions')
    required=('provenance_and_complete_record','constructor_and_static_inventory','common_scalar_domain',
        'current_projected_constraints','true_coupled_residual','selected_source_and_direct_accounting',
        'matching_full_traction_and_k_flux','strong_continuity','physical_mass_and_complete_energy',
        'common_inventory_physical_ranges','optional_native_budget_crosschecks')
    receipt['success']=not bool(receipt['failures']) and all(name in results for name in required)
    receipt.update(status='completed' if receipt['success'] else 'incomplete_or_failed',finished_utc=stamp(),
        current_phase=None,success_definition='All terminal-only phases completed with unchanged bindings; actual convergence/range/accounting values are separate observations')
    flush()
    print(json.dumps(dict(success=receipt['success'],status=receipt['status'],completed_phases=list(results),
        seconds=receipt['elapsed_seconds'],failures=receipt['failures'])),flush=True)
    return 0 if receipt['success'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
