#!/usr/bin/env python3
"""Audit one durable R3 alpha=400 snapshot on authoritative Linux."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
import time
import xml.etree.ElementTree as ET

H0 = 'b8ba7917424057a9f6f16a026d82d53ede43985871537facea170ef6e1841b71'
THREADS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
           'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')
RUNTIME_KEYS = ('source_inputs', 'runtime_artifacts', 'runtime_bindings',
                'runtime_environment', 'python_packages', 'tool_images')
TOLERANCE = 2e-12
VALUE_TOLERANCE = 2e-10
ACCOUNT_TOLERANCE = 3e-10


def require(condition, message):
    if not condition:
        raise ValueError(message)


def stamp():
    return datetime.now(timezone.utc).isoformat()


def decode(data):
    def nonfinite(value):
        raise ValueError('Nonfinite JSON constant '+value)
    return json.loads(data, parse_int=lambda value: -0.0 if value == '-0' else int(value),
                      parse_constant=nonfinite)


def load(path):
    return decode(Path(path).read_text())


def canonical(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def identity(path, limit=None):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    count = 0
    with path.open('rb') as stream:
        while limit is None or count < limit:
            block = stream.read(1024*1024 if limit is None else min(1024*1024, limit-count))
            if not block:
                break
            digest.update(block)
            count += len(block)
    require(limit is None or count == limit, 'Published progress prefix was truncated')
    return dict(path=str(path), bytes=count, sha256=digest.hexdigest())


def git(source, *words):
    return subprocess.check_output(['git', '-C', str(source), *words], text=True).strip()


def import_helper(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def compact_domain(domain):
    result = {key:value for key,value in domain.items() if key not in
              ('independent_trace', 'independent_branch_context', 'profile_branch_pattern')}
    if 'independent_trace' in domain:
        result['independent_trace'] = {key:value for key,value in domain['independent_trace'].items()
                                       if key not in ('samples', 'projected_target_coefficients')}
    if 'profile_branch_pattern' in domain:
        pattern = domain['profile_branch_pattern']
        result['profile_branch_populations'] = dict(Counter(pattern))
        result['profile_branch_sha256'] = hashlib.sha256(pattern.encode()).hexdigest()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('run', 'source', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--expected-analysis-source', required=True)
    parser.add_argument('--expected-producer-source', required=True)
    parser.add_argument('--cpu', type=int, required=True)
    parser.add_argument('--r9-local', action='store_true')
    parser.add_argument('--r10-local', action='store_true')
    args = parser.parse_args()
    local_comparison = args.r9_local or args.r10_local
    run, source, output = (path.resolve() for path in (args.run, args.source, args.output))
    require(platform.system() == 'Linux' and os.sched_getaffinity(0) == {args.cpu},
            'Use the authoritative Linux host and one explicitly allocated CPU')
    require(all(os.environ.get(key) == '1' for key in THREADS), 'Set all six numerical thread settings to one')
    require(not output.exists() and not output.is_relative_to(source), 'Use a fresh external audit output')
    sys.dont_write_bytecode = True
    started = time.monotonic()
    result = dict(schema=1, success=False, status='running', started_utc=stamp(), command_argv=sys.argv,
                  analysis_source=args.expected_analysis_source, producing_source=args.expected_producer_source,
                  original_h0_sha256=H0, alpha=400., time=.0025, tolerance=TOLERANCE, cpu=args.cpu,
                  numerical_threads={key:'1' for key in THREADS}, phases={}, failures=[], inputs={},
                  native_calls=0, native_history_advanced=False, physical_trajectory_qualification=False,
                  derivative_accuracy_assessed=False, root_absence_proven=False, M8_closed=False)
    output.parent.mkdir(parents=True, exist_ok=True)

    def save():
        result['elapsed_seconds'] = time.monotonic()-started
        with tempfile.NamedTemporaryFile(mode='w', dir=output.parent, prefix=output.name+'.', delete=False) as stream:
            json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write('\n')
            stream.flush()
            os.fsync(stream.fileno())
            temporary = Path(stream.name)
        temporary.replace(output)

    def bind(path, expected=None):
        item = identity(path)
        if expected is not None:
            require(item == expected, 'Changed bound artifact '+item['path'])
        require(result['inputs'].setdefault(item['path'], item) == item, 'Artifact changed during audit')
        return item

    def bound(item):
        return Path(bind(item['path'], item)['path'])

    def phase(name, function):
        result['current_phase'] = name
        save()
        value = function()
        result['phases'][name] = value
        save()
        return value

    save()
    try:
        require(git(source, 'rev-parse', 'HEAD') == args.expected_analysis_source and
                not git(source, 'status', '--porcelain'), 'Wrong or unclean analysis checkout')
        bind(__file__)
        tests = source/'solvers/IncNavierStokesSolver/Tests/RANS'
        sys.path.insert(0, str(tests))
        import numpy as np
        import cls_epawf_m8_r2_stage_analysis as analysis
        from cls_epawf_m8_r2_cases import parameters, solver_info
        helper_path = source/'docs/plans/cls-epawf-revised-evidence/m8-r3/prepare-predictor.py'
        prepare = import_helper(helper_path, 'r3_terminal_physics_binding')
        bind(helper_path)
        for module in list(sys.modules.values()):
            filename = getattr(module, '__file__', None)
            if filename and Path(filename).resolve().is_relative_to(source):
                bind(filename)
        result['analysis_interpreter'] = bind(sys.executable)
        result['numpy_version'] = np.__version__

        def provenance():
            nonlocal progress, prefix
            receipt_path = run/'receipt.json'
            bind(receipt_path)
            receipt = load(receipt_path)
            require(receipt.get('schema') == 1 and receipt.get('phase') == ('r5-local' if local_comparison else 'target') and receipt.get('finished_utc') and
                    receipt.get('status') in ('completed', 'native_stopped', 'failed') and receipt.get('ranks') == 1 and
                    receipt.get('runtime_unchanged') is True and receipt.get('source_clean_start') is True and
                    receipt.get('source_clean_end') is True and
                    receipt['source_commit'] == receipt['source_commit_end'] == args.expected_producer_source and
                    receipt.get('numerical_threads') == result['numerical_threads'], 'No finalized bound R3 target run')
            before_path, after_path = run/'runtime-before.json', run/'runtime-after.json'
            before_pin, after_pin = bind(before_path), bind(after_path)
            before, after = load(before_path), load(after_path)
            for fingerprint in (before, after):
                require(fingerprint.get('success') is True and fingerprint.get('clean') is True and
                        fingerprint.get('source_commit') == args.expected_producer_source, 'Invalid producing fingerprint')
            require(all(before[key] == after[key] for key in (*RUNTIME_KEYS, 'runtime_resolved_paths', 'source_root', 'build_root')) and
                    after.get('compared_fingerprint') == str(before_path) and
                    after.get('compared_fingerprint_sha256') == before_pin['sha256'], 'Producing runtime binding changed')
            executable = bound(receipt['executable'])
            key = 'build/IncNavierStokesSolver'
            require(executable == Path(before['runtime_resolved_paths'][key]).resolve() and
                    receipt['executable']['sha256'] == before['runtime_artifacts'][key], 'Wrong producing solver image')
            producer = Path(before['source_root']).resolve()
            require(git(producer, 'rev-parse', 'HEAD') == args.expected_producer_source and
                    not git(producer, 'status', '--porcelain'), 'Producing checkout changed')
            for name in ('RANSEquilibriumWallEvolution.cpp', 'RANSWallEvolutionPatch.cpp',
                         'RANSWallEvolutionFeasible.cpp', 'RANSWallEvolutionRobust.cpp',
                         'RANSWallEvolutionTransition.cpp'):
                relative = 'solvers/IncNavierStokesSolver/EquationSystems/'+name
                digest = before['source_inputs'][relative]
                require(bind(producer/relative)['sha256'] == bind(source/relative)['sha256'] == digest,
                        'Analysis/producer operator source differs '+relative)
            commands = receipt['commands']
            require([row['name'] for row in commands] == ['fingerprint-before', 'native', 'fingerprint-after'] and
                    commands[0]['exit_code'] == commands[-1]['exit_code'] == 0, 'Unexpected target command sequence')
            for row in commands:
                bound(row['log'])
            native = commands[1]
            require(type(native.get('exit_code')) is int and type(native.get('seconds')) in (int, float), 'Missing completed process observation')
            require(native['timeout_seconds'] == receipt['maximum_seconds'] and native['seconds'] >= 0,
                    'Changed native budget/elapsed observation')
            session, manifest_path = bound(receipt['session']), bound(receipt['manifest'])
            require(native['argv'] == ['taskset', '-c', str(receipt['native_cpu']), str(executable), str(session)] and
                    Path(native['cwd']).resolve() == session.parent, 'Target launch did not use the bound XML/executable')
            manifest = load(manifest_path)
            require(manifest['current_target_xml'] == receipt['session'] and manifest['current_executable'] == receipt['executable'] and
                    manifest.get('producing_source_commit') == args.expected_producer_source and
                    manifest.get('expected_h0_sha256') == H0 and manifest['origin']['converged'] is False and
                    manifest['target'] == dict(alpha=400., time=.0025, eta=0., tolerance=TOLERANCE, history_advanced=False),
                    'Predictor manifest changed its target or h0 binding')
            preflight = bound(manifest['current_preflight'])
            preflight_session = bound(manifest['current_preflight_session'])
            built_runtime = load(bound(manifest['current_runtime']))
            require(all(before[key] == built_runtime[key] for key in RUNTIME_KEYS), 'Preflight runtime differs from actual target runtime')
            require(prepare.protected_xml(session) == prepare.protected_xml(preflight_session) ==
                    manifest['protected_physical_xml_sha256'], 'Protected target physics/geometry changed')
            xml = ET.parse(session).getroot()
            info, params = solver_info(xml), parameters(xml)
            require(info.get('RANSWALLSTAGETEST') == ('DirectionExperimentR5' if local_comparison else 'ExactTargetR3') and
                    (not args.r9_local or info.get('RANSWALLR9EXPERIMENT') == 'On') and
                    (not args.r10_local or info.get('RANSWALLR10EXPERIMENT') == 'On') and
                    info.get('RANSWALLNONLINEARPOLICY') == 'ScalarFeasibleR3' and
                    info.get('RANSWALLACTIONPOLICY') == 'VerifiedR3Action' and
                    info.get('RANSWALLINNERSTRESS') == 'SaturatedR2' and float(params['TIMESTEP']) == .0025 and
                    int(params['NUMSTEPS']) == 0 and float(params.get('RANSWALLNONLINEARTOLERANCE', '2e-12')) == TOLERANCE,
                    'Changed exact-target session scope')
            require(not any(item.get('NAME', '').upper() in
                            ('BODYFORCEINNER', 'BODYFORCEOUTER', 'SCALARSOURCEINNER', 'SCALARSOURCEOUTER')
                            for item in xml.findall('.//FUNCTION')), 'Prescribed forcing is outside the audit scope')
            raw = load(preflight)
            require(raw['time'] == raw['steps'] == raw['eta'] == 0 and raw['dt'] == .0025 and
                    raw['dimension'] == 2 and raw['mapped'] is False and raw['ranks'] == 1 and
                    raw['inner_stress_model'] == 'SaturatedR2' and raw['model_variant'] == 'Regularized' and
                    raw['robust_parameters'] == [.4, 0, .2, .1, 0, 1, 1] and
                    hashlib.sha256(raw['signature'].encode()).hexdigest() == manifest['current_native_signature_sha256'] and
                    raw['initial'] == raw['state'] and analysis.vector_hash(raw['state']) == H0,
                    'Current preflight fields/physics/h0 changed')
            collection_path = bound(receipt['durable_collection']['file'])
            collection = load(collection_path)
            protocol = collection['protocol']
            require(collection['success'] is True and collection['original_h0_sha256'] == H0 and
                    collection['unknowns'] == 2165 and protocol['benchmark_only'] is False and
                    protocol['alpha'] == 400 and protocol['time'] == .0025 and protocol['tolerance'] == TOLERANCE,
                    'Wrong durable collection protocol or h0')
            progress = (session.parent/info.get('RANSWALLR3PROGRESSFILE', 'case.r3-progress.jsonl')).resolve()
            observed = collection['progress']
            prefix = identity(progress, observed['complete_prefix_bytes'])
            require(prefix['path'] == observed['path'] and prefix['sha256'] == observed['complete_prefix_sha256'], 'Changed durable event prefix')
            last = bound(collection['last_complete_snapshot'])
            snapshot_pin = collection['last_snapshot']['snapshot']
            snapshot = load(last)
            expected_snapshot = Path(str(progress)+'.solve1.iterate'+str(snapshot['iteration'])+'.snapshot.json')
            require(last == Path(str(progress)+'.last.json') and bound(snapshot_pin) == expected_snapshot and
                    identity(last)['sha256'] == snapshot_pin['sha256'] and
                    canonical(snapshot) == canonical({key:value for key,value in collection['last_snapshot'].items() if key != 'snapshot'}),
                    'Atomic last snapshot differs from immutable completed record')
            require(snapshot.get('schema') == 1 and snapshot['policy'] == 'ScalarFeasibleR3:v1' and
                    snapshot['action_policy'] == 'VerifiedR3Action:v1' and snapshot['solve_number'] == 1 and
                    snapshot['alpha'] == 400 and snapshot['time'] == .0025 and snapshot['is_complete_internal_iterate'] is True and
                    snapshot['history_advanced'] is False and snapshot['is_physical_time_sample'] is False,
                    'Invalid complete exact-target snapshot scope')
            vectors = {}
            for name in ('state', 'history', 'residual', 'constraints'):
                path = bound(snapshot[name]); expected_size = 2*snapshot['constraint_pairs'] if name == 'constraints' else 2165
                require(path.parent == progress.parent and snapshot[name]['bytes'] == 8*expected_size,
                        'Invalid durable '+name+' payload shape/path')
                vectors[name] = np.fromfile(path, dtype='<f8')
                require(len(vectors[name]) == expected_size and np.all(np.isfinite(vectors[name])), 'Invalid binary payload values')
            require(snapshot['history']['sha256'] == H0 and analysis.vector_hash(vectors['history']) == H0,
                    'Durable history differs bitwise from original h0')
            pairs = vectors['constraints'].reshape(-1, 2)
            rows = [int(row) for row in pairs[:, 0]]
            require(rows == sorted(set(rows)) and all(0 <= row < 2165 and row == value for row,value in zip(rows,pairs[:,0])),
                    'Invalid durable constraint-row encoding')
            terminal = collection['terminal']
            with progress.open('rb') as stream:
                events = [decode(line) for line in stream.read(prefix['bytes']).splitlines()]
            final_events = [event for event in events if event['kind'] == 'solve_terminal']
            require(final_events == ([] if terminal is None else [terminal]), 'Collection terminal differs from published stream')
            published = [event for event in events if event['kind'] == 'iterate_completed']
            atomic_ahead = collection['atomic_last_snapshot_precedes_completed_event']
            require(bool(published) and (published[-1]['iteration'] == snapshot['iteration']-1 if atomic_ahead else
                    published[-1]['iteration'] == snapshot['iteration'] and published[-1]['snapshot'] == snapshot_pin),
                    'Last snapshot differs from completed event sequence')
            converged = None if terminal is None else terminal['converged']
            require((terminal is None or type(converged) is bool) and
                    receipt['durable_collection']['native_converged_observation'] is converged,
                    'Conflicting durable convergence observations')
            objects.update(receipt=receipt, raw=raw, session=session, snapshot=snapshot, vectors=vectors,
                           constraints=pairs.tolist(), converged=converged, producing_checkout=producer)
            return dict(receipt=result['inputs'][str(receipt_path)], current_preflight=str(preflight),
                        runtime_before=before_pin, runtime_after=after_pin,
                        runtime_sections={key:canonical(before[key]) for key in RUNTIME_KEYS},
                        selected_iteration=snapshot['iteration'], last_snapshot=snapshot_pin,
                        payloads={key:snapshot[key] for key in vectors}, event_prefix=prefix,
                        native_terminal_observation=terminal['terminal'] if terminal else None,
                        native_converged_observation=converged, native_exit_code=native['exit_code'],
                        durable_native_scaled_block_norms=snapshot['scaled_block_norms'],
                        supervisor_wall_cap_observed=native['exit_code'] == 124,
                        last_snapshot_precedes_event=atomic_ahead)

        objects = {}; progress = prefix = None
        provenance_result = phase('bindings', provenance)
        raw, snapshot, vectors = (objects[key] for key in ('raw', 'snapshot', 'vectors'))
        state, history, residual = (vectors[key] for key in ('state', 'history', 'residual'))

        def native_report():
            item = objects['receipt'].get('native_report')
            if item is None:
                return dict(available=False, reason='Final native report was not published')
            path = bound(item)
            require(path == objects['session'].with_suffix('.evolution.json'), 'Native report is outside the actual target case')
            try:
                report = load(path)
            except (ValueError, OSError) as error:
                return dict(available=False, partial_file=item, reason=str(error),
                            stopped_process_partial_report=provenance_result['supervisor_wall_cap_observed'])
            try:
                require(report['stage_test'] == ('DirectionExperimentR5' if local_comparison else 'ExactTargetR3') and report['stage_test_native_history_advanced'] is False and
                        report['time'] == report['steps'] == 0 and analysis.vector_hash(report['state']) == H0 and
                        canonical(report['fields']) == canonical(raw['fields']), 'Native final report changes fields/h0/history scope')
                levels = [level for level in report['stage_test_levels'] if level.get('attempted') and
                          level['alpha'] == 400 and level['time'] == .0025 and
                          analysis.vector_hash(level['state']) == snapshot['state']['sha256'] and
                          analysis.vector_hash(level['history']) == H0]
                require(len(levels) == 1 and levels[0]['native_history_advanced'] is False, 'No matching native final-state level')
                start = report['stage_test_initial_native_history']
                accepted = report.get('accepted_states')
                require(start['time'] == start['steps'] == 0 and
                        analysis.vector_hash(start['state']) == analysis.vector_hash(report['state']) == H0 and
                        (accepted is None or len(accepted) == 1 and analysis.vector_hash(accepted[0]) == H0) and
                        [analysis.vector_hash(slot) for slot in start['history']] ==
                        [analysis.vector_hash(slot) for slot in report['native_history']] and
                        start['times'] == report['native_history_times'], 'Native history slots changed')
            except (KeyError, TypeError, ValueError) as error:
                return dict(available=False, file=item, binding_mismatch=True, reason=str(error))
            level = levels[0]
            objects['native_level'] = level
            return dict(available=True, file=item, native_postqualification_failure=level.get('qualification_audit_failure'),
                        native_integrator_history_slots_verified=True,
                        native_accepted_state_table_crosschecked=accepted is not None,
                        residual_matches_durable_bytes=analysis.vector_hash(level['residual']) == snapshot['residual']['sha256'],
                        native_target_converged=report['stage_test_target_converged'],
                        native_target_physical_range=report['stage_test_target_physical_range'])

        report_metadata = phase('native_report', native_report)

        class SnapshotAnalyzer(analysis.Analyzer):
            # R2 independently reconstructs all moving and fixed scalar rows.
            # Capture those rows for a timeout with no native final table.
            range_inventory = None

            def domain_result(self, admitted, phase_name, records, counts, margins, minimum, native):
                self.observed_rows = records
                self.scalar_progress = dict(admissible=bool(admitted), phase=phase_name, counts=dict(counts),
                    invalid_points=[row for row in records if not row['admissible']][:4], minima=minimum)
                return super().domain_result(admitted, phase_name, records, counts, margins, minimum, native)

            def domain(self, state, native=None, retain=False):
                if native is self.range_inventory and native is not None:
                    return self.range_domain
                return super().domain(state, native, retain)

        analyzer = SnapshotAnalyzer(dict(raw, stage_test_target_time=.0025), objects['session'], source)
        level = objects.get('native_level', {})
        native_domain = dict(snapshot['scalar_domain'])
        full_domain = level.get('scalar_domain', {})
        for key in ('observations', 'effective_state', 'trace_target_state'):
            if key in full_domain:
                native_domain[key] = full_domain[key]

        def scalar():
            try:
                domain = analyzer.domain(state, native_domain, retain=True)
                require(domain['admissible'], 'Independent scalar domain rejected the complete durable state')
            except Exception:
                result['independent_scalar_failure_observations'] = getattr(analyzer, 'scalar_progress', None)
                raise
            rows = analyzer.observed_rows
            analyzer.range_inventory = {'observations':native_domain.get('observations') or rows}
            analyzer.range_domain = domain
            return dict(compact_domain(domain), static_coordinate_coverage=analyzer.static_coverage,
                        native_final_coordinate_multiset_crosschecked=bool(native_domain.get('observations')),
                        sign_unresolved_observations=sum(not row['sign_resolved'] for row in rows))

        phase('scalar_and_projected_trace', scalar)
        effective = analyzer.external_lift(state)

        def relative_component_error(left, right):
            left, right = np.asarray(left, float), np.asarray(right, float)
            require(left.shape == right.shape and np.all(np.isfinite(left)) and np.all(np.isfinite(right)),
                    'Invalid comparison vector')
            return float(np.max(np.abs(left-right)/np.maximum(1., np.abs(right)), initial=0.))

        def coupled_residual():
            constraints = analyzer.constraints(state)
            actual = analyzer.audit.constraints(objects['constraints'])
            require(set(actual) == set(constraints), 'Actual and independent constraint rows differ')
            ids = sorted(constraints)
            constraint_error = relative_component_error([actual[i] for i in ids], [constraints[i] for i in ids])
            loads, quadrature = analysis.weak_loads(raw, analyzer.fields, effective)
            unmasked = loads + analyzer.time_mass(effective, history, 400.)
            independent = unmasked.copy()
            for row, value in constraints.items():
                independent[row] = state[row]-value
            objects.update(unmasked=unmasked, independent_constraints=constraints)
            observed = analyzer.audit.norm_record(residual, constraints)
            reconstructed = analyzer.audit.norm_record(independent, constraints)
            difference = analyzer.audit.norm_record(residual-independent, constraints)
            matrix, _ = analyzer.affine_continuity_matrix()
            exact = matrix@effective
            weak_rows = analyzer.continuity_rows()
            affine_errors = {name:analyzer.continuity_error_norms(value[weak_rows]-exact).tolist()
                             for name,value in (('native', residual), ('independent', independent))}
            component_error = relative_component_error(residual, independent)
            norm_error = relative_component_error(observed['blocks'], snapshot['scaled_block_norms'])
            return dict(native=observed, independent=reconstructed, difference=difference, quadrature=quadrature,
                        native_tolerance_met=observed['maximum'] < TOLERANCE,
                        independent_tolerance_met=reconstructed['maximum'] < TOLERANCE,
                        relative_component_value_error=component_error, native_norm_table_error=norm_error,
                        independent_constraint_target_error=constraint_error, constrained_rows=len(ids),
                        complete_value_comparison_met=max(component_error,norm_error,constraint_error) <= VALUE_TOLERANCE,
                        independent_residual_float64_LE_sha256=analysis.vector_hash(independent),
                        weak_continuity_and_gauge_absolute_errors=affine_errors,
                        weak_affine_values_gate_met=max(max(values) for values in affine_errors.values()) <= 1e-7)

        residual_result = phase('coupled_residual', coupled_residual)

        def physical():
            report = analyzer.physical(state, history, 400., unmasked=objects['unmasked'],
                                       constraints=objects['independent_constraints'], native_domain=analyzer.range_inventory)
            common = report['common_inventory_physical_qualification']
            common['coordinate_origin'] = ('matching native final point table' if native_domain.get('observations')
                                           else 'independently reconstructed finite set; durable counts/minima/branches crosschecked')
            traction = report['full_total_traction_and_k_flux']
            report['matching_observations'] = dict(zip(('area', 'reduced_pressure_jump_integral', 'normal_total_traction_jump_rms',
                'full_total_traction_jump_rms', 'max_absolute_normal_traction_jump', 'k_flux_jump_rms',
                'inner_normal_velocity_integral', 'normal_total_traction_jump_integral'), traction))
            report['matching_observations']['tangential_total_traction_jump_rms'] = report['tangential_traction_RMS']
            energy, direct = report['energy'], report['source_energy']
            _, E, D, C, *_ = direct['robust_source']
            source_scale = max(1., abs(E), abs(D), abs(C))
            report['accounting_defects'] = dict(E_minus_D_minus_C=E-D-C, integrated_redistribution=C,
                energy_weak_minus_direct=energy['direct_remainder'])
            report['independent_accounting_identities_met'] = bool(max(abs(E-D-C),abs(C))/source_scale <= ACCOUNT_TOLERANCE and
                abs(energy['direct_remainder'])/max(1.,abs(energy['actual_unmasked_weak_work']),abs(energy['direct_sum'])) <= ACCOUNT_TOLERANCE)
            native = level.get('audit', {}) if not level.get('qualification_audit_failure') else {}
            energy_order = ('actual_unmasked_weak_work', 'storage', 'stress_work_minus_k_source', 'resolved_advection',
                'roughness_force_work', 'prescribed_body_work', 'direct_sum', 'direct_remainder', 'stored_energy',
                'essential_reaction', 'anchor_work', 'free_weak_work')
            expected = {key:direct[key] for key in ('robust_source', 'wall_sources', 'constitutive_modification',
                                                   'inner_stress', 'inner_constitutive_modification')}
            expected.update(interface=traction, strong_continuity=report['strong_continuity'],
                energy=[energy[key] for key in energy_order], r2_outer_range=common['outer_range'],
                r2_covariance=common['covariance_minima'], r2_qualified=[int(common['qualified_on_common_inventory'])])
            errors = {key:relative_component_error(native[key],value) for key,value in expected.items() if key in native}
            report['native_budget_crosschecks'] = dict(available=bool(errors), complete=set(errors)==set(expected),
                errors=errors, tolerance=ACCOUNT_TOLERANCE, passed=all(error <= ACCOUNT_TOLERANCE for error in errors.values()) if errors else None)
            return report

        physical_result = phase('finite_ranges_and_accounting', physical)
        native_checks = physical_result['native_budget_crosschecks']
        finite = physical_result['common_inventory_physical_qualification']['qualified_on_common_inventory']
        both_residuals = residual_result['native_tolerance_met'] and residual_result['independent_tolerance_met']
        checked = bool(both_residuals and residual_result['complete_value_comparison_met'] and
                       residual_result['weak_affine_values_gate_met'] and finite and
                       physical_result['independent_accounting_identities_met'] and native_checks['passed'] is not False)
        converged = provenance_result['native_converged_observation']
        native_report_agreement = (report_metadata['native_target_converged'] is converged
                                   if report_metadata['available'] else None)
        checked = checked and native_report_agreement is not False and not report_metadata.get('binding_mismatch', False)
        result.update(all_checked_numerical_gates_met=checked,
                      durable_and_final_report_convergence_agree=native_report_agreement,
                      corroborated_native_convergence=bool(converged is True and checked),
                      terminal_classification=('complete_iterate_exceeds_coupled_tolerance' if not both_residuals else
                          'low_residual_without_native_convergence_observation' if converged is not True else
                          'native_convergence_corroborated_on_checked_gates' if checked else
                          'native_convergence_observed_with_unmet_terminal_gates'))
        require(identity(progress, prefix['bytes']) == prefix, 'Durable prefix changed during analysis')
        for item in result['inputs'].values():
            require(identity(item['path']) == item, 'Bound file changed during analysis '+item['path'])
        require(git(source, 'rev-parse', 'HEAD') == args.expected_analysis_source and
                not git(source, 'status', '--porcelain'), 'Analysis checkout changed')
        producer = objects['producing_checkout']
        require(git(producer, 'rev-parse', 'HEAD') == args.expected_producer_source and
                not git(producer, 'status', '--porcelain'), 'Producing checkout changed during audit')
        result['independent_checks'] = analyzer.audit.checks
        result['success'] = True
    except Exception as error:
        result['failures'].append(dict(phase=result.get('current_phase'), type=type(error).__name__, message=str(error)))
        result.update(terminal_classification='independent_audit_incomplete',
                      all_checked_numerical_gates_met=None, corroborated_native_convergence=None)
    result.update(status='completed' if result['success'] else 'incomplete_or_failed', current_phase=None, finished_utc=stamp())
    save()
    return 0 if result['success'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
