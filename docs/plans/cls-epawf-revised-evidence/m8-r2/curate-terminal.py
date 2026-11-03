#!/usr/bin/env python3
"""Compact one completed R2 terminal audit; stdlib metadata and hashes only.

No scientific library, operator reconstruction, native solver, Git operation or
remote command is invoked. Success means all terminal-audit phases completed
with unchanged bindings. Convergence and finite-range observations stay separate.
Run on the host retaining the original bound inputs and phase checkpoints.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import tempfile


PRODUCER = '9f20b3053c303e8805f5183f668a8ec468d9656b'
EXECUTABLE = 'a90a196fe554d591fa39c6556166bf58ad43afa787b6b8e8d3d5cbf3355a5668'
RUNTIME_KEYS = ('source_inputs', 'runtime_artifacts', 'runtime_bindings',
    'runtime_environment', 'python_packages', 'tool_images')
CPP_NAMES = ('RANSEquilibriumWallEvolution.cpp', 'RANSWallEvolutionPatch.cpp',
    'RANSWallEvolutionFeasible.cpp')
PHASES = ('provenance_and_complete_record', 'constructor_and_static_inventory',
    'common_scalar_domain', 'current_projected_constraints', 'true_coupled_residual',
    'selected_source_and_direct_accounting', 'matching_full_traction_and_k_flux',
    'strong_continuity', 'physical_mass_and_complete_energy',
    'common_inventory_physical_ranges', 'optional_native_budget_crosschecks')
SCOPE = ('One complete native terminal record and its completed independent audit; '
    'no Newton-path, current-action, target-search, final-history or physical-trajectory '
    'qualification, root-absence inference or M8 closure.')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def canonical(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
        allow_nan=False).encode()).hexdigest()


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for part in iter(lambda: stream.read(1024*1024), b''):
            result.update(part)
    return result.hexdigest()


def load(path):
    def invalid(value):
        raise ValueError('Nonfinite JSON constant: '+value)
    return json.loads(Path(path).read_text(), parse_constant=invalid)


def pick(value, names):
    return {name: value[name] for name in names if name in value}


def without(value, names):
    return {name: item for name, item in value.items() if name not in names}


class Bindings:
    """Hash immutable files once and reject changes while extracting metadata."""
    def __init__(self):
        self.files, self.signatures = {}, {}

    @staticmethod
    def signature(path):
        info = path.stat()
        return info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns

    def bind(self, path, expected=None):
        path = Path(path).resolve()
        name = str(path)
        if name not in self.files:
            before = self.signature(path)
            item = dict(path=name, bytes=path.stat().st_size, sha256=digest(path))
            require(before == self.signature(path), 'File changed while hashing: '+name)
            self.files[name], self.signatures[name] = item, before
        item = self.files[name]
        if expected is not None:
            require(isinstance(expected, dict) and expected.get('path') == name and
                expected.get('bytes') == item['bytes'] and expected.get('sha256') == item['sha256'],
                'Changed or incomplete file binding: '+name)
        return item

    def finish(self):
        for name, signature in self.signatures.items():
            require(self.signature(Path(name)) == signature, 'Bound file changed during curation: '+name)


def completed_phases(receipt, result):
    require(receipt.get('schema') == result.get('schema') == 1, 'Unsupported terminal-audit schema')
    require(receipt.get('success') is True and result.get('success') is True and
        result.get('analysis_success') is True and receipt.get('status') == result.get('status') == 'completed',
        'Terminal audit is incomplete or unsuccessful')
    require(receipt.get('failures') == result.get('failures') == [] and
        receipt.get('current_phase') is None, 'Terminal audit retains failures or an active phase')
    require(receipt.get('scope') == result.get('scope'), 'Receipt/result scopes differ')
    rows = receipt.get('phases', [])
    require([row.get('name') for row in rows] == list(PHASES) and
        all(row.get('status') == 'completed' for row in rows) and
        set(result.get('phases', {})) == set(PHASES), 'Missing, reordered or incomplete terminal phase')
    for name in ('native_calls',):
        require(receipt.get(name) == result.get(name) == 0, 'Not a read-only terminal audit')
    for name in ('physical_trajectory_qualification', 'full_native_case_protocol_assessed',
                 'final_native_history_slots_assessed', 'M8_closed'):
        require(receipt.get(name) is False and result.get(name) is False,
            'Unexpected broader qualification claim: '+name)
    require(receipt.get('analyzed_iteration_paths') == receipt.get('fresh_residual_endpoints') == 0 and
        receipt.get('matrix_files') == {} and receipt.get('derivative_accuracy_assessed') is False and
        receipt.get('target_search_protocol_assessed') is False and receipt.get('root_absence_proven') is False and
        result.get('iteration_path_assessed') is False and
        result.get('current_action_or_derivative_gates_assessed') is False,
        'Audit includes a path, endpoint, matrix or broader target claim')
    return rows


def phase_result_binding(name, saved, final):
    """Permit only the wrapper's two-field final prefix enrichment.

    provenance.record_snapshot aliases receipt.active_trace_prefix in the
    producing wrapper. Its final provenance pass adds these fields after the
    completed phase checkpoint has been written. Every other value, including
    all state/history/residual hashes, must retain its canonical JSON binding.
    """
    if canonical(saved) == canonical(final):
        return None
    message = 'Phase/result binding mismatch: '+name
    require(name == 'provenance_and_complete_record' and
        isinstance(saved, dict) and isinstance(final, dict), message)
    before, after = saved.get('record_snapshot'), final.get('record_snapshot')
    fields = ('prefix_sha256_end', 'trace_bytes_end')
    require(isinstance(before, dict) and isinstance(after, dict) and
        set(after)-set(before) == set(fields) and not set(before)-set(after), message)
    require(after['prefix_sha256_end'] == before['prefix_sha256'] and
        type(after['trace_bytes_end']) is int and type(before['prefix_bytes']) is int and
        after['trace_bytes_end'] >= before['prefix_bytes'] > 0, message)
    expected = dict(saved, record_snapshot=dict(before, **{key: after[key] for key in fields}))
    require(canonical(expected) == canonical(final), message)
    return dict(phase=name, location='result.record_snapshot',
        fields_added_after_completed_checkpoint={key: after[key] for key in fields},
        reason='Wrapper shared prefix dictionary receives final provenance metadata after the phase checkpoint',
        ending_prefix_hash_matches_starting_hash=True,
        ending_trace_bytes_cover_frozen_prefix=True,
        all_other_canonical_result_fields_unchanged=True)


def command_options(receipt):
    argv = receipt['command_argv']
    require(isinstance(argv, list) and len(argv) > 1 and len(argv[1:]) % 2 == 0,
        'Malformed recorded terminal command')
    options = {}
    for flag, value in zip(argv[1::2], argv[2::2]):
        require(isinstance(flag, str) and flag.startswith('--') and flag not in options,
            'Duplicate or malformed recorded terminal argument')
        options[flag] = value
    for name in ('source', 'expected-commit', 'preflight', 'preflight-session',
                 'stage-session', 'trace', 'runtime-pin', 'output'):
        require('--'+name in options, 'Missing recorded terminal argument: '+name)
    return argv[0], options


def prefix_binding(prefix):
    """Rehash only the published prefix/selected record; later appends are allowed."""
    count, offset, selected = (prefix[key] for key in ('prefix_bytes', 'selected_offset', 'selected_bytes'))
    require(all(type(value) is int for value in (count, offset, selected)) and
        count > 0 and offset >= 0 and selected > 0 and offset+selected == count,
        'Invalid complete-record byte range')
    require(prefix['prefix_sha256'] == prefix['prefix_sha256_end'] and
        prefix['later_append_permitted'] is True and prefix['iteration_entries_analyzed'] == 0,
        'Published prefix changed or iteration paths were analyzed')
    whole, record = hashlib.sha256(), hashlib.sha256()
    cursor = 0
    with Path(prefix['path']).open('rb') as stream:
        while cursor < count:
            block = stream.read(min(1024*1024, count-cursor))
            require(bool(block), 'Original published native prefix was truncated')
            whole.update(block)
            if cursor+len(block) > offset:
                record.update(block[max(0, offset-cursor):])
            cursor += len(block)
    require(whole.hexdigest() == prefix['prefix_sha256'] and
        record.hexdigest() == prefix['selected_sha256'], 'Original published native prefix/record changed')
    return dict(prefix_rehashed=True, selected_record_bytes_rehashed=True,
        later_trace_appends_permitted=True, iteration_entries_inspected=False)


def compact_results(phases):
    """Copy observations, dropping large traces and lossless branch encodings."""
    domain = phases['common_scalar_domain']
    require(domain.get('admissible') is True and domain.get('phase') == 'complete' and
        isinstance(domain.get('coordinate_correspondence'), dict), 'Common scalar audit did not complete')
    compact_domain = without(domain, ('independent_trace', 'independent_branch_context', 'profile_branch_pattern'))
    trace = domain['independent_trace']
    compact_domain['independent_trace_summary'] = without(trace, ('samples', 'projected_target_coefficients'))
    pattern = domain['profile_branch_pattern']
    compact_domain['profile_branch_summary'] = dict(count=len(pattern), population=dict(Counter(pattern)),
        pattern_sha256=hashlib.sha256(pattern.encode()).hexdigest())
    compact_domain['independent_branch_context'] = {
        name: without(value, ('lossless_RLE',)) if isinstance(value, dict) else value
        for name, value in domain['independent_branch_context'].items()}
    residual = dict(phases['true_coupled_residual'])
    weak = dict(residual['weak_residual'])
    weak['exact_affine_weak_continuity'] = without(weak['exact_affine_weak_continuity'], ('exact_values',))
    residual['weak_residual'] = weak
    constraints = phases['current_projected_constraints']
    compact_constraints = without(constraints, ('entries',))
    compact_constraints.update(entries_count=len(constraints['entries']),
        entries_canonical_sha256=canonical(constraints['entries']))
    traction = phases['matching_full_traction_and_k_flux']
    values = traction['full_total_traction_and_k_flux']
    require(len(values) == 8 and 'tangential_traction_RMS' in traction, 'Incomplete full matching-interface budget')
    matching = dict(traction)
    matching['named_observations'] = dict(zip(('area', 'reduced_pressure_jump_integral',
        'normal_total_traction_jump_RMS', 'full_total_traction_jump_RMS',
        'maximum_absolute_normal_total_traction_jump', 'k_flux_jump_RMS',
        'inner_normal_velocity_integral', 'normal_total_traction_jump_integral'), values))
    matching['named_observations']['tangential_total_traction_jump_RMS'] = traction['tangential_traction_RMS']
    continuity = dict(phases['strong_continuity'])
    require(len(continuity['independent_strong_continuity']) == 5, 'Incomplete strong-continuity budget')
    continuity['named_observations'] = dict(zip(('inner_divergence_L2', 'outer_divergence_L2',
        'inner_divergence_integral', 'outer_divergence_integral', 'external_normal_velocity_flux'),
        continuity['independent_strong_continuity']))
    return dict(common_scalar_domain=compact_domain,
        static_inventory=phases['constructor_and_static_inventory'],
        projected_constraints=compact_constraints, true_coupled_and_scalar_residual=residual,
        matching_interface=matching, strong_continuity=continuity,
        source_and_direct_accounting=phases['selected_source_and_direct_accounting'],
        mass_and_complete_energy=phases['physical_mass_and_complete_energy'],
        common_inventory_physical_ranges=phases['common_inventory_physical_ranges'],
        optional_native_budget_crosschecks=phases['optional_native_budget_crosschecks'])


def curate(run):
    run = Path(run).resolve()
    binder = Bindings()
    receipt_pin, result_pin = (binder.bind(run/name) for name in ('receipt.json', 'validation.json'))
    receipt, result = (load(item['path']) for item in (receipt_pin, result_pin))
    rows = completed_phases(receipt, result)
    phases = result['phases']
    checkpoints, enrichments = [], []
    for index, row in enumerate(rows, 1):
        expected = run/('phase-'+str(index).zfill(2)+'-'+row['name']+'.json')
        require(Path(row['checkpoint']).resolve() == expected, 'Phase checkpoint changed location')
        pin = binder.bind(expected, row['completed_checkpoint'])
        data = load(expected)
        require(data.get('schema') == 1 and data.get('scope') == receipt['scope'] and
            data.get('phase') == without(row, ('completed_checkpoint',)),
            'Phase checkpoint metadata binding mismatch: '+row['name'])
        enrichment = phase_result_binding(row['name'], data.get('result'), phases[row['name']])
        if enrichment is not None:
            enrichments.append(dict(enrichment, checkpoint=pin))
        checkpoints.append(dict(name=row['name'], seconds=row['seconds'], checkpoint=pin))
    driver, options = command_options(receipt)
    require(Path(options['--output']).resolve() == run, 'Recorded terminal output differs from --run')
    require(receipt.get('source_clean_start') is True and receipt.get('source_clean_end') is True and
        receipt['source_commit'] == receipt['source_commit_end'] == options['--expected-commit'],
        'Analysis source did not retain its clean commit binding')
    inputs = receipt['inputs']
    require(isinstance(inputs, dict) and inputs, 'No exact terminal input pins')
    for name, item in inputs.items():
        require(name == item.get('path'), 'Input key/path mismatch')
        binder.bind(name, item)
    required_paths = [driver]+[options['--'+name] for name in
        ('preflight', 'preflight-session', 'stage-session', 'runtime-pin')]
    require(all(str(Path(path).resolve()) in inputs for path in required_paths), 'Missing required input pin')
    require(inputs[str(Path(driver).resolve())]['sha256'] == receipt['wrapper_sha256'], 'Terminal wrapper hash differs')
    analysis_source = Path(options['--source']).resolve()
    require(analysis_source in Path(driver).resolve().parents and receipt['imported_source_files_start'],
        'Missing source-bound wrapper or imported analysis modules')
    for path, expected in receipt['imported_source_files_start'].items():
        require(analysis_source in Path(path).resolve().parents and path in inputs and
            inputs[path]['sha256'] == expected, 'Imported source pin differs: '+path)
    runtime = load(options['--runtime-pin'])
    require(runtime.get('success') is True and runtime.get('clean') is True and
        receipt['producing_source_commit'] == runtime['source_commit'] == PRODUCER and
        runtime['source_root'] == receipt['producing_source_path'], 'Producing source/runtime pin differs')
    require(set(receipt['producing_runtime_sections']) == set(RUNTIME_KEYS) and
        all(canonical(runtime[key]) == receipt['producing_runtime_sections'][key] for key in RUNTIME_KEYS),
        'Producing runtime section changed')
    solver = receipt['producing_solver']
    solver_key = 'build/IncNavierStokesSolver'
    require(solver['sha256'] == runtime['runtime_artifacts'][solver_key] == EXECUTABLE and
        solver['path'] == runtime['runtime_resolved_paths'][solver_key] and
        inputs.get(solver['path']) == solver, 'Producing executable binding differs')
    cpp = {}
    for name in CPP_NAMES:
        relative = 'solvers/IncNavierStokesSolver/EquationSystems/'+name
        path = str((Path(receipt['producing_source_path'])/relative).resolve())
        require(path in inputs and inputs[path]['sha256'] == runtime['source_inputs'][relative],
            'Producing C++ file binding differs: '+name)
        cpp[name] = inputs[path]
    provenance = phases['provenance_and_complete_record']
    prefix = receipt['active_trace_prefix']
    require(provenance['record_snapshot'] == prefix and
        provenance['producing_solver'] == solver and provenance['producing_source_commit'] == PRODUCER and
        provenance['matrix_files_used'] == provenance['iteration_paths_analyzed'] == 0,
        'Native record provenance differs')
    snapshot = inputs.get(prefix['snapshot'])
    require(isinstance(snapshot, dict) and snapshot['sha256'] == prefix['selected_sha256'] and
        snapshot['bytes'] == prefix['selected_bytes'] and options['--trace'] == prefix['path'] and
        int(options.get('--record-index', 0)) == receipt['record_index'], 'Native snapshot binding differs')
    prefix_check = prefix_binding(prefix)
    binding = provenance['input_binding']
    residual = phases['true_coupled_residual']
    source = phases['selected_source_and_direct_accounting']
    require(provenance['selected_history_float64_LE_sha256'] == binding['preflight_h0_float64_LE_sha256'] and
        binding['native_converged_observation'] is residual['native_converged_observation'] and
        binding['native_frozen_scalar_diagnostic'] is residual['frozen_scalar_scope'] and
        binding['time'] == .0025 and binding['alpha'] in (3200, 1600, 800, 400) and
        binding['mode'] == ('FrozenScalarR2' if residual['frozen_scalar_scope'] else 'MatchedR2') and
        (not residual['frozen_scalar_scope'] or binding['alpha'] == 400), 'Selected terminal fields/h0 differ')
    native_audit = provenance['native_audit']
    available = native_audit.get('available')
    require(type(available) is bool, 'Missing optional native-accounting availability')
    flags = [source['native_source_budget_table_available'],
        phases['matching_full_traction_and_k_flux']['native_traction_budget_table_available'],
        phases['strong_continuity']['native_strong_continuity_table_available'],
        phases['physical_mass_and_complete_energy']['native_energy_budget_table_available'],
        phases['optional_native_budget_crosschecks']['performed']]
    require(all(flag is available for flag in flags), 'Optional native-budget availability changed')
    if available:
        require('--native-report' in options and
            inputs.get(options['--native-report']) == native_audit['report'], 'Optional native report binding differs')
    else:
        require('--native-report' not in options, 'Supplied native report became unavailable')
    observations = compact_results(phases)
    binder.finish()
    return dict(schema=1, success=True,
        success_definition='All terminal-only phases completed with unchanged bindings; convergence and physical gates are separate observations',
        scope=SCOPE, curated_utc=datetime.now(timezone.utc).isoformat(),
        curator_sha256=digest(__file__), receipt=receipt_pin, validation=result_pin,
        analysis=dict(source_commit=receipt['source_commit'], source_commit_end=receipt['source_commit_end'],
            source_clean_start=receipt['source_clean_start'], source_clean_end=receipt['source_clean_end'],
            driver=inputs[str(Path(driver).resolve())], command_argv=receipt['command_argv'],
            imported_source_files_start=receipt['imported_source_files_start'],
            numerical_threads=receipt['numerical_threads'], cpu_affinity=receipt['cpu_affinity'],
            started_utc=receipt['started_utc'], finished_utc=receipt['finished_utc'],
            elapsed_seconds=receipt['elapsed_seconds'], checks=receipt['checks'], errors=receipt['errors']),
        producer=dict(source_commit=PRODUCER, cpp=cpp, executable=solver,
            runtime_pin=inputs[str(Path(options['--runtime-pin']).resolve())],
            runtime_sections=receipt['producing_runtime_sections'],
            runtime_sections_are_bound_producing_pin_not_new_runtime_capture=True),
        input_pins=list(inputs.values()), completed_phases=checkpoints,
        checkpoint_final_enrichments=enrichments,
        selected_record=dict(record_index=receipt['record_index'], snapshot=prefix, prefix_check=prefix_check,
            input_binding=binding, selected_state_float64_LE_sha256=provenance['selected_state_float64_LE_sha256'],
            selected_history_float64_LE_sha256=provenance['selected_history_float64_LE_sha256'],
            selected_native_residual_float64_LE_sha256=provenance['selected_native_residual_float64_LE_sha256'],
            native_audit=native_audit),
        terminal_observations=dict(alpha=binding['alpha'], time=binding['time'], mode=binding['mode'],
            inner_stress_model=source['independent_source_and_direct_accounting']['inner_stress_model'],
            converged=binding['native_converged_observation'], terminal=binding['native_terminal_observation'],
            frozen_scalar_diagnostic=binding['native_frozen_scalar_diagnostic'], tolerance=residual['tolerance']),
        observations=observations, native_calls=0, scientific_reconstructions=0,
        iteration_path_assessed=False, current_action_or_derivative_gates_assessed=False,
        full_native_case_protocol_assessed=False, final_native_history_slots_assessed=False,
        physical_trajectory_qualification=False, M8_closed=False)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    output = args.output.resolve()
    require(not output.exists(), 'Use a fresh compact-evidence output path')
    answer = curate(args.run)
    protected = {item['path'] for item in answer['input_pins']}
    protected.update(item['checkpoint']['path'] for item in answer['completed_phases'])
    protected.update(answer[key]['path'] for key in ('receipt', 'validation'))
    protected.add(answer['selected_record']['snapshot']['path'])
    require(str(output) not in protected and output != Path(__file__).resolve(), 'Output would overwrite bound evidence/source')
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode='w', dir=output.parent, prefix=output.name+'.',
                                    suffix='.tmp', delete=False) as stream:
        temporary = Path(stream.name)
        try:
            json.dump(answer, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write('\n')
        except Exception:
            temporary.unlink()
            raise
    try:
        temporary.replace(output)
    finally:
        if temporary.exists():
            temporary.unlink()
    print(json.dumps(dict(success=True, output=str(output), sha256=digest(output),
        alpha=answer['terminal_observations']['alpha'],
        native_converged_observation=answer['terminal_observations']['converged'])))
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (OSError, ValueError, KeyError, TypeError) as error:
        print(json.dumps(dict(success=False, type=type(error).__name__, message=str(error))), file=sys.stderr)
        raise SystemExit(1)
