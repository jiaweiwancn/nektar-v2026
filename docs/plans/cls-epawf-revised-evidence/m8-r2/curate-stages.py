"""Curate bounded native stage records without rerunning numerical operators."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import struct

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--run', type=Path, required=True)
p.add_argument('--analysis', type=Path, action='append', default=[])
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
run = a.run.resolve()
def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for part in iter(lambda: stream.read(1024*1024), b''):
            digest.update(part)
    return digest.hexdigest()
def file(path):
    path = Path(path)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=sha(path))
def vector(values):
    return hashlib.sha256(struct.pack('<'+str(len(values))+'d', *values)).hexdigest()
def domain(value):
    names = ('admissible', 'scope', 'phase', 'failure', 'minimum_dimensionless_margin',
        'limiting_inventory', 'limiting_point', 'minimum_inner_k', 'minimum_outer_k',
        'minimum_outer_epsilon', 'minimum_sampled_K', 'minimum_projected_trace_epsilon',
        'minimum_b', 'minimum_lambda', 'cutoff_range', 'inventory_counts')
    result = {key: value.get(key) for key in names}
    pattern = value.get('profile_branch_pattern', '')
    result['profile_branch_pattern_sha256'] = hashlib.sha256(pattern.encode()).hexdigest()
    result['profile_branch_counts'] = dict(Counter(pattern))
    return result

receipt_path = run/'receipt.json'
receipt = json.loads(receipt_path.read_text())
assert receipt['status'] != 'running', 'Do not curate an unfinished family as a final result'
result = dict(schema=1, milestone='M8-R2 F', created_utc=datetime.now(timezone.utc).isoformat(),
    curator=file(__file__), run_receipt=file(receipt_path),
    method='Read-only extraction/hashes of native records; no numerical reexecution',
    native_calls=0, operator_reconstructions=0, mode=receipt['mode'],
    producer_source=receipt['source_commit'], producer_source_end=receipt['source_commit_end'],
    runtime_unchanged=receipt['runtime_unchanged'], source_clean=receipt['source_clean_end'],
    worker=dict(ranks=receipt['ranks'], threads=receipt['numerical_threads'],
        maximum_concurrent_workers=receipt['maximum_concurrent_workers']),
    matrix_encoding=receipt['matrix_encoding'],
    preflight=dict(path=receipt['preflight_path'], sha256=receipt['preflight_sha256']),
    target=dict(dt=receipt['target_dt'], time=receipt['target_time'], alpha=400.,
        fixed_history=True, native_history_advance=False),
    allowed_alpha_levels=receipt['alpha_levels'], maximum_Newton_updates=12,
    runtime_before=file(run/'runtime-before.json'), runtime_after=file(run/'runtime-after.json'),
    stage_plan=file(run/'cases/stage-plan.json'), cases=[], independent_analysis=[])
pin = json.loads((run/'runtime-before.json').read_text())
result['producing_program'] = dict(path=str(Path(pin['build_root'])/'IncNavierStokesSolver'),
    sha256=pin['runtime_artifacts']['build/IncNavierStokesSolver'])
for row in receipt['cases']:
    xml = Path(row['xml_path'])
    assert sha(xml) == row['xml_sha256']
    native_command = next(item for item in receipt['commands'] if item['name'] == 'native-'+row['inner_stress'])
    item = dict(inner_stress=row['inner_stress'], xml=file(xml),
        native_seconds=row['native_seconds'], native_exit_code=row['native_exit_code'],
        command=native_command, target_attempted=row.get('target_attempted'),
        target_true_coupled_convergence=row.get('target_converged'),
        target_native_physical_range_flag=row.get('target_physical_range'),
        target_case_report_available=bool(row.get('report_path')),
        target_flags_origin='Final case-report receipt only; missing flags stay unknown. The strict stage assessment separately inspects actual logs and complete JSONL records.',
        note='Native range flag is separate from independent operator/action/accounting qualification',
        levels=[], nonlinear_solves=[], matrices=[])
    if row.get('report_path'):
        path = Path(row['report_path'])
        assert sha(path) == row['report_sha256']
        raw = json.loads(path.read_text())
        item.update(report=file(path), initial_state_float64_LE_sha256=vector(raw['initial']),
            signature_sha256=hashlib.sha256(raw['signature'].encode()).hexdigest(),
            native_time=raw['time'], native_steps=raw['steps'],
            native_history_float64_LE_sha256=[vector(value) for value in raw['native_history']],
            native_history_times=raw['native_history_times'],
            initial_native_snapshot=raw.get('stage_test_initial_native_history'),
            accepted_physical_states=len(raw['accepted_states']),
            stage_test=raw['stage_test'], failure=raw.get('failure'))
        # The initial snapshot contains one complete coefficient vector and
        # four history slots. Retain hashes rather than duplicate heavy arrays.
        snapshot = item.get('initial_native_snapshot')
        if snapshot:
            item['initial_native_snapshot'] = dict(time=snapshot['time'], steps=snapshot['steps'],
                state_sha256=vector(snapshot['state']), times=snapshot['times'],
                history_sha256=[vector(value) for value in snapshot['history']])
        for level in raw['stage_test_levels']:
            entry = {key: level.get(key) for key in ('alpha', 'time', 'attempted', 'reason',
                'converged', 'true_coupled_residual_converged', 'frozen_velocity_pressure_diagnostic',
                'iterations', 'failure', 'residual_verification_failure', 'qualification_audit_failure',
                'scaled_block_norms', 'physical_range_passed', 'native_history_advanced')}
            if level['attempted']:
                entry.update(state_sha256=vector(level['state']), history_sha256=vector(level['history']),
                    residual_sha256=vector(level['residual']), constraints=level['constraints'],
                    scalar_domain=domain(level['scalar_domain']), audit=level['audit'])
            item['levels'].append(entry)
    if row.get('trace_path'):
        path = Path(row['trace_path'])
        assert sha(path) == row['trace_sha256']
        item['trace'] = file(path)
        with path.open() as stream:
            for line in stream:
                record = json.loads(line)
                solve = {key: record.get(key) for key in ('policy', 'inner_stress_model', 'time',
                    'alpha', 'tolerance', 'maximum_iterations', 'converged', 'terminal',
                    'frozen_velocity_pressure_diagnostic', 'history_advanced', 'residual_unchanged',
                    'jacobian_rebuild_policy', 'scaled_block_norms')}
                solve.update(state_sha256=vector(record['state']), history_sha256=vector(record['history']),
                    residual_sha256=vector(record['residual']), domain=domain(record['domain']), iterations=[])
                for iteration in record['iterations']:
                    observed = {key: iteration.get(key) for key in ('iteration', 'scaled_block_norms',
                        'linear_scaled_block_norms', 'dimensionless_step', 'fixed_location_maximum_fraction',
                        'boundary_limited', 'limiting_inventory', 'limiting_point', 'limiting_value',
                        'accepted_fraction', 'accepted')}
                    trials = iteration['trials']
                    observed.update(state_sha256=vector(iteration['state']), domain=domain(iteration['domain']),
                        jacobian_file=iteration['jacobian_binary_file'],
                        column_stencils=dict(Counter(column['stencil'] for column in iteration['columns'])),
                        trial_classifications=dict(Counter(trial['classification'] for trial in trials)),
                        branch_changed_trials=sum(trial.get('profile_branches_changed', False) for trial in trials),
                        trials=[{key: trial.get(key) for key in ('fraction', 'classification',
                            'scaled_block_norms', 'type', 'message')} | dict(domain=domain(trial['domain']))
                            if 'domain' in trial else {key: trial.get(key) for key in ('fraction',
                                'classification', 'type', 'message')} for trial in trials])
                    solve['iterations'].append(observed)
                item['nonlinear_solves'].append(solve)
    # Include even a matrix written before a factorization/reporting failure.
    for path in sorted(xml.parent.glob('case.feasible.jsonl.solve*.J*.bin')):
        item['matrices'].append(file(path))
    result['cases'].append(item)
for path in a.analysis:
    path = path.resolve()
    data = json.loads(path.read_text())
    result['independent_analysis'].append(dict(file=file(path), success=data.get('success'),
        analysis_success=data.get('analysis_success'), source=data.get('source_commit'),
        producing_source=data.get('producing_source_commit'), failures=data.get('failures')))
result['scope'] = ('Native bounded frozen-stage observations, with separately linked independent checks. '
    'A failed search is not a proof of root absence. No auxiliary state is a physical trajectory sample. M8 remains OPEN.')
a.output.parent.mkdir(parents=True, exist_ok=True)
with a.output.open('x') as stream:
    json.dump(result, stream, indent=2, allow_nan=False)
    stream.write('\n')
print(json.dumps(dict(output=str(a.output), sha256=sha(a.output), cases=len(result['cases']))))
