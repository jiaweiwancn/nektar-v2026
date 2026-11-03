"""Read-only curation of stored R2 B results; never imports numerical oracles."""
import argparse
import hashlib
import json
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--run', type=Path, required=True)
p.add_argument('--generator-sha256', required=True)
a = p.parse_args()
run = a.run.resolve()
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
record = lambda path: dict(path=str(path), bytes=path.stat().st_size, sha256=sha(path))
validation_path, receipt_path = run/'results/validation.json', run/'receipt.json'
v, receipt = (json.loads(path.read_text()) for path in (validation_path, receipt_path))
assert v['success'] and v['verification_complete'] and not v['failures'] and not v['commands']
assert receipt['status'] == 'passed' and receipt['runtime_unchanged']
assert receipt['source_commit'] == receipt['source_commit_end'] == v['producing_source_commit']
assert receipt['analysis_commit'] == receipt['analysis_commit_end']
summary = dict(schema=1, milestone='M8-R2 B', status=v['status'],
    summary_method='Read-only extraction and hashes of stored Linux results; no numerical or native reexecution',
    summary_generator_sha256=a.generator_sha256,
    producing_source_commit=v['producing_source_commit'], analysis_source_commit=receipt['analysis_commit'],
    independent_checks=v['checks'], native_kernel_checks=v['kernel']['native_checks'],
    compiled_kernel_cases=len(v['kernel']['cases']), native_fixtures=len(v['cases']),
    native_calls_in_final_replay=len(v['commands']), native_kernel_reexecuted=v['kernel_native_reexecution'],
    solver_stage_launched=v['solver_stage_launched'], trajectory_launched=v['trajectory_launched'],
    cpu=receipt['cpu'], ranks=receipt['ranks'], threads=receipt['numerical_threads'],
    seconds=v['seconds'], runtime_unchanged=receipt['runtime_unchanged'],
    analysis_and_runtime_clean=(receipt['source_clean_start'] and receipt['source_clean_end']
        and receipt['analysis_clean_start'] and receipt['analysis_clean_end']),
    authoritative_inputs={name: record(path) for name, path in [
        ('validation', validation_path), ('wrapper_receipt', receipt_path),
        ('runtime_before', run/'runtime-before.json'), ('runtime_after', run/'runtime-after.json'),
        ('log', run/'suite.log')]},
    helper_sha256=v['helper_sha256'], kernel_probe_sha256=v['kernel_probe_sha256'],
    wall_probe_sha256=v['wall_probe_sha256'], kernel_report_sha256=v['kernel_sha256'],
    template_sha256=v['template_sha256'], imported_source_files=v['imported_source_files'],
    wrapper_suite_command=next(row['argv'] for row in receipt['commands'] if row['name'] == 'suite'),
    errors=v['errors'], failures=v['failures'], invalid_values=v['invalid_values'],
    kernel_cases=[{key: row[key] for key in ('name', 'selector', 'active', 'directional_only')}
                  for row in v['kernel']['cases']],
    cases=[], all_saturation_statistics_independently_certified=v['all_saturation_statistics_independently_certified'],
    unresolved_diagnostics=v['unresolved_diagnostics'], independent_outer_evaluations=v['independent_outer_evaluations'])
for row in v['cases']:
    result = row['result']
    raw_path = Path(row['input_raw'])
    assert sha(raw_path) == row['raw_sha256'] and sha(Path(row['input_xml'])) == row['xml_sha256']
    raw = json.loads(raw_path.read_text())
    audit = raw['constitutive_audit']
    summary['cases'].append(dict(name=row['name'], selector=result['inner_stress_model'],
        fixture=result['fixture'], chi=result['chi'], native_reexecution=row['native_reexecution'],
        input_raw=record(raw_path), input_xml=record(Path(row['input_xml'])),
        input_log=record(raw_path.with_name('run.log')), weak=result['weak'],
        actual_source_budgets=audit['robust_source'], actual_wall_sources=audit['wall_sources'],
        actual_source_inner_covariance_minimum=audit['source_volume_inner'][3],
        actual_inner_zero_covariance=audit['inner_zero_covariance'],
        actual_inner_statistics=audit.get('inner_stress'),
        actual_inner_statistic_maximum_point=audit.get('inner_stress_max_ratio_point'),
        actual_energy_budget=raw['energy_budget'],
        supplied_diagnostic_state=result['supplied_audit_state_key'],
        coupled_action=result['coupled_action'], lifted_traction=result['lifted_traction'],
        saturation_statistics=result['saturation_statistics'],
        physical_signature_sha256=hashlib.sha256(result['physical_signature'].encode()).hexdigest(),
        qualification=result['qualification']))
cap = dict(v['retained_pressure_cap'])
rays = cap.pop('rays')
cap['rays_count'] = len(rays)
cap['b_range'] = [min(row['b'] for row in rays), max(row['b'] for row in rays)]
cap['lambda_range'] = [min(row['lambda_value'] for row in rays), max(row['lambda_value'] for row in rays)]
cap['chi_range'] = [min(row['chi'] for row in rays), max(row['chi'] for row in rays)]
cap['first_attenuated_ray'] = next(row for row in rays if row['lambda_value'] < 1)
summary['retained_pressure_cap'] = cap
pressure_validation_path = next(Path(path) for path in cap['input_files'] if path.endswith('/results/validation.json'))
pressure = json.loads(pressure_validation_path.read_text())
summary['separate_pressure_action_verification'] = dict(input=record(pressure_validation_path),
    source_commit=pressure['source_commit'], success=pressure['success'],
    case_checks={name: row['checks'] for name, row in pressure['cases'].items()},
    physical_action_tolerance=1e-7,
    scope='Separate C gate; B mixed-action component errors and its 8e-6 threshold do not substitute for the pressure physical-norm threshold')
summary['retained_failed_attempts'] = []
for directory in sorted(set(run.parent.glob('inner-*')) | set(run.parent.parent.glob('inner-*'))):
    path = directory/'receipt.json'
    if not path.is_file():
        continue
    old = json.loads(path.read_text())
    if old.get('status') != 'failed':
        continue
    failure = dict(wrapper=record(path), producing_source_commit=old.get('source_commit'),
                   analysis_source_commit=old.get('analysis_commit'), runtime_unchanged=old.get('runtime_unchanged'))
    for name in ('suite.log', 'results/validation.json'):
        item = directory/name
        if item.is_file():
            failure[name] = record(item)
    summary['retained_failed_attempts'].append(failure)
summary['interpretation'] = dict(
    passed='Compiled kernel, active/inactive actual weak loads, coupled mixed actions, full stress/production/traction/energy/source consistency, and exact native zero-k covariance controls on declared finite fixtures',
    unresolved='Coefficient-derived maximum a/c in the zero-k fixture remains unresolved at original3e-10; actual exported primitive semantics are separately checked without changing tolerances',
    not_established='Qualified implicit-stage root, fixed outer-box coverage, polynomial-everywhere positivity, fully discrete positivity, curved MMS, long-time stability, scalability, terrain or physical validation',
    milestone_status='M8 remains OPEN')
print(json.dumps(summary, indent=2, allow_nan=False))
