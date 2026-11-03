"""Structural/hash integrity of the compact evidence; no numerical rerun."""
import argparse, ast, hashlib, json
from pathlib import Path

here = Path(__file__).resolve().parent
repo = here.parents[3]
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
load = lambda name: json.loads((here/name).read_text())
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--write-manifest', action='store_true')
args = parser.parse_args()
checks = 0
for path in sorted(here.iterdir()):
    if path.suffix == '.json': json.loads(path.read_text()); checks += 1
    if path.suffix == '.py': ast.parse(path.read_text(), filename=str(path)); checks += 1
for path in sorted(here.glob('*curation-index.json')):
    entries = json.loads(path.read_text()).get('files', {})
    if isinstance(entries, dict): entries = [dict(name=key, **value) for key, value in entries.items()]
    for entry in entries:
        name = entry.get('name', entry.get('file'))
        expected = entry.get('curated_sha256', entry.get('sha256'))
        if name and expected:
            assert (here/name).is_file(), (path.name, name)
            assert sha(here/name) == expected, (path.name, name)
            checks += 1
ledger = load('first-failure-ledger.json')
assert ledger['status'] == 'completed_bounded_diagnosis' and ledger['M8_status'] == 'OPEN'
assert not ledger['physical_model_qualified'] and not ledger['production_constitutive_adoption']
assert ledger['bounded_matrix']['used_slots'] == [1,2,3,4,5,6]
assert ledger['bounded_matrix']['additional_evolution_slots'] == 0
assert len(ledger['archived_first_failure']['accepted_epochs']) == 12
for epoch in ledger['archived_first_failure']['accepted_epochs']:
    assert epoch['minimum_lambda'] == 1 and epoch['normalization_capped_volume'] == 0
    assert epoch['minimum_inner_k'] > 0 and epoch['minimum_outer_epsilon'] > 0
    checks += 4
for name, entry in ledger['components'].items():
    assert sha(here/name) == entry['sha256'], name
    assert (here/name).stat().st_size == entry['bytes'], name
    checks += 2
for name in ('nonlinear-b1-ledger.json', 'nonlinear-b2-ledger.json'):
    value = load(name)
    assert len(value['failed_stage_state']) == len(value['failed_stage_history']) == len(value['failed_stage_residual']) == 425
    assert value['success'] and value['failures'] == []
    assert not value['fresh_jacobian_same_stage_retry'].get('history_advanced', False)
    checks += 3
binding = load('source-binding-receipt.json')
assert binding['status'] == 'passed' and binding['compiled_inputs_byte_identical']
assert not binding['all_input_manifest_identity_claimed']
for name, expected in binding['helpers'].items():
    assert sha(repo/name) == expected, name
    ast.parse((repo/name).read_text(), filename=name)
    checks += 2
runtime = load('runtime-summary.json')
for name, epoch in runtime['epochs'].items():
    assert epoch['runtime_images_count'] == 100 and epoch['source_inputs_count'] == 3921, name
    assert epoch['success'] and epoch['clean'], name
    checks += 3
assert all(all(value.values()) for value in runtime['before_after_comparisons'].values())
assert load('nonlinear-b2-ledger.json')['observed_Jacobian_accuracy_failures']
points = load('archive-points-summary.json')
assert points['success'] and points['completed_epoch_count'] == 12
assert points['points_file']['records'] == 18180 and len(points['points_file']['sha256']) == 64
assert load('focused-checks-receipt.json')['status'] == 'passed'
assert load('kernel-summary.json')['outer_proof']['original_broad_domain_status'] == 'FAILED, unchanged'
assert not load('kernel-summary.json')['production_adopted']
paths = sorted(path for path in here.iterdir() if path.is_file() and path.name != 'evidence-manifest.json')
paths += [here/'../../cls-epawf-revised-m8-r1.md', here/'../../cls-epawf-revised-m8.md']
paths += sorted((repo/'solvers/IncNavierStokesSolver/Tests/RANS').glob('cls_epawf_m8_r1_*.py'))
entries = {str(path.relative_to(here)) if path.is_relative_to(here) else '../../../../'+str(path.relative_to(repo)):
    dict(sha256=sha(path), bytes=path.stat().st_size) for path in paths}
manifest = dict(schema=1, method='Local byte/JSON/AST/receipt consistency only; no numerical work',
    checks=checks, external_large_artifacts='Not copied; exact Linux paths/hashes retained in component ledgers', files=entries)
path = here/'evidence-manifest.json'
if args.write_manifest:
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True)+'\n')
else:
    recorded = load(path.name)
    assert recorded['files'] == manifest['files'], 'Evidence changed after manifest freeze'
print(json.dumps(dict(status='passed', structural_checks=checks, indexed_files=len(entries), numerical_re_evaluation=False)))
