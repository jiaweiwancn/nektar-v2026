"""Pure provenance curation; no build, solver or numerical re-evaluation."""
import hashlib, json
from pathlib import Path

base = Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised')
r1 = base/'m8-r1'
out = r1/'runtime-curation-001'; out.mkdir(exist_ok=False)
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
canonical = lambda v: hashlib.sha256(json.dumps(v, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
summary = dict(schema=1, method='Field selection only; exact original manifests remain external',
    curator_sha256=sha(__file__), numerical_re_evaluation=False, epochs={}, receipts={}, logs={})
paths = {
    'diagnostics-build': r1/'build-001/runtime-fingerprint.json',
    'original-reproduction-before': r1/'reproduction-001/runtime-before.json',
    'original-reproduction-after': r1/'reproduction-001/runtime-after.json',
    'controlled-build': r1/'controlled/build-002/runtime-fingerprint.json',
    'admissible-reproduction-before': r1/'controlled/reproduction-001/runtime-before.json',
    'admissible-reproduction-after': r1/'controlled/reproduction-001/runtime-after.json',
    'focused-before': r1/'controlled/focused-checks-001/runtime-before.json',
    'focused-after': r1/'controlled/focused-checks-001/runtime-after.json',
    'controls-before': r1/'controlled/controls-001/runtime-before.json',
    'controls-after': r1/'controlled/controls-001/runtime-after.json',
}
for name, path in paths.items():
    data = json.loads(path.read_text()); assert data['success'] and data['clean']
    entry = dict(path=str(path), sha256=sha(path), bytes=path.stat().st_size,
        source_commit=data['source_commit'], source_root=data['source_root'],
        source_inputs_count=len(data['source_inputs']), source_inputs_canonical_sha256=canonical(data['source_inputs']),
        runtime_images_count=len(data['runtime_artifacts']), runtime_artifacts=data['runtime_artifacts'],
        runtime_bindings_canonical_sha256=canonical(data['runtime_bindings']),
        python_packages=data['python_packages'], tool_images=data['tool_images'],
        runtime_environment=data['runtime_environment'], success=data['success'], clean=data['clean'],
        source_selection=data.get('source_selection'), failures=data.get('failures'))
    entry['selected_source_inputs'] = {key: value for key, value in data['source_inputs'].items()
        if any(token in key for token in ('RANSWallEvolution', 'RANSEquilibriumWallEvolution',
            'RANSRobustVariant', 'RANSModels.', 'AnalyticalWall', 'EquilibriumAWF', 'RANSWallSolveDiagnostics'))}
    summary['epochs'][name] = entry
for label, first, second in (
    ('first_reproduction', 'original-reproduction-before', 'original-reproduction-after'),
    ('admissible_reproduction', 'admissible-reproduction-before', 'admissible-reproduction-after'),
    ('focused_checks', 'focused-before', 'focused-after'),
    ('controls', 'controls-before', 'controls-after')):
    a, b = (json.loads(paths[key].read_text()) for key in (first, second))
    compared = ('source_inputs', 'runtime_artifacts', 'runtime_bindings', 'python_packages', 'tool_images', 'runtime_environment')
    equal = {key: a[key] == b[key] for key in compared}
    assert all(equal.values()), (label, equal)
    summary.setdefault('before_after_comparisons', {})[label] = equal
for name, path in {
    'build-initial-receipt.json': r1/'build-001/receipt.json',
    'build-controlled-failed-receipt.json': r1/'controlled/build-001/receipt.json',
    'build-controlled-passed-receipt.json': r1/'controlled/build-002/receipt.json',
    'reproduction-initial-receipt.json': r1/'reproduction-001/receipt.json',
    'reproduction-admissible-receipt.json': r1/'controlled/reproduction-001/receipt.json',
    'focused-checks-receipt.json': r1/'controlled/focused-checks-001/receipt.json',
}.items():
    data = json.loads(path.read_text()); assert data['status'] != 'running'
    (out/name).write_bytes(path.read_bytes())
    summary['receipts'][name] = dict(original=str(path), sha256=sha(path), bytes=path.stat().st_size)
for name, path in {
    'build-initial.log': r1/'build-001/build.log',
    'build-controlled-failed.log': r1/'controlled/build-001/build.log',
    'build-controlled-passed.log': r1/'controlled/build-002/build.log',
    'focused-native-probes.log': r1/'controlled/focused-checks-001/native-probes.log',
}.items():
    assert path.stat().st_size < 2_000_000
    (out/name).write_bytes(path.read_bytes())
    summary['logs'][name] = dict(original=str(path), sha256=sha(path), bytes=path.stat().st_size)
(out/'runtime-summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
print(json.dumps(dict(output=str(out), manifest_epochs=len(paths),
    summary_sha256=sha(out/'runtime-summary.json'), files=len(list(out.iterdir())))))
