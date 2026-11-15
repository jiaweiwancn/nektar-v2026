#!/usr/bin/env python3
"""Close pre-publication accounting without changing any evaluated CT12 seal.

Run under ct_accepted3_command.py. The launcher's own execution.json is written
after this program exits and must accompany the exported supplement. Final
arithmetic adds that closed receipt; it never rewrites evaluated-v1.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile


def read(path):
    return json.loads(path.read_text())


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--self-control', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    a.output.mkdir(exist_ok=False)
    ev = a.root / 'evaluated-v1'
    sealed = read(ev / 'completed-artifacts.json')
    assert all(sha(ev / k) == v for k, v in sealed['files'].items())
    families = {}
    for d in sorted(a.root.glob('resources-*')):
        if d == a.self_control:
            continue
        e = read(d / 'execution.json')  # absence is an unclosed-family failure
        assert e['resource']['project_status'] == 'pass'
        assert e['resource']['aggregate_status'] == 'pass'
        high = {}
        for f in d.iterdir():
            if f.name in ('project-highwater.jsonl', 'worker-guard-highwater.jsonl',
                          'worker-guard-events.jsonl'):
                high[f.name] = [json.loads(x) for x in f.read_text().splitlines()]
        families[d.name] = dict(execution=e, sha256=sha(d / 'execution.json'),
                                highwater=high)
    roots = {v['execution'][key] for v in families.values()
             for key in ('control_parent_pid', 'guard_pid', 'owned_family_pid')}
    union_peak = 0
    witness = None
    for v in families.values():
        for row in v['highwater'].get('worker-guard-highwater.jsonl', []):
            owned = set(roots)
            while True:
                add = {x['pid'] for x in row['processes'] if x['parent'] in owned} - owned
                if not add:
                    break
                owned.update(add)
            n = sum(x['threads'] for x in row['processes'] if x['pid'] in owned)
            if n > union_peak:
                union_peak = n
                witness = dict(time_unix=row['time_unix'], threads=n,
                               processes=[x for x in row['processes'] if x['pid'] in owned])
    cpu = sum(x['execution']['process_cpu_seconds'] for x in families.values())
    wall = sum(max(0., x['execution']['wall_seconds'] -
                   x['execution']['allocation_wait_seconds']) for x in families.values())
    aggregate = max(x['execution']['resource']['maximum_aggregate'] for x in families.values())
    assert union_peak <= 16 and aggregate <= 72
    # Reserve final packaging/receipt work. Its closed receipt must fit this
    # bound; this is not a substitute for its measured process CPU afterward.
    assert cpu + 60 < 12000 and wall + 120 < 14400
    output_bytes = sum(f.stat().st_size for f in a.root.rglob('*') if f.is_file()
                       and f.relative_to(a.root).parts[0] not in ('source', 'tmp'))
    assert output_bytes + 64 * 1024**2 < 12 * 1024**3
    write(a.output / 'prepublication-accounting.json', dict(
        schema='CT12ClosedPublicationSupplementV1', status='pass',
        all_previous_families_closed=True, family_count=len(families),
        process_CPU_seconds=cpu, active_wall_upper_seconds=wall,
        output_bytes_before_packaging=output_bytes, maximum_combined_project=union_peak,
        maximum_aggregate=aggregate, combined_project_witness=witness,
        publication_CPU_reserve=60, publication_wall_reserve=120,
        scope='add the final closed publication execution receipt; elapsed minus allocation wait is a conservative active-wall upper bound, without subtracting partially resumed intervals',
        families=families))
    source = a.root / 'source'
    revision = subprocess.check_output(['git', '-c', 'core.preloadIndex=false',
                                        'rev-parse', 'HEAD'], cwd=source, text=True).strip()
    dirty = subprocess.check_output(['git', '-c', 'core.preloadIndex=false',
                                     'status', '--porcelain'], cwd=source, text=True)
    assert not dirty
    provenance = read(ev / 'provenance.json')
    for name, digest in provenance['unchanged_CPP'].items():
        assert sha(source / name) == digest
    manifest = source / 'docs/plans/cao-tamura-accepted-8-evidence/physical-reference-manifest.json'
    assert sha(manifest) == '004bc0fc82054296a68d56de77751fc1f0a41570ee16b0a174ae8087f262d296'
    runs = a.root.parent
    runtime_files = {
        'governing': runs / 'ct-accepted-9-20261114/build-v8/runtime/IncNavierStokesSolver',
        'probe': runs / 'ct-accepted-11-20261115/build-v2/runtime/RANSPreservationProbe'}
    expected = {'governing': '5830a2cd8a4b0b5b12095555171eb28bc0f63e1daa97652d63da244c3fd808c7',
                'probe': '20236b14eb4a90dea0f2a9a38b92e41135a83e9241c978fb265e185098f586e8'}
    runtimes = {k: dict(path=str(f), sha256=sha(f)) for k, f in runtime_files.items()}
    assert all(runtimes[k]['sha256'] == expected[k] for k in expected)
    mount = subprocess.check_output(['findmnt', '-T', str(a.root), '-n', '-o',
                                    'TARGET,SOURCE,FSTYPE,OPTIONS'], text=True).strip()
    assert '/media/jay/新加卷1 ' in mount
    write(a.output / 'publication-provenance.json', dict(
        schema='CT12PublicationProvenanceV1', source_revision=revision, source_clean=True,
        archive=str(a.root), mount=mount, runtime=runtimes,
        source_provenance_sha256=sha(ev / 'provenance.json'),
        evaluated_seal_sha256=sha(ev / 'completed-artifacts.json'),
        physical_manifest_sha256=sha(manifest),
        driver_sha256=sha(Path(__file__)),
        threads={k: os.environ.get(k) for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')},
        historical='CT4–CT11 files, seals, runtime installations and failure classifications unchanged',
        publication='exact evaluated bodies exported; final closed launcher receipt and additive arithmetic supplement kept separately'))
    external = []
    for name in ('pressure-v5/pressure-vectors.json', 'axes-v2/pressure-vectors.json',
                 'retention-v1/index.json', 'retention-v1/required-proof.json',
                 'face-v1/physical-face-samples.json', 'static-v1/results.json',
                 'pressure-v4/request-manifest.json', 'profile-v3/request-manifest.json'):
        f = a.root / name
        if f.is_file():
            external.append(dict(path=str(f), bytes=f.stat().st_size, sha256=sha(f)))
    write(a.output / 'external-artifact-index.json', dict(
        schema='CT12SelectedExternalArtifactsV1', files=external,
        scope='selected additional heavy witness identities; complete native requests/state/history/CSR identities are bound in each evaluated body, not duplicated here'))
    shutil.copy2(Path(__file__), a.output / 'publication-driver.py')
    archive = a.output / 'ct12-compact-v1.tar.gz'
    with tarfile.open(archive, 'w:gz', compresslevel=6) as t:
        t.add(ev, arcname='evaluated-v1')
        for f in sorted(a.output.iterdir()):
            if f != archive:
                t.add(f, arcname=f.name)
    write(a.output / 'transfer-sha256.json', dict(file=archive.name, bytes=archive.stat().st_size,
                                                sha256=sha(archive)))
    print(json.dumps(dict(status='pass', source_revision=revision, CPU_before=cpu,
                         wall_upper_before=wall, maximum_project=union_peak,
                         maximum_aggregate=aggregate, transfer=read(a.output / 'transfer-sha256.json'))), flush=True)


if __name__ == '__main__':
    main()
