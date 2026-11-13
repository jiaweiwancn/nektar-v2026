#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Publish/verify the bounded H13 preparation; never launch a native process.

The original science-v2 producer remains immutable. A descriptive manifest view
adds the explicit model constants and flat-reference/terrain-geometry roles;
its bridge verifies every inherited value, rather than regenerating the case.
"""
import argparse
import ctypes
import json
import os
from pathlib import Path
import platform
import re
import shutil
import sys
import xml.etree.ElementTree as ET

from wj_h13_cases import model_metadata, save
from wj_h13_domain import identity, preflight, topology


def load(path):
    return json.loads(Path(path).read_text())


def check(record, path=None):
    actual = identity(record['path'] if path is None else path)
    assert (actual['bytes'], actual['sha256']) == (record['bytes'], record['sha256'])


def methods(root):
    result = {}
    for folder in sorted((root/'operations').iterdir()):
        receipt = folder/'receipt.json'
        log = folder/'command.log'
        if not receipt.exists() or not log.exists():
            continue
        for line in log.read_text().splitlines():
            match = re.match(r'^(test_\w+) \(test_wj_h13_domain\..*\) \.\.\. ok$', line)
            if match:
                result.setdefault(match[1], []).append(dict(operation=folder.name, receipt=identity(receipt), log=identity(log)))
    return result


def source_records(root):
    capsules = []
    for name in ('wj-h13-source-v1', 'wj-h13-observer-v3', 'wj-h13-observer-v4', 'wj-h13-observer-v5'):
        descriptor = root/(name+'.json')
        if not descriptor.exists():
            continue
        meta = load(descriptor)
        target = root/('source' if name.endswith('source-v1') else 'observer-'+name.rsplit('-', 1)[1])
        checked = []
        # v3 was transported but never unpacked or executed. Preserve that fact.
        if target.is_dir():
            for rec in meta['files']:
                path = target/rec['path']
                check(rec, path)
                checked.append(identity(path))
        capsules.append(dict(descriptor=identity(descriptor), archive=identity(root/(name+'.tar.gz')),
                             commit=meta['commit'], verified_unpacked_files=len(checked),
                             executed_or_imported=name != 'wj-h13-observer-v3', selected_files=checked if len(checked)<20 else []))
    return dict(native_compilation='NOT_RUN', starting_commit='36115115e7d4e0ce2e98d9ed83c4c28bc2e9ec89',
                initial_source_commit='e60e8b32f46b9a0d7fd4993006e9068be99b4ae0',
                science_v2_overlay_commit='65136787e7b31b3c3649acc358bc32ee9ef9b634',
                capsules=capsules, fixture_v2=dict(commit='ed22a0b3', file=identity(root/'observer-v2/test_wj_h13_domain.py')),
                publisher=identity(__file__),
                native_producer_and_runtime='not invoked; historical evidence only, no ELF/loader investigation')


def publish(root):
    out = root/'publication-v1'
    out.mkdir(exist_ok=False)
    science = root/'science-v2'
    shutil.copytree(science/'paired-case', out/'paired-case')
    for path in sorted(science.iterdir()):
        if path.is_file():
            shutil.copyfile(path, out/path.name)
    original = load(science/'paired-case/paired-domain.json')
    view = dict(original, model=model_metadata(), reference_scope=
                'flat analytical reference/calibration; terrain artifact is geometry only and has no initialized field or imported flat cutoff')
    assert all(view[key] == value for key, value in original.items())
    save(out/'paired-domain-model-view-v1.json', view)
    save(out/'model-view-bridge.json', dict(original=identity(science/'paired-case/paired-domain.json'),
        original_publication_copy=identity(out/'paired-case/paired-domain.json'),
        descriptive_view=identity(out/'paired-domain-model-view-v1.json'),
        all_original_fields_equal=True, added_keys=sorted(set(view)-set(original)),
        scientific_case_regenerated=False, numerical_or_input_changes=False))
    save(out/'source-identities.json', source_records(root))
    import numpy
    import scipy
    save(out/'offline-runtime.json', dict(system=platform.platform(), machine=platform.machine(), Python=sys.version,
        numpy=numpy.__version__, scipy=scipy.__version__, long_double_bytes=ctypes.sizeof(ctypes.c_longdouble),
        threads={k:os.environ[k] for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS')},
        native_runtime='NOT_USED; no native binaries loaded/launched or runtime closure requalified'))
    shutil.copyfile(root/'wj-h13-literature.json', out/'literature-identities.json')
    saved_methods = methods(root)
    assert len(saved_methods) == 9, saved_methods
    save(out/'focused-tests.json', dict(unique_passing_methods=9, methods=saved_methods,
        actual_execution='verbose unittest methods; CTest registered but neither configured nor executed',
        failures_preserved=True))
    receipts = out/'operations'
    receipts.mkdir()
    failures = []
    for folder in sorted((root/'operations').iterdir()):
        if not (folder/'receipt.json').exists():
            continue  # The current publication reservation is still active.
        record = load(folder/'receipt.json')
        shutil.copytree(folder, receipts/folder.name)
        if record['status'] != 'PASS':
            failures.append(dict(label=folder.name, receipt=identity(folder/'receipt.json'), log=identity(folder/'command.log')))
    save(out/'failure-lineage.json', dict(failed_operations=failures,
        causes={'focused-tests-v1': 'three synthetic fixtures lacked CONDITIONS; external-preflight fixture also found an invalid double-hyphen XML comment',
                'paired-offline-science-v1': 'same invalid double-hyphen comment; retained malformed XML, no native launch'},
        original_failed_case=identity(root/'science/paired-case/flat/case.xml'),
        repair='fixture nodes and well-formed comment; separate observer versions, receipts and science-v2 output',
        abandoned_observer_v3='transported only, not executed', native_invocations=0))
    # Snapshot before hashing. This is explicitly not the final closed ledger.
    save(out/'accounting-at-publication.json', load(root/'accounting.json'))
    save(out/'verification-scope.json', dict(
        Linux='fresh source capsule, selected original H10 XML/manifest, generated paired artifacts and receipt hashes',
        Mac='original literature/specification bytes; compact publication transport will be independently checked',
        historical_native_states='not re-read or requalified; no checkpoint/state import',
        aliases=[dict(original=identity(p), publication=identity(out/p.relative_to(science)))
                 for p in sorted(science.rglob('*')) if p.is_file()],
        no_compressed_science_payload=True))
    clocks = load(out/'paired-case/preservation-protocol.json')['stages']['B']['clocks']
    counts = original['counts']
    geom = load(out/'paired-geometry-checks.json')
    ref = load(out/'reference-boundary-checks.json')
    summary = dict(status=load(out/'component-gates.json'), counts=counts, storage=original['storage'],
        clocks=clocks, geometry={k:{n:v[n] for n in ('owned_sample_count','minimum_sampled_determinant','physical_shared_face_error','opposite_cofactor_error','top_lateral_coordinate_error','certificate')} for k,v in geom.items()},
        reference={k:v for k,v in ref.items() if k != 'outer_reference_rows'})
    save(out/'summary.json', summary)
    manifest = dict(schema='WJH13CompactPublicationV1',
        artifacts=[identity(p) for p in sorted(out.rglob('*')) if p.is_file()],
        immutable_science=[identity(p) for p in sorted(science.rglob('*')) if p.is_file()],
        source=identity(__file__), source_capsule=identity(root/'wj-h13-source-v1.tar.gz'),
        status='preparation PASS; target capacity OBSTRUCTED; native/duration/experimental stages NOT_RUN')
    save(out/'manifest.json', manifest)
    save(out/'seal.json', dict(manifest=identity(out/'manifest.json'), component_gates=identity(out/'component-gates.json')))
    print('H13 compact publication PASS; nine focused methods; native invocations zero', flush=True)


def verify(root):
    out = root/'publication-v1'
    manifest = load(out/'manifest.json')
    for rec in manifest['artifacts'] + manifest['immutable_science']:
        check(rec)
    for key in ('manifest','component_gates'):
        check(load(out/'seal.json')[key])
    original = load(out/'paired-case/paired-domain.json')
    view = load(out/'paired-domain-model-view-v1.json')
    bridge = load(out/'model-view-bridge.json')
    assert all(view[key] == value for key, value in original.items())
    for key in ('original','original_publication_copy','descriptive_view'):
        check(bridge[key])
    xml = out/'paired-case/flat/case.xml'
    guard = preflight(xml, original, root/'source')
    assert guard['status'] == 'OBSTRUCTED' and not guard['native_launch_allowed']
    terrain = ET.parse(out/'paired-case/terrain-chart-topology.xml').getroot()
    assert terrain.find('CONDITIONS') is None
    assert topology(terrain)['sha256'] == topology(xml)['sha256']
    assert topology(xml)['sha256'] == original['shared_topology']['sha256']
    for partner in original['geometry_partners'].values():
        check(partner)
    assert original['counts']['state'] == 24688
    assert original['counts']['free_solve_blocks']['epsilon'] == 4327
    gate = load(out/'component-gates.json')
    assert gate['all_native_invocations'] == gate['native_builds'] == gate['global_dense_matrix_allocations'] == 0
    assert gate['duration_preservation'] == gate['target_domain_evolution'] == 'NOT_RUN'
    assert len(methods(root)) == 9
    closure = root/'closure-v1'
    closure.mkdir(exist_ok=False)
    save(closure/'verification.json', dict(status='PASS', manifest=identity(out/'manifest.json'),
        source=identity(__file__), all_compact_and_original_science_bytes_verified=True,
        paired_topology_and_geometry_roles_verified=True, external_capacity_preflight=guard,
        unique_passing_test_methods=9, new_native_invocations=0, new_global_dense_matrix_allocations=0,
        decoded_science='no compressed science payload; source capsules separately archive- and member-hash verified',
        independent_Mac_transport_verification='pending after this immutable receipt'))
    print('H13 compact/semantic/source-role verification PASS', flush=True)


def close(root):
    """Run under the fixed closure charge after publication verification/QA.

    It does not use an execute() reservation that it would itself close while
    running. All completed operation receipts remain immutable and are copied.
    """
    from run_wj_h7 import update_journal
    from run_wj_h4 import size
    closure = root/'closure-v1'
    assert load(closure/'verification.json')['status'] == 'PASS'
    def finish(book):
        assert not book['closed']
        assert all(v['status'] != 'RESERVED' for v in book['entries'].values())
        ops = sum(v.get('seconds',v['reserved_seconds']) for v in book['entries'].values() if v['phase'] != 'closure')
        science = sum(v.get('seconds',v['reserved_seconds']) for v in book['entries'].values() if v['phase'] == 'science')
        closure_ops = sum(v.get('seconds',v['reserved_seconds']) for v in book['entries'].values() if v['phase'] == 'closure')
        assert ops <= 480 and science <= 360 and closure_ops <= 90 and 30+ops+90 <= 600
        book.update(closed=True, charged_seconds=30+ops+90, setup_charged_seconds=30., closure_charged_seconds=90.,
            launched_operation_seconds=ops, scientific_seconds=science, recorded_closure_operation_seconds=closure_ops,
            new_native_governing_attempts=0,new_native_CP0=0,new_native_probes=0,new_strict_readbacks=0,new_native_builds=0,
            new_restart_conversions=0,new_global_dense_matrix_allocations=0,
            output_bytes_before_final_closure=size(root), peak_recorded_tree_rss_bytes=max(v.get('peak_tree_rss_bytes',0) for v in book['entries'].values()),
            target_capacity='OBSTRUCTED', scientific_status='preparation PASS; preservation NOT_RUN',
            closure_scope='fully charged 90s covers publication, visual QA, transport, checks and finalization; no allowance transfer')
    update_journal(root/'accounting.json', finish)
    save(closure/'accounting-final.json', load(root/'accounting.json'))
    receipts = closure/'operation-receipts'
    receipts.mkdir()
    for folder in sorted((root/'operations').iterdir()):
        if (folder/'receipt.json').exists():
            shutil.copytree(folder, receipts/folder.name)
    save(closure/'manifest.json', dict(artifacts=[identity(p) for p in sorted(closure.rglob('*')) if p.is_file()],
        publication_manifest=identity(root/'publication-v1/manifest.json')))
    save(closure/'seal.json', dict(manifest=identity(closure/'manifest.json'), accounting=identity(closure/'accounting-final.json')))
    print(json.dumps({k:v for k,v in load(closure/'accounting-final.json').items() if k!='entries'}, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--mode', choices=('publish','verify','close'), required=True)
    args = parser.parse_args()
    globals()[args.mode](args.root)
