#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Close CT16 receipts without converting missing qualification into a PASS.

Run after the final guarded collector has exited. This administrative exporter,
transfers and final publication consume the single existing publication reserve.
It cannot resume the allocation, launch a solver or change a resource ceiling.
"""
import argparse
import hashlib
import json
import subprocess
import time
from pathlib import Path

from ct_accepted15_collect import logical_size, max_address_space

START = 1794835825
CPU_RESERVE = 240
WALL_RESERVE = 600
BYTE_RESERVE = 256 * 1024**2


def read(path):
    return json.loads(path.read_text())


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024**2), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'snapshot', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    root, snapshot, output = (getattr(args, k).resolve()
                              for k in ('root', 'snapshot', 'output'))
    if snapshot.parent != root or output.parent != root or output.exists():
        raise ValueError('closure needs a sealed owned snapshot and new owned output')
    allocation = read(root/'allocation.json')
    if allocation['start_unix'] != START or allocation['limits'] != dict(
            soft_CPU=23760, soft_wall=28560, soft_live_bytes=15*1024**3):
        raise ValueError('original CT16 allocation identity changed')
    seal = read(snapshot/'seal.json')
    for name, digest in seal['files'].items():
        if Path(name).name != name or sha(snapshot/name) != digest:
            raise ValueError('snapshot seal mismatch: '+name)
    result = read(snapshot/'results.json')
    scientific = read(snapshot/'scientific-evidence.json')
    if result['status'] != 'UNRESOLVED' or result['physical_T0_407'] != 'UNATTEMPTED':
        raise ValueError('this exporter cannot establish admission or physical preservation')
    receipts, incomplete = {}, {}
    peak = written = 0
    for folder in sorted(root.glob('resources-*')):
        path = folder/'execution.json'
        if not path.is_file():
            raise ValueError('unclosed owned family: '+str(folder))
        receipt = read(path)
        receipts[folder.name] = dict(sha256=sha(path), **receipt)
        resource = receipt['resource']
        if any(resource.get(k) != 'pass' for k in ('aggregate_status', 'project_status')):
            incomplete[folder.name] = resource
        watch_path = folder/'package-budget-watch.json'
        if watch_path.is_file():
            sample = read(watch_path).get('last_sample', {})
            peak = max(peak, sample.get('peak_sampled_live_bytes', 0))
            written += sample.get('sampled_cumulative_process_wchar', 0)
    # Verify that every snapshot receipt is still the same closed receipt.
    old = result['engineering_ledger']['receipts']
    for name, receipt in old.items():
        if name not in receipts or receipt['sha256'] != receipts[name]['sha256']:
            raise ValueError('closed snapshot receipt changed: '+name)
    pending = Path(result['engineering_ledger']['current_collector_receipt_pending']).name
    if set(receipts) != set(old) | {pending}:
        raise ValueError('work was performed after the claimed final collector')
    output.mkdir()
    CPU = sum(q['process_cpu_seconds'] for q in receipts.values())
    elapsed = time.time()-START
    retained = logical_size(root)
    peak = max(peak, retained)
    maxima = {key: max((q['resource'].get(key, 0) for q in receipts.values()), default=0)
              for key in ('maximum_aggregate', 'maximum_project')}
    checks = dict(CPU=CPU+CPU_RESERVE <= 24000,
                  wall=elapsed+WALL_RESERVE <= 28800,
                  live_bytes=peak+BYTE_RESERVE <= 16*1024**3)
    if not all(checks.values()):
        raise ValueError('finite allocation plus existing publication reserve does not fit')
    accounting = dict(
        schema='CT16ClosedPublicationAccountingV1',
        status='closed_receipts_with_explicit_publication_reserve',
        outcome='UNRESOLVED', external_root=str(root),
        snapshot_sha256=sha(snapshot/'results.json'), snapshot_seal_sha256=sha(snapshot/'seal.json'),
        source_at_snapshot=result['source_revision'],
        source_at_export=subprocess.check_output(
            ['git', '-c', 'index.threads=1', 'rev-parse', 'HEAD'],
            cwd=root/'source', text=True).strip(),
        original_start_unix=START, closed_family_count=len(receipts),
        closed_family_CPU=CPU, publication_reserve=dict(
            CPU=CPU_RESERVE, wall=WALL_RESERVE, bytes=BYTE_RESERVE),
        CPU_upper_with_reserve=CPU+CPU_RESERVE,
        elapsed_wall_at_export=elapsed, wall_upper_with_reserve=elapsed+WALL_RESERVE,
        logical_live_bytes_at_export=retained, peak_observed_live_bytes=peak,
        new_live_upper_with_publication_reserve=peak+BYTE_RESERVE,
        sampled_cumulative_process_wchar_lower_bound=written,
        write_accounting='Sampled per-family process wchar highwaters include pipes and may miss short-lived writes; not exact physical disk I/O.',
        maximum_reported_workers=maxima, incomplete_or_failed_worker_audits=incomplete,
        unclosed_resources=[],
        maximum_recorded_process_AS_bytes=max_address_space(scientific),
        finite_CPU_wall_disk_checks=checks,
        all_receipts=receipts,
        reserve_scope='The one existing 240 CPU/600 wall/256 MiB publication reserve includes this exporter, source inspection, transfer duplication, Mac compaction, Git publication and final verification. No second administrative reserve is added.',
        separate_ledgers='Historical CT11-CT15 and prospective production excluded from new CT16 engineering; all new failures and conservative stop witnesses remain charged.',
        qualification_scope='Closed resource accounting does not complete G1-G5. Fresh candidate pair, selected controls, full workflow heldout, native retention cycles and complete admitted forecast remain missing.',
        physical_T0_407='UNATTEMPTED', full_transit='UNQUALIFIED')
    path = output/'closed-publication-accounting.json'
    path.write_text(json.dumps(accounting, indent=2, allow_nan=False)+'\n')
    (output/'seal.json').write_text(json.dumps(
        dict(files={path.name: sha(path)}), indent=2)+'\n')
    print(json.dumps({k: accounting[k] for k in (
        'outcome', 'closed_family_count', 'closed_family_CPU', 'CPU_upper_with_reserve',
        'wall_upper_with_reserve', 'new_live_upper_with_publication_reserve',
        'maximum_reported_workers', 'finite_CPU_wall_disk_checks')}))


if __name__ == '__main__':
    main()
