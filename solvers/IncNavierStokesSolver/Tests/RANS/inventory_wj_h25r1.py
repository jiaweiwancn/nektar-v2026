#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only takeover inventory. Never starts, signals or changes a science job."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


def pin(path):
    p = Path(path)
    return dict(path=str(p), bytes=p.stat().st_size,
                sha256=hashlib.sha256(p.read_bytes()).hexdigest())


def capture(root):
    root = Path(root)
    sys.path.insert(0, str(root/'source-v8/solvers/IncNavierStokesSolver/Tests/RANS'))
    from wj6_resources import scientific_workers
    ledger = json.loads((root/'accounting.json').read_text())
    phases = {k: sum(v.get('seconds', v['reserved_seconds'])
                    for v in ledger['entries'].values() if v['phase'] == k)
              for k in ('build', 'native', 'passive')}
    phases['total'] = sum(phases.values()) + ledger['setup_seconds'] + ledger['closure_seconds']
    operations = []
    for d in sorted((root/'operations').iterdir()):
        if d.is_dir():
            operations.append(dict(label=d.name, files=[pin(p) for p in sorted(d.iterdir())
                if p.name in ('request.json', 'receipt.json', 'attempt-counts.json')]))
    native = json.loads((root/'science/runtime-build-scalar-v1.json').read_text())
    inputs = json.loads((root/'inputs/baseline.json').read_text())
    checked = []

    def verify_tree(obj):
        if isinstance(obj, dict):
            if {'path', 'sha256', 'bytes'} <= obj.keys():
                actual = pin(obj['path'])
                checked.append(dict(expected=obj, actual=actual,
                                    status='PASS' if actual == {k: obj[k] for k in actual} else 'FAIL'))
            else:
                for v in obj.values(): verify_tree(v)
        elif isinstance(obj, list):
            for v in obj: verify_tree(v)

    verify_tree(native)
    verify_tree(inputs)
    runtime_deps = json.loads((root/'science/runtime-dependencies-v1.json').read_text())
    verify_tree(runtime_deps)
    source = Path('/media/jay/新加卷1/nektar-v2026')
    git = lambda *args: subprocess.check_output(['git', '-C', str(source), *args], text=True).strip()
    disk = shutil.disk_usage(root)
    return dict(schema='WJH25R1TakeoverV1', utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        hostname=os.uname().nodename, root=str(root), original_ledger=pin(root/'accounting.json'),
        original_account=ledger, original_usage=phases, allocation=scientific_workers(),
        operations=operations, unresolved_reservations=[k for k,v in ledger['entries'].items() if v['status']=='RESERVED'],
        identity_checks=checked, all_identity_checks_pass=all(v['status']=='PASS' for v in checked),
        general_checkout=dict(path=str(source), head=git('rev-parse', 'HEAD'), status=git('status','--short','--branch')),
        disk=dict(total=disk.total, used=disk.used, free=disk.free,
            filesystems=subprocess.check_output(['df','-PT',str(root),str(source)],text=True)),
        pending_paths={str(root/'science'/s): (root/'science'/s).exists()
            for s in ('h350','certificate-feasibility-v1.json','P-lookup-restart-v1')},
        inventory_pid=os.getpid(),
        h25_processes=[p for p in scientific_workers()['processes'] if str(root) in p['args']
            and p['pid'] != os.getpid() and 'inventory_wj_h25r1' not in p['args']])


if __name__ == '__main__':
    print(json.dumps(capture(sys.argv[1]), indent=2))
