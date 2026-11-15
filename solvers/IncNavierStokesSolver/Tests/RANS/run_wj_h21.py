#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H21 passive-only admission using the existing locked journal/tree helpers."""
import argparse
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
from run_wj_h4 import THREADS, load, save, sha, input_files
from run_wj_h7 import update_journal, complete_operation
from wj6_resources import scientific_workers, tree_usage, volume
from wj7_command import owned_groups

ROOT = 'wj-h21-outlet-saved-state-20261115'
FORBIDDEN = {'IncNavierStokesSolver', 'RANSWallResolutionProbe', 'WJKernelProbe',
             'cmake', 'make', 'ninja', 'mpirun', 'mpiexec'}


def passive_intent(argv):
    assert Path(argv[0]).name in ('python', 'python3'), 'passive Python only'
    assert not any(Path(str(v)).name.removesuffix('.elf') in FORBIDDEN for v in argv)
    assert not set(argv) & {'--governing', '--strict', '--cp0', '--restart', '--build', '--mass-solve'}
    assert '-c' not in argv, 'no arbitrary command bridge'
    if '-m' in argv:
        assert argv[argv.index('-m') + 1] == 'unittest'
    else:
        assert any(Path(str(v)).name == 'wj_h21_outlet.py' for v in argv)


def execute(root, label, argv, seconds, phase):
    root = Path(root)
    assert root.name == ROOT and phase in ('operation', 'science', 'closure')
    assert 0 < seconds <= (180 if phase == 'science' else 90)
    passive_intent(argv)
    allocation = scientific_workers()
    assert allocation['workers'] + 2 <= 72
    journal = root / 'accounting.json'
    folder = root / 'operations' / label
    folder.mkdir(exist_ok=False)

    def reserve(book):
        assert not book['closed'] and label not in book['entries']
        spent = lambda phases: sum(v.get('seconds', v['reserved_seconds'])
                                  for v in book['entries'].values() if v['phase'] in phases)
        if phase == 'closure':
            assert spent(('closure',)) + seconds <= 90
        else:
            assert spent(('operation', 'science')) + seconds <= 480
            if phase == 'science':
                assert spent(('science',)) + seconds <= 180
        book['entries'][label] = dict(status='RESERVED', phase=phase, reserved_seconds=seconds)

    update_journal(journal, reserve)
    request = dict(argv=list(map(str, argv)), threads=THREADS, phase=phase,
                   reserved_seconds=seconds, inputs=input_files(argv), allocation=allocation,
                   driver_sha256=sha(__file__))
    save(folder / 'request.json', request)
    env = dict(os.environ, **THREADS, TMPDIR=str(root / 'tmp'))

    def limits():
        import resource
        resource.setrlimit(resource.RLIMIT_AS, (2 * 1024**3, 2 * 1024**3))
        resource.setrlimit(resource.RLIMIT_FSIZE, (256 * 1024**2, 256 * 1024**2))

    start = time.monotonic()
    peak = 0
    groups = set()
    reason = None
    with (folder / 'command.log').open('x') as log:
        child = subprocess.Popen(['/usr/bin/time', '-v', *map(str, argv)], cwd=root,
                                 env=env, stdout=log, stderr=subprocess.STDOUT,
                                 start_new_session=True, preexec_fn=limits)
        while child.poll() is None:
            rss, _ = tree_usage(os.getpid())
            peak = max(peak, rss)
            groups |= owned_groups(child.pid)
            if rss > 2 * 1024**3 or time.monotonic() - start > seconds:
                reason = 'RESOURCE_LIMIT time/RSS'
            source_bytes = volume(root / 'source')
            if source_bytes > 128 * 1024**2 or volume(root) - source_bytes > 256 * 1024**2:
                reason = 'RESOURCE_LIMIT source/output'
            if reason:
                for group in groups:
                    try:
                        os.killpg(group, signal.SIGTERM)
                    except ProcessLookupError:
                        pass
                break
            time.sleep(.2)
        try:
            code = child.wait(timeout=5)
        except subprocess.TimeoutExpired:
            for group in groups:
                try:
                    os.killpg(group, signal.SIGKILL)
                except ProcessLookupError:
                    pass
            code = child.wait()
    source_bytes = volume(root / 'source')
    entry = dict(status='RESOURCE_LIMIT' if reason else 'PASS' if code == 0 else 'FAIL',
                 seconds=time.monotonic() - start, phase=phase, reserved_seconds=seconds,
                 exit=code, reason=reason, peak_tree_rss_bytes=peak,
                 source_bytes=source_bytes, new_output_bytes=volume(root) - source_bytes,
                 log_sha256=sha(folder / 'command.log'))
    complete_operation(journal, label, entry)
    if entry['status'] != 'PASS':
        update_journal(journal, lambda b: b.setdefault('first_failure', dict(label=label, **entry)))
    save(folder / 'receipt.json', dict(request, **entry))
    print(label, entry['status'], entry['seconds'], flush=True)
    return entry['status'] == 'PASS'


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--label', required=True)
    p.add_argument('--seconds', type=float, required=True)
    p.add_argument('--phase', choices=('operation', 'science', 'closure'), default='operation')
    p.add_argument('command', nargs=argparse.REMAINDER)
    a = p.parse_args()
    command = a.command[1:] if a.command[:1] == ['--'] else a.command
    sys.exit(0 if execute(a.root, a.label, command, a.seconds, a.phase) else 1)
