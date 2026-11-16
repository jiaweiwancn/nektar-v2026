#!/usr/bin/env python3
"""Release a verified owned HIT resource hold after existing guards finish.

This bounded scheduler never signals another experiment. Sleeping MPI service
threads cause older RANS guards to suspend themselves, so defer HIT until those
guards exit. The saved PID/start-time identities prevent signaling reused PIDs.
"""
import argparse
import json
import os
from pathlib import Path
import signal
import sys
import time


def process(pid):
    try:
        fields = (Path('/proc') / str(pid) / 'stat').read_text().rsplit(')', 1)[1].split()
        return {"state": fields[0], "process_start_ticks": fields[19]}
    except (FileNotFoundError, ProcessLookupError):
        return None


def same_process(item):
    current = process(item["pid"])
    return current is not None and current["process_start_ticks"] == item["process_start_ticks"]


def guards():
    found = []
    for path in Path('/proc').iterdir():
        if not path.name.isdigit():
            continue
        try:
            name = (path / 'comm').read_text().strip()
            if name.startswith('python') and b'guard_rans_worker_budget.py' in (path / 'cmdline').read_bytes():
                found.append(int(path.name))
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            pass
    return found


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hold-file', type=Path, required=True)
    parser.add_argument('--source-root', type=Path, required=True)
    parser.add_argument('--max-wait-seconds', type=float, default=4*3600)
    args = parser.parse_args()
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(args.source_root / 'solvers/IncNavierStokesSolver/Tests/HITInflow/pilots'))
    import run_band100 as run
    hold = json.loads(args.hold_file.read_text())
    queue_path = args.hold_file.parent / 'controls-01/queue.json'
    queue = json.loads(queue_path.read_text())
    run.require(hold['schema'] == 'NEKTAR_HIT_RESOURCE_HOLD_V1'
                and hold['queue_pid'] == queue['pid'] and hold.get('all_owned_processes_held'),
                'unverified or different resource hold')
    run.require(0 < args.max_wait_seconds <= 4*3600, 'resource wait must be bounded to at most4h')
    run.require(same_process(hold), 'study PID identity no longer matches')
    command = (Path('/proc') / str(queue['pid']) / 'cmdline').read_bytes()
    run.require(b'queue_band100_controls.py' in command and str(queue_path.parent).encode() in command,
                'queue PID command does not match this experiment')
    hold['queue_process_start_ticks'] = process(queue['pid'])['process_start_ticks']
    hold.update(helper_pid=os.getpid(), helper_sha256=run.sha(__file__),
                maximum_wait_seconds=args.max_wait_seconds,
                hard_wait_timeout_action='Kill only the saved owned HIT family; retain all partial evidence',
                helper_started_utc=run.utc())
    cpu = queue['allocation']['cpus'][-1]
    os.sched_setaffinity(0, {cpu})
    quiet_since = None
    while True:
        run.require(same_process(hold), 'held study exited or its PID identity changed')
        run.require(same_process(dict(pid=queue['pid'], process_start_ticks=hold['queue_process_start_ticks'])),
                    'queue exited or its PID identity changed')
        existing = guards()
        budget = run.external_worker_budget(queue['pid'], len(queue['allocation']['cpus']))
        quiet_since = (quiet_since or time.monotonic()) if not existing and budget['passed'] else None
        hold.update(status='WAITING_FOR_EXISTING_GUARDS', active_guard_pids=existing,
                    worker_budget=budget, updated_utc=run.utc(),
                    elapsed_wait_seconds=time.time()-hold['started_unix'])
        run.atomic_json(args.hold_file, hold)
        if time.time()-hold['started_unix'] >= args.max_wait_seconds:
            killed = []
            for item in hold['held_processes']:
                if same_process(item):
                    os.kill(item['pid'], signal.SIGKILL)
                    killed.append(item['pid'])
            hold.update(status='HARD_RESOURCE_WAIT_TIMEOUT', killed_owned_pids=killed, finished_utc=run.utc())
            run.atomic_json(args.hold_file, hold)
            return 124
        if quiet_since is not None and time.monotonic()-quiet_since >= 30:
            run.verify_files(queue['bindings'])
            if guards():
                quiet_since = None
                continue
            resumed = []
            # Restore native children before their waiting controller.
            items = [p for p in hold['held_processes'] if p['pid'] != hold['pid']] + [hold]
            for item in items:
                if same_process(item):
                    os.kill(item['pid'], signal.SIGCONT)
                    resumed.append(item['pid'])
            hold.update(status='RESUMED', resumed_owned_pids=resumed, finished_utc=run.utc())
            run.atomic_json(args.hold_file, hold)
            print('HIT resumed after existing guards exited; native binaries/inputs unchanged', flush=True)
            return 0
        time.sleep(5)


if __name__ == '__main__':
    raise SystemExit(main())
