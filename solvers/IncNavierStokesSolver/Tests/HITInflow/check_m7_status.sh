#!/usr/bin/env bash
# Read-only Milestone 7 status over Tailscale. Optional argument: refresh seconds.
set -eu
interval="${1:-0}"
case "$interval" in ''|*[!0-9]*) echo 'Usage: check_m7_status.sh [refresh_seconds; 0 = once]' >&2; exit 2;; esac
exec ssh -o BatchMode=yes -o ConnectTimeout=15 -o ServerAliveInterval=30 \
  jay@100.105.247.110 "python3 - $interval" <<'PY'
from pathlib import Path
import datetime
import json
import math
import os
import sys
import time

ROOT = Path('/media/jay/新加卷/nektar-v2026-runs/hit-inflow/milestone-7/97014ab')
INTERVAL = int(sys.argv[1])
if hasattr(os, 'sched_getaffinity') and 52 in os.sched_getaffinity(0):
    os.sched_setaffinity(0, {52})


def read(relative):
    try:
        return json.loads((ROOT / relative).read_text())
    except (FileNotFoundError, json.JSONDecodeError):
        return {}


def report():
    if not ROOT.is_dir():
        raise SystemExit('Run directory unavailable: ' + str(ROOT))
    now = datetime.datetime.now(datetime.timezone.utc)
    progress = read('live-progress.json')
    rows = {c['id']: c for c in progress.get('cases', [])}
    recipe = read('cbc-frozen-study.json').get('baseline', {})
    end, dt = recipe.get('duration', 1.6), recipe.get('dt', .000025)
    completed = 0
    failed = []
    print('\n' + now.strftime('%Y-%m-%d %H:%M:%S UTC') + ' — HIT Milestone 7', flush=True)
    print(f'{"Case":<27} {"State":<15} {"Progress":>9} {"Simulated time":>19} {"Step":>14}')
    for family, prefix in [('cbc-production', 'cbc'), ('vk-production', 'vk')]:
        study = read(family + '/study.receipt.json')
        study_cases = {c['id']: c for c in study.get('cases', [])}
        for seed in [20261028, 4315228324]:
            name = prefix + '-seed-' + str(seed)
            directory = family + '/' + name
            run = read(directory + '/run.receipt.json')
            solver = read(directory + '/solver.receipt.json')
            initialization = read(directory + '/initialization.receipt.json')
            state, physical = 'queued', None
            if run.get('solver', {}).get('exit_code') == 0:
                state, physical = 'finished', end
                completed += 1
            elif ((solver and (solver.get('exit_code') != 0 or solver.get('native_iteration_limit_exceeded')))
                  or (initialization and initialization.get('exit_code') != 0)
                  or study_cases.get(name, {}).get('status') == 'failed'):
                state = 'FAILED'
                failed.append(name)
            elif name in rows:
                physical = rows[name].get('last_history_time_s')
                state = 'running' if physical is not None else 'initializing'
            elif (ROOT / directory).is_dir():
                state = 'initializing'
            valid = physical is not None and math.isfinite(physical) and 0 <= physical <= end + 2*dt
            percent = f'{min(100, 100*physical/end):.1f}%' if valid else '--'
            clock = f'{physical:.4f} / {end:g} s' if valid else '--'
            step = f'{round(physical/dt)} / {round(end/dt)}' if valid else '--'
            print(f'{name:<27} {state:<15} {percent:>9} {clock:>19} {step:>14}')
    vk = read('vk-queue-receipt.json').get('status', 'not started')
    post = read('postprocessing-da607a9/receipt.json')
    figures = read('production-figures-receipt.json')
    analyzed = sum(c.get('status') == 'completed' for c in post.get('cases', []))
    print(f'Von Karman queue: {vk}')
    print(f'Postprocessing: {post.get("status", "not started")} ({analyzed}/4 cases)')
    print(f'Figures: {figures.get("status", "not started")}')
    print(f'Progress snapshot: {progress.get("utc", "not yet available")}')
    if progress.get('affinity_violation') or progress.get('worker_warning'):
        print('ATTENTION: the job monitor reports a CPU allocation flag.')
    if failed or vk in ('failed', 'blocked_by_CBC_failure') or post.get('status') in ('failed', 'blocked_by_pipeline_check') or figures.get('status') in ('failed', 'blocked_by_postprocessing'):
        print('PIPELINE NEEDS ATTENTION. Retain the logs and report this status.')
        return 'failed'
    if completed == 4:
        print('ALL FOUR SOLVER RUNS FINISHED.')
        if post.get('status') == 'completed' and figures.get('status') == 'completed':
            print('READY FOR MILESTONE 7 REVIEW.')
            return 'ready'
    sys.stdout.flush()
    return 'running'


try:
    while True:
        status = report()
        if not INTERVAL or status != 'running':
            raise SystemExit(1 if status == 'failed' else 0)
        time.sleep(INTERVAL)
except KeyboardInterrupt:
    pass
PY
