#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Durable scientific stop observer for an already authorized continuation.

Does not launch/retry native steps. Existing checker and full-epsilon metric
own the unchanged scientific screen. A stop request acts only after a native
block; checkpoints and the production acceptance path are never rewritten.
"""
import argparse
from datetime import datetime, timezone
import os
from pathlib import Path
import subprocess
import sys
import time

from continue_cls_epawf_pimple import CampaignLock, atomic, load_progress, process_identity
from run_cls_epawf_picard import THREADS, load, require, sha


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('campaign', 'policy', 'probe', 'output'): p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--maximum-step', type=int, default=9600); p.add_argument('--interval', type=float, default=30)
    a = p.parse_args(); os.environ.update(THREADS); a.output.mkdir(parents=True, exist_ok=False); lock = CampaignLock(a.output)
    policy = load(a.policy); require(policy['primary_physical_window'] == .5 and policy['required_consecutive_windows'] == 3, 'unchanged steady policy required')
    out = dict(status='observing', started_UTC=datetime.now(timezone.utc).isoformat(), script_sha256=sha(Path(__file__)),
        observer_identity=process_identity(os.getpid()), campaign=str(a.campaign), maximum_step=a.maximum_step,
        probe_sha256=sha(a.probe), policy_sha256=sha(a.policy), assessed_steps=[], diagnostics={},
        stopping_policy='First three passing unchanged half-unit windows, INCLUDING full transported epsilon, or maximum authorized endpoint. Stop only at a completed native block.')
    atomic(a.output/'progress.json', out)
    started=time.monotonic()
    while True:
        if not (a.campaign/'progress.json').exists():
            out.update(status='waiting_for_native_progress',updated_UTC=datetime.now(timezone.utc).isoformat())
            if time.monotonic()-started>300:
                out.update(status='infrastructure_interrupted',failure='native progress not published within startup allowance')
                atomic(a.output/'progress.json',out);return
            atomic(a.output/'progress.json',out);time.sleep(2);continue
        progress = load_progress(a.campaign/'progress.json')
        require(max(progress['requested_targets'])<=a.maximum_step,'supervised continuation exceeds the authorized endpoint')
        out.update(accepted_steps=progress['accepted_steps'], accepted_time=progress['accepted_time'], checkpoint=progress['checkpoint'],
            checkpoint_sha256=progress.get('checkpoint_sha256'), native_identity=progress.get('native_identity'),
            updated_UTC=datetime.now(timezone.utc).isoformat())
        if progress['status'] in ('targets_completed', 'numerically_rejected', 'stopped_at_checkpoint'):
            out.update(status=progress['status'], failure=progress.get('failure')); atomic(a.output/'progress.json', out); return
        stride = round(.5/progress['dt']); segments = {int(v['end_step']): v for v in progress['segments'] if not v['failure'] and v['exit_code'] == 0}
        eligible = [n for n in segments if n%stride == 0 and n not in out['assessed_steps'] and all(n-i*stride in segments for i in range(4))]
        if eligible:
            n = max(eligible); base = a.output/f'step-{n:06d}'; base.mkdir()
            common = base/'common'; epsilon = base/'full-epsilon'
            commands = [[sys.executable, '-B', str(Path(__file__).with_name('check_cls_awf_steady_windows.py')),
                '--campaign', 'Refined05='+str(a.campaign), '--policy', str(a.policy), '--output', str(common)],
                [sys.executable, '-B', str(Path(__file__).with_name('confirm_cls_awf_outer_epsilon.py')),
                '--common-panel', str(common), '--case', 'Refined05', '--probe', str(a.probe), '--output', str(epsilon)]]
            results = []
            with (base/'observer.log').open('w') as log:
                for cmd in commands:
                    code = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, timeout=1800).returncode
                    results.append(code)
                    if code: break
            out['assessed_steps'].append(n)
            entry = dict(commands=commands, exit_codes=results, folder=str(base))
            if results == [0, 0]:
                screen = load(epsilon/'summary.json'); entry.update(passing_windows=screen['passing_windows'],
                    summary_sha256=sha(epsilon/'summary.json'))
                checked=load(common/'summary.json')['cases']['Refined05']['selected_steps'][-1]
                entry['actually_assessed_endpoint_step']=checked
                if screen['passing_windows'] == 3:
                    stop = dict(reason='three_complete_steady_windows_with_full_outer_epsilon', qualified_global_step=checked,
                        scientific_summary=str(epsilon/'summary.json'), scientific_summary_sha256=sha(epsilon/'summary.json'),
                        policy_sha256=sha(a.policy), observer_identity=out['observer_identity'])
                    require(not (a.campaign/'stop-after-block.json').exists(), 'a separate stop request already exists')
                    atomic(a.campaign/'stop-after-block.json', stop); out['status'] = 'scientific_stop_requested'
            else: entry['status'] = 'diagnostic_failed_production_unchanged'
            out['diagnostics'][str(n)] = entry
        atomic(a.output/'progress.json', out)
        time.sleep(a.interval)


if __name__ == '__main__': main()
