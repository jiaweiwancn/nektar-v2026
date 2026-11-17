#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""One bounded, fail-fast original-pulse standard recovery on level5."""
import argparse, os, subprocess, time
from pathlib import Path
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import require, sha, state_hash, THREADS
from wp2_budget import wall_limit


def main():
    parser = argparse.ArgumentParser()
    for key in ('root', 'solver', 'h1', 'preflight'):
        parser.add_argument('--' + key, type=Path, required=True)
    args = parser.parse_args()
    r = args.root
    record = r / 'confirmation-recovery-supervision.json'
    require(not record.exists(), 'confirmation recovery already has an owner')
    budget = load(r / 'resource-extension-24h.json')
    require(wall_limit(r) == 86400 and sha(args.solver) == budget['solver_sha256'],
            'unqualified confirmation runtime or resource scope')
    remaining = wall_limit(r) - (time.time() - load(r / 'origin.json')['started_epoch'])
    require(remaining > budget['forecast']['native_seconds'] + budget['forecast']['analysis_reserve_seconds'],
            'insufficient measured forecast headroom')
    scripts = Path(__file__).parent
    env = dict(os.environ, **THREADS, TMPDIR=str(r / 'tmp'))
    state = dict(version='WP2ConfirmationRecoveryControlV1', status='starting', pid=os.getpid(),
        source_sha256=sha(Path(__file__)), amendment_sha256=sha(r / 'resource-extension-24h.json'),
        solver_sha256=sha(args.solver), level=5, dt=.01, physical_endpoint=1.2, jobs={})
    atomic(record, state)

    def invoke(label, argv):
        state['status'] = label
        state['jobs'][label] = dict(argv=argv, started_epoch=time.time())
        atomic(record, state)
        timeout = wall_limit(r) - (time.time() - load(r / 'origin.json')['started_epoch']) - 120
        require(timeout > 0, 'cumulative resource envelope exhausted')
        with (r / ('confirmation-recovery-' + label + '.log')).open('x') as log:
            process = subprocess.run(argv, cwd=scripts, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=timeout)
        state['jobs'][label].update(exit_code=process.returncode, ended_epoch=time.time())
        atomic(record, state)
        require(process.returncode == 0, 'confirmation control failed; no later job launched: ' + label)

    base = ['/usr/bin/python3', str(scripts / 'run_wp2_preservation.py'), '--root', str(r),
        '--solver', str(args.solver), '--h1', str(args.h1), '--preflight', str(args.preflight),
        '--model', 'standard', '--level', '5', '--dt', '.01']
    seed = r / 'recovery-seed-confirmation-standard'
    output = r / 'recovery-confirmation-standard'
    try:
        invoke('seed', base + ['--output', str(seed), '--mode', 'seed'])
        initial = load(seed / 'progress.json')
        require(initial['status'] == 'completed_bounded_trajectory' and initial['accepted_steps'] == 0,
                'fresh confirmation seed did not pass its native audit')
        seed_report = load(Path(initial['observations'][0]['report']))
        require(seed_report['diagnostic_ivp']['maximum_actual_relative_k_departure'] >= .02,
                'confirmation pulse is below the predeclared amplitude')
        state['seed_sha256'] = state_hash(seed_report['state'])
        state['fixed_reference_sha256'] = state_hash(seed_report['fixed_reference'])
        state['fixed_defect_sha256'] = state_hash(seed_report['fixed_defect'])
        invoke('physical', base + ['--output', str(output), '--mode', 'recover', '--endpoint', '1.2',
                                  '--parent-progress', str(seed / 'progress.json')])
        progress = load(output / 'progress.json')
        require(progress['status'] == 'completed_bounded_trajectory' and progress['accepted_steps'] == 120 and
                progress['actual_time'] >= progress['Treq'], 'confirmation physical horizon incomplete')
        require(progress['observations'][0]['seed_sha256'] == state['seed_sha256'],
                'physical trajectory changed its fresh confirmation seed')
        terminal = r / 'confirmation-recovery-terminal-progress.json'
        terminal.write_bytes((output / 'progress.json').read_bytes())
        state['terminal_progress_sha256'] = sha(terminal)
        invoke('assessment', ['/usr/bin/python3', str(scripts / 'audit_wp2_recovery_leg.py'), '--root', str(r),
            '--model', 'standard', '--trajectory', str(terminal), '--output', str(r / 'confirmation-recovery-assessment')])
        invoke('window', ['/usr/bin/python3', str(scripts / 'audit_wp2_window_envelope.py'), '--root', str(r),
            '--models', 'standard', '--phases', 'recovery', '--trajectory', str(terminal),
            '--output', str(r / 'confirmation-recovery-window')])
        assessment = load(r / 'confirmation-recovery-assessment/result.json')
        state.update(status='completed_bounded_control', recovery_passed=assessment['passed'],
                     accepted_steps=120, actual_time=progress['actual_time'])
    except Exception as error:
        state.update(status='terminal_stop', failure=str(error))
        raise
    finally:
        state['ended_epoch'] = time.time()
        atomic(record, state)


if __name__ == '__main__':
    main()
