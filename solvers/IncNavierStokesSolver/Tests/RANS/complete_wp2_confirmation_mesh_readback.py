#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded read-only tail of the owned confirmation control; never restarts it."""
import argparse, os, subprocess, time
from pathlib import Path
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import require, sha, THREADS
from wp2_budget import wall_limit


def qualified_terminal(root, supervision):
    require(supervision.get('status') != 'terminal_stop', 'physical supervisor stopped unexpectedly')
    if not supervision.get('terminal_progress_sha256'):
        return None
    require(supervision['jobs']['physical']['exit_code'] == 0, 'physical job did not succeed')
    path = root / 'confirmation-recovery-terminal-progress.json'
    require(sha(path) == supervision['terminal_progress_sha256'], 'terminal snapshot changed')
    result = load(path)
    require(result['status'] == 'completed_bounded_trajectory' and result['accepted_steps'] == 120 and
            result['actual_time'] >= result['Treq'], 'terminal physical horizon is incomplete')
    return path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--supervisor-pid', type=int, required=True)
    args = parser.parse_args()
    r = args.root
    receipt = r / 'confirmation-recovery-mesh-postprocess.json'
    require(not receipt.exists(), 'mesh readback already has an owner')
    plan = r / 'confirmation-recovery-mesh-comparison-command.json'
    argv = load(plan)
    program = Path(argv[1])
    require(program.name == 'analyze_wp2_recovery_mesh_control.py' and
            program.parent == Path(__file__).parent, 'readback plan is not the authorized analysis')
    command_hash, program_hash = sha(plan), sha(program)
    state = dict(version='WP2ConfirmationMeshReadbackTailV1', pid=os.getpid(), status='waiting',
        physical_supervisor_pid=args.supervisor_pid, source_sha256=sha(Path(__file__)),
        command_sha256=command_hash, program_sha256=program_hash, started_epoch=time.time(), governing_solves=0)
    atomic(receipt, state)
    try:
        while True:
            remaining = wall_limit(r) - (time.time()-load(r/'origin.json')['started_epoch'])
            require(remaining > 120, 'cumulative readback envelope exhausted')
            supervision = load(r / 'confirmation-recovery-supervision.json')
            require(supervision['pid'] == args.supervisor_pid, 'physical supervisor ownership changed')
            terminal = qualified_terminal(r, supervision)
            if terminal is not None:
                break
            command_line = (Path('/proc') / str(args.supervisor_pid) / 'cmdline').read_bytes()
            require(b'run_wp2_confirmation_recovery.py' in command_line,
                    'physical controller missing; no restart or later job is allowed')
            time.sleep(30)
        require(sha(plan) == command_hash and sha(program) == program_hash, 'readback source or plan changed')
        state.update(status='analyzing', terminal_progress_sha256=sha(terminal))
        atomic(receipt, state)
        from run_cao_tamura_budget_wp2 import charge
        charge(r, 'matched confirmation recovery mesh readback', readback_comparisons=1)
        with (r/'confirmation-recovery-mesh-comparison.log').open('x') as log:
            p = subprocess.run(argv, cwd=program.parent, env=dict(os.environ, **THREADS),
                stdout=log, stderr=subprocess.STDOUT, timeout=min(600, remaining-120))
        require(p.returncode == 0, 'matched mesh analysis failed; no subsequent numerical job')
        result = r/'confirmation-recovery-mesh-comparison.json'
        state.update(status='completed_readback', exit_code=p.returncode, result_sha256=sha(result),
                     quarter_budget_passed=load(result)['passed'])
    except Exception as error:
        state.update(status='terminal_stop', failure=str(error))
        raise
    finally:
        state['ended_epoch'] = time.time()
        atomic(receipt, state)


if __name__ == '__main__':
    main()
