#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Hand completed declared H=.5 short runs to separate durable campaigns.

This controller never advances, repairs, or resumes a short run. A published
successful 160-step result, exact native checkpoint/history validation, and no
remaining prefix process are required before the one-time handoff. The existing
supervisor owns subsequent infrastructure recovery and numerical terminal rules.
Named companions require their own saved definition, reason and roots;
there is no automatic physical-model fallback from a failed faithful case.
"""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time

from continue_cls_epawf_pimple import CampaignLock, atomic, header, process_identity
from run_cls_accepted_awf import (FAMILY, FORMULATION, FUNCTIONAL_VERSION, INVENTORY,
    SPEC_SHA256, WORK, check_report, model_choices)
from run_cls_epawf_picard import THREADS, load, require, sha, state_hash

CASES = {'Accepted-AWF-05G': ('primary-05G', 'selected-05G'),
         'Accepted-AWF-05O': ('primary-05O', 'selected-05O')}
TERMINAL = {'numerically_rejected', 'targets_completed', 'infrastructure_interrupted'}
EXECUTED_SCRIPT_SHA256 = sha(Path(__file__))


def now():
    return datetime.now(timezone.utc).isoformat()


def argv(pid):
    try:
        return (Path('/proc')/str(int(pid))/'cmdline').read_bytes().rstrip(b'\0').decode().split('\0')
    except (FileNotFoundError, ProcessLookupError, PermissionError):
        return []


def live_commands(script, output):
    """Match a complete --output argument, rather than an unqualified PID."""
    found = []
    for p in Path('/proc').iterdir():
        if not p.name.isdigit():
            continue
        command = argv(p.name)
        if not any(Path(x).name == script for x in command) or '--output' not in command:
            continue
        i = command.index('--output')
        if i+1 < len(command) and Path(command[i+1]).resolve() == output:
            identity = process_identity(p.name)
            if identity:
                found.append(dict(identity=identity, command=command))
    return found


def prefix_children(base):
    live = []
    for packet in base.glob('*/process.json'):
        record = load(packet)
        actual = process_identity(record['pid'])
        # The old packet has no start token. A reused PID is unrelated unless
        # its full native command and working directory also match.
        if actual and actual['cwd'] == str(packet.parent.resolve()) and argv(record['pid']) == record['command']:
            live.append(dict(packet=str(packet), identity=actual, command=record['command']))
    return live


def validate_definition(definition, case, companion='faithful'):
    require(definition['family'] == FAMILY and definition['case'] == case and
        definition['formulation'] == FORMULATION and definition['wall_policy'] == FORMULATION,
        'different accepted-AWF family or case')
    require(definition['appendix_sha256'] == SPEC_SHA256 and definition['H'] == .5 and
        definition['z_K'] == definition['z_m'] == .5 and definition['h'] == .1,
        'different approved physical/reference configuration')
    require(definition['K_sampling'] == 'InterfaceTrace' and
        definition['K_functional_version'] == FUNCTIONAL_VERSION and
        definition['inventory'] == INVENTORY and definition['bound_solver'] == 'RetainedPool' and
        definition['work'] == WORK and definition['k_floor_relative'] == 1e-8,
        'different sampling, inventory, floor, or projection controls')
    require(companion in ('faithful','regularized','regularized-saturated') and definition['companion'] == companion and
        (not definition['companion_reason'] if companion == 'faithful' else bool(definition['companion_reason'].strip())) and
        definition['inner_refinement'] == 'none' and definition['perturbation']['amplitude'] == 0 and
        definition['model_choices'] == model_choices(companion),
        'the short run is not the declared unperturbed primary or named companion case')
    require(definition['dt'] == .000625 and definition['steps'] in (0, 160) and
        definition['eta'] == definition['prescribed_forcing'] == 0 and definition['relaxation'] == 1.,
        'different short-run temporal or source configuration')


def validate_completed(result, native, checkpoint, case, companion='faithful'):
    """Read-only checks: never manufacture an accepted report or native history."""
    require(result['status'] == 'numerical_trajectory_completed' and
        result['steps'] == native['steps'] == checkpoint['steps'] == 160 and
        abs(result['time']-.1) < 1e-13 and abs(native['time']-.1) < 1e-13 and
        abs(checkpoint['time']-.1) < 1e-13 and not result.get('failure') and
        not native.get('failure') and result['trajectory_receipt']['exit_code'] == 0,
        'short trajectory did not complete the accepted fixed horizon')
    validate_definition(result['definition'], case, companion)
    check_report(native, True, result['definition'])
    require(len(native['pimple_attempts']) == native['implicit_callback_count'] == 160 and
        all(t['numerically_accepted'] for t in native['pimple_attempts']),
        'short native prefix does not contain exactly 160 accepted scheduled passes')
    require(checkpoint['order'] == 1 and checkpoint['dt'] == .000625 and
        checkpoint['sha256'] == result['accepted_checkpoint']['sha256'], 'checkpoint identity differs')
    require(state_hash(checkpoint['coupled_state']) == state_hash(native['state']),
        'native checkpoint coupled state differs from report')
    require([state_hash(s) for s in checkpoint['history']] ==
        [state_hash(s) for s in native['native_history']] and
        state_hash(checkpoint['history_times']) == state_hash(native['native_history_times']),
        'native checkpoint histories differ from report')
    require(state_hash(native['accepted_states'][0]) == result['actual_h0_sha256'],
        'short run did not retain its own actual native h0')
    return dict(steps=160, time=checkpoint['time'], checkpoint_sha256=checkpoint['sha256'],
        state_sha256=state_hash(checkpoint['coupled_state']),
        native_history_sha256=[state_hash(s) for s in checkpoint['history']],
        native_history_times_sha256=state_hash(checkpoint['history_times']),
        actual_h0_sha256=result['actual_h0_sha256'],
        correction_integral=result.get('correction_integral'),
        physical_steps_verified=len(native['pimple_attempts']))


def copy_exact(source, target):
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        require(sha(source) == sha(target), 'existing handoff input differs: '+str(target))
    else:
        temporary = target.with_suffix(target.suffix+'.copying')
        shutil.copy2(source, temporary)
        require(sha(source) == sha(temporary), 'handoff copy hash mismatch')
        os.replace(temporary, target)
    return dict(source=str(source), path=str(target), sha256=sha(target))


def supervisor_command(solver, selected, case, inputs, targets, timeout, companion='faithful', reason=''):
    script = Path(__file__).with_name('supervise_cls_epawf_pimple.py')
    return [sys.executable, '-B', str(script), '--solver', str(solver), '--output', str(selected),
        '--wall-policy', FORMULATION, '--accepted-awf-case', case,
        '--accepted-awf-companion',companion,'--accepted-awf-companion-reason',reason,'--inventory', INVENTORY,
        '--bound-solver', 'RetainedPool', '--checkpoint', str(inputs/'checkpoint'),
        '--prefix-report', str(inputs/'prefix.evolution.json'), '--targets', *map(str, targets),
        '--block-steps', '8', '--timeout', str(timeout),
        '--bound-max-add', '2000', '--bound-max-total', '20000',
        '--bound-max-separations', '2000', '--bound-max-seconds', '900']


def inspect_or_launch(root, solver, solver_sha, case, targets, timeout, env, companion='faithful'):
    short_name, selected_name = CASES[case]
    if companion != 'faithful':
        short_name = short_name.replace('primary-', companion+'-')
        selected_name = selected_name.replace('selected-', 'selected-'+companion+'-')
    short = root/short_name
    selected = root/selected_name
    existing = live_commands('supervise_cls_epawf_pimple.py', selected)
    require(len(existing) <= 1, 'multiple selected supervisors; do not launch a replacement')
    if existing:
        return dict(status='successor_running', output=str(selected), process=existing[0])
    if (selected/'status.json').exists():
        state = load(selected/'status.json')
        if state['status'] in TERMINAL:
            return dict(status='successor_'+state['status'], output=str(selected), terminal=state)
    if (selected/'launch-intention.json').exists():
        return dict(status='handoff_infrastructure_interrupted', output=str(selected),
            reason='A launch intention already exists without a live supervisor; inspect recorded process and lock. No outer automatic retry.')
    base = short/case
    result_file = base/'result.json'
    if not result_file.exists():
        wrappers = live_commands('run_cls_accepted_awf.py', short)
        children = prefix_children(base) if base.exists() else []
        return dict(status='waiting_for_short_publication' if wrappers or children else 'short_infrastructure_interrupted',
            short=str(short), wrapper_processes=wrappers, native_or_probe_processes=children,
            reason='Only a fully published successful native result can initiate a handoff')
    result = load(result_file)
    if result['status'] != 'numerical_trajectory_completed':
        return dict(status='short_'+result['status'], result=str(result_file),
            failure=result.get('failure'), steps=result.get('steps'), time=result.get('time'))
    wrappers = live_commands('run_cls_accepted_awf.py', short)
    children = prefix_children(base)
    if wrappers or children:
        return dict(status='waiting_for_short_process_exit', wrapper_processes=wrappers,
            native_or_probe_processes=children)
    manifest = load(short/'manifest.json')
    require(manifest['status'] == 'suite_recorded' and manifest['solver']['sha256'] == solver_sha,
        'short wrapper completion/runtime differs')
    for path, expected in manifest['solver']['libraries'].items():
        require(sha(Path(path)) == expected, 'short linked-library identity changed: '+path)
    native_file = base/'trajectory/case.evolution.json'
    native = load(native_file)
    require(sha(native_file) == result.get('final_native_report_sha256', sha(native_file)), 'native report identity differs')
    cp = Path(native['checkpoint']).resolve()
    proof = validate_completed(result, native, header(cp, True), case, companion)
    # Validate the actual trajectory definition too: the result carries the
    # fresh initialization definition, while this one includes native restart.
    trajectory_definition=load(base/'trajectory/definition.json')
    validate_definition(trajectory_definition, case, companion)
    require(trajectory_definition['companion_reason'] == result['definition']['companion_reason'],
        'native initialization and trajectory companion reasons differ')
    require(sha(solver) == solver_sha, 'master frozen solver changed')
    inputs = selected/'inputs'
    copies = [copy_exact(cp/'State.bin', inputs/'checkpoint/State.bin'),
        copy_exact(native_file, inputs/'prefix.evolution.json'),
        copy_exact(result_file, inputs/'short-result.json'),
        copy_exact(short/'manifest.json', inputs/'short-manifest.json'),
        copy_exact(base/'initial/definition.json', inputs/'fresh-definition.json'),
        copy_exact(base/'initial/case.evolution.json', inputs/'native-initial.evolution.json')]
    command = supervisor_command(solver, selected, case, inputs, targets, timeout,
        companion, result['definition']['companion_reason'])
    intention = dict(schema=2, time=now(), case=case, family=FAMILY, command=command,
        companion=companion,companion_reason=result['definition']['companion_reason'],
        executed_controller_sha256=EXECUTED_SCRIPT_SHA256,
        validation=proof, immutable_inputs=copies, solver_sha256=solver_sha,
        thread_environment=THREADS, LD_LIBRARY_PATH=env['LD_LIBRARY_PATH'],
        controller_identity=process_identity(os.getpid()),
        source_files={str(p):sha(p) for p in (Path(__file__), Path(command[2]),
            Path(__file__).with_name('continue_cls_epawf_pimple.py'),
            Path(__file__).with_name('run_cls_accepted_awf.py'))})
    atomic(selected/'launch-intention.json', intention)
    with (selected/'supervisor.log').open('a') as stream:
        child = subprocess.Popen(command, cwd=selected, env=env, stdin=subprocess.DEVNULL,
            stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
    receipt = dict(intention, process_identity=process_identity(child.pid), pid=child.pid)
    atomic(selected/'launch-receipt.json', receipt)
    return dict(status='successor_launched', output=str(selected), process=receipt['process_identity'],
        checkpoint_sha256=proof['checkpoint_sha256'], prefix_steps=160, prefix_time=proof['time'])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--solver', type=Path, required=True)
    p.add_argument('--solver-sha256', required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--cases', choices=tuple(CASES), nargs='+', default=list(CASES))
    p.add_argument('--companion', choices=('faithful','regularized','regularized-saturated'), default='faithful')
    p.add_argument('--targets', type=int, nargs='+', default=[800,1600,3200,8000])
    p.add_argument('--timeout', type=float, default=7200)
    p.add_argument('--interval', type=float, default=15)
    p.add_argument('--max-wait', type=float, default=86400)
    a = p.parse_args()
    require(platform.system() == 'Linux', 'authoritative campaign orchestration is Linux-only')
    a.root = a.root.resolve();a.solver = a.solver.resolve();a.output = a.output.resolve()
    require(a.root.is_relative_to(Path('/media/jay/新加卷/nektar-v2026-runs')) and
        a.output.is_relative_to(a.root), 'isolated heavy experiment storage required')
    require(a.interval > 0 and a.max_wait > 0 and a.timeout > 0 and min(a.targets) > 160,
        'invalid bounded handoff controls')
    require(sha(a.solver) == a.solver_sha256, 'unexpected frozen master runtime')
    a.output.mkdir(parents=True, exist_ok=True);lock = CampaignLock(a.output)
    env = os.environ.copy();env.update(THREADS)
    env.update(LD_LIBRARY_PATH=str(a.solver.parent/'lib'), TMPDIR=str(a.root/'tmp'))
    os.environ.update(THREADS)
    started = time.monotonic();states = {}
    while True:
        for case in a.cases:
            if states.get(case,{}).get('status') in ('successor_launched','successor_running','handoff_validation_failed') or \
                    states.get(case,{}).get('status','').startswith(('short_','successor_numerically_','successor_targets_','handoff_infrastructure_')):
                continue
            try:
                states[case] = inspect_or_launch(a.root,a.solver,a.solver_sha256,case,a.targets,a.timeout,env,a.companion)
            except Exception as e:
                states[case] = dict(status='handoff_validation_failed',reason=str(e))
        waiting = any(s['status'].startswith('waiting_') for s in states.values())
        timed_out = time.monotonic()-started >= a.max_wait
        atomic(a.output/'status.json',dict(schema=2,companion=a.companion,
            executed_controller_sha256=EXECUTED_SCRIPT_SHA256,status='waiting' if waiting and not timed_out else
            ('wait_limit' if waiting else 'handoffs_recorded'),updated=now(),
            controller_identity=process_identity(os.getpid()),elapsed_seconds=time.monotonic()-started,cases=states))
        if not waiting or timed_out:
            return
        time.sleep(a.interval)


if __name__ == '__main__':
    main()
