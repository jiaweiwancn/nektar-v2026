#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Continue a completed H=.5 accepted-AWF campaign from its actual T=5 state.

The completed parent is preserved. A separately locked successor keeps the
same frozen executable/library closure, coupled state and native IMEX1 history.
An isolated restart step is verification-only; the production first checkpoint
must equal it. Scientific steady-window thresholds are not acceptance gates.
"""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from benchmark_cls_awf_evaluation_reuse import replay
from continue_cls_epawf_pimple import (CampaignLock, atomic, assert_native_stopped,
    header, prefix_budget, process_identity)
from diagnose_cls_accepted_covariance import checkpoint
from follow_cls_accepted_awf_pair import live_commands
from handoff_cls_awf_runtime import launch_validated_successor, matches
from run_cls_accepted_awf import WORK, check_report, model_choices
from run_cls_epawf_picard import THREADS, load, require, sha, state_hash


def validate_extension(progress, native, saved, targets):
    """Check the real complete prefix, not a manufactured initial-value state."""
    require(progress['status'] == 'targets_completed' and not progress.get('failure') and
        not native.get('failure'), 'only a successfully completed parent can be extended')
    require(progress['accepted_steps'] == native['steps'] == saved['steps'] == 8000 and
        abs(saved['time']-5.) < 1e-10 and
        progress['accepted_time'] == native['time'] == saved['time'], 'parent is not the actual T=5 endpoint')
    require(saved['order'] == native['order'] == 1 and
        saved['dt'] == progress['dt'] == native['dt'] == .000625, 'timestep/history scheme differs')
    require(targets == [9600], 'extension target must be T=6')
    require(progress['wall_policy'] == 'CLSCalibratedBudgetAWF' and
        progress['case_family'] == 'AcceptedAWFV1' and progress['accepted_awf_case'] == 'Accepted-AWF-05G',
        'extension is restricted to the declared H=.5 graded budget family')
    definition = progress['case_definition']; companion = progress['accepted_awf_companion']
    require(companion in ('regularized', 'regularized-saturated') and
        definition['model_choices'] == model_choices(companion) and
        definition['companion'] == companion and
        definition['companion_reason'] == progress['accepted_awf_companion_reason'], 'companion identity differs')
    require(definition['H'] == definition['z_K'] == definition['z_m'] == .5 and
        definition['K_sampling'] == 'InterfaceTrace' and
        definition['inner_refinement'] == progress['accepted_awf_inner_refinement'] == 'none' and
        definition.get('outer_refinement', 'none') == 'none' and
        definition['dt'] == .000625 and definition['k_floor_relative'] == 1e-8 and
        definition['eta'] == definition['prescribed_forcing'] == 0 and definition['relaxation'] == 1,
        'physical or numerical method configuration differs')
    require(progress['projection_work_controls'] == definition['work'] == WORK and
        progress['inventory'] == 'AcceptedBudgetV1' and progress['bound_solver'] == 'RetainedPool' and
        progress['evaluation_reuse'] == 'ExactStepV1' and
        progress['audit_revision'] == native['audit_revision'] == 'NativeElementSideV2',
        'inventory, audit, cache or projection controls differ')
    require(saved['signature'] == native['signature'], 'physical checkpoint signature differs')
    require(state_hash(saved['coupled']) == state_hash(native['state']) and
        state_hash(saved['reference']) == state_hash(native['fixed_reference']) and
        state_hash(saved['defect']) == state_hash(native['fixed_defect']), 'saved state/reference/defect differs')
    require([state_hash(v) for v in saved['history']] == [state_hash(v) for v in native['native_history']] and
        state_hash(saved['history_times']) == state_hash(native['native_history_times']), 'native history differs')
    check_report(native, True, definition)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--block-steps', type=int, choices=(8, 32), default=32)
    parser.add_argument('--timeout', type=float, default=7200)
    args = parser.parse_args()
    require(sys.platform == 'linux', 'authoritative continuation is Linux-only')
    args.campaign = args.campaign.resolve(); args.output = args.output.resolve()
    args.targets = [9600]
    root = Path('/media/jay/新加卷/nektar-v2026-runs')
    require(args.campaign.is_relative_to(root) and args.output.is_relative_to(root), 'isolated run storage required')
    args.output.mkdir(parents=True, exist_ok=False)
    owner = CampaignLock(args.output)
    os.environ.update(THREADS)
    receipt = dict(status='validating_completed_parent', started_UTC=datetime.now(timezone.utc).isoformat(),
        operation='CompletedT5ToT6V1', parent_campaign=str(args.campaign),
        targets_steps=args.targets, targets_physical=[6.], observation_steps=[8800, 9600], block_steps=args.block_steps,
        production_dt=.000625, scientific_policy_unchanged=True,
        source_revision=subprocess.check_output(['git', '-C', str(Path(__file__).resolve().parents[4]),
            'rev-parse', 'HEAD'], text=True).strip(), controller_sha256=sha(Path(__file__)),
        verification_governing_sweeps=0, thread_environment=THREADS)
    def save(**changes):
        receipt.update(changes); atomic(args.output/'handoff.json', receipt)
    save()
    try:
        # Acquire positively rather than relying only on historical PID files.
        guard = CampaignLock(args.campaign)
        try:
            old = load(args.campaign/'progress.json')
            require(not (args.campaign/'runtime-handoff.json').exists(), 'parent already has a successor')
            require(not live_commands('supervise_cls_epawf_pimple.py', args.campaign.parent), 'parent supervisor is live')
            require(not matches(process_identity(old['pid']), dict(pid=old['pid'], start=old['process_start'])), 'parent driver is live')
            assert_native_stopped(old)
            args.solver = Path(old['executable'])
            pinned = load(args.campaign/'runtime/identity.json')['frozen']
            require(sha(args.solver) == old['executable_sha256'] == pinned['sha256'], 'frozen parent executable differs')
            for file, digest in pinned['libraries'].items():
                require(sha(Path(file)) == digest, 'frozen linked library differs: '+file)
            cp = Path(old['checkpoint']); saved = checkpoint(cp); h = header(cp, True)
            require(sha(cp/'State.bin') == old['checkpoint_sha256'] == h['sha256'], 'parent checkpoint hash differs')
            segment = old['segments'][-1]
            require(segment['end_step'] == 8000 and segment['exit_code'] == 0 and not segment['failure'], 'parent segment did not finish')
            report = Path(segment['checkpoint']).parent/'case.evolution.json.gz'
            with gzip.open(report, 'rb') as stream: payload = stream.read()
            require(hashlib.sha256(payload).hexdigest() == segment['report_sha256'], 'parent native report bytes differ')
            native = json.loads(payload, parse_int=float)
            validate_extension(old, native, saved, args.targets)
            budget = prefix_budget(args.campaign, h)
            inputs = args.output/'parent'; (inputs/'checkpoint').mkdir(parents=True)
            shutil.copy2(cp/'State.bin', inputs/'checkpoint/State.bin')
            require(sha(inputs/'checkpoint/State.bin') == h['sha256'], 'checkpoint copy differs')
            (inputs/'prefix.evolution.json').write_bytes(payload)
            shutil.copy2(args.campaign/'progress.json', inputs/'progress.json')
            shutil.copy2(args.campaign/'steps.jsonl', inputs/'steps.jsonl')
            shutil.copy2(Path(old['segment'])/'case.xml', inputs/'case.xml')
            save(parent_step=h['steps'], parent_time=h['time'], parent_checkpoint_sha256=h['sha256'],
                parent_state_sha256=state_hash(saved['coupled']), parent_signature=saved['signature'],
                parent_native_history_sha256=[state_hash(v) for v in saved['history']],
                parent_native_history_times_sha256=state_hash(saved['history_times']), prefix_correction_integral=budget,
                parent_runtime=pinned, solver_sha256=old['executable_sha256'],
                old_audit_revision=old['audit_revision'], new_audit_revision=old['audit_revision'],
                evaluation_reuse=old['evaluation_reuse'], status='validating_isolated_restart')
        finally: guard.close()
        fixture = args.output/'verification-input'; (fixture/'start').mkdir(parents=True)
        shutil.copy2(inputs/'checkpoint/State.bin', fixture/'start/State.bin')
        shutil.copy2(inputs/'case.xml', fixture/'production.xml')
        atomic(fixture/'capture.json', dict(segment=dict(start_step=8000, end_step=8001)))
        test = replay(fixture, args.output/'verification-restart', args.solver, 'ExactStepV1', args.timeout)
        with gzip.open(args.output/'verification-restart/case.evolution.json.gz', 'rt') as stream:
            restarted = json.load(stream, parse_int=float)
        initial = restarted['initial_k_bound_check']
        require(initial['activations'] == 0 and initial['raw_sha256'] == initial['corrected_sha256'] and
            test['packets'][0]['stages']['committed'] == state_hash(saved['coupled']), 'startup changed the accepted parent')
        save(status='validated_completed_parent', verification_governing_sweeps=1,
            expected_first_successor_checkpoint_sha256=test['checkpoint_sha256'], startup_bound_inactive=True,
            startup_committed_state_unchanged=True)
        launch_validated_successor(args, old, h, inputs, test['checkpoint_sha256'], save)
    except BaseException as error:
        # A launched successor is never replaced here. Its own supervisor
        # retains numerical-terminal and infrastructure-recovery semantics.
        save(status='extension_observation_failed' if receipt.get('status') == 'successor_launched' else 'extension_failed',
            failure=repr(error)); raise


if __name__ == '__main__': main()
