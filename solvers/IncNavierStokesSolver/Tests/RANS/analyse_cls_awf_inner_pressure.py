#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Matched physical readbacks of the declared inner-P2 pressure experiment.

This is a thin consumer of the existing strict checkpoint observer and common
physical-panel reductions. It never launches a time step or modifies a campaign.
Previously verified P1 readbacks are reused only after checking their artifacts.
"""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import time

from analyse_cls_accepted_awf import normalized, reference_errors
from analyse_cls_epawf_height_comparison import differences
from continue_cls_epawf_pimple import CampaignLock, atomic, load_progress, process_identity
from diagnose_cls_accepted_covariance import checkpoint
from qualify_cls_awf_neutral import observe
from run_cls_epawf_picard import THREADS, load, require, sha, state_hash


def verified_physical(row):
    """Bind a compact observation to its immutable native output and seed."""
    receipt = row['receipt']; output = Path(receipt['output'])
    require(sha(output) == receipt['output_sha256'], 'retained readback output differs')
    selection = row['selection']; cp = Path(selection['checkpoint'])
    require(sha(cp/'State.bin') == selection['checkpoint_sha256'], 'selected checkpoint differs')
    saved = checkpoint(cp)
    require(state_hash(saved['coupled']) == row['state_sha256'], 'selected coupled state differs')
    require(state_hash(saved['reference']) == row['reference_sha256'], 'fixed reference differs')
    require([state_hash(v) for v in saved['history']] == row['native_history_sha256'], 'native histories differ')
    d = load(output); q = d['diagnostics']['physical']
    require(q['epsilon_reporting_lower'] == .5 and q['common_reporting_queries_valid'], 'incomplete physical comparison')
    require(d['governing_solves'] == d['physical_commits'] == 0 and not d['native_history_advanced'], 'readback advanced state')
    return q


def select(progress, step):
    matches = [s for s in progress['segments'] if int(s['end_step']) == step and
               not s['failure'] and s['exit_code'] == 0]
    if not matches: return None
    require(len(matches) == 1, 'ambiguous accepted native observation')
    cp = Path(matches[0]['checkpoint']); saved = checkpoint(cp)
    require(saved['steps'] == step, 'selected native global count differs')
    folder = cp.parent; report = folder/'case.evolution.json.gz'
    require(report.exists(), 'native block has not published its compressed report')
    return dict(checkpoint=str(cp), checkpoint_sha256=sha(cp/'State.bin'), xml=str(folder/'case.xml'),
                native_step=step, native_time=saved['time'], source_report=str(report), H=.5)


def comparison(left, right):
    a, b = verified_physical(left), verified_physical(right)
    require(abs(a['time']-b['time']) < 1e-13, 'physical observation times differ')
    delta = differences(a, b)
    return dict(physical_time=a['time'], dimensional_field_difference=delta,
        common_normalized_field_difference=normalized(delta),
        P2_minus_P1_mean_traction=a['mean_wall_stress']-b['mean_wall_stress'],
        P2_reference_error=reference_errors(a), P1_reference_error=reference_errors(b),
        P2_common_normalized_reference_error=normalized(reference_errors(a)),
        P1_common_normalized_reference_error=normalized(reference_errors(b)),
        P2_strong_continuity=left['complete_audit']['strong_continuity'],
        P1_strong_continuity=right['complete_audit']['strong_continuity'],
        P2_interface=left['complete_audit']['interface'], P1_interface=right['complete_audit']['interface'],
        semantics='P2 minus P1 on common physical panels; one pressure-mean shift across both regions, retaining their relative offset; full .5<=z<=4 transported epsilon. No coefficient-array comparison across pressure spaces.')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('campaign', 'probe', 'initial', 'output'): p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--baseline', type=Path, nargs='+', required=True, help='verified qualify_cls_awf_neutral summary files')
    p.add_argument('--steps', type=int, nargs='+', default=[0, 32, 80, 160, 800])
    p.add_argument('--watch', action='store_true'); p.add_argument('--interval', type=float, default=30.)
    a = p.parse_args(); os.environ.update(THREADS); a.output.mkdir(parents=True, exist_ok=False)
    lock = CampaignLock(a.output)
    require(a.steps == sorted(set(a.steps)) and a.steps[0] == 0, 'ordered observations including native h0 required')
    baseline = {}; receipts = []
    for file in a.baseline:
        doc = load(file); require(doc['probe_sha256'] == sha(a.probe), 'baseline native probe differs')
        receipts.append(dict(path=str(file), sha256=sha(file)))
        for label, row in doc['cases'].items():
            if not label.startswith('refined-'): continue
            n = int(row['selection']['native_step'])
            if n in a.steps:
                require(n not in baseline, 'duplicate baseline observation')
                verified_physical(row); baseline[n] = row
    require(sorted(baseline) == a.steps, 'baseline lacks a requested exact observation')
    initial = load(a.initial); cp = Path(initial['checkpoint']); saved = checkpoint(cp)
    require(saved['steps'] == 0 and saved['time'] == 0., 'fresh native initialization required')
    selection0 = dict(checkpoint=str(cp), checkpoint_sha256=sha(cp/'State.bin'),
        xml=str(a.initial.parent/'case.xml'), source_report=str(a.initial), native_step=0, native_time=0., H=.5)
    out = dict(status='observing', started_UTC=datetime.now(timezone.utc).isoformat(),
        process_identity=process_identity(os.getpid()), script_sha256=sha(Path(__file__)), probe_sha256=sha(a.probe),
        campaign=str(a.campaign), baseline_summaries=receipts, requested_steps=a.steps, observations={}, comparisons={},
        governing_solves=0, physical_commits=0, native_history_advanced=False)
    try:
        while True:
            progress = load_progress(a.campaign/'progress.json')
            require(progress['accepted_awf_pressure_space'] == 'inner-p2-v1' and
                    progress['accepted_awf_inner_refinement'] == 'cutoff-v1', 'not the declared pressure experiment')
            out.update(campaign_status=progress['status'], accepted_steps=progress['accepted_steps'], accepted_time=progress['accepted_time'],
                checkpoint=progress.get('checkpoint'), checkpoint_sha256=progress.get('checkpoint_sha256'),
                updated_UTC=datetime.now(timezone.utc).isoformat())
            for n in a.steps:
                if str(n) in out['observations']: continue
                selection = selection0 if n == 0 else select(progress, n)
                if selection is None: continue
                if n == 0 and progress['accepted_steps']:
                    with (a.campaign/'steps.jsonl').open() as f: first = json.loads(next(f))
                    require(first['global_step'] == 1 and first['attempt']['incoming_sha256'] == state_hash(saved['coupled']),
                        'verified fresh h0 does not equal the production incoming seed')
                    out['production_h0_bitwise_verified'] = True
                row, physical = observe(a.probe, selection, a.output/f'p2-{n:06d}')
                row['dimensional_reference_errors'] = reference_errors(physical)
                row['common_normalized_reference_errors'] = normalized(row['dimensional_reference_errors'])
                row['drift_from_actual_h0'] = differences(physical, verified_physical(out['observations']['0'])) if n else {}
                out['observations'][str(n)] = row
                out['comparisons'][str(n)] = comparison(row, baseline[n])
                atomic(a.output/'summary.json', out)
            if len(out['observations']) == len(a.steps): out['status'] = 'completed_requested_readbacks'; break
            if progress['status'] in ('numerically_rejected', 'targets_completed', 'stopped_at_checkpoint'):
                out.update(status='campaign_terminal_before_all_observations', failure=progress.get('failure')); break
            atomic(a.output/'summary.json', out)
            if not a.watch: break
            time.sleep(a.interval)
    except BaseException as error:
        out.update(status='readback_failed_production_unchanged', failure=str(error)); atomic(a.output/'summary.json', out); raise
    out.update(completed_UTC=datetime.now(timezone.utc).isoformat()); atomic(a.output/'summary.json', out)


if __name__ == '__main__': main()
