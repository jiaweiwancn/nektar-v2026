#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Common-state local temporal sensitivity; never accumulated T=0-to-6 error."""
import argparse
from datetime import datetime, timezone
import math
import os
from pathlib import Path

from analyse_cls_accepted_awf import normalized
from analyse_cls_epawf_height_comparison import differences
from continue_cls_epawf_pimple import atomic
from diagnose_cls_accepted_covariance import checkpoint
from qualify_cls_awf_neutral import observe
from run_cls_epawf_picard import THREADS, load, require, sha


def estimate(dc, df, nc, nf, target, policy):
    tiny = min(nc, nf) <= policy['near_zero_normalized_difference']
    p = None if tiny else math.log2(dc/df)
    coherent = p is not None and policy['coherent_order_interval'][0] <= p <= policy['coherent_order_interval'][1]
    middle = 2**p*nf/(2**p-1) if coherent else None
    coarse = 2**p*middle if coherent else None
    return dict(Dc=dc, Df=df, normalized_Dc=nc, normalized_Df=nf, p=p, order_resolved=coherent,
        near_floor=tiny, normalized_fine_estimate=nf/(2**p-1) if coherent else None,
        normalized_middle_estimate=middle, normalized_original_dt_estimate=coarse,
        target=target, original_dt_within_target=coarse is not None and coarse <= target,
        observed_pairwise_sensitivity_below_target=max(nc, nf) <= target,
        interpretation='Richardson estimate only with coherent non-negligible scaling; a small observed difference alone is a finite-window sensitivity result')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('study', 'probe', 'output'): parser.add_argument('--'+name, type=Path, required=True)
    a = parser.parse_args(); os.environ.update(THREADS); a.output.mkdir(parents=True, exist_ok=False)
    source = load(a.study/'summary.json'); require(source['import_validation']['status'] == 'passed', 'unvalidated late-IVP import')
    policy = load(a.study/'declared-policy.json')['late_temporal']
    out = dict(status='running', study_summary_sha256=sha(a.study/'summary.json'), probe_sha256=sha(a.probe), policy=policy,
        started_UTC=datetime.now(timezone.utc).isoformat(), observations={}, triplets={},
        scope='Late local IVP from the same T=6 state; not accumulated initial-to-T=6 error or timestep-independent steady-root proof')
    data = {}
    parent=a.study/'parent';seed=checkpoint(parent)
    seed_selection=dict(checkpoint=str(parent),checkpoint_sha256=sha(parent/'State.bin'),
        xml=str(a.study/'coarse/initial/case.xml'),native_time=seed['time'],native_step=seed['steps'],balance=False)
    out['common_seed'],seed_fields=observe(a.probe,seed_selection,a.output/'common-seed')
    for label in ('coarse', 'middle', 'fine'):
        progress = load(a.study/label/'progress.json'); require(progress['status'] == 'targets_completed', 'late level incomplete/rejected')
        data[label] = {}; out['observations'][label] = {}
        for elapsed, row in progress['observations'].items():
            cp = Path(row['checkpoint']); saved = checkpoint(cp)
            selection = dict(checkpoint=str(cp), checkpoint_sha256=sha(cp/'State.bin'), xml=row['xml'], native_time=saved['time'],
                native_step=saved['steps'], source_report=row['report'], balance=False)
            summary, physical = observe(a.probe, selection, a.output/label/elapsed)
            summary['dimensional_change_from_common_seed']=differences(physical,seed_fields)
            summary['normalized_change_from_common_seed']=normalized(summary['dimensional_change_from_common_seed'])
            out['observations'][label][elapsed] = summary; data[label][elapsed] = physical
            atomic(a.output/'summary.json', out)
    for elapsed in policy['elapsed_observations']:
        key = str(elapsed); c, m, f = [data[label][key] for label in ('coarse', 'middle', 'fine')]
        dc = differences(c, m); df = differences(m, f); nc, nf = normalized(dc), normalized(df)
        fields = {name: estimate(dc[name], df[name], nc[name], nf[name], policy['middle_step_accuracy_targets'][name], policy)
            for name in ('u', 'v', 'k', 'epsilon', 'pressure')}
        tc, tm, tf = [d['mean_wall_stress'] for d in (c, m, f)]
        fields['traction'] = estimate(abs(tc-tm), abs(tm-tf), abs(tc-tm)/.16, abs(tm-tf)/.16,
            policy['middle_step_accuracy_targets']['relative_mean_traction'], policy)
        fields['traction']['values'] = [tc, tm, tf]
        out['triplets'][key] = fields
    resolved = all(v['order_resolved'] for t in out['triplets'].values() for v in t.values())
    small = all(v['observed_pairwise_sensitivity_below_target'] for t in out['triplets'].values() for v in t.values())
    out.update(status='passed', classification='coherent_order_and_within_budget' if resolved and small else
        'observed_late_sensitivity_below_budget_order_unresolved' if small else 'late_sensitivity_requires_further_diagnosis',
        extend_to_half_unit=False if small else None, completed_UTC=datetime.now(timezone.utc).isoformat())
    atomic(a.output/'summary.json', out)


if __name__ == '__main__': main()
