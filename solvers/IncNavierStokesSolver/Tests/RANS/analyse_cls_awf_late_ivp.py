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
from diagnose_cls_accepted_covariance import checkpoint,common_packet,run_probe
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


def assess(fields):
    """A resolved above-budget component cannot be hidden by another tiny one."""
    rows = list(fields)
    resolved = all(v['order_resolved'] for v in rows)
    adequate = all(v['original_dt_within_target'] if v['order_resolved'] else
                   v['observed_span_below_target'] for v in rows)
    if resolved:
        label = 'coherent_order_and_within_budget' if adequate else 'coherent_order_estimate_above_budget'
    elif adequate:
        label = 'observed_late_sensitivity_below_budget_order_unresolved'
    else:
        label = 'late_sensitivity_requires_further_diagnosis'
    return label, adequate


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
    # Independent full consistent-mass check includes the solved matching
    # trace, whose coefficients can differ between dt levels.
    import numpy as np
    last=str(policy['elapsed_observations'][-1]);coarse=load(a.study/'coarse/progress.json')['observations'][last]
    cp=Path(coarse['checkpoint']);saved=checkpoint(cp)
    metric,rc=run_probe(a.probe,Path(coarse['xml']),dict(common_packet(cp),operation='native_epsilon_metric',
        state=saved['coupled'],history=saved['coupled'],alpha=0.,time=saved['time'],inner_stress_model='Linear'),a.output/'epsilon-mass-crosscheck')
    mass=np.array(metric['matrix']);offset=int(metric['offset']);count=len(metric['coefficients'])
    out['full_epsilon_mass_crosscheck']=dict(receipt=rc,physical_measure=metric['physical_volume'],pairs={})
    require(abs(metric['physical_volume']-21)<1e-10,'incomplete physical outer-epsilon mass')
    for elapsed in policy['elapsed_observations']:
        key = str(elapsed); c, m, f = [data[label][key] for label in ('coarse', 'middle', 'fine')]
        dc = differences(c, m); df = differences(m, f); span=differences(c,f);nc, nf = normalized(dc), normalized(df);ns=normalized(span)
        coefficients=[]
        for label in ('coarse','middle','fine'):
            cp=Path(load(a.study/label/'progress.json')['observations'][key]['checkpoint'])
            coefficients.append(np.array(checkpoint(cp)['coupled'][offset:offset+count]))
        check=[]
        for delta,norm in [(coefficients[0]-coefficients[1],dc['epsilon']),(coefficients[1]-coefficients[2],df['epsilon'])]:
            value=math.sqrt(max(0.,float(delta@mass@delta)));require(abs(value-norm)<1e-12,'full outer epsilon common quadrature disagrees with native mass')
            check.append(dict(native_mass_norm=value,quadrature_norm=norm,absolute_agreement=abs(value-norm)))
        out['full_epsilon_mass_crosscheck']['pairs'][key]=check
        fields = {name: estimate(dc[name], df[name], nc[name], nf[name], policy['middle_step_accuracy_targets'][name], policy)
            for name in ('u', 'v', 'k', 'epsilon', 'pressure')}
        for name,value in fields.items():
            value.update(coarse_finest_difference=span[name],normalized_coarse_finest_difference=ns[name],
                observed_span_below_target=max(nc[name],nf[name],ns[name])<=value['target'])
        tc, tm, tf = [d['mean_wall_stress'] for d in (c, m, f)]
        fields['traction'] = estimate(abs(tc-tm), abs(tm-tf), abs(tc-tm)/.16, abs(tm-tf)/.16,
            policy['middle_step_accuracy_targets']['relative_mean_traction'], policy)
        fields['traction']['values'] = [tc, tm, tf]
        fields['traction'].update(coarse_finest_difference=abs(tc-tf),normalized_coarse_finest_difference=abs(tc-tf)/.16,
            observed_span_below_target=max(abs(tc-tm),abs(tm-tf),abs(tc-tf))/.16<=fields['traction']['target'])
        out['triplets'][key] = fields
    classification, adequate = assess(v for t in out['triplets'].values() for v in t.values())
    out.update(status='passed', classification=classification,
        extend_to_half_unit=False if adequate else None, completed_UTC=datetime.now(timezone.utc).isoformat())
    atomic(a.output/'summary.json', out)


if __name__ == '__main__': main()
