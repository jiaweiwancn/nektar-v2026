#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Produce additive stage costs and nested equation-solve detail from replays."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path


def summarize(directory):
    r=json.loads((directory/'result.json').read_text());c=json.loads((directory/'capture.json').read_text())
    assert r['bitwise_production_state_and_all_native_histories_equal'] and r['timing']['tree_partition_passed']
    n=r['timing']['physical_steps'];groups=r['timing']['groups']
    def value(name,metric='wall',section='step'):
        return groups.get(section+'/'+name,{}).get(metric+'_per_step',0.)
    stages={};equations={}
    for metric in ('wall','cpu'):
        sweep=value('pimple.sweep',metric);flow=value('flow.total',metric)
        eps=value('epsilon.total',metric);k=value('k.total',metric);bound=value('k.bound',metric)
        step=value('native.step',metric)
        total=r['native_wall_seconds']/n if metric=='wall' else value('process',metric,'block_overhead')
        stages[metric]={
            'Flow stage':flow,'Epsilon stage':eps,'K stage excluding bound':k-bound,'Final K bound':bound,
            'Other sweep checks':sweep-flow-eps-k,'Acceptance audit and other step work':step-sweep,
            'Block startup, reports, output':total-step}
        assert abs(sum(stages[metric].values())-total)<1e-8
        assert min(stages[metric].values())>=-1e-5
        equations[metric]={
            'Flow assembly':value('flow.assembly',metric),
            'Momentum factorization and solve':value('momentum.factorize',metric)+value('momentum.solve_validate',metric),
            'Pressure mass factor and construction':value('pressure.mass_factorize',metric)+value('pressure.build',metric),
            'Pressure factorization and solve':value('pressure.solve_validate',metric),
            'Epsilon assembly':value('epsilon.assembly',metric),
            'Epsilon factorization and solve':value('epsilon.solve_validate',metric),
            'K assembly':value('k.assembly',metric),
            'K factorization and solve':value('k.solve_validate',metric)}
    b=c['recent_completed_blocks'];production_steps=sum(x['segment']['end_step']-x['segment']['start_step'] for x in b)
    return dict(start_step=c['start']['steps'],end_step=c['expected']['steps'],physical_steps=n,
        bitwise_state_and_history_equal=True,additive_seconds_per_step=stages,nested_equation_seconds_per_step=equations,
        production_recent_seconds_per_step=sum(x['segment']['seconds'] for x in b)/production_steps,
        production_recent_steps=production_steps,
        python_report_seconds_per_step={k:r[k]/n for k in ('python_report_load_seconds','python_report_checks_seconds','python_report_gzip_seconds')},
        domain_scans_per_step=groups['step/domain.scalar']['calls']/n,
        domain_scan_seconds_per_step=value('domain.scalar'),
        complete_audit_seconds_per_step=value('native.acceptance_audit'),
        checkpoint_seconds_per_step=value('io.checkpoint'),
        snapshot_seconds_per_step=value('io.attempt_snapshot'),
        block_initialization_seconds=sum(value(name,section='block_overhead')*n for name in ('startup.object','startup.initialize')),
        block_report_seconds=value('output.report',section='block_overhead')*n,
        bound_work=r['bound_work'],capture_manifest=str(directory/'capture.json'),result_manifest=str(directory/'result.json'))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('replays',type=Path)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--plot',type=Path);a=p.parse_args()
    cases={d.name:summarize(d) for d in sorted(a.replays.iterdir()) if (d/'result.json').exists()}
    out=dict(schema=1,captured_UTC=datetime.now(timezone.utc).isoformat(),
        units='seconds per accepted physical step; eight-step immutable native block replays',
        qualification='timing study only; no change to numerical maps or live runtimes',cases=cases)
    a.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    if a.plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(11,5),layout='constrained');bottom=[0.]*len(cases)
        for stage in next(iter(cases.values()))['additive_seconds_per_step']['wall']:
            values=[c['additive_seconds_per_step']['wall'][stage] for c in cases.values()]
            ax.bar(list(cases),values,bottom=bottom,label=stage)
            bottom=[x+y for x,y in zip(bottom,values)]
        ax.set_ylabel('Wall seconds per accepted step');ax.set_title('Linux native replay: additive stage costs')
        ax.tick_params(axis='x',rotation=15);ax.legend(fontsize=8,loc='upper left',bbox_to_anchor=(1,1))
        fig.savefig(a.plot,dpi=170)


if __name__=='__main__':main()
