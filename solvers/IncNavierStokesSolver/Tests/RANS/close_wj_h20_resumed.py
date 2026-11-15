#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Close linked H20 accounting, retaining all historical and resumed failures."""
import argparse,json,tarfile
from pathlib import Path
from run_wj_h4 import load
from run_wj_h7 import update_journal
from run_wj_h20 import usage
from run_wj_h15 import size
from wj_h15_observer import identity,verify,atomic
from wj_h17_observer import native_load

def close(root):
    b=load(root/'accounting.json');assert not b['closed']and not any(v['status']=='RESERVED'for v in b['entries'].values())
    verify(b['carryforward_accounting']);old=load(b['carryforward_accounting']['path']);assert old['closed']and old['total_charged_seconds']==b['carryforward_seconds']
    assert b['target_reserved_attempts']<=22 and b['bootstrap_reserved']<=1 and len(b['strict_reservations'])<=3
    used=usage(b);assert used['total']<=7200 and used['native']<=3600
    pub=root/'publication-v1'
    for pin in load(pub/'manifest.json')['files']:verify(pin)
    closure=root/'closure-v1';closure.mkdir();counts=[]
    for level in b.get('target_levels',[]):
        p=root/f'science/{level}/case.evolution.json'
        if p.exists():
            r=native_load(p);counts.append(dict(level=level,producer=identity(p),attempted=r['implicit_callbacks'],completed=r['native_completed'],independent_checkpoint_commits=r['independently_qualified'],absolute_step=r['steps'],time=r['time'],failure=r['failure'],full_native_audits=r['full_audit_calls'],PhysicalSamples=r['physical_samples_calls']))
    output=sum(size(p)if p.is_dir()else p.stat().st_size for p in root.iterdir()if not p.name.startswith(('source','build','runtime')))
    sources=sum(size(p)if p.is_dir()else p.stat().st_size for p in root.iterdir()if p.name.startswith(('source','build','runtime')))
    prior_verification=load(verify(old['final_verification']))
    peak=max([prior_verification['aggregate_peak_RSS_bytes']]+[v.get('peak_tree_rss_bytes',0)for v in b['entries'].values()])
    assert output+prior_verification['new_output_bytes']<=2*1024**3 and sources+prior_verification['source_build_runtime_bytes']<=8*1024**3 and peak<=8*1024**3
    atomic(closure/'execution-counts.json',dict(new_branches=counts,new_governing_attempts=sum(c['attempted']for c in counts),new_completed=sum(c['completed']for c in counts),new_CP0=b['bootstrap_reserved'],new_final_strict_reads=len(b['strict_reservations']),
        historical_configuration_launches=1,historical_unknown_charged_attempt_capacity=4,historical_strict_reads=1,aggregate_reserved_attempt_capacity=4+b['target_reserved_attempts'],
        native_builds=0,spatial_initializations=0,target_native_history_initializations=b['bootstrap_reserved'],maximum_physical_horizon=.02,new_focused_distinct_methods=29,new_focused_method_executions=56,old_focused_method_executions=50))
    atomic(closure/'verification.json',dict(status='CLOSED',aggregate_total_charged_seconds=used['total'],aggregate_native_category_seconds=used['native'],ceiling=7200,native_ceiling=3600,
        carryforward_seconds=b['carryforward_seconds'],new_setup_fully_charged=60,new_closure_fully_charged=300,old_closed_accounting=b['carryforward_accounting'],
        new_output_bytes_at_snapshot=output,old_output_bytes=prior_verification['new_output_bytes'],new_source_build_runtime_bytes=sources,old_source_build_runtime_bytes=prior_verification['source_build_runtime_bytes'],aggregate_peak_science_RSS_bytes=peak,
        counts=identity(closure/'execution-counts.json'),publication=identity(pub/'manifest.json'),gates=load(pub/'component-gates.json'),
        retained_failure='old configuration rejection and new M passive postcondition receipt remain immutable; M native success independently verified without rerun',
        limits='finite checkpoint inventories and short initial relaxation; no production timestep, duration/archive or physical-experiment admission'))
    for name in ('resumption.json','active-execution-config.json','execution-config-v2.json','passive-resolution-v1.json'):
        p=root/name
        if p.exists():atomic(closure/(name+'.binding.json'),identity(p))
    def seal(book):book.update(closed=True,aggregate_total_charged_seconds=used['total'],new_charged_seconds=used['total']-book['carryforward_seconds'],new_closure_fully_charged=300,final_verification=identity(closure/'verification.json'))
    update_journal(root/'accounting.json',seal)
    atomic(closure/'accounting-final.json',load(root/'accounting.json'))
    atomic(closure/'manifest.json',dict(files=[identity(p)for p in sorted(closure.iterdir())if p.is_file()],scope='separately CLOSED linked accounting; original prefix remains closed unchanged'))
    atomic(root/'authoritative-index-v1.json',dict(publication=identity(pub/'manifest.json'),closure=identity(closure/'manifest.json'),accounting=identity(closure/'accounting-final.json'),verification=identity(closure/'verification.json'),historical_accounting=b['carryforward_accounting']))
    selected=[root/'authoritative-index-v1.json',*closure.glob('*.json'),*pub.rglob('*.json'),*pub.rglob('*.csv'),*pub.rglob('*.gz'),*pub.rglob('*.png')]
    selected.extend(root/n for n in ('resumption.json','active-execution-config.json','execution-config-v2.json','passive-resolution-v1.json'))
    selected.extend(p for p in(root/'operations').glob('*/receipt.json'))
    archive=root/'essential-v1.tar.gz'
    with tarfile.open(archive,'w:gz')as t:
        for p in sorted(set(selected)):
            if p.is_file():t.add(p,arcname=str(p.relative_to(root)))
    atomic(root/'essential-transport-intent-v1.json',dict(archive=identity(archive),files=[identity(p)for p in sorted(set(selected))if p.is_file()],scope='compact scientific results/source/runtime/closure; heavy arrays and checkpoints remain Linux'))
    print(json.dumps(dict(index=identity(root/'authoritative-index-v1.json'),archive=identity(archive),verification=load(closure/'verification.json'))),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();close(a.root)
