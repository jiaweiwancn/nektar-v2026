#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Close the obstructed H20 ledger without changing any earlier receipt."""
import argparse,json,tarfile
from pathlib import Path
from run_wj_h4 import load
from run_wj_h7 import update_journal
from run_wj_h15 import size
from wj_h15_observer import identity,verify,atomic

def close(root):
    book=load(root/'accounting.json');assert not book['closed']and book['target_halt']
    assert not any(e['status']=='RESERVED'for e in book['entries'].values())
    assert book['target_levels']==['M']and book['bootstrap_reserved']==0 and book['strict_reservations']==['M']
    pub=root/'publication-v1'
    for v in load(pub/'manifest.json')['files']:verify(v)
    old=load(root/'inputs/selected.json')['H19_closed_accounting'];verify(old)
    phase={k:sum(e.get('seconds',e['reserved_seconds'])for e in book['entries'].values()if e['phase']==k)for k in ('operation','build','qualification','target','target_readback','closure')}
    total=60+300+sum(e.get('seconds',e['reserved_seconds'])for e in book['entries'].values()if e['phase']!='closure')
    assert total<=7200 and phase['build']<=900 and sum(phase[k]for k in ('qualification','target','target_readback'))<=3600
    source=sum(size(root/n)for n in ('source','source-correction-v1','build','runtime'))
    output=sum(size(p)if p.is_dir()else p.stat().st_size for p in root.iterdir()if p.name not in ('source','source-correction-v1','build','runtime'))
    peak=max(e.get('peak_tree_rss_bytes',0)for e in book['entries'].values())
    assert source<=8*1024**3 and output<=2*1024**3 and peak<=8*1024**3
    closure=root/'closure-v1';closure.mkdir()
    for label in ('publication-v1',):
        for n in ('request.json','receipt.json','command.log'): 
            p=root/'operations'/label/n
            if p.exists():atomic(closure/(label+'-'+n+'.pin.json'),identity(p))
    counts=dict(new_CP0=0,target_native_history_initializations=0,spatial_initializations=0,native_builds=0,
        governor_process_launches=1,target_reserved_attempt_capacity=4,charged_unknown_attempt_capacity=4,observed_governing_attempts=0,observed_completed=0,observed_qualified=0,
        native_governor_full_audit_calls=0,independent_scientific_evaluations=0,final_strict_readbacks=1,final_strict_role='old H18 M step4/T=.01, genuine producer/callback',
        new_scientific_steps=0,maximum_new_physical_time=None,kernel_requests=0,Python_diagnostic_mass_solves=0,
        focused_distinct_methods=26,focused_total_executions=50,registered_tests_not_executed=False)
    atomic(closure/'execution-counts.json',counts)
    atomic(closure/'verification.json',dict(status='CLOSED_WITH_CONFIGURATION_OBSTRUCTION',phase_seconds=phase,total_charged_seconds=total,ceiling=7200,setup_fully_charged=60,closure_fully_charged=300,
        source_build_runtime_bytes=source,new_output_bytes=output,aggregate_peak_RSS_bytes=peak,old_H19_closed_accounting_unchanged=old,publication=identity(pub/'manifest.json'),counts=counts,
        allocation='unused allowance closed, never transferred; four failed/unknown reserved slots remain charged',scientific_conclusion='H20 comparisons NOT_RUN; original timeout rejection preserved'))
    def seal(b):b.update(closed=True,total_charged_seconds=total,closure_fully_charged_seconds=300,final_verification=identity(closure/'verification.json'))
    update_journal(root/'accounting.json',seal)
    atomic(closure/'accounting-final.json',load(root/'accounting.json'))
    atomic(closure/'manifest.json',dict(files=[identity(p)for p in sorted(closure.iterdir())if p.is_file()],scope='separately closed accounting; publication remains its original OPEN snapshot'))
    atomic(root/'authoritative-index-v1.json',dict(publication=identity(pub/'manifest.json'),closure=identity(closure/'manifest.json'),accounting=identity(closure/'accounting-final.json'),verification=identity(closure/'verification.json')))
    selected=[root/'authoritative-index-v1.json',*closure.glob('*.json')]
    selected.extend(p for p in pub.rglob('*')if p.is_file()and p.name!='protected-original-M-strict.json.gz')
    archive=root/'essential-v1.tar.gz'
    with tarfile.open(archive,'w:gz')as t:
        for p in sorted(set(selected)):t.add(p,arcname=str(p.relative_to(root)))
    atomic(root/'essential-transport-intent-v1.json',dict(archive=identity(archive),files=[identity(p)for p in sorted(set(selected))],remote_only=[identity(pub/'protected-original-M-strict.json.gz')]))
    print(json.dumps(dict(index=identity(root/'authoritative-index-v1.json'),archive=identity(archive),verification=load(closure/'verification.json'))),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();close(a.root)
