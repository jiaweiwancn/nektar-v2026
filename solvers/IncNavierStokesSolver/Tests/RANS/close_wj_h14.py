#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Verify selected H14 evidence, then freeze its existing ledger once."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
from run_wj_h4 import load,save
from run_wj_h7 import update_journal
from wj6_resources import scientific_workers

def identity(p):
    p=Path(p);h=hashlib.sha256()
    with p.open('rb')as f:
        for block in iter(lambda:f.read(1024**2),b''):h.update(block)
    return dict(path=str(p),bytes=p.stat().st_size,sha256=h.hexdigest())

def verify(root):
    out=root/'closure-v1';out.mkdir()
    verified=[]
    for version in ('v1','v2'):
        pub=root/('publication-'+version)
        manifest=load(pub/'manifest.json')
        for v in manifest['artifacts']:assert identity(v['path'])==v
        for v in load(pub/'heavy-native-manifest.json')['files']:assert identity(v['path'])==v
        logs=load(pub/'compressed-log-bindings.json')
        for v in logs:
            h=hashlib.sha256();count=0
            with gzip.open(v['compressed']['path'],'rb')as f:
                for block in iter(lambda:f.read(1024**2),b''):h.update(block);count+=len(block)
            assert count==v['raw_log']['bytes']and h.hexdigest()==v['raw_log']['sha256']
        verified.append(dict(publication=version,manifest=identity(pub/'manifest.json'),seal=identity(pub/'seal.json'),
            artifacts=len(manifest['artifacts']),decoded_logs=len(logs)))
    states=load(root/'publication-v2/state-bindings.json')['states']
    assert [(s['case'],s['checkpoint_absolute_step'])for s in states]==[('cp0',0),('pilot1',1)]
    for s in states:
        assert identity(s['producer']['path'])==s['producer']and identity(s['checkpoint']['path'])==s['checkpoint']
        a=load(root/f"science/target/{s['case']}/view/result.json")
        assert a['status']=='PASS'and a['absolute_step']==s['checkpoint_absolute_step']
        assert a['report']==s['producer']['sha256']and a['checkpoint']==s['checkpoint']['sha256']
        assert a['initialization_work']==[0,0]and a['initializer_calls']==a['governing_solves']==0
        assert a['scalar_domain']['admissible']
        assessment=root/f"science/target/{s['case']}/"/('assessment-v2'if s['case']=='cp0'else'assessment')/'result.json'
        a=load(assessment);assert a['status']=='PASS'and identity(a['inventory']['path'])==a['inventory']
    pilot=load(root/'science/target/pilot-result.json')
    assert pilot['completed_steps']==pilot['actual_attempts']==1 and pilot['remaining_steps_status']=='RESOURCE_LIMIT'
    assert not any((root/f'science/target/pilot{n}/case.evolution.json').exists()for n in (2,3,4))
    gates=load(root/'science/qualification-gates.json');assert gates['status']=='PASS'
    assert gates['small_governing_verification_attempts']==7
    book=load(root/'accounting.json');assert not book['closed']and not book.get('coverage_halt')and not book.get('runtime_halt')
    assert book['target_CP0_reservations']==1 and book['pilot_reserved_maximum_attempts']==4
    sizes={k:int(subprocess.check_output(['du','-sb',str(root/k)],text=True).split()[0])for k in ('science','build','source')}
    assert sizes['science']<=2*1024**3 and sizes['build']<=8*1024**3
    allocation=scientific_workers();assert allocation['workers']<=72
    save(out/'verification.json',dict(status='PASS',verified=verified,unchanged_target_bindings=states,
        scientific_output_bytes=sizes['science'],build_footprint_bytes=sizes['build'],source_footprint_bytes=sizes['source'],
        current_allocation=allocation,governing_target_attempts=1,accepted_target_steps=1,qualified_target_steps=1,
        unlaunched_target_steps=3,unlaunched_status='RESOURCE_LIMIT',small_verification_attempts=7,
        target_CP0s=1,small_CP0s=1,accepted_checkpoint_native_views=3,
        scope='remote-heavy bytes rehashed in place; compact transport independently verified separately on Mac; no new native/scientific evaluation'))
    shutil.copyfile(__file__,out/'close_wj_h14.py')
    print('closure verification PASS',sizes,flush=True)

def close(root):
    out=root/'closure-v1';assert load(out/'verification.json')['status']=='PASS'
    def finish(b):
        assert not b['closed']
        assert b['entries']['closure-verification-v1']['status']=='PASS'
        assert not any(v['status']=='RESERVED'for v in b['entries'].values())
        phases={p:sum(v.get('seconds',v['reserved_seconds'])for v in b['entries'].values()if v['phase']==p)
            for p in {v['phase']for v in b['entries'].values()}}
        total=b['setup_charged_seconds']+sum(v for p,v in phases.items()if p!='closure')+b['closure_reserve_seconds']
        assert total<=b['ceiling_seconds']and phases.get('closure',0)<=300
        b.update(closed=True,fully_charged_closure_seconds=300,phase_actual_operation_seconds=phases,total_charged_seconds=total,
            final_status='BACKEND_AND_ONE_STEP_PASS; FURTHER_PILOT_RESOURCE_LIMIT; DURATION_AND_HILL_NOT_RUN')
    update_journal(root/'accounting.json',finish)
    shutil.copyfile(root/'accounting.json',out/'accounting-final.json')
    final=load(out/'accounting-final.json')
    # Copy completed publication/closure receipts after they have been finalized.
    for label in ('compact-publication-v1','compact-publication-v2','closure-verification-v1'):
        dest=out/label;dest.mkdir()
        for name in ('request.json','receipt.json','command.log'):shutil.copyfile(root/'operations'/label/name,dest/name)
    save(out/'summary.json',dict(status='CLOSED',charged_seconds=final['total_charged_seconds'],ceiling_seconds=7200,
        setup_charged_seconds=60,closure_charged_seconds=300,closure_actual_operation_seconds=final['phase_actual_operation_seconds']['closure'],
        phase_actual_operation_seconds=final['phase_actual_operation_seconds'],
        peak_successful_native_RSS_bytes=max(final['entries'][k]['peak_tree_rss_bytes']for k in(
            'target-full-domain-CP0-v1','target-cost-gated-pilot-v1','target-CP0-strict-view-v1','target-terminal-strict-view-v1')),
        retained_failed_consumer_RSS_bytes=final['entries']['target-strict-gate-and-pilot-projection-v4']['peak_tree_rss_bytes'],
        unknown_external_precheck_charge_seconds=5,
        publication_v1='superseded compact RSS/figure summary; its original bytes/seal and science remain retained',
        publication_v2='final compact view; accounting snapshot precedes this separately frozen closure',
        closure_rule='300 seconds fully charged once; actual closure operation seconds are reported separately, not charged twice'))
    files=[identity(p)for p in sorted(out.rglob('*'))if p.is_file()]
    save(out/'manifest.json',dict(artifacts=files))
    save(out/'seal.json',dict(manifest=identity(out/'manifest.json'),accounting=identity(out/'accounting-final.json'),
        publication_seal=identity(root/'publication-v2/seal.json'),rule='snapshot saved before hashing; closed journal is never relabelled as the earlier publication snapshot'))
    print('closed charged seconds',final['total_charged_seconds'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('verify','close'),required=True)
    a=p.parse_args();assert a.root.name=='wj-h14-mapped-backend-20261113';globals()[a.mode](a.root)
