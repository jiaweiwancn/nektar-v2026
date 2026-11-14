#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only H15 publication verification followed by a separate final seal."""
import argparse
import gzip
import json
import os
from pathlib import Path
import tarfile
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from run_wj_h15 import size
from run_wj_h7 import update_journal
from wj_h15_observer import verify,identity,atomic,vector_hash
from wj3_measure import checkpoint_payload

def verify_all(root):
    out=root/'closure-v1';out.mkdir();pub=root/'publication-v2'
    seal=load(pub/'seal.json');verify(seal['manifest']);m=load(pub/'manifest.json')
    for p in m['artifacts']:verify(p)
    heavy=load(pub/'remote-heavy-manifest.json')
    for p in heavy['files']:verify(p)
    decoded=[]
    for v in load(pub/'compressed-log-bindings.json'):
        verify(v['raw']);verify(v['compressed'])
        import hashlib
        with gzip.open(v['compressed']['path'],'rb')as f:actual=hashlib.sha256(f.read()).hexdigest()
        assert actual==v['raw']['sha256'];decoded.append(v['compressed'])
    bindings=load(pub/'state-bindings.json');v=bindings['completed_but_unqualified']
    cp=checkpoint_payload(verify(v['checkpoint']).parent,include_reference=True);r=load(verify(v['producer']));cb=load(verify(v['callback']))
    for k in ('state','fixed_reference','fixed_defect'):assert vector_hash(cp[k])==vector_hash(r[k])
    assert vector_hash(cb['callback_history'])==v['callback_history_sha256']and cb['alpha_mass']==400 and cb['time']==cp['time']
    assert all(vector_hash(a)==vector_hash(b)for a,b in zip(cp['native_history'],r['native_history']))and cp['native_history_times']==r['native_history_times']
    assert cp['steps']==1 and cp['time']==.0025
    assert not(root/'science/target/monitor/qualified-cursor.json').exists()
    for i in (2,3,4):assert not(root/f'science/target/monitor/step-{i}.ready.json').exists()
    verify(load(root/'science/imports.json')['prior_ledger_unchanged'])
    t=load(root/'science/target/case.evolution.json');assert t['native_completed']==1 and t['independently_qualified']==0 and t['aborted_advances']==0
    assert t['initialization_work']==[0,0]and t['initializer_calls']==0
    assert load(root/'science/final-qualified-cp0-view/result.json')['status']=='PASS'
    assert 'Ran 25 tests'in (root/'operations/post-stop-passive-methods-v5/command.log').read_text()
    science=size(root/'science');entire=size(root)
    source_build=sum(size(p)for p in root.iterdir()if p.name.startswith(('source','build','runtime')))
    # Stronger than the science-only output limit for this small failed prefix.
    assert entire<2*1024**3 and source_build<8*1024**3
    book=load(root/'accounting.json');entries=book['entries'];known=[v.get('peak_tree_rss_bytes',0)for v in entries.values()if v['phase']!='build'];build=[v.get('peak_tree_rss_bytes',0)for v in entries.values()if v['phase']=='build']
    assert max(known)<8*1024**3 and max(build)<16*1024**3
    assert book['target_halt']['native_completed']==1 and book['target_reserved_maximum_attempts']==4 and book['target_native_blocks']==book['target_strict_readbacks']==1
    atomic(out/'verification.json',dict(status='PASS',publication=identity(pub/'manifest.json'),compact_artifacts=len(m['artifacts']),remote_heavy_files=len(heavy['files']),decoded_logs=len(decoded),
        science_bytes=science,entire_new_root_bytes=entire,source_build_runtime_bytes=source_build,peak_scientific_application_tree_RSS=max(known),peak_build_tree_RSS=max(build),
        no_following_target_stage=True,original_parent_ledger_unchanged=True,completed_replay_qualification='OBSTRUCTED_PROTOCOL; never promoted',last_qualified='original H14 CP0',scope='selected new-root bytes rehashed on Linux; not a claim that remote heavy bytes were read on Mac'))
    print('H15 compact/heavy/lineage verification PASS',flush=True)

def finalize(root):
    out=root/'closure-v1';assert load(out/'verification.json')['status']=='PASS'
    def close(b):
        assert not b['closed']
        assert all(v['status']!='RESERVED'for v in b['entries'].values())
        costs={p:sum(v.get('seconds',v['reserved_seconds'])for v in b['entries'].values()if v['phase']==p)for p in {v['phase']for v in b['entries'].values()}}
        actual_closure=costs.pop('closure',0.)
        total=60+sum(costs.values())+300;assert total<=7200 and costs['target']+costs['target_readback']<=1800 and costs['qualification']<=1800 and costs['build']<=1200
        b.update(closed=True,closure_charged_seconds=300,closure_actual_verified_seconds=actual_closure,total_charged_seconds=total,phase_charged_seconds=costs,
            final_outcome='OBSTRUCTED_PROTOCOL',target_actual_attempts=1,target_actual_completed=1,target_independently_qualified=0,target_new_scientific_progress=0,target_unlaunched=3,
            small_verification_attempts=10,fresh_CP0s=0,new_strict_readbacks=2,unknown_duration_charged_seconds=sum(v['seconds']for v in b['entries'].values()if v.get('unknown_duration_charged')))
    update_journal(root/'accounting.json',close);atomic(out/'accounting-final.json',load(root/'accounting.json'))
    atomic(out/'seal.json',dict(accounting=identity(out/'accounting-final.json'),verification=identity(out/'verification.json'),publication=identity(root/'publication-v2/manifest.json'),historical_first_partial=identity(root/'publication-v1-partial-status.json')))
    archive=root/'h15-compact-final-v2.tar.gz'
    with tarfile.open(archive,'x:gz')as t:
        for name in ('publication-v1','publication-v1-partial-status.json','publication-v2','closure-v1'):t.add(root/name,arcname=name)
    atomic(root/'compact-transport-v2.json',dict(archive=identity(archive),publication=identity(root/'publication-v2/manifest.json'),closure=identity(out/'seal.json')))
    print(json.dumps(dict(total_charged_seconds=load(root/'accounting.json')['total_charged_seconds'],transport=identity(archive))),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('verify','finalize'),required=True);a=p.parse_args();globals()[{'verify':'verify_all','finalize':'finalize'}[a.mode]](a.root)
