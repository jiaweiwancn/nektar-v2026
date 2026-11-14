#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Verify selected new H16 bytes, then close its ledger exactly once."""
import argparse
import json
import os
from pathlib import Path
import tarfile
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import verify,identity,atomic,vector_hash
from study_wj_h15 import qualified_step
from wj3_measure import checkpoint_payload
from run_wj_h15 import size
from run_wj_h7 import update_journal

def verify_all(root):
    out=root/'closure-v1';out.mkdir();pub=root/'publication-v1';verify(load(pub/'seal.json')['manifest']);manifest=load(pub/'manifest.json')
    for v in manifest['files']:verify(v)
    heavy=load(pub/'remote-heavy-manifest.json')
    for v in heavy['files']:verify(v)
    verify(load(root/'inputs/essential-index.json')['old_ledger_unchanged'])
    cfg=load(root/'execution-config.json')
    for v in cfg['files']:verify(v)
    folder=root/'science/target';q=qualified_step(folder);bound=[]
    for step in range(1,q+1):
        c=load(folder/f'monitor/step-{step}.qualification-commit.json')
        for k in ('checkpoint','callback','report','consumed_ack','prepared_cursor'):verify(c[k])
        r=load(c['report']['path']);cp=checkpoint_payload(Path(c['checkpoint']['path']).parent,include_reference=True);cb=load(c['callback']['path'])
        for k in ('state','fixed_reference','fixed_defect'):assert vector_hash(cp[k])==vector_hash(r[k])
        assert all(vector_hash(a)==vector_hash(b)for a,b in zip(cp['native_history'],r['native_history']))and cp['native_history_times']==r['native_history_times']
        assert cp['steps']==step and cb['step']==step and cb['alpha_mass']==400 and cb['time']==cp['time']
        assert vector_hash(cb['callback_history'])==c['binding']['callback_history_sha256']and c['binding']['state_sha256']==vector_hash(cp['state'])
        bound.append(dict(step=step,checkpoint=c['checkpoint'],producer=c['report'],callback=c['callback'],commit=identity(folder/f'monitor/step-{step}.qualification-commit.json')))
    b=load(root/'accounting.json');assert b['small_charged_attempts']==6 and b['target_blocks']==1 and b['target_charged_attempts']<=4 and sorted(b['strict_reservations'])==['small','target']
    build=sum(size(p)if p.is_dir()else p.stat().st_size for p in root.iterdir()if p.name.startswith(('source','build','runtime')));entire=size(root);science=size(root/'science');other=entire-build
    assert build<8*1024**3 and other<2*1024**3
    peak_build=max(v.get('peak_tree_rss_bytes',0)for v in b['entries'].values()if v['phase']=='build');peak_science=max(v.get('peak_tree_rss_bytes',0)for v in b['entries'].values()if v['phase']!='build')
    assert peak_build<16*1024**3 and peak_science<8*1024**3
    atomic(out/'verification.json',dict(status='PASS',scope_result=load(pub/'component-gates.json')['scope_result'],compact_files=len(manifest['files']),remote_heavy_files=len(heavy['files']),committed_qualified_states=bound,source_build_runtime_bytes=build,science_bytes=science,conservative_new_scientific_output_bytes=other,entire_root_bytes=entire,peak_build_tree_RSS=peak_build,peak_scientific_application_tree_RSS=peak_science,parent_H15_ledger_unchanged=True,original_native_and_independent_decisions_preserved=True))
    print('H16 compact/heavy/lineage/limits verification PASS',flush=True)

def finalize(root):
    out=root/'closure-v1';assert load(out/'verification.json')['status']=='PASS';gates=load(root/'publication-v1/component-gates.json')
    def close(b):
        assert not b['closed']and all(v['status']!='RESERVED'for v in b['entries'].values())
        phases={p:sum(v.get('seconds',v['reserved_seconds'])for v in b['entries'].values()if v['phase']==p)for p in{v['phase']for v in b['entries'].values()}}
        actual_closure=phases.pop('closure',0.);total=60+sum(phases.values())+300
        assert total<=7200 and phases.get('build',0)<=1200 and phases.get('qualification',0)<=1800 and phases.get('target',0)+phases.get('target_readback',0)<=1800 and actual_closure<=300
        b.update(closed=True,closure_charged_seconds=300,closure_actual_verified_seconds=actual_closure,total_charged_seconds=total,phase_charged_seconds=phases,final_outcome=gates['scope_result'],target_actual_completed=gates['target'].get('native_completed',0),target_committed_qualified=gates['target'].get('independently_committed',0),new_scientific_progress=gates['target'].get('new_scientific_progress',0),fresh_CP0s=0,new_strict_readbacks=len(b['strict_reservations']),unknown_duration_charged_seconds=sum(v.get('seconds',v['reserved_seconds'])for v in b['entries'].values()if v.get('unknown_attempts')))
    update_journal(root/'accounting.json',close);atomic(out/'accounting-final.json',load(root/'accounting.json'))
    atomic(out/'source-runtime-index.json',dict(native=identity(root/'science/runtime-build-v1.json'),compiled_inputs=identity(root/'science/compiled-source-objects.json'),frozen_execution=identity(root/'execution-config.json'),actual_native_maps=identity(root/'publication-v1/native-runtime-observed.json'),historical_cache_compatibility=identity(root/'science/target-cache-compatibility.json'),essential_H15_provenance=identity(root/'inputs/essential-index.json')))
    atomic(out/'seal.json',dict(accounting=identity(out/'accounting-final.json'),verification=identity(out/'verification.json'),source_runtime=identity(out/'source-runtime-index.json'),publication=identity(root/'publication-v1/manifest.json'),publication_snapshot_is_not_final_accounting=True))
    archive=root/'h16-compact-v1.tar.gz'
    with tarfile.open(archive,'x:gz')as t:
        for n in ('publication-v1','closure-v1'):t.add(root/n,arcname=n)
    atomic(root/'authoritative-index-v1.json',dict(status=gates['scope_result'],root=str(root),publication=identity(root/'publication-v1/manifest.json'),component_gates=identity(root/'publication-v1/component-gates.json'),closure=identity(out/'seal.json'),accounting=identity(out/'accounting-final.json'),archive=identity(archive),heavy_transport='NOT_REQUESTED: selected files verified in place; original coefficient/matrix bytes remain Linux'))
    print(json.dumps(dict(index=identity(root/'authoritative-index-v1.json'),closed_charge=load(root/'accounting.json')['total_charged_seconds'],archive=identity(archive))),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('verify','finalize'),required=True);a=p.parse_args();(verify_all if a.mode=='verify'else finalize)(a.root)
