#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Seal H18 accounting and essential transport without admitting another solve."""
import argparse,json,tarfile,os
from pathlib import Path
from run_wj_h4 import load
from run_wj_h7 import update_journal
from run_wj_h15 import size
from wj_h15_observer import identity,verify,atomic


def close(root):
    book=load(root/'accounting.json');assert not book['closed']
    assert not any(v['status']=='RESERVED'for v in book['entries'].values())
    pub=root/'publication-v1';manifest=load(pub/'manifest.json')
    for v in manifest['files']:verify(v)
    gates=load(pub/'component-gates.json')
    # Partial qualification remains publishable, never relabelled successful.
    if gates['scope_result']=='PASS':
        assert book['small_charged_attempts']==6 and book['target_reserved_attempts']==book['target_charged_attempts']==6
        assert book['target_cases']==['A','B']and book['strict_reservations']==['target']
        assert gates['successful_small_attempts']==5 and gates['deliberate_small_failed_attempts']==1
    old=load(root/'inputs/selected-bootstrap.json')['H16_ledger_unchanged'];verify(old)
    phase={k:sum(v.get('seconds',v['reserved_seconds'])for v in book['entries'].values()if v['phase']==k)for k in ('build','qualification','target','target_readback','operation','closure')}
    assert phase['build']<=1200 and phase['qualification']+phase['target']+phase['target_readback']<=3000 and phase['closure']<=300
    total=60+300+sum(v.get('seconds',v['reserved_seconds'])for v in book['entries'].values()if v['phase']!='closure')
    assert total<=7200
    source_names=('source','build','runtime','prequalification-build-v1','prequalification-build-v2')
    src=sum(size(root/n)for n in source_names if(root/n).exists())
    output=sum(size(p)if p.is_dir()else p.stat().st_size for p in root.iterdir()if p.name not in source_names)
    peaks={kind:max([v.get('peak_tree_rss_bytes',0)for v in book['entries'].values()if(v['phase']=='build')==(kind=='build')]+[0])for kind in ('build','science')}
    assert src<=8*1024**3 and output<=2*1024**3 and peaks['build']<=16*1024**3 and peaks['science']<=8*1024**3
    closure=root/'closure-v1';closure.mkdir()
    expected={v['ELF']['path']:v['ELF']for v in load(root/'science/runtime-build-v1.json')['executables']}
    deps={v['actual']['path']:v['actual']for v in load(root/'science/runtime-dependencies-v1.json')['libraries']}
    samples=[];used={}
    for path in sorted((root/'operations').glob('*/native-exec-samples.json')):
        for record in load(path):
            assert record['status']=='PASS'
            for native in record['native']:verify(expected[native])
            for f in record['files']:
                if '.so'not in Path(f).name and not Path(f).name.startswith('ld-linux'):continue
                assert f in deps,f
                if f not in used:used[f]=identity(verify(deps[f]))
            samples.append(dict(operation=path.parent.name,pid=record['pid'],native=record['native'],libc=record['libc']))
    atomic(closure/'runtime-final-identity.json',dict(status='PASS',actual_native_exec_samples=samples,libraries=list(used.values()),new_ELFs=list(expected.values()),pretarget_qualification=identity(root/'science/runtime-qualified-v1.json')))

    atomic(closure/'verification.json',dict(status='PASS',scientific_scope=gates['scope_result'],phase_seconds=phase,total_charged_seconds=total,ceiling=7200,setup_fully_charged=60,closure_fully_charged=300,
        source_build_runtime_bytes=src,new_output_bytes=output,aggregate_peak_RSS=peaks,H16_ledger_unchanged=old,publication=identity(pub/'manifest.json'),new_scientific_steps=0))
    def seal(b):b.update(closed=True,total_charged_seconds=total,closure_fully_charged_seconds=300,final_verification=identity(closure/'verification.json'))
    update_journal(root/'accounting.json',seal)
    atomic(closure/'accounting-final.json',load(root/'accounting.json'))
    atomic(closure/'manifest.json',dict(files=[identity(p)for p in sorted(closure.iterdir())if p.is_file()],scope='final closed accounting; publication-time snapshot retains its original OPEN identity'))
    atomic(root/'authoritative-index-v1.json',dict(publication=identity(pub/'manifest.json'),closure=identity(closure/'manifest.json'),accounting=identity(closure/'accounting-final.json'),verification=identity(closure/'verification.json')))
    selected=[root/'authoritative-index-v1.json',*closure.glob('*.json')]
    essentials={'component-gates.json','target-timeline.csv','disjoint-inclusive-timing.csv','checkpoint-storage.csv','conditional-cost-storage.json','compression-measurements.json','state-history-bindings.json','all-signed-Q1-moments.csv','compiled-source-summary.json','checkpoint-cost.png','checkpoint-storage.png','heavy-manifest.json','manifest.json','reused-context.json','final-strict-result.json.gz','audited-state-summaries.json.gz','signed-native-stage-data.json.gz'}
    for p in pub.rglob('*'):
        if p.is_file()and(p.name in essentials or p.parent.name=='controls'or p.parent.name=='small'or'operations'in p.parts and p.name in ('receipt.json','attempt-counts.json')):selected.append(p)
    archive=root/'essential-v1.tar.gz'
    with tarfile.open(archive,'w:gz')as t:
        for p in sorted(set(selected)):t.add(p,arcname=str(p.relative_to(root)))
    atomic(root/'essential-transport-intent-v1.json',dict(archive=identity(archive),files=[identity(p)for p in sorted(set(selected))],closure_fully_charged=300,scope='bounded essential subset; original full publication and all heavy arrays remain Linux'))
    print(json.dumps(dict(index=identity(root/'authoritative-index-v1.json'),archive=identity(archive),total_charged_seconds=total,resources=load(closure/'verification.json'))),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();close(a.root)
