#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Verify and separately seal H17; never admit additional numerical work."""
import argparse
import gzip
import json
import os
from pathlib import Path
import shutil
import tarfile
import time
from run_wj_h4 import load
from run_wj_h7 import update_journal
from run_wj_h15 import size
from wj_h15_observer import identity,verify,atomic


def finalize(root):
    start=time.monotonic();b=load(root/'accounting.json');assert not b['closed']
    assert not any(v['status']=='RESERVED'for v in b['entries'].values())
    assert b['full_audit_reserved']==b['full_audit_actual_calls']==12
    assert b['small_charged_attempts']==3 and b['target_charged_attempts']==b['target_reserved_attempts']==4
    assert set(b['strict_reservations'])=={'small','target'}
    assert not any(b.get(k)for k in ('coverage_halt','runtime_halt','qualification_halt','target_halt'))
    pub=root/'publication-v1';m=load(pub/'manifest.json')
    for v in m['files']:verify(v)
    gates=load(pub/'component-gates.json');assert gates['scope_result']=='PASS'and gates['new_scientific_steps']==0
    original=load(root/'inputs/selected-bootstrap.json')['H16_ledger_unchanged'];verify(original)
    for scope in ('small/A','target'):
        r=load(root/('science/'+scope+'/view/result.json'));assert r['status']=='PASS'
    target=root/'science/target';outcome=load(target/'target-outcome.json')
    assert outcome['qualified']==outcome['native_completed']==outcome['actual_attempts']==4 and outcome['new_scientific_progress']==0
    cursor=load(target/'monitor/qualified-cursor.json');assert cursor['binding']['absolute_step']=='4'
    assert cursor['binding']['state_sha256']=='4c36141ce8638d049f2526112b51d91faafe0c1b2920a875e14ed2433d901661'
    for n in range(1,5):
        c=load(target/f'monitor/step-{n}.qualification-commit.json')
        for k in ('checkpoint','callback','report','consumed_ack','prepared_cursor'):verify(c[k])
        assert load(target/f'monitor/step-{n}.native-replay-parity.json')['status']=='PASS'
        assert load(target/f'monitor/step-{n}.independent-replay-parity.json')['status']=='PASS'
    supplement=root/'publication-supplement-v1';supplement.mkdir()
    copies=[]
    def copy(p,name):
        p=Path(p)
        if not p.exists():return
        dest=supplement/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest);assert identity(p)['sha256']==identity(dest)['sha256'];copies.append(dict(original=identity(p),copy=identity(dest)))
    for name in ('after-baseline-admission.json','passive-gate-syntax-failure-v1.json','integration-entry-compatibility-v3.json','closure-code-bindings.json'):
        copy(root/'science'/name,'controls/'+name)
    for directory in ('pre-batch-config-v1','pre-monitor-config-v2'):
        for p in (root/directory).glob('*'):
            if p.is_file():copy(p,'retained-configs/'+directory+'/'+p.name)
    retained=[]
    for p in (root/'pre-batch-config-v1/runtime').glob('*.elf'):retained.append(identity(p))
    atomic(supplement/'retained-ELFs.json',dict(files=retained,scope='pre-batch native fixture binaries; distinct from the final active governor/reader'))
    signed=[]
    for n in range(1,5):
        ready=load(target/f'monitor/step-{n}.ready.json');r=load(verify(ready['report']))
        signed.append(dict(step=n,time=r['time'],native_hard_audit=r['native_hard_audit'],stage_audits=r['stage_audits'],step_audits=r['step_audits'],explicit_audits=r['explicit_audits'],actual_explicit_calls=r['actual_explicit_calls'],pimple=r['pimple_attempts'][0],immutable_rows=r['immutable_row_reuse'],report=ready['report'],callback=ready['callback']))
    atomic(supplement/'signed-stage-and-row-history.json',signed)
    copy(target/'view/result.json','terminal-strict-result.json')
    # Selected, previously sealed H16 fetch/exterior evaluations are reusable
    # because every H17 state and the independent numerical functions match.
    # Read the manifest, not the historical heavy corpus; hash only these files.
    old=root.parent/'wj-h16-repaired-monitor-20261114'
    old_index=load(old/'authoritative-index-v1.json');old_manifest=load(verify(old_index['publication']))
    record=next(v for v in old_manifest['files']if v['path']==str(old/'publication-v1/remote-heavy-manifest.json'))
    heavy=load(verify(record));pins={v['path']:v for v in heavy['files']};reused=[]
    for n in range(5):
        for role in ('fetch','exterior'):
            p=old/f'science/analysis-v1/step-{n}-{role}.json';pin=pins[str(p)];verify(pin)
            dest=supplement/f'reused-H16/step-{n}-{role}.json.gz';dest.parent.mkdir(exist_ok=True)
            with p.open('rb')as f,gzip.open(dest,'wb')as g:shutil.copyfileobj(f,g)
            reused.append(dict(original=pin,stored=identity(dest)))
    atomic(supplement/'reused-fetch-exterior-bindings.json',dict(status='REUSED_EXACT_STATE_AND_MATH',original_manifest=record,files=reused,new_scans=0,scope='copied H16 final-union fetch/exterior data; current H17 live results and bitwise state/callback parity remain separately recorded'))
    for name in ('receipt.json','command.log'):
        copy(root/'operations/compact-publication-v1'/name,'publication-operation/'+name)
    atomic(supplement/'qualification-facts.json',dict(status='PASS',focused_python_distinct_methods=13,native_catalog_units=1,read_only_small_native_fixtures=1,passive_target_processes=2,full_target_passive_complete_audits=12,small_governing_attempts=3,target_replay_attempts=4,target_completed=4,target_committed_qualified=4,new_scientific_timesteps=0,fresh_CP0=0,strict_open_checks=2,governing_initialization_calls=0,scalar_JSON_inflation=False,independent_rows_rebuilt=0,new_quadrature_campaign='NOT_RUN; unchanged native projected weak loads retained; no new universal q12 qualification',attribution='current-state contractions separated from combined accessors; lookup, validation, allocation and copy subcosts are not individually isolated'))
    atomic(supplement/'copied-evidence.json',copies)
    atomic(supplement/'manifest.json',dict(files=[identity(p)for p in sorted(supplement.rglob('*'))if p.is_file()]))
    source_names=('source','build','runtime','pre-batch-config-v1','pre-monitor-config-v2','closure-code')
    source_bytes=sum(size(root/n)for n in source_names if(root/n).exists())
    science_bytes=sum(size(p)if p.is_dir()else p.stat().st_size for p in root.iterdir()if p.name not in source_names)
    assert source_bytes<=8*1024**3 and science_bytes<=2*1024**3
    closed=root/'closure-v1';closed.mkdir()
    phases={p:sum(v.get('seconds',v['reserved_seconds'])for v in b['entries'].values()if v['phase']==p)for p in ('build','qualification','target','target_readback','closure','operation')}
    assert phases['build']<=1200 and phases['qualification']<=1800 and phases['target']+phases['target_readback']<=1800 and phases['closure']+time.monotonic()-start<=300
    total=60+300+sum(v.get('seconds',v['reserved_seconds'])for v in b['entries'].values()if v['phase']!='closure')
    assert total<=7200
    peaks={p:max([v.get('peak_tree_rss_bytes',0)for v in b['entries'].values()if(v['phase']=='build')==(p=='build')]+[0])for p in ('build','science')}
    assert peaks['build']<=16*1024**3 and peaks['science']<=8*1024**3
    atomic(closed/'verification.json',dict(status='PASS',publication=identity(pub/'manifest.json'),supplement=identity(supplement/'manifest.json'),H16_ledger_unchanged=original,phase_seconds=phases,source_build_runtime_bytes=source_bytes,new_output_bytes_before_archives=science_bytes,aggregate_RSS_sampled_peaks=peaks,counts=load(supplement/'qualification-facts.json'),closure_fully_charged=300,ceiling=7200))
    def close(book):
        assert not book['closed'];book.update(closed=True,total_charged_seconds=total,closure_charged_seconds=300,closure_scope='fully charged publication/verification/transport allowance; immutable logs only after operation completion',phase_seconds=phases)
    update_journal(root/'accounting.json',close)
    atomic(closed/'accounting-final.json',load(root/'accounting.json'))
    atomic(closed/'manifest.json',dict(files=[identity(p)for p in sorted(closed.iterdir())if p.is_file()]))
    index=dict(schema='WJH17AuthoritativeIndexV1',publication=identity(pub/'manifest.json'),supplement=identity(supplement/'manifest.json'),closure=identity(closed/'manifest.json'),accounting=identity(closed/'accounting-final.json'),verification=identity(closed/'verification.json'),new_scientific_steps=0)
    atomic(root/'authoritative-index-v2.json',index)
    archive=root/'h17-compact-v1.tar.gz'
    with tarfile.open(archive,'w:gz')as t:
        for p in (pub,supplement,closed,root/'authoritative-index-v2.json'):t.add(p,arcname=p.name)
    essential=root/'h17-essential-v1.tar.gz'
    with tarfile.open(essential,'w:gz')as t:
        for p in (root/'authoritative-index-v2.json',closed,pub/'component-gates.json',pub/'conditional-cost.json',pub/'saved-audit-cost.csv',pub/'target-timeline.csv',pub/'target-native-timing.csv',pub/'row-reuse-cost.png',supplement/'qualification-facts.json'):
            if p.exists():t.add(p,arcname=str(p.relative_to(root)))
    atomic(root/'transport-index-v1.json',dict(essential=identity(essential),complete=identity(archive),scope='closed publication/closure archive; active logs excluded; transport verification remains separate'))
    print(json.dumps(dict(status='PASS',total_charged_seconds=total,source_bytes=source_bytes,science_bytes=science_bytes,transport=load(root/'transport-index-v1.json'))),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();finalize(a.root)
