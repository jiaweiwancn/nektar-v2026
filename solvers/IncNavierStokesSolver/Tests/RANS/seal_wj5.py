#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Close WJ-5 numerical accounting and verify immutable evidence bridges."""
import argparse
import json
import os
from pathlib import Path
import re
import time
from wj_immutable import Ledger,read,write_once,artifacts,verify_artifacts,file_hash,file_identity,runtime
from wj_qualification import aggregate,component
from run_cls_epawf_picard import require


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--compact',type=Path,required=True);p.add_argument('--figures',type=Path,required=True);a=p.parse_args();r=a.root.resolve();out=a.compact.resolve();start=time.time()
    commands=json.loads((r/'command-receipts.json').read_text());inventory=json.loads((r/'inventory.log').read_text());registered={t['name']:t for t in inventory['tests']};executed={};suites=[]
    for label,count in (('wj5-regression-controls',24),('wj5-capacity-retry-native',3),('wj5-capacity-retry-readonly',2),('wj5-capacity-retry-metadata',1),('wj5-capacity-retry-composite',2)):
        inv=json.loads((r/(label+'-inventory.json')).read_text());account=read(r/(label+'-accounting.json'));require(len(inv['tests'])==count and account['exit_code']==0,'missing/failed required regression selection: '+label)
        receipt=next(c for c in commands if c['label']=='ctest-'+label);require(receipt['exit_code']==0 and file_hash(receipt['log'])==receipt['log_sha256'],'CTest receipt changed')
        names=[t['name'] for t in inv['tests']];executed.update({t['name']:t for t in inv['tests']});suites.append(dict(label=label,count=count,names=names,receipt=receipt,accounting=file_identity(r/(label+'-accounting.json'))))
    require(len(executed)==32,'required WJ-4 regression union is not32 tests')
    focused=next(c for c in commands if c['label']=='wj5-focused-final');require(focused['exit_code']==0 and file_hash(focused['log'])==focused['log_sha256'],'focused native test receipt changed')
    for name in ('WJLocalNormalMetadata','WJLocalNormalObservationEP','WJLocalNormalObservationBudget'):executed[name]=registered[name]
    metadata=next(c for c in commands if c['label']=='measurement-metadata-fixed');require(metadata['exit_code']==0,'corrected metadata selection did not pass')
    closures={}
    for test in executed.values():
        path=Path(test['command'][0]).resolve()
        if path.is_file() and str(path) not in closures:closures[str(path)]=runtime(path)
    write_once(out/'regression-verification.json',dict(schema='WJ5RegressionSealV1',status='PASS',registered=len(registered),selected_unique=len(executed),
        names=sorted(executed),required32=suites,focused=focused,metadata_supplement=metadata,runtime_closures=closures,
        allocation_rejections=[file_identity(f) for f in r.glob('*allocation*.json') if json.loads(f.read_text()).get('allowed') is False]))
    # Verify every original retained-case artifact against the sealed WJ-4
    # index, including fields/options as well as checkpoints and reports.
    source=Path(__file__).resolve().parents[4];old=read(source/'docs/plans/wj4-linux/manifest.json');index_meta=old['artifact_index'];require(file_hash(index_meta['path'])==index_meta['sha256'],'WJ-4 sealed index changed')
    old_index=read(index_meta['path']);oldroot=Path('/media/jay/新加卷/nektar-v2026-runs/wj4-20261110')
    prefixes=[str(oldroot/'study'/p/n)+'/' for p in ('EP','Budget') for n in ('horizontal-Nx1-cp0','horizontal-Nx1-32','horizontal-Nx2-cp0','horizontal-Nx2-32','longer-coarse-32')]
    checked={}
    for name,identity in old_index.items():
        if any(name.startswith(prefix) for prefix in prefixes):
            require(Path(name).is_file() and file_hash(name)==identity['sha256'],'retained WJ-4 bytes changed: '+name);checked[name]=identity
    publication=read(source/'docs/plans/wj4-linux/publication.json');verify_artifacts(publication['final_accounting'])
    require(file_hash(oldroot/'compact/receipt.json')==publication['inner_receipt_sha256'],'WJ-4 inner receipt changed')
    write_once(out/'wj4-integrity.json',dict(schema='WJ5PriorEvidenceIntegrityV1',status='PASS',original_index=index_meta,checked_count=len(checked),checked=checked,
        final_wj4_accounting=publication['final_accounting'],original_wj4_publication=file_identity(source/'docs/plans/wj4-linux/publication.json')))
    # The two original development-only observation files were made before the
    # source path was frozen. Recover their exact code bytes by Git identity;
    # they are not reused for qualification. All governing case-code paths and
    # qualified observer code paths remain frozen and independently hashed.
    development=[]
    import subprocess
    for revision in ('da9518a','7c46393'):
        archive=r/'development-code'/revision;archive.mkdir(parents=True,exist_ok=True)
        for name in ('run_wj5_native.py','wj5_measure.py','test_wj5_metadata.py'):
            relative='solvers/IncNavierStokesSolver/Tests/RANS/'+name
            data=subprocess.check_output(['git','-C',str(source),'show',revision+':'+relative]);path=archive/name
            if path.exists():require(path.read_bytes()==data,'development code archive differs')
            else:
                with path.open('xb') as stream:stream.write(data)
            development.append(dict(revision=revision,relative=relative,archive=file_identity(path),scope='retained development code; not qualified immutable reuse'))
    write_once(out/'development-code-bridge.json',development)
    ledgers={p:Ledger(r,p,cap=800) for p in ('EP','Budget')};accounting={}
    shared=[];specific={p:[] for p in ledgers}
    for c in commands:
        if c['label'].startswith('build'):continue
        policy=next((p for p in ledgers if c['label'].endswith('-'+p)),None)
        if policy:specific[policy].append(c)
        else:shared.append(c)
    # Explicitly cover direct passive inspections and denied allocations which
    # did not enter a driver session. This conservative charge is not a solve.
    for policy,l in ledgers.items():
        d=read(l.file);require(all(e['status']!='RUNNING' for e in d['entries'].values()) and all(s['status']!='RUNNING' for s in d.get('sessions',{}).values()),'unreconciled numerical reservation/session')
        desired=sum(c['seconds'] for c in shared+specific[policy])+120.
        before=l.totals();extra=max(0.,desired-before['seconds'])+max(1.,time.time()-start)
        label='final/outer-process-and-passive-inspection-coverage'
        l.reserve(label,dict(command_receipts=file_identity(r/'command-receipts.json'),policy_commands=[c['label'] for c in specific[policy]],
            shared_commands=[c['label'] for c in shared],unwrapped_passive_and_allocation_reserve_seconds=120.,definition='full outer process time; shared/control work conservatively charged to each policy; already charged sessions are never reduced; native builds recorded separately'),0,max(extra,1.))
        l.finish(label,0,extra,str(out/'publication.json'));totals=l.totals();require(totals['attempts']<=800 and totals['seconds']<=10800,'RESOURCE_LIMIT: final package accounting exceeds caps')
        snapshot=out/('final-ledger-'+policy+'.json');write_once(snapshot,read(l.file));accounting[policy]=dict(totals=totals,ledger=file_identity(snapshot),outer_process_seconds=desired,added_seconds=extra)
    write_once(out/'command-receipts-before-seal.json',commands)
    compact_files=artifacts([f for f in out.iterdir() if f.is_file()]);figure_files=artifacts([f for f in a.figures.iterdir() if f.is_file()])
    policies={}
    for policy in ledgers:
        d=read(out/(policy+'.json'));parts=d['components'];parts['actual_resource']=component('PASS',totals=accounting[policy]['totals']);policies[policy]=aggregate(parts)
    write_once(out/'publication.json',dict(schema='WJ5FinalPublicationV1',inner_receipt=file_identity(out/'receipt.json'),files=compact_files,figures=figure_files,
        policies=policies,accounting=accounting,complete_command_receipts=artifacts([out/'command-receipts-before-seal.json']),governing_attempts=0,
        native_build_seconds=sum(c['seconds'] for c in commands if c['label'].startswith('build')),publication_code=artifacts([Path(__file__)]),
        scope='bounded serial flat2D finite reference-Dirichlet strip; no exact finite-space equilibrium, global stability/realizability, stationarity or long-time qualification',
        closing_receipt='outer command receipt is copied and hash-bound after this seal exits; snapshots are not rewritten'))
    print(out/'publication.json')

if __name__=='__main__':main()
