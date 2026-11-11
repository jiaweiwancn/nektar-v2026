#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Close new WJ-7 ledgers and export exact/versioned compact evidence."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import time
from wj_immutable import Ledger,read,write_once,artifacts,file_identity,file_hash,verify_artifacts
from wj6_resources import DEFAULT,volume
from wj_qualification import aggregate,component
from run_cls_epawf_picard import require


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('root','parent','compact','figures','output'):p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();root=a.root.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.time()
    ledgers={p:Ledger(root,p,cap=800) for p in ('EP','Budget')}
    for l in ledgers.values():l.reserve('final-seal-transport-reserve',dict(code=file_identity(__file__)),0,300)
    verify_artifacts(read(a.compact/'receipt.json')['files'])
    require(file_hash(a.parent/'export-v1/transport.json')=='74649d94ba88f1ddb20e5f2ca4642a7a0bfb36da47290b195cf862cda79017df','parent seal changed')
    parent_index=read(a.parent/'export-v1/transport.json')['files']
    for name in ('publication.json.gz','final-ledger-EP.json','final-ledger-Budget.json','accepted-stage-lineage.json.gz'):
        r=parent_index[name]['original'];verify_artifacts({r['path']:r})
    heavy=read(a.compact/'heavy-index.json');verify_artifacts(heavy)
    resource={};counts={};phases={}
    for policy in ledgers:
        study=root/'study-v2'/policy;s=read(study/'summary.json');packets=[]
        for e in s['execution'].values():
            if not e.get('report'):continue
            verify_artifacts({e['report']['path']:e['report']});d=read(e['report']['path']);packets.extend(d['pimple_attempts'])
        keys=('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls',
            'bound_activations','additional_governing_equation_solves','additional_correction_or_restoration_solves')
        counts[policy]=dict(attempts=len(packets),**{k:sum(r.get(k,0) for r in packets) for k in keys})
        require(counts[policy]['bound_activations']==0,'correction is not inactive')
        for k in keys[:5]:require(counts[policy][k]==len(packets),'solve/bound count mismatch')
        require(counts[policy][keys[6]]==counts[policy][keys[7]]==0,'additional governing solve')
        peaks=[c['whole_process_peak_rss_bytes'] for c in read(root/'command-receipts.json')]
        resource[policy]=dict(whole_process_peak_rss_bytes=max(peaks),conservatively_counted_all_output_roots_bytes=volume(root),limits=DEFAULT)
        require(max(peaks)<=DEFAULT['rss_bytes'] and volume(root)<=DEFAULT['policy_output_bytes'],'RESOURCE_LIMIT: actual memory/storage')
        phases[policy]=read(a.compact/(policy+'-phases.json'))
    for policy,l in ledgers.items():
        d=read(l.file)
        require(all(s['status']=='FINISHED' for s in d.get('sessions',{}).values()),'unclosed/unknown policy session')
        require(all(e['status']=='FINISHED' or k=='final-seal-transport-reserve' for k,e in d['entries'].items()),'unreconciled attempt remains charged; cannot seal PASS')
        l.finish('final-seal-transport-reserve',0,time.time()-start+180,str(out/'finalization.json'))
        totals=l.totals();require(totals['attempts']<=800 and totals['seconds']<=10800,'RESOURCE_LIMIT: final ledger')
        write_once(a.compact/('final-ledger-'+policy+'.json'),read(l.file));resource[policy]['totals']=totals
    barriers={f.stem:read(f) for f in (root/'joint').glob('*.json')}
    terminal=barriers.get('barrier-128',component('NOT_RUN','T3 unlaunched/unreached'))
    first_limit=next((dict(name=k,gate=v) for k,v in sorted(barriers.items()) if v['status']!='PASS'),None)
    combined={p:dict(wj6_seconds=5030.391243219376 if p=='EP' else 5875.989612579346,
        wj7_seconds=resource[p]['totals']['seconds'],combined_seconds=(5030.391243219376 if p=='EP' else 5875.989612579346)+resource[p]['totals']['seconds']) for p in ledgers}
    final=dict(schema='WJ7FinalSealV1',primary_target=terminal,first_limiting_gate=first_limit,phases=phases,
        counts=counts,resources=resource,combined_historical_costs=combined,
        parent_outcome='WJ6 Tx PASS; T3 RESOURCE_LIMIT/NOT_RUN, unchanged',closure_transport_reserve_seconds=180,
        closed_ledgers=artifacts([l.file for l in ledgers.values()]),code=file_identity(__file__),governing_attempts=0)
    write_once(a.compact/'final-seal.json',final)
    with (a.compact/'final-command-receipts.json').open('xb') as stream:stream.write((root/'command-receipts.json').read_bytes())
    write_once(a.compact/'publication.json',dict(inner_receipt=file_identity(a.compact/'receipt.json'),final_seal=file_identity(a.compact/'final-seal.json'),
        command_receipts=file_identity(a.compact/'final-command-receipts.json'),figures=artifacts(a.figures.rglob('*')),
        files=artifacts(f for f in a.compact.rglob('*') if f.suffix!='.lock'),parent_transport=file_identity(a.parent/'export-v1/transport.json')))
    files={}
    def export(f,name):
        if f.suffix=='.lock' or not f.is_file():return
        raw=f.read_bytes();original=file_identity(f);encoding='exact bytes'
        if len(raw)>180000 and f.suffix in ('.json','.csv'):
            raw=gzip.compress(raw,compresslevel=9,mtime=0);name+='.gz';encoding='gzip mtime0'
        target=out/name;target.parent.mkdir(parents=True,exist_ok=True)
        with target.open('xb') as stream:stream.write(raw)
        files[name]=dict(original=original,encoding=encoding,transport_bytes=len(raw),transport_sha256=hashlib.sha256(raw).hexdigest(),uncompressed_sha256=original['sha256'])
    for f in a.compact.rglob('*'):export(f,str(f.relative_to(a.compact)))
    for f in a.figures.rglob('*'):export(f,'figures/'+str(f.relative_to(a.figures)))
    for f in (root/'operations').rglob('*'):
        if f.is_file():export(f,'operations/'+str(f.relative_to(root/'operations')))
    for policy in ('EP','Budget'):
        study=root/'study-v2'/policy
        for f in study.rglob('*'):
            if f.name in ('case.xml','requested-intent.json','xml-identity.json','resource-design.json','design.json','execution-intent.json','receipt.json'):
                export(f,'intents/'+policy+'/'+str(f.relative_to(study)))
    for f in (root/'prefix-row-cache-index-EP.json',root/'prefix-row-cache-index-Budget.json'):export(f,'identity/'+f.name)
    for f in root.glob('production-profile/*.json'):export(f,'production-profile/'+f.name)
    write_once(out/'transport.json',dict(schema='WJ7TransportV1',files=files,code=file_identity(__file__),
        parent_transport=file_identity(a.parent/'export-v1/transport.json'),governing_attempts=0,
        verification_scope='exported encoded bytes and original decoded hashes; heavy index hashes verified remotely, heavy fields not transported'))
    print(json.dumps(dict(output=str(out),files=len(files),seconds=time.time()-start,counts=counts,resources=resource,terminal=terminal['status'])))


if __name__=='__main__':main()
