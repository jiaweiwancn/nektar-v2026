#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Transport compact R2/C evidence; native payloads stay at indexed Linux paths."""
import argparse,gzip,hashlib,json,tarfile
from pathlib import Path


def load(p):return json.loads(Path(p).read_text())


def identity(p):
    p=Path(p);h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return dict(path=str(p),bytes=p.stat().st_size,sha256=h.hexdigest())


def publish(root,tag):
    base=root/'resumption-r1';out=base/('publication-'+tag);out.mkdir()
    files=[];seen=set()
    def copy(p):
        p=Path(p)
        if p in seen:return
        seen.add(p);b=p.read_bytes();name=str(p.relative_to(base)).replace('/','__')
        packed=gzip.compress(b,mtime=0)if len(b)>60000 else b
        if len(b)>60000:name+='.gz'
        (out/name).write_bytes(packed)
        files.append(dict(source=identity(p),local=name,stored_bytes=len(packed),stored_sha256=hashlib.sha256(packed).hexdigest()))
    operations={}
    for d in sorted((base/'operations').iterdir()):
        if not(d/'receipt.json').exists():continue
        r=load(d/'receipt.json');operations[d.name]={k:r.get(k)for k in ('status','seconds','phase','reason','reserved_attempts','actual_attempts','peak_tree_rss_bytes')}
        for n in ('request.json','receipt.json','native-exec-maps.json','attempt-counts.json','timing.jsonl'):
            if(d/n).exists():copy(d/n)
    for n in ('accounting.json','tools-manifest-v4.json','tools-manifest-v5.json','source-manifest-v1.json','source-manifest-v2.json','source-manifest-v3.json'):
        if(base/n).exists():copy(base/n)
    compact={'qualification.json','execution.json','execution-intent.json','intent.json','request.json','launch-result.json','outcome.json','design.json','initial-source.json','rotation.json','decision.json','case.xml','temporal-qualification.json','native.log.retention.json'}
    large={'comparison.json','temporal-10.json','temporal-12.json'}
    scientific={};heavy=[]
    for p in sorted((base/'science').rglob('*')):
        if not p.is_file():continue
        rel=p.relative_to(base/'science');n=p.name
        if n in compact or n in large or (len(rel.parts)==1 and n.endswith('.json')):
            copy(p)
        elif p.parent.name in ('memo-component-v1','positive-prototype-v1')and n=='native.json':copy(p)
        elif n in ('independent.json',)and p.parent.name=='fresh-C-v1':copy(p)
        elif p.parent.name=='monitor'and n.endswith(('.qualification-commit.json','.completion.json','.inventory.json')):copy(p)
        if n=='qualification.json':scientific[str(rel.parent)]=load(p)
        if n in ('State.bin','case.evolution.json','native.json')or p.suffix=='.npz':heavy.append(identity(p))
    for n in ('original-account-before-close.json','continuation-account-before-close.json','closure-accounts.json','final-validation.json'):
        if(base/n).exists():copy(base/n)
    costfile=base/'science/complete-cost-v1.json'
    if not costfile.exists():costfile=base/'science/provisional-cost-r2-v1.json'
    costs=load(costfile)if costfile.exists()else None
    summary=dict(schema='WJH25R1ClosureTransportV1',tag=tag,
        status='TARGET_NOT_MET'if tag.startswith('c-')else'IMPLEMENTATION_GAIN_WITH_TRANSIENT_QUALIFICATION_OPEN',
        operations=operations,scientific=scientific,cost=costs,
        heavy_data=heavy,root=str(base),
        limitations=['Prescribed finite perturbation fails frozen initial principal quality; it was not advanced.',
            'Unperturbed h/q initial screens and short C controls do not replace transient h/q/dt or production-cadence qualification.',
            'No full-duration flat or hill validation was executed.'])
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (out/'transport-manifest.json').write_text(json.dumps(dict(schema='WJH25R1TransportV1',tag=tag,root=str(base),files=files,summary=identity(out/'summary.json')),indent=2)+'\n')
    bundle=base/('publication-'+tag+'.tar.gz')
    with tarfile.open(bundle,'x:gz')as t:
        for p in sorted(out.iterdir()):t.add(p,arcname=p.name)
    print(json.dumps(dict(bundle=identity(bundle),summary=identity(out/'summary.json'),files=len(files))))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--tag',choices=('r2-v1','c-v1'),required=True)
    a=p.parse_args();publish(a.root,a.tag)
