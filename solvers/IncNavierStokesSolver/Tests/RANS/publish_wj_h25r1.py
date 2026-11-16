#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact immutable checkpoint transport; heavy arrays remain indexed on Linux."""
import argparse,gzip,hashlib,json,os,tarfile
from pathlib import Path


def identity(p):
    p=Path(p);b=p.read_bytes()
    return dict(path=str(p),bytes=len(b),sha256=hashlib.sha256(b).hexdigest())


def load(p):return json.loads(Path(p).read_text())


def publish(root,tag):
    base=root/'resumption-r1';out=base/('publication-'+tag);out.mkdir()
    operations=[p for p in sorted((base/'operations').iterdir()) if (p/'receipt.json').exists()]
    files=[]
    def copy(p,name=None):
        p=Path(p);b=p.read_bytes();name=name or str(p.relative_to(base)).replace('/','__')
        stored=gzip.compress(b,mtime=0) if len(b)>60000 else b
        if len(b)>60000:name+='.gz'
        dest=out/name;dest.write_bytes(stored)
        files.append(dict(source=identity(p),local=name,stored_bytes=len(stored),
            stored_sha256=hashlib.sha256(stored).hexdigest()))
    for d in operations:
        for name in ('request.json','receipt.json','native-exec-maps.json'):
            if (d/name).exists():copy(d/name)
    for name in ('accounting.json','science/runtime-build-seed-v1.json'):
        if (base/name).exists():copy(base/name)
    for d in ('seed700-pert-dt0025-r1-v1','seed700-dt005-r1-v1','seed350-pert-dt0025-r1-v1'):
        for n in ('request.json','execution-intent.json','launch-result.json','independent.json','qualification.json'):
            if (base/'science'/d/n).exists():copy(base/'science'/d/n)
    for name in ('certificate-feasibility-r1-v1.json','h350-diagnosis-r1-v2.json'):
        if (base/'science'/name).exists():copy(base/'science'/name)
    h=load(base/'science/h350-r1-v1/comparison.json')
    f=load(base/'science/seed700-pert-dt0025-r1-v1/independent.json')
    summary=dict(schema='WJH25R1CheckpointR1V1',status='STAGE1_EXECUTED_WITH_REJECTIONS_NO_PRODUCTION_ADMISSION',
        h350=dict(status=h['status'],source=identity(base/'science/h350-r1-v1/comparison.json'),
            numerical_status=h['numerical_screen']['status'],principal=h['numerical_screen']['principal'],
            field_utilization={region:{name:{metric:v[metric]['utilization'] for metric in ('rms','maximum')}
                for name,v in fields.items()}for region,fields in h['numerical_screen']['temporal'].items()},
            physical=h['coarse_physical']['summary'],wall_stress_scaled_max=h['wall_stress_scaled_max']),
        fine_perturbation=dict(status=f['status'],physical=f['physical']['summary'],finite=f['finite_perturbation'],
            signed_Q1_maxima={k:v['maximum']for k,v in f['physical']['all_regional_Q1_moments'].items()},
            source=identity(base/'science/seed700-pert-dt0025-r1-v1/independent.json')),
        own_dt_seed=load(base/'science/seed700-dt005-r1-v1/qualification.json'),
        target350=dict(status='NATIVE_HARD_AUDIT_REJECTED',governing_attempts=0,
            source=identity(base/'science/seed350-pert-dt0025-r1-v1/native.log'),
            distinction='source fine strict import in compatible fine context succeeded; target construction reached native hard audit; no ordinary import was weakened'),
        coefficient=identity(base/'science/certificate-feasibility-r1-v1.json'),
        source_capsule=identity(base/'source-manifest-v1.json'),
        original_account=identity(root/'accounting.json'),continuation_account=identity(base/'accounting.json'),
        operations={p.name:load(p/'receipt.json')for p in operations},
        decision='Keep frozen gates. Do not advance the rejected finite initial state. Continue bounded native assembly/bound prototypes, localized spatial successor, quadrature-only initial screen and retained-trajectory operational qualification.',
        native_source_revision='955600e1569feea0317814f6b15fea4d5d3ce9aa',
        screen_tool_revision='a285c10f',
        future_source='prototype source in Git is separate from this executed source-v1/build-seed-v1/runtime-seed-v1 provenance')
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    manifest=dict(schema='WJH25R1TransportV1',root=str(base),tag=tag,files=files,
        summary=identity(out/'summary.json'))
    (out/'transport-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    bundle=base/('publication-'+tag+'.tar.gz')
    with tarfile.open(bundle,'x:gz')as t:
        for p in sorted(out.iterdir()):t.add(p,arcname=p.name)
    print(json.dumps(dict(bundle=identity(bundle),summary=identity(out/'summary.json'))))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--tag',required=True)
    a=p.parse_args();publish(a.root,a.tag)
