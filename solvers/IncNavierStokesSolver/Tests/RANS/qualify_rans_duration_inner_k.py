#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Seal the inner-only common-k checks alongside the unchanged global-k norms."""
import argparse
import gzip
import json
import math
import os
from pathlib import Path
os.environ.update(dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS',
    'MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'),'1'))
from rans_duration_contract import require, identity, required_times
from rans_request_manifest import RequestManifest, sha
from rans_preservation_oracle import Reference
from audit_rans_common_scalar_partitions import partitions
from rans_qualification import write_json


def reduction(coarse,fine):
    require(all(math.isfinite(v) and v>=0 for v in (coarse,fine)),'invalid_common_inner_norm')
    status='unresolved_floor' if max(coarse,fine)<=1e-8 else (
        'pass' if fine<=.85*coarse and coarse-fine>5e-10 else 'fail')
    return dict(status=status,coarse=coarse,fine=fine,ratio=fine/coarse if coarse else None,
                maximum_ratio=.85,residual_floor=1e-8,uncertainty=1e-10)


def qualify(root,admission_path,cohort,times,resume=False):
    ad=json.loads(admission_path.read_text());cfg=ad['configuration'];common=Path(ad['common_directory'])/'packet.json.gz'
    require(sha(common)==ad['common_packet_sha256'],'common_packet_changed')
    source={str(t):{m:root/'evolved-audits'/(cohort+'-'+str(t)+'-'+m)/'results.json'
        for m in ('B5','B11')} for t in times}
    request=dict(schema='DurationCommonInnerKQualificationV1',configuration=cfg,
        qualifier_source_sha256=sha(Path(__file__)),
        admission_sha256=sha(admission_path),common_packet_sha256=sha(common),required_times=times,
        source_audit_sha256={t:{m:sha(p) for m,p in paths.items()} for t,paths in source.items()},
        partition_source_sha256=sha(Path(__file__).with_name('audit_rans_common_scalar_partitions.py')),
        gates=dict(maximum_ratio=.85,residual_floor=1e-8,uncertainty=1e-10),
        scope='additional inner-only fixed-common free tests; global-k and full-native norms remain separate')
    out=root/'common-inner-qualification'/cohort;existed=out.exists();seal=RequestManifest(out,request,existed)
    if seal.completed(out):return json.loads((out/'results.json').read_text())
    require(not existed or resume,'incomplete_inner_qualification_requires_explicit_resume')
    with gzip.open(common,'rt') as f:cp=json.load(f)
    ref=Reference(cfg['model'],cfg['policy'],cfg['wall']);checks=[];short=[]
    for t in (.0016,.016):
        groups={}
        for m in ('B5','B11'):
            rows=[v for v in ad['short_spatial_audits'][m]['checks'] if abs(v['time']-t)<1e-13]
            require(len(rows)==1,'missing_required_short_time')
            if 'configuration' in rows[0]:identity(rows[0]['configuration'],cfg)
            groups[m]=partitions(cp,rows[0],ref)
        gate=reduction(*[groups[m]['inner']['mass_dual'] for m in ('B5','B11')])
        short.append(dict(time=t,status='pass' if gate['status'] in ('pass','unresolved_floor') else 'fail',reduction=gate,groups=groups))
    for t in times:
        groups={}
        for m,path in source[str(t)].items():
            audit=json.loads(path.read_text());require(audit['status']=='pass','source_analysis_failed')
            rows=[v for v in audit['checks'] if abs(v['time']-t)<1e-13];require(len(rows)==1,'missing_required_time')
            identity(rows[0]['configuration'],cfg);groups[m]=partitions(cp,rows[0],ref)
        gate=reduction(*[groups[m]['inner']['mass_dual'] for m in ('B5','B11')])
        checks.append(dict(time=t,status='pass' if gate['status'] in ('pass','unresolved_floor') else 'fail',
                           reduction=gate,groups=groups))
    required_times(checks,times);required_times(short,(.0016,.016))
    result=dict(status='pass' if all(q['status']=='pass' for q in checks+short) else 'fail',
        configuration=cfg,checks=checks,short_checks=short,source_audit_sha256=request['source_audit_sha256'],
        common_packet_sha256=request['common_packet_sha256'],scope=request['scope'])
    write_json(out/'results.json',result);seal.seal(out,['results.json']);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--admission',type=Path,required=True);p.add_argument('--cohort',required=True)
    p.add_argument('--times',type=float,nargs='+',default=[.016,.1,.5]);p.add_argument('--resume',action='store_true');a=p.parse_args()
    result=qualify(a.root,a.admission,a.cohort,a.times,a.resume);print(result['status'])
    return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
