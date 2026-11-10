#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Supplement global common-k norms with inner/outer/shared fixed-test norms."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
os.environ.update(dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
                                'BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'),'1'))
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from rans_duration_contract import require, identity
from rans_request_manifest import sha
from rans_preservation_oracle import Reference
from qualify_rans_preservation import load
from rans_qualification import write_json


def partitions(common, audit, ref):
    observation=common['observation'];mass=np.asarray(common['resolution']['mass'][2])
    support=[set(),set()]
    for q in observation['volume']:support[q['region']].update(q['velocity_k']['ids'])
    shared=support[0]&support[1]
    row=audit['common_test']['fields']['k'];free=np.asarray(row['common_free_ids'],int)
    require(hashlib.sha256(np.asarray(mass[np.ix_(free,free)],dtype='<f8').tobytes()).hexdigest()==row['common_mass_sha256'],
            'common_scalar_mass_or_orientation_mismatch')
    vector=np.asarray(row['restricted_load']);require(vector.shape==free.shape,'common_scalar_load_shape')
    result={}
    for name,region in (('inner',0),('outer',1),('shared',0)):
        mask=np.asarray([i in shared if name=='shared' else i in support[region] and i not in shared for i in free])
        ids=free[mask];x=vector[mask];M=mass[np.ix_(ids,ids)]
        value=float(np.sqrt(max(0.,x@cho_solve(cho_factor(M,lower=True),x)))/(ref.K*ref.us/ref.H*np.sqrt(8.)))
        result[name]=dict(mass_dual=value,common_ids=ids.tolist(),test_count=len(ids))
    require(sum(v['test_count'] for v in result.values())==len(free),'common_scalar_partition_coverage')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--admission',type=Path,required=True)
    p.add_argument('--audits',type=Path,nargs=2,required=True);p.add_argument('--time',type=float,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();ad=load(a.admission);cfg=ad['configuration']
    common=Path(ad['common_directory'])/'packet.json.gz';require(sha(common)==ad['common_packet_sha256'],'common_packet_changed')
    with gzip.open(common,'rt') as f:cp=json.load(f)
    ref=Reference(cfg['model'],cfg['policy'],cfg['wall']);groups=[]
    for path in a.audits:
        audit=load(path);require(audit['status']=='pass','source_spatial_audit_failed')
        rows=[r for r in audit['checks'] if abs(r['time']-a.time)<1e-13];require(len(rows)==1,'missing_required_time')
        identity(rows[0]['configuration'],cfg);groups.append(partitions(cp,rows[0],ref))
    c,f=[v['inner']['mass_dual'] for v in groups]
    status='unresolved_floor' if max(c,f)<=1e-8 else 'pass' if f<=.85*c and c-f>5e-10 else 'failed_or_unresolved'
    result=dict(status='pass' if status in ('pass','unresolved_floor') else 'failed_or_unresolved',configuration=cfg,
        time=a.time,common_packet_sha256=sha(common),source_audit_sha256=[sha(p) for p in a.audits],
        groups=dict(zip(('B5','B11'),groups,strict=True)),inner_reduction=dict(status=status,
        fine_over_coarse=f/c if c else None,maximum_ratio=.85,residual_floor=1e-8,uncertainty=1e-10),
        scope='additional fixed common inner-k test norm, same physical global normalization; shared rows counted once separately; existing full/global norms remain unchanged')
    require(not a.output.exists(),'partition_output_exists');write_json(a.output,result)
    return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
