#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Absolute-pressure, equal-time CT6 degree and timestep comparisons."""
import argparse,json,math,os
from pathlib import Path
from rans_qualification import THREADS,write_json
os.environ.update(THREADS)
import numpy as np
from ct_accepted6_audit import common_grid
from ct_accepted_feasibility import parameters
from ct_accepted_open import require
from rans_request_manifest import RequestManifest,sha

NAMES=('U','w','k','outer_epsilon','Pi')


def observations(member,audit):
    r=json.loads((member/'results.json').read_text());a=json.loads((audit/'results.json').read_text())
    require(r['physical_completion']=='pass','incomplete or physically rejected required member')
    p=Path(a['common_observations']['path']);require(sha(p)==a['common_observations']['sha256'],'changed common physical packet')
    packet=json.loads(p.read_text());o=packet['observation']
    require(abs(o['time']-r['target'])<1e-13 and abs(a['time']-r['target'])<1e-13,'required endpoint not reached')
    points,sides,weights,n=common_grid()
    require(o['grid']==points and o['grid_sides']==sides and a['common_observations']['weights']==weights,'common observation identities differ')
    v=np.asarray(o['grid_values'][:n]);truth=np.asarray([q['values'] for q in packet['manufactured']['samples'][:n]])
    require(v.shape==truth.shape==(n,5) and np.all(np.isfinite(v)) and np.all(np.isfinite(truth)),'invalid physical field arrays')
    return r,a,o,v,truth,np.asarray(weights),np.asarray(o['grid_regions'][:n])


def norms(delta,weights,regions,model):
    p=parameters(model);scale=[.2,.2,p['q']*.16,.4**3/(p['kappa']*.5),.1];out={}
    for j,name in enumerate(NAMES):
        mask=regions==1 if j==3 else np.ones(len(weights),bool);w=weights[mask];d=delta[mask,j]
        square=float(w@(d*d));rms=math.sqrt(square/float(w.sum()))
        out[name]=dict(weighted_RMS=rms,true_L2=math.sqrt(square),scale=scale[j],scaled_weighted_RMS=rms/scale[j])
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--member',action='append',type=Path,required=True);p.add_argument('--audit',action='append',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();require(len(a.member)==len(a.audit),'member/audit coverage mismatch')
    records=[observations(x,y) for x,y in zip(a.member,a.audit,strict=True)]
    model=records[0][0]['model'];target=records[0][0]['target']
    require(all(r[0]['model']==model and r[0]['target']==target for r in records),'foreign model or horizon')
    indexed={(r[0]['degree'],r[0]['dt']):r for r in records};require(len(indexed)==len(records),'duplicate member')
    m=RequestManifest(a.output,dict(schema='CT6MMSComparisonsV1',model=model,target=target,
        members=[dict(receipt_sha256=sha(x/'results.json'),audit_sha256=sha(y/'results.json')) for x,y in zip(a.member,a.audit,strict=True)],
        source_sha256=sha(Path(__file__)),pressure='absolute, no mean subtraction',floor=1e-7,
        degree_ratio=.85,first_order_band=[.8,1.2]))
    result=dict(schema='CT6MMSComparisonsV1',model=model,time=target,reference_errors={},temporal={},degree={},statuses={})
    for key,r in indexed.items():
        result['reference_errors'][str(key)]=norms(r[3]-r[4],r[5],r[6],model)
        for other in records:
            require(r[2]['grid']==other[2]['grid'] and r[2]['grid_regions']==other[2]['grid_regions'],'sample regions differ')
    for degree in sorted({k[0] for k in indexed}):
        pair=[indexed.get((degree,t)) for t in (.002,.001)]
        if not all(pair):continue
        require(pair[0][0]['zero']['state_sha256']==pair[1][0]['zero']['state_sha256'],'same-space actual initial coefficients differ')
        d=norms(pair[0][3]-pair[1][3],pair[0][5],pair[0][6],model)
        row=dict(difference_002_001=d,initial_coefficients_identical=True)
        if (degree,.004) in indexed:
            q=indexed[degree,.004];require(q[0]['zero']['state_sha256']==pair[0][0]['zero']['state_sha256'],'coarse temporal actual initialization differs')
            d0=norms(q[3]-pair[0][3],q[5],q[6],model);orders={}
            for name in NAMES:
                x=d0[name]['scaled_weighted_RMS'];y=d[name]['scaled_weighted_RMS']
                resolved=min(x,y)>10e-7
                order=math.log2(x/y) if resolved else None
                orders[name]=dict(status=('pass' if .8<=order<=1.2 else 'fail') if resolved else 'unresolved_floor',
                    order=order,D004_002=x,D002_001=y,signal_floor=1e-6)
            row.update(difference_004_002=d0,orders=orders)
        result['temporal'][str(degree)]=row
    if (3,.001) in indexed and (4,.001) in indexed:
        coarse,fine=indexed[3,.001],indexed[4,.001]
        delta=norms(coarse[3]-fine[3],coarse[5],coarse[6],model)
        ce=result['reference_errors'][str((3,.001))];fe=result['reference_errors'][str((4,.001))]
        for name in NAMES:
            x,y=ce[name]['scaled_weighted_RMS'],fe[name]['scaled_weighted_RMS'];resolved=x>1e-7 and x-y>2e-8
            result['degree'][name]=dict(status=('pass' if y<=.85*x else 'fail') if resolved else 'unresolved_floor',
                coarse_error=x,fine_error=y,ratio=y/x if x>1e-7 else None,ratio_limit=.85,
                spatial_difference=delta[name],floor=1e-7,scope='fixed physical partition P3/P1 versus P4/P2, not a measured convergence order')
    orders=[v for r in result['temporal'].values() for v in r.get('orders',{}).values()]
    result['statuses']['temporal_order']='pass' if any(q['status']=='pass' for q in orders) and all(q['status']!='fail' for q in orders) else 'unresolved' if not any(q['status']=='fail' for q in orders) else 'fail'
    result['statuses']['degree_reduction']='pass' if result['degree'] and all(q['status']!='fail' for q in result['degree'].values()) and any(q['status']=='pass' for q in result['degree'].values()) else 'fail_or_unresolved'
    result['statuses']['physical_budgets']='pass' if all(r[1]['status']=='pass' for r in records) else 'fail'
    result['status']='pass' if all(q=='pass' for q in result['statuses'].values()) else 'fail_or_unresolved'
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json'])
    return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
