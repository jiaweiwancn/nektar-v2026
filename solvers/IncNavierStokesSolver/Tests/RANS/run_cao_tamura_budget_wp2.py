#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Durable bounded WP2 native work; no hill launch or implicit gate promotion."""
import argparse,os,resource,time,fcntl
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from diagnose_cls_accepted_covariance import run_probe
from run_cls_epawf_picard import require,THREADS,sha
from wp2_neutral_cases import case


def charge(root,event,**work):
    with (root/'.ledger.lock').open('a+') as guard:
        fcntl.flock(guard,fcntl.LOCK_EX)
        return _charge_locked(root,event,**work)


def _charge_locked(root,event,**work):
    origin=load(root/'origin.json');p=root/'aggregate.json'
    d=load(p) if p.exists() else dict(version='WP2AggregateWorkV1',totals={},events=[])
    for k,v in work.items():d['totals'][k]=d['totals'].get(k,0)+v
    d['events'].append(dict(event=event,time=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),work=work))
    atomic(p,d)
    for k,v in dict(spatial_cases=32,boundary_calls=24,physical_attempts=6000,stationary_iterations=120,timestep_trials=16).items():require(d['totals'].get(k,0)<=v,'WP2 aggregate '+k)
    require(time.time()-origin['started_epoch']<43200,'WP2 wall envelope')
    return d


def main():
    p=argparse.ArgumentParser()
    for key in ('root','probe','preflight','h1','output'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--models',nargs='+',choices=('CLS','standard'),default=['CLS','standard'])
    p.add_argument('--level',type=int,default=0);p.add_argument('--order',type=int,default=16)
    p.add_argument('--audit-order',type=int,default=24);p.add_argument('--boundary',action='store_true');p.add_argument('--pout',type=float,default=0.)
    p.add_argument('--projection',choices=('RadialL2V1','RadialEnergyV1'),default='RadialL2V1')
    a=p.parse_args();os.environ.update(THREADS);a.output.mkdir(exist_ok=False);pre=load(a.preflight)
    result=dict(version='WP2NativeBatchV1',status='running',probe=str(a.probe),probe_sha256=sha(a.probe),results={},governing_solves=0,physical_steps=0)
    atomic(a.output/'progress.json',result)
    try:
        for model in a.models:
            charge(a.root,'native '+model+' level'+str(a.level),spatial_cases=1,boundary_calls=int(a.boundary))
            xml=case(a.output/(model+'-case'),a.h1,pre,model,a.level,a.order,pout=a.pout,projection=a.projection)
            packet=dict(operation='wp2_equilibrium',inner_stress_model='Linear',xml_sha256=sha(xml),
                        expected_model='CLSHighRe' if model=='CLS' else 'StandardKEpsilon',order=a.audit_order,boundary_checks=a.boundary)
            d,receipt=run_probe(a.probe,xml,packet,a.output/(model+'-native'))
            receipt['child_peak_RSS_KiB_cumulative_max']=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
            spatial=d['spatial'];definition=load(xml.parent/'wp2-definition.json');require(spatial['size']==definition['estimated_slots'],'native dimensions differ from case estimate')
            result['results'][model]=dict(receipt=receipt,size=spatial['size'],state_sha256=d['state_sha256'],
                boundary_passed=d['boundary']['passed'] if d['boundary'] else None,
                free_block_norms=spatial['unanchored_free_block_norms'],analytical_free_block_norms=spatial['analytical_free_block_norms'])
            atomic(a.output/'progress.json',result)
        result['status']='completed_native_batch'
    except Exception as error:result.update(status='failed',failure=str(error));raise
    finally:atomic(a.output/'progress.json',result)


if __name__=='__main__':main()
