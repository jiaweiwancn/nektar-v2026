#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Actual shared completed .05/.06 states; unperturbed window only."""
import argparse,os
from pathlib import Path
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
import numpy as np
from wj_h15_observer import identity,verify,atomic,assess
from wj_h17_observer import load_compatible,native_load
from wj_h20_observer import common_points
from wj_h19_analysis import compare


def run(root):
    d=root/'resumption-r1/science/C-wide-v1';assert load(d/'qualification.json')['status']=='PASS_NATIVE_AND_PHYSICAL_PENDING_TEMPORAL'
    cache,_=load_compatible(d/'cache-compatible.json')
    fixed=load(verify(load(root/'inputs/baseline.json')['support']))['points']
    records=[]
    for nc,nm,name in ((10,20,'W1-lookup-v1'),(12,24,'W0')):
        cp=d/f'monitor/step-{nc}.report.json';mp=root/f'science/{name}/monitor/step-{nm}.report.json'
        c=native_load(cp);m=native_load(mp)
        assert c['steps']==nc and m['steps']==nm and c['dt']==.005 and m['dt']==.0025
        assert c['last_stage_alpha']==200 and m['last_stage_alpha']==400
        assert abs(c['time']-m['time'])<=128*np.finfo(float).eps*max(1.,abs(c['time']),abs(m['time']))
        points=common_points(cache,fixed,{'C':c,'M':m},{'C':c,'M':m})
        rows,initial,stats=cache.view(points)
        pc,dc=assess(rows,c['state'],initial,nc);pm,dm=assess(rows,m['state'],initial,nm)
        result,actions=compare(rows,dc,dm,initial)
        wall=[i for i,(p,o)in enumerate(rows.rows[:rows.point_count])if o['chart'][2]==0]
        wall_error=float(np.max(np.linalg.norm(dc['T'][wall,:,2]-dm['T'][wall,:,2],axis=1))/rows.config.u_tau**2)
        status='PASS'if result['status']=='PASS'and pc['summary']['status']==pm['summary']['status']=='PASS'and wall_error<=1e-3 else'REJECTED'
        arrays={prefix+'_'+k:v for prefix,data in (('C',dc),('M',dm))for k,v in data.items()}
        arrays.update(physical=np.asarray([p for p,o in rows.rows]),geometry=np.asarray([o['geometry']for p,o in rows.rows]),weights=rows.weights)
        arrays.update({'action_'+k:v for k,v in actions.items()})
        payload=d/f'temporal-{nc}-signed-fields.npz';np.savez_compressed(payload,**arrays)
        record=dict(status=status,coarse=identity(cp),fine=identity(mp),time_C=c['time'],time_M=m['time'],
            coarse_physical=pc,fine_physical=pm,comparison=result,wall_traction_scaled_MAX=wall_error,
            signed_fields=identity(payload),cache=stats,current_points=len(points),owned_rows=rows.point_count,
            scope='short evolving unperturbed flat control; no substituted finite perturbation, temporal order or production-duration claim')
        atomic(d/f'temporal-{nc}.json',record);records.append(record)
    status='PASS_SHORT_UNPERTURBED_WINDOW'if all(r['status']=='PASS'for r in records)else'REJECTED'
    atomic(d/'temporal-qualification.json',dict(status=status,records=[identity(d/f'temporal-{n}.json')for n in (10,12)],
        earlier_window='H24 C/M/F .02 to .04 reused after new exact coefficient/history parity',
        finite_perturbation='REJECTED_AT_INITIAL_QUALITY_GATE; transient h/q/dt qualification remains open'))
    print(status,flush=True)
    assert status=='PASS_SHORT_UNPERTURBED_WINDOW'


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);run(p.parse_args().root)
