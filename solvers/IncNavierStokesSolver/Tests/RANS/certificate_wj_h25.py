#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Coefficient-only feasibility screen, never an online acceptance waiver."""
import argparse,hashlib,os,time
from pathlib import Path
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
import numpy as np
from wj_h15_observer import identity,atomic,verify
from wj_h17_observer import native_load
from validate_cls_epawf_m3 import basis

def run(root):
    started=time.monotonic()
    T=np.asarray([[3,0,0,0],[2,1,1,-2],[1,2,1,2],[0,3,0,0]],dtype=np.longdouble)
    x=np.asarray([-1.,-.75,0.,.5,1.]);t=(x+1)/2
    B=np.stack(((1-t)**3,3*t*(1-t)**2,3*t*t*(1-t),t**3),axis=1)
    assert np.max(abs(B@np.asarray(T,dtype=float)/3-basis(3,x)[0]))<64*np.finfo(float).eps
    monitor=root/'science/W0/monitor';ready=load(monitor/'step-20.completion.json')
    assert ready['schema']=='WJNativeCompletionV1'
    report=verify(ready['report']);sp=monitor/'structure.json'
    assert identity(sp)['sha256']==ready['binding']['structure_sha256']
    r=native_load(report);structure=load(sp)
    assert r['failure']=='' and int(r['steps'])==20 and len(r['state'])==24688
    baseline=load(root/'inputs/baseline.json')
    original=native_load(verify(baseline['own_cp0']['report']))
    assert structure['fields']==original['fields'] and structure['offsets']==original['offsets']
    assert r['signature']==original['signature']
    assert structure['physical_signature_sha256']==hashlib.sha256(r['signature'].encode()).hexdigest()
    assert ready['binding']['state_sha256']==hashlib.sha256(np.asarray(r['state'],dtype='<f8').tobytes()).hexdigest()
    state=np.asarray(r['state'],dtype=np.longdouble);assert np.all(np.isfinite(state))
    # K0 is the immutable reference scale, not a surrogate for current endpoint K.
    K=.008809015145515192;threshold=np.longdouble(K)*(np.longdouble('1e-8')+4096*np.finfo(float).eps)
    precision=dict(dtype=str(np.dtype(np.longdouble)),mantissa_bits=int(np.finfo(np.longdouble).nmant),
        epsilon=str(np.finfo(np.longdouble).eps),extended=bool(np.finfo(np.longdouble).nmant>52))
    records=[]
    for block in (2,3):
        cells=[]
        for el in structure['fields'][block]['mesh']:
            assert el['degree']==3
            assert len(el['ids'])==len(el['sign'])==16 and set(el['sign'])<={-1,1}
            assert all(0<=i<structure['fields'][block]['size'] for i in el['ids'])
            c=(state[np.asarray(el['ids'],dtype=int)+structure['offsets'][block]]*np.asarray(el['sign'],dtype=np.longdouble)).reshape(4,4)
            coefficients=T@c@T.T/9;absolute=abs(T)@abs(c)@abs(T.T)/9
            reserve=4096*np.finfo(float).eps*(K+absolute)
            lower=np.min(coefficients-reserve);limit=threshold if block==2 else np.longdouble(0.)
            cells.append(dict(geometry=el['geometry'],heuristic_lower=str(lower),limit=str(limit),screen_positive=bool(lower>limit)))
        records.append(dict(block=block,cells=cells,screen_positive=sum(v['screen_positive']for v in cells),total=len(cells),
            minimum=str(min(np.longdouble(v['heuristic_lower']) for v in cells))))
    out=root/'resumption-r1/science/certificate-feasibility-r1-v1.json'
    atomic(out,dict(status='READ_ONLY_FEASIBILITY_NOT_ACCEPTANCE',source=identity(report),
        structure=identity(sp),completion=identity(monitor/'step-20.completion.json'),precision=precision,
        basis='native Modified_A cubic to Bernstein; i fast, physical global orientation signs applied; integer transform/9 in long double',
        K_reference=K,k_floor_relative=1e-8,scope_of_floor='whole-cell positive reference floor is stronger than nonnegative moving-ray acceptance',
        reserve='heuristic 4096*binary64 epsilon*(K0+sum absolute transformed terms); NOT a proved outward enclosure',
        records=records,seconds=time.monotonic()-started,
        moving_ray_runtime_weighted_coverage=dict(status='NOT_MEASURED',fraction=None,
            reason='saved state/structure and aggregate native timings do not resolve moving-ray cost per cell; cell counts are not runtime coverage'),
        decision='Retain every exact guard; no sufficient-positivity prototype is justified by cell counts without weighted cost evidence and an outward proof',
        scope='saved coefficient-only feasibility; no solve/factorization, no guard skipped; current-K/ray/root/ground/trace gates still required'))
    print([{k:v[k]for k in ('block','screen_positive','total','minimum')}for v in records],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True,type=Path);run(p.parse_args().root)
