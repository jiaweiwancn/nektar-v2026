#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Coefficient-only feasibility screen, never an online acceptance waiver."""
import argparse,os
from pathlib import Path
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
import numpy as np
from wj_h15_observer import identity,atomic
from wj_h17_observer import native_load
from validate_cls_epawf_m3 import basis

def run(root):
    T=np.asarray([[3,0,0,0],[2,1,1,-2],[1,2,1,2],[0,3,0,0]],dtype=np.longdouble)
    x=np.asarray([-1.,-.75,0.,.5,1.]);t=(x+1)/2
    B=np.stack(((1-t)**3,3*t*(1-t)**2,3*t*t*(1-t),t**3),axis=1)
    assert np.max(abs(B@np.asarray(T,dtype=float)/3-basis(3,x)[0]))<64*np.finfo(float).eps
    r=native_load(root/'science/W0/monitor/step-20.report.json');structure=load(root/'science/W0/monitor/structure.json')
    state=np.asarray(r['state'],dtype=np.longdouble);K=.008809015145515192;threshold=K*(1e-8+4096*np.finfo(float).eps)
    records=[]
    for block in (2,3):
        cells=[]
        for el in structure['fields'][block]['mesh']:
            assert el['degree']==3
            c=(state[np.asarray(el['ids'],dtype=int)+structure['offsets'][block]]*np.asarray(el['sign'],dtype=np.longdouble)).reshape(4,4)
            coefficients=T@c@T.T/9;absolute=abs(T)@abs(c)@abs(T.T)/9
            reserve=4096*np.finfo(float).eps*(K+absolute)
            lower=float(np.min(coefficients-reserve));limit=threshold if block==2 else 0.
            cells.append(dict(geometry=el['geometry'],lower=lower,limit=limit,certifiable=lower>limit))
        records.append(dict(block=block,cells=cells,certifiable=sum(v['certifiable']for v in cells),total=len(cells),minimum=min(v['lower']for v in cells)))
    atomic(root/'science/certificate-feasibility-v1.json',dict(status='READ_ONLY_FEASIBILITY_NOT_ACCEPTANCE',source=identity(root/'science/W0/monitor/step-20.report.json'),
        structure=identity(root/'science/W0/monitor/structure.json'),basis='native Modified_A cubic to Bernstein, exact integer transform/9 in long double',
        reserve='4096*binary64 epsilon*(K+sum absolute transformed terms); includes cancellation reserve',records=records,
        scope='saved coefficient-only feasibility; no solve/factorization, no guard skipped; current-K/ray/root/ground/trace gates still required'))
    print([{k:v[k]for k in ('block','certifiable','total','minimum')}for v in records],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True,type=Path);run(p.parse_args().root)
