#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Select a local h successor from saved failed cells, never relax a gate."""
import argparse
import os
from pathlib import Path
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
import numpy as np
from wj_h15_observer import identity,atomic
from study_wj_h14 import stored_reference
from wj_h19_analysis import GATE,LIMITS,principal_screen


def run(root):
    d=root/'resumption-r1/science/h350-r1-v1'
    arrays=d/'common-signed-fields.npz';a=np.load(arrays)
    design=load(d/'design.json');c=stored_reference()
    p=a['physical'];r=p[:,2]-c.d0;outer=r>c.H
    # All one-sided physical owners are retained; geometry determines which
    # coarse parent interval carries a face value, not a nearest-height guess.
    coarse=load(d/'case.evolution.json');mesh={e['geometry']:e for e in coarse['fields'][0]['mesh']}
    z=design['parent_knots'];ids=a['coarse_geometry']
    pairs=np.asarray([z.index(mesh[int(g)]['lower'][2]) for g in ids])
    gradient=c.u_tau**2/(c.nu+c.w['a']*np.maximum(r-c.w['zv'],0.))
    gradient[outer]=c.u_tau/(.3819178635467331*r[outer])
    scales=dict(U=np.ones(len(p)),k=np.full(len(p),c.K),R=np.full(len(p),c.u_tau**2),
        pi=np.full(len(p),c.u_tau**2),G=np.maximum(abs(gradient),c.u_tau/19.9125),epsilon=c.A/np.maximum(r,c.H))
    errors={}
    for name,limit in LIMITS.items():
        v=a['coarse_'+name]-a['fine_'+name]
        normalized=np.linalg.norm(v.reshape(len(p),-1),axis=1)/scales[name]/limit
        if name=='epsilon':normalized[~outer]=0.
        errors[name]=normalized
    budget=.1*(a['fine_margin']-GATE)
    quality=np.zeros(len(p));quality[outer]=abs(a['coarse_margin'][outer]-a['fine_margin'][outer])/budget[outer]
    errors['principal_minimized']=quality
    records=[]
    for pair in range(25):
        mask=pairs==pair
        maxima={n:float(np.max(v[mask])) for n,v in errors.items()}
        records.append(dict(pair=pair,chart_interval=z[pair:pair+2],utilization=maxima,restore_fine=max(maxima.values())>1))
    restored=[v['pair'] for v in records if v['restore_fine']]
    # Current comparisons already demand retaining every inner interval on
    # stress/velocity/gradient grounds. Global CG projection can couple adjacent
    # intervals; this is a successor proposal requiring fresh qualification.
    proposed=sorted(set(restored)|set(range(8)))
    atomic(root/'resumption-r1/science/h350-diagnosis-r1-v1.json',dict(status='SUCCESSOR_PROPOSAL_NOT_ACCEPTANCE',
        source=identity(arrays),comparison=identity(d/'comparison.json'),records=records,
        restore_pairs=proposed,cells=14*(25+len(proposed)),
        rationale='restore both original children in each pair failing a frozen pointwise budget; preserve fine inner layer',
        limitation='principal minimized action is a localization diagnostic; all directions and complete support remain mandatory in the new comparison'))
    print('restore pairs',proposed,'cells',14*(25+len(proposed)),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);run(p.parse_args().root)
