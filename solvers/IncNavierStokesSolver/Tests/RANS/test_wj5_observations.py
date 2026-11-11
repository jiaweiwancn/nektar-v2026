#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compiled custom-z sizes, all-face ownership and actual constraint maps."""
import argparse
from pathlib import Path
import uuid
import numpy as np
from numpy.testing import assert_allclose
from run_wj5_native import Study
from wj5_measure import constraint_estimate
from wj_immutable import write_once


def main():
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--probe',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--policy',choices=('EP','Budget'),required=True);a=p.parse_args();a.retained=a.output;a.output=a.output/('run-'+uuid.uuid4().hex[:12]);a.solver=a.solver.resolve();a.probe=a.probe.resolve()
    s=Study(a)
    try:
        records=[]
        for mesh,nx,count in (('M1',1,1409),('M1',2,2447),('M2',1,1713),('M2',2,2975)):
            d,r,f=s.native('cp0-'+mesh+'-'+str(nx),mesh=mesh,nx=nx,role='focused_control');s.valid(d,r);assert len(d['state'])==count
            records.append(dict(mesh=mesh,nx=nx,size=count,report=str(f/'case.evolution.json')))
            if mesh=='M1' and nx==2:initial=d;producer=f
        m,c=s.measure('actual-maps',initial,producer,mesh='M1',nx=2,assembly=True)
        corner=[r for r in m['observations'] if r['position']==[.5,0.,.1]]
        assert len(corner)==4 and len({r['geometry'] for r in corner})==4
        assert {(r['requested_x_side'],r['requested_normal_side']) for r in corner}=={(-1,-1),(-1,1),(1,-1),(1,1)}
        rng=np.random.default_rng(541);free=rng.normal(size=len(initial['state']));constraints={int(i):v for i,v in m['constraints']}
        for i in constraints:free[i]=0.
        fixed_count=free_count=pressure_count=dependent_count=0;maximum=0.
        for r in m['observations']:
            x,_,z=r['position'];boundary=x in (0.,1.) or z==2.
            for j,desc in enumerate(r['constraint_rows']):
                if desc is None:continue
                native=sum(float(p)*free[int(i)] for i,p in desc['row']);mapped=sum(float(p)*free[int(i)] for i,p in desc['free_row'])
                # On dependent H epsilon the map includes free K, rather than
                # treating its constrained coefficient as independent data.
                if j==3 and z==1.:
                    dependent_count+=1;continue
                maximum=max(maximum,abs(native-mapped));assert_allclose(native,mapped,rtol=0,atol=3e-13)
                if boundary and j<4:
                    fixed_count+=1;assert desc['free_row_L2']<1e-12
                    estimate=constraint_estimate(desc,1e-3)
                    assert 0<estimate['total']<1e-8
                    assert constraint_estimate(desc,1e-3,prescribed_uncertainty=1e-6)['total']>1e-6
                if not boundary and desc['free_row_L2']>.1:free_count+=1
                if j==4:pressure_count+=1;assert desc['free_row_L2']>.1
        assert min(fixed_count,free_count,pressure_count,dependent_count)>0
        write_once(s.root/'focused-result.json',dict(status='PASS',meshes=records,corner=corner,maximum_free_action_defect=maximum,
            fixed_rows=fixed_count,free_rows=free_count,pressure_rows=pressure_count,dependent_rows=dependent_count,
            diagnosis=c['diagnosis'],independent_loads=c['independent_loads'],governing_diagnostic_solves=0))
        print(s.root/'focused-result.json')
    finally:s.ledger.finish_session(s.session);s.save()

if __name__=='__main__':main()
