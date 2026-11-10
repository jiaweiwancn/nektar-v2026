#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compiled multi-column immutable readback and composite-rule checks."""
import argparse
from pathlib import Path
import uuid
import numpy as np
from numpy.testing import assert_allclose
from run_wj4_native import Study
from wj4_measure import rule,boundary_rule
from wj_immutable import write_once


def main():
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--policy',choices=('EP','Budget'),required=True);a=p.parse_args()
    a.output=a.output.resolve()/('run-'+uuid.uuid4().hex[:12]);a.solver=a.solver.resolve();a.probe=a.probe.resolve();s=Study(a)
    physical=rule();weighted=np.array([r for r in physical if r[3]>0]);x,z,w=weighted[:,0],weighted[:,2],weighted[:,3]
    boundary=boundary_rule();assert_allclose(sum(boundary.values()),5.,rtol=0,atol=3e-14)
    assert_allclose(sum(w for (x,z),w in boundary.items() if z>1),3.,rtol=0,atol=3e-14)
    assert_allclose([sum(w),w@(x**6),w@(z**6)],[2,2/7,128/7],rtol=2e-14,atol=2e-14)
    z,r,f=s.native('cp0',level=0,nx=2);s.valid(z,r)
    assert len(z['state'])==1109
    m,c=s.measure('cp0',z,f,level=0,nx=2,assembly=True)
    evolved,r,producer=s.native('step',1,z['checkpoint'],level=0,nx=2);s.valid(evolved,r)
    m,c=s.measure('step',evolved,producer,level=0,nx=2,assembly=True)
    # Four distinct native elements meet at the fitted H/vertical-face crossing.
    corner=[r for r in m['observations'] if r['position']==[.5,0.,1.]]
    assert len(corner)==4 and len({r['geometry'] for r in corner})==4
    assert {(r['side'],r['requested_x_side']) for r in corner}=={(-1,-1),(-1,1),(1,-1),(1,1)}
    assert all(r['provenance']=='NativeElementSideV2' for r in corner)
    assert m['completed_stage']['alpha_mass']==1/3.2e-5
    q,qc=s.measure('step-q6',evolved,producer,level=0,nx=2,xorder=6)
    for kind in ('continuum_error','initial_drift'):
        for field,row in c['errors'][kind]['full']['fields'].items():
            assert_allclose(row['L2'],qc['errors'][kind]['full']['fields'][field]['L2'],rtol=1e-8,atol=1e-10)
    write_once(s.root/'focused-result.json',dict(status='PASS',unknowns=1109,corner=corner,independent_loads=c['independent_loads'],
        weak_modes=c['independent_weak_modes'],governing_diagnostic_solves=0,composite_orders=[4,6],actual_completed_stage=True))
    s.ledger.finish_session(s.session);s.save()
    print(s.root/'focused-result.json')

if __name__=='__main__':main()
