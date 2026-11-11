#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""One predeclared refined half-dt diagnostic; no continuation of failed state."""
import argparse
import json
import os
from pathlib import Path
import numpy as np
from wj_h2_cases import case_xml
from run_wj_h2 import execute,load,save,sha


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--solver',type=Path,required=True);a=p.parse_args()
    root=a.root.resolve();out=root/'diagnosis-v1';out.mkdir(exist_ok=False)
    closed=root/'study-v1/progress.json';study=load(closed)
    assert study['status']=='OBSTRUCTED' and study['first_limiting_gate']=='refined physical pilot'
    coarse=load(study['runs']['refined-pilot']['report']);assert coarse['steps']==2 and coarse['time']==.04
    assert coarse['physical_range_events'][-1]['reasons']==['WJ_principal_in_plane_exit']
    save(out/'intent.json',dict(schema='WJH2RefinedHalfDtDiagnosticV1',parent=str(closed),parent_sha256=sha(closed),
        native_report_sha256=sha(study['runs']['refined-pilot']['report']),solver=str(a.solver),solver_sha256=sha(a.solver),
        dt=.01,maximum_steps=8,target=.08,stop='first hard/correction exit',seconds_limit=300,
        role='timestep diagnosis; original T=8 matched study remains OBSTRUCTED/NOT_RUN'))
    env=dict(LD_LIBRARY_PATH=os.environ['LD_LIBRARY_PATH'])
    zero=case_xml(out/'cp0','physical',0,.01,True,'ExactStepV1')
    assert execute(root,'diagnosis-refined-cp0',[a.solver,zero],zero.parent,0,180,extra_env=env)
    cp0=load(zero.with_suffix('.evolution.json'));prior=load(study['runs']['refined-cp0']['report'])
    assert np.asarray(cp0['state'],dtype='<f8').tobytes()==np.asarray(prior['state'],dtype='<f8').tobytes(),'changed physical initialization'
    case=case_xml(out/'half-dt','physical',8,.01,True,'ExactStepV1',restart=cp0['checkpoint'])
    okay=execute(root,'diagnosis-refined-half-dt',[a.solver,case],case.parent,8,300,extra_env=env)
    report=load(case.with_suffix('.evolution.json'))
    result=dict(execution='PASS' if okay else 'HARD_EXIT',qualification='PASS' if okay and report['steps']==8 else 'OBSTRUCTED',
        accepted_steps=report['steps'],time=report['time'],failure=report['failure'],
        xml=str(case),xml_sha256=sha(case),report=str(case.with_suffix('.evolution.json')),
        report_sha256=sha(case.with_suffix('.evolution.json')),original_study_unchanged=sha(closed)==load(out/'intent.json')['parent_sha256'],
        attempt_count=len(report['pimple_attempts']),complete_attempts=sum(v['scheduled_complete'] for v in report['pimple_attempts']),
        bound_activations=sum(v['bound_activations'] for v in report['pimple_attempts']),
        accepted_initial_physical_fields_equal=True,retains_native_dt_labels=True)
    save(out/'result.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':main()
