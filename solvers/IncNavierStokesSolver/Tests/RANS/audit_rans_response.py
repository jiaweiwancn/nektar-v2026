#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Actual response-state residuals, scalar coverage and quadrature; no steady gate."""
import argparse,os,subprocess,sys,time
from pathlib import Path
from qualify_rans_preservation import load
from rans_qualification import THREADS,write_json
from rans_request_manifest import RequestManifest,sha
from rans_duration_contract import MEMBERS,identity,require,scalar_coverage
from admit_rans_configuration import validate_audit
from rans_worker_budget import spawn


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--admission',type=Path,required=True);p.add_argument('--cohort-prefix',required=True);p.add_argument('--time',type=float,required=True);a=p.parse_args();r=a.root.resolve();ad=load(a.admission);cfg=ad['configuration'];source=Path(__file__).parent;out=r/'response-audits'/(a.cohort_prefix+'-'+str(a.time));out.mkdir(parents=True,exist_ok=False);env=dict(os.environ,**THREADS);checks=[];spatial={};scalars={};labels=[a.cohort_prefix+'-'+m for m in MEMBERS]
    def run(name,cmd):
        with (out/(name+'.log')).open('x') as f:
            proc=None
            while proc is None:
                proc=spawn(r,[sys.executable,'-B',*map(str,cmd)],env=env,stdout=f,stderr=subprocess.STDOUT)
                if proc is None:time.sleep(5)
            code=proc.wait()
        checks.append(dict(name=name,exit_code=code,command=list(map(str,cmd))));write_json(out/'progress.json',dict(status='running',execution=checks));require(code==0,'response_audit_failed:'+name)
    try:
        for m,label in zip(MEMBERS,labels,strict=True):
            folder=r/'campaigns'/label;request=load(folder/'request-manifest.json')['request'];identity(request,cfg);manifest=RequestManifest(folder,request,True);reports=[]
            for T in request['horizons']:
                if T>a.time+1e-13:continue
                endpoint=folder/('T-'+format(T,'.4f'));require(manifest.completed(endpoint),'response_observation_endpoint_missing');reports.append(load(endpoint/'case.evolution.json'))
            if m in MEMBERS[:2]:
                dest=r/'evolved-audits'/('response-'+a.cohort_prefix+'-'+str(a.time)+'-'+m)
                run('spatial-'+m,[source/'audit_rans_balanced_evolved.py','--root',r,'--analysis-runtime',r/'runtime-resolution-v3','--campaign-root',r,'--label',label,'--output-label',dest.name,'--times',str(a.time),'--common-directory',ad['common_directory']])
                spatial[m]=validate_audit(dest/'results.json',folder,cfg,[a.time],Path(ad['common_directory']))
            leaf='response-scalar-observations-'+str(a.time);run('scalars-'+m,[source/'audit_rans_pilot_samples.py','--root',r,'--label',label,'--through-time',str(a.time),'--output-label',leaf]);obs=load(folder/leaf/'results.json');coverage=scalar_coverage(obs,reports,cfg,a.time);scalars[m]=dict(observer_sha256=sha(folder/leaf/'results.json'),coverage=coverage,receipt={k:v for k,v in obs.items() if k!='observations'})
        qlabel='response-quadrature-'+a.cohort_prefix+'-'+str(a.time);run('quadrature',[source/'check_rans_balanced_state_quadrature.py','--root',r,'--analysis-runtime',r/'runtime-resolution-v3','--model',cfg['model'],'--policy',cfg['policy'],'--walls',cfg['wall'],'--label-prefix',a.cohort_prefix,'--time',str(a.time),'--output-label',qlabel]);quad=load(r/qlabel/'results.json');require(quad['status']=='pass','response_quadrature_failed')
        result=dict(status='pass',configuration=cfg,time=a.time,spatial=spatial,scalars=scalars,quadrature=quad,execution=checks,neutral_residual_reduction='not_applicable_to_departure_from_equilibrium',scope='actual alpha-zero native/common residuals and nonsteady response checks; no steady-error reduction imposed')
    except Exception as e:result=dict(status='fail',configuration=cfg,time=a.time,failure=repr(e),execution=checks)
    write_json(out/'results.json',result);return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
