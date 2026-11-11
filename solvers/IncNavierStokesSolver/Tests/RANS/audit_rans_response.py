#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Gather independent response audits; retain failures and return nonzero."""
import argparse,os,re,subprocess,sys,time
from pathlib import Path
from qualify_rans_preservation import load
from rans_qualification import THREADS,write_json
from rans_request_manifest import RequestManifest,sha,runtime_identity
from rans_duration_contract import MEMBERS,identity,require,scalar_coverage
from admit_rans_configuration import validate_audit
from rans_worker_budget import spawn


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--admission',type=Path,required=True);p.add_argument('--cohort-prefix',required=True)
    p.add_argument('--time',type=float,required=True);p.add_argument('--resume',action='store_true')
    p.add_argument('--receipt-tag');p.add_argument('--supplement-from',type=Path);a=p.parse_args()
    r=a.root.resolve();ad=load(a.admission);cfg=ad['configuration'];source=Path(__file__).parent
    suffix=''
    if a.receipt_tag:
        require(re.fullmatch('[A-Za-z0-9_-]+',a.receipt_tag) is not None and a.supplement_from,'invalid_response_supplement')
        previous=load(a.supplement_from);identity(previous['configuration'],cfg)
        require(previous['status']!='pass' and abs(previous['time']-a.time)<1e-13,'supplement_requires_matching_preserved_failure')
        suffix='-'+a.receipt_tag
    else:require(a.supplement_from is None,'supplement_requires_new_receipt_tag')
    out=r/'response-audits'/(a.cohort_prefix+'-'+str(a.time)+suffix)
    env=dict(os.environ,**THREADS);checks=[];spatial={};scalars={};failures=[];files=['results.json']
    labels=[a.cohort_prefix+'-'+m for m in MEMBERS]
    request=dict(schema='ResponseFinalAuditV2',configuration=cfg,time=a.time,
        admission_sha256=sha(a.admission),meshes=ad['meshes'],common_packet_sha256=ad['common_packet_sha256'],
        diagnostic_runtime=ad['duration_diagnostic_runtime'],policy='GatherIndependentChecksV1',
        trajectory_requests={m:sha(r/'campaigns'/label/'request-manifest.json') for m,label in zip(MEMBERS,labels,strict=True)},
        previous_failure=dict(path=str(a.supplement_from),sha256=sha(a.supplement_from)) if a.supplement_from else None)
    if a.supplement_from:
        previous_request=load(a.supplement_from.parent/'request-manifest.json')['request']
        for key in ('trajectory_requests','common_packet_sha256','diagnostic_runtime','admission_sha256'):
            require(previous_request[key]==request[key],'response_supplement_input_mismatch:'+key)
    seal=RequestManifest(out,request,a.resume)
    if seal.completed(out):return 0 if load(out/'results.json')['status']=='pass' else 1
    os.environ['LD_LIBRARY_PATH']=str(r/'runtime-resolution-v3/lib');env['LD_LIBRARY_PATH']=os.environ['LD_LIBRARY_PATH']
    require(runtime_identity(r/'runtime-resolution-v3/RANSPreservationProbe')==ad['duration_diagnostic_runtime'],'admitted_response_diagnostic_runtime_changed')
    require(sha(Path(ad['common_directory'])/'packet.json.gz')==ad['common_packet_sha256'],'response_common_packet_changed')
    def run(name,cmd):
        with (out/(name+'.log')).open('x') as f:
            proc=None
            while proc is None:
                proc=spawn(r,[sys.executable,'-B',*map(str,cmd)],env=env,stdout=f,stderr=subprocess.STDOUT)
                if proc is None:time.sleep(5)
            code=proc.wait()
        checks.append(dict(name=name,exit_code=code,command=list(map(str,cmd))))
        if code:failures.append('response_audit_failed:'+name)
        write_json(out/'progress.json',dict(status='running',execution=checks,failures=failures))
    for m,label in zip(MEMBERS,labels,strict=True):
        folder=r/'campaigns'/label
        try:
            q=load(folder/'request-manifest.json')['request'];identity(q,cfg);manifest=RequestManifest(folder,q,True);reports=[]
            for T in q['horizons']:
                if T>a.time+1e-13:continue
                endpoint=folder/('T-'+format(T,'.4f'));require(manifest.completed(endpoint),'response_observation_endpoint_missing')
                reports.append(load(endpoint/'case.evolution.json'))
        except Exception as error:
            failures.append(m+':'+repr(error));continue
        if m in MEMBERS[:2]:
            try:
                dest=r/'evolved-audits'/('response-'+a.cohort_prefix+'-'+str(a.time)+'-'+m)
                if not dest.exists():run('spatial-'+m,[source/'audit_rans_balanced_evolved.py','--root',r,'--analysis-runtime',r/'runtime-resolution-v3','--campaign-root',r,'--label',label,'--output-label',dest.name,'--times',str(a.time),'--common-directory',ad['common_directory']])
                spatial[m]=validate_audit(dest/'results.json',folder,cfg,[a.time],Path(ad['common_directory']),allow_failed=True,expected_runtime=ad['duration_diagnostic_runtime'])
                if spatial[m]['status']!='pass':failures.append('spatial_numerical_check_failed:'+m)
                files.extend([str(dest/'results.json'),str(dest/'immutable-request.json'),str(dest/('T-'+format(a.time,'.4f'))/'packet.json.gz')])
            except Exception as error:failures.append('spatial-'+m+':'+repr(error))
        try:
            leaf=folder/('response-scalar-observations-'+str(a.time))
            if not leaf.exists():run('scalars-'+m,[source/'audit_rans_pilot_samples.py','--root',r,'--label',label,'--through-time',str(a.time),'--output-label',leaf.name])
            obs=load(leaf/'results.json')
            scalars[m]=dict(observer_sha256=sha(leaf/'results.json'),coverage=dict(status='not_verified'),receipt={k:v for k,v in obs.items() if k!='observations'})
            files.extend([str(leaf/'results.json'),str(leaf/'native-basis.json.gz')])
            scalars[m]['coverage']=scalar_coverage(obs,reports,cfg,a.time)
        except Exception as error:failures.append('scalars-'+m+':'+repr(error))
    quad={}
    try:
        qlabel='response-quadrature-'+a.cohort_prefix+'-'+str(a.time);path=r/qlabel/'results.json'
        if not path.exists():run('quadrature',[source/'check_rans_balanced_state_quadrature.py','--root',r,'--analysis-runtime',r/'runtime-resolution-v3','--model',cfg['model'],'--policy',cfg['policy'],'--walls',cfg['wall'],'--label-prefix',a.cohort_prefix,'--time',str(a.time),'--output-label',qlabel]+(['--space-names',*ad['space_names']] if 'space_names' in ad else []))
        quad=load(path);files.append(str(path));require(quad['status']=='pass','response_quadrature_failed')
    except Exception as error:failures.append('quadrature:'+repr(error))
    require(set(spatial).issubset(MEMBERS[:2]) and set(scalars).issubset(MEMBERS),'response_audit_member_identity')
    complete=set(spatial)==set(MEMBERS[:2]) and set(scalars)==set(MEMBERS) and bool(quad)
    if not complete:failures.append('required_independent_audit_coverage_incomplete')
    result=dict(status='pass' if not failures and complete else 'fail',execution_status='completed',configuration=cfg,
        time=a.time,spatial=spatial,scalars=scalars,quadrature=quad,execution=checks,failures=failures,
        independent_coverage='complete' if complete else 'incomplete',
        neutral_residual_reduction='not_applicable_to_departure_from_equilibrium',
        scope='actual native/common alpha-zero loads, scalar coverage and frozen quadrature; failures do not suppress independent diagnostics or qualify the trajectory')
    write_json(out/'results.json',result);seal.seal(out,files)
    print(result['status'],result['independent_coverage']);return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
