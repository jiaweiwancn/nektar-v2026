#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Configuration-specific sealed milestone audits, including true time coverage."""
import argparse,gzip,json,os,subprocess,sys,time
from pathlib import Path
from qualify_rans_preservation import load
from rans_request_manifest import RequestManifest,sha,runtime_identity
from rans_qualification import THREADS,write_json
from rans_duration_contract import MEMBERS,identity,require,required_times,scalar_coverage,residual_reduction,aggregate,prefix
from admit_rans_configuration import validate_audit
from analyze_rans_duration_constitutive import analyze as constitutive
from schedule_rans_balanced_evolution import workers
from rans_worker_budget import spawn


def finish(r, admission_path, dt, T, resume=False):
    source=Path(__file__).parent;ad=load(admission_path);config=ad['configuration'];tag=prefix(config['model'],config['policy'],config['wall'])+'-duration-'+str(dt);runtime=r/'runtime-resolution-v3';require(runtime.exists(),'forward_probe_runtime_not_frozen')
    os.environ.update(THREADS,LD_LIBRARY_PATH=str(runtime/'lib'),TMPDIR=str(r/'tmp'))
    common=Path(ad['common_directory']);require(sha(common/'packet.json.gz')==ad['common_packet_sha256'],'common_packet_changed');require(runtime_identity(runtime/'RANSPreservationProbe')==ad['duration_diagnostic_runtime'],'admitted_diagnostic_runtime_changed')
    times=[t for t in ad['milestones'] if t<=T+1e-13];out=r/'duration-audits'/(tag+'-'+str(T));labels={m:tag+'-'+m for m in MEMBERS}
    request=dict(schema='ConfigurationMilestoneAuditV2',configuration=config,timesteps=[dt,dt/2],required_times=times,meshes=ad['meshes'],common_packet_sha256=ad['common_packet_sha256'],admission_sha256=sha(admission_path),analysis_runtime=runtime_identity(runtime/'RANSPreservationProbe'),trajectory_requests={m:sha(r/'campaigns'/v/'request-manifest.json') for m,v in labels.items()})
    manifest=RequestManifest(out,request,resume)
    if manifest.completed(out):return load(out/'results.json')
    env=dict(os.environ,**THREADS,LD_LIBRARY_PATH=str(runtime/'lib'),TMPDIR=str(r/'tmp'));execution=[];files=[]
    def run(name,cmd):
        # Read-only probes also consume MPI helper threads; share the same limit.
        with (out/(name+'.log')).open('x') as f:
            proc=None
            while proc is None:
                proc=spawn(r,[sys.executable,'-B',*map(str,cmd)],env=env,stdout=f,stderr=subprocess.STDOUT)
                if proc is None:time.sleep(5)
            code=proc.wait()
        execution.append(dict(name=name,exit_code=code,command=list(map(str,cmd))));write_json(out/'progress.json',dict(status='running',execution=execution))
        require(code==0,'audit_subprocess_failed:'+name)
    try:
        audits={};companions={};stress={}
        for m,label in labels.items():
            folder=r/'campaigns'/label;req=load(folder/'request-manifest.json')['request'];identity(req,config);reports=[]
            for t in times:
                endpoint=folder/('T-'+format(t,'.4f'));require(RequestManifest(folder,req,True).completed(endpoint),'unsealed_companion_endpoint:'+m);reports.append(load(endpoint/'case.evolution.json'))
            if m in MEMBERS[:2]:
                dest=r/'evolved-audits'/(tag+'-'+str(T)+'-'+m)
                if not dest.exists():run('spatial-'+m,[source/'audit_rans_balanced_evolved.py','--root',r,'--analysis-runtime',runtime,'--campaign-root',r,'--label',label,'--output-label',dest.name,'--times',str(T),'--common-directory',common,'--constitutive-samples'])
                current=validate_audit(dest/'results.json',folder,config,[T],common)
                combined=[]
                for t in times:
                    path=r/'evolved-audits'/(tag+'-'+str(t)+'-'+m)/'results.json';a=validate_audit(path,folder,config,[t],common);combined.extend(a['checks']);files.append(path)
                required_times(combined,times);audits[m]=dict(status='pass',configuration=config,checks=combined)
                packet_path=dest/('T-'+format(T,'.4f'))/'packet.json.gz'
                with gzip.open(packet_path,'rt') as f:packet=json.load(f)
                stress[m]=constitutive(packet,config,req['evaluation']['weights']);files.append(packet_path)
            obsdir=folder/('scalar-observations-'+str(T))
            if not obsdir.exists():run('scalars-'+m,[source/'audit_rans_pilot_samples.py','--root',r,'--label',label,'--through-time',str(T),'--output-label',obsdir.name])
            obs=load(obsdir/'results.json');coverage=scalar_coverage(obs,reports,config,T);require(coverage['actual_state_hashes_checked'],'scalar_state_hashes_missing')
            companions[m]=dict(status='pass',configuration=config,request_sha256=sha(folder/'request-manifest.json'),scalar_coverage=coverage,scalar_observer={k:v for k,v in obs.items() if k!='observations'},observer_sha256=sha(obsdir/'results.json'));files.append(obsdir/'results.json')
        qlabel='duration-quadrature-'+tag+'-'+str(T);qpath=r/qlabel/'results.json'
        if not qpath.exists():run('frozen-quadrature',[source/'check_rans_balanced_state_quadrature.py','--root',r,'--analysis-runtime',runtime,'--model',config['model'],'--policy',config['policy'],'--walls',config['wall'],'--label-prefix',tag,'--time',str(T),'--output-label',qlabel])
        q=load(qpath);require(q['status']=='pass' and len(q['checks'])==2,'frozen_quadrature_failed');files.append(qpath)
        result=dict(configuration=config,status='running',companions=companions,spatial_audits=audits,constitutive=stress,quadrature=q,execution=execution)
        gate=aggregate(dict(config,audit_times=times),result);result.update(status=gate['status'],aggregate=gate,evolved_reductions=[dict(time=t,**residual_reduction(*[next(x for x in audits[m]['checks'] if abs(x['time']-t)<1e-13) for m in MEMBERS[:2]],config['model'])) for t in times])
    except Exception as e:result=dict(status='fail',configuration=config,failure=repr(e),execution=execution)
    write_json(out/'results.json',result);manifest.seal(out,['results.json',*map(str,files)]);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--admission',type=Path,required=True);p.add_argument('--dt',type=float,required=True);p.add_argument('--time',type=float,required=True);p.add_argument('--resume',action='store_true');a=p.parse_args();result=finish(a.root.resolve(),a.admission.resolve(),a.dt,a.time,a.resume);print(result['status']);return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
