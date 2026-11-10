#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Four companions, immutable configuration admission, gates before continuation."""
import argparse,os,subprocess,sys,time
from pathlib import Path
from qualify_rans_preservation import load
from rans_qualification import THREADS,write_json
from rans_request_manifest import RequestManifest,sha,runtime_identity
from rans_duration_contract import require,prefix,MEMBERS
from compare_rans_balanced_evolution import compare
from schedule_rans_balanced_evolution import workers
from rans_worker_budget import spawn


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--admission',type=Path,required=True);p.add_argument('--dt',type=float,choices=(.002,.001),default=.002);p.add_argument('--max-time',type=float,choices=(.016,.1,.5),default=.016);p.add_argument('--concurrency',type=int,choices=(1,2,4),default=2);p.add_argument('--resume',action='store_true');a=p.parse_args();r=a.root.resolve();ad=load(a.admission);require(ad['status']=='pass','admission_failed');cfg=ad['configuration'];require([a.dt,a.dt/2] in ad['timestep_trials'],'timestep_trial_not_admitted');tag=prefix(cfg['model'],cfg['policy'],cfg['wall'])+'-duration-'+str(a.dt);source=Path(__file__).parent
    os.environ.update(THREADS,LD_LIBRARY_PATH=str(r/'runtime/lib'),TMPDIR=str(r/'tmp'))
    request=dict(schema='ConfigurationDurationControllerV2',configuration=cfg,admission_sha256=sha(a.admission),timesteps=[a.dt,a.dt/2],milestones=ad['milestones'],meshes=ad['meshes'],common_packet_sha256=ad['common_packet_sha256'],runtime={n:runtime_identity(r/'runtime'/n) for n in ad['runtime']},history=ad['history'],gates=ad['neutral_gates'])
    for name,value in request['runtime'].items():require(value==ad['runtime'][name],'admitted_runtime_changed:'+name)
    cohort=RequestManifest(r/'duration-cohorts'/tag,request,a.resume);out=r/'duration-schedules'/(tag+'-'+str(a.max_time));out.mkdir(parents=True,exist_ok=a.resume);env=dict(os.environ,**THREADS);events=[]
    for T in ad['milestones']:
        if T>a.max_time+1e-13:break
        labels=[tag+'-'+m for m in MEMBERS];queue=[];running={};done=[]
        for member,label in zip(MEMBERS,labels,strict=True):
            mesh='B5' if member.startswith('B5') else 'B11';dt=a.dt/(2 if member.endswith('dt2') else 1)
            cmd=[sys.executable,'-B',str(source/'run_rans_neutral_pilot.py'),'--root',str(r),'--label',label,'--model',cfg['model'],'--policy',cfg['policy'],'--wall',cfg['wall'],'--mesh-file',ad['meshes'][mesh]['path'],'--dt',str(dt),'--milestones',*map(str,ad['milestones']),'--max-time',str(T),'--timeout','14400']
            require(sha(ad['meshes'][mesh]['path'])==ad['meshes'][mesh]['sha256'],'admitted_mesh_changed')
            if (r/'campaigns'/label).exists():cmd+=['--resume']
            queue.append((label,cmd))
        while queue or running:
            for label,(proc,stream,cmd) in list(running.items()):
                code=proc.poll()
                if code is None:continue
                stream.close();done.append(dict(label=label,exit_code=code,command=cmd));del running[label]
            if queue and len(running)<a.concurrency:
                label,cmd=queue[0];stream=(out/(label+'-'+str(T)+'-'+str(time.time_ns())+'.log')).open('x');proc=spawn(r,cmd,stdout=stream,stderr=subprocess.STDOUT,env=env)
                if proc is None:stream.close()
                else:queue.pop(0);running[label]=(proc,stream,cmd)
            write_json(out/'progress.json',dict(status='running',time=T,workers=workers(),running=list(running),queued=[v[0] for v in queue],completed=done))
            if queue or running:time.sleep(5)
        dest=r/'comparisons'/(tag+'-'+str(T)+'.json');dest.parent.mkdir(exist_ok=True)
        try:
            result=compare(r,labels,T)
            if dest.exists():require(load(dest)==result,'completed_comparison_changed')
            else:write_json(dest,result)
        except Exception as e:
            result=dict(status='fail',failure=repr(e))
            if not dest.exists():write_json(dest,result)
        field_status=result['status'];audit_status='not_attempted';audit_file=r/'duration-audits'/(tag+'-'+str(T))/'results.json'
        if field_status=='pass' and all(x['exit_code']==0 for x in done):
            cmd=[sys.executable,'-B',str(source/'finish_rans_balanced_duration.py'),'--root',str(r),'--admission',str(a.admission),'--dt',str(a.dt),'--time',str(T)]
            if audit_file.parent.exists():cmd+=['--resume']
            with (out/('audit-'+str(T)+'-'+str(time.time_ns())+'.log')).open('x') as f:code=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=env).returncode
            audit_status=load(audit_file)['status'] if audit_file.exists() else 'missing';require(code==0 or audit_status!='pass','audit_process_status_disagrees')
        status='pass' if field_status=='pass' and audit_status=='pass' and all(x['exit_code']==0 for x in done) else 'failed_or_unresolved'
        events.append(dict(time=T,status=status,process_completion='pass' if all(x['exit_code']==0 for x in done) else 'fail',field_and_temporal=field_status,spatial_scalar_quadrature_forward=audit_status,comparison=str(dest),audit=str(audit_file),execution=done));write_json(out/'results.json',dict(status=status,configuration=cfg,cohort_request_sha256=cohort.payload['sha256'],milestones=events));print(T,status,flush=True)
        if status!='pass':return 1
    return 0


if __name__=='__main__':raise SystemExit(main())
