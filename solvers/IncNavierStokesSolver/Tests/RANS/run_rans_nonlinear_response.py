#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Prepare actual initial states, then advance admitted matched response cohorts."""
import argparse,os,subprocess,sys,time
from pathlib import Path
from qualify_rans_preservation import load
from rans_duration_contract import MEMBERS,identity,require,prefix
from rans_request_manifest import RequestManifest,sha
from rans_qualification import THREADS,write_json
from rans_worker_budget import spawn


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--admission',type=Path,required=True);p.add_argument('--accepted-neutral-proof',type=Path,required=True);p.add_argument('--dt',type=float,choices=(.002,.001),default=.002);p.add_argument('--max-time',type=float,choices=(.016,.1),default=.016);p.add_argument('--amplitude',type=float,choices=(2e-5,1e-5),default=2e-5);p.add_argument('--concurrency',type=int,choices=(1,2,4),default=2);p.add_argument('--resume',action='store_true');a=p.parse_args();r=a.root.resolve();source=Path(__file__).parent;ad=load(a.admission);cfg=ad['configuration'];proof=load(a.accepted_neutral_proof);require(proof['status']=='pass','neutral_duration_not_qualified');identity(proof['configuration'],dict(model=cfg['model'],policy=cfg['policy'],wall='Accepted'))
    for m in MEMBERS:require(proof['companions'][m]['scalar_coverage']['last_time']>=.5-1e-13,'neutral_duration_horizon_incomplete')
    tag=prefix(cfg['model'],cfg['policy'],cfg['wall']);baseline=[tag+'-duration-'+str(a.dt)+'-'+m for m in MEMBERS];name=tag+'-response-'+str(a.amplitude)+'-dt-'+str(a.dt);labels=[name+'-'+m for m in MEMBERS];out=r/'response-cohorts'/name
    for label in baseline:
        request=load(r/'campaigns'/label/'request-manifest.json')['request'];identity(request,cfg);require(RequestManifest(r/'campaigns'/label,request,True).completed(r/'campaigns'/label/'T-0.0160'),'response_neutral_baseline_missing')
    request=dict(schema='NonlinearResponseControllerV1',configuration=cfg,amplitude=a.amplitude,timesteps=[a.dt,a.dt/2],milestones=[.016,.1],admission_sha256=sha(a.admission),accepted_neutral_proof_sha256=sha(a.accepted_neutral_proof),meshes=ad['meshes'],runtime=ad['runtime'],baseline_requests=[sha(r/'campaigns'/label/'request-manifest.json') for label in baseline],gates=dict(spatial_absolute_floor=5e-7,spatial_relative=.1,temporal_absolute_floor=1e-7,temporal_relative=.05),observables=['u','w','k','epsilon_outer','Pi','pi_mechanical','Tzz','Txz','Rxx','Ryy','Rzz','Rxz','Qxx','Qyy','Qzz','Qxz','P'],initialization='actual zero-step coefficients/native histories sealed before any advancement')
    manifest=RequestManifest(out,request,a.resume);env=dict(os.environ,**THREADS);events=[]
    def command(i,T,initial=False):
        m=MEMBERS[i];mesh='B5' if m.startswith('B5') else 'B11';label=labels[i];cmd=[sys.executable,'-B',str(source/'run_rans_neutral_pilot.py'),'--root',str(r),'--label',label,'--model',cfg['model'],'--policy',cfg['policy'],'--wall',cfg['wall'],'--mesh-file',ad['meshes'][mesh]['path'],'--dt',str(a.dt/(2 if m.endswith('dt2') else 1)),'--velocity-amplitude',str(a.amplitude),'--milestones','.016','.1','--max-time',str(T),'--timeout','14400']
        if (r/'campaigns'/label).exists():cmd+=['--resume']
        if initial:cmd+=['--initialize-only']
        return cmd
    def batch(commands,phase):
        queued=list(enumerate(commands));running={};done=[]
        while queued or running:
            for i,(proc,stream,cmd) in list(running.items()):
                code=proc.poll()
                if code is None:continue
                stream.close();done.append(dict(member=MEMBERS[i],exit_code=code,command=cmd));del running[i]
            if queued and len(running)<a.concurrency:
                i,cmd=queued[0];f=(out/(phase+'-'+str(i)+'-'+str(time.time_ns())+'.log')).open('x');proc=spawn(r,cmd,env=env,stdout=f,stderr=subprocess.STDOUT)
                if proc is None:f.close()
                else:queued.pop(0);running[i]=(proc,f,cmd)
            write_json(out/'progress.json',dict(status='running',phase=phase,running=[MEMBERS[i] for i in running],queued=[MEMBERS[i] for i,c in queued],completed=done))
            if queued or running:time.sleep(5)
        return done
    prepared=batch([command(i,.016,True) for i in range(4)],'initialize')
    require(all(x['exit_code']==0 for x in prepared),'actual_initialization_failed')
    initial=[load(r/'campaigns'/label/'request-manifest.json')['request'] for label in labels]
    for i,j in ((0,2),(1,3)):
        require(initial[i]['actual_initial']['state_sha256']==initial[j]['actual_initial']['state_sha256'],'same_space_actual_initial_coefficients_differ')
        for key in ('reference_sha256','raw_neutral_sha256','boundary'):require(initial[i][key]==initial[j][key],'same_space_reference_boundary_differ:'+key)
    initial_receipt=dict(status='pass',cases=[dict(label=label,actual_initial=q['actual_initial'],neutral_reference_sha256=q['reference_sha256']) for label,q in zip(labels,initial,strict=True)],execution=prepared)
    initial_path=out/'actual-initial-comparison.json'
    if initial_path.exists():require(load(initial_path)['cases']==initial_receipt['cases'],'sealed_actual_initial_comparison_changed')
    else:write_json(initial_path,initial_receipt)
    for T in (.016,.1):
        if T>a.max_time+1e-13:break
        for label in baseline:
            q=load(r/'campaigns'/label/'request-manifest.json')['request'];require(RequestManifest(r/'campaigns'/label,q,True).completed(r/'campaigns'/label/('T-'+format(T,'.4f'))),'matching_baseline_endpoint_missing')
        execution=batch([command(i,T) for i in range(4)],'T-'+str(T));destination=out/('T-'+str(T));status='failed_execution'
        if all(x['exit_code']==0 for x in execution):
            if destination.exists():require(manifest.completed(destination),'unsealed_response_comparison')
            else:
                cmd=[sys.executable,'-B',str(source/'compare_rans_nonlinear_response.py'),'--root',str(r),'--perturbed',*labels,'--baseline',*baseline,'--time',str(T),'--output-directory',str(destination)];proc=None
                with (out/('compare-'+str(T)+'.log')).open('x') as f:
                    while proc is None:
                        proc=spawn(r,cmd,env=env,stdout=f,stderr=subprocess.STDOUT)
                        if proc is None:time.sleep(5)
                    code=proc.wait()
                require((destination/'results.json').exists(),'response_comparison_missing');manifest.seal(destination,['results.json',*(['response-vectors.npz'] if (destination/'response-vectors.npz').exists() else []),*map(str,destination.glob('*/packet.json.gz'))])
            result=load(destination/'results.json');status=result['status']
            if status=='pass' and not result['eligible_for_extension']:status='unresolved_response_signal'
            if status=='pass':
                audit=r/'response-audits'/(name+'-'+str(T))/'results.json'
                cmd=[sys.executable,'-B',str(source/'audit_rans_response.py'),'--root',str(r),'--admission',str(a.admission),'--cohort-prefix',name,'--time',str(T)]
                if audit.parent.exists():cmd+=['--resume']
                with (out/('audit-'+str(T)+'-'+str(time.time_ns())+'.log')).open('x') as f:
                    code=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
                status='pass' if code==0 and audit.exists() and load(audit)['status']=='pass' else 'failed_response_audit'
        events.append(dict(time=T,status=status,execution=execution,comparison=str(destination/'results.json')));write_json(out/'results.json',dict(status=status,configuration=cfg,milestones=events))
        if status!='pass':return 1
    return 0


if __name__=='__main__':raise SystemExit(main())
