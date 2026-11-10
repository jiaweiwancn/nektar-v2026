#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Conditional fresh-history duration cohort; every milestone has four companions."""
import argparse,os,subprocess,sys,time
from pathlib import Path
from rans_qualification import THREADS,write_json
from qualify_rans_preservation import load
from schedule_rans_balanced_evolution import workers

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--model',choices=('E94E95','A95','A93'),required=True);p.add_argument('--wall',choices=('EP','Accepted'),required=True);p.add_argument('--dt',type=float,choices=(.004,.002),default=.004);p.add_argument('--max-time',type=float,choices=(.016,.1,.5),default=.016);p.add_argument('--concurrency',type=int,choices=(1,2,4),default=4);p.add_argument('--resume',action='store_true');a=p.parse_args();r=a.root.resolve();source=Path(__file__).parent;env=dict(os.environ,**THREADS)
 admission=load(r/'duration-admission.json');assert admission['status']=='pass' and [a.model,'Natural',a.wall] in admission['configurations']
 prefix=('E' if a.model=='E94E95' else a.model)+'-N-'+a.wall;tag=prefix+'-duration-'+str(a.dt);out=r/'duration-schedules'/(tag+'-'+str(a.max_time));out.mkdir(parents=True,exist_ok=False)
 labels=[tag+'-B5',tag+'-B11',tag+'-B5-dt2',tag+'-B11-dt2'];events=[]
 for T in (.016,.1,.5):
  if T>a.max_time+1e-14:break
  receipt=r/'comparisons'/(tag+'-'+str(T)+'.json')
  if receipt.exists():
   assert a.resume and load(receipt)['status']=='pass';continue
  queue=[];running={};done=[]
  for label,mesh,dt in zip(labels,('Balanced5V1','Balanced11V1','Balanced5V1','Balanced11V1'),(a.dt,a.dt,a.dt/2,a.dt/2),strict=True):
   target=r/'campaigns'/label/('T-'+format(T,'.4f'))
   if (target/'completed-artifacts.json').exists():assert a.resume;done.append(dict(label=label,exit_code=0,reused=True));continue
   cmd=[sys.executable,'-B',str(source/'run_rans_neutral_pilot.py'),'--root',str(r),'--label',label,'--model',a.model,'--wall',a.wall,'--mesh-file',str(r/'meshes'/(prefix)/(mesh+'.json')),'--dt',str(dt),'--milestones','.016','.1','.5','--max-time',str(T),'--timeout','14400']
   if (r/'campaigns'/label).exists():cmd+=['--resume']
   queue.append((label,cmd))
  write_json(out/('request-'+str(T)+'.json'),dict(target=T,commands=[c for n,c in queue],admission=admission))
  while queue or running:
   for name,(proc,f,cmd) in list(running.items()):
    code=proc.poll()
    if code is None:continue
    f.close();done.append(dict(label=name,exit_code=code,command=cmd));del running[name]
   if queue and len(running)<a.concurrency and workers()+10<=72:
    name,cmd=queue.pop(0);f=(out/(name+'-'+str(T)+'.log')).open('x');proc=subprocess.Popen(cmd,stdout=f,stderr=subprocess.STDOUT,env=env);running[name]=(proc,f,cmd)
   write_json(out/'progress.json',dict(status='running',target=T,workers=workers(),running=list(running),queued=[n for n,c in queue],completed=done))
   if queue or running:time.sleep(5)
  cmd=[sys.executable,'-B',str(source/'compare_rans_balanced_evolution.py'),'--root',str(r),'--labels',*labels,'--time',str(T),'--output',str(receipt)]
  with (out/('comparison-'+str(T)+'.log')).open('x') as f:code=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=env).returncode
  events.append(dict(time=T,execution=done,comparison=str(receipt),exit_code=code,status=load(receipt)['status']));write_json(out/'results.json',dict(status='pass' if code==0 and all(x['exit_code']==0 for x in done) else 'failed_or_unresolved',milestones=events))
  if code or any(x['exit_code'] for x in done):return 1
 return 0
if __name__=='__main__':raise SystemExit(main())
