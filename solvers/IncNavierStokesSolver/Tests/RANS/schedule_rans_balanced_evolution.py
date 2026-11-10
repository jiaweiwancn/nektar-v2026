#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""One bounded stage of four-member comparisons, with a shared worker budget."""
import argparse,json,os,subprocess,sys,time
from pathlib import Path
from rans_qualification import THREADS,write_json
from qualify_rans_preservation import load

def workers():
 total=0
 for line in subprocess.check_output(['ps','-eo','nlwp,comm,args'],text=True).splitlines()[1:]:
  p=line.split(None,2)
  if len(p)==3 and (p[1].startswith(('IncNavier','RANS','Nek','cc1plus','f951','orted','mpiexec')) or (p[1].startswith('python') and 'nektar' in p[2])):total+=int(p[0])
 return total

def labels(model,policy,wall):
 prefix=('E' if model=='E94E95' else model)+'-'+('N' if policy=='Natural' else 'P')+'-'+wall
 return [prefix+'-B5',prefix+'-B11',prefix+'-B5-dt2',prefix+'-B11-dt2']

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--model',choices=('E94E95','A95','A93'),default='E94E95');p.add_argument('--policy',choices=('Natural','Prescribed'),default='Natural');p.add_argument('--walls',nargs='+',choices=('EP','Accepted'),default=['EP','Accepted']);p.add_argument('--time',type=float,choices=(.0016,.016,.1,.5),required=True);p.add_argument('--concurrency',type=int,default=4);a=p.parse_args();assert 1<=a.concurrency<=4
 r=a.root.resolve();source=Path(__file__).parent;out=r/'evolution-schedules'/(a.model+'-'+a.policy+'-'+str(a.time));out.mkdir(parents=True,exist_ok=False);env=dict(os.environ,**THREADS);queue=[];running={};completed=[]
 for wall in a.walls:
  if a.model=='E94E95' and a.policy=='Natural':assert load(r/('initial-gate-'+wall+'.json'))['status']=='pass'
  if a.time>.0016:
   previous=.0016 if a.time==.016 else .016 if a.time==.1 else .1
   assert load(r/'comparisons'/(a.model+'-'+a.policy+'-'+wall+'-'+str(previous)+'.json'))['status']=='pass','preceding evolved gate not passed'
  if a.time>=.1:
   admission=load(r/'duration-admission.json');assert [a.model,a.policy,wall] in admission['configurations'] and admission['status']=='pass'
  for label,mesh,dt in zip(labels(a.model,a.policy,wall),('Balanced5V1','Balanced11V1','Balanced5V1','Balanced11V1'),(.0004,.0004,.0002,.0002),strict=True):
   m=r/'meshes'/((('E' if a.model=='E94E95' else a.model)+'-'+('N' if a.policy=='Natural' else 'P')+'-'+wall))/(mesh+'.json')
   cmd=[sys.executable,'-B',str(source/'run_rans_neutral_pilot.py'),'--root',str(r),'--label',label,'--model',a.model,'--policy',a.policy,'--wall',wall,'--degree','5','--mesh-file',str(m),'--dt',str(dt),'--max-time',str(a.time),'--milestones','.0016','.016','.1','.5','--timeout','14400']
   if (r/'campaigns'/label).exists():cmd+=['--resume']
   queue.append((label,cmd))
 write_json(out/'immutable-request.json',dict(model=a.model,policy=a.policy,walls=a.walls,target=a.time,commands=[c for _,c in queue],worker_limit=72,reserve=10,concurrency=a.concurrency))
 while queue or running:
  for label,(process,stream,cmd) in list(running.items()):
   code=process.poll()
   if code is None:continue
   stream.close();completed.append(dict(label=label,exit_code=code,command=cmd));del running[label]
  if queue and len(running)<a.concurrency and workers()+10<=72:
   label,cmd=queue.pop(0);stream=(out/(label+'.log')).open('x');process=subprocess.Popen(cmd,env=env,stdout=stream,stderr=subprocess.STDOUT);running[label]=(process,stream,cmd)
  write_json(out/'progress.json',dict(status='running',workers=workers(),running=list(running),queued=[n for n,c in queue],completed=completed))
  if queue or running:time.sleep(5)
 comparisons=[];(r/'comparisons').mkdir(exist_ok=True)
 for wall in a.walls:
  dest=r/'comparisons'/(a.model+'-'+a.policy+'-'+wall+'-'+str(a.time)+'.json');cmd=[sys.executable,'-B',str(source/'compare_rans_balanced_evolution.py'),'--root',str(r),'--labels',*labels(a.model,a.policy,wall),'--time',str(a.time),'--output',str(dest)]
  with (out/(wall+'-comparison.log')).open('x') as f:code=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
  comparisons.append(dict(wall=wall,exit_code=code,receipt=str(dest),status=load(dest)['status']))
 status='pass' if all(x['exit_code']==0 for x in completed+comparisons) else 'failed_or_unresolved';write_json(out/'results.json',dict(status=status,completed=completed,comparisons=comparisons,worker_limit=72));print(status);return 0 if status=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
