#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded natural-kappa DAG; at most four pilot controllers, one extra slot.

Existing controllers/requests are observed, never interrupted. Every numerical
advance uses the immutable-request runner and native accuracy/domain guards.
"""
import argparse,json,os,subprocess,sys,time
from pathlib import Path
from rans_qualification import THREADS,write_json
from qualify_rans_preservation import load

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--concurrency',type=int,default=4);a=p.parse_args();r=a.root.resolve();out=r/'natural-schedule';out.mkdir(exist_ok=False);source=Path(__file__).parent;env=dict(os.environ,**THREADS);children={};done={};blocked={};events=[]
 def log(kind,**fields):events.append(dict(event=kind,time=time.time(),**fields));write_json(out/'events.json',events)
 def workers():
  total=0
  for line in subprocess.check_output(['ps','-eo','nlwp,comm,args'],text=True).splitlines()[1:]:
   parts=line.split(None,2)
   if len(parts)<3:continue
   n,comm,cmd=parts
   if comm.startswith(('IncNavier','RANS','Nek','cc1plus','f951','orted','mpiexec')) or (comm.startswith('python') and 'nektar' in cmd):total+=int(n)
  return total
 def running_labels():
  labels=set()
  for line in subprocess.check_output(['ps','-eo','args'],text=True).splitlines():
   if 'run_rans_neutral_pilot.py' not in line or '--label ' not in line or str(r) not in line:continue
   labels.add(line.split('--label ',1)[1].split()[0])
  return labels
 def endpoint(label,T):
  f=r/'campaigns'/label/('T-'+format(T,'.4f'))
  return (f/'completed-artifacts.json').exists() and (f/'metrics.json').exists() and load(f/'metrics.json')['status']=='pass'
 def failed(label):
  f=r/'campaigns'/label/'results.json'
  return f.exists() and bool(load(f).get('first_failure'))
 def start(label,model,wall,dt,T,refinement=False):
  key=label+'@'+str(T)
  if endpoint(label,T) or failed(label) or label in running_labels() or key in children or key in blocked:return False
  if len(running_labels())>=a.concurrency or workers()+10>72:return False
  command=[sys.executable,'-B',str(source/'run_rans_neutral_pilot.py'),'--root',str(r),'--label',label,'--model',model,'--wall',wall,'--dt',str(dt),'--max-time',str(T)]
  if refinement:command+=['--refinement','Interface8V1']
  if (r/'campaigns'/label).exists():command+=['--resume']
  f=(out/(key.replace('@','-')+'.log')).open('x');process=subprocess.Popen(command,stdout=f,stderr=subprocess.STDOUT,env=env);children[key]=(process,f,command);log('start',key=key,command=command,pid=process.pid);return True
 def gate(model,wall,prefix,T):
  labels=[prefix+'-h0p5',prefix+'-Interface8p5',prefix+'-h0p5-dt2'];key=prefix+'-space-'+str(T)
  if key in done:return done[key]
  if key in blocked:return False
  if not all(endpoint(label,T) for label in labels):return False
  target=out/(key+'.json');cmd=[sys.executable,'-B',str(source/'compare_rans_neutral_pilot.py'),'--root',str(r),'--coarse',labels[0],'--fine',labels[1],'--control',labels[2],'--time',str(T),'--output',str(target)]
  with (out/(key+'.log')).open('x') as f:code=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=env).returncode
  status=load(target)['status'];log('gate',key=key,exit_code=code,status=status)
  if code==0 and status=='pass':done[key]=True;return True
  blocked[key]=dict(status=status,receipt=str(target));return False
 while True:
  for key,(process,f,cmd) in list(children.items()):
   code=process.poll()
   if code is None:continue
   f.close();log('complete',key=key,exit_code=code)
   if code:blocked[key]=dict(exit_code=code,command=cmd)
   del children[key]
  for model in ('E94E95','A95','A93'):
   for wall in ('EP','Accepted'):
    prefix=('E' if model=='E94E95' else model)+'-N-'+wall
    for suffix,dt,refine in [('-h0p5',.0004,False),('-Interface8p5',.0004,True),('-h0p5-dt2',.0002,False)]:start(prefix+suffix,model,wall,dt,.0016,refine)
    if gate(model,wall,prefix,.0016):
     for suffix,dt,refine in [('-h0p5',.0004,False),('-Interface8p5',.0004,True),('-h0p5-dt2',.0002,False)]:start(prefix+suffix,model,wall,dt,.016,refine)
    if gate(model,wall,prefix,.016):
     label=prefix+'-h0p5'
     if endpoint(label,.05):start(label,model,wall,.0004,.1)
     else:start(label,model,wall,.0004,.05)
  required=[('E' if m=='E94E95' else m)+'-N-'+w+'-h0p5' for m in ('E94E95','A95','A93') for w in ('EP','Accepted')]
  status='pass' if all(endpoint(x,.1) for x in required) else 'running'
  progress=dict(status=status,accounted_worker_threads=workers(),worker_limit=72,children=list(children),external_running=sorted(running_labels()),gates=done,blocked=blocked,completed_pilots=[x for x in required if endpoint(x,.1)],failed_pilots=[x for x in required if failed(x)])
  write_json(out/'progress.json',progress)
  if status=='pass':break
  if not children and not running_labels():progress['status']='blocked_by_failed_or_unresolved_gates';write_json(out/'progress.json',progress);return 1
  time.sleep(5)
 return 0
if __name__=='__main__':sys.exit(main())
