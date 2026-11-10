#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only final-horizon audits after the bounded four-member comparison."""
import argparse,os,subprocess,sys,time
from pathlib import Path
from rans_qualification import THREADS,write_json
from qualify_rans_preservation import load

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--wait',action='store_true');a=p.parse_args();r=a.root.resolve();source=Path(__file__).parent;target=r/'comparisons/E-N-Accepted-duration-0.002-0.5.json';env=dict(os.environ,**THREADS)
 deadline=time.monotonic()+14400
 while a.wait and not target.exists():
  if time.monotonic()>deadline:raise TimeoutError('bounded duration comparison wait expired')
  time.sleep(5)
 if not target.exists():raise FileNotFoundError('final comparison missing')
 out=r/'duration-final-audits';out.mkdir(exist_ok=False);comparison=load(target);jobs=[]
 labels=['E-N-Accepted-duration-0.002-'+s for s in ('B5','B11','B5-dt2','B11-dt2')]
 for index,label in enumerate(labels):
  commands=[]
  endpoint=r/'campaigns'/label/'T-0.5000/metrics.json'
  if not endpoint.exists() or load(endpoint)['status']!='pass':continue
  if index<2:commands.append([sys.executable,'-B',str(source/'audit_rans_balanced_evolved.py'),'--root',str(r),'--analysis-runtime',str(r/'runtime-resolution-v2'),'--campaign-root',str(r),'--label',label,'--output-label','duration-at-0.5-'+('B5' if index==0 else 'B11'),'--times','.5'])
  commands.append([sys.executable,'-B',str(source/'audit_rans_pilot_samples.py'),'--root',str(r),'--label',label]);jobs.append((label,commands))
 running={};checks=[]
 while jobs or running:
  for label,(process,f,commands,position) in list(running.items()):
   code=process.poll()
   if code is None:continue
   f.close();checks.append(dict(label=label,command=commands[position],exit_code=code))
   if position+1<len(commands):
    position+=1;f=(out/(label+'-'+str(position)+'.log')).open('x');process=subprocess.Popen(commands[position],stdout=f,stderr=subprocess.STDOUT,env=env);running[label]=(process,f,commands,position)
   else:del running[label]
  if jobs and len(running)<2:
   label,commands=jobs.pop(0);f=(out/(label+'-0.log')).open('x');process=subprocess.Popen(commands[0],stdout=f,stderr=subprocess.STDOUT,env=env);running[label]=(process,f,commands,0)
  write_json(out/'progress.json',dict(status='running',checks=checks,running=list(running),queued=[x[0] for x in jobs]))
  if jobs or running:time.sleep(5)
 if comparison['status']=='pass':
  cmd=[sys.executable,'-B',str(source/'check_rans_balanced_state_quadrature.py'),'--root',str(r),'--time','.5','--walls','Accepted','--label-prefix','E-N-Accepted-duration-0.002']
  with (out/'final-quadrature.log').open('x') as f:code=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=env).returncode
  checks.append(dict(label='final-quadrature',command=cmd,exit_code=code))
 status='pass' if comparison['status']=='pass' and len(checks)==7 and all(x['exit_code']==0 for x in checks) else 'failed_or_incomplete';write_json(out/'results.json',dict(status=status,checks=checks,comparison_status=comparison['status']));return 0 if status=='pass' else 1
if __name__=='__main__':
 if '--admission' in sys.argv:
  from finish_rans_configuration_duration import main
 raise SystemExit(main())
