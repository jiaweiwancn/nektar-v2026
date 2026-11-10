#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Apply the screened E94E95 strategy to A95/A93, preserving independent failures."""
import argparse,os,subprocess,sys
from pathlib import Path
from rans_qualification import THREADS,write_json
from rans_preservation_oracle import Reference
from make_rans_balanced_meshes import specifications
from qualify_rans_preservation import load

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--models',nargs='+',default=['A95','A93']);a=p.parse_args();r=a.root.resolve();out=r/'model-qualification';out.mkdir(exist_ok=False);source=Path(__file__).parent;env=dict(os.environ,**THREADS);checks=[]
 def call(name,cmd):
  with (out/(name+'.log')).open('x') as f:code=subprocess.run([sys.executable,'-B',*map(str,cmd)],env=env,stdout=f,stderr=subprocess.STDOUT).returncode
  checks.append(dict(name=name,exit_code=code,status='pass' if code==0 else 'fail',command=list(map(str,cmd))));write_json(out/'results.json',dict(status='running',checks=checks));return code
 for model in a.models:
  walls=[]
  for wall in ('EP','Accepted'):
   label=model+'-Natural-'+wall;folder=r/'meshes'/(model+'-N-'+wall);folder.mkdir(parents=True,exist_ok=False);ref=Reference(model,'Natural',wall);spaces=specifications(ref)
   for n,s in spaces.items():write_json(folder/(n+'.json'),s)
   spaces['BaseP5']=dict(degree=5,knots=ref.knots());write_json(folder/'screen-spaces.json',spaces)
   if call(label+'-initial',[source/'diagnose_rans_balanced_spatial.py','--root',r,'--analysis-runtime',r/'runtime-resolution-v2','--model',model,'--wall',wall,'--label',label,'--specifications',folder/'screen-spaces.json']):continue
   gate=r/('initial-gate-'+label+'.json')
   if call(label+'-gate',[source/'classify_rans_balanced_spatial.py','--root',r,'--model',model,'--wall',wall,'--screen',label,'--output',gate])==0:walls.append(wall)
  if not walls:continue
  code=call(model+'-evolved-0.0016',[source/'schedule_rans_balanced_evolution.py','--root',r,'--model',model,'--walls',*walls,'--time','.0016','--concurrency','4'])
  # Continue only each independently passing model/wall comparison.
  passed=[w for w in walls if (r/'comparisons'/(model+'-Natural-'+w+'-0.0016.json')).exists() and load(r/'comparisons'/(model+'-Natural-'+w+'-0.0016.json'))['status']=='pass']
  if passed:call(model+'-evolved-0.016',[source/'schedule_rans_balanced_evolution.py','--root',r,'--model',model,'--walls',*passed,'--time','.016','--concurrency','4'])
 status='pass' if len(checks)==12 and all(x['status']=='pass' for x in checks) else 'failed_or_incomplete';write_json(out/'results.json',dict(status=status,checks=checks));return 0 if status=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
