#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Two selected prescribed-kappa controls; no full prescribed campaign claim."""
import argparse,os,subprocess,sys
from pathlib import Path
from make_rans_balanced_meshes import specifications
from rans_preservation_oracle import Reference
from rans_qualification import THREADS,write_json
from qualify_rans_preservation import load

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args();r=a.root.resolve();out=r/'prescribed-controls';out.mkdir(exist_ok=False);source=Path(__file__).parent;env=dict(os.environ,**THREADS);rows=[]
 for model,wall in [('A95','Accepted'),('A93','EP')]:
  prefix=model+'-P-'+wall;folder=r/'meshes'/prefix;folder.mkdir(parents=True,exist_ok=False);ref=Reference(model,'Prescribed',wall);spec=specifications(ref)['Balanced11V1'];write_json(folder/'Balanced11V1.json',spec);label=prefix+'-B11-short'
  cmd=[sys.executable,'-B',str(source/'run_rans_neutral_pilot.py'),'--root',str(r),'--label',label,'--model',model,'--policy','Prescribed','--wall',wall,'--mesh-file',str(folder/'Balanced11V1.json'),'--dt','.0004','--max-time','.0016']
  with (out/(prefix+'.log')).open('w') as f:code=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
  receipt=r/'campaigns'/label/'T-0.0016/metrics.json';metric=load(receipt) if receipt.exists() else None
  rows.append(dict(model=model,wall=wall,policy='Prescribed',kappa=.41,label=label,exit_code=code,status='pass' if code==0 and metric and metric['status']=='pass' else 'fail',source=str(receipt)))
 status='pass' if len(rows)==2 and all(x['status']=='pass' for x in rows) else 'fail';write_json(out/'results.json',dict(status=status,checks=rows,scope='two prescribed .41 fine-space endpoint controls through .0016; no new prescribed temporal/spatial order'));return 0 if status=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
