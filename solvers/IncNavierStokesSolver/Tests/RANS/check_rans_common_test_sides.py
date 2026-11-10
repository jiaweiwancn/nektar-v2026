#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Native explicit-parent test-basis selectors, including intended rejections."""
import argparse,gzip,json,os,shutil,subprocess
from pathlib import Path
os.environ.update(dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'),'1'))
from rans_qualification import THREADS,write_json
from qualify_rans_preservation import replay_request,load

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args();r=a.root.resolve();runtime=r/'runtime-resolution-v2';out=r/'common-test-selector-regression';out.mkdir(exist_ok=False);os.environ.update(THREADS,LD_LIBRARY_PATH=str(runtime/'lib'),TMPDIR=str(r/'tmp'));common=r/'diagnosis/initial-Accepted/common'
 with gzip.open(common/'packet.json.gz','rt') as f:seed=json.load(f)
 elements=seed['observation']['elements'];inner=next(e for e in elements if e['lower'][2]==0);outer=next(e for e in elements if e['lower'][2]==2);q=[0.,0.,(inner['lower'][2]+inner['upper'][2])/2];base=replay_request(seed,[q]);base.update(export_test_basis=True,grid_test_geometries=[inner['geometry']]);checks=[]
 controls=[('valid',base,None),('unknown',dict(base,grid_test_geometries=[999999]),'unknown explicit native test geometry'),('count',dict(base,grid_test_geometries=[inner['geometry'],inner['geometry']]),'test geometry selector count differs'),('side-conflict',dict(base,grid=[[0.,0.,2.]],grid_sides=[-1],grid_test_geometries=[outer['geometry']]),'test geometry conflicts with declared interface side')]
 for name,request,reason in controls:
  folder=out/name;folder.mkdir();shutil.copy2(common/'case.xml',folder/'case.xml');write_json(folder/'request.json',request)
  with (folder/'probe.log').open('w') as f:c=subprocess.run([str(runtime/'RANSPreservationProbe'),'case.xml','request.json','packet.json'],cwd=folder,stdout=f,stderr=subprocess.STDOUT,timeout=120)
  d=load(folder/'packet.json')
  if reason:assert c.returncode!=0 and reason in d.get('failure','')
  else:assert c.returncode==0 and d['test_basis'][0]['velocity_k']['geometry']==inner['geometry']
  checks.append(dict(name=name,status='pass',exit_code=c.returncode,expected_reason=reason,observed_reason=d.get('failure')))
 write_json(out/'results.json',dict(status='pass',checks=checks));return 0
if __name__=='__main__':raise SystemExit(main())
