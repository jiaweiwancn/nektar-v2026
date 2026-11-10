#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Production preservation packet: explicit one-sided interface provenance."""
import argparse,copy,os,subprocess,sys
from pathlib import Path
import numpy as np
from qualify_rans_preservation import probe,replay_request,load,Q
from rans_qualification import THREADS,write_json
from rans_preservation_oracle import Reference
from run_a93_a95_e94e95_matrix import case_xml
from validate_a93_a95_e94e95 import oracle

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path);a=p.parse_args();r=a.root.resolve();out=a.output or r/'interface-side-regression';out.mkdir(parents=True,exist_ok=False)
 os.environ.update(THREADS,LD_LIBRARY_PATH=str(r/'runtime/lib'),TMPDIR=str(r/'tmp'))
 result=dict(status='running',checks=[])
 try:
  ref=Reference('A95','Natural','Accepted');case_xml(out/'case','A95','Natural','Accepted',steps=0,descriptor_probe=r/'runtime/RANSNewModelsProbe',radial_knots=ref.knots(),quadrature=Q[0],boundary_policy='FrozenReferenceV1')
  folder=out/'case';exe=r/'runtime/RANSPreservationProbe';base=probe(exe,folder,dict(operation='neutral',grid=[],export_basis=True));o=base['observation'];off=o['offsets'];state=list(o['state'])
  faces=[f for f in o['faces'] if f['interface']];outer_ids=set(i for f in faces if f['region']==1 for i in f['velocity_k']['ids']);boundary=dict(base['fixed_boundary']);chosen=None
  for f in faces:
   if f['region']!=0:continue
   row=f['velocity_k']
   for j,i in enumerate(row['ids']):
    if i not in outer_ids and off[2]+i not in boundary and abs(row['phi'][j])<1e-12 and abs(row['derivative'][2][j])>1e-4:
     chosen=(f['coordinate'],i);break
   if chosen:break
  assert chosen,'no interior one-sided derivative mode'
  x,i=chosen;state[off[0]+i]+=1e-7;state[off[2]+i]+=1e-7
  request=replay_request(base,[x,x],state=state);request['grid_sides']=[-1,1];packet=probe(exe,folder,request)['observation']
  v=np.array(packet['grid_values']);g=np.array(packet['grid_gradients']);geom=packet['grid_geometry_by_field']
  assert np.max(np.abs(v[0,:3]-v[1,:3]))<1e-12,'common velocity/k values'
  assert abs(g[0,0,2]-g[1,0,2])>1e-10 and abs(g[0,2,2]-g[1,2,2])>1e-10,'one-sided derivatives collapsed'
  assert all(geom[0][j]!=geom[1][j] for j in (0,1,2,4)) and geom[0][3]==-1
  assert abs(v[0,4]-v[1,4])>1e-6 and packet['grid_regions']==[0,1]
  expected=[]
  for side in (0,1):
   G=np.zeros((3,3));G[0]=g[side,0];G[2]=g[side,1];k=v[side,2];pi=v[side,4]
   if side==0:stress=(ref.nu+ref.kernel.nt(2.))*(G+G.T)-pi*np.eye(3)
   else:stress=-oracle('A95',k,v[side,3],G)['R']+(2*k/3-pi)*np.eye(3)
   expected.append(stress[:,2].tolist())
  error=float(np.max(np.abs(np.array(expected)-packet['grid_traction'])));assert error<1e-12,error
  result['checks'].append(dict(name='one-sided-values-gradients-pressure-traction',status='pass',packet={k:packet[k] for k in ('grid','grid_sides','grid_regions','grid_geometry_by_field','grid_values','grid_gradients','grid_traction','grid_inner_prescribed_epsilon')},independent_traction_max_error=error))
  bad=replay_request(base,[x],state=state);write_json(folder/'ambiguous.json',bad)
  with (folder/'ambiguous.log').open('w') as f:c=subprocess.run([str(exe),'case.xml','ambiguous.json','ambiguous-result.json'],cwd=folder,stdout=f,stderr=subprocess.STDOUT)
  reason=load(folder/'ambiguous-result.json').get('failure','');assert c.returncode!=0 and 'ambiguous interface sample' in reason,reason
  result['checks'].append(dict(name='ambiguous-interface-rejected',status='pass',failure=reason,exit_code=c.returncode));result['status']='pass'
 except Exception as e:result['checks'].append(dict(name='interface-regression',status='fail',failure=repr(e)));result['status']='fail'
 write_json(out/'results.json',result);print(result['status']);return 0 if result['status']=='pass' else 1
if __name__=='__main__':sys.exit(main())
