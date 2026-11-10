#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Frozen evolved coefficients: actual Q48 versus Q64 operators on both spaces."""
import argparse,os,shutil,xml.etree.ElementTree as ET
from pathlib import Path
os.environ.update(dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'),'1'))
from rans_qualification import THREADS,write_json,exact_vector_hash
from qualify_rans_preservation import load,replay_request,Q
from diagnose_rans_balanced_spatial import run_probe
from run_cls_epawf_picard import set_parameter,write

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--time',type=float,default=.016);p.add_argument('--walls',choices=('EP','Accepted'),nargs='+',default=['EP','Accepted']);p.add_argument('--label-prefix');a=p.parse_args();r=a.root.resolve();runtime=r/'runtime-resolution-v2';out=r/('evolved-quadrature-'+str(a.time));out.mkdir(exist_ok=False);os.environ.update(THREADS,LD_LIBRARY_PATH=str(runtime/'lib'),TMPDIR=str(r/'tmp'));checks=[]
 for wall in a.walls:
  folders=[r/'campaigns'/((a.label_prefix or ('E-N-'+wall))+'-'+suffix) for suffix in ('B5','B11')];base=[load(f/('T-'+format(a.time,'.4f'))/'metrics.json')['spatial_alpha_zero_full'] for f in folders]
  for index,folder in enumerate(folders):
   endpoint=folder/('T-'+format(a.time,'.4f'));native=load(endpoint/'case.evolution.json');seed=load(folder/'initial-packet.json');target=out/(wall+'-'+str(index));target.mkdir();root=ET.parse(endpoint/'case.xml').getroot()
   for key,value in zip(('RANSWallNormalOrder','RANSWallTangentialOrder','RANSWallMassOrder','RANSWallVolumeOrder','RANSWallInterfaceOrder'),Q[2],strict=True):set_parameter(root,key,value)
   write(root,target/'case.xml');request=replay_request(seed,[],state=native['state'],time=native['time']);request['fixed_boundary']=native['initialization_receipt']['live_fixed_boundary'];packet,_=run_probe(runtime,target,request);q=packet['observation']['audit_alpha_zero'];receipt=packet['initialization_receipt'];assert exact_vector_hash(native['state'])==receipt['raw_initial_sha256']==receipt['effective_initial_sha256'];assert receipt['explicit_boundary_sha256']==native['initialization_receipt']['explicit_boundary_sha256']
   rows=[]
   for key in ('coupled_norms','native_k_inner_outer_shared_interface_dual'):
    for j,(x,y,z) in enumerate(zip(base[index][key],q[key],base[1-index][key],strict=True)):
     change=abs(x-y);difference=abs(x-z);limit=max(1e-10,.01*difference);rows.append(dict(metric=key,index=j,quadrature_change=change,spatial_difference=difference,limit=limit,status='pass' if change<=limit else 'fail'))
   checks.append(dict(wall=wall,space='Balanced5V1' if index==0 else 'Balanced11V1',time=native['time'],status='pass' if all(x['status']=='pass' for x in rows) else 'fail',same_raw_effective_reference_boundary=True,state_sha256=receipt['raw_initial_sha256'],source_orders=Q[1],target_orders=Q[2],actual_target_quadrature=packet['observation']['actual_quadrature'],checks=rows))
 expected={(wall,space) for wall in a.walls for space in ('Balanced5V1','Balanced11V1')};observed={(x['wall'],x['space']) for x in checks}
 status='pass' if len(checks)==len(expected) and observed==expected and all(x['status']=='pass' for x in checks) else 'fail';write_json(out/'results.json',dict(status=status,checks=checks,required_coverage=sorted(expected),observed_coverage=sorted(observed),scope='frozen actual evolved-state spatial residuals; no initial regeneration or time advance'));return 0 if status=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
