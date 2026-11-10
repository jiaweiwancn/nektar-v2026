#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Same-coefficient p5 frozen quadrature check with identical effective boundaries."""
import argparse,os,shutil,sys,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS,write_json,exact_vector_hash
from qualify_rans_preservation import load,probe,replay_request,Q
from run_cls_epawf_picard import set_parameter,write

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args();r=a.root.resolve();out=r/'selected-quadrature';out.mkdir(exist_ok=False);os.environ.update(THREADS,LD_LIBRARY_PATH=str(r/'runtime/lib'),TMPDIR=str(r/'tmp'));result=dict(status='running',checks=[])
 for wall in ('EP','Accepted'):
  try:
   coarse=r/'campaigns'/('E-N-'+wall+'-h0p5');fine=r/'campaigns'/('E-N-'+wall+'-Interface8p5');seed=load(coarse/'initial-packet.json');c=load(coarse/'initial-metrics.json')['residual'];f=load(fine/'initial-metrics.json')['residual'];folder=out/wall;folder.mkdir();root=ET.parse(coarse/'T-0.0016/case.xml').getroot()
   for key,value in zip(('RANSWallNormalOrder','RANSWallTangentialOrder','RANSWallMassOrder','RANSWallVolumeOrder','RANSWallInterfaceOrder'),Q[2],strict=True):set_parameter(root,key,value)
   write(root,folder/'case.xml');full=probe(r/'runtime/RANSPreservationProbe',folder,replay_request(seed,[],state=seed['reference']));q=full['observation']['audit_alpha_zero'];receipt=full['initialization_receipt'];assert exact_vector_hash(seed['reference'])==receipt['raw_initial_sha256']==receipt['effective_initial_sha256'];assert receipt['explicit_boundary_sha256']==seed['initialization_receipt']['explicit_boundary_sha256']
   gates=[]
   for key in ('coupled_norms','native_k_inner_outer_shared_interface_dual'):
    for j,(v,z,s) in enumerate(zip(c[key],q[key],f[key],strict=True)):
     delta=abs(z-v);bound=max(1e-10,.01*abs(v-s));gates.append(dict(metric=key,index=j,quadrature_change=delta,spatial_difference=abs(v-s),limit=bound,status='pass' if delta<=bound else 'unresolved'))
   result['checks'].append(dict(wall=wall,status='pass' if all(x['status']=='pass' for x in gates) else 'unresolved',same_coefficient_and_effective_boundary_inputs=True,source_quadrature=Q[1],target_quadrature=Q[2],actual_target_quadrature=full['observation']['actual_quadrature'],gates=gates,initialization_receipt={k:v for k,v in receipt.items() if k not in ('live_fixed_boundary','effective_initial')}))
  except Exception as e:result['checks'].append(dict(wall=wall,status='fail',failure=repr(e)))
  write_json(out/'results.json',result)
 result['status']='pass' if len(result['checks'])==2 and all(x['status']=='pass' for x in result['checks']) else 'fail';write_json(out/'results.json',result);print(result['status']);return 0 if result['status']=='pass' else 1
if __name__=='__main__':sys.exit(main())
