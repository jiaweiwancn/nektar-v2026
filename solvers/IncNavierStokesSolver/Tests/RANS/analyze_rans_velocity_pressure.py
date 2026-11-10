#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded divergence-free response, compared to unperturbed numerical states."""
import argparse,os,sys
from pathlib import Path
from rans_qualification import THREADS,compare_samples,write_json,audited_extrema
os.environ.update(THREADS)
import numpy as np
from rans_request_manifest import RequestManifest
from rans_preservation_oracle import Reference
from run_rans_neutral_pilot import scaled
from qualify_rans_preservation import load
from validate_a93_a95_e94e95 import oracle

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=a.root;labels=['A95-N-EP-h0p5','A95-N-EP-velocity'];folders=[r/'campaigns'/s for s in labels];requests=[load(f/'request-manifest.json')['request'] for f in folders]
 for key in ('descriptor','knots','degree','pressure_degree','quadrature','dt','boundary','runtime','raw_neutral_sha256'):
  assert requests[0][key]==requests[1][key],key
 actual=requests[1]['actual_initial']['receipt'];assert actual['raw_initial_sha256']==actual['effective_initial_sha256'] and actual['explicit_boundary_sha256']==actual['implicit_boundary_sha256'] and actual['raw_effective_max_difference']==0.
 ref=Reference('A95','Natural','EP');initial=load(folders[1]/'initial-metrics.json');div=initial['residual']['strong_continuity'][:2];assert max(div)<1e-12
 result=dict(status='running',amplitude=requests[1]['velocity_amplitude'],initial_divergence_L2=div,actual_initial_identity={k:v for k,v in actual.items() if k not in ('live_fixed_boundary','effective_initial')},initial_native_history_sha256=requests[1]['actual_initial']['history_sha256'],observations=[])
 for T in (.0016,.016):
  packets=[];reports=[]
  for f,request in zip(folders,requests,strict=True):
   d=f/('T-'+format(T,'.4f'));assert RequestManifest(f,request,True).completed(d)
   m=load(d/'metrics.json');assert m['status']=='pass';packets.append(load(d/'endpoint-samples.json'));reports.append(load(d/'case.evolution.json'))
  fields=scaled(compare_samples(*packets,packets[0]['grid_weights']),ref);w=np.array(packets[0]['grid_weights']);regions=np.array(packets[0]['grid_regions']);outer=regions==1
  normals=[];traction=[];maximum_oracle_error=0.
  for packet in packets:
   values=np.array(packet['grid_values']);g=np.array(packet['grid_gradients']);t=np.array(packet['grid_traction']);Rzz=2*values[:,2]/3-values[:,4]-t[:,2];Rzz[~outer]+=2*ref.nu*g[~outer,1,2]
   diag=[]
   for j in np.flatnonzero(outer):
    G=np.zeros((3,3));G[0]=g[j,0];G[2]=g[j,1];R=oracle('A95',values[j,2],values[j,3],G)['R'];diag.append(np.diag(R));maximum_oracle_error=max(maximum_oracle_error,abs(R[2,2]-Rzz[j]))
   normals.append(np.array(diag));traction.append(Rzz)
  assert maximum_oracle_error<1e-12,maximum_oracle_error
  delta=normals[1]-normals[0];stress=np.sqrt(np.sum(w[outer,None]*delta**2,axis=0)/np.sum(w[outer]))/ref.us**2
  item=dict(time=T,theta=T*ref.us/ref.H,perturbed_minus_unperturbed_weighted_RMS=fields,outer_normal_stress_response_Rxx_Ryy_Rzz_scaled=stress.tolist(),whole_patch_Rzz_response_scaled=float(np.sqrt(w@((traction[1]-traction[0])**2)/w.sum())/ref.us**2),independent_stress_vs_native_traction_max_error=maximum_oracle_error,physical_extrema=audited_extrema(reports[1]))
  item['physical_extrema'].pop('final_state',None);result['observations'].append(item)
 assert all(x['perturbed_minus_unperturbed_weighted_RMS']['Pi']['weighted_RMS']>1e-10 for x in result['observations'])
 result['status']='pass';result['scope']='bounded response, no monotone decay or global stability assertion';assert not a.output.exists();write_json(a.output,result);print(result['status']);return 0
if __name__=='__main__':sys.exit(main())
