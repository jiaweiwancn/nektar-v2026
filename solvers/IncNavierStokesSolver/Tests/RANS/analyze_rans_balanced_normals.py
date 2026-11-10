#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Nonlinear-model normal-stress and pressure receipts on unchanged samples."""
import argparse,os
from pathlib import Path
os.environ.update(dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'),'1'))
import numpy as np
from qualify_rans_preservation import load
from validate_a93_a95_e94e95 import oracle
from rans_preservation_oracle import Reference
from rans_qualification import write_json

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();rows=[]
 for model in ('A95','A93'):
  for wall in ('EP','Accepted'):
   ref=Reference(model,'Natural',wall)
   for space in ('B5','B11'):
    folder=a.root/'campaigns'/(model+'-N-'+wall+'-'+space)
    for T in (.0016,.016):
     endpoint=folder/('T-'+format(T,'.4f'));m=load(endpoint/'metrics.json');q=load(endpoint/'endpoint-samples.json');v=np.array(q['grid_values']);g=np.array(q['grid_gradients']);t=np.array(q['grid_traction']);w=np.array(q['grid_weights']);outer=np.array(q['grid_regions'])==1;assert len(v)==len(g)==len(t)==len(w)==len(outer) and abs(q['time']-T)<1e-13
     Rzz=2*v[:,2]/3-v[:,4]-t[:,2];diags=[];exact=[];error=0.;div=[]
     for j in np.flatnonzero(outer):
      G=np.zeros((3,3));G[0]=g[j,0];G[2]=g[j,1];c=oracle(model,v[j,2],v[j,3],G);R=c['R'];diags.append(np.diag(R));exact.append(np.diag(ref.at(q['grid'][j][2],1)['R']));error=max(error,abs(R[2,2]-Rzz[j]));div.append(np.trace(G))
     assert error<1e-12,error
     delta=np.array(diags)-np.array(exact);normal=np.sqrt(np.sum(w[outer,None]*delta*delta,axis=0)/w[outer].sum())/ref.us**2
     rows.append(dict(model=model,wall=wall,space=space,time=q['time'],status='pass',outer_Rxx_Ryy_Rzz_reference_RMS_scaled=normal.tolist(),independent_outer_Rzz_vs_native_full_traction_maximum=error,raw_outer_divergence_maximum=float(np.max(np.abs(div))),pressure_reference_RMS=m['reference_error']['Pi'],vertical_velocity_reference_RMS=m['reference_error']['v'],native_spatial_residual=m['spatial_alpha_zero'],native_extrema=m['extrema'],source=str(endpoint)))
 assert not a.output.exists();write_json(a.output,dict(status='pass',checks=rows,scope='independent full normal-stress evaluations; native full normal traction validates Rzz; full 3D kernels/covariance evidence reused within recorded source identity'));return 0
if __name__=='__main__':raise SystemExit(main())
