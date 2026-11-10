#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Declared evolved-space and temporal-uncertainty gates on common samples."""
import argparse,json,sys
from pathlib import Path
from rans_qualification import compare_samples,write_json
from rans_preservation_oracle import Reference
from run_rans_neutral_pilot import scaled,FIELDS
from qualify_rans_preservation import load
from rans_request_manifest import RequestManifest

def compare(root,coarse,fine,control,T):
 folders=[root/'campaigns'/x for x in (coarse,fine,control)]
 requests=[load(f/'request-manifest.json')['request'] for f in folders]
 for key in ('model','policy','wall','descriptor','boundary_policy','runtime','velocity_amplitude'):
  if any(r[key]!=requests[0][key] for r in requests[1:]):raise ValueError('comparison identity differs: '+key)
 if requests[0]['dt']!=requests[1]['dt'] or requests[2]['dt']!=requests[0]['dt']/2:raise ValueError('comparison timesteps differ')
 if any(requests[0][k]!=requests[2][k] for k in ('knots','degree','pressure_degree','quadrature','raw_neutral_sha256','boundary')):raise ValueError('same-space temporal initial/boundary identity differs')
 metrics=[];samples=[]
 for f,r in zip(folders,requests,strict=True):
  m=RequestManifest(f,r,True);p=f/('T-'+format(T,'.4f'))
  if not m.completed(p):raise ValueError('incomplete comparison endpoint '+str(p))
  metrics.append(load(p/'metrics.json'));samples.append(load(p/'endpoint-samples.json'))
 ref=Reference(requests[0]['model'],requests[0]['policy'],requests[0]['wall']);w=samples[0]['grid_weights']
 spatial=scaled(compare_samples(samples[0],samples[1],w),ref);temporal=scaled(compare_samples(samples[0],samples[2],w),ref);fields={}
 for f in FIELDS:
  ec,ef=[m['reference_error'][f]['weighted_RMS'] for m in metrics[:2]];delta=spatial[f]['weighted_RMS'];td=temporal[f]['weighted_RMS'];tol=max(1e-7,.1*delta)
  reduction=ef/ec if ec>1e-10 else None
  fields[f]=dict(coarse_reference_RMS=ec,fine_reference_RMS=ef,fine_over_coarse=reduction,spatial_difference_RMS=delta,temporal_difference_RMS=td,temporal_limit=tol,
   spatial_status='unresolved_floor' if reduction is None else 'pass' if reduction<=.85 else 'fail',temporal_uncertainty_status='pass' if td<=tol else 'unresolved',triangle_status='pass' if delta<=ec+ef+1e-12 else 'fail')
 status='pass' if all(f['spatial_status'] in ('pass','unresolved_floor') and f['temporal_uncertainty_status']=='pass' and f['triangle_status']=='pass' for f in fields.values()) else 'unresolved'
 return dict(status=status,time=T,theta=T*ref.us/ref.H,coarse=coarse,fine=fine,temporal_control=control,fields=fields,scope='evolved spatial reduction and dt/2 uncertainty; no inferred temporal order')

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--coarse',required=True);p.add_argument('--fine',required=True);p.add_argument('--control',required=True);p.add_argument('--time',type=float,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 try:r=compare(a.root,a.coarse,a.fine,a.control,a.time)
 except Exception as e:r=dict(status='fail',failure=repr(e))
 if a.output.exists():raise FileExistsError('comparison receipt already exists')
 write_json(a.output,r);print(json.dumps(r));return 0 if r['status']=='pass' else 2 if r['status']=='unresolved' else 1
if __name__=='__main__':sys.exit(main())
