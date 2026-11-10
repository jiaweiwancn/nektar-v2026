#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Matched four-trajectory space/dt uncertainty with actual consumed identities."""
import argparse,hashlib,json
from pathlib import Path
from qualify_rans_preservation import load
from rans_qualification import compare_samples,write_json
from rans_request_manifest import RequestManifest
from rans_preservation_oracle import Reference
from run_rans_neutral_pilot import scaled,FIELDS

def compare(root,labels,T):
    folders=[root/'campaigns'/label for label in labels];reqs=[load(f/'request-manifest.json')['request'] for f in folders];metrics=[];samples=[];provenance=[]
    for key in ('model','policy','wall','descriptor','boundary_policy','runtime','velocity_amplitude'):
        if any(q[key]!=reqs[0][key] for q in reqs[1:]):raise ValueError('physics/runtime identity differs '+key)
    if reqs[0]['dt']!=reqs[1]['dt'] or reqs[2]['dt']!=reqs[0]['dt']/2 or reqs[3]['dt']!=reqs[1]['dt']/2:raise ValueError('invalid dt pair')
    for a,b in ((0,2),(1,3)):
        for key in ('knots','degree','pressure_degree','columns','quadrature','raw_neutral_sha256','boundary','declared_mesh'):
            if reqs[a].get(key)!=reqs[b].get(key):raise ValueError('same-space initial request differs '+key)
        receipt=[load(folders[i]/('T-'+format(reqs[i]['horizons'][0],'.4f'))/'metrics.json')['initialization_receipt'] for i in (a,b)]
        for key in ('raw_initial_sha256','effective_initial_sha256','fixed_reference_sha256','explicit_boundary_sha256','implicit_boundary_sha256'):
            if receipt[0][key]!=receipt[1][key]:raise ValueError('actual consumed inputs differ '+key)
        if any(x['raw_effective_max_difference']!=0 for x in receipt):raise ValueError('initial effective lift differs')
        for i,x in zip((a,b),receipt,strict=True):
            expected=reqs[i].get('actual_initial',{}).get('state_sha256',reqs[i]['raw_neutral_sha256'])
            if x['raw_initial_sha256']!=expected or x['fixed_reference_sha256']!=reqs[i]['reference_sha256'] or x['explicit_boundary_sha256']!=x['implicit_boundary_sha256']:raise ValueError('actual inputs do not match declared initialization')
    for folder,request in zip(folders,reqs,strict=True):
        admission=RequestManifest(folder,request,True);endpoint=folder/('T-'+format(T,'.4f'))
        if not admission.completed(endpoint):raise ValueError('incomplete endpoint '+str(endpoint))
        m=load(endpoint/'metrics.json');s=load(endpoint/'endpoint-samples.json');assert m['status']=='pass' and abs(m['time']-T)<1e-13
        metrics.append(m);samples.append(s);provenance.append(dict(label=folder.name,request_sha256=admission.payload['sha256'],space={k:request[k] for k in ('knots','degree','pressure_degree','quadrature','columns')},dt=request['dt'],time=m['time'],history=m['history_identity'],state_sha256=m['final_state_sha256'],metrics_sha256=hashlib.sha256((endpoint/'metrics.json').read_bytes()).hexdigest(),source=str(endpoint)))
    ref=Reference(reqs[0]['model'],reqs[0]['policy'],reqs[0]['wall']);weights=samples[0]['grid_weights'];dh=scaled(compare_samples(samples[0],samples[1],weights),ref);dc=scaled(compare_samples(samples[0],samples[2],weights),ref);df=scaled(compare_samples(samples[1],samples[3],weights),ref);fields={}
    for f in FIELDS:
        ec,ef=[m['reference_error'][f]['weighted_RMS'] for m in metrics[:2]];h=dh[f]['weighted_RMS'];limit=max(1e-7,.1*h);ratio=ef/ec if ec>1e-10 else None
        fields[f]=dict(coarse_reference_RMS=ec,fine_reference_RMS=ef,fine_over_coarse=ratio,spatial_difference=h,coarse_dt_difference=dc[f]['weighted_RMS'],fine_dt_difference=df[f]['weighted_RMS'],temporal_limit=limit,
            spatial_status='unresolved_floor' if ratio is None else 'pass' if ratio<=.85 else 'fail',temporal_status='pass' if max(dc[f]['weighted_RMS'],df[f]['weighted_RMS'])<=limit else 'unresolved',triangle_status='pass' if h<=ec+ef+1e-12 else 'fail',change_from_actual_initial=[m['change_from_initial'][f] for m in metrics])
    conservation=[]
    for m in metrics:
        defects=m['source_and_constraint_maxima'];checks={k:abs(v)<=(0 if k=='bound_calls' else 1e-10 if k=='corrected_continuity_gauge' else 1e-9) for k,v in defects.items() if k!='refreshed_trace_defect'}
        conservation.append(dict(status='pass' if all(checks.values()) else 'fail',checks=checks,defects=defects,refreshed_trace_is_unrepaired_splitting_error=True))
    status='pass' if all(v['spatial_status'] in ('pass','unresolved_floor') and v['temporal_status']=='pass' and v['triangle_status']=='pass' for v in fields.values()) and all(x['status']=='pass' for x in conservation) else 'unresolved_or_failed'
    return dict(status=status,time=T,theta=T*ref.us/ref.H,actual_consumed_temporal_inputs_identical_for_both_spaces=True,provenance=provenance,fields=fields,conservation=conservation,scope='field accuracy and temporal uncertainty, not a temporal-order measurement; alpha-zero spatial term audits separate')

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--labels',nargs=4,required=True);p.add_argument('--time',type=float,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();assert not a.output.exists()
 try:r=compare(a.root,a.labels,a.time)
 except Exception as e:r=dict(status='fail',failure=repr(e))
 write_json(a.output,r);print(r['status']);return 0 if r['status']=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
