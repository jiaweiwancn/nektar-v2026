#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only CT9 .001/.0005 reaggregation at the declared endpoint tolerance."""
import argparse,copy,json,os
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from ct_accepted9_pressure_run import sealed
from ct_accepted9_gate import pair
from ct_accepted8_native import norms,observations,read
from ct_accepted_feasibility import require,write_json
from ct_accepted5_traction import difference
from rans_request_manifest import RequestManifest,sha

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('ct8','finer','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();r=sealed(a.finer);pre=sealed(a.ct8/'preflight-v6');c=pre['case']
    m=RequestManifest(a.output,dict(schema='CT9PressureComparisonV2',finer=sha(a.finer/'results.json'),ct8=sha(a.ct8/'strip-T008-v1/results.json'),source=sha(Path(__file__)),endpoint_tolerance=1e-13,correction='read-only replacement of an exact floating-time equality; no state, history, threshold or old receipt edited'))
    ini=read(Path(pre['initial']['0.001']['directory'])/'initial.json');n=pre['global_count'];weights=np.asarray(pre['weights']);regions=np.asarray(ini['observation']['grid_regions'][:n]);scale=c['ustar']**2
    reports={k:read(Path(v['directory'])/'case.evolution.json') for k,v in r['members'].items()}
    for k,v in r['members'].items():require(sha(Path(v['directory'])/'case.evolution.json')==v['native_report_sha256'] and abs(reports[k]['time']-.008)<1e-13,'changed/incomplete actual report')
    left,right=[reports[k] for k in ('0.001','0.0005')]
    obs=[observations(ini,q['state'],n) for q in (left,right)];delta=norms(*obs,weights,regions,c);r=copy.deepcopy(r);r['temporal']={};r['traction_temporal']={};r['common_time_pressure']=[]
    for f,v in delta.items():
        drift=max(x['accepted_metrics'][-1]['drift'][f]['scaled'] for x in r['members'].values());r['temporal'][f]=dict(difference=v,measured_drift=drift,limit=max(1e-7,.1*drift))
    rows=list(r['members'].values());spread=difference(*(x['traction_history'][-1]['profile'] for x in rows),c);drifts=[difference(x['actual_initial_traction_profile'],x['traction_history'][-1]['profile'],c) for x in rows]
    for f in ('ground','interface','signed_mean'):
        d=abs(spread['signed_mean_difference']) if f=='signed_mean' else spread['profiles'][f]['normalized_weighted_RMS'];v=max(abs(x['signed_mean_difference']) if f=='signed_mean' else x['profiles'][f]['normalized_weighted_RMS'] for x in drifts);r['traction_temporal'][f]=dict(spread=d,measured_drift=v,limit=max(1e-7,.1*v))
    histories=[{round(t,12):s for t,s in zip(q['audit_times'],q['accepted_states'],strict=True)} for q in (left,right)]
    for t in sorted(histories[0].keys()&histories[1].keys()):
        x,y=[observations(ini,q[t],n)[:,4] for q in histories];d=x-y;entry={}
        for label,mask in [('whole',np.ones(n,bool)),('inner',regions==0),('outer',regions==1)]:
            w=weights[mask];v=d[mask];mean=float(w@v/w.sum());entry[label]=dict(signed_mean=mean,weighted_RMS=float(np.sqrt(w@(v*v)/w.sum())),spatial_variation_RMS=float(np.sqrt(w@((v-mean)**2)/w.sum())),scale=scale)
        r['common_time_pressure'].append(dict(time=t,states=[exact_vector_hash(q[t]) for q in histories],pressure=entry))
    r['comparison_version']='CT9PressureComparisonV2';r['qualification']=pair(r);r['retained_incomplete_comparison']=dict(path=str(a.finer/'results.json'),sha256=sha(a.finer/'results.json'),reason='exact float endpoint check omitted temporal fields; physical member complete; preserved unchanged')
    write_json(a.output/'results.json',r);m.seal(a.output,['results.json']);print('CT9 actual .001/.0005 comparison',r['qualification'],flush=True);return 0 if r['qualification']['status']=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
