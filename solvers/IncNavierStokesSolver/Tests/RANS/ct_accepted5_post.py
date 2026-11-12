#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Additional equal-time CT5 contrasts and CT4 scalar/outlet stage inventories."""
import argparse,json,os
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from ct_accepted4_reference import write_json
from ct_accepted4_run import expected,traction
from ct_accepted2_audit import values
from ct_accepted3_pressure import pressure_action
from ct_accepted2_storage import artifact
from ct_accepted_open import require
from ct_accepted_pilot import compare
from rans_request_manifest import RequestManifest,sha

def read(p):return json.loads(Path(p).read_text())

def scalar_inventory(packet,state,n,weights):
    o=values(packet,state,n)['observation'];v=np.array(o['grid_values']);regions=np.array(o['grid_regions']);w=np.array(weights);out={}
    for label,j,mask in [('inner_k',2,regions==0),('outer_k',2,regions==1),('outer_epsilon',3,regions==1)]:
        out[label]=dict(weighted_mean=float(w[mask]@v[mask,j]/w[mask].sum()),sampled_minimum=float(np.min(v[mask,j])),sampled_maximum=float(np.max(v[mask,j])))
    return out

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','ct4','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();newprep=read(a.root/'initial/A95/results.json');c=newprep['case'];n=newprep['global_count'];weights=newprep['weights'];baseline=a.ct4/'evolution-v1/A95-6.0-P4/T0.008'
    common=a.root/'recovery-analysis-v1/spatial/0.001/common-basis.json';bp=read(common)
    m=RequestManifest(a.output,dict(schema='CT5EqualTimeAndStageInventoriesV1',common_basis_sha256=sha(common),new_initial_sha256=sha(a.root/'initial/A95/results.json'),source_sha256=sha(Path(__file__)),
        consumed_native_reports={str(f):sha(f) for dt in ('0.002','0.001') for f in (baseline/dt/'case.evolution.json',a.root/'runs/A95'/dt/'T0.008/native/case.evolution.json')},
        baseline_initial_sha256=sha(a.ct4/'initial-v3/A95-6.0-P4/0.001/initial.json')))
    out=dict(status='pass',equal_time={},baseline_stage_scalars=[],baseline_outlet={},scope='read-only stored coefficients, unchanged physical reference, absolute pressure; no residual correction or additional governing solve')
    for dt in ('0.002','0.001'):
        old=read(baseline/dt/'case.evolution.json');newfolder=a.root/'runs/A95'/dt/'T0.008/native';new=read(newfolder/'case.evolution.json');npkt=read(a.root/'initial/A95'/dt/'initial.json');times=sorted(set(old['audit_times'])&set(new['audit_times']))
        out['equal_time'][dt]={}
        for t in times:
            if t not in (0.,.001,.002,.004,.008):continue
            ss=[]
            for report in (old,new):ss.append(report['accepted_states'][next(i for i,v in enumerate(report['audit_times']) if abs(v-t)<1e-13)])
            obs=[]
            for packet,s in zip((bp,npkt),ss,strict=True):v=values(packet,s,n);v['observation']['time']=t;obs.append(v)
            truth,_=expected(npkt,c);truth['observation']={k:v[:n] if isinstance(v,list) else v for k,v in truth['observation'].items()};truth['observation']['time']=t
            gradients=[]
            for packet in obs:
                grad=np.array(packet['observation']['grid_gradients']);d=grad[:,0,0]+grad[:,1,2];w=np.array(weights)
                gradients.append(dict(divergence_weighted_RMS=float(np.sqrt(w@(d*d)/w.sum())),signed_mean=float(w@d/w.sum()),maximum=float(np.max(abs(d)))))
            out['equal_time'][dt][str(t)]=dict(spatial_fields=compare(*obs,weights,c),reference_errors=[compare(v,truth,weights,c) for v in obs],
                gradients=gradients,states=[exact_vector_hash(s) for s in ss],scope='initial projection at0 separate from subsequent evolved contrasts; common physical weights/sides; no order inferred')
    oldprep=read(a.ct4/'initial-v3/A95-6.0-P4/results.json');ini=read(a.ct4/'initial-v3/A95-6.0-P4/0.001/initial.json');native=read(baseline/'0.001/case.evolution.json');oc=oldprep['case']
    for t in (0.,.001,.002,.004,.008):
        i=next(i for i,x in enumerate(native['audit_times']) if abs(x-t)<1e-13);state=native['accepted_states'][i];tr=traction(ini,state,oc)
        out['baseline_outlet'][str(t)]=dict(normalized_reconstructed_outlet_RMS=tr['normalized_outlet_reconstruction_RMS'],
            prescribed_outer_total_traction=oc['reference_column']['outer_equilibrium']['outlet_traction'],prescribed_inner_total_traction=[0.,oc['ustar']**2],
            scope='actual reconstructed outlet mismatch, distinct from the prescribed natural datum; not added to ground/interface accuracy norms')
    for attempt in native['pimple_attempts']:
        if not any(abs(attempt['time']-t)<1e-13 for t in (.001,.002,.004,.008)):continue
        fpath=artifact(attempt['flow_snapshot'],attempt['flow_snapshot_sha256']);kpath=artifact(attempt['k_snapshot'],attempt['k_snapshot_sha256']);f=read(fpath);k=read(kpath)
        states=[('input',f['input']),('predictor',f['momentum_predictor']),('corrected_flow',f['corrected_flow']),('post_flow',k['post_flow']),('post_epsilon',k['pre_k']),('post_k',k['corrected'])]
        pa,rp,owners=pressure_action(ini['observation'],np.array(f['corrected_flow'])-f['momentum_predictor'],oc['roughness_input'],oc['pressure_degree']);ids=np.array(f['rows'],int);nv=f['velocity_size']
        defect=float(np.max(abs(pa[ids[:nv]]-np.array(f['full_pressure_action'])[ids[:nv]])));require(defect<2e-9,'full_rough_pressure_correction_action')
        out['baseline_stage_scalars'].append(dict(time=attempt['time'],flow_snapshot_sha256=sha(fpath),k_snapshot_sha256=sha(kpath),
            scalar_inventories={name:dict(state_sha256=exact_vector_hash(state),**scalar_inventory(ini,state,oldprep['global_count'],oldprep['weights'])) for name,state in states},
            independent_full_pressure_action_defect=defect,rough_pressure_action_maximum=float(np.max(abs(rp))),trace_geometries=sorted({p['trace_geometry'] for p in owners})))
    require(len(out['baseline_stage_scalars'])==4,'missing_required_stage_time')
    write_json(a.output/'results.json',out);m.seal(a.output,['results.json']);print('CT5 equal-time/scalar/outlet postchecks pass');return 0
if __name__=='__main__':raise SystemExit(main())
