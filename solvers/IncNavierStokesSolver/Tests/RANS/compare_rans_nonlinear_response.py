#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Matched baseline-subtracted responses, never a neutral-error ratio gate."""
import argparse,gzip,json,os,shutil
from pathlib import Path
os.environ.update(dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'),'1'))
import numpy as np
from compare_rans_balanced_evolution import compare
from rans_duration_contract import require,identity
from rans_request_manifest import sha,runtime_identity
from rans_qualification import write_json
from rans_preservation_oracle import Reference
from qualify_rans_preservation import load,replay_request
from diagnose_rans_balanced_spatial import run_probe
from analyze_rans_duration_constitutive import analyze


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--perturbed',nargs=4,required=True);p.add_argument('--baseline',nargs=4,required=True);p.add_argument('--time',type=float,required=True);p.add_argument('--output-directory',type=Path,required=True);a=p.parse_args();r=a.root.resolve();out=a.output_directory;out.mkdir(parents=True,exist_ok=False);runtime=r/'runtime-resolution-v3';os.environ.update(LD_LIBRARY_PATH=str(runtime/'lib'),TMPDIR=str(r/'tmp'))
    pert=compare(r,a.perturbed,a.time,inputs_only=True);base=compare(r,a.baseline,a.time,inputs_only=True);requests=pert['requests'];config={k:requests[0][k] for k in ('model','policy','wall','descriptor')};ref=Reference(config['model'],config['policy'],config['wall']);weights=np.array(requests[0]['evaluation']['weights']);arrays=[];initial=[];checks=[];provenance=[]
    for i,(pl,bl) in enumerate(zip(a.perturbed,a.baseline,strict=True)):
        pq,bq=requests[i],base['requests'][i]
        for key in ('model','policy','wall','descriptor','knots','degree','pressure_degree','columns','quadrature','dt','runtime','boundary_policy','boundary','reference_sha256','cache'):
            require(pq[key]==bq[key],'response_baseline_identity_mismatch:'+key)
        require(pq['velocity_amplitude']!=0 and bq['velocity_amplitude']==0,'response_baseline_perturbation_invalid')
        pv,bv=[np.array(bundle['samples'][i]['grid_values']) for bundle in (pert,base)];w=np.array(pert['samples'][i]['grid_weights']);require(np.array_equal(w,weights) and pert['samples'][i]['grid']==base['samples'][i]['grid'],'response_grid_mismatch')
        delta=pv-bv;delta[:,4]-=np.sum(w*delta[:,4])/w.sum();physical=delta/np.array([ref.us,ref.us,ref.K,ref.A/ref.H,ref.us**2]);outer=np.array(pert['samples'][i]['grid_regions'])==1;responses={name:(physical[outer,j],w[outer]) if j==3 else (physical[:,j],w) for j,name in enumerate(('u','w','k','epsilon_outer','Pi'))}
        native=[]
        for label in (pl,bl):
            folder=r/'campaigns'/label;endpoint=folder/('T-'+format(a.time,'.4f'));state=load(endpoint/'case.evolution.json');target=out/(str(i)+'-'+('perturbed' if label==pl else 'baseline'));target.mkdir();shutil.copy2(endpoint/'case.xml',target/'case.xml');seed=load(folder/'initial-packet.json');q=replay_request(seed,pq['evaluation']['grid'],state=state['state'],time=state['time']);q['fixed_boundary']=state['initialization_receipt']['live_fixed_boundary'];q['export_constitutive']=True;packet,_=run_probe(runtime,target,q);packet['observation']['grid_weights']=w.tolist();checks.append(analyze(packet,config,w));native.append(packet)
            provenance.append(dict(label=label,request_sha256=sha(folder/'request-manifest.json'),native_report_sha256=sha(endpoint/'case.evolution.json'),history=load(endpoint/'metrics.json')['history_identity'],consumed_restart_state_sha256=state['initialization_receipt'].get('raw_initial_sha256'),final_state_sha256=load(endpoint/'metrics.json')['final_state_sha256'],packet_sha256=sha(target/'packet.json.gz')))
        pp,bb=[v['native_constitutive']['samples'] for v in native];require(len(pp)==len(bb),'response_constitutive_shape')
        stress=[];sw=[]
        for ps,bs in zip(pp,bb,strict=True):
            require((ps['sample_index'],ps['coordinate'],ps['side'],ps['geometry'])==(bs['sample_index'],bs['coordinate'],bs['side'],bs['geometry']),'response_constitutive_identity')
            R=np.array(ps['R'])-np.array(bs['R']);Q=np.array(ps['Q'])-np.array(bs['Q']);stress.append([R[0,0],R[1,1],R[2,2],R[0,2],Q[0,2],ps['P']-bs['P']]);sw.append(w[ps['sample_index']])
        stress=np.array(stress);sw=np.array(sw)
        for j,name in enumerate(('Rxx','Ryy','Rzz','Rxz','Qxz','P')):responses[name]=(stress[:,j]/(ref.A/ref.H if name=='P' else ref.us**2),sw)
        arrays.append(responses)
        init=load(r/'campaigns'/pl/'actual-initial/observation.json');g=np.array(init['grid_gradients']);div=g[:,0,0]+g[:,1,2];initial.append(dict(label=pl,state_sha256=pq['actual_initial']['state_sha256'],native_history=pq['actual_initial']['history_header'],raw_divergence_RMS=float(np.sqrt(np.sum(w*div**2)/w.sum())),raw_divergence_maximum=float(np.max(np.abs(div))),identity=pq['actual_initial']['receipt']))
    fields={}
    def rms(v,w):return float(np.sqrt(np.sum(w*v*v)/np.sum(w)))
    for name in arrays[0]:
        v=[q[name][0] for q in arrays];w=arrays[0][name][1]
        require(all(np.array_equal(w,q[name][1]) for q in arrays),'response_weights_differ')
        signal=[rms(x,w) for x in v];dh=rms(v[0]-v[1],w);dt=[rms(v[0]-v[2],w),rms(v[1]-v[3],w)];amplitude=max(signal);hlimit=max(5e-7,.1*amplitude);tlimit=max(1e-7,.05*amplitude)
        fields[name]=dict(response_RMS=signal,response_L2=[v*float(np.sqrt(w.sum())) for v in signal],spatial_difference=dh,timestep_differences=dt,spatial_limit=hlimit,temporal_limit=tlimit,signal_status='resolved' if amplitude>5e-7 else 'unresolved_floor',spatial_status='pass' if dh<=hlimit else 'unresolved',temporal_status='pass' if max(dt)<=tlimit else 'unresolved')
    conservation=[]
    for bundle in (pert,base):
        for m in bundle['metrics']:
            d=m['source_and_constraint_maxima'];conservation.append(dict(defects=d,status='pass' if all(abs(v)<=(0 if k=='bound_calls' else 1e-10 if k=='corrected_continuity_gauge' else 1e-9) for k,v in d.items() if k!='refreshed_trace_defect') else 'fail'))
    status='pass' if all(x['spatial_status']=='pass' and x['temporal_status']=='pass' for x in fields.values()) and all(x['status']=='pass' for x in conservation) else 'failed_or_unresolved'
    resolved=all(fields[n]['signal_status']=='resolved' for n in ('Pi','w','Qxz','P'))
    np.savez_compressed(out/'response-vectors.npz',**{str(i)+'_'+name:values[0] for i,row in enumerate(arrays) for name,values in row.items()},**{name+'_weights':value[1] for name,value in arrays[0].items()})
    result=dict(status=status,configuration=config,time=a.time,theta=a.time*ref.us/ref.H,fields=fields,core_pressure_velocity_nonlinear_production_response_resolved=resolved,eligible_for_extension=status=='pass' and resolved,initialization=initial,provenance=provenance,independent_constitutive=checks,conservation=conservation,analytical_reference_errors=[m['reference_error'] for m in pert['metrics']],actual_initial_drifts=[m['change_from_initial'] for m in pert['metrics']],diagnostic_runtime=runtime_identity(runtime/'RANSPreservationProbe'),response_vectors_sha256=sha(out/'response-vectors.npz'),scope='departure from matched numerical baselines; no neutral spatial-error ratio, monotone-decay or temporal-order claim')
    write_json(out/'results.json',result);print(status);return 0 if status=='pass' else 1


if __name__=='__main__':
 import sys
 try:raise SystemExit(main())
 except Exception as e:
  destination=Path(sys.argv[sys.argv.index('--output-directory')+1]);destination.mkdir(parents=True,exist_ok=True)
  if not (destination/'results.json').exists():write_json(destination/'results.json',dict(status='fail',failure=repr(e),scope='preserved incomplete response analysis; native packets and logs retained'))
  raise SystemExit(1)
