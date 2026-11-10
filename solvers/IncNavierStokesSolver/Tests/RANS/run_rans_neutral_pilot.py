#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded, identity-checked neutral pilot; stop at the first native accuracy exit.

--max-time limits work on an immutable .1-horizon request; --resume validates
all completed artifacts and advances only missing milestones with native history.
No field, checkpoint signature, timestep or boundary vector is rewritten.
"""
import argparse,gzip,json,math,os,subprocess,sys,tempfile,time,xml.etree.ElementTree as ET
from pathlib import Path
from contextlib import nullcontext
import numpy as np
from rans_qualification import THREADS,write_json,exact_vector_hash,compare_samples,audited_extrema
from rans_request_manifest import RequestManifest,runtime_identity
from rans_preservation_oracle import Reference,physical_grid,representation
from qualify_rans_preservation import load,probe,replay_request,nested,Q,residual_metrics
from run_a93_a95_e94e95_matrix import case_xml
from run_cls_epawf_picard import set_info,write,sha
from run_cls_epawf_pimple import identity
from continue_cls_epawf_pimple import header

FIELDS=('u','v','k','epsilon_outer','Pi')
TARGETS=(.0016,.016,.05,.1)

def configure(folder,a,ref,steps,restart=None):
 case_xml(folder,a.model,a.policy,a.wall,steps=steps,dt=a.dt,restart=restart,descriptor_probe=a.root/'runtime/RANSNewModelsProbe',radial_knots=nested(ref.knots(),a.level),degree=a.degree,quadrature=Q[a.quadrature_index],checkpoint_interval=max(steps,1),cache=a.cache,boundary_policy='FrozenReferenceV1',xknots=[-1.,0.,1.] if a.columns==2 else [-1.,1.])
 root=ET.parse(folder/'case.xml').getroot();set_info(root,'RANSWallPilotAccuracy','NeutralV1')
 if a.velocity_amplitude:
  # Streamfunction a*(1-x*x)^2*f(z). Both f,f' vanish at each
  # wall/interface/top and g,g' at x boundaries. p>=5 represents it.
  for region,f,fp in [('Inner','y^2*(2-y)^2','2*y*(2-y)*(2-2*y)'),('Outer','(y-2)^2*(4-y)^2','2*(y-2)*(4-y)*(6-2*y)')]:
   fun=ET.SubElement(root.find('CONDITIONS'),'FUNCTION',{'NAME':'InitialPerturbation'+region})
   ET.SubElement(fun,'E',{'VAR':'u','VALUE':str(a.velocity_amplitude)+'*(1-x*x)^2*('+fp+')'})
   ET.SubElement(fun,'E',{'VAR':'v','VALUE':str(4*a.velocity_amplitude)+'*x*(1-x*x)*('+f+')'})
 write(root,folder/'case.xml')

def scaled(comparison,ref):
 scales=(ref.us,ref.us,ref.K,ref.A/ref.H,ref.us**2)
 return {f:dict(weighted_RMS=comparison['fields'][f]['weighted_RMS']/s,L2=comparison['fields'][f]['L2']/s) for f,s in zip(FIELDS,scales,strict=True)}

def accuracy(ref,packet,initial,weights):
 rep=representation(ref,packet,weights);base=dict(initial,time=packet['time'])
 return dict(reference_error=scaled(rep,ref),change_from_initial=scaled(compare_samples(packet,base,weights),ref),gradient_weighted_RMS=rep['gradient_weighted_RMS'],global_pressure_gauge=rep['global_pressure_gauge_shifts'])

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--label',required=True);p.add_argument('--model',choices=('E94E95','A95','A93'),default='E94E95');p.add_argument('--policy',choices=('Natural','Prescribed'),default='Natural');p.add_argument('--wall',choices=('EP','Accepted'),required=True);p.add_argument('--degree',type=int,default=5);p.add_argument('--level',type=int,default=0);p.add_argument('--columns',type=int,choices=(1,2),default=1);p.add_argument('--dt',type=float,default=.0001);p.add_argument('--max-time',type=float,default=.0016);p.add_argument('--quadrature-index',type=int,choices=(0,1,2),default=1);p.add_argument('--velocity-amplitude',type=float,default=0.);p.add_argument('--cache',choices=('Off','ExactStepV1'),default='ExactStepV1');p.add_argument('--resume',action='store_true');p.add_argument('--timeout',type=int,default=3600);a=p.parse_args();a.root=a.root.resolve()
 if a.degree<4 or a.level<0:raise ValueError('declared pilot candidate requires degree >=4, nonnegative refinement')
 r=a.root;out=r/'campaigns'/a.label;ref=Reference(a.model,a.policy,a.wall)
 os.environ.update(THREADS,LD_LIBRARY_PATH=str(r/'runtime/lib'),TMPDIR=str(r/'tmp'))
 # One common physical grid resolves both declared candidates and profile breaks.
 grid,weights=physical_grid(nested(ref.knots(),max(1,a.level)),10,8)
 sampler=r/'runtime/RANSPreservationProbe';solver=r/'runtime/IncNavierStokesSolver'
 with nullcontext(tempfile.mkdtemp(prefix='pilot-request-',dir=r/'tmp')) as tmp:
  pre=Path(tmp)/'seed';configure(pre,a,ref,0);seed=probe(sampler,pre,dict(operation='neutral',grid=grid,export_basis=False))
  initial=dict(seed['observation'],grid_weights=weights)
  request=dict(schema='NeutralPilotV1',model=a.model,policy=a.policy,wall=a.wall,descriptor=seed['observation']['resolved_descriptor'],
   knots=nested(ref.knots(),a.level),degree=a.degree,pressure_degree=a.degree-2,columns=a.columns,quadrature=seed['observation']['actual_quadrature'],
   dt=a.dt,horizons=TARGETS,integration='native IMEX1 PimpleSingleSweep',history='native same-method restart; independently initialized for each dt',
   initialization='constrained neutral projection plus declared divergence-free streamfunction',velocity_amplitude=a.velocity_amplitude,perturbation='StreamfunctionPolyV1',
   raw_neutral_sha256=exact_vector_hash(initial['state']),reference_sha256=exact_vector_hash(seed['reference']),boundary_policy='FrozenReferenceV1',boundary=seed['fixed_boundary'],
   runtime={n:runtime_identity(r/'runtime'/n) for n in ('IncNavierStokesSolver','RANSPreservationProbe','RANSNewModelsProbe')},
   evaluation=dict(grid=grid,weights=weights,sides=[0]*len(grid),export_basis=False,metric='physical weighted RMS and true L2; common global pressure gauge'),
   gates=dict(field_RMS=1e-3,vertical_RMS=1e-4,wall_traction=1e-3,interface_traction=1e-3,spatial_ratio=.85,temporal_fraction=.1,comparison_floor=1e-7),cache=a.cache,accuracy_policy='NeutralV1')
  admission=RequestManifest(out,request,a.resume)
  if not a.resume:
   write_json(out/'initial-packet.json',seed);write_json(out/'evaluation-grid.json',dict(grid=grid,weights=weights));write_json(out/'initial-metrics.json',dict(accuracy=accuracy(ref,initial,initial,weights),residual=initial['audit_alpha_zero'],identity=seed['initialization_receipt'],dof=initial['sizes'],analytical=ref.self_check()))
  result=load(out/'results.json') if a.resume and (out/'results.json').exists() else dict(status='running',execution_status='running',numerical_accuracy_status='not_assessed',milestones=[],request_sha256=admission.payload['sha256'])
  if result.get('first_failure'):
   print(json.dumps(dict(status='preserved_failure',failure=result['first_failure'])));return 1
  restart=None;start=0.
  for T in TARGETS:
   if T>a.max_time+1e-14:break
   label='T-'+format(T,'.4f');folder=out/label
   if a.resume and admission.completed(folder):
    m=load(folder/'metrics.json');restart=Path(m['checkpoint']);start=T;continue
   steps=round((T-start)/a.dt)
   if steps<=0 or abs(steps*a.dt-(T-start))>1e-13:raise ValueError('milestone is not an integer number of unchanged timesteps')
   configure(folder,a,ref,steps,restart);env=dict(os.environ,NEKTAR_RANS_TIMING_FILE=str(folder/'timing.jsonl'));t0=time.monotonic()
   with (folder/'solver.log').open('w') as f:
    try:c=subprocess.run([str(solver),'case.xml'],cwd=folder,stdout=f,stderr=subprocess.STDOUT,env=env,timeout=a.timeout);code=c.returncode
    except subprocess.TimeoutExpired:code=124
   elapsed=time.monotonic()-t0;report=load(folder/'case.evolution.json') if (folder/'case.evolution.json').exists() else {}
   complete=code==0 and not report.get('failure') and abs(report.get('time',-1)-T)<1e-13 and report.get('steps')==round(T/a.dt)
   metric=dict(time=report.get('time'),theta=report.get('time',0)*ref.us/ref.H,target=T,exit_code=code,wall_seconds=elapsed,seconds_per_step=elapsed/max(1,steps),steps=report.get('steps'),execution_status='completed' if complete else 'failed',
    initialization_receipt=report.get('initialization_receipt'),extrema=audited_extrema(report),failure=report.get('failure'),configuration_sha256=sha(folder/'case.xml'))
   audits=report.get('state_audits',[]);metric['all_audited_accuracy']=[dict(time=t,values=x['pilot_accuracy_E_D_wall_gauge']) for t,x in zip(report.get('audit_times',[]),audits,strict=True) if 'pilot_accuracy_E_D_wall_gauge' in x]
   metric['final_stage_residual']=report.get('stage_audits',[])[-1:];metric['final_audit']=audits[-1] if audits else None
   attempts=report.get('pimple_attempts',[]);metric['source_and_constraint_maxima']={k:max((float(np.max(np.abs(x.get(k,0)))) for x in attempts),default=0.) for k in ('refreshed_trace_defect','momentum_splitting_identity_defect','k_same_state_source_storage_defect','epsilon_same_state_source_storage_defect','corrected_continuity_gauge','bound_calls')}
   metric['energy_budget']=report.get('energy_budget');metric['roughness_budget']=report.get('roughness_budget')
   if complete:
    req=replay_request(seed,grid,state=report['state'],time=T);req['fixed_boundary']=report['initialization_receipt']['live_fixed_boundary']
    final=probe(sampler,folder,req);packet=dict(final['observation'],grid_weights=weights)
    # Perturbed initial state must be sampled from the actual initialized vector.
    if a.velocity_amplitude:
     first=(out/'T-0.0016'/'case.evolution.json');ir=load(first) if first.exists() else report
     initreq=replay_request(seed,grid,state=ir['accepted_states'][0],time=0.);ip=probe(sampler,folder,initreq);initial=dict(ip['observation'],grid_weights=weights)
    metric.update(accuracy(ref,packet,initial,weights));metric['spatial_alpha_zero']=residual_metrics(packet['audit_alpha_zero']);metric['spatial_alpha_zero_full']=packet['audit_alpha_zero'];write_json(folder/'endpoint-samples.json',packet)
    checks=[metric['reference_error'][f]['weighted_RMS']<=(1e-4 if f=='v' else 1e-3) for f in FIELDS]
    metric['accuracy_status']='pass' if all(checks) else 'fail';metric['checkpoint']=str(folder/('checkpoint.'+str(round(T/a.dt))))
    # The native checkpoint naming is authoritative, not guessed from IO steps.
    checkpoint=Path(report['checkpoint']);assert abs(header(checkpoint)['time']-T)<1e-13
    metric['checkpoint']=str(checkpoint);metric['history_identity']=header(checkpoint);metric['final_state_sha256']=exact_vector_hash(report['state'])
    if a.dt!=metric['history_identity']['dt']:raise AssertionError('history timestep changed')
    restart=checkpoint;start=T
   else:metric['accuracy_status']='failed_or_not_reached';metric['first_failure_witness']=dict(rejected_state=report.get('rejected_state'),range_events=report.get('physical_range_events'),native_history=report.get('native_history'))
   metric['status']='pass' if complete and metric['accuracy_status']=='pass' else 'fail';write_json(folder/'metrics.json',metric)
   files=['case.xml','solver.log','metrics.json']+(['case.evolution.json'] if report else [])+(['endpoint-samples.json'] if complete else [])
   if complete:files+=[str(Path(metric['checkpoint']).relative_to(folder)/'State.bin')]
   admission.seal(folder,files)
   result['milestones'].append(metric);result['execution_status']=metric['execution_status'];result['numerical_accuracy_status']=metric['accuracy_status'];result['status']=metric['status'];write_json(out/'results.json',result)
   print(json.dumps({k:metric[k] for k in ('target','status','seconds_per_step','failure')}),flush=True)
   if metric['status']!='pass':result['first_failure']=dict(folder=str(folder),failure=metric.get('failure'),time=metric.get('time'));write_json(out/'results.json',result);return 1
  result['status']='pass' if result['milestones'] else 'unresolved';result['persistence_scope']=dict(last_time=start,theta=start*ref.us/ref.H,requested_final_time=.1,bounded=True);write_json(out/'results.json',result);return 0 if result['status']=='pass' else 2
if __name__=='__main__':sys.exit(main())
