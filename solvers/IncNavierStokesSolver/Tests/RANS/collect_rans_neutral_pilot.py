#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact immutable pilot export. Missing required evidence returns nonzero."""
import argparse,gzip,hashlib,json,os,platform,subprocess,sys
from pathlib import Path
from rans_qualification import THREADS,write_json,coverage
from rans_request_manifest import RequestManifest
from qualify_rans_preservation import load
from rans_preservation_oracle import Reference
from rans_pilot_evidence import required_comparisons,epsilon_budget_gate
import numpy,scipy

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def compact_init(r):return {k:v for k,v in (r or {}).items() if k not in ('effective_initial','live_fixed_boundary')}
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--comparison-directory',type=Path);a=p.parse_args();base=a.root.resolve();r=base/'p5';out=a.output;out.mkdir(parents=True,exist_ok=False)
 result=dict(schema='NeutralPilotEvidenceV1',source_start='1297301156c703fd1e4fe3b7776dc3676a5f3267',source_report=subprocess.check_output(['git','-c','core.filemode=false','-C',str(base/'source'),'rev-parse','HEAD'],text=True).strip(),environment=dict(python=sys.version,numpy=numpy.__version__,scipy=scipy.__version__,platform=platform.platform(),threads=THREADS,ranks=1),cases=[],verification={},failures=[])
 for runtime in (base/'runtime',r/'runtime'):
  for name in ('source-build.json','actual-loaded-libraries.json'):
   if (runtime/name).exists():
    dest=out/(("p5-" if runtime.parent==r else "first-")+name);dest.write_bytes((runtime/name).read_bytes())
 for key,path in [('frozen-boundary',base/'frozen-boundary-regression/results.json'),('interface-side',base/'interface-side-regression-v2/results.json'),('native-controls',r/'native-regression/qualification.json'),('validated-resume',r/'resume-cli-verification.json'),('ground-spectrum-domain',base/'ground-polynomial-unit.json'),('mathematical-reuse',base/'mathematical-evidence-reuse.json'),('extended-ctest',base/'extended-ctest.json'),('velocity-pressure',r/'velocity-pressure-response.json'),('streamwise',r/'streamwise-sensitivity-v2/results.json'),('selected-controls',r/'selected-controls/results.json'),('selected-quadrature',r/'selected-quadrature/results.json'),('nearwall-readonly',r/'nearwall-frozen-diagnosis/results.json'),('failure-contract',base/'collector-failure-contract.json')]:
  if path.exists():
   d=load(path)
   if key=='frozen-boundary':
    for c in d['checks']:
     if 'initialization_receipts' in c:c['initialization_receipts']=[compact_init(x) for x in c['initialization_receipts']]
   if key=='native-controls':d={k:d[k] for k in ('status','checks','required_check_names','coverage') if k in d}
   if key=='selected-controls':
    for c in d['checks']:
     v=c.get('result',{})
     if 'metrics' in v:
      m=v['metrics'];request=v['request'];c['result']=dict(label=v['label'],request_sha256=request['sha256'],descriptor=request['request']['descriptor'],reference_error=m.get('reference_error'),initialization=compact_init(m.get('initialization_receipt')),history=m.get('history_identity'),source_and_constraint_maxima=m.get('source_and_constraint_maxima'))
   result['verification'][key]=dict(status=d.get('status'),source=str(path),sha256=sha(path));write_json(out/(key+'.json'),d)
  else:result['verification'][key]=dict(status='not_attempted',source=str(path))
 for folder in sorted((r/'campaigns').iterdir()):
  if not (folder/'request-manifest.json').exists():continue
  payload=load(folder/'request-manifest.json');req=payload['request'];admission=RequestManifest(folder,req,True);ref=Reference(req['model'],req['policy'],req['wall'])
  row=dict(label=folder.name,model=req['model'],policy=req['policy'],wall=req['wall'],request_sha256=payload['sha256'],configuration={k:req[k] for k in ('descriptor','knots','degree','pressure_degree','columns','quadrature','dt','velocity_amplitude','boundary_policy','raw_neutral_sha256','reference_sha256')},runtime_sha256=req['runtime']['IncNavierStokesSolver']['sha256'],boundary_identity_sha256=hashlib.sha256(json.dumps(req['boundary'],separators=(',',':')).encode()).hexdigest(),source=str(folder),milestones=[])
  ini=load(folder/'initial-metrics.json') if (folder/'initial-metrics.json').exists() else {};row['initial']=dict(dof=ini.get('dof'),accuracy=ini.get('accuracy'),alpha_zero=ini.get('residual'),identity=compact_init(ini.get('identity')),aerodynamic_roughness=ref.H*__import__('math').exp(-ref.kappa*ref.UH/ref.us),wall_matching_velocity=ref.UH,analytical=ini.get('analytical'))
  for f in sorted(folder.glob('T-*/metrics.json')):
   if not admission.completed(f.parent):continue
   m=load(f);rep=load(f.parent/'case.evolution.json') if (f.parent/'case.evolution.json').exists() else {};attempts=rep.get('pimple_attempts',[]);audits=rep.get('state_audits',[]);times=rep.get('audit_times',[])
   item={k:m.get(k) for k in ('target','time','theta','status','failure','execution_status','accuracy_status','reference_error','change_from_initial','gradient_weighted_RMS','spatial_alpha_zero','final_stage_residual','source_and_constraint_maxima','history_identity','final_state_sha256','seconds_per_step','wall_seconds','steps')};item['initialization']=compact_init(m.get('initialization_receipt'));item['native_report_sha256']=sha(f.parent/'case.evolution.json') if rep else None
   calls=sum(x.get('bound_calls',0) for x in attempts);item['corrections']=dict(calls=calls,activations=sum(x.get('bound_activations',0) for x in attempts),signed_integral=0. if calls==0 else None,absolute_integral=0. if calls==0 else None)
   item['extrema']={k:v for k,v in m.get('extrema',{}).items() if k!='final_state'};item['final_extrema']={k:v for k,v in (m.get('extrema',{}).get('final_state') or {}).items() if any(t in k for t in ('domain','covariance','radicand','interface','trace','cutoff','continuity'))};item['accuracy_maxima']={};item['cutoff_motion']=[]
   for t,d in zip(times,audits,strict=True):
    vals=d.get('pilot_accuracy_E_D_wall_gauge',[])
    for j,name in enumerate(('u','v','k','epsilon_outer','Pi')):
     if vals and vals[j]>item['accuracy_maxima'].get(name,{}).get('value',-1):item['accuracy_maxima'][name]=dict(value=vals[j],time=t,stage='accepted')
   for step in attempts:
    for stage,d in step.items():
     if stage.endswith('_domain') and isinstance(d,dict) and 'cutoff_range' in d:item['cutoff_motion'].append(dict(time=step['time'],stage=stage,range=d['cutoff_range'],required=d.get('required_cutoff_count'),budget_minimum=d.get('minimum_required_budget_ratio')))
   motion=item.pop('cutoff_motion');item['cutoff_motion_summary']=dict(samples=len(motion),minimum=min(motion,key=lambda x:x['range'][0]) if motion else None,maximum=max(motion,key=lambda x:x['range'][1]) if motion else None,first=motion[:1],last=motion[-1:])
   item['momentum_k_weak_load_reaction']=audits[-1].get('momentum_k_load_and_boundary_reaction') if audits else None
   item['energy_budget']=rep.get('energy_budget');item['roughness_budget']=rep.get('roughness_budget')
   row['milestones'].append(item)
   if m['status']=='fail':
    witness=dict(label=folder.name,target=m['target'],failure=m.get('failure'),committed_time=rep.get('time'),state=rep.get('state'),history=rep.get('native_history'),history_times=rep.get('native_history_times'),rejected_state=rep.get('rejected_state'),events=rep.get('physical_range_events'),signature=rep.get('signature'),configuration_sha256=sha(f.parent/'case.xml'))
    with gzip.open(out/(folder.name+'-failure.json.gz'),'wt') as stream:json.dump(witness,stream)
    (out/(folder.name+'-failure.xml')).write_bytes((f.parent/'case.xml').read_bytes());result['failures'].append({k:witness[k] for k in ('label','target','failure','committed_time')})
  obs=folder/'scalar-observations/results.json'
  if obs.exists():
   d=load(obs);row['read_only_scalar_stage_observations']={k:d[k] for k in ('status','scope','minimum','maximum','epsilon_weak_budgets','analytical_weak_balance','source_packet_sha256')};row['observer_receipt_sha256']=sha(obs);row['read_only_scalar_stage_observations']['epsilon_budget_gate']=epsilon_budget_gate(d['epsilon_weak_budgets'],ref)
  row['status']=load(folder/'results.json').get('status') if (folder/'results.json').exists() else 'incomplete';result['cases'].append(row)
 first=base/'campaigns/E-N-EP-h0p5/T-0.0016/case.evolution.json'
 if first.exists():
  d=load(first);witness=dict(label='first-p5-cubic-guard',source=str(first),failure=d.get('failure'),state=d.get('state'),history=d.get('native_history'),time=d.get('time'),steps=d.get('steps'),signature=d.get('signature'))
  with gzip.open(out/'first-p5-cubic-guard.json.gz','wt') as stream:json.dump(witness,stream)
  result['failures'].append({k:witness[k] for k in ('label','source','failure','time','steps')})
 comparisons=[]
 for f in sorted(r.glob('*-space-*.json')):
  d=load(f);comparisons.append(dict(path=str(f),sha256=sha(f),**d))
 result['evolved_space_comparisons']=comparisons
 result['exploratory_comparison_failures']=[x for x in comparisons if x.get('status')!='pass']
 expected={('E' if model=='E94E95' else model)+'-N-'+wall+'-h0p5' for model in ('E94E95','A95','A93') for wall in ('EP','Accepted')};actual={x['label']:x for x in result['cases']};gates={}
 for label in sorted(expected):
  c=actual.get(label);milestones={m['target']:m for m in c['milestones']} if c else {};gates[label]=dict(status='pass' if all(t in milestones and milestones[t]['status']=='pass' for t in (.0016,.016,.05,.1)) else 'failed' if c and any(m['status']=='fail' for m in c['milestones']) else 'not_completed',completed_times=sorted(milestones))
 result['required_evolved_space_comparisons']=required_comparisons(r,a.comparison_directory)
 result['required_scalar_stage_observations']={label:actual.get(label,{}).get('read_only_scalar_stage_observations',{}).get('status','not_completed') if actual.get(label,{}).get('read_only_scalar_stage_observations',{}).get('epsilon_budget_gate',{}).get('status')=='pass' else 'failed_or_missing_epsilon_budget' for label in expected}
 result['required_natural_pilots']=gates;result['status']='pass' if all(v['status']=='pass' for v in result['verification'].values()) and all(v['status']=='pass' for v in gates.values()) and all(v.get('status')=='pass' and 'cases' in v for v in result['required_evolved_space_comparisons'].values()) and all(v=='pass' for v in result['required_scalar_stage_observations'].values()) else 'incomplete_or_failed'
 result['full_spatial_residual_consistency_status']='not fully demonstrated in evolved spaces; accepted inner-k and momentum dual norms require balanced near-wall/interface refinement; global h1/p5 frozen controls improve them but 3801 unknowns exceed the 3000 evolution guard'
 result['temporal_scope']='all six natural configurations have dt/2 uncertainty checks at .0016 and .016; .1 is bounded integration with a finer E94E95 accepted-wall comparator, not a measured .1 temporal order'
 result['historical_epsilon_order_status']='eight previous observations remain unresolved; no retrospective threshold change';result['scope']='serial flat2D prescribed-reference bounded patch; no long-time, downstream ABL, curved/3D/MPI or physical validation claim'
 write_json(out/'summary.json',result);print(result['status']);return 0 if result['status']=='pass' else 1
if __name__=='__main__':sys.exit(main())
