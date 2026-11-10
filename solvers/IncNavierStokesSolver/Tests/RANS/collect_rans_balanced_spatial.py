#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Complete compact receipts for the balanced spatial study; missing gates fail."""
import argparse,gzip,hashlib,json,subprocess,sys,platform
import numpy,scipy
from pathlib import Path
from qualify_rans_preservation import load
from rans_qualification import write_json,exact_vector_hash,THREADS
from rans_request_manifest import RequestManifest

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def receipt(path):
 if not path.exists():return dict(status='not_completed',source=str(path))
 return dict(load(path),source=str(path),source_file_sha256=sha(path))
def compact_spatial(s):
 d={k:v for k,v in s.items() if k not in ('one_sided_jumps','constrained_reactions')};byz={}
 for j in s['one_sided_jumps']:byz.setdefault(j['coordinate'][2],[]).append(j)
 d['element_localization']=[dict(row,terms={field:{name:term for name,term in terms.items() if any(value!=0 for value in term.values())} for field,terms in row['terms'].items() if any(any(value!=0 for value in term.values()) for term in terms.values())}) for row in d['element_localization']]
 d['dominant_native_rows_with_basis_support']={k:v[:20] for k,v in d['dominant_native_rows_with_basis_support'].items()}
 d['one_sided_jump_summary']=[]
 for z,items in sorted(byz.items()):
  worst_t=max(items,key=lambda j:sum(x*x for x in j['traction']));worst_k=max(items,key=lambda j:abs(j['k_flux']));w=[j['one_sided'][0]['weight'] for j in items];total=sum(w)
  small=lambda j:{k:v for k,v in j.items() if k!='one_sided'} | {'one_sided':[{k:v for k,v in a.items() if k!='test'} for a in j['one_sided']]}
  d['one_sided_jump_summary'].append(dict(z=z,interface=items[0]['interface'],samples=len(items),traction_RMS=(sum(v*sum(x*x for x in j['traction']) for v,j in zip(w,items,strict=True))/total)**.5,k_flux_RMS=(sum(v*j['k_flux']**2 for v,j in zip(w,items,strict=True))/total)**.5,worst_traction=small(worst_t),worst_k_flux=small(worst_k)))
 d['prescribed_reaction_summary']={k:dict(rows=len(v['rows']),signed_sum=sum(v['load']),absolute_sum=sum(abs(x) for x in v['load']),maximum=max(map(abs,v['load']),default=0)) for k,v in s['constrained_reactions'].items()}
 return d

def compact_analysis(d):
 d=dict(d)
 if 'spatial' in d:d['spatial']=compact_spatial(d['spatial'])
 if 'common_test' in d:
  d['common_test']=dict(d['common_test']);d['common_test']['fields']={name:{k:v for k,v in row.items() if k not in ('restricted_load','common_free_ids')} | dict(restricted_load_sha256=exact_vector_hash(row.get('restricted_load',[])),free_ids_sha256=exact_vector_hash(row.get('common_free_ids',[]))) for name,row in d['common_test']['fields'].items()}
 return d

def audit_summaries(native):
    maxima={};cutoffs=[]
    for t,audit in zip(native.get('audit_times',[]),native.get('state_audits',[]),strict=True):
        for key,values in audit.items():
            if not any(token in key for token in ('interface','traction','flux','trace','gauge','divergence','budget','pilot_accuracy')):continue
            rows=maxima.setdefault(key,[None]*len(values))
            if len(rows)!=len(values):raise ValueError('native diagnostic shape changed')
            for j,v in enumerate(values):
                if v is not None and (rows[j] is None or abs(v)>abs(rows[j]['value'])):rows[j]=dict(value=v,time=t,stage='accepted',component=j)
    for attempt in native.get('pimple_attempts',[]):
        for stage,d in attempt.items():
            if stage.endswith('_domain') and isinstance(d,dict) and 'cutoff_range' in d:cutoffs.append(dict(time=attempt['time'],stage=stage,range=d['cutoff_range'],minimum_required_ratio=d.get('minimum_required_budget_ratio'),required=d.get('required_cutoff_count'),limiting_point=d.get('cutoff_limiting_point')))
    return dict(accepted_componentwise_absolute_maxima=maxima,cutoff_motion=dict(samples=len(cutoffs),minimum=min(cutoffs,key=lambda x:x['range'][0]) if cutoffs else None,maximum=max(cutoffs,key=lambda x:x['range'][1]) if cutoffs else None,first=cutoffs[:1],last=cutoffs[-1:]),scope='all retained accepted-state diagnostic components; all audited stage cutoff ranges; full named vector schema and records remain in the hashed native report')

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=a.root.resolve();out=a.output;out.mkdir(parents=True,exist_ok=False)
 result=dict(schema='BalancedSpatialEvidenceV1',environment=dict(python=sys.version,numpy=numpy.__version__,scipy=scipy.__version__,platform=platform.platform(),threads=THREADS,ranks=1,aggregate_worker_limit=72),starting_revision='3b0961e12142d7eca56d75962feebfffa34630cc',source_report=subprocess.check_output(['git','-C',str(r/'source'),'rev-parse','HEAD'],text=True).strip(),verification={},initial_screens={},initial_gates={},evolved_audits={},required_short_comparisons={},cases=[],duration={},failures=[],scope='serial flat2D prescribed-reference patch; unchanged equations; full native and fixed-common tests')
 for dirname in ('runtime','runtime-resolution-v2'):
  for n in ('source-build.json','actual-loaded-libraries.json','lineage.json'):
   p=r/dirname/n
   if p.exists():(out/(dirname+'-'+n)).write_bytes(p.read_bytes())
 for key,path in {'previous-comparison-contract':r/'evidence-contract/results.json','native-observer':r/'diagnostic-regression/results.json','actual-perturbed-initial':r/'observer-perturbed-baseline-v2/results.json','explicit-parent-selector':r/'common-test-selector-regression/results.json','coarse-quadrature':r/'balanced-quadrature/results.json','fine-quadrature':r/'balanced-fine-quadrature/results.json','evolved-quadrature':r/'evolved-quadrature-0.016/results.json','evaluation-uncertainty':r/'evaluation-uncertainty/results.json','prescribed-controls':r/'prescribed-controls/results.json','mathematical-reuse':r/'mathematical-source-reuse.json','nonlinear-normal-stresses':r/'nonlinear-normal-stresses.json','targeted-ctest':r/'targeted-ctest.json','portable-sampling-witness':r/'portable-sampling-witness/results.json','incomplete-collector-contract':r/'collector-incomplete-contract.json','exclusive-timing':r/'exclusive-timing.json'}.items():
  d=receipt(path)
  # The full evaluation grid remains in the immutable request, not duplicated.
  for row in d.get('checks',[]):row.pop('source_grid',None)
  result['verification'][key]=dict(status=d['status'],source=str(path),sha256=sha(path) if path.exists() else None);write_json(out/(key+'.json'),d)
 previous=r/'evidence-contract/supplement/summary.json'
 if previous.exists():
  d=load(previous);body=d['required_evolved_space_comparisons'];assert len(body)==12
  write_json(out/'previous-twelve-complete-comparisons.json',dict(status='pass' if all(x['status']=='pass' and len(x.get('cases',[]))==3 for x in body.values()) else 'fail',comparisons=body,source=str(previous),source_sha256=sha(previous),exploratory_failures=d.get('exploratory_comparison_failures',[x for x in d.get('evolved_space_comparisons',[]) if x.get('status')!='pass'])))
 else:result['verification']['previous-bodies']=dict(status='not_completed')
 for p in sorted((r/'diagnosis').glob('*/results.json')):
  d=load(p);result['initial_screens'][p.parent.name]=dict(status=d['status'],source=str(p),source_sha256=sha(p),checks=[compact_analysis(x) for x in d['checks']])
 for p in sorted(r.glob('initial-gate-*.json')):result['initial_gates'][p.stem]=receipt(p)
 for p in sorted((r/'evolved-audits').glob('*/results.json')):
  d=load(p);result['evolved_audits'][p.parent.name]=dict(status=d['status'],source=str(p),source_sha256=sha(p),checks=[compact_analysis(x) for x in d['checks']])
 for model in ('E94E95','A95','A93'):
  for wall in ('EP','Accepted'):
   for t in (.0016,.016):
    name=model+'-Natural-'+wall+'-'+str(t);result['required_short_comparisons'][name]=receipt(r/'comparisons'/(name+'.json'))
 for folder in sorted((r/'campaigns').iterdir()):
  p=folder/'request-manifest.json'
  if not p.exists():continue
  manifest=load(p);req=manifest['request'];admission=RequestManifest(folder,req,True);row=dict(label=folder.name,request_sha256=manifest['sha256'],source=str(folder),model=req['model'],policy=req['policy'],wall=req['wall'],descriptor=req['descriptor'],space={k:req[k] for k in ('knots','degree','pressure_degree','columns','quadrature')},dt=req['dt'],horizons=req['horizons'],runtime=req['runtime'],initial_reference_boundary=dict(raw=req['raw_neutral_sha256'],reference=req['reference_sha256'],boundary_sha256=hashlib.sha256(json.dumps(req['boundary'],separators=(',',':')).encode()).hexdigest()),milestones=[])
  for p in sorted(folder.glob('T-*/metrics.json')):
   if not admission.completed(p.parent):continue
   m=load(p);native=load(p.parent/'case.evolution.json')
   receipt_init={k:v for k,v in m.get('initialization_receipt',{}).items() if k not in ('live_fixed_boundary','effective_initial')}
   item={k:v for k,v in m.items() if k in ('status','target','time','theta','exit_code','failure','execution_status','accuracy_status','seconds_per_step','wall_seconds','steps','reference_error','change_from_initial','gradient_weighted_RMS','spatial_alpha_zero','final_stage_residual','source_and_constraint_maxima','history_identity','final_state_sha256')};item['initialization']=receipt_init;item['extrema']={k:v for k,v in m.get('extrema',{}).items() if k!='final_state'};item['final_state_audit']=m.get('final_audit');item['configuration_sha256']=sha(p.parent/'case.xml');item['native_report_sha256']=sha(p.parent/'case.evolution.json');item['energy_budget']=native.get('energy_budget');item['roughness_budget']=native.get('roughness_budget');item['corrections']=dict(calls=sum(x.get('bound_calls',0) for x in native.get('pimple_attempts',[])),activations=sum(x.get('bound_activations',0) for x in native.get('pimple_attempts',[])))
   item['all_state_and_stage_summaries']=audit_summaries(native)
   item['corrections'].update(signed_integral=0. if item['corrections']['calls']==0 else None,absolute_integral=0. if item['corrections']['calls']==0 else None)
   row['milestones'].append(item)
   if m['status']!='pass':
    witness=dict(label=folder.name,metric=m,failure=native.get('failure'),state=native.get('state'),rejected_state=native.get('rejected_state'),history=native.get('native_history'),events=native.get('physical_range_events'))
    with gzip.open(out/(folder.name+'-'+str(m['target'])+'-failure.json.gz'),'wt') as f:json.dump(witness,f)
    (out/(folder.name+'-failure.xml')).write_bytes((p.parent/'case.xml').read_bytes());result['failures'].append(dict(label=folder.name,target=m['target'],failure=m.get('failure')))
  observer=folder/'scalar-observations/results.json'
  if observer.exists():
   obs=load(observer);row['scalar_observations']={k:v for k,v in obs.items() if k not in ('observations',)};row['scalar_observations_source_sha256']=sha(observer)
  result['cases'].append(row)
 for p in sorted((r/'comparisons').glob('*duration*.json')):result['duration'][p.stem]=receipt(p)
 if (r/'duration-admission.json').exists():result['duration_admission']=receipt(r/'duration-admission.json')
 # Preserve first diagnostic/harness failures, independently of later corrections.
 for name in ('replay-previous.log','replay-previous-v2.log','check_rans_observer_baseline.py.log'):
  p=r/name
  if p.exists():(out/name).write_bytes(p.read_bytes())
 result['historical_epsilon_orders']='eight previous observations remain unresolved under their original criterion'
 required_audits=['balanced-'+model+'-N-'+wall+'-'+suffix for model in ('E','A95','A93') for wall in ('EP','Accepted') for suffix in ('B5','B11')]
 result['coverage']=dict(required_short_comparisons=12,present_short_comparisons=sum(x.get('status')=='pass' for x in result['required_short_comparisons'].values()),required_evolved_term_audits=required_audits,complete_evolved_term_audits=[n for n in required_audits if result['evolved_audits'].get(n,{}).get('status')=='pass'])
 required_initial=['initial-gate-'+w for w in ('EP','Accepted')]+['initial-gate-'+m+'-Natural-'+w for m in ('A95','A93') for w in ('EP','Accepted')]
 result['coverage']['required_initial_gates']=required_initial
 required_observers=[m+'-N-'+w+'-'+v for m in ('E','A95','A93') for w in ('EP','Accepted') for v in ('B5','B11')]
 result['coverage']['required_natural_scalar_observers']=required_observers
 bylabel={c['label']:c for c in result['cases']}
 result['coverage']['complete_natural_scalar_observers']=[n for n in required_observers if bylabel.get(n,{}).get('scalar_observations',{}).get('status')=='pass']
 duration_final=[x for k,x in result['duration'].items() if k.endswith('-0.5')]
 result['duration_final_audits']=receipt(r/'duration-final-audits/results.json')
 result['gates']=dict(natural_scalar_observations='pass' if len(result['coverage']['complete_natural_scalar_observers'])==12 else 'failed_or_missing',initial_spatial='pass' if all(result['initial_gates'].get(k,{}).get('status')=='pass' for k in required_initial) else 'failed_or_missing',verification='pass' if all(x['status']=='pass' for x in result['verification'].values()) else 'failed_or_missing',short_spatial_temporal='pass' if result['coverage']['present_short_comparisons']==12 else 'failed_or_missing',all_model_term_replay='pass' if len(result['coverage']['complete_evolved_term_audits'])==12 else 'failed_or_missing',targeted_duration='pass' if duration_final and all(x['status']=='pass' for x in duration_final) and result['duration_final_audits']['status']=='pass' else 'failed_or_unattempted')
 result['status']='pass' if all(x=='pass' for x in result['gates'].values()) else 'failed_unresolved_or_incomplete';write_json(out/'summary.json',result);print(result['status']);return 0 if result['status']=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
