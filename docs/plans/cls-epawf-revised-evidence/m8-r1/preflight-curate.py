"""Read-only preflight synthesis and compact exact-byte curation; no solver."""
import fcntl, hashlib, json, math, os
from pathlib import Path
import shutil, struct, subprocess, sys
from datetime import datetime, timezone

M8=Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8')
R1=M8.parent/'m8-r1'; SOURCE=R1/'archive-source'; HIST=M8/'initialization'
INITIAL=R1/'preflight-001'; COMPLETE=R1/'preflight-completion-001'
OUT=R1/'preflight-postprocess-001'; OUT.mkdir(exist_ok=False)
CURATED=OUT/'curated'; CURATED.mkdir()
lock=(HIST/'validation.lock').open('a');fcntl.flock(lock,fcntl.LOCK_SH|fcntl.LOCK_NB)
runner=json.loads((COMPLETE/'runner.json').read_text());assert runner['status']=='passed'
for key,value in runner['environment'].items():os.environ[key]=value
os.environ['PYTHONDONTWRITEBYTECODE']='1';os.sched_setaffinity(0,{49})
fingerprint=M8.parent/'transfers/fingerprint-runtime.py'

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def vector_hash(value):return hashlib.sha256(struct.pack('<'+str(len(value))+'d',*value)).hexdigest()
def save(name,value):(CURATED/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
subprocess.run(['python3','-B',str(fingerprint),'--root',str(HIST),'--output',str(OUT/'runtime-before.json'),
    '--compare',str(COMPLETE/'runtime-after.json')],check=True)
sys.path.insert(0,str(SOURCE/'solvers/IncNavierStokesSolver/Tests/RANS'))
import numpy as np
from cls_epawf_m8_state_coverage import reconstruction_rules,point_sets,local_values
from validate_cls_epawf_m4 import spaces
FILES={
 'preflight-initial-attempt.json':INITIAL/'runner.json',
 'preflight-runner.json':COMPLETE/'runner.json',
 'preflight-plan.json':INITIAL/'cases/preflight-plan.json',
 'preflight-casegen-001.py':R1/'preflight-tools-001/cls_epawf_m8_r1_cases.py',
 'preflight-casegen-002.py':R1/'preflight-tools-002/cls_epawf_m8_r1_cases.py',
 'preflight-audit.py':R1/'preflight-tools-002/cls_epawf_m8_r1_preflight.py',
 'preflight-run-001.py':R1/'preflight-tools-001/m8-r1-run-preflight-001.py',
 'preflight-complete-001.py':R1/'preflight-tools-002/m8-r1-run-preflight-completion-001.py',
 'preflight-original-metadata-failure.json':INITIAL/'cases/original-base/preflight.json',
 'preflight-curate.py':Path(__file__).resolve(),
}
audits,raws,rows,oracle_pins={}, {}, [], {}
for case in runner['cases']:
 name=case['name'];audit_path=Path(case['audit']);raw_path=Path(case['report'])
 assert sha(audit_path)==case['audit_sha256'] and sha(raw_path)==case['report_sha256']
 audit=json.loads(audit_path.read_text());raw=json.loads(raw_path.read_text());answer=audit['result']
 assert audit['status']=='passed' and audit['input_hashes_unchanged']
 for module in audit['imported_oracles'].values():
  assert sha(module['path'])==module['sha256']
  oracle_pins[module['path']]=module['sha256']
 audits[name],raws[name]=audit,raw
 FILES['preflight-'+name+'.json']=audit_path;FILES['preflight-'+name+'.xml']=Path(case['xml'])
 source=answer['independent_source_energy']['robust_source'];ranges=answer['finite_ranges']
 commands=runner['commands']
 seconds=case.get('original_solver_seconds')
 if seconds is None:seconds=next(row['seconds'] for row in commands if row['name']=='solve-'+name)
 rows.append(dict(name=name,unknowns=answer['unknowns'],source_quadrature=case['source_quadrature'],
  native_initialization_seconds=seconds,independent_checks=answer['independent_checks'],
  finite_admissibility=answer['initial_finite_admissibility'],
  outer_gradient_entry_deviation=ranges['outer_maximum_gradient_entry_deviation'],
  outer_k_over_K0=ranges['outer_k_over_K0'],outer_epsilon_over_epsilon0=ranges['outer_epsilon_over_epsilon0'],
  inner_maximum_relative_AWF_gradient_error=max(row['maximum_relative_AWF_gradient_error'] for row in ranges['inner'].values()),
  inner_maximum_r_R=max(row['maximum_r_R'] for row in ranges['inner'].values()),
  chi0_trace_epsilon_ratio=answer['independent_trace']['chi_zero_epsilon_ratio_range'],
  source_layout=['preE','E','D','C','limiter_active_volume','cap_active_volume','minimum_lambda','maximum_rate'],
  source=source,trace_sample_b_range=[min(row['b'] for row in answer['independent_trace']['samples']),max(row['b'] for row in answer['independent_trace']['samples'])],
  trace_error_over_E0=answer['independent_trace']['relative_to_E0'],
  raw_unanchored_block_norms=answer['raw_unanchored_residual']['independent_block_norms'],
  complete_weak_load_error=answer['raw_unanchored_residual']['maximum_complete_weak_load_error'],
  interface=answer['independent_interface'],
  same_state_quadrature_sensitivity=answer['same_state_quadrature_sensitivity'],
  initial_state_float64_LE_sha256=answer['coefficient_float64_LE_sha256'],
  native_history_slot_float64_LE_sha256=[vector_hash(slot) for slot in raw['native_history']],
  native_history_times=raw['native_history_times'],order=raw['order'],dt=raw['dt'],
  report_sha256=sha(raw_path),XML_sha256=sha(case['xml']),audit_sha256=sha(audit_path),
  historical_radius001_kEPS99_101_covered=ranges['historical_radius001_kEPS99_101_covered']))
comparisons={}
for representation in ('original','resolved'):
 a,b=raws[representation+'-base'],raws[representation+'-high']
 first,second=a['accepted_states'][0],b['accepted_states'][0]
 same=vector_hash(first)==vector_hash(second)
 comparisons[representation]=dict(initial_coefficients_bitwise_identical=same,
  maximum_initial_coefficient_difference=max(abs(x-y) for x,y in zip(first,second)),
  native_initial_raw_load_maximum_difference=max(abs(x-y) for x,y in zip(a['raw_reference_loads'],b['raw_reference_loads'])),
  native_initial_block_norm_differences=[y-x for x,y in zip(a['raw_reference_norms'],b['raw_reference_norms'])],
  interpretation='Quadrature changes without changing space or physical inputs; observed integration sensitivity, no order or convergence claim')
original_raw=M8/'flow-001/qualification/control-coarse/case.evolution.json'
old=json.loads(original_raw.read_text());new=raws['original-base']
original_match=dict(archived_report_sha256=sha(original_raw),new_report_sha256=sha(INITIAL/'cases/original-base/case.evolution.json'),
 initial_coefficients_bitwise_identical=vector_hash(old['accepted_states'][0])==vector_hash(new['accepted_states'][0]),
 maximum_initial_coefficient_difference=max(abs(x-y) for x,y in zip(old['accepted_states'][0],new['accepted_states'][0])),
 exact_physics_signature_preserved=old['signature']==new['signature'],
 native_history_interpretation='Only t=0 coefficients are compared. No old evolving history or rejected trajectory is reused.')
witnesses={}
for name in ('original-base','resolved-base'):
 raw=raws[name];fields=spaces(raw);state=np.asarray(raw['accepted_states'][0]);rules,issue=reconstruction_rules(raw,INITIAL/'cases'/name/'case.xml');assert issue is None
 flow,normal=np.asarray(raw['robust_reference_flow']),np.asarray(raw['robust_reference_normal'])
 frame=np.column_stack([flow,np.cross(normal,flow),normal]);center=np.zeros((3,3));center[0,2]=1.
 best=None;us,kappa=raw['ustar'],raw['parameters'][7]
 for kind,identity,points in point_sets(raw,fields,rules):
  G,k,epsilon,distance=local_values(raw,fields,state,identity,points)
  A=np.einsum('ai,qab,bj->qij',frame,G,frame)/(us/(kappa*distance))[:,None,None]
  deviations=np.max(np.abs(A-center),axis=(1,2));index=int(np.argmax(deviations))
  value=float(deviations[index])
  if best is None or value>best['gradient_entry_deviation']:
   best=dict(set=kind,geometry=identity,chart=points[index].tolist(),gradient=G[index].tolist(),
    normalized_gradient=A[index].tolist(),gradient_entry_deviation=value,
    k_over_K0=float(k[index]/raw['parameters'][5]),epsilon_over_epsilon0=float(epsilon[index]/(us**3/(kappa*distance[index]))))
 assert abs(best['gradient_entry_deviation']-audits[name]['result']['finite_ranges']['outer_maximum_gradient_entry_deviation'])<1e-14
 witnesses[name]=best
runtime_summaries={}
for prefix in (INITIAL,COMPLETE):
 for point in ('before','after'):
  p=prefix/('runtime-'+point+'.json');raw=json.loads(p.read_text())
  runtime_summaries[str(p)]=dict(sha256=sha(p),success=raw['success'],source_commit=raw['source_commit'],
   source_inputs=len(raw['source_inputs']),runtime_images=len(raw['runtime_artifacts']),
   runtime_environment=raw['runtime_environment'],identical_qualification_inputs_and_runtime=raw.get('identical_qualification_inputs_and_runtime'))
copied=[]
for name,path in FILES.items():
 data=path.read_bytes();(CURATED/name).write_bytes(data)
 assert sha(CURATED/name)==sha(path)
 copied.append(dict(file=name,original=str(path),bytes=len(data),sha256=sha(path),transformation='byte-for-byte copy'))
summary=dict(schema=1,phase='C',status='implementation_audits_passed',generated_utc=datetime.now(timezone.utc).isoformat(),
 source=runner['reviewed_source'],actual_build_source=runner['actual_build_source'],solver=runner['solver'],
 Cxx_source_equivalence=runner['reviewed_ab5_Cxx_identical_to_actual_compiled_ed34'],
 scope='Exactly four zero-step native initializations; no evolution, no field floors, stress clipping or physical/model retuning',
 runtime=dict(numerical_threads=runner['numerical_threads'],affinity=[49],maximum_workers=1),
 mesh=dict(tangential=[-4,-3,-2,-1,1,2],original_radial=[0,2,4],
  resolved_radial=[0,.01,.1,.25,.5,1,2,8/3,10/3,4],degree=3,pressure_degree=1,native_GLL_points=[8,8],
  dense_guard=3000,preserved=True),cases=rows,quadrature_comparisons=comparisons,
 archived_initial_correspondence=original_match,outer_gradient_witnesses=witnesses,
 decision=dict(original='The actual archived mesh has initial negative inner covariance with undepleted k; exact zero-step/archive coefficient correspondence is recorded separately.',
  resolved='The one approved refinement yields independently positive finite inner/outer covariance, positive k/EPS and inactive source limiter/cap.',
  certified_coverage='FAILED: resolved gradient-entry deviation 0.0034123386947804146 exceeds the already frozen radius .001 region; no mesh or box adjustment.',
  next_control='Any parent-authorized resolved trajectory must be diagnostic/unqualified from t=0 and use fresh valid IMEX initialization.',
  excluded=['polynomial-everywhere admissibility','continuum invariant-region proof','adopted model qualification','physical terrain validation']),
 metadata_attempt=dict(preserved=str(INITIAL/'runner.json'),reason='Tuple/list serialization mismatch in time-integration plan comparison; source/physics strings unchanged. Successful original/base native initialization retained and independently audited after the metadata-only fix; only three missing initializations executed thereafter.'),
 copied_input_hashes={row['original']:row['sha256'] for row in copied},curation_recipe_sha256=sha(__file__))
save('preflight-summary.json',summary)
save('preflight-runtime-summary.json',dict(schema=1,full_manifests_retained_on_Linux=True,
 manifests=runtime_summaries,build_pin_comparison=runner['build_pin_comparison']))
save('preflight-curation-index.json',dict(schema=1,files=copied,curation_recipe_sha256=sha(__file__),
 excluded='Native coefficient fields/histories and full runtime/source manifests remain in the source-bound Linux run directories.'))
subprocess.run(['python3','-B',str(fingerprint),'--root',str(HIST),'--output',str(OUT/'runtime-after.json'),
    '--compare',str(OUT/'runtime-before.json')],check=True)
after=json.loads((OUT/'runtime-after.json').read_text());assert after['success'] and after['identical_qualification_inputs_and_runtime']
assert all(sha(row['original'])==row['sha256'] for row in copied)
assert all(sha(path)==digest for path,digest in oracle_pins.items())
save('preflight-postprocess-receipt.json',dict(schema=1,status='passed',driver_sha256=sha(__file__),
 numerical_threads={key:os.environ[key] for key in runner['numerical_threads']},affinity=[49],
 actual_preflight_runtime_identical_before_postprocess=True,postprocess_runtime_identical_before_after=True,
 before=dict(path=str(OUT/'runtime-before.json'),sha256=sha(OUT/'runtime-before.json')),
 after=dict(path=str(OUT/'runtime-after.json'),sha256=sha(OUT/'runtime-after.json')),
 source_inputs_unchanged=True,solver_invoked=False,summary_sha256=sha(CURATED/'preflight-summary.json')))
print(json.dumps(dict(status='passed',copied_files=len(copied),copied_bytes=sum(row['bytes'] for row in copied),
 output=str(CURATED),quadrature_comparisons=comparisons,archived_initial_correspondence=original_match,
 outer_gradient_witnesses=witnesses)),flush=True)
