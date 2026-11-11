import os,sys,json,time,math,subprocess,pathlib,hashlib
from pathlib import Path
r=Path(__file__).parent;s=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(s))
from rans_qualification import THREADS,write_json,exact_vector_hash
os.environ.update(THREADS,LD_LIBRARY_PATH=str(r/'runtime/lib'),TMPDIR=str(r/'tmp'),NEKTAR_RANS_PROBE_TIMEOUT='14400')
from qualify_rans_preservation import load,replay_request,nested
from rans_request_manifest import RequestManifest,sha,runtime_identity
from rans_duration_contract import require,identity
from rans_worker_budget import spawn
from rans_preservation_oracle import Reference,physical_grid,weak_balance
from analyze_rans_physical_budgets import whole_patch
ref=Reference('A95','Natural','Accepted');cfg=dict(model='A95',policy='Natural',wall='Accepted')
cohort='A95-N-Accepted-response-2e-05-dt-0.002'
binding=r/'physical-budget-bindings'/(cohort+'-0.016')/'results.json'
b=load(binding);require(b['accuracy_status']=='failed_or_unresolved' and b['reconstruction_status']=='pass','no_failed_accuracy_to_diagnose')
analysis_path=r/'physical-budget-analyses'/(cohort+'-0.016')/'results.json'
require(sha(analysis_path)==b['binding']['analysis_sha256'],'main_analysis_changed')
main_analysis=load(analysis_path)
main_inputs=load(analysis_path.parent/'request-manifest.json')['request']['inputs']
def run(out,name,cmd):
 start=time.monotonic();proc=None
 with (out/(name+'.log')).open('x') as stream:
  while proc is None:
   proc=spawn(r,['python3','-B',*map(str,cmd)],stdout=stream,stderr=subprocess.STDOUT,env=dict(os.environ))
   if proc is None:time.sleep(3)
  pid,status,usage=os.wait4(proc.pid,0);code=os.waitstatus_to_exitcode(status);proc.returncode=code
 return dict(name=name,command=list(map(str,cmd)),exit_code=code,wall_seconds=time.monotonic()-start,
  process_cpu_seconds=dict(user=usage.ru_utime,system=usage.ru_stime,total=usage.ru_utime+usage.ru_stime))
def budget_terms(d):return dict(d['scaled'],whole_patch_k=whole_patch(d,ref))
def faces(d):
 groups={}
 for f in d['physical_budget']['physical_faces']:
  key='region-'+str(f['region'])+'/axis-'+str(f['axis'])+'/side-'+str(f['outward_side']);v=groups.setdefault(key,{})
  for name in ('k_frozen_flux','k_current_flux','epsilon_frozen_flux','epsilon_current_flux'):
   if f[name] is not None:v[name]=v.get(name,0.)+f['weight']*f[name]
 return groups
results=[]
for control in ('NormalEndpointsV1','StreamwiseP6V1'):
 parent=r/'budget-axis-controls'/control;out=parent/'execution';h=load(parent/'hypothesis.json')
 request=dict(schema='PhysicalBudgetAxisExecutionV1',configuration=cfg,control=control,hypothesis_sha256=sha(parent/'hypothesis.json'),mesh_sha256=sha(parent/'mesh.json'),failed_binding_sha256=sha(binding),driver_sha256=sha(s/'run_rans_budget_axis_control.py'),controller_sha256=sha(Path(__file__)),runtime=runtime_identity(r/'runtime/IncNavierStokesSolver'),time=.016,timesteps=[.002,.001],amplitudes=[0.,2e-5],scope='single predeclared directional control; original required accuracy decision retained')
 seal=RequestManifest(out,request,False);events=[];labels={};first_failure=None
 def execute(name,cmd):
  value=run(out,name,cmd);events.append(value);write_json(out/'progress.json',dict(status='running',events=events));return value['exit_code']
 try:
  for kind,amp in [('baseline',0.),('perturbed',2e-5)]:
   for timing,dt in [('primary',.002),('half',.001)]:
    label='A95-N-Accepted-budget-'+control+'-'+kind+'-'+timing;labels[kind,timing]=label
    cmd=[s/'run_rans_budget_axis_control.py','--root',r,'--control',control,'--hypothesis',parent/'hypothesis.json','--binding',binding,'--mesh-file',parent/'mesh.json','--label',label,'--dt',dt,'--amplitude',amp]
    require(execute('initialize-'+kind+'-'+timing,[*cmd,'--initialize-only'])==0,'initialization_screen_failed:'+kind+'/'+timing)
  initial={}
  for kind in ('baseline','perturbed'):
   q=[load(r/'campaigns'/labels[kind,t]/'request-manifest.json')['request'] for t in ('primary','half')]
   for key in ('raw_neutral_sha256','reference_sha256','boundary','declared_mesh','spatial_control'):
    require(q[0][key]==q[1][key],'same_space_initial_input_mismatch:'+key)
   if kind=='perturbed':require(q[0]['actual_initial']['state_sha256']==q[1]['actual_initial']['state_sha256'],'same_space_actual_initial_changed')
   initial[kind]=[dict(request_sha256=sha(r/'campaigns'/labels[kind,t]/'request-manifest.json'),actual_initial=x.get('actual_initial'),raw_neutral_sha256=x['raw_neutral_sha256'],reference_sha256=x['reference_sha256'],boundary=x['boundary']) for t,x in zip(('primary','half'),q,strict=True)]
  write_json(out/'initial-identities.json',dict(status='pass',cases=initial))
  for kind,amp in [('baseline',0.),('perturbed',2e-5)]:
   for timing,dt in [('primary',.002),('half',.001)]:
    label=labels[kind,timing];cmd=[s/'run_rans_budget_axis_control.py','--root',r,'--control',control,'--hypothesis',parent/'hypothesis.json','--binding',binding,'--mesh-file',parent/'mesh.json','--label',label,'--dt',dt,'--amplitude',amp,'--resume']
    code=execute('evolve-'+kind+'-'+timing,cmd)
    if code and first_failure is None:first_failure=dict(phase='native_evolution',label=label,exit_code=code)
  require(first_failure is None,'native_control_failed')
  receipts={};leaf_paths={}
  for kind in ('baseline','perturbed'):
   for timing in ('primary','half'):
    for qi,face in (((1,16),(2,16),(2,32)) if timing=='primary' else ((2,32),)):
     name='axis-'+control+'-'+kind+'-'+timing+'-Q'+str(qi)+'-F'+str(face);label=labels[kind,timing]
     require(execute('budget-'+name,[s/'audit_rans_physical_scalar_budget.py','--root',r,'--label',label,'--time','.016','--output-label',name,'--quadrature-index',qi,'--face-order',face])==0,'physical_budget_reconstruction_failed:'+name)
     path=r/'physical-budgets'/name;req=load(path/'request-manifest.json')['request'];require(RequestManifest(path,req,True).completed(path),'unsealed_axis_budget')
     d=load(path/'results.json');identity(d['configuration'],cfg);require(abs(d['time']-.016)<1e-13,'wrong_control_time')
     receipts[kind,timing,qi,face]=d;leaf_paths[name]=dict(path=str(path),sha256=sha(path/'results.json'))
  comparisons=[];quadrature=[];decomp={}
  for kind in ('baseline','perturbed'):
   main={}
   for timing,member in [('primary','B5'),('half','B5-dt2')]:
    name=cohort+'-'+kind+'-'+member+'-0.016-Q2-F32';path=r/'physical-budgets'/name
    require(sha(path/'results.json')==main_inputs[name]['result_sha256'],'main_budget_changed')
    mr=load(path/'request-manifest.json')['request'];require(RequestManifest(path,mr,True).completed(path),'main_budget_unsealed')
    main[timing]=load(path/'results.json');leaf_paths[name]=dict(path=str(path),sha256=sha(path/'results.json'))
   for timing in ('primary','half'):
    d=receipts[kind,timing,2,32];decomp[kind+'-'+timing]=dict(scaled=budget_terms(d),equations=d['physical_budget']['equations'],physical_face_integrals=faces(d),native_report_sha256=d['source_report_sha256'],state_sha256=d['state_sha256'],previous_state_sha256=d['previous_state_sha256'])
   versions=[receipts[kind,'primary',q,f] for q,f in ((1,16),(2,16),(2,32))]
   require(len({x['state_sha256'] for x in versions})==1 and len({x['previous_state_sha256'] for x in versions})==1,'control_quadrature_state_changed')
   for eq in ('inner_k','outer_k','outer_epsilon','whole_patch_k'):
    terms=[budget_terms(x)[eq] for x in versions];key='frozen_stage_physical_defect'
    changes={k:max(abs(terms[0][k]-terms[1][k]),abs(terms[1][k]-terms[2][k])) for k in terms[2]}
    for term,value in changes.items():
     limit=max(1e-9,.01*abs(terms[2][term]-budget_terms(main['primary'])[eq][term]));quadrature.append(dict(kind=kind,equation=eq,term=term,change=value,limit=limit,status='pass' if value<=limit else 'unresolved'))
    c=budget_terms(main['primary'])[eq][key];f=terms[2][key];ct=budget_terms(main['half'])[eq][key];ft=budget_terms(receipts[kind,'half',2,32])[eq][key]
    original_uncertainty=next(x['quadrature_uncertainty'] for x in main_analysis['finite_defects'] if x['kind']==kind and x['equation']==eq)
    uncertainty=max(changes[key],original_uncertainty);floor=max(5e-9,5*uncertainty);delta=max(abs(c-ct),abs(f-ft));limit=max(1e-9,.1*abs(c-f))
    spatial='unresolved_floor' if max(abs(c),abs(f))<=floor else 'pass' if abs(f)<=.85*abs(c) and abs(c)-abs(f)>5*uncertainty else 'failed_or_unresolved'
    comparisons.append(dict(kind=kind,equation=eq,coarse_B5=c,control=f,coarse_half_dt=ct,control_half_dt=ft,absolute_ratio=abs(f/c) if c else None,spatial_status=spatial,floor=floor,quadrature_uncertainty=uncertainty,timestep_difference=delta,temporal_limit=limit,temporal_status='pass' if delta<=limit else 'unresolved'))
  native={}
  for (kind,timing),label in labels.items():
   folder=r/'campaigns'/label;metrics=load(folder/'T-0.0160/metrics.json');rq=load(folder/'request-manifest.json')['request']
   native[kind+'-'+timing]=dict(label=label,initial=load(folder/'initial-metrics.json'),metrics=metrics,request_sha256=sha(folder/'request-manifest.json'),request=rq,native_report_sha256=sha(folder/'T-0.0160/case.evolution.json'))
  result=dict(status='diagnosis_completed',numerical_accuracy_status='pass' if all(x['spatial_status'] in ('pass','unresolved_floor') and x['temporal_status']=='pass' for x in comparisons) and all(x['status']=='pass' for x in quadrature) else 'failed_or_unresolved',configuration=cfg,control=control,numerical_operator='GeneratedVolumeOwnerV1',time=.016,unknowns=h['expected_unknowns'],initialization=initial,comparisons=comparisons,quadrature=quadrature,decomposition=decomp,native=native,inputs=leaf_paths,execution=events,first_failure=None,original_required_accuracy_reclassified=False,scope='directional physical-integral diagnostic; no old history, mesh transplant, source correction, or extension')
 except Exception as error:
  result=dict(status='diagnosis_failed',control=control,first_failure=first_failure or dict(error=repr(error)),execution=events,labels={str(k):v for k,v in labels.items()},original_required_accuracy_reclassified=False)
 write_json(out/'results.json',result);seal.seal(out,['results.json',*(['initial-identities.json'] if (out/'initial-identities.json').exists() else [])]);write_json(parent/'results.json',dict(result,execution_request_sha256=seal.payload['sha256']));results.append(dict(control=control,status=result['status'],accuracy_status=result.get('numerical_accuracy_status'),receipt=str(out/'results.json')))
complete=all(x['status']=='diagnosis_completed' for x in results);accurate=all(x.get('accuracy_status')=='pass' for x in results)
write_json(r/'budget-axis-control-results.json',dict(status='completed' if complete else 'failed_or_incomplete',numerical_accuracy_status='pass' if accurate else 'failed_or_unresolved',controls=results,original_accuracy_decision_unchanged=True))
raise SystemExit(0 if complete and accurate else 1)
