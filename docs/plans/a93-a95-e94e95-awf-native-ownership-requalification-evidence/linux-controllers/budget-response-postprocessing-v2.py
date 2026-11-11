import pathlib,sys,os,time,subprocess,json
from concurrent.futures import ThreadPoolExecutor
r=pathlib.Path(__file__).parent;s=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(s))
from qualify_rans_preservation import load
from rans_qualification import THREADS,write_json
from rans_worker_budget import spawn
from rans_duration_contract import MEMBERS
os.environ.update(THREADS,LD_LIBRARY_PATH=str(r/'runtime/lib'),TMPDIR=str(r/'tmp'),NEKTAR_RANS_PROBE_TIMEOUT='14400')
out=r/'physical-postprocessing';out.mkdir(exist_ok=True)
neutral=r/'corrective-duration-response-v2/A95-N-Accepted-neutral.json'
while not neutral.exists():time.sleep(5)
q=load(neutral)
if q['status']!='pass':
 write_json(out/'results.json',dict(status='not_admitted',reason='corrected A95 neutral .5 did not pass'));sys.exit(1)
dt=q['dt'];cohort='A95-N-Accepted-response-2e-05-dt-'+str(dt)
full=r/'response-cohorts'/cohort/'T-0.016/results.json'
while not full.exists():
 controller=r/'corrective-duration-response-v2/response-results.json'
 if controller.exists() and load(controller)['status']=='not_admitted':sys.exit(1)
 time.sleep(5)
events=[]
def run(name,cmd,native=True):
 start=time.monotonic()
 with (out/(name+'.log')).open('x') as f:
  if native:
   proc=None
   while proc is None:
    proc=spawn(r,['python3','-B',*map(str,cmd)],env=dict(os.environ),stdout=f,stderr=subprocess.STDOUT)
    if proc is None:time.sleep(3)
   pid,status,usage=os.wait4(proc.pid,0)
   code=os.waitstatus_to_exitcode(status);proc.returncode=code
   cpu=dict(user=usage.ru_utime,system=usage.ru_stime,total=usage.ru_utime+usage.ru_stime,scope='wait4 for this worker, including its reaped descendants; independent of the other parallel worker')
  else:
   import resource
   before=resource.getrusage(resource.RUSAGE_CHILDREN)
   code=subprocess.run(['python3','-B',*map(str,cmd)],env=dict(os.environ),stdout=f,stderr=subprocess.STDOUT).returncode
   after=resource.getrusage(resource.RUSAGE_CHILDREN)
   cpu=dict(user=after.ru_utime-before.ru_utime,system=after.ru_stime-before.ru_stime,total=after.ru_utime-before.ru_utime+after.ru_stime-before.ru_stime,scope='sequential postprocessing child and reaped descendants')
 return dict(name=name,command=list(map(str,cmd)),exit_code=code,process_cpu_seconds=cpu,wall_seconds=time.monotonic()-start)
jobs=[]
for kind in ('perturbed','baseline'):
 for member in MEMBERS:
  for qi,face in (((1,16),(2,16),(2,32)) if member in MEMBERS[:2] else ((2,32),)):
   label=(cohort if kind=='perturbed' else 'A95-N-Accepted-duration-'+str(dt))+'-'+member
   name=cohort+'-'+kind+'-'+member+'-0.016-Q'+str(qi)+('-F32' if face==32 else '')
   jobs.append((name,[s/'audit_rans_physical_scalar_budget.py','--root',r,'--label',label,'--time','.016','--output-label',name,'--quadrature-index',qi,'--face-order',face]))
with ThreadPoolExecutor(max_workers=2) as pool:events=list(pool.map(lambda job:run(*job),jobs))
write_json(out/'execution.json',dict(checks=events))
if any(x['exit_code'] for x in events):
 write_json(out/'results.json',dict(status='reconstruction_failed',execution=events));sys.exit(1)
analysis=r/'physical-budget-analyses'/(cohort+'-0.016');analysis.parent.mkdir(exist_ok=True)
accuracy=r/'physical-budget-accuracy'/(cohort+'-0.016');accuracy.parent.mkdir(exist_ok=True)
binding=r/'physical-budget-bindings'/(cohort+'-0.016');binding.parent.mkdir(exist_ok=True)
events.append(run('analysis',[s/'analyze_rans_physical_budgets.py','--root',r,'--cohort',cohort,'--time','.016','--output',analysis],False))
events.append(run('accuracy',[s/'classify_rans_physical_budget_accuracy.py','--source',analysis/'results.json','--output',accuracy],False))
events.append(run('binding',[s/'bind_rans_ownership_budget_accuracy.py','--root',r,'--admission',q['admission'],'--cohort',cohort,'--dt',dt,'--time','.016','--analysis',analysis/'results.json','--accuracy',accuracy/'results.json','--output',binding],False))
# Independent decomposition retains raw dimensional terms and physical faces.
from rans_preservation_oracle import Reference
ref=Reference('A95','Natural','Accepted');decomposition={}
for kind in ('perturbed','baseline'):
 for member in MEMBERS:
  leaf=r/'physical-budgets'/(cohort+'-'+kind+'-'+member+'-0.016-Q2-F32');receipt=load(leaf/'results.json')
  b=receipt['physical_budget'];face_groups={}
  for face in b['physical_faces']:
   key='region-'+str(face['region'])+'/axis-'+str(face['axis'])+'/side-'+str(face['outward_side'])
   group=face_groups.setdefault(key,dict(k_frozen_flux=0.,k_current_flux=0.,epsilon_frozen_flux=0.,epsilon_current_flux=0.,count=0))
   group['count']+=1
   for name in ('k_frozen_flux','k_current_flux','epsilon_frozen_flux','epsilon_current_flux'):
    if face[name] is not None:group[name]+=face['weight']*face[name]
  decomposition[kind+'-'+member]=dict(configuration=receipt['configuration'],dt=receipt['dt'],equations=b['equations'],physical_face_integrals=face_groups,scaled_terms=receipt['scaled'],source_report_sha256=receipt['source_report_sha256'],state_sha256=receipt['state_sha256'],previous_state_sha256=receipt['previous_state_sha256'])
dest=r/'physical-budget-diagnosis'/cohort;dest.mkdir(parents=True)
write_json(dest/'results.json',dict(status='diagnosis_completed',time=.016,members=decomposition,
 scope='signed storage/source/physical face flux contributions; interface retained by region but excluded from whole-patch external total in the separately sealed analysis; no prescribed-row reaction substituted'))
write_json(out/'results.json',dict(status='pass' if events[-1]['exit_code']==0 else 'accuracy_failed_or_unresolved',execution=events,analysis=str(analysis/'results.json'),accuracy=str(accuracy/'results.json'),binding=str(binding/'results.json')))
responses=r/'corrective-duration-response-v2/response-results.json'
while not responses.exists():time.sleep(5)
half=r/'response-cohorts'/('A95-N-Accepted-response-1e-05-dt-'+str(dt))/'T-0.016'
if half.exists():
 events.append(run('half-amplitude',[s/'compare_rans_half_response.py','--full',full.parent,'--half',half,'--output',out/'half-amplitude.json'],False))
ep=load(responses).get('EP_neutral',{})
if ep.get('status')=='pass' and ep['dt']==dt:
 events.append(run('wall-comparison',[s/'compare_rans_wall_responses.py','--root',r,'--model','A95','--accepted-prefix',cohort,'--ep-prefix','A95-N-EP-response-2e-05-dt-'+str(dt),'--dt',dt,'--time','.016','--output',r/'common-wall-response'],True))
write_json(out/'complete-execution.json',dict(checks=events))
