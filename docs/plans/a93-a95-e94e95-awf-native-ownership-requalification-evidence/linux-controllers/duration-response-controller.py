import pathlib,sys,os,subprocess,time,json
from concurrent.futures import ThreadPoolExecutor
r=pathlib.Path(__file__).parent;s=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(s))
from qualify_rans_preservation import load
from rans_qualification import write_json,THREADS
from rans_duration_contract import MEMBERS,prefix,require
from rans_request_manifest import sha
os.environ.update(THREADS,LD_LIBRARY_PATH=str(r/'runtime/lib'),TMPDIR=str(r/'tmp'),NEKTAR_RANS_PROBE_TIMEOUT='14400')
out=r/'corrective-duration-response';out.mkdir()
while not (r/'ownership-short-qualification/results.json').exists():
 if (r/'short-controller-completion.json').exists():raise RuntimeError('short controller ended without complete receipt')
 time.sleep(5)
def call(name,cmd):
 with (out/(name+'.log')).open('x') as f:c=subprocess.run(['python3','-B',*map(str,cmd)],stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ))
 write_json(out/(name+'.execution.json'),dict(exit_code=c.returncode,command=list(map(str,cmd))))
 return c.returncode
def duration(model,wall='Accepted',limit=.5):
 tag=prefix(model,'Natural',wall);ad=r/'admissions'/(tag+'.json')
 if not ad.exists():
  code=call('admit-'+tag,[s/'admit_rans_ownership_duration.py','--root',r,'--model',model,'--wall',wall])
  if code:return dict(model=model,wall=wall,status='not_admitted',admission_log=str(out/('admit-'+tag+'.log')))
 trials=[]
 for dt in (.002,.001):
  code=call('duration-'+tag+'-'+str(dt),[s/'run_rans_configuration_duration.py','--root',r,'--admission',ad,'--dt',dt,'--max-time',limit,'--concurrency','2'])
  schedule=r/'duration-schedules'/(tag+'-duration-'+str(dt)+'-'+str(limit))/'results.json'
  data=load(schedule) if schedule.exists() else dict(status='fail',failure='controller ended without schedule')
  trials.append(dict(dt=dt,exit_code=code,schedule=str(schedule),body=data))
  if code==0:
   value=dict(model=model,wall=wall,status='pass',dt=dt,admission=str(ad),trials=trials)
   write_json(out/(tag+'-neutral.json'),value);return value
  # A new trial is permitted only for an actual timestep-uncertainty failure.
  last=(data.get('milestones') or [{}])[-1];comparison=Path(last.get('comparison','/nonexistent'))
  c=load(comparison) if comparison.exists() else {}
  temporal=any(v['temporal_status']!='pass' for v in c.get('fields',{}).values())
  otherwise=bool(c.get('fields')) and all(v['spatial_status'] in ('pass','unresolved_floor') and v['triangle_status']=='pass' for v in c['fields'].values()) and all(v['status']=='pass' for v in c.get('conservation',[]))
  if not (temporal and otherwise):break
 value=dict(model=model,wall=wall,status='failed_or_unresolved',trials=trials,admission=str(ad))
 write_json(out/(tag+'-neutral.json'),value);return value
with ThreadPoolExecutor(max_workers=2) as pool:neutral=list(pool.map(duration,('A95','A93')))
write_json(out/'neutral-results.json',dict(configurations=neutral))
a95=next(x for x in neutral if x['model']=='A95')
if a95['status']!='pass':
 write_json(out/'response-results.json',dict(status='not_admitted',reason='corrected A95 Accepted complete neutral contract did not pass'));sys.exit(1)
dt=a95['dt'];proof=r/'duration-audits'/('A95-N-Accepted-duration-'+str(dt)+'-0.5')/'results.json'
responses=[]
for amplitude in (2e-5,1e-5):
 name='A95-N-Accepted-response-'+str(amplitude)+'-dt-'+str(dt)
 code=call('response-'+name,[s/'run_rans_nonlinear_response.py','--root',r,'--admission',a95['admission'],'--accepted-neutral-proof',proof,'--dt',dt,'--max-time','.016','--amplitude',amplitude,'--concurrency','2'])
 responses.append(dict(cohort=name,exit_code=code,receipt=str(r/'response-cohorts'/name/'results.json')))
ep=duration('A95','EP',.016)
if ep['status']=='pass':
 epdt=ep['dt'];name='A95-N-EP-response-2e-05-dt-'+str(epdt)
 code=call('response-'+name,[s/'run_rans_nonlinear_response.py','--root',r,'--admission',ep['admission'],'--accepted-neutral-proof',proof,'--dt',epdt,'--max-time','.016','--amplitude','2e-5','--concurrency','2'])
 responses.append(dict(cohort=name,exit_code=code,receipt=str(r/'response-cohorts'/name/'results.json')))
write_json(out/'response-results.json',dict(status='pass' if len(responses)==3 and all(x['exit_code']==0 for x in responses) else 'failed_unresolved_or_incomplete',checks=responses,EP_neutral=ep))
