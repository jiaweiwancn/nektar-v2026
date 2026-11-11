import pathlib,sys,os,subprocess,time,json
from concurrent.futures import ThreadPoolExecutor
Path=pathlib.Path
r=Path(__file__).parent;s=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(s))
from qualify_rans_preservation import load
from rans_qualification import write_json,THREADS
from rans_duration_contract import MEMBERS,prefix,require,identity
from rans_request_manifest import sha,RequestManifest
os.environ.update(THREADS,LD_LIBRARY_PATH=str(r/'runtime/lib'),TMPDIR=str(r/'tmp'),NEKTAR_RANS_PROBE_TIMEOUT='14400')
out=r/'corrective-duration-response-v2';out.mkdir(exist_ok=True)
initial_process={'A95':606051,'A93':616396}
def call(name,cmd):
 with (out/(name+'.log')).open('x') as f:
  c=subprocess.run(['python3','-B',*map(str,cmd)],stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ))
 value=dict(exit_code=c.returncode,command=list(map(str,cmd)))
 write_json(out/(name+'.execution.json'),value)
 return c.returncode
def schedule_result(model,dt,pid=None,wall='Accepted',horizon=.5):
 tag=prefix(model,'Natural',wall)
 schedule=r/'duration-schedules'/(tag+'-duration-'+str(dt)+'-'+str(horizon))/'results.json'
 while True:
  d=load(schedule) if schedule.exists() else {}
  milestones=d.get('milestones',[])
  if milestones and milestones[-1]['status']!='pass':return schedule,d
  if milestones and abs(milestones[-1]['time']-horizon)<1e-13:return schedule,d
  if pid is not None:
   p=Path('/proc')/str(pid)/'cmdline'
   cmd=p.read_bytes().replace(b'\0',b' ').decode() if p.exists() else ''
   if str(r) not in cmd or 'run_rans_configuration_duration.py' not in cmd:
    return schedule,dict(d,status='fail',failure='existing exact duration controller ended before required final milestone')
  time.sleep(5)
def classify(model,dt,schedule,d,wall='Accepted',horizon=.5):
 tag=prefix(model,'Natural',wall);ad=r/'admissions'/(tag+'.json');a=load(ad)
 identity(d.get('configuration',{}),a['configuration'])
 expected=[t for t in a['milestones'] if t<=horizon+1e-13]
 m=d.get('milestones',[])
 passed=d['status']=='pass' and len(m)==len(expected) and all(abs(x['time']-t)<1e-13 and x['status']=='pass' for x,t in zip(m,expected,strict=True))
 return dict(model=model,wall=wall,status='pass' if passed else 'failed_or_unresolved',dt=dt,admission=str(ad),
  trials=[dict(dt=dt,exit_code=0 if passed else 1,schedule=str(schedule),schedule_sha256=sha(schedule) if schedule.exists() else None,body=d)],
  scheduling_policy='existing native histories unchanged; independently consumed complete configuration-specific duration receipts')
def neutral(model):
 dt=.002;path,d=schedule_result(model,dt,initial_process[model]);v=classify(model,dt,path,d)
 if v['status']!='pass':
  last=(d.get('milestones') or [{}])[-1];comparison=Path(last.get('comparison','/nonexistent'));c=load(comparison) if comparison.exists() else {}
  temporal=any(x['temporal_status']!='pass' for x in c.get('fields',{}).values())
  otherwise=bool(c.get('fields')) and all(x['spatial_status'] in ('pass','unresolved_floor') and x['triangle_status']=='pass' for x in c['fields'].values()) and all(x['status']=='pass' for x in c.get('conservation',[]))
  if temporal and otherwise:
   tag=prefix(model,'Natural','Accepted');ad=r/'admissions'/(tag+'.json');dt=.001
   call('duration-'+tag+'-'+str(dt),[s/'run_rans_configuration_duration.py','--root',r,'--admission',ad,'--dt',dt,'--max-time','.5','--concurrency','2'])
   path=r/'duration-schedules'/(tag+'-duration-'+str(dt)+'-0.5')/'results.json';d=load(path)
   refined=classify(model,dt,path,d);refined['trials']=v['trials']+refined['trials'];v=refined
 write_json(out/(prefix(model,'Natural','Accepted')+'-neutral.json'),v);return v
def response(a95):
 if a95['status']!='pass':
  q=dict(status='not_admitted',reason='corrected A95 complete .5 neutral contract did not pass');write_json(out/'response-results.json',q);return q
 dt=a95['dt'];proof=r/'duration-audits'/('A95-N-Accepted-duration-'+str(dt)+'-0.5')/'results.json'
 rows=[]
 for amp in (2e-5,1e-5):
  name='A95-N-Accepted-response-'+str(amp)+'-dt-'+str(dt)
  code=call('response-'+name,[s/'run_rans_nonlinear_response.py','--root',r,'--admission',a95['admission'],'--accepted-neutral-proof',proof,'--dt',dt,'--max-time','.016','--amplitude',amp,'--concurrency','2'])
  rows.append(dict(cohort=name,exit_code=code,receipt=str(r/'response-cohorts'/name/'results.json')))
 epad=r/'admissions/A95-N-EP.json'
 code=call('admit-A95-N-EP',[s/'admit_rans_ownership_duration.py','--root',r,'--model','A95','--wall','EP']) if not epad.exists() else 0
 ep=dict(model='A95',wall='EP',status='not_admitted',admission=str(epad))
 if code==0:
  trials=[]
  for epdt in (.002,.001):
   code=call('duration-A95-N-EP-'+str(epdt),[s/'run_rans_configuration_duration.py','--root',r,'--admission',epad,'--dt',epdt,'--max-time','.016','--concurrency','2'])
   path=r/'duration-schedules'/('A95-N-EP-duration-'+str(epdt)+'-0.016')/'results.json';d=load(path);trial=classify('A95',epdt,path,d,wall='EP',horizon=.016);trials+=trial['trials'];ep=dict(trial,trials=trials)
   if code==0:break
   last=(d.get('milestones') or [{}])[-1];c=load(Path(last['comparison'])) if last.get('comparison') else {}
   if not (any(x['temporal_status']!='pass' for x in c.get('fields',{}).values()) and all(x['spatial_status'] in ('pass','unresolved_floor') and x['triangle_status']=='pass' for x in c.get('fields',{}).values())):break
  write_json(out/'A95-N-EP-neutral.json',ep)
 if ep['status']=='pass':
  epdt=ep['dt'];name='A95-N-EP-response-2e-05-dt-'+str(epdt)
  code=call('response-'+name,[s/'run_rans_nonlinear_response.py','--root',r,'--admission',epad,'--accepted-neutral-proof',proof,'--dt',epdt,'--max-time','.016','--amplitude','2e-5','--concurrency','2'])
  rows.append(dict(cohort=name,exit_code=code,receipt=str(r/'response-cohorts'/name/'results.json')))
 result=dict(status='pass' if len(rows)==3 and all(x['exit_code']==0 for x in rows) else 'failed_unresolved_or_incomplete',checks=rows,EP_neutral=ep)
 write_json(out/'response-results.json',result);return result
with ThreadPoolExecutor(max_workers=2) as pool:
 a93=pool.submit(neutral,'A93');a95=pool.submit(neutral,'A95')
 answer95=a95.result();responses=response(answer95);answer93=a93.result()
write_json(out/'neutral-results.json',dict(configurations=[answer95,answer93]))
