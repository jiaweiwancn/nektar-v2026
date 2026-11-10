from pathlib import Path
import os,sys,subprocess,json,time
r=Path(__file__).resolve().parent;old=r.parent/'a93-a95-e94e95-balanced-20261110';t=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(t));from rans_worker_budget import spawn
from rans_request_manifest import sha
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(r/'tmp'),NEKTAR_RANS_PROBE_TIMEOUT='14400');events=[]
commands=[]
for m in ('B5','B11'):
 label='A95-N-EP-'+m;commands.append((label,['python3','-B',str(t/'audit_rans_pilot_samples.py'),'--root',str(old),'--label',label,'--through-time','.016','--output-directory',str(r/'reused-short-observers'/label),'--tmp-directory',str(r/'tmp')]))
commands.append(('A95-EP-admission',['python3','-B',str(t/'admit_rans_balanced_duration.py'),'--root',str(r),'--previous-root',str(old),'--model','A95','--wall','EP','--scalar-observer-root',str(r/'reused-short-observers'),'--output',str(r/'admissions/A95-N-EP-loader-v3.json')]))
for name,cmd in commands:
 with (r/(name+'-preparation-v2.log')).open('x') as f:
  p=None
  while p is None:
   p=spawn(r,cmd,env=env,stdout=f,stderr=subprocess.STDOUT)
   if p is None:time.sleep(5)
  code=p.wait()
 events.append(dict(name=name,command=cmd,exit_code=code));(r/'A95-EP-preparation-v2.json').write_text(json.dumps(dict(status='running' if code==0 and len(events)<len(commands) else 'pass' if code==0 else 'fail',events=events,scope='read-only short-evidence supplementation; no new EP trajectory admitted yet'),indent=2)+'\n')
 if code:raise SystemExit(code)
