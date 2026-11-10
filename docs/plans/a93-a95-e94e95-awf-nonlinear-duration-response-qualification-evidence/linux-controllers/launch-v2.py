from pathlib import Path
import os,subprocess,json
r=Path(__file__).resolve().parent;t=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';old=r.parent/'a93-a95-e94e95-balanced-20261110';env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(r/'tmp'));rows=[]
cmd=['python3','-B',str(t/'check_rans_duration_failures.py'),'--output',str(r/'aggregate-negative-controls')]
with (r/'aggregate-negative-controls.log').open('x') as f:code=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
assert code==0
for model in ('A93','A95'):
 tag=model+'-N-Accepted';ad=r/'admissions'/(tag+'-loader-v2.json');cmd=['python3','-B',str(t/'admit_rans_balanced_duration.py'),'--root',str(r),'--previous-root',str(old),'--model',model,'--wall','Accepted','--output',str(ad)]
 with (r/('admit-v4-'+tag+'.log')).open('x') as f:code=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
 print(tag,'admission',code,flush=True);rows.append(dict(model=model,admission=str(ad),admission_exit=code))
 if code:continue
 cmd=['python3','-B',str(t/'run_rans_balanced_duration.py'),'--root',str(r),'--admission',str(ad),'--dt','.002','--max-time','.5','--concurrency','2']
 with (r/(tag+'-duration-v2.log')).open('x') as f:p=subprocess.Popen(cmd,env=env,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
 rows[-1].update(pid=p.pid,command=cmd);print('COHORT',model,p.pid,flush=True)
 (r/'launch-v2-results.json').write_text(json.dumps(dict(status='launched' if all(x['admission_exit']==0 for x in rows) else 'preparation_failure',checks=rows),indent=2)+'\n')
