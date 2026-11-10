from pathlib import Path
import os,subprocess,json,time
r=Path(__file__).resolve().parent;old=r.parent/'a93-a95-e94e95-balanced-20261110';t=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(r/'tmp'));rows=[]
for model in ('A93','A95'):
 tag=model+'-N-Accepted'
 for m in ('B5-dt2','B11-dt2'):
  label=tag+'-'+m
  cmd=['python3','-B',str(t/'audit_rans_pilot_samples.py'),'--root',str(old),'--label',label,'--through-time','.016','--output-directory',str(r/'reused-short-observers'/label),'--tmp-directory',str(r/'tmp')]
  with (r/('repair-'+label+'.log')).open('x') as f:code=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
  rows.append(dict(label=label,exit_code=code,command=cmd));print(label,code,flush=True)
 cmd=['python3','-B',str(t/'admit_rans_balanced_duration.py'),'--root',str(r),'--previous-root',str(old),'--model',model,'--wall','Accepted']
 with (r/('admit-v2-'+tag+'.log')).open('x') as f:code=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
 rows.append(dict(label=tag+'-admission',exit_code=code,command=cmd));print(tag+'-admission',code,flush=True)
 (r/'short-coverage-repair-results.json').write_text(json.dumps(dict(status='pass' if all(x['exit_code']==0 for x in rows) else 'fail',checks=rows),indent=2)+'\n')
