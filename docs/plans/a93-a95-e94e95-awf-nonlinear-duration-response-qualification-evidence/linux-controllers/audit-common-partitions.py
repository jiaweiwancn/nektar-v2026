from pathlib import Path
import os,subprocess,json
r=Path(__file__).resolve().parent;t=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1');checks=[]
for model,times in [('A93',[.016,.1]),('A95',[.016])]:
 for T in times:
  tag=model+'-N-Accepted-duration-0.002-'+str(T);out=r/'common-scalar-partitions'/tag;out.mkdir(parents=True)
  cmd=['python3','-B',str(t/'audit_rans_common_scalar_partitions.py'),'--admission',str(r/'admissions'/(model+'-N-Accepted-loader-v2.json')),'--audits',*[str(r/'evolved-audits'/(tag+'-'+m)/'results.json') for m in ('B5','B11')],'--time',str(T),'--output',str(out/'results.json')]
  with (out/'process.log').open('x') as f:code=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
  checks.append(dict(model=model,time=T,exit_code=code))
(r/'common-partition-regression.json').write_text(json.dumps(dict(status='pass' if all(x['exit_code']==0 for x in checks) else 'fail',checks=checks),indent=2)+'\n')
