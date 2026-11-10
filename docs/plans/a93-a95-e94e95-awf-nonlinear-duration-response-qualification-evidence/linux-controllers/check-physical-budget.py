from pathlib import Path
import os,subprocess,json,time,sys
r=Path(__file__).resolve().parent;t=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(t));from rans_worker_budget import spawn
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(r/'tmp'),NEKTAR_RANS_PROBE_TIMEOUT='14400');cmd=['python3','-B',str(t/'audit_rans_physical_scalar_budget.py'),'--root',str(r),'--label','A93-N-Accepted-duration-0.002-B5','--time','.016','--output-label','test-neutral-A93-B5-0.016-Q64']
with (r/'check-physical-budget-child.log').open('x') as f:
 p=None
 while p is None:
  p=spawn(r,cmd,env=env,stdout=f,stderr=subprocess.STDOUT)
  if p is None:time.sleep(5)
 code=p.wait()
(r/'physical-budget-regression.json').write_text(json.dumps(dict(status='pass' if code==0 else 'fail',exit_code=code,command=cmd,scope='optional physical budget reconstruction, not a conservation accuracy certificate'))+'\n')
