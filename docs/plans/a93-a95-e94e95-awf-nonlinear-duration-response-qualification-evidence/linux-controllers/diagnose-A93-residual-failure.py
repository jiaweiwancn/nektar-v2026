from pathlib import Path
import os,sys,subprocess,time,json
r=Path(__file__).resolve().parent;t=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(t));from rans_worker_budget import spawn
while not (r/'evolved-audits/A93-N-Accepted-duration-0.002-0.5-B11-dt2/results.json').exists():time.sleep(10)
cmd=['python3','-B',str(t/'diagnose_rans_duration_residual.py'),'--root',str(r),'--admission',str(r/'admissions/A93-N-Accepted-loader-v2.json'),'--cohort','A93-N-Accepted-duration-0.002','--time','0.5','--output',str(r/'residual-diagnosis/A93-N-Accepted-0.5')]
env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(r/'tmp'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
with (r/'A93-residual-diagnosis.log').open('x') as f:
 p=None
 while p is None:
  p=spawn(r,cmd,env=env,stdout=f,stderr=subprocess.STDOUT)
  if p is None:time.sleep(5)
 code=p.wait()
(r/'A93-residual-diagnosis-process.json').write_text(json.dumps(dict(exit_code=code,command=cmd),indent=2)+'\n')
