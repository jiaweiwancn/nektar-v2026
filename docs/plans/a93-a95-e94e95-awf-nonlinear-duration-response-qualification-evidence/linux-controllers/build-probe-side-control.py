from pathlib import Path
import subprocess,os,json,time
r=Path(__file__).resolve().parent;t=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(r/'tmp'));cmd=['cmake','--build',str(r/'build'),'--parallel','1','--target','RANSPreservationProbe'];start=time.monotonic()
with (r/'build-probe-side-control.log').open('x') as f:code=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
result=dict(status='pass' if code==0 else 'fail',command=cmd,exit_code=code,seconds=time.monotonic()-start)
if code==0:
 freeze=['python3','-B',str(t/'freeze_rans_diagnostic_runtime.py'),'--root',str(r),'--output-name','runtime-probe-v6','--parent-name','runtime-budget-v5','--build-workers','1']
 with (r/'freeze-probe-v6.log').open('x') as f:result['freeze_exit_code']=subprocess.run(freeze,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
 if result['freeze_exit_code']!=0:result['status']='fail'
with (r/'test-physical-accounting.log').open('x') as f:result['accounting_unit_exit_code']=subprocess.run(['python3','-B',str(t/'test_rans_physical_budget_analysis.py')],env=env,stdout=f,stderr=subprocess.STDOUT).returncode
(r/'build-probe-side-control.json').write_text(json.dumps(result,indent=2)+'\n')
