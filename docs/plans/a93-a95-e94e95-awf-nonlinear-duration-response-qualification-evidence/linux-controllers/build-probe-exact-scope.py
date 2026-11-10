from pathlib import Path
import subprocess,os,time,json,sys
r=Path(__file__).resolve().parent;t=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(t));from rans_worker_budget import spawn
cmd=['cmake','--build',str(r/'build'),'--parallel','1','--target','RANSPreservationProbe'];env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(r/'tmp'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1');start=time.monotonic()
with (r/'build-probe-exact-scope.log').open('x') as f:
 p=None
 while p is None:
  p=spawn(r,cmd,env=env,stdout=f,stderr=subprocess.STDOUT)
  if p is None:time.sleep(5)
 rc=p.wait()
 d=dict(status='pass' if rc==0 else 'fail',command=cmd,exit_code=rc,seconds=time.monotonic()-start)
 if rc==0:
  cmd=['python3','-B',str(t/'freeze_rans_diagnostic_runtime.py'),'--root',str(r),'--output-name','runtime-probe-v7','--parent-name','runtime-probe-v6','--build-workers','1']
  p=subprocess.run(cmd,env=env,capture_output=True,text=True);d.update(freeze_exit_code=p.returncode,freeze_output=p.stdout,freeze_error=p.stderr)
 (r/'build-probe-exact-scope-control.json').write_text(json.dumps(d,indent=2)+'\n')
 if d.get('freeze_exit_code')==0:subprocess.Popen(['python3','-B',str(r/'run-side-exact-scope-regression.py')],env=env,start_new_session=True,stdout=(r/'run-side-exact-scope-coordinator.log').open('x'),stderr=subprocess.STDOUT)
