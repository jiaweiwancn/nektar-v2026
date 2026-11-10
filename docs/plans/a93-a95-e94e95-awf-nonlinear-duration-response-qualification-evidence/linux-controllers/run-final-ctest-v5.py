from pathlib import Path
import os,sys,subprocess,time,json
r=Path(__file__).resolve().parent;t=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(t));from rans_worker_budget import spawn
from rans_request_manifest import sha
cmd=['ctest','--test-dir',str(r/'build'),'-j','1','--output-on-failure','-R','^RANSForwardInterfaceTolerance$'];env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(r/'tmp'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',LD_LIBRARY_PATH=str(r/'runtime-probe-v7/lib'))
with (r/'ctest-final-v5.log').open('x') as f:
 p=None
 while p is None:
  p=spawn(r,cmd,env=env,stdout=f,stderr=subprocess.STDOUT)
  if p is None:time.sleep(5)
 rc=p.wait()
(r/'ctest-final-v5.json').write_text(json.dumps(dict(status='pass' if rc==0 else 'fail',command=cmd,exit_code=rc,earlier_four_passing_checks=dict(source=str(r/'ctest-final-v4.json'),sha256=sha(r/'ctest-final-v4.json')),preserved_missing_target_failure=str(r/'ctest-final-v4-1.log'),descriptor_sha256=sha(r/'build/RANSNewModelsProbe'),probe_sha256=sha(r/'build/RANSPreservationProbe'),registration_sha256=sha(r/'build/CTestTestfile.cmake'),source_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=r/'source',text=True).strip()),indent=2)+'\n')
