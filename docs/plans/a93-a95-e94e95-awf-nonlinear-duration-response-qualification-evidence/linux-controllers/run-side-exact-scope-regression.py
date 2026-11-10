from pathlib import Path
import os,sys,subprocess,json,time
r=Path(__file__).resolve().parent;t=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(t));from rans_worker_budget import spawn
from rans_request_manifest import sha
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(r/'tmp'),NEKTAR_RANS_PROBE_TIMEOUT='14400',LD_LIBRARY_PATH=str(r/'runtime-probe-v7/lib'));cmd=['python3','-B',str(t/'check_rans_forward_interface_tolerance.py'),'--probe',str(r/'runtime-probe-v7/RANSPreservationProbe'),'--legacy-probe',str(r/'runtime-resolution-v3/RANSPreservationProbe'),'--descriptor-probe',str(r/'runtime/RANSNewModelsProbe'),'--output',str(r/'forward-interface-exact-scope-regression')]
with (r/'forward-interface-exact-scope-process.log').open('x') as f:
 p=None
 while p is None:
  p=spawn(r,cmd,env=env,stdout=f,stderr=subprocess.STDOUT)
  if p is None:time.sleep(5)
 code=p.wait()
(r/'forward-interface-exact-scope-status.json').write_text(json.dumps(dict(status='pass' if code==0 else 'fail',exit_code=code,command=cmd,script_sha256=sha(t/'check_rans_forward_interface_tolerance.py')),indent=2)+'\n')
