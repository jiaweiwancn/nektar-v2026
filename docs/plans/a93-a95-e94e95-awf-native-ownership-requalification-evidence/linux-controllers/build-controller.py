import pathlib,sys,subprocess,os,time,json
r=pathlib.Path(__file__).parent;s=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(s))
from rans_worker_budget import spawn
from rans_qualification import THREADS,write_json
os.environ.update(THREADS,LD_LIBRARY_PATH=str(r/'dependency'),TMPDIR=str(r/'tmp'))
cmd=['cmake','--build',str(r/'build'),'--parallel','4','--target','IncNavierStokesSolver','RANSPreservationProbe','RANSNewModelsProbe','RANSNativeDiagnosticChecks','RANSWallResolutionProbe']
with (r/'build.log').open('x') as f:
 proc=None
 while proc is None:
  proc=spawn(r,cmd,env=dict(os.environ),stdout=f,stderr=subprocess.STDOUT)
  if proc is None:time.sleep(3)
 code=proc.wait()
write_json(r/'build-completion.json',dict(exit_code=code,command=cmd))
