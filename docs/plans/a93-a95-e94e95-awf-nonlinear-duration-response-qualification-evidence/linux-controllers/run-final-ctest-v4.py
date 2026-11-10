from pathlib import Path
import os,sys,subprocess,time,json,shutil,hashlib
r=Path(__file__).resolve().parent;t=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(t));from rans_worker_budget import spawn
from rans_request_manifest import sha
env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(r/'tmp'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',LD_LIBRARY_PATH=str(r/'runtime-probe-v7/lib'));checks=[]
commands=[['cmake','-S',str(t),'-B',str(r/'build')],['ctest','--test-dir',str(r/'build'),'-j','1','--output-on-failure','-R','^(RANSDurationCoverageContract|RANSNonlinearEvidenceContract|RANSPilotEvidenceContract|RANSQualificationFailureContract|RANSForwardInterfaceTolerance)$']]
for j,cmd in enumerate(commands):
 with (r/('ctest-final-v4-'+str(j)+'.log')).open('x') as f:
  p=None
  while p is None:
   p=spawn(r,cmd,env=env,stdout=f,stderr=subprocess.STDOUT)
   if p is None:time.sleep(5)
  rc=p.wait();checks.append(dict(command=cmd,exit_code=rc))
 if rc:break
(r/'ctest-final-v4.json').write_text(json.dumps(dict(status='pass' if len(checks)==2 and all(c['exit_code']==0 for c in checks) else 'fail',checks=checks),indent=2)+'\n')
folder=r/'loader-probe-v7';folder.mkdir();shutil.copy2(r/'forward-interface-exact-scope-regression/A95-Accepted/case.xml',folder/'case.xml');shutil.copy2(r/'forward-interface-exact-scope-regression/A95-Accepted/request.json',folder/'request.json');cmd=[str(r/'runtime-probe-v7/RANSPreservationProbe'),'case.xml','request.json','packet.json']
with (folder/'probe.log').open('x') as f:
 p=None
 while p is None:
  p=spawn(r,cmd,cwd=folder,env=env,stdout=f,stderr=subprocess.STDOUT)
  if p is None:time.sleep(5)
 libraries={}
 while p.poll() is None:
  maps=Path('/proc')/str(p.pid)/'maps'
  if maps.exists():
   for line in maps.read_text().splitlines():
    fields=line.split(maxsplit=5)
    if len(fields)==6 and fields[5].startswith('/'):
     v=Path(fields[5])
     if v.is_file() and '.so' in v.name:libraries[str(v)]=sha(v)
  time.sleep(.02)
 rc=p.returncode
 assert rc==0 and libraries
 (r/'runtime-probe-v7/actual-loaded-libraries-current.json').write_text(json.dumps(dict(status='pass',executable=str(r/'runtime-probe-v7/RANSPreservationProbe'),executable_sha256=sha(r/'runtime-probe-v7/RANSPreservationProbe'),libraries=libraries,command=cmd,exit_code=rc),indent=2)+'\n')
