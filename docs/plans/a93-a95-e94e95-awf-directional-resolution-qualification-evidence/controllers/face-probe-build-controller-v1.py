import os,sys,pathlib,subprocess,resource,time,json,hashlib,shutil
r=pathlib.Path(__file__).parent
s=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(s))
from rans_request_manifest import runtime_identity
from rans_qualification import THREADS,write_json
os.environ.update(THREADS,TMPDIR=str(r/'tmp'),PYTHONDONTWRITEBYTECODE='1')
start=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN)
with (r/'logs/face-probe-build-v1.log').open('x') as f:
 code=subprocess.run(['cmake','--build',str(r/'build'),'--parallel','1','--target','RANSPreservationProbe'],stdout=f,stderr=subprocess.STDOUT).returncode
after=resource.getrusage(resource.RUSAGE_CHILDREN)
result=dict(exit_code=code,wall_seconds=time.monotonic()-start,process_cpu_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,source_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=r/'source',text=True).strip(),scope='probe-only additional face-component output; frozen scientific solver and probe unchanged')
if code==0:
 out=r/'runtime-face-decomposition-v1';out.mkdir();(out/'lib').symlink_to(r/'runtime/lib');shutil.copy2(r/'build/bin/RANSPreservationProbe',out/'RANSPreservationProbe')
 result['runtime']=runtime_identity(out/'RANSPreservationProbe')
 result['scientific_runtime_unchanged']=json.loads((r/'runtime-identity.json').read_text())
write_json(r/'face-probe-build-execution.json',result)
raise SystemExit(code)
