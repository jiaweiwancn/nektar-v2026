import pathlib,subprocess,sys,os,time,json,shutil
r=pathlib.Path(__file__).parent;s=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(s))
from rans_worker_budget import spawn
from rans_request_manifest import runtime_identity,sha
from rans_qualification import THREADS,write_json
os.environ.update(THREADS,TMPDIR=str(r/'tmp'),LD_LIBRARY_PATH=str(r/'dependency'))
while not (r/'build-completion.json').exists():time.sleep(5)
assert json.loads((r/'build-completion.json').read_text())['exit_code']==0
with (r/'freeze-runtime.log').open('x') as f:subprocess.run(['python3','-B',str(s/'freeze_rans_pilot_runtime.py'),'--source',str(r/'source'),'--build',str(r/'build'),'--output',str(r/'runtime')],stdout=f,stderr=subprocess.STDOUT,check=True)
shutil.copy2(r/'build/RANSWallResolutionProbe',r/'runtime/RANSWallResolutionProbe')
manifest=json.loads((r/'runtime/source-build.json').read_text())
manifest.update(numerical_operator='GeneratedVolumeOwnerV1',build_started_cpp_revision='7cb6757',compiled_cpp_revision='7cb6757')
write_json(r/'runtime/source-build.json',manifest)
for name in ('runtime-resolution-v2','runtime-resolution-v3','runtime-budget-v5'):
 dest=r/name;shutil.copytree(r/'runtime',dest,symlinks=True)
 os.environ['LD_LIBRARY_PATH']=str(dest/'lib')
 copy=dict(manifest,runtime=runtime_identity(dest/'RANSPreservationProbe'),scope='same corrected compiled snapshot; alias name is not a legacy operator')
 write_json(dest/'source-build.json',copy)
os.environ['LD_LIBRARY_PATH']=str(r/'runtime/lib')
old=r.with_name('a93-a95-e94e95-nonlinear-20261110')
cmd=['python3','-B',str(s/'check_rans_native_ownership.py'),'--probe',str(r/'runtime/RANSPreservationProbe'),'--descriptor-probe',str(r/'runtime/RANSNewModelsProbe'),'--old-probe',str(old/'runtime-resolution-v3/RANSPreservationProbe'),'--previous-root',str(old),'--output',str(r/'geometric-regressions-v1')]
with (r/'geometric-regressions-v1.log').open('x') as f:
 proc=None
 while proc is None:
  proc=spawn(r,cmd,stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ))
  if proc is None:time.sleep(3)
 code=proc.wait()
write_json(r/'geometry-execution.json',dict(exit_code=code,command=cmd))
assert code==0
cmd=['python3','-B',str(s/'run_rans_ownership_short.py'),'--root',str(r),'--geometry-receipt',str(r/'geometric-regressions-v1/results.json'),'--concurrency','2']
with (r/'short-controller-v1.log').open('x') as f:
 proc=subprocess.Popen(cmd,stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ))
 write_json(r/'short-controller-launch.json',dict(pid=proc.pid,command=cmd))
 code=proc.wait()
write_json(r/'short-controller-completion.json',dict(exit_code=code))
