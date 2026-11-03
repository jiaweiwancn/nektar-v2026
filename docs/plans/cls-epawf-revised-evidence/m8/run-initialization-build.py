"""Build the corrected capped initializer without changing the AA campaign."""
from datetime import datetime,timezone
import hashlib,json,os,subprocess,time
from pathlib import Path
base=Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8')
root=base/'initialization';source=root/'source';build=root/'build';out=root/'build-001';out.mkdir(exist_ok=False)
def git(path,*args):return subprocess.check_output(['git','-C',str(path),*args],text=True).strip()
assert git(source,'rev-parse','HEAD')=='ed34052fe4dd967e963c22272a0ec0888cf35b70' and not git(source,'status','--porcelain')
assert git(base/'source','rev-parse','HEAD')=='aa854f4c5b05451a81a4e1dd60eaaca98eb547e4'
old=json.loads((base/'campaign-001/runtime-before.json').read_text())
def old_images():
 for key,value in old['runtime_artifacts'].items():
  assert hashlib.sha256(Path(old['runtime_resolved_paths'][key]).read_bytes()).hexdigest()==value,key
old_images();(root/'tmp').mkdir(exist_ok=True)
env=os.environ.copy();env.update({k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')})
env['PYTHONDONTWRITEBYTECODE']='1';env['TMPDIR']=str(root/'tmp');env['LD_LIBRARY_PATH']=str(build/'geometry-lib')+':'+env.get('LD_LIBRARY_PATH','')
receipt=dict(schema=1,status='running',source_commit=git(source,'rev-parse','HEAD'),compiler_workers=16,
 connection='tailscale:100.105.247.110',unchanged_primary_commit=git(base/'source','rev-parse','HEAD'),
 original_runtime_images=len(old['runtime_artifacts']),
 runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 started_utc=datetime.now(timezone.utc).isoformat(),commands=[])
def save():(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
save();commands=[('configure',['cmake','-S',str(source/'solvers/IncNavierStokesSolver/Tests/RANS'),'-B',str(build),
 '-DRANS_BUILD_SOLVER=ON','-DCMAKE_BUILD_TYPE=Release','-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake'],240),
 ('build',['cmake','--build',str(build),'--target','IncNavierStokesSolver','--parallel','16'],1800),
 ('fingerprint',['python3','-B',str(source/'docs/plans/cls-epawf-revised-evidence/m7/fingerprint-runtime.py'),
  '--root',str(root),'--output',str(out/'runtime-fingerprint.json')],240)]
code=0
for name,command,timeout in commands:
 start=time.monotonic()
 try:
  with (out/(name+'.log')).open('w') as log:code=subprocess.run(command,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=timeout).returncode
 except subprocess.TimeoutExpired:code=124
 receipt['commands'].append(dict(name=name,command=command,exit_code=code,seconds=time.monotonic()-start,
  timeout=timeout,log_sha256=hashlib.sha256((out/(name+'.log')).read_bytes()).hexdigest()));save()
 print(name,code,flush=True);print((out/(name+'.log')).read_text()[-1500:],flush=True)
 if code:break
old_images();assert not git(source,'status','--porcelain') and not git(base/'source','status','--porcelain')
assert git(base/'source','rev-parse','HEAD')==receipt['unchanged_primary_commit']
receipt.update(status='passed' if code==0 else 'failed',exit_code=code,original_runtime_images_unchanged=True,
 finished_utc=datetime.now(timezone.utc).isoformat());save();raise SystemExit(code)
