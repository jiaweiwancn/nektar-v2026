"""Complete native gates after a receipt-only fix; replays are explicit."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime,timezone
import fcntl,hashlib,json,os,subprocess,time
from pathlib import Path
root=Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8')
source=root/'source';audit_source=root/'verification/source';build=root/'build';out=root/'native-completion-001'
out.mkdir(exist_ok=False)
lock=(root/'validation.lock').open('a');fcntl.flock(lock,fcntl.LOCK_SH|fcntl.LOCK_NB)
def git(path,*args):return subprocess.check_output(['git','-C',str(path),*args],text=True).strip()
assert git(source,'rev-parse','HEAD')=='aa854f4c5b05451a81a4e1dd60eaaca98eb547e4'
assert git(audit_source,'rev-parse','HEAD')=='6284787bf4c436203f4404498ba3dfb68aa144ef'
assert not git(source,'status','--porcelain') and not git(audit_source,'status','--porcelain')
prior=json.loads((root/'campaign-001/receipt.json').read_text())
for name in ('scaled-native','regularized-native'):
 assert prior['jobs'][name]['status']=='failed', 'Wait until the producing native cases stop'
 assert 'PosixPath is not JSON serializable' in (root/'campaign-001'/(name+'.log')).read_text()
env=os.environ.copy();env.update({k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')})
env['PYTHONDONTWRITEBYTECODE']='1';env['TMPDIR']=str(root/'tmp');env['LD_LIBRARY_PATH']=str(build/'geometry-lib')+':'+env.get('LD_LIBRARY_PATH','')
fingerprint=source/'docs/plans/cls-epawf-revised-evidence/m7/fingerprint-runtime.py'
def pin(name,compare=None):
 command=['python3','-B',str(fingerprint),'--root',str(root),'--output',str(out/(name+'.json'))]
 if compare:command+=['--compare',str(out/(compare+'.json'))]
 with (out/(name+'.log')).open('w') as log:subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=240,check=True)
def audit_inputs():
 files=git(audit_source,'ls-files').splitlines()
 return {name:hashlib.sha256((audit_source/name).read_bytes()).hexdigest() for name in files
  if Path(name).suffix in ('.py','.xml','.cpp','.h','.hpp','.cmake') or Path(name).name=='CMakeLists.txt'}
inputs=audit_inputs();(out/'audit-source-inputs.json').write_text(json.dumps(inputs,sort_keys=True,indent=2)+'\n')
pin('runtime-before')
receipt=dict(schema=1,status='running',compiled_source_commit=git(source,'rev-parse','HEAD'),
 audit_source_commit=git(audit_source,'rev-parse','HEAD'),maximum_concurrent_ranks=4,numerical_threads=1,
 native_execution_vs_replay='Only verified successful byte-identical prior cases reused, every independent audit repeated; checkpoint/restart and missing cases execute afresh',
 runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),audit_input_count=len(inputs),
 started_utc=datetime.now(timezone.utc).isoformat(),jobs={})
def save():(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
driver=audit_source/'solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_epawf_m8.py'
jobs=[]
for name,variant,gamma,cpus in [('scaled-native','CubicScaled',0,'0-3'),('regularized-native','Regularized',1,'4-7')]:
 command=['taskset','-c',cpus,'python3','-B',str(driver),'--build',str(build),'--output',str(out/name),
  '--mode','native','--variant',variant,'--beta','.4','--gamma-scale',str(gamma),
  '--reuse',str(root/'campaign-001'/name),'--reuse-runtime-pin',str(root/'campaign-001/runtime-before.json')]
 receipt['jobs'][name]=dict(status='running',command=command,maximum_ranks=2,timeout=14400)
 jobs.append((name,command))
save()
def execute(job):
 name,command=job;start=time.monotonic()
 try:
  with (out/(name+'.log')).open('w') as log:code=subprocess.run(command,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=14400).returncode
 except subprocess.TimeoutExpired:code=124
 return name,dict(status='passed' if code==0 else 'failed',exit_code=code,seconds=time.monotonic()-start,
  log_sha256=hashlib.sha256((out/(name+'.log')).read_bytes()).hexdigest())
with ThreadPoolExecutor(max_workers=2) as pool:
 for future in as_completed([pool.submit(execute,job) for job in jobs]):
  name,result=future.result();receipt['jobs'][name].update(result);save();print(name,json.dumps(result),flush=True)
  print((out/(name+'.log')).read_text()[-1500:],flush=True)
pin('runtime-after','runtime-before')
assert audit_inputs()==inputs and not git(audit_source,'status','--porcelain')
assert git(source,'rev-parse','HEAD')==receipt['compiled_source_commit'] and not git(source,'status','--porcelain')
receipt.update(status='passed' if all(j['exit_code']==0 for j in receipt['jobs'].values()) else 'failed',
 runtime_unchanged=True,audit_inputs_unchanged=True,finished_utc=datetime.now(timezone.utc).isoformat())
save();print(json.dumps(receipt,indent=2),flush=True);raise SystemExit(0 if receipt['status']=='passed' else 1)
