"""Execute independent M8 candidate gates against one immutable Linux build."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

p=argparse.ArgumentParser();p.add_argument('--output',required=True);p.add_argument('--revision',required=True)
a=p.parse_args()
root=Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8')
source,build=root/'source',root/'build';out=root/a.output;out.mkdir(exist_ok=False)
lock=(root/'validation.lock').open('a');fcntl.flock(lock,fcntl.LOCK_SH|fcntl.LOCK_NB)
def git(*args):return subprocess.check_output(['git','-C',str(source),*args],text=True).strip()
assert git('rev-parse','HEAD')==a.revision and not git('status','--porcelain')
env=os.environ.copy();env.update({k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')})
env['PYTHONDONTWRITEBYTECODE']='1';env['TMPDIR']=str(root/'tmp')
env['LD_LIBRARY_PATH']=str(build/'geometry-lib')+':'+env.get('LD_LIBRARY_PATH','')
driver=source/'solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_epawf_m8.py'
fingerprint=source/'docs/plans/cls-epawf-revised-evidence/m7/fingerprint-runtime.py'
def pin(name,compare=None):
 command=['python3','-B',str(fingerprint),'--root',str(root),'--output',str(out/(name+'.json'))]
 if compare:command+=['--compare',str(out/(compare+'.json'))]
 with (out/(name+'.log')).open('w') as log:
  subprocess.run(command,env=env,cwd=out,stdout=log,stderr=subprocess.STDOUT,timeout=240,check=True)
pin('runtime-before')
jobs=[]
for name,variant,mode,gamma,cpus in [('scaled-native','CubicScaled','native',0,'0-3'),
    ('regularized-native','Regularized','native',1,'4-7'),
    ('regularized-mms','Regularized','mms',1,'8-11')]:
 jobs.append((name,['taskset','-c',cpus,'python3','-B',str(driver),'--build',str(build),
  '--output',str(out/name),'--mode',mode,'--variant',variant,'--beta','.4','--gamma-scale',str(gamma)],21600,2))
jobs.append(('legacy-regression',['taskset','-c','12-35','ctest','--test-dir',str(build),'--parallel','12',
 '-E','^RANSCLSEPAWFM8','--output-on-failure','--output-junit',str(out/'ctest.xml')],21600,12))
receipt=dict(schema=1,status='running',source_commit=a.revision,compiler_workers=0,
 maximum_concurrent_test_ranks=18,independent_certificate_workers=1,total_worker_bound=19,
 numerical_threads={k:env[k] for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')},
 connection='tailscale:100.105.247.110',started_utc=datetime.now(timezone.utc).isoformat(),
 runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),jobs={})
def save():(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
for name,command,timeout,ranks in jobs:receipt['jobs'][name]=dict(status='pending',command=command,timeout=timeout,maximum_ranks=ranks)
save()
def execute(job):
 name,command,timeout,ranks=job;start=time.monotonic()
 try:
  with (out/(name+'.log')).open('w') as log:
   code=subprocess.run(command,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=timeout).returncode
 except subprocess.TimeoutExpired:code=124
 return name,dict(status='passed' if code==0 else 'failed',exit_code=code,seconds=time.monotonic()-start,
  finished_utc=datetime.now(timezone.utc).isoformat(),log_sha256=hashlib.sha256((out/(name+'.log')).read_bytes()).hexdigest())
with ThreadPoolExecutor(max_workers=4) as executor:
 futures={executor.submit(execute,job):job[0] for job in jobs}
 for future in as_completed(futures):
  name,result=future.result();receipt['jobs'][name].update(result);save()
  print(name,json.dumps(result),flush=True);print((out/(name+'.log')).read_text()[-1200:],flush=True)
pin('runtime-after','runtime-before')
assert git('rev-parse','HEAD')==a.revision and not git('status','--porcelain')
receipt.update(status='passed' if all(j['exit_code']==0 for j in receipt['jobs'].values()) else 'failed',
 runtime_unchanged=True,source_clean_end=True,finished_utc=datetime.now(timezone.utc).isoformat())
save();print(json.dumps(receipt,indent=2),flush=True)
raise SystemExit(0 if receipt['status']=='passed' else 1)
