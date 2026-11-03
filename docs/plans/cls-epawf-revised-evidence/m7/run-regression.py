import fcntl,hashlib,json,os,subprocess,time
from pathlib import Path
from datetime import datetime,timezone
root=Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised')
source=root/'source';main=Path('/media/jay/新加卷/nektar-v2026');out=root/'m7-regression-001'
# The already-running M7 campaign owns the canonical source lock. This
# additional lease prohibits subsequent source synchronization until CTest ends.
lock=(root/'m7-regression.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
def git(path,*args):return subprocess.check_output(['git','-C',str(path),*args],text=True).strip()
expected='bf202be017bab7aceeedf9a48d86bcfeeca0580f'
assert git(source,'rev-parse','HEAD')==expected and not git(source,'status','--porcelain')
assert not git(main,'status','--porcelain')
main_revision=git(main,'rev-parse','HEAD');out.mkdir(exist_ok=False)
env=os.environ.copy();threads={k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')}
env.update(threads);env['PYTHONDONTWRITEBYTECODE']='1';env['TMPDIR']=str(root/'tmp')
env['LD_LIBRARY_PATH']=str(root/'build/geometry-lib')+':'+env.get('LD_LIBRARY_PATH','')
command=['taskset','-c','2-23','ctest','--test-dir',str(root/'build'),'--parallel','12','-E','^RANSCLSEPAWFM7','--output-on-failure','--output-junit',str(out/'ctest.xml')]
receipt=dict(schema=1,source_commit=expected,main_commit=main_revision,main_status='',connection='tailscale:100.105.247.110',started_utc=datetime.now(timezone.utc).isoformat(),command=command,thread_environment=threads,ctest_allocated_ranks=12,overlapping_M7_campaign_maximum_ranks=3,combined_maximum_allocated_ranks=15,compiler_workers_during_tests=0,worker_limit=72,source_write_guard='canonical M7 campaign lock followed by m7-regression.lock; parent waits for both before any source mutation',runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
(out/'processes-before.txt').write_text(subprocess.check_output(['ps','-eo','pid,nlwp,pcpu,comm'],text=True))
start=time.monotonic()
try:
 with (out/'ctest.log').open('w') as log:r=subprocess.run(command,env=env,cwd=out,stdout=log,stderr=subprocess.STDOUT,timeout=14400)
 code=r.returncode
except subprocess.TimeoutExpired:code=124
receipt.update(status='passed' if code==0 else 'failed',exit_code=code,seconds=time.monotonic()-start,finished_utc=datetime.now(timezone.utc).isoformat(),source_end_commit=git(source,'rev-parse','HEAD'),source_end_status=git(source,'status','--porcelain'),main_end_commit=git(main,'rev-parse','HEAD'),main_end_status=git(main,'status','--porcelain'),log_sha256=hashlib.sha256((out/'ctest.log').read_bytes()).hexdigest())
assert receipt['source_end_commit']==expected and not receipt['source_end_status']
assert receipt['main_end_commit']==main_revision and not receipt['main_end_status']
(out/'processes-after.txt').write_text(subprocess.check_output(['ps','-eo','pid,nlwp,pcpu,comm'],text=True))
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2),flush=True)
print((out/'ctest.log').read_text()[-3000:],flush=True)
raise SystemExit(code)
