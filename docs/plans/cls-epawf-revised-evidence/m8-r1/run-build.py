"""Isolated Linux build of the committed R1 diagnostics; no trajectories."""
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

base=Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised')
root=base/'m8-r1';main=Path('/media/jay/新加卷/nektar-v2026')
source,build,out=root/'source',root/'build',root/'build-001'
out.mkdir(exist_ok=False)
lock=(root/'validation.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX)
def git(path,*argv):return subprocess.check_output(['git','-C',str(path),*argv],text=True).strip()
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
assert git(main,'rev-parse','HEAD')=='f4b70baec10016dd098005bc5a5507b87e4a57bf' and not git(main,'status','--porcelain')
assert git(root/'archive-source','rev-parse','HEAD')=='ab5c3c6fee19761e8d88692e8e6c331e5889234e' and not git(root/'archive-source','status','--porcelain')
subprocess.run(['git','-C',str(main),'fetch','--no-tags',str(root/'diagnostics.bundle'),'HEAD'],check=True)
revision=git(main,'rev-parse','FETCH_HEAD');assert revision.startswith('426d48c')
subprocess.run(['git','-C',str(main),'worktree','add','--detach',str(source),revision],check=True)
env=dict(os.environ)
env.update({k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')})
env.update(PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(root/'tmp'),LD_LIBRARY_PATH=str(build/'geometry-lib'))
receipt=dict(schema=1,status='running',source_commit=revision,source_clean_start=True,
    compiler_workers=16,maximum_concurrent_analysis_workers=3,total_worker_bound=19,hard_limit=72,
    numerical_threads={k:env[k] for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')},
    runner_sha256=sha(__file__),bundle_sha256=sha(root/'diagnostics.bundle'),connection='tailscale:100.105.247.110',
    commands=[],started_utc=datetime.now(timezone.utc).isoformat())
def save():(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
commands=[('configure',['cmake','-S',str(source/'solvers/IncNavierStokesSolver/Tests/RANS'),'-B',str(build),
    '-DRANS_BUILD_SOLVER=ON','-DCMAKE_BUILD_TYPE=Release','-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake']),
    ('build',['taskset','-c','0-31','cmake','--build',str(build),'--parallel','16']),
    ('fingerprint',['python3','-B',str(source/'docs/plans/cls-epawf-revised-evidence/m7/fingerprint-runtime.py'),
                   '--root',str(root),'--output',str(out/'runtime-fingerprint.json')])]
save();code=0
for name,argv in commands:
    started=time.monotonic();log=out/(name+'.log')
    try:
        with log.open('w') as stream:code=subprocess.run(argv,cwd=out,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=1800 if name=='build' else 240).returncode
    except subprocess.TimeoutExpired:code=124
    receipt['commands'].append(dict(name=name,argv=argv,exit_code=code,seconds=time.monotonic()-started,log_sha256=sha(log)));save()
    print(name,code,log.read_text()[-2000:],flush=True)
    if code:break
assert git(source,'rev-parse','HEAD')==revision and not git(source,'status','--porcelain')
assert git(main,'rev-parse','HEAD')=='f4b70baec10016dd098005bc5a5507b87e4a57bf' and not git(main,'status','--porcelain')
receipt.update(status='passed' if code==0 else 'failed',source_clean_end=True,
               finished_utc=datetime.now(timezone.utc).isoformat());save()
raise SystemExit(code)
