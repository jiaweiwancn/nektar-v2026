"""Build the second diagnostic epoch, preserving the first R1 binary/source."""
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

base=Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised')
r1=base/'m8-r1';root=r1/'controlled';root.mkdir(exist_ok=False);(root/'tmp').mkdir()
main=Path('/media/jay/新加卷/nektar-v2026');source,build,out=root/'source',root/'build',root/'build-001';out.mkdir()
lock=(root/'validation.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX)
git=lambda path,*a:subprocess.check_output(['git','-C',str(path),*a],text=True).strip()
sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
old=json.loads((r1/'build-001/runtime-fingerprint.json').read_text())
def preserve():
    assert git(r1/'source','rev-parse','HEAD')=='426d48cd44f9f4af017244f6be4c1a6f62d36c01' and not git(r1/'source','status','--porcelain')
    for name,digest in old['runtime_artifacts'].items():assert sha(old['runtime_resolved_paths'][name])==digest,name
    assert git(main,'rev-parse','HEAD')=='f4b70baec10016dd098005bc5a5507b87e4a57bf' and not git(main,'status','--porcelain')
preserve()
subprocess.run(['git','-C',str(main),'fetch','--no-tags',str(r1/'controlled.bundle'),'HEAD'],check=True)
revision=git(main,'rev-parse','FETCH_HEAD');assert revision.startswith('23cbbd5')
subprocess.run(['git','-C',str(main),'worktree','add','--detach',str(source),revision],check=True)
env=dict(os.environ);env.update({key:'1' for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')})
env.update(PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(root/'tmp'),LD_LIBRARY_PATH=str(build/'geometry-lib'))
receipt=dict(schema=1,status='running',source_commit=revision,source_clean_start=True,
    compiler_workers=16,maximum_concurrent_analysis_workers=3,total_worker_bound=19,hard_limit=72,
    runner_sha256=sha(__file__),bundle_sha256=sha(r1/'controlled.bundle'),connection='tailscale:100.105.247.110',
    preserved_first_R1_images=len(old['runtime_artifacts']),commands=[],started_utc=datetime.now(timezone.utc).isoformat())
def save():(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
commands=[('configure',['cmake','-S',str(source/'solvers/IncNavierStokesSolver/Tests/RANS'),'-B',str(build),
    '-DRANS_BUILD_SOLVER=ON','-DCMAKE_BUILD_TYPE=Release','-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake']),
    ('build',['taskset','-c','0-31','cmake','--build',str(build),'--parallel','16']),
    ('fingerprint',['python3','-B',str(source/'docs/plans/cls-epawf-revised-evidence/m7/fingerprint-runtime.py'),'--root',str(root),'--output',str(out/'runtime-fingerprint.json')])]
save();code=0
for name,argv in commands:
    started=time.monotonic();log=out/(name+'.log')
    try:
        with log.open('w') as stream:code=subprocess.run(argv,cwd=out,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=1800 if name=='build' else 240).returncode
    except subprocess.TimeoutExpired:code=124
    receipt['commands'].append(dict(name=name,argv=argv,exit_code=code,seconds=time.monotonic()-started,log_sha256=sha(log)));save()
    print(name,code,log.read_text()[-1800:],flush=True)
    if code:break
preserve();assert git(source,'rev-parse','HEAD')==revision and not git(source,'status','--porcelain')
receipt.update(status='passed' if code==0 else 'failed',first_R1_source_and_images_preserved=True,source_clean_end=True,
               finished_utc=datetime.now(timezone.utc).isoformat());save();raise SystemExit(code)
