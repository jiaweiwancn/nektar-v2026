"""Complete the failed strict build after its formatting-only source fix."""
from datetime import datetime, timezone
import fcntl,hashlib,json,os,subprocess,time
from pathlib import Path
r1=Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8-r1')
root=r1/'controlled';source,build,out=root/'source',root/'build',root/'build-002';out.mkdir(exist_ok=False)
lock=(root/'validation.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX)
git=lambda *a:subprocess.check_output(['git','-C',str(source),*a],text=True).strip()
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert git('rev-parse','HEAD')=='23cbbd5c25136d2a611db0ee5a01589e65f9781f' and not git('status','--porcelain')
subprocess.run(['git','-C',str(source),'fetch','--no-tags',str(r1/'controlled-format.bundle'),'HEAD'],check=True)
revision=git('rev-parse','FETCH_HEAD');assert revision.startswith('360cbb7')
assert git('diff','--name-only','HEAD',revision)=='solvers/IncNavierStokesSolver/EquationSystems/RANSEquilibriumWallEvolution.cpp'
subprocess.run(['git','-C',str(source),'checkout','--detach',revision],check=True)
env=dict(os.environ);env.update({key:'1' for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')})
env.update(PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(root/'tmp'),LD_LIBRARY_PATH=str(build/'geometry-lib'))
receipt=dict(schema=1,status='running',source_commit=revision,compiler_workers=16,total_worker_bound=19,hard_limit=72,
    prior_attempt_preserved=str(root/'build-001'),reason='Split an if/output and declaration to satisfy unchanged -Werror=misleading-indentation; no equations, tolerances or flags changed',
    runner_sha256=sha(__file__),commands=[],started_utc=datetime.now(timezone.utc).isoformat())
def save():(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
commands=[('build',['taskset','-c','0-31','cmake','--build',str(build),'--parallel','16']),
    ('fingerprint',['python3','-B',str(source/'docs/plans/cls-epawf-revised-evidence/m7/fingerprint-runtime.py'),'--root',str(root),'--output',str(out/'runtime-fingerprint.json')])]
save();code=0
for name,argv in commands:
    start=time.monotonic();log=out/(name+'.log')
    try:
        with log.open('w') as stream:code=subprocess.run(argv,cwd=out,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=1800 if name=='build' else 240).returncode
    except subprocess.TimeoutExpired:code=124
    receipt['commands'].append(dict(name=name,argv=argv,exit_code=code,seconds=time.monotonic()-start,log_sha256=sha(log)));save()
    print(name,code,log.read_text()[-2000:],flush=True)
    if code:break
assert git('rev-parse','HEAD')==revision and not git('status','--porcelain')
old=json.loads((r1/'build-001/runtime-fingerprint.json').read_text())
for key,value in old['runtime_artifacts'].items():assert sha(old['runtime_resolved_paths'][key])==value
receipt.update(status='passed' if code==0 else 'failed',first_R1_images_preserved=True,source_clean_end=True,
               finished_utc=datetime.now(timezone.utc).isoformat());save();raise SystemExit(code)
