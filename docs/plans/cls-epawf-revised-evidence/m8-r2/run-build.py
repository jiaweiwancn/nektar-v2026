"""First isolated R2 Release build; no solver/trajectory execution."""
from datetime import datetime, timezone
import hashlib, json, os, subprocess, time
from pathlib import Path

base=Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised')
root=base/'m8-r2';main=Path('/media/jay/新加卷/nektar-v2026')
source,build,out=root/'source',root/'build',root/'build-001'
out.mkdir(exist_ok=False);(root/'tmp').mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
git=lambda p,*a:subprocess.check_output(['git','-C',str(p),*a],text=True).strip()
main_head=git(main,'rev-parse','HEAD');assert main_head=='f4b70baec10016dd098005bc5a5507b87e4a57bf'
assert not git(main,'status','--porcelain')
historical={}
for p in (base/'m8/source',base/'m8-r1/source',base/'m8-r1/controlled/source',base/'m8-r1/finalization-001/source'):
    historical[str(p)]=dict(commit=git(p,'rev-parse','HEAD'),status=git(p,'status','--porcelain'))
    assert not historical[str(p)]['status']
old_images={str(p):sha(p) for p in (base/'m8/build/IncNavierStokesSolver',
    base/'m8-r1/build/IncNavierStokesSolver',base/'m8-r1/controlled/build/IncNavierStokesSolver')}
subprocess.run(['git','-C',str(main),'fetch','--no-tags',str(root/'initial.bundle'),'HEAD'],check=True)
revision=git(main,'rev-parse','FETCH_HEAD');assert revision=='26d137f296bf42c39f1ef0e10ac99f8a2e2ef7cb'
subprocess.run(['git','-C',str(main),'worktree','add','--detach',str(source),revision],check=True)
env=dict(os.environ)
keys=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')
env.update({key:'1' for key in keys})
env.update(PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(root/'tmp'),LD_LIBRARY_PATH=str(build/'geometry-lib'))
receipt=dict(schema=1,status='running',source_commit=revision,source_clean_start=True,
    main_checkout=dict(commit=main_head,branch=git(main,'branch','--show-current'),clean=True),
    historical_sources=historical,historical_solver_images=old_images,
    compiler_workers=16,maximum_concurrent_analysis_workers=3,total_worker_ceiling=19,hard_worker_limit=72,
    numerical_threads={key:env[key] for key in keys},connection='tailscale:100.105.247.110',
    runner_sha256=sha(__file__),bundle_sha256=sha(root/'initial.bundle'),commands=[],
    started_utc=datetime.now(timezone.utc).isoformat())
def save():(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
commands=[('configure',['cmake','-S',str(source/'solvers/IncNavierStokesSolver/Tests/RANS'),'-B',str(build),
    '-DRANS_BUILD_SOLVER=ON','-DCMAKE_BUILD_TYPE=Release','-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake']),
    ('build',['taskset','-c','0-31','cmake','--build',str(build),'--parallel','16']),
    ('fingerprint',['python3','-B',str(source/'docs/plans/cls-epawf-revised-evidence/m7/fingerprint-runtime.py'),
                    '--root',str(root),'--output',str(out/'runtime-fingerprint.json')])]
save();code=0
for name,argv in commands:
    start=time.monotonic();log=out/(name+'.log')
    try:
        with log.open('w') as stream:code=subprocess.run(argv,cwd=out,env=env,stdout=stream,stderr=subprocess.STDOUT,
            timeout=1800 if name=='build' else 240).returncode
    except subprocess.TimeoutExpired:code=124
    receipt['commands'].append(dict(name=name,argv=argv,exit_code=code,seconds=time.monotonic()-start,log_sha256=sha(log)));save()
    print(name,code,log.read_text()[-3000:],flush=True)
    if code:break
receipt.update(status='passed' if code==0 else 'failed',source_commit_end=git(source,'rev-parse','HEAD'),
    source_clean_end=not bool(git(source,'status','--porcelain')),main_preserved=git(main,'rev-parse','HEAD')==main_head and not git(main,'status','--porcelain'),
    historical_sources_preserved=all(git(Path(p),'rev-parse','HEAD')==value['commit'] and not git(Path(p),'status','--porcelain') for p,value in historical.items()),
    historical_solver_images_preserved=all(sha(p)==value for p,value in old_images.items()),finished_utc=datetime.now(timezone.utc).isoformat())
save();assert receipt['main_preserved'] and receipt['historical_sources_preserved'] and receipt['historical_solver_images_preserved']
raise SystemExit(code)
