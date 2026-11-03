"""Retained, explicit completion attempt after a diagnosed R2 build error."""
import argparse,fcntl,hashlib,json,os,subprocess,time
from datetime import datetime,timezone
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--bundle',type=Path,required=True);p.add_argument('--expected',required=True)
p.add_argument('--previous-source',required=True);p.add_argument('--attempt',required=True)
p.add_argument('--allowed-path',action='append',required=True);a=p.parse_args()
root=Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8-r2')
source,build,out=root/'source',root/'build',root/a.attempt;out.mkdir(exist_ok=False)
lock=(root/'validation.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX)
sha=lambda x:hashlib.sha256(Path(x).read_bytes()).hexdigest()
git=lambda *x:subprocess.check_output(['git','-C',str(source),*x],text=True).strip()
assert git('rev-parse','HEAD')==a.previous_source and not git('status','--porcelain')
subprocess.run(['git','-C',str(source),'fetch','--no-tags',str(a.bundle),'HEAD'],check=True)
assert git('rev-parse','FETCH_HEAD')==a.expected
changed=git('diff','--name-only',a.previous_source,a.expected).splitlines();assert sorted(changed)==sorted(a.allowed_path)
subprocess.run(['git','-C',str(source),'checkout','--detach',a.expected],check=True)
env=dict(os.environ);keys=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')
env.update({k:'1' for k in keys});env.update(PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(root/'tmp'),LD_LIBRARY_PATH=str(build/'geometry-lib'))
receipt=dict(schema=1,status='running',source_commit=a.expected,previous_source=a.previous_source,
    explicit_changed_paths=changed,runner_sha256=sha(__file__),bundle_sha256=sha(a.bundle),
    compiler_workers=16,maximum_analysis_workers=3,total_worker_ceiling=19,hard_worker_limit=72,
    numerical_threads={k:env[k] for k in keys},commands=[],started_utc=datetime.now(timezone.utc).isoformat())
def save():(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
save();code=0
for name,argv in [('build',['taskset','-c','0-31','cmake','--build',str(build),'--parallel','16']),
    ('fingerprint',['python3','-B',str(source/'docs/plans/cls-epawf-revised-evidence/m7/fingerprint-runtime.py'),
       '--root',str(root),'--output',str(out/'runtime-fingerprint.json')])]:
    start=time.monotonic();log=out/(name+'.log')
    try:
        with log.open('w') as stream:code=subprocess.run(argv,cwd=out,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=1800 if name=='build' else 240).returncode
    except subprocess.TimeoutExpired:code=124
    receipt['commands'].append(dict(name=name,argv=argv,exit_code=code,seconds=time.monotonic()-start,log_sha256=sha(log)));save()
    print(name,code,log.read_text()[-4000:],flush=True)
    if code:break
receipt.update(status='passed' if code==0 else 'failed',source_clean_end=not bool(git('status','--porcelain')),
    source_commit_end=git('rev-parse','HEAD'),finished_utc=datetime.now(timezone.utc).isoformat());save()
raise SystemExit(code)
