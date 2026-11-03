import hashlib,json,os,subprocess,time
from pathlib import Path
from datetime import datetime,timezone
root=Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised');source=root/'m7-audit-source';out=root/'m7-native-final-001';out.mkdir(exist_ok=False)
def git(*args):return subprocess.check_output(['git','-C',str(source),*args],text=True).strip()
revision=git('rev-parse','HEAD');assert not git('status','--porcelain')
changed=git('diff','--name-only','bf202be','--','library','solvers/IncNavierStokesSolver').splitlines()
assert not [p for p in changed if Path(p).suffix in ('.cpp','.h','.hpp','.cc','.cxx','.c')]
pin=json.loads((root/'m7-attempt-004/runtime-fingerprint.json').read_text())
def runtime_check():
 for key,expected in pin['runtime_artifacts'].items():assert hashlib.sha256(Path(pin['runtime_resolved_paths'][key]).read_bytes()).hexdigest()==expected,key
runtime_check()
env=os.environ.copy();env.update({k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')});env['PYTHONDONTWRITEBYTECODE']='1';env['TMPDIR']=str(root/'tmp');env['LD_LIBRARY_PATH']=str(root/'build/geometry-lib')+':'+env.get('LD_LIBRARY_PATH','')
cmd=['taskset','-c','24-27','python3','-B',str(source/'solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_epawf_m7.py'),'--build',str(root/'build'),'--output',str(out/'qualification')]
receipt=dict(schema=1,source_commit=revision,compiled_physics_commit=pin['source_commit'],compiled_source_differences=[],unchanged_runtime_images=len(pin['runtime_artifacts']),thread_environment={k:env[k] for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')},maximum_allocated_ranks=2,connection='tailscale:100.105.247.110',started_utc=datetime.now(timezone.utc).isoformat(),command=cmd,runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');start=time.monotonic()
try:
 with (out/'qualification.log').open('w') as log:r=subprocess.run(cmd,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=14400)
 code=r.returncode
except subprocess.TimeoutExpired:code=124
runtime_check();assert git('rev-parse','HEAD')==revision and not git('status','--porcelain')
receipt.update(status='passed' if code==0 else 'failed',exit_code=code,seconds=time.monotonic()-start,finished_utc=datetime.now(timezone.utc).isoformat(),source_clean_end=True,log_sha256=hashlib.sha256((out/'qualification.log').read_bytes()).hexdigest())
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2),flush=True);print((out/'qualification.log').read_text()[-3000:],flush=True)
raise SystemExit(code)
