"""Bounded unqualified-flow diagnosis; never a developed terrain campaign."""
import fcntl,hashlib,json,os,subprocess,time
from datetime import datetime,timezone
from pathlib import Path
root=Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8')
source,build=root/'source',root/'build';out=root/'flow-001';out.mkdir(exist_ok=False)
lock=(root/'validation.lock').open('a');fcntl.flock(lock,fcntl.LOCK_SH|fcntl.LOCK_NB)
revision=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()
assert revision=='aa854f4c5b05451a81a4e1dd60eaaca98eb547e4'
env=os.environ.copy();env.update({k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')})
env['PYTHONDONTWRITEBYTECODE']='1';env['TMPDIR']=str(root/'tmp');env['LD_LIBRARY_PATH']=str(build/'geometry-lib')+':'+env.get('LD_LIBRARY_PATH','')
fingerprint=source/'docs/plans/cls-epawf-revised-evidence/m7/fingerprint-runtime.py'
def pin(name,compare=None):
 cmd=['python3','-B',str(fingerprint),'--root',str(root),'--output',str(out/(name+'.json'))]
 if compare:cmd+=['--compare',str(out/(compare+'.json'))]
 with (out/(name+'.log')).open('w') as log:subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=240,check=True)
pin('runtime-before')
command=['taskset','-c','36','python3','-B',str(source/'solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_epawf_m8.py'),
 '--build',str(build),'--output',str(out/'qualification'),'--mode','flow','--diagnostic',
 '--variant','Regularized','--beta','.4','--gamma-scale','1','--timeout','3600']
receipt=dict(schema=1,status='running',source_commit=revision,command=command,maximum_ranks=1,
 numerical_threads=1,diagnostic_only=True,candidate_qualified=False,experimental_accuracy_claimed=False,
 runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 started_utc=datetime.now(timezone.utc).isoformat())
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');start=time.monotonic()
try:
 with (out/'flow.log').open('w') as log:code=subprocess.run(command,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=18000).returncode
except subprocess.TimeoutExpired:code=124
pin('runtime-after','runtime-before')
assert subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()==revision
assert not subprocess.check_output(['git','-C',str(source),'status','--porcelain'],text=True).strip()
receipt.update(status='passed' if code==0 else 'failed',exit_code=code,seconds=time.monotonic()-start,
 runtime_unchanged=True,finished_utc=datetime.now(timezone.utc).isoformat())
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2));print((out/'flow.log').read_text()[-2500:]);raise SystemExit(code)
