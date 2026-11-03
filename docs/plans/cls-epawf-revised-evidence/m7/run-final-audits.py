import hashlib,json,os,subprocess,time
from pathlib import Path
from datetime import datetime,timezone
root=Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised');source=root/'m7-audit-source';build=root/'build';out=root/'m7-final-audits-001'
out.mkdir(exist_ok=False)
def git(*args):return subprocess.check_output(['git','-C',str(source),*args],text=True).strip()
revision=git('rev-parse','HEAD');assert not git('status','--porcelain')
env=os.environ.copy();env.update({k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')});env['PYTHONDONTWRITEBYTECODE']='1';env['TMPDIR']=str(root/'tmp');env['LD_LIBRARY_PATH']=str(build/'geometry-lib')+':'+env.get('LD_LIBRARY_PATH','')
tests=source/'solvers/IncNavierStokesSolver/Tests/RANS'
commands=[('native-replay',[str(tests/'validate_cls_epawf_m7.py'),'--build',str(build),'--output',str(root/'m7-replay-003'),'--replay',str(root/'m7-attempt-004/native')]),
 ('dense-guard',[str(tests/'cls_epawf_m7_dense_guard.py'),'--build',str(build),'--output',str(root/'m7-dense-guard-003')]),
 ('transverse',[str(tests/'validate_cls_epawf_m7_transverse.py'),'--build',str(build),'--output',str(root/'m7-strict-transverse-002')]),
 ('active-energy',[str(root/'transfers/cls_epawf_m7_energy_replay_last.py')]),
 ('flow-replay',[str(source/'docs/plans/cls-epawf-revised-evidence/m7/replay-flow-diagnosis.py'),'--source',str(source),'--input',str(root/'m7-flow-diagnosis-001/validation.json'),'--output',str(root/'m7-flow-replay-001')])]
receipt=dict(schema=1,source_commit=revision,connection='tailscale:100.105.247.110',started_utc=datetime.now(timezone.utc).isoformat(),runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),commands=[],allocated_workers=1,thread_environment={k:env[k] for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')})
for name,args in commands:
 cmd=['taskset','-c','24','python3','-B',*args];start=time.monotonic()
 try:
  with (out/(name+'.log')).open('w') as log:r=subprocess.run(cmd,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=1800)
  code=r.returncode
 except subprocess.TimeoutExpired:code=124
 receipt['commands'].append(dict(name=name,argv=cmd,exit_code=code,seconds=time.monotonic()-start,log_sha256=hashlib.sha256((out/(name+'.log')).read_bytes()).hexdigest()))
 receipt['status']='passed' if all(x['exit_code']==0 for x in receipt['commands']) else 'failed'
 (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(name,code,flush=True);print((out/(name+'.log')).read_text()[-2000:],flush=True)
 if code:break
assert git('rev-parse','HEAD')==revision and not git('status','--porcelain')
receipt['finished_utc']=datetime.now(timezone.utc).isoformat();receipt['source_clean_end']=True
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
raise SystemExit(0 if receipt['status']=='passed' else 1)
