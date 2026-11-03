#!/usr/bin/env python3
"""After all jobs end, sync the primary worktree without changing qualified physics."""
import argparse
from datetime import datetime,timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--revision',required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--fingerprint-tool',type=Path,required=True)
a=p.parse_args();root=a.root.resolve();source=root/'source';audit=root/'m7-audit-source';main=Path('/media/jay/新加卷/nektar-v2026');out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
def git(path,*args):return subprocess.check_output(['git','-c','gc.auto=0','-c','pack.threads=1','-c','checkout.workers=1','-C',str(path),*args],text=True).strip()
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
locks=[]
for name in ('m3-validation.lock','m7-regression.lock'):
 lock=(root/name).open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);locks.append(lock)
native=json.loads((root/'m7-native-final-001/receipt.json').read_text());legacy=json.loads((root/'m7-regression-001/receipt.json').read_text())
assert native['status']==legacy['status']=='passed' and native.get('finished_utc') and legacy.get('finished_utc')
assert git(audit,'rev-parse','HEAD')==a.revision==native['source_commit']
for path in (source,audit,main):assert not git(path,'status','--porcelain'),str(path)
main_revision=git(main,'rev-parse','HEAD');pin=json.loads((root/'m7-attempt-004/runtime-fingerprint.json').read_text())
receipt=dict(schema=1,source_before=git(source,'rev-parse','HEAD'),source_after=a.revision,main_commit=main_revision,started_utc=datetime.now(timezone.utc).isoformat(),connection='tailscale:100.105.247.110',compiler_workers=16,numerical_threads=1,commands=[],runner_sha256=sha(__file__))
subprocess.run(['git','-c','gc.auto=0','-c','checkout.workers=1','-C',str(source),'checkout','--detach',a.revision],check=True)
allowed={str(p.relative_to(source)) for p in (source/'solvers/IncNavierStokesSolver/Tests/RANS').glob('*m7*.py')}
for name,expected in pin['source_inputs'].items():
 if name not in allowed:assert sha(source/name)==expected,name
changed=git(source,'diff','--name-only',pin['source_commit'],a.revision,'--','library','solvers/IncNavierStokesSolver').splitlines()
assert not [name for name in changed if Path(name).suffix in ('.cpp','.h','.hpp','.cc','.cxx','.c')]
receipt['changed_paths']=changed
for key,expected in pin['runtime_artifacts'].items():assert sha(pin['runtime_resolved_paths'][key])==expected,key
env=os.environ.copy();env.update({k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')});env['PYTHONDONTWRITEBYTECODE']='1';env['TMPDIR']=str(root/'tmp');env['LD_LIBRARY_PATH']=str(root/'build/geometry-lib')+':'+env.get('LD_LIBRARY_PATH','')
commands=[('configure',['cmake','-S',str(source/'solvers/IncNavierStokesSolver/Tests/RANS'),'-B',str(root/'build'),'-DRANS_BUILD_SOLVER=ON','-DCMAKE_BUILD_TYPE=Release','-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake']),
 ('incremental-build',['cmake','--build',str(root/'build'),'--parallel','16']),
 ('fingerprint',['python3','-B',str(a.fingerprint_tool),'--root',str(root),'--output',str(out/'qualified-runtime-fingerprint.json')])]
for name,cmd in commands:
 start=time.monotonic()
 with (out/(name+'.log')).open('w') as log:result=subprocess.run(cmd,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=900)
 receipt['commands'].append(dict(name=name,argv=cmd,exit_code=result.returncode,seconds=time.monotonic()-start,log_sha256=sha(out/(name+'.log'))))
 assert result.returncode==0,name
for key,expected in pin['runtime_artifacts'].items():assert sha(pin['runtime_resolved_paths'][key])==expected,key
qualified=json.loads((out/'qualified-runtime-fingerprint.json').read_text())
assert qualified['success'] and qualified['source_commit']==a.revision
for name in ('runtime_artifacts','runtime_bindings','runtime_environment','python_packages','tool_images'):
 assert qualified[name]==pin[name],'Campaign-to-qualified runtime context changed: '+name
assert git(main,'rev-parse','HEAD')==main_revision and not git(main,'status','--porcelain')
assert git(source,'rev-parse','HEAD')==a.revision and not git(source,'status','--porcelain')
receipt.update(status='passed',unchanged_runtime_images=len(pin['runtime_artifacts']),identical_campaign_runtime_context=True,
               compared_runtime_fields=['runtime_artifacts','runtime_bindings','runtime_environment','python_packages','tool_images'],
               main_preserved=True,source_clean=True,finished_utc=datetime.now(timezone.utc).isoformat())
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2),flush=True)
