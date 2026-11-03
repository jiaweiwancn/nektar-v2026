#!/usr/bin/env python3
"""Publish a local documentation closure to the isolated Linux worktree."""
import argparse
from datetime import datetime,timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--bundle',type=Path,required=True);p.add_argument('--revision',required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--qualified',type=Path,required=True);p.add_argument('--fingerprint-tool',type=Path,required=True)
a=p.parse_args();root=a.root.resolve();source=root/'source';main=Path('/media/jay/新加卷/nektar-v2026');out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
def git(path,*args):return subprocess.check_output(['git','-c','gc.auto=0','-c','pack.threads=1','-c','checkout.workers=1','-C',str(path),*args],text=True).strip()
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
locks=[]
for name in ('m3-validation.lock','m7-regression.lock'):
 lock=(root/name).open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);locks.append(lock)
qualified=json.loads(a.qualified.read_text());assert qualified['success']
assert git(source,'rev-parse','HEAD')==qualified['source_commit']
for path in (source,main):assert not git(path,'status','--porcelain')
main_revision=git(main,'rev-parse','HEAD')
subprocess.run(['git','-c','gc.auto=0','-c','pack.threads=1','-C',str(main),'fetch',str(a.bundle),'refs/heads/codex/cls-epawf-revised:refs/remotes/cls-epawf/m7-closed'],check=True)
assert git(main,'rev-parse','refs/remotes/cls-epawf/m7-closed')==a.revision
changed=git(main,'diff','--name-only',qualified['source_commit'],a.revision).splitlines()
assert all(path.startswith('docs/') or path=='README.md' or (path.startswith('solvers/IncNavierStokesSolver/Examples/') and path.endswith('.md')) for path in changed),changed
subprocess.run(['git','-c','checkout.workers=1','-C',str(source),'checkout','--detach',a.revision],check=True)
env=os.environ.copy();env.update({k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')});env['PYTHONDONTWRITEBYTECODE']='1';env['TMPDIR']=str(root/'tmp');env['LD_LIBRARY_PATH']=str(root/'build/geometry-lib')+':'+env.get('LD_LIBRARY_PATH','')
cmd=['python3','-B',str(a.fingerprint_tool),'--root',str(root),'--output',str(out/'final-runtime-fingerprint.json'),'--compare',str(a.qualified)]
with (out/'fingerprint.log').open('w') as log:result=subprocess.run(cmd,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=600)
assert result.returncode==0,(out/'fingerprint.log').read_text()
fingerprint=json.loads((out/'final-runtime-fingerprint.json').read_text());assert fingerprint['success'] and fingerprint['identical_qualification_inputs_and_runtime']
assert git(source,'rev-parse','HEAD')==a.revision and not git(source,'status','--porcelain')
assert git(main,'rev-parse','HEAD')==main_revision and not git(main,'status','--porcelain')
receipt=dict(schema=1,status='passed',final_commit=a.revision,qualified_commit=qualified['source_commit'],changed_paths=changed,identical_qualification_inputs_and_runtime=True,numerical_source_inputs=len(fingerprint['source_inputs']),runtime_images=len(fingerprint['runtime_artifacts']),main_commit=main_revision,main_preserved=True,source_clean=True,connection='tailscale:100.105.247.110',command=cmd,fingerprint_sha256=sha(out/'final-runtime-fingerprint.json'),runner_sha256=sha(__file__),finished_utc=datetime.now(timezone.utc).isoformat())
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2),flush=True)
