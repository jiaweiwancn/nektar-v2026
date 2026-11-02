import argparse, fcntl, hashlib, json, os, shutil, subprocess, time
from pathlib import Path
from datetime import datetime, timezone

p=argparse.ArgumentParser()
p.add_argument('bundle');p.add_argument('reference');p.add_argument('attempt')
p.add_argument('--mode', choices=['native','solver','development','oracle','all','suite'], default='native')
p.add_argument('--wait-for-lock',action='store_true')
args=p.parse_args()
root=Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised')
main=Path('/media/jay/新加卷/nektar-v2026')
source=root/'source';build=root/'build'
lock=(root/'m3-validation.lock').open('w')
if args.wait_for_lock:print('Waiting for the previous source-bound validation to finish',flush=True)
fcntl.flock(lock, fcntl.LOCK_EX|(0 if args.wait_for_lock else fcntl.LOCK_NB))
def git(path,*args): return subprocess.check_output(['git','-C',str(path),*args],text=True).strip()
assert not git(source,'status','--porcelain') and not git(main,'status','--porcelain')
subprocess.run(['git','-C',str(main),'fetch',str(root/'transfers'/args.bundle),
                'refs/heads/codex/cls-epawf-revised:refs/remotes/cls-epawf/'+args.reference],check=True)
revision=git(main,'rev-parse','refs/remotes/cls-epawf/'+args.reference)
subprocess.run(['git','-C',str(source),'checkout','--detach',revision],check=True)
out=root/args.attempt;out.mkdir(exist_ok=False)
env=os.environ.copy()
threads={k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')}
env.update(threads);env['TMPDIR']=str(root/'tmp');env['PYTHONDONTWRITEBYTECODE']='1'
env['LD_LIBRARY_PATH']=str(build/'geometry-lib')+':'+env.get('LD_LIBRARY_PATH','')
receipt={'schema':1,'source_commit':revision,'source_status':git(source,'status','--porcelain'),
         'main_commit':git(main,'rev-parse','HEAD'),'main_status':git(main,'status','--porcelain'),
         'started_utc':datetime.now(timezone.utc).isoformat(),'compiler_workers':16,
         'test_workers':1,'mpi_ranks':1,'thread_environment':threads,
         'runner_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'commands':[]}
(out/'processes-before.txt').write_text(subprocess.check_output(['ps','-eo','pid,nlwp,pcpu,comm'],text=True))
commands=[('configure',['cmake','-S',str(source/'solvers/IncNavierStokesSolver/Tests/RANS'),'-B',str(build),
                       '-DRANS_BUILD_SOLVER=ON','-DCMAKE_BUILD_TYPE=Release',
                       '-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake'],180)]
cmd=['cmake','--build',str(build),'--parallel','16']
if args.mode=='native':cmd+=['--target','RANSEquilibriumPatchProbe']
commands.append(('build',cmd,900))
if args.mode in ('native','solver'):
 shutil.copyfile(source/'solvers/IncNavierStokesSolver/Tests/RANS/EquilibriumWallFlat.xml',out/'EquilibriumWallFlat.xml')
 if args.mode=='native':
  commands.append(('native',['taskset','-c','0',str(build/'RANSEquilibriumPatchProbe'),
                  str(out/'EquilibriumWallFlat.xml'),str(out/'patch.json')],600))
 else:
  commands.append(('solver',['taskset','-c','0',str(build/'IncNavierStokesSolver'),str(out/'EquilibriumWallFlat.xml')],900))
elif args.mode=='development':
 commands.append(('oracle-residuals',['taskset','-c','0','python3',str(source/'solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_epawf_m3.py'),
                           '--build',str(build),'--output',str(out/'oracle'),'--suite','residuals'],900))
 shutil.copyfile(source/'solvers/IncNavierStokesSolver/Tests/RANS/EquilibriumWallFlat.xml',out/'EquilibriumWallFlat.xml')
 commands.append(('solver',['taskset','-c','0',str(build/'IncNavierStokesSolver'),str(out/'EquilibriumWallFlat.xml')],900))
elif args.mode!='suite':
 commands.append(('oracle',['taskset','-c','0','python3',str(source/'solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_epawf_m3.py'),
                           '--build',str(build),'--output',str(out/'oracle')],1800))
if args.mode in ('all','suite'):
 commands.append(('suite',['taskset','-c','0','ctest','--test-dir',str(build),'--parallel','1','--output-on-failure'],1800))
for name,cmd,timeout in commands:
 start=time.monotonic()
 try:
  with (out/(name+'.log')).open('w') as f:
   result=subprocess.run(cmd,cwd=out,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=timeout)
  code=result.returncode
 except subprocess.TimeoutExpired:code=124
 receipt['commands'].append({'command':cmd,'exit_code':code,'seconds':time.monotonic()-start,
                             'timeout_seconds':timeout,'log_sha256':hashlib.sha256((out/(name+'.log')).read_bytes()).hexdigest()})
 receipt['status']='passed' if code==0 else 'failed'
 (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(name,code,flush=True)
 log=(out/(name+'.log')).read_text()
 print('\n'.join(l for l in log.splitlines() if 'error:' in l)[-4000:] if code else log[-1500:],flush=True)
 if code:
  print(log[-3000:],flush=True)
  break
receipt['finished_utc']=datetime.now(timezone.utc).isoformat()
if args.mode=='suite' and receipt['status']=='passed':
 reports=[]
 for path in build.glob('cls-epawf-m3-*/validation.json'):
  data=json.loads(path.read_text())
  if data['source_commit']==revision and data['verification_status']=='passed' and data['milestone_closure']:
   reports.append(path)
 assert len(reports)==1
 (out/'oracle').symlink_to(reports[0].parent,target_is_directory=True)
 receipt['qualification_report']=str(reports[0])
 receipt['qualification_execution']='complete registered CTest; no duplicate direct oracle run'
receipt['source_end_commit']=git(source,'rev-parse','HEAD')
receipt['source_end_status']=git(source,'status','--porcelain')
receipt['main_end_commit']=git(main,'rev-parse','HEAD')
receipt['main_end_status']=git(main,'status','--porcelain')
assert receipt['source_end_commit']==revision and not receipt['source_end_status']
assert receipt['main_end_commit']==receipt['main_commit'] and not receipt['main_end_status']
(out/'processes-after.txt').write_text(subprocess.check_output(['ps','-eo','pid,nlwp,pcpu,comm'],text=True))
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
