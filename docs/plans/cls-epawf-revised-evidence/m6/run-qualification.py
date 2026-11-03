import argparse, fcntl, hashlib, json, os, shutil, subprocess, time
from pathlib import Path
from datetime import datetime, timezone

p=argparse.ArgumentParser()
p.add_argument('bundle');p.add_argument('reference');p.add_argument('attempt')
p.add_argument('--mode', choices=['build','geometry','suite','oracle','regression','note','replay','m6-pilot','m6-diagnose','m6-kernel','m6-all','m6-probes'], default='build')
p.add_argument('--wait-for-lock',action='store_true')
args=p.parse_args()
root=Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised')
main=Path('/media/jay/新加卷/nektar-v2026')
source=root/'source';build=root/'build'
lock=(root/'m3-validation.lock').open('w')
if args.wait_for_lock:print('Waiting for the previous source-bound validation to finish',flush=True)
fcntl.flock(lock, fcntl.LOCK_EX|(0 if args.wait_for_lock else fcntl.LOCK_NB))
def git(path,*args): return subprocess.check_output(['git','-c','gc.auto=0','-c','pack.threads=1','-C',str(path),*args],text=True).strip()
assert not git(source,'status','--porcelain') and not git(main,'status','--porcelain')
subprocess.run(['git','-c','gc.auto=0','-c','pack.threads=1','-C',str(main),'fetch',str(root/'transfers'/args.bundle),
                'refs/heads/codex/cls-epawf-revised:refs/remotes/cls-epawf/'+args.reference],check=True)
revision=git(main,'rev-parse','refs/remotes/cls-epawf/'+args.reference)
subprocess.run(['git','-c','gc.auto=0','-c','pack.threads=1','-C',str(source),'checkout','--detach',revision],check=True)
out=root/args.attempt;out.mkdir(exist_ok=False)
env=os.environ.copy()
threads={k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')}
env.update(threads);env['TMPDIR']=str(root/'tmp');env['PYTHONDONTWRITEBYTECODE']='1'
env['LD_LIBRARY_PATH']=str(build/'geometry-lib')+':'+env.get('LD_LIBRARY_PATH','')
receipt={'schema':1,'source_commit':revision,'source_status':git(source,'status','--porcelain'),
         'main_commit':git(main,'rev-parse','HEAD'),'main_status':git(main,'status','--porcelain'),
         'started_utc':datetime.now(timezone.utc).isoformat(),'compiler_workers':16,
         'test_workers':1,'mpi_ranks':[1,2,4],'maximum_concurrent_test_ranks':4,'thread_environment':threads,
         'connection':'tailscale:100.105.247.110','runner_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'commands':[]}
(out/'processes-before.txt').write_text(subprocess.check_output(['ps','-eo','pid,nlwp,pcpu,comm'],text=True))
commands=[('configure',['cmake','-S',str(source/'solvers/IncNavierStokesSolver/Tests/RANS'),'-B',str(build),
                       '-DRANS_BUILD_SOLVER=ON','-DCMAKE_BUILD_TYPE=Release',
                       '-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake'],180)]
cmd=['cmake','--build',str(build),'--parallel','16']

commands.append(('build',cmd,900))
if args.mode.startswith('m6-'):
 mode=args.mode[3:]
 if mode=='diagnose':
  commands.append(('generate-initial',['python3',str(source/'solvers/IncNavierStokesSolver/Tests/RANS/cls_epawf_m6_cases.py'),str(out/'initial.xml'),'--family','flat','--steps','0'],120))
  commands.append(('expressions',[str(build/'RANSWallExpressionProbe'),str(out/'initial.xml'),str(out/'expressions.json')],120))
  commands.append(('initial',[str(build/'IncNavierStokesSolver'),str(out/'initial.xml')],180))
  commands.append(('kernel',['taskset','-c','0','python3',str(source/'solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_epawf_m6_kernel.py'),'--probe',str(build/'RANSEquilibriumTransitionProbe'),'--output',str(out/'kernel')],1800))
 elif mode=='pilot':
  commands.append(('generate-pilot',['python3',str(source/'solvers/IncNavierStokesSolver/Tests/RANS/cls_epawf_m6_cases.py'),str(out/'pilot.xml'),'--family','flat','--steps','1'],120))
  commands.append(('pilot',[str(build/'IncNavierStokesSolver'),str(out/'pilot.xml')],1800))
 else:
  script='validate_cls_epawf_m6'+('_kernel' if mode=='kernel' else '')+'.py'
  options=['--probe',str(build/'RANSEquilibriumTransitionProbe')] if mode=='kernel' else ['--build',str(build),'--mode',mode]
  commands.append(('m6-'+mode,['taskset','-c','0','python3',str(source/'solvers/IncNavierStokesSolver/Tests/RANS'/script),*options,'--output',str(out/'qualification')],10800))
elif args.mode=='geometry':
 commands.append(('geometry',[str(build/'RANSCurvedGeometryProbe'),'--output',str(out/'geometry.json')],180))
elif args.mode=='replay':
 commands.append(('replay',['taskset','-c','0','python3',str(source/'solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_epawf_m5_session_replay.py'),
                           '--build',str(build),'--qualification',str(root/'m5-attempt-005/oracle'),'--output',str(out/'replay')],1800))
elif args.mode=='note':
 commands.append(('profiles-native',[str(build/'RANSEquilibriumProfilesProbe'),'--output',str(out/'profiles.json')],180))
 commands.append(('note',['taskset','-c','0','python3',str(source/'solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_epawf_m5_note.py'),
                         '--probe',str(build/'RANSCurvedNoteProbe'),'--profiles',str(out/'profiles.json'),'--output',str(out/'note')],180))
elif args.mode=='oracle':
 commands.append(('oracle',['taskset','-c','0','python3',str(source/'solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_epawf_m5.py'),
                           '--build',str(build),'--output',str(out/'oracle')],10800))
if args.mode in ('all','suite','regression'):
 suite=['taskset','-c','0','ctest','--test-dir',str(build),'--parallel','1','--output-on-failure','--output-junit',str(out/'ctest.xml')]
 if args.mode=='regression':suite+=['-E','^RANSCLSEPAWFM6']
 commands.append(('suite',suite,10800))
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
receipt['source_end_commit']=git(source,'rev-parse','HEAD')
receipt['source_end_status']=git(source,'status','--porcelain')
receipt['main_end_commit']=git(main,'rev-parse','HEAD')
receipt['main_end_status']=git(main,'status','--porcelain')
assert receipt['source_end_commit']==revision and not receipt['source_end_status']
assert receipt['main_end_commit']==receipt['main_commit'] and not receipt['main_end_status']
(out/'processes-after.txt').write_text(subprocess.check_output(['ps','-eo','pid,nlwp,pcpu,comm'],text=True))
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
