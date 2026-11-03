set -eo pipefail
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export LD_LIBRARY_PATH="/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8/initialization/build/geometry-lib:${LD_LIBRARY_PATH:-}"
export TMPDIR='/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8/initialization/tmp'
taskset -c 48 python3 -B - <<'PY'
import datetime,fcntl,hashlib,json,os,pathlib,subprocess,sys,time
root=pathlib.Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8/initialization')
source=root/'source';build=root/'build';out=root/'capped-native-001';work=root/'capped-native-001-runner';work.mkdir(exist_ok=False)
expected='ed34052fe4dd967e963c22272a0ec0888cf35b70'
lock=(root/'validation.lock').open('a');fcntl.flock(lock,fcntl.LOCK_SH|fcntl.LOCK_NB)
def git(*args):return subprocess.check_output(['git','-C',str(source),*args],text=True).strip()
def sha(path):
 h=hashlib.sha256()
 with pathlib.Path(path).open('rb') as stream:
  for block in iter(lambda:stream.read(1048576),b''):h.update(block)
 return h.hexdigest()
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
receipt=dict(schema=1,milestone='M8',status='running',started_utc=stamp(),source_commit=git('rev-parse','HEAD'),source_clean_start=not git('status','--porcelain'),connection='tailscale:100.105.247.110',allocated_workers=1,cpu_affinity=sorted(os.sched_getaffinity(0)),numerical_threads={k:os.environ[k] for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')},maximum_reserved_server_workers=38,lock=dict(path=str(root/'validation.lock'),mode='shared',held_throughout=True),commands=[],failures=[],physical_robustness_qualified=False)
receipt_path=work/'receipt.json'
def save():receipt_path.write_text(json.dumps(receipt,indent=2)+'\n')
def command(name,argv,timeout):
 log=work/(name+'.log');start=time.monotonic()
 try:
  with log.open('w') as stream:code=subprocess.run(list(map(str,argv)),cwd=work,env=os.environ.copy(),stdout=stream,stderr=subprocess.STDOUT,timeout=timeout).returncode
 except subprocess.TimeoutExpired:code=124
 receipt['commands'].append(dict(name=name,argv=list(map(str,argv)),exit_code=code,seconds=time.monotonic()-start,timeout=timeout,log=str(log),log_sha256=sha(log)));save()
 return code
fingerprint=source/'docs/plans/cls-epawf-revised-evidence/m7/fingerprint-runtime.py'
def pin(name,compare=None):
 argv=[sys.executable,'-B',fingerprint,'--root',root,'--output',work/(name+'.json')]
 if compare is not None:argv+=['--compare',work/(compare+'.json')]
 assert command(name,argv,300)==0,name
 return json.loads((work/(name+'.json')).read_text())
code=1;before=None;started=time.monotonic();save()
try:
 assert receipt['source_commit']==expected and receipt['source_clean_start']
 assert not out.exists(),'Fresh capped output exists already'
 baseline_path=root/'build-002/runtime-fingerprint.json';baseline=json.loads(baseline_path.read_text())
 assert baseline['success'] and baseline['source_commit']==expected
 receipt['build002_fingerprint_sha256']=sha(baseline_path)
 before=pin('runtime-before');assert before['success'] and before['source_commit']==expected
 for key in ('source_inputs','runtime_artifacts','runtime_bindings','python_packages','tool_images'):
  assert before[key]==baseline[key],key+' differs from the successful full build'
 receipt['runtime_image_count']=len(before['runtime_artifacts']);receipt['source_input_count']=len(before['source_inputs'])
 receipt['environment_differences_from_build002']={key:dict(build002=value,current=before['runtime_environment'].get(key)) for key,value in baseline['runtime_environment'].items() if value!=before['runtime_environment'].get(key)}
 driver=source/'solvers/IncNavierStokesSolver/Tests/RANS/cls_epawf_m8_capped_native.py'
 receipt['driver_sha256']=sha(driver);receipt['solver_sha256']=sha(build/'IncNavierStokesSolver')
 argv=[sys.executable,'-B',driver,'--build',build,'--output',out,'--timeout','1800']
 code=command('capped-native',argv,4200)
 if (out/'validation.json').exists():
  result=json.loads((out/'validation.json').read_text());receipt.update(qualification_receipt=str(out/'validation.json'),qualification_sha256=sha(out/'validation.json'),mechanism_verified=result.get('mechanism_verified'),evolution_verified=result.get('evolution_verified'),case_count=len(result.get('cases',[])),qualification_status=result.get('status'),qualification_failures=result.get('failures'))
 receipt['status']='passed' if code==0 else 'failed'
except Exception as error:
 receipt['failures'].append(str(error));code=1;receipt['status']='failed'
try:
 if before is not None:
  after=pin('runtime-after','runtime-before');assert after['success'] and after['identical_qualification_inputs_and_runtime'];receipt['runtime_unchanged']=True
 receipt.update(source_commit_end=git('rev-parse','HEAD'),source_clean_end=not git('status','--porcelain'))
 assert receipt['source_commit_end']==expected and receipt['source_clean_end']
except Exception as error:
 receipt['failures'].append('End provenance: '+str(error));code=1;receipt['status']='failed'
receipt.update(exit_code=code,finished_utc=stamp(),seconds=time.monotonic()-started);save()
print(json.dumps(receipt,indent=2),flush=True)
if (work/'capped-native.log').exists():print((work/'capped-native.log').read_text()[-6000:],flush=True)
fcntl.flock(lock,fcntl.LOCK_UN);lock.close()
raise SystemExit(code)
PY
