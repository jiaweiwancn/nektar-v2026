"""One frozen focused operator suite, with actual source/runtime pins."""
import argparse,fcntl,hashlib,json,os,subprocess,time
from datetime import datetime,timezone
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--suite',choices=('inner','pressure'),required=True);a=p.parse_args()
base=Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised');root=base/'m8-r2'
source,build=root/'source',root/'build';out=root/(a.suite+'-verification-001');out.mkdir(exist_ok=False)
lock=(root/'validation.lock').open('a');fcntl.flock(lock,fcntl.LOCK_SH)
sha=lambda x:hashlib.sha256(Path(x).read_bytes()).hexdigest()
git=lambda *x:subprocess.check_output(['git','-C',str(source),*x],text=True).strip()
assert git('rev-parse','HEAD')=='1ee3a8454c857c807bcb6911024186afb9d22721' and not git('status','--porcelain')
assert json.loads((root/'build-003/receipt.json').read_text())['status']=='passed'
env=dict(os.environ);keys=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')
env.update({k:'1' for k in keys});env.update(PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(root/'tmp'),LD_LIBRARY_PATH=str(build/'geometry-lib'))
tests=source/'solvers/IncNavierStokesSolver/Tests/RANS';cpu=50 if a.suite=='inner' else 48
if a.suite=='inner':
    driver=tests/'validate_cls_epawf_m8_r2_inner.py'
    command=['taskset','-c',str(cpu),'python3','-B',str(driver),'--kernel-probe',str(build/'RANSWallInnerStressProbe'),
        '--wall-probe',str(build/'RANSCurvedWallProbe'),'--template',str(tests/'EquilibriumWallEvolution.xml'),
        '--output',str(out/'results')]
else:
    driver=tests/'validate_cls_epawf_m8_r2_pressure.py'
    command=['taskset','-c',str(cpu),'python3','-B',str(driver),'--build',str(build),
        '--r1-trace',str(base/'m8-r1/controlled/reproduction-001/case/case.nonlinear.json'),
        '--output',str(out/'results')]
receipt=dict(schema=1,status='running',suite=a.suite,source_commit=git('rev-parse','HEAD'),source_clean_start=True,
    runner_sha256=sha(__file__),driver_sha256=sha(driver),ranks=1,cpu=cpu,
    numerical_threads={k:env[k] for k in keys},trajectory_calls=0,commands=[],started_utc=datetime.now(timezone.utc).isoformat())
def save():(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
def run(name,argv,timeout=240):
    started=time.monotonic();log=out/(name+'.log')
    try:
        with log.open('w') as stream:code=subprocess.run(argv,cwd=out,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=timeout).returncode
    except subprocess.TimeoutExpired:code=124
    receipt['commands'].append(dict(name=name,argv=argv,cwd=str(out),exit_code=code,seconds=time.monotonic()-started,log_sha256=sha(log)));save();return code
fp=source/'docs/plans/cls-epawf-revised-evidence/m7/fingerprint-runtime.py';save()
code=run('fingerprint-before',['python3','-B',str(fp),'--root',str(root),'--output',str(out/'runtime-before.json')])
assert code==0
before=json.loads((out/'runtime-before.json').read_text());built=json.loads((root/'build-003/runtime-fingerprint.json').read_text())
for k in ('source_inputs','runtime_artifacts','runtime_bindings','python_packages','tool_images','runtime_environment'):assert before[k]==built[k],k
code=run('suite',command,1800)
after=run('fingerprint-after',['python3','-B',str(fp),'--root',str(root),'--output',str(out/'runtime-after.json'),
    '--compare',str(out/'runtime-before.json')])
receipt.update(status='passed' if code==0 and after==0 else 'failed',suite_exit_code=code,runtime_unchanged=after==0,
    source_clean_end=not bool(git('status','--porcelain')),source_commit_end=git('rev-parse','HEAD'),finished_utc=datetime.now(timezone.utc).isoformat());save()
print(json.dumps({k:receipt[k] for k in ('suite','status','suite_exit_code','runtime_unchanged')}),flush=True)
print((out/'suite.log').read_text()[-4000:],flush=True);raise SystemExit(0 if receipt['status']=='passed' else 1)
