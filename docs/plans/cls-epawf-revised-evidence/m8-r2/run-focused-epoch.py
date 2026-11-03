"""One focused suite bound to a clean, immutable R2 source/build epoch."""
import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--root',type=Path,required=True)
p.add_argument('--expected-source',required=True)
p.add_argument('--suite',choices=('inner','pressure','legacy'),required=True)
p.add_argument('--attempt',required=True)
p.add_argument('--analysis-source',type=Path)
p.add_argument('--analysis-commit')
p.add_argument('--reuse-inner-run',type=Path)
p.add_argument('--reuse-extra-run',type=Path,action='append',default=[])
p.add_argument('--capped-pressure-run',type=Path)
p.add_argument('--analysis-only',action='store_true')
a=p.parse_args()
root=a.root.resolve();source,build=root/'source',root/'build'
out=root/a.attempt;out.mkdir(exist_ok=False)
lock=(root/'validation.lock').open('a');fcntl.flock(lock,fcntl.LOCK_SH)
sha=lambda x:hashlib.sha256(Path(x).read_bytes()).hexdigest()
git=lambda *x:subprocess.check_output(['git','-C',str(source),*x],text=True).strip()
assert git('rev-parse','HEAD')==a.expected_source and not git('status','--porcelain')
assert json.loads((root/'build-receipt.json').read_text())['status']=='passed'
keys=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')
env=dict(os.environ,**{k:'1' for k in keys})
env.update(PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(root/'tmp'),LD_LIBRARY_PATH=str(build/'geometry-lib'))
tests=source/'solvers/IncNavierStokesSolver/Tests/RANS'
analysis_source=a.analysis_source.resolve() if a.analysis_source else source
analysis_git=lambda *x:subprocess.check_output(['git','-C',str(analysis_source),*x],text=True).strip()
if a.analysis_source:
    assert a.analysis_commit and analysis_git('rev-parse','HEAD')==a.analysis_commit
    assert not analysis_git('status','--porcelain')
analysis_tests=analysis_source/'solvers/IncNavierStokesSolver/Tests/RANS'
if a.suite=='inner':
    driver=analysis_tests/'validate_cls_epawf_m8_r2_inner.py';cpu=50
    command=['taskset','-c',str(cpu),'python3','-B',str(driver),'--kernel-probe',str(build/'RANSWallInnerStressProbe'),
        '--wall-probe',str(build/'RANSCurvedWallProbe'),'--template',str(tests/'EquilibriumWallEvolution.xml'),
        '--producing-source-commit',a.expected_source,'--output',str(out/'results')]
    if a.reuse_inner_run:
        command+=['--reuse-run',str(a.reuse_inner_run)]
        if not a.analysis_only:command+=['--continue-missing']
    if a.analysis_only:assert a.reuse_inner_run
    for retained in a.reuse_extra_run:command+=['--reuse-extra-run',str(retained)]
    if a.capped_pressure_run:command+=['--capped-pressure-run',str(a.capped_pressure_run)]
elif a.suite=='pressure':
    driver=tests/'validate_cls_epawf_m8_r2_pressure.py';cpu=48
    base=Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised')
    command=['taskset','-c',str(cpu),'python3','-B',str(driver),'--build',str(build),
        '--r1-trace',str(base/'m8-r1/controlled/reproduction-001/case/case.nonlinear.json'),'--output',str(out/'results')]
else:
    driver=source/'solvers/IncNavierStokesSolver/Tests/RANS/CMakeLists.txt';cpu=49
    command=['taskset','-c',str(cpu),'ctest','--test-dir',str(build),'--parallel','1','--output-on-failure',
        '-R','^(RANSEquilibriumTimeNative|RANSEquilibriumEvolutionNative|RANSEquilibriumWallEvolutionSmoke)$']
receipt=dict(schema=2,status='running',suite=a.suite,source_commit=git('rev-parse','HEAD'),source_clean_start=True,
    analysis_source=str(analysis_source),analysis_commit=analysis_git('rev-parse','HEAD'),analysis_clean_start=True,
    runner_sha256=sha(__file__),driver_sha256=sha(driver),ranks=1,cpu=cpu,
    numerical_threads={k:env[k] for k in keys},native_physical_trajectories=0,
    scope='focused operator/action or existing isolated legacy verification; no R2 qualification evolution',
    commands=[],started_utc=datetime.now(timezone.utc).isoformat())
def save():(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
def run(name,argv,timeout):
    started=time.monotonic();log=out/(name+'.log')
    try:
        with log.open('w') as stream:
            code=subprocess.run(argv,cwd=out,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=timeout).returncode
    except subprocess.TimeoutExpired:code=124
    receipt['commands'].append(dict(name=name,argv=argv,cwd=str(out),exit_code=code,
        seconds=time.monotonic()-started,log_sha256=sha(log)));save()
    return code
fp=source/'docs/plans/cls-epawf-revised-evidence/m7/fingerprint-runtime.py';save()
code=run('fingerprint-before',['python3','-B',str(fp),'--root',str(root),'--output',str(out/'runtime-before.json')],240)
assert code==0
before=json.loads((out/'runtime-before.json').read_text());built=json.loads((root/'runtime-fingerprint.json').read_text())
for k in ('source_inputs','runtime_artifacts','runtime_bindings','python_packages','tool_images','runtime_environment'):
    assert before[k]==built[k],k
code=run('suite',command,1800)
after=run('fingerprint-after',['python3','-B',str(fp),'--root',str(root),'--output',str(out/'runtime-after.json'),
    '--compare',str(out/'runtime-before.json')],240)
receipt.update(status='passed' if code==0 and after==0 else 'failed',suite_exit_code=code,runtime_unchanged=after==0,
    source_clean_end=not bool(git('status','--porcelain')),source_commit_end=git('rev-parse','HEAD'),
    analysis_clean_end=not bool(analysis_git('status','--porcelain')),analysis_commit_end=analysis_git('rev-parse','HEAD'),
    finished_utc=datetime.now(timezone.utc).isoformat())
if not (receipt['source_clean_end'] and receipt['source_commit_end']==a.expected_source and
        receipt['analysis_clean_end'] and receipt['analysis_commit_end']==receipt['analysis_commit']):
    receipt['status']='failed'
save()
print(json.dumps({k:receipt[k] for k in ('suite','status','suite_exit_code','runtime_unchanged')}),flush=True)
print((out/'suite.log').read_text()[-4000:],flush=True)
raise SystemExit(0 if receipt['status']=='passed' else 1)
