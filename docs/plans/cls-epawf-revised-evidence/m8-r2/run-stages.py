"""Run only the predeclared matched stage pair, or one conditional scalar test."""
import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import struct
import subprocess
import sys
import time

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--root',type=Path,required=True)
p.add_argument('--expected-source',required=True)
p.add_argument('--preflight',type=Path,required=True)
p.add_argument('--original',type=Path,required=True)
p.add_argument('--representation',choices=('outer6','outer9'),required=True)
p.add_argument('--mode',choices=('matched','frozen-scalar'),default='matched')
p.add_argument('--matched-receipt',type=Path)
p.add_argument('--matched-assessment',type=Path,
    help='Strict metadata assessment of the original matched receipt; never permits a coupled retry')
p.add_argument('--output',type=Path,required=True)
a=p.parse_args()
root=a.root.resolve();source,build=root/'source',root/'build';out=a.output.resolve()
assert out.is_relative_to(root) and not out.is_relative_to(source)
out.mkdir(exist_ok=False)
lock=(root/'validation.lock').open('a');fcntl.flock(lock,fcntl.LOCK_SH)
sha=lambda x:hashlib.sha256(Path(x).read_bytes()).hexdigest()
vector_sha=lambda x:hashlib.sha256(struct.pack('<'+str(len(x))+'d',*x)).hexdigest()
git=lambda *x:subprocess.check_output(['git','-C',str(source),*x],text=True).strip()
assert git('rev-parse','HEAD')==a.expected_source and not git('status','--porcelain')
keys=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')
env=dict(os.environ,**{k:'1' for k in keys})
env.update(PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(root/'tmp'),LD_LIBRARY_PATH=str(build/'geometry-lib'))
tests=source/'solvers/IncNavierStokesSolver/Tests/RANS'
fp=source/'docs/plans/cls-epawf-revised-evidence/m7/fingerprint-runtime.py'
preflight=json.loads(a.preflight.read_text())
assert preflight['status']=='passed' and preflight['full_runtime_source_inputs_unchanged']
assert preflight['source_commit']==a.expected_source
assert preflight['matched_initial_state_identity_all_quadratures']
initial={row['inner_stress']:row for row in preflight['cases'] if row['quadrature']=='base'}
assert set(initial)=={'Linear','SaturatedR2'}
assert initial['SaturatedR2']['finite_observed_joint_preflight']
matched_binding=None
if a.mode=='matched':
    if a.matched_receipt is not None or a.matched_assessment is not None:
        raise ValueError('Matched runs cannot consume fallback/retry evidence')
else:
    assert a.matched_receipt is not None
    # Future conditional recipe only. Recheck the original logs, complete
    # JSONL records, h0 and runtime instead of treating a missing final report
    # as a failed target. Positive convergence evidence blocks this branch.
    sys.dont_write_bytecode=True
    spec=importlib.util.spec_from_file_location('r2_stage_evidence',Path(__file__).with_name('assess-stage-evidence.py'))
    evidence=importlib.util.module_from_spec(spec);spec.loader.exec_module(evidence)
    matched_binding=evidence.check_frozen_scalar(a.matched_assessment,a.matched_receipt,
        expected_source=a.expected_source,preflight=a.preflight,original=a.original,
        representation=a.representation,runtime_pin=root/'runtime-fingerprint.json')
receipt=dict(schema=2,status='running',mode=a.mode,source_commit=a.expected_source,source_clean_start=True,
    runner_sha256=sha(__file__),preflight_path=str(a.preflight),preflight_sha256=sha(a.preflight),
    original_path=str(a.original),original_sha256=sha(a.original),representation=a.representation,
    target_dt=.0025,target_time=.0025,alpha_levels=[400.] if a.mode=='frozen-scalar' else [3200.,1600.,800.,400.],
    maximum_Newton_updates_per_level=12,maximum_seconds_per_native_case=14400,
    continuation='auxiliary mass parameter only; original target/history/time/forcing fixed; no Advance',
    ranks=1,numerical_threads={k:env[k] for k in keys},maximum_concurrent_workers=1,
    matrix_encoding=dict(dtype='IEEE754 binary64',byte_order=sys.byteorder,layout='row-major'),
    started_utc=datetime.now(timezone.utc).isoformat(),commands=[],cases=[])
if a.matched_receipt:receipt['matched_receipt']=dict(path=str(a.matched_receipt),sha256=sha(a.matched_receipt))
if matched_binding is not None:receipt['matched_evidence_gate']=matched_binding
def save():(out/'receipt.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
def command(name,argv,cwd,timeout,cpu=0):
    log=out/(name+'.log');started=time.monotonic()
    argv=['taskset','-c',str(cpu),*map(str,argv)]
    item=dict(name=name,argv=argv,cwd=str(cwd),timeout_seconds=timeout,cpu=cpu)
    receipt['commands'].append(item);save()
    with log.open('w') as stream:
        process=subprocess.Popen(argv,cwd=cwd,env=env,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
        item['pid']=process.pid;save();last=started
        while process.poll() is None:
            now=time.monotonic()
            if now-started>timeout:
                os.killpg(process.pid,signal.SIGTERM)
                try:process.wait(10)
                except subprocess.TimeoutExpired:os.killpg(process.pid,signal.SIGKILL);process.wait()
                item['exit_code']=124;break
            if now-last>=30:
                print(name+' alive '+str(round(now-started,1))+'s; '+str(log.stat().st_size)+' log bytes',flush=True);last=now
            time.sleep(1)
        else:item['exit_code']=process.returncode
    item.update(seconds=time.monotonic()-started,log_path=str(log),log_sha256=sha(log));save()
    print(name+' exit '+str(item['exit_code'])+' '+str(round(item['seconds'],3))+'s',flush=True)
    if item['exit_code']:print(log.read_text()[-3000:],flush=True)
    return item
save();before=False
try:
    item=command('fingerprint-before',['python3','-B',fp,'--root',root,'--output',out/'runtime-before.json'],out,240,49)
    assert item['exit_code']==0;before=True
    captured=json.loads((out/'runtime-before.json').read_text());built=json.loads((root/'runtime-fingerprint.json').read_text())
    for k in ('source_inputs','runtime_artifacts','runtime_bindings','runtime_environment','python_packages','tool_images'):
        assert captured[k]==built[k],k
    family='frozen-scalar' if a.mode=='frozen-scalar' else 'stage'
    item=command('prepare',['python3','-B',tests/'cls_epawf_m8_r2_cases.py','--original',a.original,
        '--output',out/'cases','--representation',a.representation,'--family',family],out,60,49)
    assert item['exit_code']==0
    plan=out/'cases/stage-plan.json';receipt['stage_plan_sha256']=sha(plan)
    for case in json.loads(plan.read_text())['cases']:
        xml=Path(case['xml']);name=case['name'];law=case['inner_stress']
        row=dict(name=name,inner_stress=law,xml_path=str(xml),xml_sha256=sha(xml),status='started')
        receipt['cases'].append(row);save()
        # Only actual predeclared levels inside this native entry point.
        item=command('native-'+law,[build/'IncNavierStokesSolver',xml],xml.parent,14400)
        raw_path=xml.with_name('case.evolution.json');trace=xml.with_name('case.feasible.jsonl')
        row.update(native_exit_code=item['exit_code'],native_seconds=item['seconds'])
        if raw_path.exists():
            raw=json.loads(raw_path.read_text())
            row.update(report_path=str(raw_path),report_sha256=sha(raw_path),
                actual_initial_state_sha256=vector_sha(raw['initial']))
            assert row['actual_initial_state_sha256']==initial[law]['initial_state_float64_LE_sha256']
            levels=raw.get('stage_test_levels',[])
            row.update(target_attempted=any(level.get('attempted') and level['alpha']==400. for level in levels),
                target_converged=raw.get('stage_test_target_converged',False),
                target_physical_range=raw.get('stage_test_target_physical_range',False),
                levels=[{k:level.get(k) for k in ('alpha','attempted','converged','iterations','failure',
                    'true_coupled_residual_converged','physical_range_passed','residual_verification_failure',
                    'qualification_audit_failure','scaled_block_norms')} for level in levels])
        if trace.exists():row.update(trace_path=str(trace),trace_sha256=sha(trace))
        row['status']='native_observed' if raw_path.exists() and trace.exists() else 'native_incomplete'
        # Analysis remains separate from obtaining a stage root; a helper
        # failure never silently removes the other matched target.
        save()
    receipt['status']='completed' if all(row['status']=='native_observed' for row in receipt['cases']) else 'incomplete'
except Exception as error:
    receipt.update(status='failed',failure=dict(type=type(error).__name__,message=str(error)))
finally:
    if before:
        item=command('fingerprint-after',['python3','-B',fp,'--root',root,'--output',out/'runtime-after.json',
            '--compare',out/'runtime-before.json'],out,240,49)
        receipt['runtime_unchanged']=item['exit_code']==0
    else:receipt['runtime_unchanged']=False
    receipt.update(source_commit_end=git('rev-parse','HEAD'),source_clean_end=not bool(git('status','--porcelain')),
        finished_utc=datetime.now(timezone.utc).isoformat())
    save();print(json.dumps(dict(status=receipt['status'],cases=receipt['cases'])),flush=True)
raise SystemExit(0 if receipt['status']=='completed' and receipt['runtime_unchanged'] else 1)
