"""One unchanged-policy diagnostic reproduction and same-stage fresh retry.

This is unqualified from t=0 because the archived initial state already loses
inner realizability and leaves the local outer certificates. It is not a
qualification trajectory. The counterfactual retry never advances IMEX history.
"""
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
import xml.etree.ElementTree as ET

base=Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised')
root=base/'m8-r1';source,build=root/'source',root/'build';out=root/'reproduction-001'
out.mkdir(exist_ok=False);case=out/'case';case.mkdir()
lock=(root/'validation.lock').open('a');fcntl.flock(lock,fcntl.LOCK_SH)
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def git(*a):return subprocess.check_output(['git','-C',str(source),*a],text=True).strip()
old=base/'m8/flow-001/qualification/control-coarse';archive=root/'archive-attempt-003/archive-ledger.json'
assert json.loads(archive.read_text())['success'],'Complete Phase A is required before a new trajectory'
revision=git('rev-parse','HEAD');assert revision=='426d48cd44f9f4af017244f6be4c1a6f62d36c01' and not git('status','--porcelain')
tree=ET.parse(old/'case.xml');conditions=tree.getroot().find('CONDITIONS');infos=conditions.find('SOLVERINFO')
flags=dict(RANSWallSolveDiagnostics='Failure',RANSWallSolveDiagnosticReplay='Fresh',
           RANSWallSolveDiagnosticIterations='12',RANSWallSolveDiagnosticFile='case.nonlinear.json')
for key,value in flags.items():
    assert not any(n.get('PROPERTY')==key for n in infos)
    ET.SubElement(infos,'I',PROPERTY=key,VALUE=value)
xml=case/'case.xml';tree.write(xml,encoding='utf-8',xml_declaration=True)
def semantic(node):
    children=[n for n in node if not (n.tag=='I' and n.get('PROPERTY') in flags)]
    return [node.tag,sorted(node.attrib.items()),(node.text or '').strip(),[semantic(n) for n in children]]
assert semantic(ET.parse(old/'case.xml').getroot())==semantic(ET.parse(xml).getroot())
env=dict(os.environ);env.update({key:'1' for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
    'BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')})
env.update(PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(root/'tmp'),LD_LIBRARY_PATH=str(build/'geometry-lib'))
receipt=dict(schema=1,status='running',source_commit=revision,runner_sha256=sha(__file__),
    source_clean_start=True,old_XML_sha256=sha(old/'case.xml'),XML_sha256=sha(xml),
    original_report_sha256=sha(old/'case.evolution.json'),complete_archive_ledger_sha256=sha(archive),
    physical_session_semantics_unchanged=True,diagnostic_solver_info=flags,
    executable=str(build/'IncNavierStokesSolver'),executable_sha256=sha(build/'IncNavierStokesSolver'),
    connection='tailscale:100.105.247.110',ranks=1,numerical_threads=1,cpu=0,
    diagnostic_unqualified=True,qualification=False,commands=[],
    run_budget=dict(maximum_native_evolution_and_same_stage_retries=6,
        allocated_here=[dict(slot=1,kind='original unforced control, unchanged solver policy, diagnostic from t0'),
                        dict(slot=2,kind='counterfactual same failed stage; fresh Jacobian every iteration, max12, no history advancement')],
        remaining_slots=[3,4,5,6],reproduction_duration='Original XML NumSteps and dt retained; expected failure at .060; hard wall timeout3600s'),
    started_utc=datetime.now(timezone.utc).isoformat())
def save():(out/'receipt.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
def run(name,argv,cwd=out,timeout=240):
    start=time.monotonic();log=out/(name+'.log')
    try:
        with log.open('w') as stream:code=subprocess.run(argv,cwd=cwd,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=timeout).returncode
    except subprocess.TimeoutExpired:code=124
    receipt['commands'].append(dict(name=name,argv=argv,cwd=str(cwd),exit_code=code,seconds=time.monotonic()-start,
        timeout=timeout,log=str(log),log_sha256=sha(log)));save();return code
fp=source/'docs/plans/cls-epawf-revised-evidence/m7/fingerprint-runtime.py'
save()
try:
    assert run('fingerprint-before',['python3','-B',str(fp),'--root',str(root),'--output',str(out/'runtime-before.json')])==0
    pin=json.loads((root/'build-001/runtime-fingerprint.json').read_text());before=json.loads((out/'runtime-before.json').read_text())
    for key in ('source_inputs','runtime_artifacts','runtime_bindings','python_packages','tool_images'):
        assert pin[key]==before[key],key
    receipt['environment_changes_from_build']={key:[pin['runtime_environment'].get(key),v] for key,v in before['runtime_environment'].items() if pin['runtime_environment'].get(key)!=v}
    save();code=run('solver',['taskset','-c','0',str(build/'IncNavierStokesSolver'),str(xml)],cwd=case,timeout=3600)
    report=case/'case.evolution.json';trace=case/'case.nonlinear.json'
    receipt.update(solver_exit_code=code,report_sha256=sha(report) if report.exists() else None,
        trace_sha256=sha(trace) if trace.exists() else None)
    if not report.exists() or not trace.exists():raise RuntimeError('Expected reports missing; do not retry automatically')
    raw=json.loads(report.read_text());diagnostic=json.loads(trace.read_text())
    expected='evolution stage rejected every admissible Newton trial'
    receipt.update(original_failure=raw.get('failure'),failed_stage_time=diagnostic['time'],
        failed_stage_alpha=diagnostic['alpha'],last_accepted_time=max(raw['audit_times']),
        accepted_state_count=len(raw['accepted_states']),same_stage_diagnostic=diagnostic.get('fresh_jacobian_retry'))
    # Keep large retry iterates only in the trace; retain its terminal metrics here.
    value=receipt['same_stage_diagnostic']
    if value:receipt['same_stage_diagnostic']={k:value[k] for k in ('converged','terminal','maximum_iterations','history_advanced','scaled_block_norms') if k in value}
    assert code!=0 and expected in raw.get('failure','') and expected in diagnostic['failure']
    assert abs(diagnostic['time']-.06)<1e-12 and abs(max(raw['audit_times'])-.055)<1e-12
    receipt['status']='original_failure_reproduced'
except Exception as error:receipt.update(status='failed_or_unexpected_result',error=str(error))
finally:
    after=run('fingerprint-after',['python3','-B',str(fp),'--root',str(root),'--output',str(out/'runtime-after.json'),'--compare',str(out/'runtime-before.json')])
    receipt.update(runtime_unchanged=after==0,source_clean_end=not bool(git('status','--porcelain')),source_commit_end=git('rev-parse','HEAD'),
        finished_utc=datetime.now(timezone.utc).isoformat());save()
print(json.dumps({key:receipt[key] for key in ('status','failed_stage_time','last_accepted_time','same_stage_diagnostic','error') if key in receipt}),flush=True)
raise SystemExit(0 if receipt['status']=='original_failure_reproduced' and receipt['runtime_unchanged'] and receipt['source_clean_end'] else 1)
