"""One justified repeated prefix for finite-domain Jacobian diagnostics.

The first unchanged-policy reproduction showed that a fresh original central
column leaves the normalization domain. This comparison keeps the native
policy/residual unchanged and changes only its post-failure counterfactual:
scaled admissible central/one-sided columns, with the same trial fractions.
It is unqualified throughout and consumes conservative budget slots 3/4.
"""
from datetime import datetime, timezone
import fcntl,hashlib,json,os,subprocess,time
from pathlib import Path
import xml.etree.ElementTree as ET
base=Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised');r1=base/'m8-r1';root=r1/'controlled'
source,build=root/'source',root/'build';out=root/'reproduction-001';out.mkdir(exist_ok=False);case=out/'case';case.mkdir()
lock=(root/'validation.lock').open('a');fcntl.flock(lock,fcntl.LOCK_SH)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
git=lambda *a:subprocess.check_output(['git','-C',str(source),*a],text=True).strip()
assert git('rev-parse','HEAD')=='360cbb7ac728b575cdf9af21e6d0963a1d88acb9' and not git('status','--porcelain')
assert json.loads((r1/'archive-attempt-003/archive-ledger.json').read_text())['success']
previous=json.loads((r1/'reproduction-001/receipt.json').read_text())
assert previous['status']=='original_failure_reproduced' and previous['same_stage_diagnostic']['terminal']=='diagnostic_evaluation_or_linear_failure'
old=base/'m8/flow-001/qualification/control-coarse';tree=ET.parse(old/'case.xml');infos=tree.getroot().find('CONDITIONS/SOLVERINFO')
flags=dict(RANSWallSolveDiagnostics='Failure',RANSWallSolveDiagnosticReplay='AdmissibleFresh',
           RANSWallSolveDiagnosticIterations='12',RANSWallSolveDiagnosticFile='case.nonlinear.json')
for key,value in flags.items():ET.SubElement(infos,'I',PROPERTY=key,VALUE=value)
xml=case/'case.xml';tree.write(xml,encoding='utf-8',xml_declaration=True)
def semantic(node):return [node.tag,sorted(node.attrib.items()),(node.text or '').strip(),
    [semantic(n) for n in node if not (n.tag=='I' and n.get('PROPERTY') in flags)]]
assert semantic(ET.parse(old/'case.xml').getroot())==semantic(ET.parse(xml).getroot())
env=dict(os.environ);env.update({key:'1' for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
    'BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')})
env.update(PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(root/'tmp'),LD_LIBRARY_PATH=str(build/'geometry-lib'))
receipt=dict(schema=1,status='running',source_commit=git('rev-parse','HEAD'),source_clean_start=True,
    runner_sha256=sha(__file__),executable_sha256=sha(build/'IncNavierStokesSolver'),
    XML_sha256=sha(xml),old_XML_sha256=sha(old/'case.xml'),original_report_sha256=sha(old/'case.evolution.json'),
    previous_reproduction_receipt_sha256=sha(r1/'reproduction-001/receipt.json'),physical_semantics_unchanged=True,
    diagnostic_unqualified=True,physical_guard='Off for exact numerical-failure reproduction; already physically/range-invalid at t0',
    diagnostic_solver_info=flags,ranks=1,numerical_threads=1,cpu=0,commands=[],
    bounded_matrix_slots=[3,4],remaining_native_slots=[5,6],
    justification='Original fresh central assembly was outside finite evaluation domain; independently scaled/admissible differences are required, without changing original line-search minimum or tolerances.',
    started_utc=datetime.now(timezone.utc).isoformat())
def save():(out/'receipt.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
def run(name,argv,cwd=out,timeout=240):
    start=time.monotonic();log=out/(name+'.log')
    try:
        with log.open('w') as stream:code=subprocess.run(argv,cwd=cwd,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=timeout).returncode
    except subprocess.TimeoutExpired:code=124
    receipt['commands'].append(dict(name=name,argv=argv,cwd=str(cwd),exit_code=code,seconds=time.monotonic()-start,
        log_sha256=sha(log),timeout=timeout));save();return code
fp=source/'docs/plans/cls-epawf-revised-evidence/m7/fingerprint-runtime.py';save()
try:
    assert run('fingerprint-before',['python3','-B',str(fp),'--root',str(root),'--output',str(out/'runtime-before.json')])==0
    a=json.loads((out/'runtime-before.json').read_text());b=json.loads((root/'build-002/runtime-fingerprint.json').read_text())
    for key in ('source_inputs','runtime_artifacts','runtime_bindings','python_packages','tool_images'):assert a[key]==b[key],key
    receipt['environment_changes_from_build']={key:[b['runtime_environment'].get(key),v] for key,v in a['runtime_environment'].items() if b['runtime_environment'].get(key)!=v}
    code=run('solver',['taskset','-c','0',str(build/'IncNavierStokesSolver'),str(xml)],case,3600)
    path=case/'case.evolution.json';trace=case/'case.nonlinear.json'
    receipt.update(solver_exit_code=code,report_sha256=sha(path) if path.exists() else None,trace_sha256=sha(trace) if trace.exists() else None)
    assert path.exists() and trace.exists(),'Missing preserved diagnostic result'
    raw=json.loads(path.read_text());detail=json.loads(trace.read_text())
    assert code!=0 and 'evolution stage rejected every admissible Newton trial' in raw.get('failure','')
    assert abs(detail['time']-.06)<1e-12 and abs(max(raw['audit_times'])-.055)<1e-12
    receipt.update(status='original_failure_reproduced',failed_stage_time=detail['time'],last_accepted_time=max(raw['audit_times']),
        fresh_matrix_available=detail['fresh_jacobian'] is not None,directional_ladder_available=detail['directional_checks'] is not None,
        diagnostic_retry={key:detail['fresh_jacobian_retry'][key] for key in ('terminal','converged','history_advanced','maximum_iterations','scaled_block_norms') if key in detail['fresh_jacobian_retry']})
except Exception as error:receipt.update(status='failed_or_unexpected_result',error=str(error))
finally:
    result=run('fingerprint-after',['python3','-B',str(fp),'--root',str(root),'--output',str(out/'runtime-after.json'),'--compare',str(out/'runtime-before.json')])
    receipt.update(runtime_unchanged=result==0,source_clean_end=not bool(git('status','--porcelain')),source_commit_end=git('rev-parse','HEAD'),finished_utc=datetime.now(timezone.utc).isoformat());save()
print(json.dumps({key:receipt[key] for key in ('status','failed_stage_time','fresh_matrix_available','directional_ladder_available','diagnostic_retry','error') if key in receipt}),flush=True)
raise SystemExit(0 if receipt['status']=='original_failure_reproduced' and receipt['runtime_unchanged'] and receipt['source_clean_end'] else 1)
