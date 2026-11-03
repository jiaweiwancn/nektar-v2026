"""Audit the successful Original MMS solve and run the Neutral smoke afresh.

The preserved 002 attempt completed the solver and exact-field/traction audit,
then requested an unsupported full manufactured-forcing energy oracle. Reuse
only that successful solve under its unchanged full runtime pin. This check
does not claim a new independent full forcing/energy reconstruction.
"""
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
root = Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8/initialization')
source, build, old, out = root/'source', root/'build', root/'original-controls-002', root/'original-controls-003'
out.mkdir(exist_ok=False)
lock = (root/'validation.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_SH)
env = dict(os.environ)
env.update({k: '1' for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                           'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')})
env.update(PYTHONDONTWRITEBYTECODE='1', TMPDIR=str(root/'tmp'), LD_LIBRARY_PATH=str(build/'geometry-lib'))
os.environ.update(env); os.sched_setaffinity(0, {49})
tests = source/'solvers/IncNavierStokesSolver/Tests/RANS'; sys.path.insert(0, str(tests))
from cls_epawf_m8_mms import flat_session, verify
from cls_epawf_m8_native_oracle import verify as verify_energy
from validate_cls_epawf_m7 import audits
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
git = lambda *a: subprocess.check_output(['git', '-C', str(source), *a], text=True).strip()
receipt = dict(schema=1, status='running', source_commit=git('rev-parse', 'HEAD'),
               runner_sha256=sha(__file__), numerical_workers=1, numerical_threads=1,
               commands=[], started_utc=datetime.now(timezone.utc).isoformat(),
               limitation='Independent all-field and traction MMS plus selected accounting; '
               'full manufactured-forcing volume-energy reconstruction is not claimed')
def save(): (out/'receipt.json').write_text(json.dumps(receipt, indent=2, allow_nan=False)+'\n')
def run(name, argv):
    start=time.monotonic(); log=out/(name+'.log')
    with log.open('w') as stream: p=subprocess.run(argv, env=env, cwd=out, stdout=stream, stderr=subprocess.STDOUT, timeout=600)
    receipt['commands'].append(dict(name=name, argv=argv, exit_code=p.returncode, seconds=time.monotonic()-start, log_sha256=sha(log))); save()
    if p.returncode: raise RuntimeError(name+' failed')
fp=source/'docs/plans/cls-epawf-revised-evidence/m7/fingerprint-runtime.py'
save()
try:
    assert receipt['source_commit']=='ed34052fe4dd967e963c22272a0ec0888cf35b70' and not git('status','--porcelain')
    run('fingerprint-before',['python3','-B',str(fp),'--root',str(root),'--output',str(out/'runtime-before.json'),'--compare',str(old/'runtime-after.json')])
    previous=json.loads((old/'receipt.json').read_text())
    command=next(c for c in previous['commands'] if c['name']=='original-functions')
    assert command['exit_code']==0 and previous['source_and_runtime_unchanged']
    assert previous['failure']=='Full independent M8 energy oracle does not infer manufactured forcing'
    xml=flat_session(tests/'EquilibriumWallEvolution.xml',out/'regenerated.xml',variant='Original',beta=1.,gamma=0.,gamma_scale=0.,steps=1)
    assert sha(xml)==sha(old/'flat-functions/case.xml')
    path=old/'flat-functions/case.evolution.json'; raw=json.loads(path.read_text())
    assert not raw.get('failure') and raw['model_variant']=='Original' and raw['robust_parameters'][6]==0
    assert not raw['neutral_reference'] and not raw['mapped'] and not raw['upstream_transition']
    assert ':capped_initial_trace=' not in raw['signature']
    exact=verify(raw); exact['scope']='Original flat Functions MMS, gamma=0; independent exact fields/traction'
    receipt['original_functions']=dict(success=True, original_solver_command=command, original_report=str(path),
        report_sha256=sha(path), xml_sha256=sha(xml), previous_receipt_sha256=sha(old/'receipt.json'),
        manufactured=exact, native_audit=audits(raw,'original-functions'),
        selected_accounting=verify_energy(raw,full_volume=False),
        full_forcing_energy='NOT ESTABLISHED by this fixture; unsupported oracle request retained in 002')
    run('original-neutral-smoke',['ctest','--test-dir',str(build),'--parallel','1','-R','^RANSEquilibriumWallEvolutionSmoke$',
        '--output-on-failure','--output-junit',str(out/'ctest.xml')])
    cases=ET.parse(out/'ctest.xml').getroot().findall('testcase')
    assert len(cases)==1 and cases[0].get('name')=='RANSEquilibriumWallEvolutionSmoke' and not cases[0].findall('failure')
    receipt['original_neutral_smoke']=dict(success=True,junit_sha256=sha(out/'ctest.xml'))
    receipt['status']='passed'
except Exception as error: receipt.update(status='failed',failure=str(error))
finally:
    try:
        run('fingerprint-after',['python3','-B',str(fp),'--root',str(root),'--output',str(out/'runtime-after.json'),'--compare',str(out/'runtime-before.json')])
        assert git('rev-parse','HEAD')==receipt['source_commit'] and not git('status','--porcelain')
        receipt['source_and_runtime_unchanged']=True
    except Exception as error: receipt.update(status='failed',fingerprint_failure=str(error))
    receipt['finished_utc']=datetime.now(timezone.utc).isoformat();save()
print(json.dumps({k:receipt[k] for k in ('status','failure','source_and_runtime_unchanged') if k in receipt}))
raise SystemExit(0 if receipt['status']=='passed' else 1)
