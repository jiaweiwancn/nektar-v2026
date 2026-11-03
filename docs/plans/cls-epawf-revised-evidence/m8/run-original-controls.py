"""Check both unchanged Original startup routes in the corrected solver."""
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
source, build, out = root/'source', root/'build', root/'original-controls-002'
out.mkdir(exist_ok=False)
lock = (root/'validation.lock').open('a')
fcntl.flock(lock, fcntl.LOCK_SH)
env = dict(os.environ)
env.update({name: '1' for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                                  'MKL_NUM_THREADS', 'BLIS_NUM_THREADS',
                                  'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')})
env.update(PYTHONDONTWRITEBYTECODE='1', TMPDIR=str(root/'tmp'),
           LD_LIBRARY_PATH=str(build/'geometry-lib'))
os.environ.update(env)
os.sched_setaffinity(0, {49})
tests = source/'solvers/IncNavierStokesSolver/Tests/RANS'
sys.path.insert(0, str(tests))
from cls_epawf_m8_mms import flat_session, verify
from cls_epawf_m8_native_oracle import verify as verify_energy
from validate_cls_epawf_m7 import audits


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def git(*args):
    return subprocess.check_output(['git', '-C', str(source), *args], text=True).strip()


receipt = dict(schema=1, status='running', source_commit=git('rev-parse', 'HEAD'),
               runner_sha256=sha(__file__), commands=[], numerical_workers=1,
               numerical_threads=1, connection='tailscale:100.105.247.110',
               scope='Focused Original/Original-source regression on corrected startup solver',
               started_utc=datetime.now(timezone.utc).isoformat())
assert receipt['source_commit'] == 'ed34052fe4dd967e963c22272a0ec0888cf35b70'
assert not git('status', '--porcelain')


def save():
    (out/'receipt.json').write_text(json.dumps(receipt, indent=2, allow_nan=False)+'\n')


def command(name, argv, cwd=out, timeout=600):
    started = time.monotonic()
    log = out/(name+'.log')
    with log.open('w') as stream:
        result = subprocess.run(argv, cwd=cwd, env=env, stdout=stream,
                                stderr=subprocess.STDOUT, timeout=timeout)
    receipt['commands'].append(dict(name=name, argv=list(map(str, argv)), cwd=str(cwd),
        exit_code=result.returncode, seconds=time.monotonic()-started,
        log_sha256=sha(log), timeout=timeout))
    save()
    if result.returncode:
        raise RuntimeError(name+' failed; inspect preserved log')


fingerprint = source/'docs/plans/cls-epawf-revised-evidence/m7/fingerprint-runtime.py'
save()
try:
    command('fingerprint-before', ['python3', '-B', str(fingerprint), '--root', str(root),
        '--output', str(out/'runtime-before.json')])
    before = json.loads((out/'runtime-before.json').read_text())
    built = json.loads((root/'build-002/runtime-fingerprint.json').read_text())
    for key in ('source_inputs', 'runtime_artifacts', 'runtime_bindings', 'python_packages', 'tool_images'):
        assert before[key] == built[key], key
    receipt['build_environment_differences'] = {key: [built['runtime_environment'].get(key), value]
        for key, value in before['runtime_environment'].items()
        if value != built['runtime_environment'].get(key)}
    directory = out/'flat-functions'
    directory.mkdir()
    xml = flat_session(tests/'EquilibriumWallEvolution.xml', directory/'case.xml',
                       variant='Original', beta=1., gamma=0., gamma_scale=0., steps=1)
    command('original-functions', [str(build/'IncNavierStokesSolver'), str(xml)], directory)
    path = directory/'case.evolution.json'
    raw = json.loads(path.read_text())
    assert not raw.get('failure') and raw['model_variant'] == 'Original'
    assert not raw['neutral_reference'] and not raw['mapped'] and not raw['upstream_transition']
    assert ':capped_initial_trace=' not in raw['signature'] and raw['robust_parameters'][6] == 0
    assert raw['robust_parameters'][0] == 1 and raw['robust_parameters'][1] == raw['robust_parameters'][5] == 0
    exact = verify(raw)
    exact['scope'] = 'Original flat Functions MMS, zero gamma; no new model selection'
    receipt['original_functions'] = dict(success=True, xml_sha256=sha(xml),
        report_sha256=sha(path), manufactured=exact, native_audit=audits(raw, 'original-functions'),
        energy=verify_energy(raw, full_volume=True))
    command('original-neutral-smoke', ['ctest', '--test-dir', str(build), '--parallel', '1',
        '-R', '^RANSEquilibriumWallEvolutionSmoke$', '--output-on-failure',
        '--output-junit', str(out/'ctest.xml')])
    cases = ET.parse(out/'ctest.xml').getroot().findall('testcase')
    assert len(cases) == 1 and not cases[0].findall('failure') and cases[0].get('status') == 'run'
    receipt['original_neutral_smoke'] = dict(success=True, junit_sha256=sha(out/'ctest.xml'))
    receipt['status'] = 'passed'
except Exception as error:
    receipt.update(status='failed', failure=str(error))
finally:
    try:
        command('fingerprint-after', ['python3', '-B', str(fingerprint), '--root', str(root),
            '--output', str(out/'runtime-after.json'), '--compare', str(out/'runtime-before.json')])
        assert git('rev-parse', 'HEAD') == receipt['source_commit'] and not git('status', '--porcelain')
        receipt['source_and_runtime_unchanged'] = True
    except Exception as error:
        receipt.update(status='failed', fingerprint_failure=str(error))
    receipt['finished_utc'] = datetime.now(timezone.utc).isoformat()
    save()
print(json.dumps({k: receipt[k] for k in ('status', 'source_commit', 'source_and_runtime_unchanged') if k in receipt}))
raise SystemExit(0 if receipt['status'] == 'passed' else 1)
