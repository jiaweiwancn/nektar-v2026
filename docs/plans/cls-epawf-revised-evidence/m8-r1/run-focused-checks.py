"""Focused existing unit probes and two zero-step range-guard checks.

No additional PDE evolution trajectory is permitted by this runner. Existing
time/patch/evolution probes exercise toy IMEX history, field/constraint loads,
and dense algebra; both solver cases request NumSteps=0 and must stop at t=0.
"""
from datetime import datetime, timezone
import fcntl, hashlib, json, os, subprocess, time
from pathlib import Path
import xml.etree.ElementTree as ET

base = Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised')
r1 = base/'m8-r1'; root = r1/'controlled'
source, build = root/'source', root/'build'
out = root/'focused-checks-001'; out.mkdir(exist_ok=False)
lock = (root/'validation.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_SH)
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
git = lambda *a: subprocess.check_output(['git', '-C', str(source), *a], text=True).strip()
assert git('rev-parse', 'HEAD') == '360cbb7ac728b575cdf9af21e6d0963a1d88acb9'
assert not git('status', '--porcelain')
env = dict(os.environ)
keys = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
        'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')
env.update({key: '1' for key in keys})
env.update(PYTHONDONTWRITEBYTECODE='1', TMPDIR=str(root/'tmp'), LD_LIBRARY_PATH=str(build/'geometry-lib'))
receipt = dict(schema=1, status='running', source_commit=git('rev-parse', 'HEAD'),
    source_clean_start=True, runner_sha256=sha(__file__), ranks=1, cpu=1,
    numerical_threads={key: env[key] for key in keys}, commands=[], zero_step_checks=[],
    additional_PDE_trajectories=0, started_utc=datetime.now(timezone.utc).isoformat())
def save(): (out/'receipt.json').write_text(json.dumps(receipt, indent=2, allow_nan=False)+'\n')
def run(name, argv, cwd=out, timeout=300):
    start = time.monotonic(); log = out/(name+'.log')
    try:
        with log.open('w') as stream:
            code = subprocess.run(argv, cwd=cwd, env=env, stdout=stream,
                stderr=subprocess.STDOUT, timeout=timeout).returncode
    except subprocess.TimeoutExpired: code = 124
    receipt['commands'].append(dict(name=name, argv=argv, cwd=str(cwd), exit_code=code,
        seconds=time.monotonic()-start, log_sha256=sha(log), timeout=timeout))
    save(); return code
fp = source/'docs/plans/cls-epawf-revised-evidence/m7/fingerprint-runtime.py'
save()
try:
    assert run('fingerprint-before', ['python3', '-B', str(fp), '--root', str(root),
        '--output', str(out/'runtime-before.json')]) == 0
    a = json.loads((out/'runtime-before.json').read_text())
    b = json.loads((root/'build-002/runtime-fingerprint.json').read_text())
    for key in ('source_inputs', 'runtime_artifacts', 'runtime_bindings', 'python_packages', 'tool_images'):
        assert a[key] == b[key], key
    assert run('native-probes', ['taskset', '-c', '1', 'ctest', '--test-dir', str(build),
        '-R', '^RANSEquilibrium(Time|Evolution|Patch)Native$', '--parallel', '1',
        '--output-on-failure'], timeout=900) == 0
    archived = base/'m8/flow-001/qualification/control-coarse/case.xml'
    resolved = root/'controls-001/cases/resolved-dt005/case.xml'
    for name, template, guard, reason, size in (
        ('original-physical-stop', archived, 'R1DiagnosticV1', 'covariance:', 425),
        ('resolved-local-stop', resolved, 'R1LocalV1', 'outer_gradient_box_exit', 1571)):
        directory = out/name; directory.mkdir(); tree = ET.parse(template); xml = directory/'case.xml'
        p = next(node for node in tree.getroot().findall('CONDITIONS/PARAMETERS/P')
                 if node.text.split('=', 1)[0].strip().upper() == 'NUMSTEPS')
        p.text = 'NumSteps = 0'
        info = tree.getroot().find('CONDITIONS/SOLVERINFO')
        for key, value in dict(RANSWallPhysicalRangeGuard=guard,
            RANSWallSolveDiagnostics='Failure', RANSWallSolveDiagnosticReplay='Off',
            RANSWallSolveDiagnosticFile='case.nonlinear.json').items():
            nodes = [n for n in info if n.get('PROPERTY', '').upper() == key.upper()]
            assert len(nodes) <= 1
            if nodes: nodes[0].set('VALUE', value)
            else: ET.SubElement(info, 'I', PROPERTY=key, VALUE=value)
        tree.write(xml, encoding='utf-8', xml_declaration=True)
        code = run(name, ['taskset', '-c', '1', str(build/'IncNavierStokesSolver'), str(xml)], directory)
        report = directory/'case.evolution.json'; raw = json.loads(report.read_text())
        events = raw['physical_range_events']; event = events[0]
        assert code != 0 and len(events) == 1 and event['time'] == 0
        assert event['stopped_before_history_advance'] and any(reason in v for v in event['reasons'])
        assert len(event['state']) == size and len(event['history']) == size
        assert not raw['accepted_states'] and not raw['audit_times']
        receipt['zero_step_checks'].append(dict(name=name, guard=guard, XML_sha256=sha(xml),
            report_sha256=sha(report), exit_code=code, state_unknowns=size,
            time=event['time'], reasons=event['reasons'], stopped_before_history_advance=True))
        save()
    receipt['status'] = 'passed'
except Exception as error:
    receipt.update(status='failed', error=repr(error))
finally:
    code = run('fingerprint-after', ['python3', '-B', str(fp), '--root', str(root),
        '--output', str(out/'runtime-after.json'), '--compare', str(out/'runtime-before.json')])
    receipt.update(runtime_unchanged=code == 0, source_clean_end=not bool(git('status', '--porcelain')),
        source_commit_end=git('rev-parse', 'HEAD'), finished_utc=datetime.now(timezone.utc).isoformat())
    save()
print(json.dumps({key: receipt[key] for key in ('status', 'error', 'zero_step_checks') if key in receipt}), flush=True)
raise SystemExit(0 if receipt['status'] == 'passed' and receipt['runtime_unchanged'] else 1)
