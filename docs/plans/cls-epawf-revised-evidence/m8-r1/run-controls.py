"""The final two sequential short M8-R1 native evolution allocations.

Both use a fresh IMEX initialization. The frozen local certificate has already
been exited at t=0: these are explicitly unqualified diagnostics. The guard
records that exit, then stops before advancing history at the first audited
scalar/covariance/lambda violation. No production solver policy is changed.
"""
from datetime import datetime, timezone
import fcntl, hashlib, json, os, subprocess, sys, time
from pathlib import Path

base = Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised')
r1 = base/'m8-r1'
root = r1/'controlled'
source, build = root/'source', root/'build'
out = root/'controls-001'
out.mkdir(exist_ok=False)
lock = (root/'validation.lock').open('a')
fcntl.flock(lock, fcntl.LOCK_SH)
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
git = lambda *a: subprocess.check_output(['git', '-C', str(source), *a], text=True).strip()
assert git('rev-parse', 'HEAD') == '360cbb7ac728b575cdf9af21e6d0963a1d88acb9'
assert not git('status', '--porcelain')
prior = json.loads((root/'reproduction-001/receipt.json').read_text())
assert prior['status'] == 'original_failure_reproduced' and prior['runtime_unchanged']
tools = Path(__file__).parent
sys.path[:0] = [str(tools), str(r1/'archive-source/solvers/IncNavierStokesSolver/Tests/RANS')]
from cls_epawf_m8_r1_cases import prepare_trajectories
plan = prepare_trajectories(base/'m8/flow-001/qualification/control-coarse/case.xml', out/'cases')
env = dict(os.environ)
threads = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
           'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')
env.update({key: '1' for key in threads})
env.update(PYTHONDONTWRITEBYTECODE='1', TMPDIR=str(root/'tmp'),
           LD_LIBRARY_PATH=str(build/'geometry-lib'))
receipt = dict(schema=1, status='running', source_commit=git('rev-parse', 'HEAD'),
    source_clean_start=True, runner_sha256=sha(__file__),
    generator_sha256=sha(tools/'cls_epawf_m8_r1_cases.py'),
    executable_sha256=sha(build/'IncNavierStokesSolver'),
    plan_sha256=sha(out/'cases/trajectory-plan.json'),
    ranks=1, numerical_threads={key: env[key] for key in threads}, cpu=0,
    bounded_matrix_slots=[5, 6], additional_evolution_allocations_permitted=0,
    controls_sequential=True, diagnostic_unqualified_from_t0=True,
    commands=[], results=[], started_utc=datetime.now(timezone.utc).isoformat())
def save():
    (out/'receipt.json').write_text(json.dumps(receipt, indent=2, allow_nan=False)+'\n')
def run(name, argv, cwd=out, timeout=240):
    start = time.monotonic()
    log = out/(name+'.log')
    try:
        with log.open('w') as stream:
            code = subprocess.run(argv, cwd=cwd, env=env, stdout=stream,
                                  stderr=subprocess.STDOUT, timeout=timeout).returncode
    except subprocess.TimeoutExpired:
        code = 124
    receipt['commands'].append(dict(name=name, argv=argv, cwd=str(cwd),
        exit_code=code, seconds=time.monotonic()-start, log_sha256=sha(log), timeout=timeout))
    save()
    return code
fp = source/'docs/plans/cls-epawf-revised-evidence/m7/fingerprint-runtime.py'
save()
try:
    assert run('fingerprint-before', ['python3', '-B', str(fp), '--root', str(root),
                                    '--output', str(out/'runtime-before.json')]) == 0
    a = json.loads((out/'runtime-before.json').read_text())
    b = json.loads((root/'build-002/runtime-fingerprint.json').read_text())
    for key in ('source_inputs', 'runtime_artifacts', 'runtime_bindings', 'python_packages', 'tool_images'):
        assert a[key] == b[key], key
    receipt['environment_changes_from_build'] = {key: [b['runtime_environment'].get(key), value]
        for key, value in a['runtime_environment'].items() if b['runtime_environment'].get(key) != value}
    for item in plan['cases']:
        xml = Path(item['xml'])
        code = run(item['name'], ['taskset', '-c', '0', str(build/'IncNavierStokesSolver'), str(xml)],
                   xml.parent, 3600)
        result = dict(name=item['name'], XML_sha256=sha(xml), exit_code=code, fresh_history=True)
        report = xml.parent/'case.evolution.json'
        if report.exists():
            raw = json.loads(report.read_text())
            result.update(report_sha256=sha(report), last_audit_time=max(raw.get('audit_times', [0])),
                          success=raw.get('success'), failure=raw.get('failure'),
                          physical_range_events=raw.get('physical_range_events', []))
        trace = xml.parent/'case.nonlinear.json'
        result['trace_sha256'] = sha(trace) if trace.exists() else None
        receipt['results'].append(result)
        save()
    receipt['status'] = 'bounded_controls_completed'
except Exception as error:
    receipt.update(status='failed_or_unexpected_result', error=str(error))
finally:
    result = run('fingerprint-after', ['python3', '-B', str(fp), '--root', str(root),
        '--output', str(out/'runtime-after.json'), '--compare', str(out/'runtime-before.json')])
    receipt.update(runtime_unchanged=result == 0, source_clean_end=not bool(git('status', '--porcelain')),
        source_commit_end=git('rev-parse', 'HEAD'), finished_utc=datetime.now(timezone.utc).isoformat())
    save()
print(json.dumps(dict(status=receipt['status'],
    results=[{key: value for key, value in row.items() if key != 'physical_range_events'}
             for row in receipt['results']])), flush=True)
raise SystemExit(0 if receipt['status'] == 'bounded_controls_completed' and receipt['runtime_unchanged'] else 1)
