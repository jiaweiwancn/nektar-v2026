"""Export all offending archived inner points; never launch a solver."""
import hashlib, json, os, subprocess, time
from pathlib import Path

r1 = Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8-r1')
tools = r1/'point-tools-001'; out = r1/'archive-points-runner-001'; out.mkdir(exist_ok=False)
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
driver = tools/'cls_epawf_m8_r1_points.py'
assert sha(driver) == '87fbb11cf86163035499bd55da355ad568a63854439b086745144c169e453501'
argv = ['taskset', '-c', '50', 'python3', '-B', str(driver), '--source', str(r1/'archive-source'),
    '--archive', str(r1.parent/'m8/flow-001/qualification/control-coarse'),
    '--ledger', str(r1/'archive-attempt-003/archive-ledger.json'),
    '--reconstruction-driver', str(r1/'cls_epawf_m8_r1_archive.py'),
    '--output', str(r1/'archive-points-001')]
env = dict(os.environ)
env.update({key: '1' for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
    'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')})
env.update(PYTHONDONTWRITEBYTECODE='1', TMPDIR=str(r1/'tmp'))
start = time.monotonic(); log = out/'run.log'
with log.open('w') as stream:
    code = subprocess.run(argv, cwd=out, env=env, stdout=stream,
                          stderr=subprocess.STDOUT, timeout=900).returncode
receipt = dict(schema=1, source='archive-only postprocessing; no new trajectory',
    runner_sha256=sha(__file__), driver_sha256=sha(driver), command=argv,
    cwd=str(out), exit_code=code, seconds=time.monotonic()-start, log_sha256=sha(log))
(out/'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
print(log.read_text()[-2000:]); raise SystemExit(code)
