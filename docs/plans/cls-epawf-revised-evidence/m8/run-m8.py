"""Source-bound Linux M8 runner. All generated data stay outside the checkout."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
from datetime import datetime, timezone

p = argparse.ArgumentParser()
p.add_argument('--bundle', required=True)
p.add_argument('--ref', required=True)
p.add_argument('--attempt', required=True)
p.add_argument('--mode', choices=['build', 'kernel', 'pilot', 'native', 'transverse', 'regression'], default='build')
p.add_argument('--gamma-scale', type=float, default=0)
p.add_argument('--variant', choices=['CubicScaled','Regularized'], default='CubicScaled')
p.add_argument('--beta', type=float, default=.4)
a = p.parse_args()
root = Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8')
main = Path('/media/jay/新加卷/nektar-v2026')
source, build = root/'source', root/'build'
root.mkdir(exist_ok=True)
lock = (root/'validation.lock').open('a')
fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
def git(path, *args):
    return subprocess.check_output(['git', '-c', 'gc.auto=0', '-c', 'pack.threads=1', '-C', str(path), *args], text=True).strip()
assert not git(main, 'status', '--porcelain')
main_revision = git(main, 'rev-parse', 'HEAD')
if source.exists():
    assert not git(source, 'status', '--porcelain')
remote_ref = 'refs/remotes/cls-epawf/m8-'+a.ref
subprocess.run(['git', '-c', 'gc.auto=0', '-c', 'pack.threads=1', '-C', str(main), 'fetch',
                str(root/'transfers'/a.bundle), 'refs/heads/codex/cls-epawf-revised:'+remote_ref], check=True)
revision = git(main, 'rev-parse', remote_ref)
if source.exists():
    subprocess.run(['git', '-C', str(source), 'checkout', '--detach', revision], check=True)
else:
    subprocess.run(['git', '-C', str(main), 'worktree', 'add', '--detach', str(source), revision], check=True)
out = root/a.attempt
out.mkdir(exist_ok=False)
(root/'tmp').mkdir(exist_ok=True)
env = os.environ.copy()
threads = {k: '1' for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS')}
env.update(threads)
env['PYTHONDONTWRITEBYTECODE'] = '1'
env['TMPDIR'] = str(root/'tmp')
env['LD_LIBRARY_PATH'] = str(build/'geometry-lib')+':'+env.get('LD_LIBRARY_PATH', '')
receipt = dict(schema=1, status='running', source_commit=revision, main_commit=main_revision,
               connection='tailscale:100.105.247.110', compiler_workers=16,
               numerical_threads=threads, mode=a.mode, commands=[],
               maximum_test_ranks=12 if a.mode=='regression' else 2 if a.mode=='native' else 1,
               runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               started_utc=datetime.now(timezone.utc).isoformat())
def save():
    (out/'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
save()
commands = [('configure', ['cmake', '-S', str(source/'solvers/IncNavierStokesSolver/Tests/RANS'),
    '-B', str(build), '-DRANS_BUILD_SOLVER=ON', '-DCMAKE_BUILD_TYPE=Release',
    '-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake'], 240),
    ('build', ['cmake', '--build', str(build), '--parallel', '16'], 1800)]
fingerprint=source/'docs/plans/cls-epawf-revised-evidence/m7/fingerprint-runtime.py'
commands += [('fingerprint-before',['python3','-B',str(fingerprint),'--root',str(root),
    '--output',str(out/'runtime-before.json')],240)]
tests = source/'solvers/IncNavierStokesSolver/Tests/RANS'
if a.mode == 'kernel':
    commands += [('kernel', ['python3', '-B', str(tests/'validate_cls_epawf_m8_kernel.py'),
        '--probe', str(build/'RANSRobustVariantProbe'), '--output', str(out/'qualification')], 3600)]
elif a.mode in ('pilot', 'native'):
    commands += [('native', ['python3', '-B', str(tests/'validate_cls_epawf_m8.py'),
        '--build', str(build), '--output', str(out/'qualification'), '--mode', a.mode,
        '--gamma-scale', str(a.gamma_scale), '--variant', a.variant, '--beta', str(a.beta)], 21600)]
elif a.mode == 'transverse':
    commands += [('transverse', ['python3', '-B', str(tests/'validate_cls_epawf_m8_transverse.py'),
        '--build', str(build), '--output', str(out/'qualification'),
        '--gamma-scale', str(a.gamma_scale)], 7200)]
elif a.mode == 'regression':
    commands += [('regression', ['taskset', '-c', '0-23', 'ctest', '--test-dir', str(build),
        '--parallel', '12', '-E', '^RANSCLSEPAWFM8', '--output-on-failure',
        '--output-junit', str(out/'ctest.xml')], 21600)]
code = 0
for name, command, timeout in commands:
    started = time.monotonic()
    try:
        with (out/(name+'.log')).open('w') as log:
            proc = subprocess.run(command, cwd=out, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=timeout)
        code = proc.returncode
    except subprocess.TimeoutExpired:
        code = 124
    receipt['commands'].append(dict(name=name, command=command, exit_code=code,
        seconds=time.monotonic()-started, timeout_seconds=timeout,
        log_sha256=hashlib.sha256((out/(name+'.log')).read_bytes()).hexdigest()))
    save()
    print(name, code, flush=True)
    print((out/(name+'.log')).read_text()[-2500:], flush=True)
    if code:
        break
if (out/'runtime-before.json').is_file():
    after_command=['python3','-B',str(fingerprint),'--root',str(root),'--output',str(out/'runtime-after.json'),
        '--compare',str(out/'runtime-before.json')]
    with (out/'fingerprint-after.log').open('w') as log:
        after=subprocess.run(after_command,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=240)
    receipt['runtime_unchanged']=after.returncode==0
    receipt['runtime_after_command']=after_command
    if after.returncode and code==0:code=after.returncode
receipt.update(status='passed' if code == 0 else 'failed', exit_code=code,
    finished_utc=datetime.now(timezone.utc).isoformat(), source_end_commit=git(source, 'rev-parse', 'HEAD'),
    source_end_status=git(source, 'status', '--porcelain'), main_end_commit=git(main, 'rev-parse', 'HEAD'),
    main_end_status=git(main, 'status', '--porcelain'))
assert receipt['source_end_commit'] == revision and not receipt['source_end_status']
assert receipt['main_end_commit'] == main_revision and not receipt['main_end_status']
save()
raise SystemExit(code)
