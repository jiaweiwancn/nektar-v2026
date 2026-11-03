"""Build an additional immutable R2 epoch without replacing earlier binaries."""
import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--bundle', type=Path, required=True)
p.add_argument('--revision', required=True)
p.add_argument('--epoch', required=True)
a = p.parse_args()
base = Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised')
root = base/'m8-r2'
assert a.epoch.startswith('epoch-') and '/' not in a.epoch
epoch = root/a.epoch
epoch.mkdir(exist_ok=False)
source, build = epoch/'source', epoch/'build'
(epoch/'tmp').mkdir()
lock = (root/'validation.lock').open('a')
fcntl.flock(lock, fcntl.LOCK_EX)
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
git = lambda path, *argv: subprocess.check_output(['git', '-C', str(path), *argv], text=True).strip()
preserved_sources = {str(path): dict(commit=git(path, 'rev-parse', 'HEAD'),
    status=git(path, 'status', '--porcelain')) for path in (
        Path('/media/jay/新加卷/nektar-v2026'), base/'m8/source', base/'m8-r1/source',
        base/'m8-r1/controlled/source', base/'m8-r1/finalization-001/source', root/'source')}
assert all(not item['status'] for item in preserved_sources.values())
preserved_images = {str(path): sha(path) for path in (
    base/'m8/build/IncNavierStokesSolver', base/'m8-r1/build/IncNavierStokesSolver',
    base/'m8-r1/controlled/build/IncNavierStokesSolver', root/'build/IncNavierStokesSolver',
    root/'build/RANSCurvedWallProbe', root/'build/RANSPressureActionProbe')}
subprocess.run(['git', '-C', str(root/'source'), 'fetch', '--no-tags', str(a.bundle), 'HEAD'], check=True)
assert git(root/'source', 'rev-parse', 'FETCH_HEAD') == a.revision
subprocess.run(['git', '-C', str(root/'source'), 'worktree', 'add', '--detach', str(source), a.revision], check=True)
keys = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS',
        'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')
env = dict(os.environ, **{key: '1' for key in keys})
env.update(PYTHONDONTWRITEBYTECODE='1', TMPDIR=str(epoch/'tmp'), LD_LIBRARY_PATH=str(build/'geometry-lib'))
receipt = dict(schema=1, status='running', source_commit=a.revision,
    source_clean_start=True, runner_sha256=sha(__file__), bundle_sha256=sha(a.bundle),
    preserved_sources=preserved_sources, preserved_images=preserved_images,
    compiler_workers=16, maximum_concurrent_analysis_workers=3, allocation_ceiling=19,
    server_worker_limit=72, numerical_threads={key: env[key] for key in keys},
    started_utc=datetime.now(timezone.utc).isoformat(), commands=[])
def save():
    (epoch/'build-receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
save()
commands = [('configure', ['cmake', '-S', str(source/'solvers/IncNavierStokesSolver/Tests/RANS'),
    '-B', str(build), '-DRANS_BUILD_SOLVER=ON', '-DCMAKE_BUILD_TYPE=Release',
    '-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake']),
    ('build', ['taskset', '-c', '0-31', 'cmake', '--build', str(build), '--parallel', '16']),
    ('fingerprint', ['python3', '-B', str(source/'docs/plans/cls-epawf-revised-evidence/m7/fingerprint-runtime.py'),
        '--root', str(epoch), '--output', str(epoch/'runtime-fingerprint.json')])]
code = 0
for name, argv in commands:
    started = time.monotonic()
    log = epoch/(name+'.log')
    try:
        with log.open('w') as stream:
            code = subprocess.run(argv, cwd=epoch, env=env, stdout=stream,
                stderr=subprocess.STDOUT, timeout=1800 if name=='build' else 240).returncode
    except subprocess.TimeoutExpired:
        code = 124
    receipt['commands'].append(dict(name=name, argv=argv, cwd=str(epoch), exit_code=code,
        seconds=time.monotonic()-started, log_sha256=sha(log)))
    save()
    print(name, code, log.read_text()[-2500:], flush=True)
    if code:
        break
receipt.update(status='passed' if code==0 else 'failed',
    source_clean_end=not bool(git(source, 'status', '--porcelain')),
    source_commit_end=git(source, 'rev-parse', 'HEAD'),
    sources_preserved=all(git(Path(path), 'rev-parse', 'HEAD')==item['commit'] and
        git(Path(path), 'status', '--porcelain')==item['status'] for path, item in preserved_sources.items()),
    images_preserved=all(sha(path)==digest for path, digest in preserved_images.items()),
    finished_utc=datetime.now(timezone.utc).isoformat())
save()
assert receipt['sources_preserved'] and receipt['images_preserved']
raise SystemExit(code)
