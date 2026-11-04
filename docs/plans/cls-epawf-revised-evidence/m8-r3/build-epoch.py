#!/usr/bin/env python3
"""Build a fresh isolated R3 epoch on allocated Linux workers only."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import time


def pin(path):
    path = Path(path).resolve(); digest = hashlib.sha256()
    with path.open('rb') as stream:
        for part in iter(lambda: stream.read(1024*1024), b''): digest.update(part)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=digest.hexdigest())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('source', 'build', 'prefix', 'output'): p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--expected-source', required=True); p.add_argument('--jobs', type=int, required=True)
    p.add_argument('--target', action='append', default=[])
    a = p.parse_args()
    if platform.system() != 'Linux' or not 1 <= a.jobs <= 16:
        raise ValueError('Authoritative Linux and an explicit allocation of 1..16 compiler workers required')
    source, build, output = a.source.resolve(), a.build.resolve(), a.output.resolve()
    if build.exists() or output.exists() or build.is_relative_to(source) or output.is_relative_to(source):
        raise ValueError('Fresh external build and receipt paths required')
    git = lambda *args: subprocess.check_output(['git', '-C', str(source), *args], text=True).strip()
    if git('rev-parse', 'HEAD') != a.expected_source or git('status', '--porcelain'):
        raise ValueError('Wrong or unclean producing source')
    output.mkdir(parents=True)
    environment = dict(os.environ, **{key:'1' for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
        'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')})
    receipt = dict(schema=1, status='running', source_commit=a.expected_source, source_clean_start=True,
        builder=pin(__file__), started_utc=datetime.now(timezone.utc).isoformat(), compiler_workers=a.jobs,
        source=str(source), build=str(build), prefix=str(a.prefix.resolve()), commands=[])
    def save():
        (output/'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    commands = [('configure', ['cmake', '-S', str(source/'solvers/IncNavierStokesSolver/Tests/RANS'),
        '-B', str(build), '-DRANS_BUILD_SOLVER=ON', '-DCMAKE_BUILD_TYPE=Release',
        '-DCMAKE_PREFIX_PATH='+str(a.prefix.resolve()), '-DCMAKE_EXPORT_COMPILE_COMMANDS=ON',
        '-DCMAKE_CXX_FLAGS=-fno-fast-math -ffp-contract=off']),
        ('build', ['cmake', '--build', str(build), '--parallel', str(a.jobs), '--target',
            *(a.target or ['IncNavierStokesSolver', 'RANSPressureActionProbe'])])]
    save()
    try:
        for name, argv in commands:
            row = dict(name=name, argv=argv); receipt['commands'].append(row); save()
            start = time.monotonic(); log = output/(name+'.log')
            with log.open('w') as stream:
                completed = subprocess.run(argv, cwd=output, env=environment, stdout=stream, stderr=subprocess.STDOUT)
            row.update(exit_code=completed.returncode, seconds=time.monotonic()-start, log=pin(log)); save()
            if completed.returncode: raise ValueError('Build command failed: '+name)
        receipt['status'] = 'compiled'
        if (build/'IncNavierStokesSolver').exists(): receipt['executable'] = pin(build/'IncNavierStokesSolver')
        receipt['compile_commands'] = pin(build/'compile_commands.json')
    except Exception as error:
        receipt.update(status='failed', failure=dict(type=type(error).__name__, message=str(error)))
    finally:
        receipt.update(source_commit_end=git('rev-parse', 'HEAD'), source_clean_end=not bool(git('status', '--porcelain')),
            finished_utc=datetime.now(timezone.utc).isoformat()); save()
    return 0 if receipt['status'] == 'compiled' and receipt['source_clean_end'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
