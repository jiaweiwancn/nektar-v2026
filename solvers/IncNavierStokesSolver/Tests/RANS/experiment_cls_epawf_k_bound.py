#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent fixed-horizon post-k experiments, with bounded concurrency.

Numerical rejection is retained as a result. A failed case never shortens
another case's horizon. Compilation must finish before this driver starts.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import platform
import subprocess

from run_cls_epawf_picard import load, require, write_json
from run_cls_epawf_pimple import ROOT, identity


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--solver', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--checkpoint', type=Path)
    p.add_argument('--jobs', type=int, default=3)
    a = p.parse_args()
    require(platform.system() == 'Linux' and 1 <= a.jobs <= 3, 'Linux; at most three one-thread runs')
    a.output = a.output.resolve(); a.output.mkdir(parents=True, exist_ok=False)
    a.solver = a.solver.resolve(strict=True)
    # Include helper threads in the conservative worker audit. Do not
    # interrupt any existing solver, compiler or unrelated experiment.
    processes = subprocess.check_output(['ps', '-eo', 'comm=,nlwp='], text=True)
    workers = sum(int(line.split()[-1]) for line in processes.splitlines()
                  if line.split()[0] in ('IncNavierStokes', 'IncNavierStokesSolver', 'cc1plus', 'RANSWallResolut'))
    require(workers + 3*a.jobs <= 72, 'insufficient worker allocation')
    runner = Path(__file__).with_name('run_cls_epawf_pimple.py')
    cases = []
    for i, dt in enumerate((.000625, .0003125, .00015625)):
        cases.append((f'bound-{i}', dt, .005, 'MassProjection', ['--diagnostic'] +
                      (['--budget-time', '.0025'] if i == 0 else [])))
    for i, dt in enumerate((.000625, .0003125, .00015625)):
        cases.append((f'off-{i}', dt, .0025, 'Off', ['--diagnostic'] +
                      (['--budget-time', '.0025'] if i == 0 else [])))
    for name, floor in (('floor-low', '1e-10'), ('floor-high', '1e-6')):
        cases.append((name, .000625, .005, 'MassProjection', ['--diagnostic', '--floor-relative', floor]))
    cases += [('strict', .000625, .0025, 'MassProjection', []),
              ('hard-reject', .000625, .0025, 'MassProjection', ['--diagnostic', '--hard-reject-step', '4']),
              ('inactive-on', .00015625, .000625, 'MassProjection', ['--diagnostic']),
              ('inactive-off', .00015625, .000625, 'Off', ['--diagnostic'])]
    if a.checkpoint:
        cases.append(('restart', .000625, .005, 'MassProjection',
                      ['--diagnostic', '--restart', str(a.checkpoint.resolve())]))
    summary = dict(schema=1, source=subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip(),
                   executable=identity(a.solver), workers_existing=workers, concurrent_jobs=a.jobs,
                   ranks_per_job=1, threads_per_job=1, cases={}, status='running')
    write_json(a.output/'summary.json', summary)

    def execute(case):
        name, dt, final, mode, extra = case
        command = ['python3', '-B', str(runner), '--solver', str(a.solver), '--output', str(a.output/name),
                   '--dt', str(dt), '--final-time', str(final), '--k-bound', mode, '--snapshots', *extra]
        with (a.output/(name+'.driver.log')).open('w') as log:
            code = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT).returncode
        return name, command, code

    with ThreadPoolExecutor(max_workers=a.jobs) as pool:
        for future in as_completed([pool.submit(execute, case) for case in cases]):
            name, command, code = future.result()
            row = dict(command=command, driver_exit_code=code)
            if (a.output/name/'manifest.json').exists():
                m = load(a.output/name/'manifest.json')
                row.update(status=m['status'], accepted_steps=m['accepted_steps'], final_time=m['final_time'],
                           seconds=m['execution']['seconds'], failure=m['failure'])
            summary['cases'][name] = row
            write_json(a.output/'summary.json', summary)
            print(json.dumps(dict(name=name, **row)), flush=True)
    summary['status'] = 'recorded' if all(x['driver_exit_code'] == 0 for x in summary['cases'].values()) else 'driver_failure'
    write_json(a.output/'summary.json', summary)
    require(summary['status'] == 'recorded', 'inspect driver failures; numerical rejections remain valid results')


if __name__ == '__main__':
    main()
