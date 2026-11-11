# SPDX-License-Identifier: MIT
"""Finite process-tree monitoring for the bounded serial WJ-6 package."""
import math
import os
from pathlib import Path
import signal
import subprocess
import time
import re

from run_cls_epawf_picard import THREADS

GIB = 1024**3
DEFAULT = dict(native_seconds=600., observation_seconds=480., arithmetic_seconds=600.,
               rss_bytes=8*GIB, invocation_output_bytes=8*GIB,
               policy_output_bytes=32*GIB, server_workers=72)


def validate(limits):
    if set(limits) != set(DEFAULT):
        raise ValueError('incomplete finite WJ6 resource design')
    for key, value in limits.items():
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
            raise ValueError('invalid finite WJ6 resource limit: '+key)
        if key.endswith('_bytes') or key == 'server_workers':
            if int(value) != value:
                raise ValueError('resource count must be an exact integer: '+key)
    if limits['server_workers'] != 72:
        raise ValueError('server worker ceiling must remain72')


def scientific_workers():
    rows = []
    text = subprocess.check_output(['ps', '-eo', 'pid=,ppid=,nlwp=,args='], text=True)
    for line in text.splitlines():
        pid, parent, threads, args = line.strip().split(None, 3)
        command = args.split()[0].split('/')[-1]
        if any(command.startswith(s) for s in ('IncNavierStokes', 'RANS', 'python', 'cc1',
                                               'orted', 'mpirun', 'mpiexec', 'lto1', 'f951')):
            rows.append(dict(pid=int(pid), parent=int(parent), threads=int(threads), args=args))
    return dict(workers=sum(r['threads'] for r in rows), processes=rows)


def tree_usage(pid):
    queue, seen, rss, threads = [pid], set(), 0, 0
    while queue:
        current = queue.pop()
        if current in seen:
            continue
        seen.add(current)
        try:
            status = (Path('/proc')/str(current)/'status').read_text()
            values = dict(line.split(':', 1) for line in status.splitlines() if ':' in line)
            rss += int(values.get('VmRSS', '0 kB').split()[0])*1024
            threads += int(values.get('Threads', '0'))
            children = Path('/proc')/str(current)/'task'/str(current)/'children'
            queue.extend(int(p) for p in children.read_text().split())
        except (FileNotFoundError, ProcessLookupError):
            pass
    return rss, threads


def volume(folder):
    total=0
    for p in Path(folder).rglob('*'):
        try:
            if p.is_file():total+=p.stat().st_size
        except FileNotFoundError:pass # Concurrent atomic rename, not lost data.
    return total


def run_limited(argv, folder, log, limits, kind='native', seconds=None):
    """Only the launched process group may be stopped by a ceiling."""
    validate(limits)
    timeout = limits[kind+'_seconds'] if seconds is None else seconds
    if not math.isfinite(timeout) or timeout <= 0:
        raise ValueError('invalid finite process timeout')
    allocation = scientific_workers()
    if allocation['workers']+4 > 72:
        raise RuntimeError('RESOURCE_LIMIT: current scientific workers plus4 exceed72')
    start = time.monotonic()
    peak_rss = peak_threads = 0
    reason = None
    last_volume = last_workers = start
    output = volume(folder)
    time_file=Path(folder)/'resource-time.txt'
    if time_file.exists():raise ValueError('unreconciled resource-time output exists')
    launcher=['/usr/bin/time','-v','-o',str(time_file.resolve()),'--',*list(map(str,argv))]
    with Path(log).open('x') as stream:
        proc = subprocess.Popen(launcher, cwd=folder, stdout=stream,
                                stderr=subprocess.STDOUT, start_new_session=True,
                                env=dict(os.environ, **THREADS))
        while proc.poll() is None:
            rss, threads = tree_usage(proc.pid)
            peak_rss, peak_threads = max(peak_rss, rss), max(peak_threads, threads)
            now = time.monotonic()
            if now-last_volume >= 2:
                output, last_volume = volume(folder), now
            if now-start > timeout:
                reason = 'process_wall_seconds'
            elif peak_rss > limits['rss_bytes']:
                reason = 'process_tree_rss_bytes'
            elif output > limits['invocation_output_bytes']:
                reason = 'invocation_output_bytes'
            if now-last_workers >= 3:
                if scientific_workers()['workers'] > 72:
                    reason = 'server_worker_ceiling'
                last_workers = now
            if reason:
                os.killpg(proc.pid, signal.SIGTERM)
                try:
                    proc.wait(timeout=3)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                break
            time.sleep(.1)
        exit_code = proc.wait()
    match=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',time_file.read_text()) if time_file.exists() else None
    measured_peak=int(match.group(1))*1024 if match else None
    return dict(exit_code=exit_code if not reason else 124, limit_reason=reason,
        seconds=time.monotonic()-start, peak_rss_bytes=max(peak_rss,measured_peak or 0),
        sampled_process_tree_peak_rss_bytes=peak_rss,wait4_maximum_rss_bytes=measured_peak,
        peak_process_threads=peak_threads, output_bytes=volume(folder),
        timeout_seconds=timeout, limits=limits, allocation=allocation, threads=THREADS,
        launcher=launcher,resource_time_file=str(time_file.resolve()))
