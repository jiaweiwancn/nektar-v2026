#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Same guarded serial family, with separately measured monitor/launcher CPU."""
import argparse, os, resource, subprocess, sys, time
from pathlib import Path
from rans_qualification import THREADS, write_json
os.environ.update(THREADS)
from rans_worker_budget import spawn


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--control',type=Path,required=True)
    p.add_argument('--owned-child',action='store_true')
    p.add_argument('--queue-deadline-unix',type=float,default=1794805033+28560,
        help='Latest worker-lease admission, using the original CT15 wall allocation and closure reserve')
    p.add_argument('command',nargs=argparse.REMAINDER)
    a=p.parse_args();cmd=a.command[1:] if a.command and a.command[0]=='--' else a.command
    if a.owned_child:return subprocess.run(cmd,env=dict(os.environ)).returncode
    root=a.control.resolve();root.mkdir(parents=True,exist_ok=False)
    before=resource.getrusage(resource.RUSAGE_CHILDREN);start=time.monotonic();self_start=time.process_time();wait=0.
    guard=subprocess.Popen([sys.executable,'-B',str(Path(__file__).with_name('guard_rans_worker_budget.py')),'--root',str(root),'--project-limit','16'])
    try:
        child=None
        while child is None and time.time()<a.queue_deadline_unix:
            child=spawn(root,[sys.executable,'-B',str(Path(__file__).resolve()),'--owned-child','--control',str(root),'--',*cmd],env=dict(os.environ))
            if child is None:time.sleep(2);wait+=2
        code=child.wait() if child is not None else 124
    finally:
        before_guard=resource.getrusage(resource.RUSAGE_CHILDREN)
        (root/'worker-guard.stop').write_text('complete\n');guard.wait(timeout=10)
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    r=__import__('json').loads((root/'resource-final.json').read_text())
    launcher_CPU=time.process_time()-self_start
    children_CPU=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime
    guard_CPU=after.ru_utime+after.ru_stime-before_guard.ru_utime-before_guard.ru_stime
    write_json(root/'execution.json',dict(command=cmd,exit_code=code,threads=THREADS,control_parent_pid=os.getpid(),guard_pid=guard.pid,
        owned_family_pid=child.pid if child is not None else None,wall_seconds=time.monotonic()-start,allocation_wait_seconds=wait,
        status='completed' if child is not None else 'declined_before_spawn',queue_deadline_unix=a.queue_deadline_unix,
        process_cpu_seconds=children_CPU+launcher_CPU,launcher_self_CPU_seconds=launcher_CPU,worker_monitor_CPU_seconds=guard_CPU,owned_descendant_family_CPU_seconds=children_CPU-guard_CPU,
        resource=r,process_cpu_scope='owned descendants plus worker monitor plus launcher self; wall remains elapsed including pauses and allocation wait'))
    return code if code else (0 if r['aggregate_status']==r['project_status']=='pass' else 2)
if __name__=='__main__':raise SystemExit(main())
