#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""One owned serial command family, explicit thread limits and peak evidence.

Use a NEW control directory outside source. The numerical command may have its
own immutable output directory. Its exit status is returned unchanged.
"""
import argparse,os,resource,subprocess,sys,time
from pathlib import Path
from rans_qualification import THREADS,write_json
os.environ.update(THREADS)
from rans_worker_budget import spawn


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--control',type=Path,required=True);p.add_argument('--owned-child',action='store_true');p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args()
    command=a.command[1:] if a.command and a.command[0]=='--' else a.command
    if not command:raise ValueError('a command after -- is required')
    root=a.control.resolve()
    if a.owned_child:
        # Keep this waiting parent visible in the lease: its root argument
        # identifies all native/MPI descendants even if their argv differs.
        return subprocess.run(command,env=dict(os.environ)).returncode
    root.mkdir(parents=True,exist_ok=False)
    start=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN);wait=0.
    guard=subprocess.Popen([sys.executable,'-B',str(Path(__file__).with_name('guard_rans_worker_budget.py')),'--root',str(root)],env=dict(os.environ))
    try:
        child=None
        while child is None:
            child=spawn(root,[sys.executable,'-B',str(Path(__file__).resolve()),'--owned-child','--control',str(root),'--',*command],env=dict(os.environ))
            if child is None:time.sleep(2);wait+=2
        code=child.wait()
    finally:
        (root/'worker-guard.stop').write_text('owned command complete\n');guard.wait(timeout=5)
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    write_json(root/'execution.json',dict(command=command,exit_code=code,threads=THREADS,
        wall_seconds=time.monotonic()-start,allocation_wait_seconds=wait,
        process_cpu_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
        process_cpu_scope='owned command family plus its short-lived resource monitor',
        family_reservation=8,project_limit=16,aggregate_limit=72,startup_reserve=16))
    return code
if __name__=='__main__':raise SystemExit(main())
