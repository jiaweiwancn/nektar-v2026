#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Close an old CT15 launcher that is still waiting before its first spawn.

Read-only by default. This narrowly repairs accounting for a launcher started
before the queue deadline was added. Never signal a numerical child or another
campaign. Preserve the stopped shell parent until its child becomes a zombie,
so /proc supplies the launcher's final own plus reaped-child CPU counters.
"""
import argparse, hashlib, json, os, signal, time
from pathlib import Path


def read(path):return json.loads(path.read_text())
def write(path,value):
    with path.open('x') as f:json.dump(value,f,indent=2);f.write('\n')
def proc(pid):
    p=Path('/proc')/str(pid)
    fields=(p/'stat').read_text().rsplit(')',1)[1].split()
    return dict(pid=pid,parent=int(fields[1]),state=fields[0],start=int(fields[19]),
        ticks=[int(fields[i]) for i in (11,12,13,14)],
        command=(p/'cmdline').read_bytes().decode().split('\0')[:-1],uid=p.stat().st_uid)
def require(condition,message):
    if not condition:raise ValueError(message)
def stopped(pid,identity,timeout=5):
    until=time.monotonic()+timeout
    while time.monotonic()<until:
        p=proc(pid);require(p['start']==identity['start'],'PID identity changed')
        if p['state'] in ('T','t'):return p
        time.sleep(.02)
    raise ValueError('owned process did not stop')
def children(root):
    rows={}
    for p in Path('/proc').iterdir():
        if not p.name.isdigit():continue
        try:rows[int(p.name)]=proc(int(p.name))
        except (FileNotFoundError,ProcessLookupError,PermissionError):pass
    family={root}
    while True:
        more={pid for pid,p in rows.items() if p['parent'] in family}-family
        if not more:return [rows[p] for p in family-{root} if p in rows]
        family.update(more)


def main():
    own_CPU=time.process_time();started=time.monotonic()
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--control',type=Path,required=True)
    for key in ('launcher','parent','guard'):p.add_argument('--'+key,type=int,required=True)
    mode=p.add_mutually_exclusive_group()
    mode.add_argument('--cancel',action='store_true')
    mode.add_argument('--cancel-at-unix',type=float,help='Wait for admission; close only a still-unadmitted old launcher at this absolute deadline')
    a=p.parse_args();root=a.control.resolve()
    require(root.name.startswith('resources-') and root.parent.name=='ct-accepted-15-20261116','foreign CT15 control root')
    require(not (root/'execution.json').exists(),'family already closed')
    launcher,parent,guard=[proc(i) for i in (a.launcher,a.parent,a.guard)]
    require(all(q['uid']==os.getuid() for q in (launcher,parent,guard)),'foreign process owner')
    require(launcher['parent']==a.parent and guard['parent']==a.launcher,'foreign process ancestry')
    require(any(Path(s).name=='ct_accepted15_command.py' for s in launcher['command']) and str(root) in launcher['command'] and '--owned-child' not in launcher['command'],'foreign launcher')
    require(any(Path(s).name=='guard_rans_worker_budget.py' for s in guard['command']) and str(root) in guard['command'],'foreign worker guard')
    require(Path(parent['command'][0]).name in ('bash','sh','zsh'),'parent is not the owned command shell')
    require(not read(root/'worker-leases.json') and not (root/'live-budget.json').exists(),'numerical family already admitted')
    status=read(root/'worker-guard-status.json')
    witness=dict(schema='CT15PreSpawnQueueCancellationV1',control=str(root),launcher=launcher,parent=parent,guard=guard,worker_status=status,action='read_only_preflight')
    def watch_receipt(status):
        write(root/'queue-deadline-watch.json',dict(schema='CT15QueueDeadlineWatchV1',status=status,
            deadline_unix=a.cancel_at_unix,self_CPU=time.process_time()-own_CPU,wall=time.monotonic()-started,
            CPU_accounting='administrative helper charged against the fixed closure/transport/publication reserve, not added again to the launcher receipt'))
    if a.cancel_at_unix is not None:
        require(a.cancel_at_unix<=1794805033+28200,'deadline must leave the original 600-second publication reserve')
        require(not (root/'queue-deadline-watch.json').exists(),'deadline watcher already closed')
        while time.time()<a.cancel_at_unix:
            if (root/'execution.json').exists() or read(root/'worker-leases.json') or (root/'live-budget.json').exists():
                watch_receipt('family_admitted_or_closed');print('family admitted or closed; existing package guard remains authoritative');return 0
            time.sleep(min(10,max(0,a.cancel_at_unix-time.time())))
        a.cancel=True
    if not a.cancel:print(json.dumps(witness,indent=2));return 0
    # A stopped launcher cannot race the admission check. The guard remains
    # running until the launcher's ordinary finally block closes it.
    parent_stopped=launcher_stopped=False
    try:
        os.kill(a.parent,signal.SIGSTOP);parent_stopped=True;stopped(a.parent,parent)
        os.kill(a.launcher,signal.SIGSTOP);launcher_stopped=True;stopped(a.launcher,launcher)
        descendants=children(a.launcher)
        require(not read(root/'worker-leases.json') and not (root/'live-budget.json').exists(),'admission raced cancellation')
        require(all(q['pid']==a.guard or Path(q['command'][0]).name=='ps' for q in descendants if q['command']),'unexpected owned child; cancellation refused')
        ignored=int(next(s.split()[1] for s in (Path('/proc')/str(a.launcher)/'status').read_text().splitlines() if s.startswith('SigIgn:')),16)
        require(not ignored & (1<<(signal.SIGINT-1)),'launcher ignores SIGINT')
        os.kill(a.launcher,signal.SIGINT);os.kill(a.launcher,signal.SIGCONT);launcher_stopped=False
        until=time.monotonic()+15
        while True:
            final=proc(a.launcher);require(final['start']==launcher['start'],'launcher PID changed')
            if final['state']=='Z':break
            require(time.monotonic()<until,'launcher failed to close; no receipt manufactured')
            time.sleep(.05)
        require((root/'resource-final.json').is_file(),'worker audit did not close')
        require(not read(root/'worker-leases.json'),'child was spawned during cancellation')
        hz=os.sysconf('SC_CLK_TCK');boot=int(next(s.split()[1] for s in Path('/proc/stat').read_text().splitlines() if s.startswith('btime ')))
        elapsed=time.time()-(boot+launcher['start']/hz)+1.
        CPU=sum(final['ticks'])/hz;rounding=4/hz
        witness.update(action='cancelled_before_spawn',final_zombie_counters=final,final_CPU_ticks_per_second=hz,
            final_CPU_seconds=CPU,CPU_rounding_reserve=rounding,wall_upper_seconds=elapsed,
            scientific_execution='UNATTEMPTED',resource_final_sha256=hashlib.sha256((root/'resource-final.json').read_bytes()).hexdigest())
        write(root/'queue-cancellation.json',witness)
        command=launcher['command'];command=command[command.index('--')+1:]
        write(root/'execution.json',dict(command=command,exit_code=130,status='cancelled_before_spawn',
            control_parent_pid=a.launcher,guard_pid=a.guard,owned_family_pid=None,
            wall_seconds=elapsed,allocation_wait_seconds=elapsed,process_cpu_seconds=CPU+rounding,
            resource=read(root/'resource-final.json'),queue_cancellation_sha256=hashlib.sha256((root/'queue-cancellation.json').read_bytes()).hexdigest(),
            process_cpu_scope='final /proc launcher own plus reaped-child ticks, including its guard and all lease probes; four-tick rounding reserve; helper/publication separately reserved',
            wall_scope='conservative process-birth-to-receipt upper bound; no scientific child was admitted'))
        print(json.dumps(dict(status='cancelled_before_spawn',CPU_upper=CPU+rounding,wall_upper=elapsed)))
    finally:
        if launcher_stopped:os.kill(a.launcher,signal.SIGCONT)
        if parent_stopped:os.kill(a.parent,signal.SIGCONT)
    if a.cancel_at_unix is not None:watch_receipt('cancelled_before_spawn')
    return 0


if __name__=='__main__':raise SystemExit(main())
