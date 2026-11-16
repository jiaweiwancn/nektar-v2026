#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""One CT16 family, with the finite package watch outside worker suspension."""
import argparse,json,os,resource,signal,subprocess,sys,time
from pathlib import Path
from rans_qualification import THREADS,write_json
os.environ.update(THREADS)
from rans_worker_budget import spawn


def read(path):return json.loads(path.read_text())
def stat(pid):
    p=Path('/proc')/str(pid);v=(p/'stat').read_text().rsplit(')',1)[1].split()
    return dict(pid=pid,parent=int(v[1]),state=v[0],start=int(v[19]),ticks=sum(int(v[k]) for k in (11,12,13,14)))
def tree(root):
    rows={};todo=[root]
    while todo:
        pid=todo.pop()
        if pid in rows:continue
        try:
            rows[pid]=stat(pid)
            for thread in (Path('/proc')/str(pid)/'task').iterdir():
                try:todo.extend(map(int,(thread/'children').read_text().split()))
                except (FileNotFoundError,ProcessLookupError,PermissionError):pass
        except (FileNotFoundError,ProcessLookupError,PermissionError):pass
    return rows
def logical_bytes(root):
    n=0
    for folder,dirs,files in os.walk(root):
        for name in files:
            try:n+=(Path(folder)/name).lstat().st_size
            except FileNotFoundError:pass
    return n
def own_and_reaped_CPU():
    a=resource.getrusage(resource.RUSAGE_SELF);b=resource.getrusage(resource.RUSAGE_CHILDREN)
    return a.ru_utime+a.ru_stime+b.ru_utime+b.ru_stime
def stop_owned(root,guard,control):
    """Stable handles, frozen CPU witness, then kill only this owned subtree.

    No numerical acceptance is supplied by an interrupted family. Its receipt
    charges a conservative live/reaped CPU upper bound even if an interrupted
    intermediate controller cannot reap its children.
    """
    handles={};before=own_and_reaped_CPU();own=time.process_time();hz=os.sysconf('SC_CLK_TCK')
    def freeze(pid):
        if pid in handles:return
        try:
            first=stat(pid);fd=os.pidfd_open(pid);second=stat(pid)
            if first['start']!=second['start']:os.close(fd);raise RuntimeError('owned PID identity changed')
            signal.pidfd_send_signal(fd,signal.SIGSTOP);handles[pid]=fd
        except ProcessLookupError:pass
        except FileNotFoundError:pass
    # Stop the worker supervisor first so it cannot resume the family while
    # this separate package supervisor freezes its live accounting witness.
    freeze(guard.pid)
    for _ in range(4):
        rows=tree(os.getpid())
        for pid in rows:
            if pid!=os.getpid():freeze(pid)
        time.sleep(.05)
    rows=tree(os.getpid());witness=sum(r['ticks'] for r in rows.values())/hz
    for pid,fd in sorted(handles.items(),reverse=True):
        try:signal.pidfd_send_signal(fd,signal.SIGKILL)
        except ProcessLookupError:pass
        finally:os.close(fd)
    result=dict(schema='CT16OwnedBudgetStopV1',processes=rows,CPU_before_stop=before,
        frozen_tree_CPU=witness,rounding_and_cleanup_CPU_reserve=2.,
        CPU_upper=max(before,witness)+max(0.,time.process_time()-own)+2.,
        scope='owned descendants/worker guard only; all handles opened after ancestry discovery and start-identity checks; no external process signalled')
    write_json(control/'budget-stop-witness.json',result);return result


def main():
    resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3))
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--control',type=Path,required=True)
    p.add_argument('--owned-child',action='store_true')
    p.add_argument('--expected-CPU',type=float,default=0);p.add_argument('--expected-wall',type=float,default=0)
    p.add_argument('--expected-new-bytes',type=int,default=0)
    p.add_argument('--family-wall-limit',type=float,default=0,help='Optional smaller per-family control; does not enlarge package limits')
    p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args();cmd=a.command[1:] if a.command[:1]==['--'] else a.command
    if a.owned_child:return subprocess.run(cmd,env=dict(os.environ)).returncode
    root=a.root.resolve();control=a.control.resolve();allocation=read(root/'allocation.json')
    if control.parent!=root or not control.name.startswith('resources-'):raise ValueError('foreign control directory')
    start_unix=allocation['start_unix'];limits=allocation['limits'];closed=0.;prior={}
    if start_unix!=1794835825 or limits!=dict(soft_CPU=23760,soft_wall=28560,soft_live_bytes=15*1024**3):raise ValueError('CT16 allocation identity/limits changed')
    if not cmd or a.expected_CPU<=0 or a.expected_wall<=0 or a.expected_new_bytes<0:raise ValueError('missing whole-scope prelaunch estimate')
    for folder in sorted(root.glob('resources-*')):
        if not (folder/'execution.json').is_file():raise ValueError('unclosed prior family:'+str(folder))
        r=read(folder/'execution.json');closed+=r['process_cpu_seconds'];prior[folder.name]=r['process_cpu_seconds']
    control.mkdir(exist_ok=False);started=time.monotonic();wait=0.;maximum=closed;peak_disk=0;last_disk=-1e9;last_write=-1e9;written={};child=None;stop=None;code=124;stop_witness=None;sample={}
    guard=subprocess.Popen([sys.executable,'-B',str(Path(__file__).with_name('guard_rans_worker_budget.py')),'--root',str(control),'--project-limit','16'])
    try:
        while True:
            now=time.monotonic();rows=tree(os.getpid());live=sum(r['ticks'] for r in rows.values())/os.sysconf('SC_CLK_TCK');CPU=closed+live;maximum=max(maximum,CPU);wall=time.time()-start_unix
            for pid,r in rows.items():
                try:
                    io=dict(s.split(':',1) for s in (Path('/proc')/str(pid)/'io').read_text().splitlines());key=(pid,r['start']);written[key]=max(written.get(key,0),int(io['wchar']))
                except (FileNotFoundError,ProcessLookupError,PermissionError,KeyError):pass
            if now-last_disk>=30:
                retained=logical_bytes(root);peak_disk=max(peak_disk,retained);last_disk=now
            sample=dict(closed_CPU=closed,live_owned_family_CPU=live,total_CPU_sample=CPU,whole_elapsed_wall=wall,retained_bytes=retained,peak_sampled_live_bytes=peak_disk,sampled_cumulative_process_wchar=sum(written.values()))
            if child is not None and child.poll() is not None:code=child.returncode;break
            if CPU>=limits['soft_CPU']:stop='package CPU reserve'
            elif wall>=limits['soft_wall']:stop='original package wall reserve'
            elif retained>=limits['soft_live_bytes']:stop='package live-byte reserve'
            elif a.family_wall_limit and now-started>=a.family_wall_limit:stop='declared family wall control'
            elif child is None and (CPU+a.expected_CPU>limits['soft_CPU'] or wall+a.expected_wall>limits['soft_wall'] or retained+a.expected_new_bytes>limits['soft_live_bytes']):stop='whole-scope prelaunch admission declined'
            if stop:
                stop_witness=stop_owned(root,guard,control);break
            if child is None:
                child=spawn(control,[sys.executable,'-B',str(Path(__file__).resolve()),'--owned-child','--root',str(root),'--control',str(control),'--',*cmd],env=dict(os.environ))
                if child is None:wait+=.5
            if now-last_write>=2:
                write_json(control/'live-budget.json',dict(maximum_sampled_total_CPU=maximum,last_sample=sample));last_write=now
            time.sleep(.5)
    except BaseException as error:
        stop='launcher exception:'+type(error).__name__+':'+str(error)
        stop_witness=stop_owned(root,guard,control)
    finally:
        if stop_witness is None:
            (control/'worker-guard.stop').write_text('complete\n')
            try:guard.wait(timeout=10)
            except subprocess.TimeoutExpired:
                stop=stop or 'worker guard failed to close';stop_witness=stop_owned(root,guard,control)
        if child is not None:
            try:child.wait(timeout=10)
            except subprocess.TimeoutExpired:
                stop=stop or 'owned family failed to reap';stop_witness=stop_owned(root,guard,control);child.wait(timeout=10)
        guard.wait(timeout=10)
    r=read(control/'resource-final.json') if (control/'resource-final.json').is_file() else dict(aggregate_status='incomplete',project_status='incomplete',aggregate_limit=72,project_limit=16)
    if not (control/'resource-final.json').exists():write_json(control/'resource-final.json',r)
    CPU=max(own_and_reaped_CPU(),stop_witness['CPU_upper'] if stop_witness else 0.)
    peak_disk=max(peak_disk,logical_bytes(root));sample['peak_sampled_live_bytes']=peak_disk
    write_json(control/'package-budget-watch.json',dict(schema='CT16ExternalPackageTreeV1',status='stopped_before_allocation' if stop else 'completed',reason=stop,closed_families_CPU=prior,maximum_sampled_total_CPU=maximum,last_sample=sample,scope='outside worker-paused family; /proc own+waited-child CPU across owned process tree; final reaped receipt authoritative unless explicit conservative stop witness'))
    write_json(control/'execution.json',dict(command=cmd,exit_code=124 if stop else code,status='stopped_by_package_guard' if stop else 'completed',threads=THREADS,wall_seconds=time.monotonic()-started,allocation_wait_seconds=wait,process_cpu_seconds=CPU,resource=r,stop_reason=stop,allocation_start_unix=start_unix,expected=dict(CPU=a.expected_CPU,wall=a.expected_wall,new_bytes=a.expected_new_bytes),process_cpu_scope='launcher self plus reaped descendants including worker monitor; conservative frozen-tree upper bound on forced-stop exception; no historical CT15 charge'))
    return 124 if stop else code if code else 0 if r['aggregate_status']==r['project_status']=='pass' else 2


if __name__=='__main__':raise SystemExit(main())
