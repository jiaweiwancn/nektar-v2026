#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Stop the owned CT14 retry before its total package allocation is exhausted.

Run inside ct_accepted3_command.py. Include both outer launcher and worker
monitor in live CPU; leave them alive to close their receipt after a stop.
"""
import argparse,json,os,resource,signal,subprocess,time
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
from ct_accepted_open import require
from ct_accepted_feasibility import write_json

def processes():
    out={}
    for path in Path('/proc').iterdir():
        if not path.name.isdigit():continue
        try:
            # comm may contain spaces or parentheses. Fields after its last
            # closing parenthesis begin with state,ppid,... at stat field3.
            fields=(path/'stat').read_text().rsplit(')',1)[1].split()
            out[int(path.name)]=dict(parent=int(fields[1]),start=int(fields[19]),ticks=sum(int(fields[k]) for k in (11,12,13,14)))
        except (FileNotFoundError,ProcessLookupError,PermissionError):pass
    return out

def family(rows,root):
    owned={root}
    while True:
        more={pid for pid,r in rows.items() if r['parent'] in owned}-owned
        if not more:return owned
        owned.update(more)

def main():
    resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3))
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--control',type=Path,required=True);p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args();cmd=a.command[1:] if a.command[:1]==['--'] else a.command
    rows=processes();parent=rows[os.getpid()]['parent'];outer=rows[parent]['parent'];require(outer in rows,'missing owned outer launcher')
    closed=0.;receipts={}
    for folder in sorted(a.root.glob('resources-*')):
        if folder==a.control:continue
        f=folder/'execution.json';require(f.is_file(),'unclosed previous family:'+str(folder));r=json.loads(f.read_text());closed+=r['process_cpu_seconds'];receipts[folder.name]=r['process_cpu_seconds']
    hz=os.sysconf('SC_CLK_TCK');child=subprocess.Popen(cmd,env=dict(os.environ));maximum=closed;reason=None;last_scan=0.;retained=0;witness={}
    while child.poll() is None:
        rows=processes();owned=family(rows,outer);live=sum(rows[pid]['ticks'] for pid in owned if pid in rows)/hz;cpu=closed+live;maximum=max(maximum,cpu);whole=time.time()-1794751358
        if time.monotonic()-last_scan>30:
            retained=sum(f.lstat().st_size for f in a.root.rglob('*') if f.is_file() or f.is_symlink());last_scan=time.monotonic()
        witness=dict(closed_CPU=closed,live_owned_family_CPU=live,total_CPU_sample=cpu,whole_elapsed_wall=whole,retained_bytes=retained)
        if cpu>=23760:reason='CPU allocation reserve'
        elif whole>=28560:reason='wall allocation reserve'
        elif retained>=16*1024**3-256*1024**2:reason='output allocation reserve'
        if reason:
            victims=family(rows,child.pid)
            for pid in sorted(victims,reverse=True):
                try:os.kill(pid,signal.SIGTERM)
                except ProcessLookupError:pass
            time.sleep(2)
            for pid in victims:
                try:os.kill(pid,signal.SIGKILL)
                except ProcessLookupError:pass
            break
        time.sleep(.5)
    code=child.wait();write_json(a.control/'package-budget-watch.json',dict(schema='CT14PackageBudgetWatchV1',status='stopped_before_allocation' if reason else 'completed',reason=reason,command=cmd,exit_code=code,closed_families_CPU=receipts,maximum_sampled_total_CPU=maximum,last_sample=witness,reserves=dict(CPU=240,wall=240,output_bytes=256*1024**2),scope='live /proc own plus waited-child CPU over outer launcher/guard/owned family; original closed launcher receipts are the final accounting authority; no scientific acceptance supplied'))
    return 124 if reason else code

if __name__=='__main__':raise SystemExit(main())
