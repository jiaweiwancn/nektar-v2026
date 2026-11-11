#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Whole-process WJ-7 receipts/RSS ceiling, including parent algebra and I/O."""
import argparse
import json
import os
from pathlib import Path
import signal
import subprocess
import time
from wj_immutable import Ledger,read,write_once,atomic,file_identity,digest
from wj6_resources import tree_usage,scientific_workers,DEFAULT
from run_cls_epawf_picard import THREADS,require


def owned_groups(pid):
    todo=[pid];seen=set();groups=set()
    while todo:
        q=todo.pop()
        if q in seen:continue
        seen.add(q)
        try:
            groups.add(os.getpgid(q))
            todo.extend(map(int,(Path('/proc')/str(q)/'task'/str(q)/'children').read_text().split()))
        except (OSError,ProcessLookupError):pass
    return groups


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--label',required=True)
    p.add_argument('--seconds',type=float,required=True);p.add_argument('--control-units',type=int,default=0)
    p.add_argument('argv',nargs=argparse.REMAINDER);a=p.parse_args();argv=a.argv[1:] if a.argv[:1]==['--'] else a.argv
    require(argv and a.seconds>0,'finite command required');os.environ.update(THREADS)
    root=a.root.resolve();out=root/'operations'/a.label;out.mkdir(parents=True,exist_ok=False)
    ledgers={q:Ledger(root,q,cap=800) for q in ('EP','Budget')}
    def session_seconds(l):
        return sum(s.get('seconds',max(0,time.time()-s['started_epoch'])) for s in read(l.file).get('sessions',{}).values())
    before={q:session_seconds(l) for q,l in ledgers.items()}
    allocation=scientific_workers();require(allocation['workers']+8<=72,'RESOURCE_LIMIT: outer controller allocation')
    intent=dict(argv=argv,limits=DEFAULT,timeout=a.seconds,control_units=a.control_units,allocation=allocation,
        source=file_identity(__file__),threads=THREADS)
    write_once(out/'intent.json',intent)
    # The child accounts its own sessions/attempts; this reservation covers
    # only unaccounted wrapper/launch/I/O time. Unknown work remains reserved.
    label='outer/'+a.label
    for l in ledgers.values():l.reserve(label,intent,a.control_units,60)
    start=time.monotonic();peak=threads=0;failure=None
    with (out/'command.log').open('x') as log:
        child=subprocess.Popen(['/usr/bin/time','-v',*argv],stdout=log,stderr=subprocess.STDOUT,start_new_session=True,env=dict(os.environ,**THREADS))
        while child.poll() is None:
            rss,n=tree_usage(child.pid);peak=max(peak,rss);threads=max(threads,n)
            if rss>DEFAULT['rss_bytes']:failure='RESOURCE_LIMIT: whole-process RSS'
            if time.monotonic()-start>a.seconds:failure='RESOURCE_LIMIT: whole-process timeout'
            if scientific_workers()['workers']>72:failure='RESOURCE_LIMIT: current total workers'
            if failure:
                for pg in owned_groups(child.pid):
                    try:os.killpg(pg,signal.SIGTERM)
                    except ProcessLookupError:pass
                break
            time.sleep(.5)
        try:code=child.wait(timeout=10)
        except subprocess.TimeoutExpired:
            for pg in owned_groups(child.pid):
                try:os.killpg(pg,signal.SIGKILL)
                except ProcessLookupError:pass
            code=child.wait()
    seconds=time.monotonic()-start
    # Sum only sessions begun by this command, so the other policy's serial
    # governing work is not charged again as this policy's parent overhead.
    session_delta=sum(max(0,session_seconds(l)-before[q]) for q,l in ledgers.items())
    overhead=max(0,seconds-session_delta)
    receipt=dict(schema='WJ7WholeProcessReceiptV1',label=a.label,argv=argv,exit_code=code,
        failure=failure,seconds=seconds,whole_process_peak_rss_bytes=peak,peak_tree_threads=threads,
        unaccounted_parent_seconds=overhead,log=file_identity(out/'command.log'),intent=file_identity(out/'intent.json'),threads=THREADS)
    write_once(out/'receipt.json',receipt)
    for l in ledgers.values():l.finish(label,a.control_units,overhead,str(out/'receipt.json'))
    index=root/'command-receipts.json';receipts=read(index) if index.exists() else []
    receipts.append(receipt);atomic(index,receipts)
    print(json.dumps(receipt));raise SystemExit(code if code else 0)


if __name__=='__main__':main()
