#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H14 scope adapter over the qualified atomic-journal/resource helpers."""
import argparse
import os
from pathlib import Path
import signal
import subprocess
import time
from run_wj_h4 import load, save, sha, THREADS, input_files, size
os.environ.update(THREADS)  # before any optional admission module imports NumPy
from run_wj_h7 import update_journal, complete_operation
from wj6_resources import scientific_workers, tree_usage
from wj7_command import owned_groups
from run_wj_h10 import runtime_snapshot

ROOT='wj-h14-mapped-backend-20261113'
CAPS={'operation':1200,'build':2400,'qualification':1800,'target_cp0':600,'pilot':600,'closure':300}

def halt_admission(book,phase):
    assert not book['closed'],'closed H14 allocation'
    if phase in ('qualification','target_cp0','pilot'):
        assert not book.get('coverage_halt')and not book.get('runtime_halt'),'monotone H14 native/qualification halt'

def execute(root,label,argv,seconds,phase='operation',workers=1):
    root=Path(root);assert root.name==ROOT and phase in CAPS and 0<seconds<=CAPS[phase]
    assert 1<=workers<=(16 if phase=='build'else 1)
    book=root/'accounting.json';allocation=scientific_workers();assert allocation['workers']+workers+6<=72
    folder=root/'operations'/label;folder.mkdir(exist_ok=False)
    request=dict(argv=list(map(str,argv)),phase=phase,reserved_seconds=seconds,workers=workers,allocation=allocation,
        threads=THREADS,inputs=input_files(argv),driver_sha256=sha(__file__))
    def reserve(b):
        halt_admission(b,phase);assert label not in b['entries']
        charged=lambda pp:sum(v.get('seconds',v['reserved_seconds'])for v in b['entries'].values()if v['phase']in pp)
        assert 60+charged(set(CAPS)-{'closure'})+(0 if phase=='closure'else seconds)+300<=7200
        if phase!='operation':assert charged({phase})+seconds<=CAPS[phase]
        if phase in ('target_cp0','pilot'):
            assert (root/'target-admission.json').exists(),'target preflight remains closed'
            admission=load(root/'target-admission.json');assert admission['status']=='PASS'
            assert admission['height_mm']==6 and admission['state_size']==24688
            assert not any(v['phase']==phase for v in b['entries'].values()),'native reservation cannot be repeated'
            from wj_h14_execution import native_command
            native_command(admission,phase,argv)
            if phase=='target_cp0':b['target_CP0_reservations']=b.get('target_CP0_reservations',0)+1
            else:b['pilot_reserved_maximum_attempts']=4
        b['entries'][label]=dict(status='RESERVED',phase=phase,reserved_seconds=seconds)
    update_journal(book,reserve);save(folder/'request.json',request)
    env=dict(os.environ,**THREADS,TMPDIR=str(root/'tmp'),MPLCONFIGDIR=str(root/'tmp/matplotlib'),
        NEKTAR_RANS_TIMING_FILE=str(folder/'timing.jsonl'))
    rssLimit=(16 if phase=='build'else 8)*1024**3
    def limits():
        import resource
        resource.setrlimit(resource.RLIMIT_FSIZE,(2*1024**3,2*1024**3))
    start=time.monotonic();peak=0;reason=None;groups=set();samples=[];seen=set()
    def mappings(pid):
        nonlocal reason
        try:
            files=sorted({str(Path(s.split(maxsplit=5)[5]).resolve())for s in Path(f'/proc/{pid}/maps').read_text().splitlines()
                if len(s.split(maxsplit=5))==6 and s.split(maxsplit=5)[5].startswith('/')and not s.endswith(' (deleted)')})
            snapshot=runtime_snapshot(root,files)
            if snapshot:
                key=(pid,tuple(files))
                if key not in seen:snapshot.update(pid=pid,elapsed=time.monotonic()-start);samples.append(snapshot);seen.add(key)
                if snapshot['status']!='PASS':reason='RUNTIME_CLOSURE_MISMATCH actual native-exec maps'
            for child in Path(f'/proc/{pid}/task/{pid}/children').read_text().split():mappings(int(child))
        except (FileNotFoundError,ProcessLookupError):pass
    with(folder/'command.log').open('x')as log:
        p=subprocess.Popen(['/usr/bin/time','-v',*map(str,argv)],cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,preexec_fn=limits)
        while p.poll()is None:
            rss,n=tree_usage(p.pid);peak=max(peak,rss);groups|=owned_groups(p.pid);mappings(p.pid)
            if rss>rssLimit or time.monotonic()-start>seconds:reason='RESOURCE_LIMIT RSS/time'
            if (root/'science').exists()and size(root/'science')>2*1024**3:reason='RESOURCE_LIMIT scientific output'
            if (root/'build').exists()and size(root/'build')>8*1024**3:reason='RESOURCE_LIMIT build footprint'
            if scientific_workers()['workers']>72:reason='RESOURCE_LIMIT server workers'
            if reason:
                for g in groups:
                    try:os.killpg(g,signal.SIGTERM)
                    except ProcessLookupError:pass
                break
            time.sleep(.5)
        try:code=p.wait(timeout=5)
        except subprocess.TimeoutExpired:
            for g in groups:
                try:os.killpg(g,signal.SIGKILL)
                except ProcessLookupError:pass
            code=p.wait()
    entry=dict(status=('OBSTRUCTED'if reason.startswith('RUNTIME_')else'RESOURCE_LIMIT')if reason else'PASS'if code==0 else'FAIL',phase=phase,reserved_seconds=seconds,
        seconds=time.monotonic()-start,exit=code,reason=reason,peak_tree_rss_bytes=peak,log_sha256=sha(folder/'command.log'))
    if reason and reason.startswith('RUNTIME_'):
        def halt(b):b.setdefault('runtime_halt',dict(operation=label,reason=reason))
        update_journal(book,halt)
    save(folder/'native-exec-samples.json',samples)
    complete_operation(book,label,entry);save(folder/'receipt.json',dict(request,**entry));print(label,entry['status'],entry['seconds'],flush=True)
    return code==0 and not reason

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--label',required=True)
    p.add_argument('--seconds',type=float,required=True);p.add_argument('--phase',choices=CAPS,default='operation');p.add_argument('--workers',type=int,default=1)
    p.add_argument('argv',nargs=argparse.REMAINDER);a=p.parse_args();argv=a.argv[1:]if a.argv[:1]==['--']else a.argv
    raise SystemExit(0 if execute(a.root,a.label,argv,a.seconds,a.phase,a.workers)else 1)
