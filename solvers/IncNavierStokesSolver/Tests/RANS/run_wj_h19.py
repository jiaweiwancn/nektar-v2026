#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H19 admission adapter over the existing atomic journal/process helpers."""
import argparse
import os
from pathlib import Path
import signal
import subprocess
import time
from run_wj_h4 import THREADS,load,save,sha,input_files
os.environ.update(THREADS)
from run_wj_h7 import update_journal,complete_operation
from run_wj_h10 import runtime_snapshot
from run_wj_h15 import size
from wj6_resources import scientific_workers,tree_usage
from wj7_command import owned_groups

ROOT='wj-h19-target-timestep-20261115'
CAPS=dict(operation=1200,build=1200,qualification=1800,target=1800,target_readback=120,closure=300)

def native_selection(root,argv,kind,attempts):
    import sys
    assert kind in ('bootstrap','target','component')
    assert attempts==(2 if kind=='target'else 0)
    expected=[sys.executable,'-B',str(Path(__file__).with_name('study_wj_h19.py')),'--root',str(root),'--mode',kind]
    assert list(map(str,argv))==expected,'undeclared native command or matrix change'


def reserve_strict(root,scope):
    assert scope in ('small','target')
    def reserve(b):
        assert not b['closed'];used=b.setdefault('strict_reservations',[]);assert scope not in used,'strict slot already consumed';used.append(scope)
    update_journal(Path(root)/'accounting.json',reserve)

def halt(root,key,record):
    update_journal(Path(root)/'accounting.json',lambda b:b.setdefault(key,record))

def frozen(root):
    from wj_h15_observer import verify
    cfg=load(Path(root)/'execution-config.json')
    for v in cfg['files']:verify(v)
    return cfg

def execute(root,label,argv,seconds,phase='operation',workers=1,kind='passive',attempts=0):
    root=Path(root);assert root.name==ROOT and phase in CAPS and 0<seconds<=CAPS[phase]
    assert kind in ('passive','bootstrap','component','target')and 1<=workers<=(8 if phase=='build'else 2 if kind=='target'else 1)
    book=root/'accounting.json';folder=root/'operations'/label;folder.mkdir(exist_ok=False)
    allocation=scientific_workers();assert allocation['workers']+workers+6<=72
    def reserve(b):
        assert not b['closed']and label not in b['entries']
        charged=lambda phases:sum(v.get('seconds',v['reserved_seconds'])for v in b['entries'].values()if v['phase']in phases)
        assert 60+charged(set(CAPS)-{'closure'})+(0 if phase=='closure'else seconds)+300<=5400
        if phase in ('qualification','target','target_readback'):assert charged({'qualification','target','target_readback'})+seconds<=1800
        elif phase!='operation':assert charged({phase})+seconds<=CAPS[phase]
        if kind!='passive':
            native_selection(root,argv,kind,attempts)
            assert not any(b.get(k)for k in ('runtime_halt','qualification_halt','target_halt','coverage_halt','screen_halt')),'monotone scientific/execution halt'
            frozen(root)
        if kind=='bootstrap':
            assert attempts==0 and not b['bootstrap_reserved'],'coarse native CP0 slot consumed'
            b['bootstrap_reserved']=1
        if kind=='target':
            intent=load(root/'target-admission.json');assert intent['status']=='PASS'
            assert attempts==2 and not b.get('target_blocks'), 'one block/two attempts only; no retry'
            assert list(map(str,argv))==intent['command']
            from wj_h15_observer import verify
            for v in intent['inputs']:verify(v)
            b['target_blocks']=1;b['target_reserved_attempts']=2
        b['entries'][label]=dict(status='RESERVED',phase=phase,reserved_seconds=seconds,kind=kind,reserved_attempts=attempts,charged_attempts=attempts)
    update_journal(book,reserve)
    request=dict(argv=list(map(str,argv)),phase=phase,kind=kind,reserved_seconds=seconds,maximum_attempts=attempts,workers=workers,allocation=allocation,threads=THREADS,inputs=input_files(argv),driver=sha(__file__))
    save(folder/'request.json',request)
    env=dict(os.environ,**THREADS,PYTHONPATH=str(Path(__file__).resolve().parent),TMPDIR=str(root/'tmp'),MPLCONFIGDIR=str(root/'tmp/matplotlib'),NEKTAR_RANS_TIMING_FILE=str(folder/'timing.jsonl'),WJ_H19_OPERATION=str(folder),WJ_H14_PARENT='/home/jay/build/wj-h13-fullheight-paired-domain-20261113/publication-v1/paired-case')
    start=time.monotonic();peak=0;reason=None;groups=set();maps=[];seen=set();rss_limit=(16 if phase=='build'else 8)*1024**3
    def mappings(pid):
        nonlocal reason
        try:
            files=sorted({str(Path(s.split(maxsplit=5)[5]).resolve())for s in Path(f'/proc/{pid}/maps').read_text().splitlines()if len(s.split(maxsplit=5))==6 and s.split(maxsplit=5)[5].startswith('/')and not s.endswith(' (deleted)')})
            v=runtime_snapshot(root,files)
            if v:
                key=(pid,tuple(files))
                if key not in seen:seen.add(key);maps.append(dict(v,pid=pid,elapsed=time.monotonic()-start))
                if v['status']!='PASS':reason='RUNTIME_CLOSURE_MISMATCH actual native-exec maps'
            for child in Path(f'/proc/{pid}/task/{pid}/children').read_text().split():mappings(int(child))
        except (FileNotFoundError,ProcessLookupError):pass
    def limits():
        import resource
        resource.setrlimit(resource.RLIMIT_FSIZE,(2*1024**3,2*1024**3))
    with(folder/'command.log').open('x')as log:
        p=subprocess.Popen(['/usr/bin/time','-v',*map(str,argv)],cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,preexec_fn=limits)
        while p.poll()is None:
            try:
                rss,_=tree_usage(os.getpid());peak=max(peak,rss);groups|=owned_groups(p.pid);mappings(p.pid)
                if rss>rss_limit or time.monotonic()-start>seconds:reason='RESOURCE_LIMIT RSS/time'
                build_bytes=sum(size(q)if q.is_dir()else q.lstat().st_size for q in root.iterdir()if q.name.startswith(('source','build','runtime')))
                if build_bytes>8*1024**3 or size(root/'science')>2*1024**3:reason='RESOURCE_LIMIT storage'
                if scientific_workers()['workers']>72:reason='RESOURCE_LIMIT workers'
            except Exception as error:reason='OPERATIONAL_MONITOR_ERROR '+repr(error);groups.add(p.pid)
            if reason:
                for group in groups:
                    try:os.killpg(group,signal.SIGTERM)
                    except ProcessLookupError:pass
                break
            time.sleep(.5)
        try:code=p.wait(timeout=5)
        except subprocess.TimeoutExpired:
            for group in groups:
                try:os.killpg(group,signal.SIGKILL)
                except ProcessLookupError:pass
            code=p.wait()
    count=attempts
    if kind=='target'and(folder/'attempt-counts.json').exists():
        c=load(folder/'attempt-counts.json');assert 0<=c['actual_attempts']<=attempts
        from wj_h15_observer import verify
        r=load(verify(c['producer']));assert c['actual_attempts']==r['implicit_callbacks'];count=c['actual_attempts']
    entry=dict(status=('OBSTRUCTED'if reason and reason.startswith(('OPERATIONAL','RUNTIME'))else'RESOURCE_LIMIT')if reason else'PASS'if code==0 else'FAIL',phase=phase,kind=kind,reserved_seconds=seconds,seconds=time.monotonic()-start,exit=code,reason=reason,peak_tree_rss_bytes=peak,reserved_attempts=attempts,charged_attempts=count,unknown_attempts=kind=='target'and not(folder/'attempt-counts.json').exists(),log_sha256=sha(folder/'command.log'))
    def charge(b):
        if kind=='target':b['target_charged_attempts']=b.get('target_charged_attempts',0)+count
        if reason and reason.startswith('RUNTIME'):b.setdefault('runtime_halt',dict(operation=label,reason=reason))
        if kind!='passive'and(code!=0 or reason):b.setdefault('qualification_halt'if kind in ('bootstrap','component')else'target_halt',dict(operation=label,reason=reason or'nonzero coordinator/native result',native_log=entry['log_sha256']))
    update_journal(book,charge);save(folder/'native-exec-samples.json',maps);complete_operation(book,label,entry);save(folder/'receipt.json',dict(request,**entry))
    print(label,entry['status'],entry['seconds'],flush=True);return code==0 and not reason

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--label',required=True);p.add_argument('--seconds',type=float,required=True);p.add_argument('--phase',choices=CAPS,default='operation');p.add_argument('--workers',type=int,default=1);p.add_argument('--kind',choices=('passive','bootstrap','component','target'),default='passive');p.add_argument('--attempts',type=int,default=0);p.add_argument('argv',nargs=argparse.REMAINDER);a=p.parse_args();argv=a.argv[1:]if a.argv[:1]==['--']else a.argv
    raise SystemExit(0 if execute(a.root,a.label,argv,a.seconds,a.phase,a.workers,a.kind,a.attempts)else 1)
