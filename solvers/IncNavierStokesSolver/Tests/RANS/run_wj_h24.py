#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Finite H24 M/C/F continuation with coherent imported-source admission."""
import argparse,os,signal,subprocess,sys,time
from pathlib import Path
from run_wj_h4 import THREADS,load,save,sha,input_files
os.environ.update(THREADS)
from run_wj_h7 import update_journal,complete_operation
from run_wj_h10 import runtime_snapshot
from run_wj_h15 import size
from wj6_resources import scientific_workers,tree_usage
from wj7_command import owned_groups

ROOT='wj-h24-flat-continuation-20261115'
PARENT=Path('/home/jay/build/wj-h22-early-phase-replay-20261115')
NATIVE={'M':1050,'strict_M':120,'C':900,'strict_C':120,'F':1650,'strict_F':120}

def usage(b):
    charged=lambda phases:sum(v.get('seconds',v['reserved_seconds'])for v in b['entries'].values()if v['phase']in phases)
    return dict(total=60+300+charged({'passive','native'}),native=charged({'native'}),passive=charged({'passive'}))

def halt(root,key,event):update_journal(Path(root)/'accounting.json',lambda b:b.setdefault(key,event))

def frozen(root):
    from wj_h15_observer import verify
    c=load(root/'execution-config.json')
    for v in c['files']:verify(v)
    from wj_h23_continuation import active_execution_bindings
    active_execution_bindings(c['compatibility'])
    return c

def command_kind(root,argv):
    base=[sys.executable,'-B',str(Path(__file__).with_name('study_wj_h24.py')),'--root',str(root),'--mode']
    for k in NATIVE:
        if list(map(str,argv))==base+[k]:return k
    raise AssertionError('outside the fixed H24 M/C/F governor and three strict-read matrix')

def reserve_native(b,kind):
    order=list(NATIVE);used=b['native_reservations']
    assert kind not in used,'native reservation is not repeatable'
    if kind.startswith('strict_'):
        z=kind.split('_')[1]
        assert used and used[-1]==z,'strict read only after that branch; no other governor'
    else:
        assert kind==order[len(used)],'frozen M -> strict -> C -> strict -> F -> strict order'
        assert not any(b.get(k)for k in ('native_halt','coverage_halt','integrity_halt','runtime_halt','screen_halt','resource_halt'))
    assert usage(b)['native']+NATIVE[kind]<=4200
    used.append(kind)
    if not kind.startswith('strict_'):
        from wj_h23_continuation import ADDITIONAL
        b['reserved_attempts']+=ADDITIONAL[kind]
        assert b['reserved_attempts']<=28

def execute(root,label,argv,seconds,phase='passive',kind='passive'):
    root=Path(root);assert root.name==ROOT and phase in ('passive','native')
    assert 0<seconds<={'passive':300,'native':1650}[phase]
    assert kind=='passive'or kind in NATIVE
    if kind!='passive':assert phase=='native'and seconds<=NATIVE[kind]and command_kind(root,argv)==kind
    assert phase!='native'or kind!='passive','no unnamed native operation'
    workers=2 if kind in ('M','C','F')else 1
    allocation=scientific_workers();assert allocation['workers']+workers+6<=72
    folder=root/'operations'/label;folder.mkdir(exist_ok=False);journal=root/'accounting.json'
    def reserve(b):
        assert not b['closed']and label not in b['entries'];u=usage(b)
        assert u['total']+seconds<=6000
        if phase=='passive':assert u['passive']+seconds<=900
        if kind!='passive':
            frozen(root)
            if not kind.startswith('strict_'):
                from study_wj_h24 import admission
                admission(root,kind)
            reserve_native(b,kind)
        b['entries'][label]=dict(status='RESERVED',phase=phase,kind=kind,reserved_seconds=seconds,
            reserved_attempts={'M':8,'C':4,'F':16}.get(kind,0))
    update_journal(journal,reserve)
    request=dict(argv=list(map(str,argv)),phase=phase,kind=kind,reserved_seconds=seconds,workers=workers,
        allocation=allocation,threads=THREADS,inputs=input_files(argv),driver_sha256=sha(__file__))
    save(folder/'request.json',request)
    env=dict(os.environ,**THREADS,PYTHONPATH=str(Path(__file__).resolve().parent),TMPDIR=str(root/'tmp'),
        MPLCONFIGDIR=str(root/'tmp/matplotlib'),NEKTAR_RANS_TIMING_FILE=str(folder/'timing.jsonl'),WJ_H24_OPERATION=str(folder),WJ_H23_OPERATION=str(folder),WJ_H23_DEADLINE=str(time.monotonic()+seconds),WJ_H23_OPERATION_START=str(time.monotonic()),
        WJ_H14_PARENT='/home/jay/build/wj-h13-fullheight-paired-domain-20261113/publication-v1/paired-case')
    start=time.monotonic();peak=0;reason=None;groups=set();maps=[];seen=set()
    def mappings(pid):
        nonlocal reason
        try:
            files=sorted({str(Path(s.split(maxsplit=5)[5]).resolve())for s in Path(f'/proc/{pid}/maps').read_text().splitlines()
                if len(s.split(maxsplit=5))==6 and s.split(maxsplit=5)[5].startswith('/')and not s.endswith(' (deleted)')})
            v=runtime_snapshot(PARENT,files) # deliberately reuse the qualified closure IN PLACE
            if v:
                key=pid,tuple(files)
                if key not in seen:seen.add(key);maps.append(dict(v,pid=pid,elapsed=time.monotonic()-start))
                if v['status']!='PASS':reason='RUNTIME_CLOSURE_MISMATCH actual native-exec maps'
            for child in Path(f'/proc/{pid}/task/{pid}/children').read_text().split():mappings(int(child))
        except (FileNotFoundError,ProcessLookupError):pass
    with(folder/'command.log').open('x')as log:
        child=subprocess.Popen(['/usr/bin/time','-v',*map(str,argv)],cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        while child.poll()is None:
            try:
                rss,_=tree_usage(os.getpid());peak=max(peak,rss);groups|=owned_groups(child.pid);mappings(child.pid)
                if rss>8*1024**3 or time.monotonic()-start>seconds:reason='RESOURCE_LIMIT RSS/time'
                source_bytes=sum(size(p)for p in root.iterdir()if p.is_dir()and p.name.startswith('source'))
                output_bytes=sum(size(root/n)for n in ('science','publication-v1')if(root/n).exists())
                if source_bytes>8*1024**3 or output_bytes>2*1024**3:reason='RESOURCE_LIMIT new storage'
                if scientific_workers()['workers']>72:reason='RESOURCE_LIMIT workers'
            except Exception as e:reason='OPERATIONAL_MONITOR_ERROR '+repr(e);groups.add(child.pid)
            if reason:
                for group in groups:
                    try:os.killpg(group,signal.SIGTERM)
                    except ProcessLookupError:pass
                break
            time.sleep(.5)
        try:code=child.wait(timeout=5)
        except subprocess.TimeoutExpired:
            for group in groups:
                try:os.killpg(group,signal.SIGKILL)
                except ProcessLookupError:pass
            code=child.wait()
    count={'M':8,'C':4,'F':16}.get(kind,0);unknown=kind in ('M','C','F')
    if(folder/'attempt-counts.json').exists():
        c=load(folder/'attempt-counts.json');count=c['actual_attempts'];assert 0<=count<={'M':8,'C':4,'F':16}.get(kind,0);unknown=False
    entry=dict(status=('OBSTRUCTED'if reason.startswith(('OPERATIONAL','RUNTIME'))else'RESOURCE_LIMIT')if reason else'PASS'if code==0 else'FAIL',phase=phase,kind=kind,
        seconds=time.monotonic()-start,reserved_seconds=seconds,exit=code,reason=reason,peak_tree_rss_bytes=peak,
        charged_attempts=count,unknown_attempts=unknown,log_sha256=sha(folder/'command.log'))
    save(folder/'native-exec-samples.json',maps);complete_operation(journal,label,entry);save(folder/'receipt.json',dict(request,**entry))
    if kind!='passive'and(code!=0 or reason):halt(root,'native_halt',dict(operation=label,kind=kind,reason=reason or'nonzero child',receipt=sha(folder/'receipt.json')))
    print(label,entry['status'],entry['seconds'],flush=True);return code==0 and not reason

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--label',required=True)
    p.add_argument('--seconds',type=float,required=True);p.add_argument('--phase',choices=('passive','native'),default='passive')
    p.add_argument('--kind',choices=('passive',*NATIVE),default='passive');p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args()
    v=a.command[1:]if a.command[:1]==['--']else a.command
    raise SystemExit(0 if execute(a.root,a.label,v,a.seconds,a.phase,a.kind)else 1)
