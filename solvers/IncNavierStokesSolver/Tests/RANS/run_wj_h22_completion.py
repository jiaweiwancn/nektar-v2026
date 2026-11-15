#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Only the linked H22 M2-verification/F-replay admission, using existing guards."""
import argparse,os,signal,subprocess,sys,time
from pathlib import Path
from run_wj_h4 import THREADS,load,save,sha,input_files
os.environ.update(THREADS)
from run_wj_h7 import update_journal,complete_operation
from run_wj_h10 import runtime_snapshot
from run_wj_h15 import size
from wj6_resources import scientific_workers,tree_usage
from wj7_command import owned_groups

ROOT='wj-h22-completion-20261115'
PARENT=Path('/home/jay/build/wj-h22-early-phase-replay-20261115')
NATIVE={'verify_M2':180,'strict_M2':120,'F':600,'strict_F':120}

def usage(b):
    charged=lambda phases:sum(v.get('seconds',v['reserved_seconds'])for v in b['entries'].values()if v['phase']in phases)
    return dict(total=b['carry']['total']+60+300+charged({'operation','analysis','native'}),
        native=b['carry']['native']+charged({'native'}),analysis=b['carry']['analysis']+charged({'analysis'}),build=b['carry']['build'])

def halt(root,key,event):update_journal(Path(root)/'accounting.json',lambda b:b.setdefault(key,event))

def frozen(root):
    from wj_h15_observer import verify
    c=load(root/'execution-config.json')
    for v in c['files']:verify(v)
    return c

def command_kind(root,argv):
    base=[sys.executable,'-B',str(Path(__file__).with_name('study_wj_h22_completion.py')),'--root',str(root),'--mode']
    for k in NATIVE:
        if list(map(str,argv))==base+[k]:return k
    raise AssertionError('outside one zero-advance verification, four F advances, two strict reads')

def reserve_native(b,kind):
    order=list(NATIVE);used=b['new_native_reservations']
    assert kind not in used and kind==order[len(used)],'native reservation is not repeatable or reorderable'
    assert not any(b.get(k)for k in ('native_halt','coverage_halt','integrity_halt','runtime_halt'))
    remaining=sum(NATIVE[k]for k in order[len(used):])
    u=usage(b);assert u['native']+1.25*remaining<=1800
    assert u['total']+1.25*remaining<=4500
    used.append(kind)
    if kind=='F':assert b['new_attempt_capacity']==0;b['new_attempt_capacity']=4

def execute(root,label,argv,seconds,phase='operation',kind='passive'):
    root=Path(root);assert root.name==ROOT and phase in ('operation','analysis','native')
    assert 0<seconds<={'operation':180,'analysis':180,'native':600}[phase]
    assert kind=='passive'or kind in NATIVE
    if kind!='passive':assert phase=='native'and seconds<=NATIVE[kind]and command_kind(root,argv)==kind
    assert phase!='native'or kind!='passive','no unnamed native operation'
    workers=2 if kind in ('verify_M2','F')else 1
    allocation=scientific_workers();assert allocation['workers']+workers+6<=72
    folder=root/'operations'/label;folder.mkdir(exist_ok=False);journal=root/'accounting.json'
    def reserve(b):
        assert not b['closed']and label not in b['entries'];u=usage(b)
        assert u['total']+seconds<=4500
        if phase=='analysis':assert u['analysis']+seconds<=600
        if kind!='passive':
            frozen(root)
            if kind=='F':
                from wj_h22_completion import require_M2_gate
                require_M2_gate(root)
            reserve_native(b,kind)
        b['entries'][label]=dict(status='RESERVED',phase=phase,kind=kind,reserved_seconds=seconds,
            reserved_attempts=4 if kind=='F'else 0)
    update_journal(journal,reserve)
    request=dict(argv=list(map(str,argv)),phase=phase,kind=kind,reserved_seconds=seconds,workers=workers,
        allocation=allocation,threads=THREADS,inputs=input_files(argv),driver_sha256=sha(__file__))
    save(folder/'request.json',request)
    env=dict(os.environ,**THREADS,PYTHONPATH=str(Path(__file__).resolve().parent),TMPDIR=str(root/'tmp'),
        MPLCONFIGDIR=str(root/'tmp/matplotlib'),NEKTAR_RANS_TIMING_FILE=str(folder/'timing.jsonl'),WJ_H22C_OPERATION=str(folder),
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
                if size(root/'source')+156584056>8*1024**3 or size(root/'science')+202392692>2*1024**3:reason='RESOURCE_LIMIT retained+new storage'
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
    count=4 if kind=='F'else 0;unknown=kind=='F'
    if(folder/'attempt-counts.json').exists():
        c=load(folder/'attempt-counts.json');count=c['actual_attempts'];assert 0<=count<=4 if kind=='F'else count==0;unknown=False
    entry=dict(status='RESOURCE_LIMIT'if reason else'PASS'if code==0 else'FAIL',phase=phase,kind=kind,
        seconds=time.monotonic()-start,reserved_seconds=seconds,exit=code,reason=reason,peak_tree_rss_bytes=peak,
        charged_attempts=count,unknown_attempts=unknown,log_sha256=sha(folder/'command.log'))
    save(folder/'native-exec-samples.json',maps);complete_operation(journal,label,entry);save(folder/'receipt.json',dict(request,**entry))
    if kind!='passive'and(code!=0 or reason):halt(root,'native_halt',dict(operation=label,kind=kind,reason=reason or'nonzero child',receipt=sha(folder/'receipt.json')))
    print(label,entry['status'],entry['seconds'],flush=True);return code==0 and not reason

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--label',required=True)
    p.add_argument('--seconds',type=float,required=True);p.add_argument('--phase',choices=('operation','analysis','native'),default='operation')
    p.add_argument('--kind',choices=('passive',*NATIVE),default='passive');p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args()
    v=a.command[1:]if a.command[:1]==['--']else a.command
    raise SystemExit(0 if execute(a.root,a.label,v,a.seconds,a.phase,a.kind)else 1)
