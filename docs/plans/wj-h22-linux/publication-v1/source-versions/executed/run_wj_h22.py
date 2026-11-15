#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H22 finite admission over existing journal, process and runtime primitives."""
import argparse,os,signal,subprocess,sys,time
from pathlib import Path
from run_wj_h4 import THREADS,load,save,sha,input_files
os.environ.update(THREADS)
from run_wj_h7 import update_journal,complete_operation
from run_wj_h10 import runtime_snapshot
from run_wj_h15 import size
from wj6_resources import scientific_workers,tree_usage
from wj7_command import owned_groups

ROOT='wj-h22-early-phase-replay-20261115'
ATTEMPTS=dict(C=1,M=2,F=4)
CAPS=dict(operation=180,build=900,qualification=120,target=600,target_readback=120,analysis=600,closure=300)

def usage(b):
    spent=lambda phases:sum(v.get('seconds',v['reserved_seconds'])for v in b['entries'].values()if v['phase']in phases)
    return dict(total=60+300+spent(set(CAPS)-{'closure'}),native=spent({'qualification','target','target_readback'}),
                build=spent({'build'}),analysis=spent({'analysis'}))

def halt(root,key,event):update_journal(Path(root)/'accounting.json',lambda b:b.setdefault(key,event))

def frozen(root):
    from wj_h15_observer import verify
    cfg=load(Path(root)/'execution-config.json')
    for v in cfg['files']:verify(v)
    return cfg

def command_level(root,argv,kind,attempts):
    base=[sys.executable,'-B',str(Path(__file__).with_name('study_wj_h22.py')),'--root',str(root),'--mode']
    assert kind in ('target','strict')
    for z,n in ATTEMPTS.items():
        if list(map(str,argv))==base+[kind,'--level',z]:
            assert attempts==(n if kind=='target'else 0);return z
    raise AssertionError('outside frozen seven-attempt/three-read matrix')

def execute(root,label,argv,seconds,phase='operation',workers=1,kind='passive',attempts=0):
    root=Path(root);assert root.name==ROOT and phase in CAPS and 0<seconds<=CAPS[phase]
    assert kind in ('passive','target','strict')and 1<=workers<=(8 if phase=='build'else 2 if kind=='target'else 1)
    journal=root/'accounting.json';folder=root/'operations'/label;folder.mkdir(exist_ok=False)
    allocation=scientific_workers();assert allocation['workers']+workers+6<=72
    def reserve(b):
        assert not b['closed']and label not in b['entries'];u=usage(b)
        assert u['total']+(0 if phase=='closure'else seconds)<=4500
        if phase in ('qualification','target','target_readback'):assert u['native']+seconds<=1800
        if phase in ('build','analysis'):assert u[phase]+seconds<=CAPS[phase]
        if kind!='passive':
            z=command_level(root,argv,kind,attempts);frozen(root)
            if kind=='target':
                assert not any(b.get(k)for k in ('target_halt','runtime_halt','coverage_halt','integrity_halt'))
                assert z==('C','M','F')[len(b['target_levels'])]and z not in b['target_levels']
                admission=load(root/f'admission-{z}.json');assert admission['status']=='PASS'and admission['command']==list(map(str,argv))
                assert b['target_reserved_attempts']+attempts<=7
                b['target_levels'].append(z);b['target_reserved_attempts']+=attempts
            else:
                assert z not in b['strict_reservations']and len(b['strict_reservations'])<3
                b['strict_reservations'].append(z)
        b['entries'][label]=dict(status='RESERVED',phase=phase,reserved_seconds=seconds,kind=kind,reserved_attempts=attempts)
    update_journal(journal,reserve)
    request=dict(argv=list(map(str,argv)),phase=phase,kind=kind,reserved_seconds=seconds,maximum_attempts=attempts,
        workers=workers,allocation=allocation,threads=THREADS,inputs=input_files(argv),driver_sha256=sha(__file__))
    save(folder/'request.json',request)
    env=dict(os.environ,**THREADS,PYTHONPATH=str(Path(__file__).resolve().parent),TMPDIR=str(root/'tmp'),MPLCONFIGDIR=str(root/'tmp/matplotlib'),
        NEKTAR_RANS_TIMING_FILE=str(folder/'timing.jsonl'),WJ_H22_OPERATION=str(folder),
        WJ_H14_PARENT='/home/jay/build/wj-h13-fullheight-paired-domain-20261113/publication-v1/paired-case')
    start=time.monotonic();peak=0;reason=None;groups=set();maps=[];seen=set()
    def mappings(pid):
        nonlocal reason
        try:
            files=sorted({str(Path(s.split(maxsplit=5)[5]).resolve())for s in Path(f'/proc/{pid}/maps').read_text().splitlines()
                if len(s.split(maxsplit=5))==6 and s.split(maxsplit=5)[5].startswith('/')and not s.endswith(' (deleted)')})
            v=runtime_snapshot(root,files)
            if v:
                key=(pid,tuple(files))
                if key not in seen:seen.add(key);maps.append(dict(v,pid=pid,elapsed=time.monotonic()-start))
                if v['status']!='PASS':reason='RUNTIME_CLOSURE_MISMATCH actual native-exec maps'
            for child in Path(f'/proc/{pid}/task/{pid}/children').read_text().split():mappings(int(child))
        except (FileNotFoundError,ProcessLookupError):pass
    with(folder/'command.log').open('x')as log:
        child=subprocess.Popen(['/usr/bin/time','-v',*map(str,argv)],cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        while child.poll()is None:
            try:
                rss,_=tree_usage(os.getpid());peak=max(peak,rss);groups|=owned_groups(child.pid);mappings(child.pid)
                if rss>(16 if phase=='build'else 8)*1024**3 or time.monotonic()-start>seconds:reason='RESOURCE_LIMIT RSS/time'
                source_bytes=sum(size(p)for p in root.iterdir()if p.is_dir()and p.name.startswith(('source','build','runtime')))
                if source_bytes>8*1024**3 or size(root/'science')>2*1024**3:reason='RESOURCE_LIMIT storage'
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
    count=attempts;unknown=kind=='target'
    if kind=='target'and(folder/'attempt-counts.json').exists():
        c=load(folder/'attempt-counts.json');assert 0<=c['actual_attempts']<=attempts
        count=c['actual_attempts'];unknown=False
    entry=dict(status='RESOURCE_LIMIT'if reason else'PASS'if code==0 else'FAIL',phase=phase,kind=kind,seconds=time.monotonic()-start,
        reserved_seconds=seconds,exit=code,reason=reason,peak_tree_rss_bytes=peak,reserved_attempts=attempts,charged_attempts=count,unknown_attempts=unknown,
        log_sha256=sha(folder/'command.log'))
    save(folder/'native-exec-samples.json',maps);complete_operation(journal,label,entry);save(folder/'receipt.json',dict(request,**entry))
    if kind=='target'and(code!=0 or reason):halt(root,'target_halt',dict(operation=label,reason=reason or'nonzero native/coordinator result',receipt=sha(folder/'receipt.json')))
    if reason and reason.startswith('RUNTIME'):halt(root,'runtime_halt',dict(operation=label,reason=reason))
    print(label,entry['status'],entry['seconds'],flush=True);return code==0 and not reason

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--label',required=True)
    p.add_argument('--seconds',type=float,required=True);p.add_argument('--phase',choices=CAPS,default='operation');p.add_argument('--workers',type=int,default=1)
    p.add_argument('--kind',choices=('passive','target','strict'),default='passive');p.add_argument('--attempts',type=int,default=0)
    p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args();v=a.command[1:]if a.command[:1]==['--']else a.command
    raise SystemExit(0 if execute(a.root,a.label,v,a.seconds,a.phase,a.workers,a.kind,a.attempts)else 1)
