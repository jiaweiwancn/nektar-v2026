#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H25-only finite operation journal; serial science, owned-group limits."""
import argparse,os,signal,subprocess,sys,time
from pathlib import Path
from run_wj_h4 import THREADS,load,save,sha,input_files
os.environ.update(THREADS)
from run_wj_h7 import update_journal,complete_operation
from wj6_resources import scientific_workers,tree_usage
from wj7_command import owned_groups
from run_wj_h15 import size

def usage(b):
    result={k:sum(v.get('seconds',v['reserved_seconds'])for v in b['entries'].values()if v['phase']==k)
            for k in ('build','native','passive')}
    result['total']=b['setup_seconds']+b['closure_seconds']+sum(result.values())
    return result

def execute(root,label,argv,seconds,phase,attempts=0):
    root=Path(root);assert root.name=='wj-h25-production-20261116'
    assert phase in ('build','native','passive')and 0<seconds<={'build':900,'native':900,'passive':300}[phase]
    assert 0<=attempts<=64 and(not attempts or phase=='native')
    allocation=scientific_workers();assert allocation['workers']+(8 if phase=='build'else 2)+6<=72
    folder=root/'operations'/label;folder.mkdir(exist_ok=False)
    def reserve(b):
        assert not b['closed']and label not in b['entries']
        u=usage(b);assert u['total']+seconds<=b['ceilings']['total']
        assert u[phase]+seconds<=b['ceilings'][phase]
        if attempts:
            assert b['reserved_attempts']+attempts<=b['maximum_attempts']
            assert b['governor_blocks']<b['maximum_governor_blocks']
            b['reserved_attempts']+=attempts;b['governor_blocks']+=1
        b['entries'][label]=dict(status='RESERVED',phase=phase,reserved_seconds=seconds,reserved_attempts=attempts)
    update_journal(root/'accounting.json',reserve)
    req=dict(argv=list(map(str,argv)),phase=phase,reserved_seconds=seconds,attempts=attempts,
             threads=THREADS,allocation=allocation,inputs=input_files(argv),runner_sha256=sha(__file__))
    save(folder/'request.json',req)
    env=dict(os.environ,**THREADS,PYTHONPATH=str(Path(__file__).resolve().parent),TMPDIR=str(root/'tmp'),
        MPLCONFIGDIR=str(root/'tmp/matplotlib'),NEKTAR_RANS_TIMING_FILE=str(folder/'timing.jsonl'),
        WJ_H25_OPERATION=str(folder),WJ_H25_DEADLINE=str(time.monotonic()+seconds),
        WJ_H14_PARENT='/home/jay/build/wj-h13-fullheight-paired-domain-20261113/publication-v1/paired-case')
    start=time.monotonic();peak=0;groups=set();reason=None;last_disk=-100;maps=[];seen_maps=set()
    def native_maps(pid):
        nonlocal reason
        try:
            paths=sorted({str(Path(s.split(maxsplit=5)[5]).resolve())for s in Path(f'/proc/{pid}/maps').read_text().splitlines()
                if len(s.split(maxsplit=5))==6 and s.split(maxsplit=5)[5].startswith('/')and not s.endswith(' (deleted)')})
            if any(p.endswith('.elf')for p in paths)and any(p.endswith('/libc.so.6')for p in paths):
                key=(pid,tuple(paths))
                if key not in seen_maps:
                    assert len(maps)<16,'bounded native runtime-map catalog'
                    seen_maps.add(key);libc=[p for p in paths if p.endswith('/libc.so.6')]
                    okay=libc==['/home/jay/build/wj-h22-early-phase-replay-20261115/runtime/glibc/libc.so.6']
                    maps.append(dict(pid=pid,status='PASS'if okay else'FAIL',files=paths,elapsed=time.monotonic()-start))
                    if not okay:reason='MIXED_NATIVE_RUNTIME'
            for p in Path(f'/proc/{pid}/task/{pid}/children').read_text().split():native_maps(int(p))
        except (FileNotFoundError,ProcessLookupError):pass
    with(folder/'command.log').open('x')as log:
        child=subprocess.Popen(['/usr/bin/time','-v',*map(str,argv)],cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        while child.poll()is None:
            try:
                rss,_=tree_usage(os.getpid());peak=max(peak,rss);groups|=owned_groups(child.pid);native_maps(child.pid)
                if rss>(16 if phase=='build'else 8)*1024**3:reason='RSS_LIMIT'
                if time.monotonic()-start>seconds:reason='TIME_LIMIT'
                if time.monotonic()-start-last_disk>10:
                    last_disk=time.monotonic()-start
                    code=sum(size(p)for p in root.iterdir()if p.is_dir()and p.name.startswith(('source','build','runtime')))
                    output=sum(size(root/n)for n in ('science','publication-v1')if(root/n).exists())
                    if code>8*1024**3 or output>8*1024**3:reason='STORAGE_LIMIT'
                    if scientific_workers()['workers']>72:reason='HOST_WORKER_LIMIT'
            except Exception as e:reason='MONITOR_ERROR '+repr(e)
            if reason:
                for g in groups|{child.pid}:
                    try:os.killpg(g,signal.SIGTERM)
                    except ProcessLookupError:pass
                break
            time.sleep(.25)
        try:code=child.wait(timeout=5)
        except subprocess.TimeoutExpired:
            for g in groups|{child.pid}:
                try:os.killpg(g,signal.SIGKILL)
                except ProcessLookupError:pass
            code=child.wait()
    seconds_used=time.monotonic()-start
    event=dict(status='PASS'if code==0 and reason is None else'FAIL',phase=phase,
        reserved_seconds=seconds,seconds=seconds_used,exit=code,reason=reason,
        peak_tree_rss_bytes=peak,reserved_attempts=attempts,unknown_attempts=bool(attempts),
        log_sha256=sha(folder/'command.log'))
    counts=folder/'attempt-counts.json'
    if counts.exists():
        actual=load(counts);assert type(actual['actual_attempts'])is int and 0<=actual['actual_attempts']<=attempts
        event.update(actual);event['unknown_attempts']=not actual.get('known',True)
    save(folder/'native-exec-maps.json',maps)
    complete_operation(root/'accounting.json',label,event);save(folder/'receipt.json',dict(req,**event))
    print(label,event['status'],seconds_used,flush=True)
    return code==0 and reason is None

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--label',required=True)
    p.add_argument('--seconds',type=float,required=True);p.add_argument('--phase',choices=('build','native','passive'),required=True)
    p.add_argument('--attempts',type=int,default=0);p.add_argument('command',nargs=argparse.REMAINDER)
    a=p.parse_args();cmd=a.command[1:]if a.command[:1]==['--']else a.command
    raise SystemExit(0 if execute(a.root,a.label,cmd,a.seconds,a.phase,a.attempts)else 1)
