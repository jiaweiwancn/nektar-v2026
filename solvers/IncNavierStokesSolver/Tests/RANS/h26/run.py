#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Separate finite H26 journal, bounded owned processes and live host admission."""
import argparse,datetime,json,os,signal,subprocess,sys,time
from pathlib import Path

THREADS={k:'1'for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS')}
LEGACY=Path('/home/jay/build/wj-h25-production-20261116/resumption-r1/source-v3/solvers/IncNavierStokesSolver/Tests/RANS')
sys.path.insert(0,str(LEGACY))
from wj6_resources import scientific_workers,tree_usage
from wj7_command import owned_groups
from run_wj_h15 import size
from wj_h15_observer import identity,atomic


def run(root,label,argv,phase,seconds,evidence=False):
    assert root.name=='wj-h26-native-flow-20261117'and phase in ('build','science','analysis','evidence')
    assert 0<seconds<=dict(build=900,science=1800,analysis=300,evidence=300)[phase]
    book=root/'accounting.json';a=json.loads(book.read_text());assert not a['closed']and label not in a['entries']
    workers=scientific_workers();reservation=16 if phase=='build'else 10 if phase=='science'else 2
    free=os.statvfs(root).f_bavail*os.statvfs(root).f_frsize
    d=root/'operations'/label;d.mkdir()
    req=dict(argv=argv,phase=phase,seconds_cap=seconds,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),allocation=workers,reservation=reservation,free_bytes=free,threads=THREADS,runner=identity(__file__))
    atomic(d/'request.json',req)
    if workers['workers']+reservation>72 or free<15*1024**3:
        event=dict(status='RESOURCE_NOT_ADMITTED',seconds=1.,phase=phase);a['entries'][label]=event;atomic(book,a);atomic(d/'receipt.json',dict(req,**event));print(label,event['status'],flush=True);return 1
    used=sum(e.get('seconds',e.get('reserved_seconds',0))for e in a['entries'].values())
    assert used+seconds+(0 if evidence or phase=='evidence'else 1800)<=21600
    a['entries'][label]=dict(status='RESERVED',phase=phase,reserved_seconds=seconds);atomic(book,a)
    start=time.monotonic();peak=0;reason=None;groups=set()
    env=dict(os.environ,**THREADS,TMPDIR=str(root/'tmp'),H26_ROOT=str(root),H26_OPERATION=str(d),PYTHONPATH=str(Path(__file__).resolve().parent)+':'+str(LEGACY))
    with(d/'command.log').open('xb')as log:
        p=subprocess.Popen(argv,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        while p.poll()is None:
            rss,_=tree_usage(os.getpid());peak=max(peak,rss);groups|=owned_groups(p.pid)
            if rss>(16 if phase=='build'else 8)*1024**3:reason='RSS_LIMIT'
            if time.monotonic()-start>seconds:reason='TIME_LIMIT'
            if scientific_workers()['workers']>72:reason='HOST_WORKER_LIMIT'
            if os.statvfs(root).f_bavail*os.statvfs(root).f_frsize<15*1024**3:reason='FREE_LIMIT'
            if size(root/'science')>8*1024**3:reason='SCIENCE_STORAGE_LIMIT'
            if sum(size(v)for v in root.iterdir()if v.is_dir()and v.name.startswith(('source','build','runtime')))>8*1024**3:reason='CODE_STORAGE_LIMIT'
            if(d/'command.log').stat().st_size>32*1024**2:reason='LOG_LIMIT'
            if reason:
                for g in groups|{p.pid}:
                    try:os.killpg(g,signal.SIGTERM)
                    except ProcessLookupError:pass
                break
            time.sleep(.3)
        try:code=p.wait(timeout=5)
        except subprocess.TimeoutExpired:
            for g in groups|{p.pid}:
                try:os.killpg(g,signal.SIGKILL)
                except ProcessLookupError:pass
            code=p.wait()
    event=dict(status='PASS'if code==0 and reason is None else'FAIL',phase=phase,seconds=time.monotonic()-start,exit=code,reason=reason,peak_tree_rss_bytes=peak,log=identity(d/'command.log'))
    a=json.loads(book.read_text());a['entries'][label]=event;atomic(book,a);atomic(d/'receipt.json',dict(req,**event));print(label,event['status'],event['seconds'],flush=True);return int(event['status']!='PASS')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--label',required=True);p.add_argument('--phase',required=True);p.add_argument('--seconds',type=float,required=True);p.add_argument('--evidence',action='store_true');p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args();raise SystemExit(run(a.root,a.label,a.command[1:]if a.command[:1]==['--']else a.command,a.phase,a.seconds,a.evidence))
