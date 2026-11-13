#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H13-only offline admission using the established journal/tree-resource helpers."""
import argparse
import os
from pathlib import Path
import signal
import subprocess
import time
from run_wj_h4 import load,save,sha,THREADS,input_files,size
from run_wj_h7 import update_journal,complete_operation
from wj6_resources import scientific_workers,tree_usage
from wj7_command import owned_groups

ROOT='wj-h13-fullheight-paired-domain-20261113'
FORBIDDEN={'IncNavierStokesSolver','RANSWallResolutionProbe','WJKernelProbe','WJTerrainCoupledProbe','cmake','make','ninja','mpirun'}

def offline_intent(argv):
    assert Path(argv[0]).name in ('python3','python'),'H13 admits only offline Python operations'
    assert not any(Path(str(v)).name.removesuffix('.elf')in FORBIDDEN for v in argv),'zero native/build scope'
    assert not set(argv)&{'--execute','--native','--governing','--cp0','--restart','--strict','--build'},'H13 native work is NOT_RUN'

def execute(root,label,argv,seconds=60,phase='operation'):
    root=Path(root);assert root.name==ROOT and phase in ('operation','science','closure')
    assert 0<seconds<=(360 if phase=='science'else 90);offline_intent(argv)
    book=root/'accounting.json';allocation=scientific_workers();assert allocation['workers']+6<=72
    folder=root/'operations'/label;folder.mkdir(exist_ok=False)
    request=dict(argv=list(map(str,argv)),cwd=str(root),threads=THREADS,phase=phase,reserved_seconds=seconds,
                 inputs=input_files(argv),allocation=allocation,driver_sha256=sha(__file__))
    def reserve(b):
        assert not b['closed']and label not in b['entries']and not b.get('coverage_halt')
        used=lambda pp:sum(v.get('seconds',v['reserved_seconds'])for v in b['entries'].values()if v['phase']in pp)
        if phase=='closure':assert used(('closure',))+seconds<=90
        else:
            assert used(('operation','science'))+seconds<=480
            assert 30+used(('operation','science'))+seconds+90<=600
            if phase=='science':assert used(('science',))+seconds<=360
        b['entries'][label]=dict(status='RESERVED',phase=phase,reserved_seconds=seconds)
    update_journal(book,reserve);save(folder/'request.json',request)
    env=dict(os.environ,**THREADS,TMPDIR=str(root/'tmp'),MPLCONFIGDIR=str(root/'tmp/matplotlib'))
    def limits():
        import resource
        resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));resource.setrlimit(resource.RLIMIT_FSIZE,(256*1024**2,256*1024**2))
    start=time.monotonic();peak=0;reason=None;groups=set()
    with(folder/'command.log').open('x')as f:
        p=subprocess.Popen(['/usr/bin/time','-v',*map(str,argv)],cwd=root,env=env,stdout=f,stderr=subprocess.STDOUT,start_new_session=True,preexec_fn=limits)
        while p.poll()is None:
            rss,n=tree_usage(os.getpid());peak=max(peak,rss);groups|=owned_groups(p.pid)
            if rss>2*1024**3 or time.monotonic()-start>seconds or size(root)>256*1024**2:reason='RESOURCE_LIMIT RSS/time/output'
            if reason:
                for g in groups:
                    try:os.killpg(g,signal.SIGTERM)
                    except ProcessLookupError:pass
                break
            time.sleep(.25)
        try:code=p.wait(timeout=5)
        except subprocess.TimeoutExpired:
            for g in groups:
                try:os.killpg(g,signal.SIGKILL)
                except ProcessLookupError:pass
            code=p.wait()
    entry=dict(status='RESOURCE_LIMIT'if reason else 'PASS'if code==0 else 'FAIL',phase=phase,reserved_seconds=seconds,
        seconds=time.monotonic()-start,exit=code,reason=reason,peak_tree_rss_bytes=peak,output_bytes=size(root),log_sha256=sha(folder/'command.log'))
    complete_operation(book,label,entry);save(folder/'receipt.json',dict(request,**entry));print(label,entry['status'],entry['seconds'],flush=True)
    return code==0 and not reason

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--label',required=True);p.add_argument('--seconds',type=float,default=60)
    p.add_argument('--phase',choices=('operation','science','closure'),default='operation');p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args()
    argv=a.command[1:]if a.command[:1]==['--']else a.command
    raise SystemExit(0 if execute(a.root,a.label,argv,a.seconds,a.phase)else 1)
