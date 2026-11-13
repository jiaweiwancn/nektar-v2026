#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H12 passive-only adapter to the existing locked resource helpers."""
import os
from pathlib import Path
import signal
import subprocess
import time
from run_wj_h4 import load,save,sha,THREADS,input_files
from run_wj_h7 import update_journal,complete_operation
from run_wj_h10 import size
from wj6_resources import scientific_workers,tree_usage
from wj7_command import owned_groups

ROOT='wj-h12-saved-phase-diagnosis-20261113'
FORBIDDEN={'IncNavierStokesSolver','RANSWallResolutionProbe','WJKernelProbe','cmake','make','ninja','mpirun'}

def passive_intent(argv):
    assert Path(argv[0]).name in ('python3','python'),'H12 accepts passive Python operations only'
    assert not any(Path(str(v)).name.removesuffix('.elf')in FORBIDDEN for v in argv),'forbidden native/build entrypoint'
    assert not set(argv)&{'--governing','--strict','--cp0','--restart','--build'},'H12 has zero native/strict/build scope'

def execute(root,label,argv,seconds=60,phase='operation'):
    root=Path(root);assert root.name==ROOT and phase in ('operation','science','closure')
    assert 0<seconds<=(360 if phase=='science'else 90)
    passive_intent(argv);journal=root/'accounting.json';allocation=scientific_workers();assert allocation['workers']+6<=72
    folder=root/'operations'/label;folder.mkdir(exist_ok=False)
    request=dict(argv=list(map(str,argv)),cwd=str(root),threads=THREADS,phase=phase,reserved_seconds=seconds,
                 allocation=allocation,inputs=input_files(argv),driver_sha256=sha(__file__))
    def reserve(b):
        assert not b.get('closed')and label not in b['entries']
        spent=lambda phases:sum(v.get('seconds',v['reserved_seconds'])for v in b['entries'].values()if v['phase']in phases)
        if phase!='closure':
            if b.get('h12_halt'):
                # H12 permits a charged purely passive correction, never a
                # native relaunch or erasure of the original stop receipt.
                assert b['h12_halt']['reason']=="KeyError('principal')",'unresolved scientific/provenance halt'
                if '--resume-from'in argv:
                    admission=b['passive_repair_admission'];assert not admission['used']
                    record=admission['receipt'];assert sha(record['path'])==record['sha256']
                    assert load(record['path'])['status']=='PASS'
                    admission['used']=True
                else:assert '--mode'in argv and argv[argv.index('--mode')+1]=='repair','monotone H12 halt'
            assert spent(('operation','science'))+seconds<=480
            if phase=='science':assert spent(('science',))+seconds<=360
        else:assert spent(('closure',))+seconds<=90
        b['entries'][label]=dict(status='RESERVED',reserved_seconds=seconds,phase=phase)
    update_journal(journal,reserve);save(folder/'request.json',request)
    env=dict(os.environ,**THREADS,TMPDIR=str(root/'tmp'))
    def limits():
        import resource
        resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));resource.setrlimit(resource.RLIMIT_FSIZE,(2*1024**3,2*1024**3))
    start=time.monotonic();peak=0;groups=set();reason=None
    with(folder/'command.log').open('x')as f:
        p=subprocess.Popen(['/usr/bin/time','-v',*map(str,argv)],cwd=root,env=env,stdout=f,stderr=subprocess.STDOUT,
                           start_new_session=True,preexec_fn=limits)
        while p.poll()is None:
            rss,n=tree_usage(os.getpid());peak=max(peak,rss);groups|=owned_groups(p.pid)
            if rss>4*1024**3 or time.monotonic()-start>seconds:reason='RESOURCE_LIMIT time/RSS'
            if size(root/'source')>1024**3 or size(root)-size(root/'source')>2*1024**3:reason='RESOURCE_LIMIT output'
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
    entry=dict(status='RESOURCE_LIMIT'if reason else 'PASS'if code==0 else 'FAIL',seconds=time.monotonic()-start,
               reserved_seconds=seconds,phase=phase,exit=code,reason=reason,peak_tree_rss_bytes=peak,
               source_bytes=size(root/'source'),new_output_bytes=size(root)-size(root/'source'),log_sha256=sha(folder/'command.log'))
    complete_operation(journal,label,entry);save(folder/'receipt.json',dict(request,**entry))
    print(label,entry['status'],entry['seconds'],flush=True)
    return code==0 and reason is None

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--label',required=True)
    p.add_argument('--seconds',type=float,default=60);p.add_argument('--phase',choices=('operation','science','closure'),default='operation')
    p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args()
    command=a.command[1:]if a.command[:1]==['--']else a.command
    raise SystemExit(0 if execute(a.root,a.label,command,a.seconds,a.phase)else 1)
