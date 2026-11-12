#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded independent H3 resource journal and bounded native command launcher."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import time
from wj6_resources import tree_usage, scientific_workers
from wj7_command import owned_groups

THREADS={k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
    'BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS')}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):
    # Native streams may emit binary64 negative zero as the JSON number -0.
    # Preserve its sign without converting integer step/geometry indices to float.
    return json.loads(Path(p).read_text(),parse_int=lambda v:-0.0 if v=='-0' else int(v))
def save(p,v):
    tmp=Path(str(p)+'.tmp');tmp.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n');os.replace(tmp,p)
def size(p):return int(subprocess.check_output(['du','-sb',str(p)],text=True).split()[0])
def input_files(argv):
    result={}
    for v in argv:
        try:is_file=Path(v).is_file()
        except OSError:is_file=False  # e.g. a long CTest regex, retained in argv
        if is_file:result[v]=dict(sha256=sha(v),bytes=Path(v).stat().st_size)
    return result


def execute(root,label,argv,cwd=None,units=0,seconds=600,expected=None,extra_env=None):
    root=Path(root).resolve();assert root.name=='wj-h3-20261112'
    assert 0<=units<=8 and 0<seconds<=900
    argv=list(map(str,argv));cwd=Path(cwd or root).resolve()
    out=root/'operations'/label;out.mkdir(parents=True,exist_ok=False)
    env=dict(os.environ,**THREADS);env.update(extra_env or {})
    env['NEKTAR_RANS_TIMING_FILE']=str(out/'timing.jsonl')
    request=dict(argv=argv,cwd=str(cwd),units=units,seconds_limit=seconds,expected_reason=expected,
        threads=THREADS,driver_sha256=sha(__file__),inputs=input_files(argv),loader=env.get('LD_LIBRARY_PATH',''))
    save(out/'request.json',request)
    journal=root/'accounting.json'
    with (root/'accounting.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        book=load(journal) if journal.exists() else dict(schema='WJH3ResourceJournalV1',unit_ceiling=44,
            seconds_ceiling=3600,closure_seconds=300,prior_profile_setup_seconds=120.,entries={})
        assert label not in book['entries'],'immutable H3 reservation exists'
        consumed=sum(v.get('seconds',v['reserved_seconds']) for v in book['entries'].values())
        charged=sum(v['units'] for v in book['entries'].values())
        assert consumed+book['closure_seconds']+book['prior_profile_setup_seconds']+seconds<=3600,'RESOURCE_LIMIT time before launch'
        assert charged+units<=44,'RESOURCE_LIMIT units before launch'
        allocation=scientific_workers();assert allocation['workers']+10<=72,'RESOURCE_LIMIT current worker allocation'
        entry=dict(units=units,reserved_seconds=seconds,status='RESERVED',allocation=allocation)
        book['entries'][label]=entry;save(journal,book)
    before=size(root);start=time.monotonic();groups=set();peak=threads=0;reason=None;lastSize=start
    def limits():
        import resource
        resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
        resource.setrlimit(resource.RLIMIT_FSIZE,(8*1024**3,8*1024**3))
    with (out/'command.log').open('x') as log:
        child=subprocess.Popen(['/usr/bin/time','-v',*argv],cwd=cwd,env=env,stdout=log,stderr=subprocess.STDOUT,
            start_new_session=True,preexec_fn=limits)
        while child.poll() is None:
            rss,n=tree_usage(child.pid);peak=max(peak,rss);threads=max(threads,n);groups|=owned_groups(child.pid)
            if rss>8*1024**3:reason='RESOURCE_LIMIT whole-process RSS'
            if time.monotonic()-start>seconds:reason='RESOURCE_LIMIT invocation timeout'
            if time.monotonic()-lastSize>=5:
                lastSize=time.monotonic();current=size(root)
                if current>16*1024**3 or current-before>8*1024**3:reason='RESOURCE_LIMIT output'
                if scientific_workers()['workers']>72:reason='RESOURCE_LIMIT global workers'
            if reason:
                for g in groups:
                    try:os.killpg(g,signal.SIGTERM)
                    except ProcessLookupError:pass
                break
            time.sleep(.5)
        try:code=child.wait(timeout=10)
        except subprocess.TimeoutExpired:
            for g in groups:
                try:os.killpg(g,signal.SIGKILL)
                except ProcessLookupError:pass
            code=child.wait()
    elapsed=time.monotonic()-start;log=out/'command.log';text=log.read_text(errors='replace')
    passed=(code==0 if expected is None else code not in (0,124) and expected in text)
    entry.update(status='RESOURCE_LIMIT' if reason else 'PASS' if passed else 'FAIL',seconds=elapsed,exit=code,
        reason=reason,peak_tree_rss_bytes=peak,peak_tree_threads=threads,root_bytes=size(root),output_bytes=size(root)-before,
        log_sha256=sha(log))
    with (root/'accounting.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX);book=load(journal);book['entries'][label]=entry;save(journal,book)
    save(out/'receipt.json',dict(request,**entry));print(json.dumps(dict(label=label,**entry)),flush=True)
    return passed and not reason


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--label',required=True)
    p.add_argument('--cwd',type=Path);p.add_argument('--units',type=int,default=0);p.add_argument('--seconds',type=int,default=600)
    p.add_argument('--expected-reason');p.add_argument('argv',nargs=argparse.REMAINDER);a=p.parse_args()
    argv=a.argv[1:] if a.argv[:1]==['--'] else a.argv
    if not argv:p.error('command required')
    raise SystemExit(0 if execute(a.root,a.label,argv,a.cwd,a.units,a.seconds,a.expected_reason) else 1)
if __name__=='__main__':main()
