#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H10 bounded native build and three-first-step admission; independent ledger."""
import json
import fcntl
import os
from pathlib import Path
import signal
import subprocess
import time
import xml.etree.ElementTree as ET
from run_wj_h4 import load,save,sha,input_files,size as _size,THREADS
from wj6_resources import scientific_workers,tree_usage
from wj7_command import owned_groups


from run_wj_h7 import complete_operation,record_coverage_halt,coverage_quality

def size(p):return _size(p)if Path(p).exists()else 0

def execute(root,label,argv,*,cwd=None,seconds=120,phase="science",kind="passive"):
    root=Path(root);assert root.name=='wj-h10-flat-inlet-20261113'
    assert 0<seconds<=2000 and phase in ('build','science') and kind in ('passive','fixture','cp0','governing','strict')
    root.mkdir(exist_ok=True);journal=root/'accounting.json'
    book=load(journal)if journal.exists()else dict(seconds_ceiling=3600,setup_seconds=60.,closure_reserve=300.,
        governing_attempts=0,cp0_reservations=0,strict_readbacks=0,builds=0,entries={})
    assert not book.get('closed'),'closed H10 allocation'
    assert not book.get('coverage_halt') and not book.get('h10_halt'),'H10 monotone scientific halt'
    assert label not in book['entries'],'write-once command label'
    charged=lambda v:v.get('seconds',v['reserved_seconds'])
    used=sum(charged(v)for v in book['entries'].values()if v['phase']==phase)+(60 if phase=='build'else 0)
    assert used+seconds<=(2400 if phase=='build'else 900),'RESOURCE_LIMIT phase allowance'
    key={'cp0':'cp0_reservations','governing':'governing_attempts','strict':'strict_readbacks'}.get(kind)
    if key:
        assert book[key]<{'cp0':3,'governing':3,'strict':6}[kind],'H10 invocation ceiling'
        book[key]+=1
    if phase=='build'and '--build'in argv:book['builds']+=1
    spent=book['setup_seconds']+sum(v.get('seconds',v['reserved_seconds'])for v in book['entries'].values())
    assert spent+seconds+book['closure_reserve']<=3600,'RESOURCE_LIMIT: H10 remaining time'
    units=0
    allocation=scientific_workers();assert allocation['workers']+6<=72,'RESOURCE_LIMIT workers'
    folder=root/'operations'/label;folder.mkdir(parents=True)
    request=dict(argv=list(map(str,argv)),cwd=str(cwd or root),threads=THREADS,reserved_seconds=seconds,units=int(kind=='governing'),kind=kind,phase=phase,
        inputs=input_files(argv),allocation=allocation,driver_sha256=sha(__file__))
    save(folder/'request.json',request);entry=dict(status='RESERVED',reserved_seconds=seconds,units=int(kind=='governing'),kind=kind,phase=phase)
    book['entries'][label]=entry;save(journal,book)
    env=dict(os.environ,**THREADS,TMPDIR=str(root/'tmp'),NEKTAR_RANS_TIMING_FILE=str(folder/'timing.jsonl'));(root/'tmp').mkdir(exist_ok=True)
    def limits():
        import resource
        resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));resource.setrlimit(resource.RLIMIT_FSIZE,(2*1024**3,2*1024**3))
    start=time.monotonic();peak=0;reason=None;groups=set();last=start;mapped=set();process_maps={}
    def mappings(pid):
        try:
            files=set()
            for line in Path(f'/proc/{pid}/maps').read_text().splitlines():
                fields=line.split(maxsplit=5)
                if len(fields)==6 and fields[5].startswith('/')and not fields[5].endswith(' (deleted)'):
                    file=str(Path(fields[5]).resolve());mapped.add(file);files.add(file)
            if files:
                old=process_maps.setdefault(str(pid),dict(command=Path(f'/proc/{pid}/cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace'),files=[]))
                old['files']=sorted(set(old['files'])|files)
            for ch in Path(f'/proc/{pid}/task/{pid}/children').read_text().split():mappings(int(ch))
        except (FileNotFoundError,ProcessLookupError):pass
    with(folder/'command.log').open('x')as f:
        child=subprocess.Popen(['/usr/bin/time','-v',*map(str,argv)],cwd=cwd or root,env=env,stdout=f,stderr=subprocess.STDOUT,
            start_new_session=True,preexec_fn=limits)
        while child.poll()is None:
            rss,n=tree_usage(child.pid);peak=max(peak,rss);groups|=owned_groups(child.pid);mappings(child.pid)
            if rss>4*1024**3 or time.monotonic()-start>seconds:reason='RESOURCE_LIMIT process memory/time'
            if time.monotonic()-last>=5:
                last=time.monotonic()
                if (sum(size(root/p)for p in ('source','build','runtime'))>8*1024**3 or
                    size(root/'science')>2*1024**3) or scientific_workers()['workers']>72:reason='RESOURCE_LIMIT output/workers'
            if reason:
                for g in groups:
                    try:os.killpg(g,signal.SIGTERM)
                    except ProcessLookupError:pass
                break
            time.sleep(.5)
        try:code=child.wait(timeout=5)
        except subprocess.TimeoutExpired:
            for g in groups:
                try:os.killpg(g,signal.SIGKILL)
                except ProcessLookupError:pass
            code=child.wait()
    loaded=[dict(path=p,bytes=Path(p).stat().st_size,sha256=sha(p))for p in sorted(mapped)if Path(p).is_file()]
    save(folder/'loaded-maps.json',loaded);save(folder/'process-maps.json',process_maps)
    entry.update(status='RESOURCE_LIMIT'if reason else 'PASS'if code==0 else 'FAIL',seconds=time.monotonic()-start,exit=code,
        reason=reason,peak_tree_rss_bytes=peak,root_bytes=size(root),build_source_runtime_bytes=sum(size(root/p)for p in ('source','build','runtime')),science_bytes=size(root/'science'),log_sha256=sha(folder/'command.log'))
    complete_operation(journal,label,entry);save(folder/'receipt.json',dict(request,**entry));print(label,entry['status'],entry['seconds'],flush=True)
    return code==0 and reason is None
