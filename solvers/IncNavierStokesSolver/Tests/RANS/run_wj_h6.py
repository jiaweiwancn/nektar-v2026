#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Small H6-only resource journal. No governing timestep launch is allowed."""
import json
import os
from pathlib import Path
import signal
import subprocess
import time
import xml.etree.ElementTree as ET
from run_wj_h4 import load,save,sha,input_files,size,THREADS
from wj6_resources import scientific_workers,tree_usage
from wj7_command import owned_groups


def execute(root,label,argv,*,cwd=None,seconds=120,cp0=False):
    root=Path(root);assert root.name=='wj-h6-neutral-projection-20261112'
    assert 0<seconds<=180
    root.mkdir(exist_ok=True);journal=root/'accounting.json'
    book=load(journal)if journal.exists()else dict(seconds_ceiling=1200,setup_seconds=40.,closure_reserve=180.,
        governing_attempts=0,cp0_reservations=0,entries={})
    assert label not in book['entries'],'write-once command label'
    spent=book['setup_seconds']+sum(v.get('seconds',v['reserved_seconds'])for v in book['entries'].values())
    assert spent+seconds+book['closure_reserve']<=1200,'RESOURCE_LIMIT: H6 remaining time'
    if cp0:
        assert book['cp0_reservations']<3,'H6 allows exactly three CP0 reservations, no retry'
        xml=ET.parse(argv[1]).getroot();params={e.text.split('=')[0].strip():e.text.split('=')[1].strip()for e in xml.findall('./CONDITIONS/PARAMETERS/P')}
        assert params['NumSteps']=='0' and float(params['Time'])==0,'H6 forbids governing timesteps'
        assert xml.find('./CONDITIONS/SOLVERINFO/I[@PROPERTY="RANSWallRestart"]')is None,'H6 forbids continuation'
        height=float(params['RANSWallHeight']);seen=book.setdefault('reserved_heights',[])
        assert height in (.125,.1375,.15)and height not in seen,'one H6 CP0 per prescribed height, no retry'
        seen.append(height)
        book['cp0_reservations']+=1
    elif any(Path(str(v)).name=='IncNavierStokesSolver'for v in argv):raise ValueError('solver requires CP0-only admission')
    allocation=scientific_workers();assert allocation['workers']+10<=72,'RESOURCE_LIMIT workers'
    folder=root/'operations'/label;folder.mkdir(parents=True)
    request=dict(argv=list(map(str,argv)),cwd=str(cwd or root),threads=THREADS,reserved_seconds=seconds,cp0=cp0,
        inputs=input_files(argv),allocation=allocation,driver_sha256=sha(__file__))
    save(folder/'request.json',request);entry=dict(status='RESERVED',reserved_seconds=seconds,cp0=cp0)
    book['entries'][label]=entry;save(journal,book)
    env=dict(os.environ,**THREADS,TMPDIR=str(root/'tmp'),NEKTAR_RANS_TIMING_FILE=str(folder/'timing.jsonl'));(root/'tmp').mkdir(exist_ok=True)
    def limits():
        import resource
        resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3));resource.setrlimit(resource.RLIMIT_FSIZE,(8*1024**3,8*1024**3))
    start=time.monotonic();peak=0;reason=None;groups=set();last=start;mapped=set()
    def mappings(pid):
        try:
            for line in Path(f'/proc/{pid}/maps').read_text().splitlines():
                fields=line.split(maxsplit=5)
                if len(fields)==6 and fields[5].startswith('/')and not fields[5].endswith(' (deleted)'):
                    mapped.add(str(Path(fields[5]).resolve()))
            for ch in Path(f'/proc/{pid}/task/{pid}/children').read_text().split():mappings(int(ch))
        except (FileNotFoundError,ProcessLookupError):pass
    with(folder/'command.log').open('x')as f:
        child=subprocess.Popen(['/usr/bin/time','-v',*map(str,argv)],cwd=cwd or root,env=env,stdout=f,stderr=subprocess.STDOUT,
            start_new_session=True,preexec_fn=limits)
        while child.poll()is None:
            rss,n=tree_usage(child.pid);peak=max(peak,rss);groups|=owned_groups(child.pid);mappings(child.pid)
            if rss>8*1024**3 or time.monotonic()-start>seconds:reason='RESOURCE_LIMIT process memory/time'
            if time.monotonic()-last>=5:
                last=time.monotonic()
                if size(root)>8*1024**3 or scientific_workers()['workers']>72:reason='RESOURCE_LIMIT output/workers'
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
    save(folder/'loaded-maps.json',loaded)
    entry.update(status='RESOURCE_LIMIT'if reason else 'PASS'if code==0 else 'FAIL',seconds=time.monotonic()-start,exit=code,
        reason=reason,peak_tree_rss_bytes=peak,root_bytes=size(root),log_sha256=sha(folder/'command.log'))
    book['entries'][label]=entry;save(journal,book);save(folder/'receipt.json',dict(request,**entry));print(label,entry['status'],entry['seconds'],flush=True)
    return code==0 and reason is None
