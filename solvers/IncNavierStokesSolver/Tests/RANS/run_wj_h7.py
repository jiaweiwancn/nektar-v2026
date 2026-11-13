#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Small H7-only bounded admission; six fixed blocks, three fresh CP0s, no retries."""
import json
import fcntl
import os
from pathlib import Path
import signal
import subprocess
import time
import xml.etree.ElementTree as ET
from run_wj_h4 import load,save,sha,input_files,size,THREADS
from wj6_resources import scientific_workers,tree_usage
from wj7_command import owned_groups


def update_journal(path,change):
    """Small shared journal mutation; children and parent use the same lock."""
    path=Path(path)
    with Path(str(path)+'.lock').open('a')as lock:
        fcntl.flock(lock,fcntl.LOCK_EX);latest=load(path);change(latest);save(path,latest)
    return latest


def complete_operation(path,label,entry):
    # Never publish the pre-spawn book: a scan child may have added a halt,
    # reservations or charges. Only this operation belongs to the parent.
    def change(latest):
        assert label in latest['entries'],'missing original reservation'
        latest['entries'][label].update(entry)
    return update_journal(path,change)


def record_coverage_halt(path,event):
    if event['status']not in ('OBSTRUCTED','UNRESOLVED'):raise ValueError('coverage halt requires a physical/sign failure')
    if event.get('role','accepted')!='accepted':raise ValueError('rejected candidate cannot become an accepted coverage halt')
    trajectory=event.get('trajectory',str(event.get('height',''))+'/'+str(event.get('branch','')))
    def change(latest):
        latest.setdefault('coverage_halt',event) # first DETECTED witness is immutable
        events=latest.setdefault('coverage_events',[])
        if event not in events:events.append(event)
        earliest=latest.setdefault('coverage_earliest',{})
        if trajectory not in earliest or event['time']<earliest[trajectory]['time']:earliest[trajectory]=event
    return update_journal(path,change)


def coverage_quality(records):
    accepted=[r['status']for r in records if r.get('role')=='accepted']
    if not accepted:return dict(coverage='NOT_RUN',quality='NOT_RUN')
    if any(s not in ('PASS','QUALITY_LOST','OBSTRUCTED','UNRESOLVED')for s in accepted):raise ValueError('unknown scan status')
    coverage='OBSTRUCTED'if 'OBSTRUCTED'in accepted else 'UNRESOLVED'if 'UNRESOLVED'in accepted else 'PASS'
    quality=coverage if coverage!='PASS'else 'QUALITY_LOST'if 'QUALITY_LOST'in accepted else 'PASS'
    return dict(coverage=coverage,quality=quality)


MATRIX=[(h,dt,n)for h in (5.5,6.,5.)for dt,n in ((.01,4),(.005,8))]
def admit_xml(path,cp0=False):
    xml=ET.parse(path).getroot()
    params={e.text.split('=')[0].strip():e.text.split('=')[1].strip()for e in xml.findall('./CONDITIONS/PARAMETERS/P')}
    h=float(params['RANSWallHeight'])*40;dt=float(params['TimeStep']);n=float(params['NumSteps'])
    assert n.is_integer()and (h,dt,int(n)if not cp0 else 8)in MATRIX,'outside fixed H7 matrix'
    assert float(params['Time'])==0 and params['IO_CheckSteps']=='1'
    infos={e.get('PROPERTY'):e.get('VALUE')for e in xml.findall('./CONDITIONS/SOLVERINFO/I')}
    assert infos['RANSWallPimpleVerify']=='Off'and infos['RANSWallEvaluationReuse']=='ExactStepV1'
    assert infos['RANSWallPimpleStopOnBound']=='On'
    if cp0:assert dt==.005 and n==0 and 'RANSWallRestart'not in infos,'fresh compatible fine CP0 only'
    else:assert infos.get('RANSWallRestart')not in (None,'None'),'H7 evolution restarts a compatible CP0'
    return f'{h:g}:{dt:g}:'+('CP0'if cp0 else 'block'),0 if cp0 else int(n)


def execute(root,label,argv,*,cwd=None,seconds=120,cp0=False,governing=False):
    root=Path(root);assert root.name=='wj-h7-neutral-transport-20261113'
    assert 0<seconds<=180
    root.mkdir(exist_ok=True);journal=root/'accounting.json'
    book=load(journal)if journal.exists()else dict(seconds_ceiling=1800,setup_seconds=60.,closure_reserve=240.,
        governing_attempts=0,cp0_reservations=0,entries={})
    assert label not in book['entries'],'write-once command label'
    spent=book['setup_seconds']+sum(v.get('seconds',v['reserved_seconds'])for v in book['entries'].values())
    assert spent+seconds+book['closure_reserve']<=1800,'RESOURCE_LIMIT: H7 remaining time'
    if cp0 or governing:
        assert not book.get('coverage_halt'),'expanded coverage obstruction stops new governing blocks'
        key,units=admit_xml(argv[1],cp0=cp0)
        seen=book.setdefault('native_reservations',[])
        assert key not in seen,'one fixed native reservation per branch, no retry'
        assert sum(v.get('units',0)for v in book['entries'].values())+units<=36,'36 governing-attempt ceiling'
        if cp0:
            assert book['cp0_reservations']<3,'only three fresh .005 CP0s'
            book['cp0_reservations']+=1
        seen.append(key)
    else:
        units=0
        if any(Path(str(v)).name=='IncNavierStokesSolver'for v in argv):raise ValueError('native solver requires H7 admission')
    allocation=scientific_workers();assert allocation['workers']+10<=72,'RESOURCE_LIMIT workers'
    folder=root/'operations'/label;folder.mkdir(parents=True)
    request=dict(argv=list(map(str,argv)),cwd=str(cwd or root),threads=THREADS,reserved_seconds=seconds,cp0=cp0,units=units,
        inputs=input_files(argv),allocation=allocation,driver_sha256=sha(__file__))
    save(folder/'request.json',request);entry=dict(status='RESERVED',reserved_seconds=seconds,cp0=cp0,units=units)
    book['entries'][label]=entry;save(journal,book)
    env=dict(os.environ,**THREADS,TMPDIR=str(root/'tmp'),NEKTAR_RANS_TIMING_FILE=str(folder/'timing.jsonl'));(root/'tmp').mkdir(exist_ok=True)
    def limits():
        import resource
        # A bounded caller may have imposed a smaller hard limit. Respect it;
        # a child must never try to raise the enclosing allocation's ceiling.
        for kind in (resource.RLIMIT_AS,resource.RLIMIT_FSIZE):
            hard=resource.getrlimit(kind)[1]
            limit=8*1024**3 if hard==resource.RLIM_INFINITY else min(8*1024**3,hard)
            resource.setrlimit(kind,(limit,limit))
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
    complete_operation(journal,label,entry);save(folder/'receipt.json',dict(request,**entry));print(label,entry['status'],entry['seconds'],flush=True)
    return code==0 and reason is None
