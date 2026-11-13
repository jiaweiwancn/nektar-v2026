#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H8 bounded admission: one CP0, one16-step block, seven strict views."""
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


from run_wj_h7 import complete_operation,record_coverage_halt,coverage_quality

def admit_xml(path,cp0=False):
    xml=ET.parse(path).getroot()
    params={e.text.split('=')[0].strip():e.text.split('=')[1].strip()for e in xml.findall('./CONDITIONS/PARAMETERS/P')}
    assert float(params['RANSWallHeight'])==.15 and float(params['TimeStep'])==.0025,'only6mm/.0025'
    n=float(params['NumSteps']);assert n.is_integer()and n==(0 if cp0 else 16),'one fixed16-step trajectory'
    assert float(params['Time'])==0 and params['IO_CheckSteps']=='1'
    assert params['RANSWallTangentialOrder']=='12'and params['RANSWallNormalOrder']=='16'
    infos={e.get('PROPERTY'):e.get('VALUE')for e in xml.findall('./CONDITIONS/SOLVERINFO/I')}
    assert infos['RANSWallPimpleVerify']=='Off'and infos['RANSWallEvaluationReuse']=='ExactStepV1'and infos['RANSWallPimpleStopOnBound']=='On'
    if cp0:assert 'RANSWallRestart'not in infos
    else:assert infos.get('RANSWallRestart')not in (None,'None')
    return 'H8-CP0'if cp0 else 'H8-block',0 if cp0 else 16


def execute(root,label,argv,*,cwd=None,seconds=120,cp0=False,governing=False,strict_view=False):
    root=Path(root);assert root.name=='wj-h8-local-timestep-20261113'
    assert 0<seconds<=180
    root.mkdir(exist_ok=True);journal=root/'accounting.json'
    book=load(journal)if journal.exists()else dict(seconds_ceiling=1200,setup_seconds=60.,closure_reserve=180.,
        governing_attempts=0,cp0_reservations=0,entries={})
    assert label not in book['entries'],'write-once command label'
    spent=book['setup_seconds']+sum(v.get('seconds',v['reserved_seconds'])for v in book['entries'].values())
    assert spent+seconds+book['closure_reserve']<=1200,'RESOURCE_LIMIT: H8 remaining time'
    if cp0 or governing:
        assert not book.get('coverage_halt'),'expanded coverage obstruction stops new governing blocks'
        key,units=admit_xml(argv[1],cp0=cp0)
        seen=book.setdefault('native_reservations',[])
        assert key not in seen,'one fixed native reservation per branch, no retry'
        assert sum(v.get('units',0)for v in book['entries'].values())+units<=16,'16 governing-attempt ceiling'
        if cp0:
            assert book['cp0_reservations']<1,'only one fresh .0025 CP0'
            book['cp0_reservations']+=1
        seen.append(key)
    else:
        units=0
        if any(Path(str(v)).name=='IncNavierStokesSolver'for v in argv):raise ValueError('native solver requires H8 admission')
    if strict_view:
        assert book.get('strict_view_reservations',0)<7,'seven strict views maximum'
        book['strict_view_reservations']=book.get('strict_view_reservations',0)+1
    allocation=scientific_workers();assert allocation['workers']+10<=72,'RESOURCE_LIMIT workers'
    folder=root/'operations'/label;folder.mkdir(parents=True)
    request=dict(argv=list(map(str,argv)),cwd=str(cwd or root),threads=THREADS,reserved_seconds=seconds,cp0=cp0,units=units,strict_view=strict_view,
        inputs=input_files(argv),allocation=allocation,driver_sha256=sha(__file__))
    save(folder/'request.json',request);entry=dict(status='RESERVED',reserved_seconds=seconds,cp0=cp0,units=units,strict_view=strict_view)
    book['entries'][label]=entry;save(journal,book)
    env=dict(os.environ,**THREADS,TMPDIR=str(root/'tmp'),NEKTAR_RANS_TIMING_FILE=str(folder/'timing.jsonl'));(root/'tmp').mkdir(exist_ok=True)
    def limits():
        import resource
        resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3));resource.setrlimit(resource.RLIMIT_FSIZE,(8*1024**3,8*1024**3))
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
    save(folder/'loaded-maps.json',loaded);save(folder/'process-maps.json',process_maps)
    entry.update(status='RESOURCE_LIMIT'if reason else 'PASS'if code==0 else 'FAIL',seconds=time.monotonic()-start,exit=code,
        reason=reason,peak_tree_rss_bytes=peak,root_bytes=size(root),log_sha256=sha(folder/'command.log'))
    complete_operation(journal,label,entry);save(folder/'receipt.json',dict(request,**entry));print(label,entry['status'],entry['seconds'],flush=True)
    return code==0 and reason is None
