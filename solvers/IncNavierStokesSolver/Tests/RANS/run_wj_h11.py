#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H11-only reservations on the existing locked journal/resource primitives."""
import os
from pathlib import Path
import signal
import subprocess
import time
from run_wj_h4 import load,save,sha,input_files,THREADS
from run_wj_h10 import size,runtime_snapshot
from run_wj_h7 import update_journal,complete_operation,record_coverage_halt
from wj6_resources import scientific_workers,tree_usage
from wj7_command import owned_groups
from wj_h11_analysis import admit_xml,native_reserve,VIEWS

H10=Path('/home/jay/build/wj-h10-flat-inlet-20261113')
ROOT='wj-h11-compatible-common-time-20261113'

def intent(root,argv,kind):
    if kind=='passive':return None,0
    assert Path(argv[0])==H10/'runtime'/('IncNavierStokesSolver'if kind=='governing'else 'RANSWallResolutionProbe')
    assert sha(root/'science/input-intent.json')==load(root/'accounting.json')['intent_sha256'],'frozen requested intent changed'
    plan=load(root/'science/input-intent.json');case=next(v for v in plan['cases']if v['xml']['path']==str(argv[1]))
    for v in [case['parent']['producer'],case['parent']['checkpoint'],*plan['runtime_files']]:assert Path(v['path']).stat().st_size==v['bytes']and sha(v['path'])==v['sha256'],'stale checkpoint/producer/runtime'
    for v in plan['runtime']['executables']:
        for key in ('ELF','wrapper'):assert sha(v[key]['path'])==v[key]['sha256']
    assert sha(argv[1])==case['xml']['sha256'];units=admit_xml(argv[1],case['branch'],case['original_xml']['path'])
    if kind=='governing':return 'block:'+case['branch'],units
    request=load(argv[2]);assert request['operation']=='open_checks'and request['native_report']==case['report']
    step=request['selected_absolute_step'];assert step in VIEWS[case['branch']]
    return 'strict:'+case['branch']+':'+str(step),0

def execute(root,label,argv,*,cwd=None,seconds=120,kind='passive'):
    root=Path(root);assert root.name==ROOT and kind in ('passive','governing','strict')and 0<seconds<=240
    journal=root/'accounting.json';assert journal.is_file(),'H11 requires frozen setup/ledger'
    key,units=intent(root,argv,kind);allocation=scientific_workers();assert allocation['workers']+6<=72
    folder=root/'operations'/label;folder.mkdir(parents=True,exist_ok=False)
    request=dict(argv=list(map(str,argv)),cwd=str(cwd or root),threads=THREADS,reserved_seconds=seconds,units=units,kind=kind,
        inputs=input_files(argv),allocation=allocation,driver_sha256=sha(__file__))
    entry=dict(status='RESERVED',reserved_seconds=seconds,units=units,kind=kind)
    def reserve(book):
        assert not book.get('closed')and label not in book['entries']
        spent=sum(v.get('seconds',v['reserved_seconds'])for v in book['entries'].values())
        assert spent+seconds<=1200 and 60+spent+seconds+240<=1500,'RESOURCE_LIMIT H11 allocation'
        if key:
            native_reserve(book,key,units)
            field='governing_blocks'if kind=='governing'else 'strict_readbacks';limit=3 if kind=='governing'else 11
            assert book[field]<limit;book[field]+=1
        book['entries'][label]=entry
    update_journal(journal,reserve);save(folder/'request.json',request)
    env=dict(os.environ,**THREADS,TMPDIR=str(root/'tmp'),NEKTAR_RANS_TIMING_FILE=str(folder/'timing.jsonl'));(root/'tmp').mkdir(exist_ok=True)
    def limits():
        import resource
        resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));resource.setrlimit(resource.RLIMIT_FSIZE,(3*1024**3,3*1024**3))
    start=time.monotonic();peak=0;reason=None;groups=set();last=start;mapped=set();process_maps={};samples=[]
    allowed={v['path']:v for v in load(root/'science/input-intent.json')['runtime_files']}if(root/'science/input-intent.json').exists()else{}
    def mappings(pid):
        nonlocal reason
        try:
            files=set()
            for line in Path(f'/proc/{pid}/maps').read_text().splitlines():
                fields=line.split(maxsplit=5)
                if len(fields)==6 and fields[5].startswith('/')and not fields[5].endswith(' (deleted)'):files.add(str(Path(fields[5]).resolve()))
            mapped.update(files);sample=runtime_snapshot(H10,sorted(files))
            if sample:
                sample.update(pid=pid,elapsed=time.monotonic()-start);samples.append(sample)
                if sample['status']!='PASS':reason='RUNTIME_CLOSURE_MISMATCH'
                for p in files:
                    if '.so'in Path(p).name or p.endswith('.elf'):
                        if p not in allowed or sha(p)!=allowed[p]['sha256']:reason='RUNTIME_CLOSURE_MISMATCH';sample['status']='OBSTRUCTED'
            if files:
                old=process_maps.setdefault(str(pid),dict(command=Path(f'/proc/{pid}/cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace'),files=[]));old['files']=sorted(set(old['files'])|files)
            for ch in Path(f'/proc/{pid}/task/{pid}/children').read_text().split():mappings(int(ch))
        except (FileNotFoundError,ProcessLookupError):pass
    with(folder/'command.log').open('x')as f:
        child=subprocess.Popen(['/usr/bin/time','-v',*map(str,argv)],cwd=cwd or root,env=env,stdout=f,stderr=subprocess.STDOUT,start_new_session=True,preexec_fn=limits)
        while child.poll()is None:
            rss,n=tree_usage(os.getpid());peak=max(peak,rss);groups|=owned_groups(child.pid);mappings(child.pid)
            if rss>4*1024**3 or time.monotonic()-start>seconds:reason='RESOURCE_LIMIT process memory/time'
            if time.monotonic()-last>=5:
                last=time.monotonic()
                if sum(size(root/p)for p in ('source','runtime'))>1024**3 or size(root/'science')>3*1024**3 or scientific_workers()['workers']>72:reason='RESOURCE_LIMIT output/workers'
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
    save(folder/'loaded-maps.json',[dict(path=p,bytes=Path(p).stat().st_size,sha256=sha(p))for p in sorted(mapped)if Path(p).is_file()])
    save(folder/'process-maps.json',process_maps);save(folder/'native-exec-samples.json',samples)
    entry.update(status=('OBSTRUCTED'if reason=='RUNTIME_CLOSURE_MISMATCH'else 'RESOURCE_LIMIT')if reason else 'PASS'if code==0 else 'FAIL',
        seconds=time.monotonic()-start,exit=code,reason=reason,peak_tree_rss_bytes=peak,source_runtime_bytes=sum(size(root/p)for p in ('source','runtime')),science_bytes=size(root/'science'),log_sha256=sha(folder/'command.log'))
    complete_operation(journal,label,entry)
    if kind!='passive'and(code!=0 or reason):
        update_journal(journal,lambda b:b.setdefault('h11_halt',dict(status=entry['status'],operation=label,reason=reason,exit=code,log_sha256=entry['log_sha256'])))
    save(folder/'receipt.json',dict(request,**entry));print(label,entry['status'],entry['seconds'],flush=True)
    return code==0 and not reason

def command(root,label,argv,seconds=120,**kwargs):
    start=time.monotonic();ok=execute(root,label,argv,seconds=seconds,**kwargs)
    def finish(b):b['entries'][label]['seconds']=max(b['entries'][label]['seconds'],time.monotonic()-start)
    update_journal(Path(root)/'accounting.json',finish);return ok
