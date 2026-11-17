#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""New H27-only finite journal and bounded process-tree execution."""
import argparse,datetime,json,os,signal,subprocess,time
from pathlib import Path
from resources import THREADS,identity,atomic,update,workers,tree,size

def run(root,label,argv,phase,seconds,evidence=False):
    root=Path(root).resolve();assert root.name.startswith('wj-h27-')and Path('/home/jay/build')in root.parents
    assert phase in ('build','science','analysis','evidence')and 0<seconds<=dict(build=1800,science=7200,analysis=900,evidence=900)[phase]
    book=root/'accounting.json';b=json.loads(book.read_text());assert not b['closed']and label not in b['entries']
    allocation=workers();reservation=16 if phase=='build'else 10 if phase=='science'else 2;operation=root/'operations'/label;operation.mkdir()
    inputs=[]
    for argument in argv:
        try:
            if '\n'not in argument and Path(argument).is_file():inputs.append(identity(argument))
        except OSError:pass
    request=dict(argv=argv,phase=phase,seconds_cap=seconds,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),allocation=allocation,reservation=reservation,threads=THREADS,runner=identity(__file__),input_files=inputs)
    atomic(operation/'request.json',request)
    admitted=allocation['workers']+reservation<=72 and os.statvfs(root).f_bavail*os.statvfs(root).f_frsize>=15*1024**3
    if not admitted:
        event=dict(status='RESOURCE_NOT_ADMITTED',phase=phase,seconds=1.)
        update(book,lambda x:x['entries'].update({label:event}));atomic(operation/'receipt.json',dict(request,**event));return 1
    def reserve(x):
        assert not x['closed']and label not in x['entries'];used=sum(v.get('seconds',v.get('reserved_seconds',0))for v in x['entries'].values())
        assert used+seconds+(0 if evidence or phase=='evidence'else 3600)<=43200,'H27 allowance/evidence reserve'
        x['entries'][label]=dict(status='RESERVED',phase=phase,reserved_seconds=seconds)
    update(book,reserve);began=time.monotonic();peak=0;peak_workers=allocation['workers'];reason=None;groups=set();env=dict(os.environ,**THREADS,PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(root/'tmp'),H27_ROOT=str(root),H27_OPERATION=str(operation),PYTHONPATH=str(Path(__file__).resolve().parent))
    with(operation/'command.log').open('x')as log:
        child=subprocess.Popen(argv,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        while child.poll()is None:
            rss,newgroups=tree(child.pid);peak=max(peak,rss);groups|=newgroups;now=time.monotonic();current=workers()['workers'];peak_workers=max(peak_workers,current)
            if rss>(16 if phase=='build'else 8)*1024**3:reason='RSS_LIMIT'
            if now-began>seconds:reason='TIME_LIMIT'
            if current>72:reason='HOST_WORKER_LIMIT'
            if os.statvfs(root).f_bavail*os.statvfs(root).f_frsize<15*1024**3:reason='FREE_LIMIT'
            if size(root/'science')>8*1024**3:reason='SCIENCE_STORAGE_LIMIT'
            if sum(size(p)for p in root.iterdir()if p.is_dir()and p.name.startswith(('source','build','runtime')))>8*1024**3:reason='SOURCE_BUILD_STORAGE_LIMIT'
            if(operation/'command.log').stat().st_size>32*1024**2:reason='LOG_LIMIT'
            if reason:
                for group in groups|{child.pid}:
                    try:os.killpg(group,signal.SIGTERM)
                    except ProcessLookupError:pass
                break
            time.sleep(.5)
        try:code=child.wait(timeout=5)
        except subprocess.TimeoutExpired:
            for group in groups|{child.pid}:
                try:os.killpg(group,signal.SIGKILL)
                except ProcessLookupError:pass
            code=child.wait()
    event=dict(status='PASS'if code==0 and reason is None else'FAIL',phase=phase,seconds=time.monotonic()-began,exit=code,reason=reason,peak_tree_rss_bytes=peak,peak_host_workers=peak_workers,log=identity(operation/'command.log'))
    update(book,lambda x:x['entries'].update({label:event}));atomic(operation/'receipt.json',dict(request,**event));print(label,event['status'],event['seconds'],flush=True);return int(event['status']!='PASS')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--label',required=True);p.add_argument('--phase',required=True);p.add_argument('--seconds',type=float,required=True);p.add_argument('--evidence',action='store_true');p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args();raise SystemExit(run(a.root,a.label,a.command[1:]if a.command[:1]==['--']else a.command,a.phase,a.seconds,a.evidence))
