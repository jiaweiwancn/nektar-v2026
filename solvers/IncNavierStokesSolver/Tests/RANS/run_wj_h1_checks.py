#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Small H1 CTest receipts/resource journal; no old package ledger is opened."""
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
from run_cls_epawf_picard import THREADS

def load(p):return json.loads(Path(p).read_text())
def save(p,v):
    tmp=Path(str(p)+'.tmp');tmp.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n');os.replace(tmp,p)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def size(p):return int(subprocess.check_output(['du','-sb',str(p)],text=True).split()[0])
def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--build',type=Path,required=True)
    p.add_argument('--label',required=True);p.add_argument('--select',required=True);p.add_argument('--units',type=int,default=0)
    p.add_argument('--seconds',type=int,default=1200);a=p.parse_args()
    root=a.root.resolve();assert root.name.startswith('wj-h1-') and a.units>=0 and 0<a.seconds<=3600
    out=root/'operations'/a.label;out.mkdir(parents=True,exist_ok=False)
    journal=root/'checks-accounting.json'
    with (root/'checks-accounting.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if journal.exists():book=load(journal)
        else:
            old=root/'operations/geometry-v1.json'
            book=dict(schema='WJH1BoundedChecksV1',unit_ceiling=40,seconds_ceiling=10800,
                closure_reserve_seconds=600,setup_allowance_seconds=120,
                prior_geometry_seconds=load(old)['seconds'] if old.exists() else 0.,entries={})
        assert a.label not in book['entries'],'immutable H1 check label already reserved'
        units=sum(v['units'] for v in book['entries'].values())
        seconds=sum(v.get('seconds',v['reserved_seconds']) for v in book['entries'].values())
        seconds+=book['closure_reserve_seconds']+book['setup_allowance_seconds']+book['prior_geometry_seconds']
        assert units+a.units<=40 and seconds+a.seconds<=10800,'H1 RESOURCE_LIMIT before launch'
        allocation=scientific_workers();assert allocation['workers']+8<=72,'current worker allocation unavailable'
        entry=dict(units=a.units,reserved_seconds=a.seconds,status='RESERVED',allocation=allocation)
        book['entries'][a.label]=entry;save(journal,book)
        argv=['ctest','--test-dir',str(a.build.resolve()),'-j','1','--output-on-failure','-R',a.select]
        env=dict(os.environ,**THREADS,NUMEXPR_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',
            NEKTAR_RANS_TIMING_FILE=str(out/'native-timing.jsonl'))
        inventory=subprocess.run(['ctest','--test-dir',str(a.build.resolve()),'--show-only=json-v1'],env=env,
            capture_output=True,text=True,check=True)
        (out/'inventory.json').write_text(inventory.stdout)
        import re
        selected=[v['name'] for v in json.loads(inventory.stdout)['tests'] if re.search(a.select,v['name'])]
        assert selected,'zero selected H1 tests is not validation'
        before=size(root);groups=set();peak=threads=0;reason=None;start=time.monotonic();lastSize=start
        request=dict(argv=argv,selected=selected,registered=len(json.loads(inventory.stdout)['tests']),
            timing_file=env['NEKTAR_RANS_TIMING_FILE'],driver_sha256=sha(__file__),allocation=allocation)
        save(out/'request.json',request)
        with (out/'command.log').open('x') as log:
            child=subprocess.Popen(['/usr/bin/time','-v',*argv],cwd=root,env=env,
                stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            while child.poll() is None:
                rss,n=tree_usage(child.pid);peak=max(peak,rss);threads=max(threads,n);groups|=owned_groups(child.pid)
                if rss>8*1024**3:reason='RESOURCE_LIMIT: whole-process RSS'
                if time.monotonic()-start>a.seconds:reason='RESOURCE_LIMIT: check timeout'
                if scientific_workers()['workers']>72:reason='RESOURCE_LIMIT: global workers'
                if time.monotonic()-lastSize>10:
                    current=size(root);lastSize=time.monotonic()
                    if current>16*1024**3 or current-before>4*1024**3:reason='RESOURCE_LIMIT: output'
                if reason:
                    for group in groups:
                        try:os.killpg(group,signal.SIGTERM)
                        except ProcessLookupError:pass
                    break
                time.sleep(.5)
            try:code=child.wait(timeout=10)
            except subprocess.TimeoutExpired:
                for group in groups:
                    try:os.killpg(group,signal.SIGKILL)
                    except ProcessLookupError:pass
                code=child.wait()
        elapsed=time.monotonic()-start
        entry.update(status='PASS' if code==0 and not reason else 'RESOURCE_LIMIT' if reason else 'FAIL',
            seconds=elapsed,exit=code,reason=reason,peak_tree_rss_bytes=peak,peak_tree_threads=threads,
            output_bytes=size(root)-before,root_bytes=size(root),log_sha256=sha(out/'command.log'))
        save(out/'receipt.json',dict(request,**entry));save(journal,book)
        print(json.dumps(dict(label=a.label,**entry)))
        raise SystemExit(0 if code==0 and not reason else 1)
if __name__=='__main__':main()
