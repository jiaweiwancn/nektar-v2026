#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Fresh finite H25 account and bounded storage preflight, no numerical imports."""
import argparse,hashlib,json,os,shutil,subprocess,time
from pathlib import Path

def initialize(root,revision):
    assert root.name=='wj-h25-production-20261116'
    root.mkdir()
    for name in ('inputs','science','operations','tmp'): (root/name).mkdir()
    def put(path,value):
        with path.open('x') as f:
            json.dump(value,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
    put(root/'accounting.json',dict(schema='WJH25FiniteDevelopmentV1',closed=False,
        ceilings=dict(total=14400,build=2400,native=7200,passive=3600),
        setup_seconds=120,closure_seconds=600,entries={},reserved_attempts=0,
        maximum_attempts=200,governor_blocks=0,maximum_governor_blocks=12,
        source_revision=revision,scope='new H25 only; all earlier accounts stay CLOSED'))
    assert shutil.disk_usage(root).free>=15*1024**3
    raw=bytes(range(256))*4096;q=root/'tmp/fsync-readback'
    with q.open('xb') as f:f.write(raw);f.flush();os.fsync(f.fileno())
    assert q.read_bytes()==raw
    for d in (q.parent,root,root.parent):
        fd=os.open(d,os.O_RDONLY|os.O_DIRECTORY)
        try:os.fsync(fd)
        finally:os.close(fd)
    put(root/'preflight.json',dict(revision=revision,free_bytes=shutil.disk_usage(root).free,
        filesystem=subprocess.check_output(['findmnt','-n','-o','FSTYPE','-T',str(root)],text=True).strip(),
        cpu=json.loads(subprocess.check_output(['lscpu','-J'],text=True)),
        processes=subprocess.check_output(['ps','-eo','pid,ppid,nlwp,rss,etimes,args'],text=True),
        readback_sha256=hashlib.sha256(raw).hexdigest(),setup_fully_charged=120))
    print(root)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--revision',required=True)
    a=p.parse_args();initialize(a.root,a.revision)
