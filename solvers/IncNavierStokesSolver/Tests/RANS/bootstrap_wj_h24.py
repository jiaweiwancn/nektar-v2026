#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""One fresh ext4 H24 ledger; no scientific/native initialization."""
import argparse,hashlib,json,os,shutil,subprocess,time
from pathlib import Path

def initialize(root,revision):
    start=time.monotonic();assert root.name=='wj-h24-flat-continuation-20261115'
    root.mkdir();[(root/n).mkdir()for n in ('inputs','science','operations','tmp')]
    def put(path,value):
        with path.open('x')as f:json.dump(value,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
    book=dict(schema='WJH24LedgerV1',closed=False,ceiling_seconds=6000,native_ceiling=4200,passive_ceiling=900,
        setup_charged_seconds=60,closure_reserve_seconds=300,native_reservations=[],reserved_attempts=0,entries={
            'source-transfer-and-initial-inspection-v1':dict(status='FULLY_CHARGED',phase='passive',kind='passive',seconds=30,reserved_seconds=30),
            'coordination-and-transport-v1':dict(status='FULLY_CHARGED',phase='passive',kind='passive',seconds=120,reserved_seconds=120)})
    put(root/'accounting.json',book)
    assert shutil.disk_usage(root).free>10*1024**3
    data=bytes(range(256))*4096;p=root/'tmp/fsync-readback'
    with p.open('xb')as f:f.write(data);f.flush();os.fsync(f.fileno())
    assert p.read_bytes()==data
    for d in (p.parent,root,root.parent):
        fd=os.open(d,os.O_RDONLY|os.O_DIRECTORY)
        try:os.fsync(fd)
        finally:os.close(fd)
    workers=subprocess.check_output(['ps','-eo','pid,ppid,nlwp,rss,etimes,args'],text=True)
    put(root/'preflight-host.json',dict(starting_revision=revision,filesystem=subprocess.check_output(['findmnt','-n','-o','FSTYPE','-T',str(root)],text=True).strip(),
        free_bytes=shutil.disk_usage(root).free,write_flush_readback_sha256=hashlib.sha256(data).hexdigest(),
        process_listing=workers,seconds=time.monotonic()-start,scope='initial read-only inspection and bounded preflight inside fully charged setup/transfer allowances'))
    print(root,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--revision',required=True);a=p.parse_args();initialize(a.root,a.revision)
