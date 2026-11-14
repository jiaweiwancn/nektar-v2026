#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""One write-once H16 bootstrap; small H15 essentials, no scientific replay."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[name]='1'
ROOT=Path('/home/jay/build/wj-h16-repaired-monitor-20261114')
OLD=Path('/home/jay/build/wj-h15-persistent-neutral-20261114')

def pin(p):
    p=Path(p);h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1024**2),b''):h.update(b)
    return dict(path=str(p),bytes=p.stat().st_size,sha256=h.hexdigest())

def save(p,v):
    tmp=Path(str(p)+'.tmp');tmp.write_text(json.dumps(v,indent=2)+'\n')
    with tmp.open('rb')as f:os.fsync(f.fileno())
    os.replace(tmp,p);fd=os.open(p.parent,os.O_RDONLY|os.O_DIRECTORY)
    try:os.fsync(fd)
    finally:os.close(fd)

def main():
    ROOT.mkdir();(ROOT/'operations').mkdir();(ROOT/'science').mkdir();(ROOT/'tmp').mkdir();(ROOT/'inputs').mkdir()
    folder=ROOT/'operations/essential-bootstrap-v1';folder.mkdir()
    book=dict(schema='WJH16LedgerV1',root=str(ROOT),ceiling_seconds=7200,setup_charged_seconds=60,closure_reserve_seconds=300,closed=False,
        starting_revision='b0afaa6cc9aa22ed85c2a0afa5eef8b8cc187f05',small_charged_attempts=0,target_reserved_attempts=0,entries={'essential-bootstrap-v1':dict(status='RESERVED',phase='operation',reserved_seconds=120)})
    save(ROOT/'accounting.json',book);start=time.monotonic()
    save(folder/'request.json',dict(source=pin(__file__),reserved_seconds=120,threads={n:os.environ[n]for n in os.environ if n.endswith('_NUM_THREADS')or n=='VECLIB_MAXIMUM_THREADS'},native_work=0))
    try:
        fs=subprocess.check_output(['findmnt','-n','-o','FSTYPE,SOURCE,TARGET','-T',str(ROOT)],text=True).strip();assert fs.startswith('ext4 ')
        assert shutil.disk_usage(ROOT).free>12*1024**3
        p=ROOT/'tmp/preflight.bin';data=bytes(range(256))*4096
        with p.open('xb')as f:f.write(data);f.flush();os.fsync(f.fileno())
        assert p.read_bytes()==data
        dest=ROOT/'inputs/h15-final';dest.mkdir();names=['publication-v3/manifest.json','publication-v3/seal.json','closure-v3/seal.json','closure-v3/accounting-final.json','closure-v3/verification.json','closure-v3/post-stop-build-binaries.json']
        original={};copied=[]
        for name in names:
            assert time.monotonic()-start<120
            q=OLD/name;rec=pin(q);original[name]=rec;to=dest/name;to.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(q,to);assert pin(to)['sha256']==rec['sha256'];copied.append(dict(original=rec,copy=pin(to)))
        manifest=json.loads((OLD/names[0]).read_text());assert original[names[0]]['bytes']==52591 and original[names[0]]['sha256']=='f7e60a290093d7dfba631bc68a578a8afdaf595a61febf7c4ba18aba3525913b'
        oldbook=json.loads((OLD/'closure-v3/accounting-final.json').read_text());assert oldbook['closed']and oldbook['total_charged_seconds']==2265.5647236939985
        binaries=json.loads((OLD/'closure-v3/post-stop-build-binaries.json').read_text())
        selected=[]
        for build in binaries['builds']:
            for rec in build['binaries']+build['link_commands']:assert pin(rec['path'])==rec;selected.append(rec)
        cache=pin(OLD/'science/observer-cache/immutable-cache.pickle');assert cache['bytes']==211582940 and cache['sha256']=='800ac28c73ecb6498e6504b4e14e774db2cff893f2c0aaf418cb7036d50d00b0'
        small=OLD/'science/small/legacy/case.evolution.json';heavy=json.loads((OLD/'publication-v3/remote-heavy-manifest.json').read_text());legacy=next(v for v in heavy['files']if v['path']==str(small));assert pin(small)==legacy
        save(ROOT/'inputs/essential-index.json',dict(status='PASS',filesystem=fs,free_bytes=shutil.disk_usage(ROOT).free,preflight=pin(p),historical_files=copied,selected_repair_build_files=selected,legacy_producer=legacy,cache=cache,
            historical_cache_manifest=pin(OLD/'science/observer-cache/manifest.json'),small_cache_manifest=pin(OLD/'science/small-cache-v2/manifest.json'),runtime_manifest=pin(OLD/'science/runtime-build-v1.json'),runtime_qualification=pin(OLD/'science/runtime-qualified-v1.json'),compiled_sources=pin(OLD/'publication-v3/compiled-source-objects.json'),
            old_scanner_proof=pin(OLD/'science/protocol-repair-checks-v2/result.json'),old_ledger_unchanged=original['closure-v3/accounting-final.json'],full_archive_not_required=True))
        entry=dict(status='PASS',phase='operation',reserved_seconds=120,seconds=time.monotonic()-start,source=pin(__file__),native_work=0)
    except Exception as error:
        entry=dict(status='OBSTRUCTED',phase='operation',reserved_seconds=120,seconds=time.monotonic()-start,reason=repr(error),native_work=0)
        save(folder/'receipt.json',entry);book['entries']['essential-bootstrap-v1']=entry;save(ROOT/'accounting.json',book);raise
    book['entries']['essential-bootstrap-v1']=entry;save(ROOT/'accounting.json',book);save(folder/'receipt.json',entry)
    print(json.dumps(dict(status='PASS',root=str(ROOT),seconds=entry['seconds'],index=pin(ROOT/'inputs/essential-index.json'))),flush=True)

if __name__=='__main__':main()
