#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Finish H7 archival closure after storage repair; never run native programs.

Restores the two exact sealed transport archives, verifies decoded and retained
heavy bytes, and emits a new bridge. Closed H7 ledgers and source artifacts are
read-only. Every destination is new; historical paths inside packets stay intact.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tarfile
import time

ARCHIVES={
    'v1':('transport.tar.gz',4726782,'84eb376e256ee5b28536c49632ef64cbd4f3f378ad33b64904a2389536ee766f',68),
    'v2':('transport-owner-v2.tar.gz',4721706,'aa6672b6cfee0e20b7c81d86c9a12065109b898e65a741628382975f688144f0',69)}


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb')as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()

def identity(p):return dict(path=str(p),bytes=Path(p).stat().st_size,sha256=sha(p))
def load(p):return json.loads(Path(p).read_text())
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def verified(path,record):
    result=identity(path)
    if result['bytes']!=record['bytes']or result['sha256']!=record['sha256']:raise ValueError('sealed bytes differ: '+str(path))
    return result


def recover(root,folder,reuse_decoded=None):
    start=time.monotonic();root=root.resolve();folder=folder.resolve()
    if root.name!='wj-h7-neutral-transport-20261113'or folder.parent!=root:raise ValueError('H7 recovery subdirectory required')
    if folder.name not in ('storage-recovery-v1','storage-recovery-v2'):raise ValueError('fixed H7 recovery version required')
    if reuse_decoded is not None:
        reuse_decoded=reuse_decoded.resolve()
        if folder.name!='storage-recovery-v2'or reuse_decoded!=root/'storage-recovery-v1':raise ValueError('only v1 decoded bytes may supply v2')
    inputs=reuse_decoded or folder
    if(folder/'receipt.json').exists()or any(folder.glob('decoded-*')):raise ValueError('immutable recovery already executed')
    mount=subprocess.check_output(['findmnt','-no','SOURCE,TARGET,FSTYPE,UUID','-T',str(root)],text=True,timeout=5).strip()
    if mount.split()[-1]!='5AB8A669B8A642FD':raise ValueError('wrong storage UUID')
    intent=dict(source=identity(__file__),primary_root=str(root),wall_limit_seconds=180,extra_native_attempts=0,
        extra_cp0=0,new_output_ceiling_bytes=128*1024**2,rank='no MPI/native invocation',threads={k:os.environ.get(k)for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS')},
        original_study_allocation='closed1800s preserved; this user-requested metadata recovery is recorded separately')
    if any(v!='1'for v in intent['threads'].values()):raise ValueError('thread controls must all be1')
    save(folder/'intent.json',intent)
    def deadline():
        if time.monotonic()-start>180:raise TimeoutError('metadata recovery time limit')
    parent_before={p:identity(p)for p in (root/'accounting.json',root/'imports.json')};decoded={};seals={}
    try:
        for key,(name,n,h,count)in ARCHIVES.items():
            archive=inputs/name;verified(archive,dict(bytes=n,sha256=h));dest=inputs/('decoded-'+key)
            if reuse_decoded is None:dest.mkdir()
            with tarfile.open(archive,'r:gz')as t:
                members=t.getmembers()
                if sum(v.size for v in members)>64*1024**2:raise ValueError('decoded archive too large')
                for v in members:
                    path=Path(v.name)
                    if path.is_absolute()or '..'in path.parts or not(v.isfile()or v.isdir()):raise ValueError('unsafe archive member')
                if reuse_decoded is None:t.extractall(dest)
            seal=load(dest/'publication/seal.json')
            if len(seal['files'])!=count:raise ValueError('unexpected publication coverage')
            decoded[key]=[verified(dest/'publication'/v['name'],v)for v in seal['files']];seals[key]=identity(dest/'publication/seal.json');deadline()
        pub=inputs/'decoded-v2/publication';manifest=load(pub/'manifest.json')
        heavy=[verified(Path(v['path']),v)for v in manifest['heavy']];deadline()
        if len(heavy)!=105:raise ValueError('heavy coverage differs')
        h6=[]
        for row in load(pub/'imports.json'):
            for name in ('report','xml','checkpoint'):h6.append(verified(Path(row[name]['path']),row[name]))
        libraries=[verified(Path(v['path']),v)for v in load(pub/'runtime-equivalence.json')['qualified_ELFs']];deadline()
        maps=load(inputs/'decoded-v1/closure/loaded-native-closure.json')
        map_records=[verified(Path(v['maps']['path']),v['maps'])for v in maps]
        for p,v in parent_before.items():verified(p,v)
        original=inputs/'decoded-v1/closure/accounting-final.json';verified(original,load(inputs/'decoded-v1/closure/receipt.json')['accounting'])
        verified(root/'accounting.json',load(original)['primary_accounting_parent'])
        # Current H7 processes, not stale PIDs that may have been reused after reboot.
        live=[]
        for d in Path('/proc').iterdir():
            if not d.name.isdigit()or int(d.name)==os.getpid():continue
            try:args=(d/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace')
            except (FileNotFoundError,PermissionError,ProcessLookupError):continue
            if str(root/'observer-10a4078')in args and ('git 'in args or 'du -sb'in args):live.append(dict(pid=int(d.name),args=args))
        if live:raise ValueError('old H7 metadata operation still active')
        size=sum(p.stat().st_size for directory in ({folder,inputs})for p in directory.rglob('*')if p.is_file())
        if reuse_decoded is None:
            total=int(subprocess.check_output(['du','-sb',str(root)],text=True,timeout=10).split()[0]);size_scope='fresh whole-root measurement';storage_status='RESOLVED'
        else:
            # The failed v1 proved all sealed bytes before du hit EIO in the
            # obsolete, unused checkout. Preserve that error; do not repeat the
            # failing traversal or claim it is a fresh whole-root measurement.
            problem=load(inputs/'failure.json')
            prior=load(inputs/'decoded-v1/closure/receipt.json')['conservative_primary_plus_side_output_bound_bytes']
            total=prior+size;size_scope='sealed pre-repair bound plus measured recovery files; obsolete checkout EIO remains unmeasured';storage_status='ARCHIVAL_RECOVERED_AUXILIARY_EIO_UNRESOLVED'
        if size>128*1024**2 or total>8*1024**3:raise ValueError('output ceiling exceeded')
        receipt=dict(status='PASS',storage_status=storage_status,scientific_result='unchanged PASS,36/36 accepted atT=.04',mount=mount,source=intent['source'],
            elapsed_seconds=time.monotonic()-start,recovery_bytes=size,primary_root_bytes_or_recorded_bound=total,size_scope=size_scope,prior_failure=identity(inputs/'failure.json')if reuse_decoded else None,decoded_counts={k:len(v)for k,v in decoded.items()},publication_seals=seals,
            heavy_verified=heavy,H6_parent_verified=h6,current_qualified_ELFs_verified=libraries,original_native_map_receipts_verified=map_records,
            original_accounting=identity(original),primary_ledgers_unchanged=list(parent_before.values()),old_metadata_operations=live,
            path_aliases={str(Path('/tmp/wj-h7-compact-20261113/publication')):str(inputs/'decoded-v1/publication'),'/tmp/wj-h7-compact-20261113/closure':str(inputs/'decoded-v1/closure'),
                '/tmp/wj-h7-owner-repair/publication-owner-v2':str(pub)},native_attempts=0,cp0_initializations=0,
            preserved='old seals, paths, numerical decisions, failed operations and1800s charge unchanged; obsolete partial observer checkout retained read-only')
        # Flush only this recovery's files, never issue a volume-wide sync.
        for path in folder.rglob('*'):
            if path.is_file():
                with path.open('rb')as f:os.fsync(f.fileno())
        save(folder/'receipt.json',receipt)
        save(folder/'verification-index.json',dict(decoded=decoded,closure_files=[identity(p)for p in sorted((inputs/'decoded-v1/closure').rglob('*'))if p.is_file()]))
        print(json.dumps({k:v for k,v in receipt.items()if k not in ('heavy_verified','H6_parent_verified','current_qualified_ELFs_verified','original_native_map_receipts_verified')},indent=2),flush=True)
    except Exception as e:
        save(folder/'failure.json',dict(status='UNRESOLVED',reason=repr(e),elapsed_seconds=time.monotonic()-start,source=intent['source'],extra_native_attempts=0));raise

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--folder',type=Path,required=True);p.add_argument('--reuse-decoded',type=Path);a=p.parse_args();recover(a.root,a.folder,a.reuse_decoded)
