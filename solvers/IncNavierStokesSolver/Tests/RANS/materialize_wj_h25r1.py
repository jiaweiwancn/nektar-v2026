#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Write-once source overlay transport, preserving every historical capsule."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import tarfile
import time


def identity(p):
    p=Path(p)
    return dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())


def materialize(parent,revision,expected,version=1):
    assert version in (1,2,3)
    root=parent/'resumption-r1';archive=root/('overlay-v'+str(version)+'.tar.gz')
    assert identity(archive)['sha256']==expected
    source=root/('source-v'+str(version));assert not source.exists()
    source_base=parent/'source-v8' if version==1 else root/('source-v'+str(version-1))
    base_manifest=parent/'source-manifest-v8.json' if version==1 else root/('source-manifest-v'+str(version-1)+'.json')
    started=time.monotonic();shutil.copytree(source_base,source,symlinks=True)
    with tarfile.open(archive) as a:
        members=a.getmembers()
        for m in members:
            assert not Path(m.name).is_absolute() and '..' not in Path(m.name).parts and (m.isfile() or m.isdir())
            p=source/m.name
            if p.exists() and p.is_file(): p.chmod(p.stat().st_mode | 0o200)
        a.extractall(source)
    files=[identity(p) for p in sorted(source.rglob('*')) if p.is_file()]
    manifest=dict(schema='WJH25R1SourceOverlayV1',source_revision=revision,
        overlay=identity(archive),base=identity(base_manifest),
        source=str(source),files=files,seconds=time.monotonic()-started,
        distinction='base native capsule plus listed committed overlay; documentation Git HEAD is not executable provenance')
    manifest_path=root/('source-manifest-v'+str(version)+'.json')
    with manifest_path.open('x') as f:json.dump(manifest,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
    print(json.dumps(dict(source=str(source),files=len(files),seconds=manifest['seconds'],manifest=identity(manifest_path))))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--revision',required=True);p.add_argument('--sha256',required=True)
    p.add_argument('--version',type=int,default=1)
    a=p.parse_args();materialize(a.root,a.revision,a.sha256,a.version)
