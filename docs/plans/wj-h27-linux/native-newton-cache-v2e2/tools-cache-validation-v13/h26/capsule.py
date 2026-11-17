#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Materialize a new H26 source capsule from a pinned base and Git archive."""
import argparse,hashlib,json,shutil,tarfile
from pathlib import Path

def identity(p):return dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def materialize(root,version,revision):
    source=root/f'source-v{version}';overlay=root/f'overlay-v{version}.tar';assert not source.exists()
    shutil.copytree(root/'source-v1',source,ignore=shutil.ignore_patterns('__pycache__'))
    with tarfile.open(overlay)as archive:
        for member in archive.getmembers():
            assert not member.name.startswith('/')and'..'not in Path(member.name).parts and not member.issym()and not member.islnk()
        archive.extractall(source)
    (root/f'source-v{version}.json').write_text(json.dumps(dict(revision=revision,base=identity(root/'source-v1.json'),overlay=identity(overlay),files=[identity(p)for p in sorted(source.rglob('*'))if p.is_file()]),indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--version',type=int,required=True);p.add_argument('--revision',required=True);a=p.parse_args();materialize(a.root,a.version,a.revision)
