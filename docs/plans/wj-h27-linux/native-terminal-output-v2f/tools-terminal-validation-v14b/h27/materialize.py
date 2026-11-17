#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Verify and materialize a fresh Mac-owned H27 source snapshot on Linux."""
import argparse,json,tarfile
from pathlib import Path,PurePosixPath
from resources import atomic,identity

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    p.add_argument('--archive',type=Path,required=True);p.add_argument('--manifest',type=Path,required=True)
    p.add_argument('--version',type=int,required=True);a=p.parse_args()
    plan=json.loads(a.manifest.read_text())
    assert identity(a.archive)['sha256']==plan['archive']['sha256']
    source=a.root/f'source-v{a.version}';assert not source.exists()
    with tarfile.open(a.archive)as archive:
        members=archive.getmembers()
        assert all(m.isfile()and not PurePosixPath(m.name).is_absolute()
                   and '..'not in PurePosixPath(m.name).parts for m in members)
        assert {m.name for m in members}=={r['path']for r in plan['files']}
        source.mkdir();archive.extractall(source)
    records=[]
    for row in plan['files']:
        measured=identity(source/row['path'])
        assert measured['sha256']==row['sha256']and measured['bytes']==row['bytes'],row['path']
        records.append(measured)
    atomic(a.root/f'source-v{a.version}.json',dict(
        scope='Mac-owned pre-checkpoint scientific source snapshot; no shared Linux checkout mutation',
        reviewed_base=plan['reviewed_base'],archive=identity(a.archive),
        manifest=identity(a.manifest),files=records))
    print('materialized',len(records),'files',source)

if __name__=='__main__':main()
