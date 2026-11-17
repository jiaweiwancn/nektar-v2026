#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Materialize current native inputs over a read-only qualified H26 capsule."""
import argparse,json,shutil,tarfile
from pathlib import Path,PurePosixPath
from resources import atomic,identity

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    p.add_argument('--archive',type=Path,required=True);p.add_argument('--sha256',required=True)
    p.add_argument('--version',type=int,required=True);a=p.parse_args()
    assert identity(a.archive)['sha256']==a.sha256
    source=a.root/f'source-v{a.version}';assert not source.exists()
    base=Path('/home/jay/build/wj-h26-native-flow-20261117/source-v8')
    with tarfile.open(a.archive)as archive:
        members=archive.getmembers()
        assert all(m.isfile()and not PurePosixPath(m.name).is_absolute()
                   and '..'not in PurePosixPath(m.name).parts for m in members)
        plan=json.load(archive.extractfile('mac-manifest.json'))
        expected={v['path']:v for v in plan['files']}
        checked=[]
        # Every source/header available to the native target is byte-checked
        # against the Mac checkout; historical scripts are not claimed current.
        for f in base.rglob('*'):
            rel=str(f.relative_to(base))
            native=rel.startswith(('library/','solvers/IncNavierStokesSolver/'))
            if f.is_file()and native and(f.suffix in('.h','.hpp','.cpp','.c','.inl')or f.name=='CMakeLists.txt'):
                record=identity(f);assert rel in expected,rel
                assert record['sha256']==expected[rel]['sha256'],rel
                checked.append(dict(relative_path=rel,**record))
        shutil.copytree(base,source,ignore=shutil.ignore_patterns('__pycache__'))
        archive.extractall(source,members=[m for m in members if m.name!='mac-manifest.json'])
    for row in plan['overlay_files']:
        got=identity(source/row['path']);assert got['sha256']==row['sha256'],row['path']
    atomic(a.root/f'source-v{a.version}.json',dict(
        scope='Verified current Mac native inputs and H27/inlet overlay over immutable historical capsule',
        reviewed_base=plan['reviewed_base'],historical_source_manifest=identity(base.parent/'source-v8.json'),
        archive=identity(a.archive),verified_native_inputs=checked,
        current_overlay=plan['overlay_files'],
        files=[identity(f)for f in sorted(source.rglob('*'))if f.is_file()]))
    print('materialized',source,'verified native inputs',len(checked))

if __name__=='__main__':main()
