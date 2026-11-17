#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Build the H27 constitutive probe in its own immutable source namespace."""
import argparse,json,os,shutil,subprocess
from pathlib import Path
from resources import THREADS,atomic,identity

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    p.add_argument('--version',type=int,required=True);a=p.parse_args()
    root=a.root.resolve();source=root/f'source-v{a.version}'
    build=root/f'build-v{a.version}';runtime=root/f'runtime-v{a.version}'
    assert source.is_dir()and not build.exists()and not runtime.exists()
    os.environ.update(THREADS)
    configure=['cmake','-S',str(source/'solvers/IncNavierStokesSolver/Tests/RANS/h27'),
        '-B',str(build),'-DCMAKE_BUILD_TYPE=Release',
        '-DCMAKE_CXX_FLAGS=-ffp-contract=off -fno-fast-math',
        '-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake',
        '-DCMAKE_EXPORT_COMPILE_COMMANDS=ON']
    subprocess.run(configure,check=True)
    command=['cmake','--build',str(build),'--target','WJH27PrincipalProbe','--parallel','4']
    subprocess.run(command,check=True)
    runtime.mkdir();elf=runtime/'WJH27PrincipalProbe.elf'
    shutil.copy2(build/'WJH27PrincipalProbe',elf)
    old=Path('/home/jay/build/wj-h22-early-phase-replay-20261115/runtime')
    loader=old/'glibc/ld-linux-x86-64.so.2'
    libs=str(old/'glibc')+':'+str(old/'geometry')+':/home/jay/.local/opt/nektar-v2026-mpi/lib64'
    wrapper=runtime/'WJH27PrincipalProbe'
    wrapper.write_text('#!/bin/sh\nexec '+str(loader)+' --library-path '+libs+' '+str(elf)+' "$@"\n')
    wrapper.chmod(0o755)
    atomic(root/'science'/f'runtime-v{a.version}.json',dict(
        scope='constitutive feasibility, no coupled PDE',
        source=identity(root/f'source-v{a.version}.json'),
        configure=configure,build=command,executable=identity(elf),
        wrapper=identity(wrapper),loader=identity(loader),libraries=libs,
        compiler=subprocess.check_output(['c++','--version'],text=True),
        compile_commands=identity(build/'compile_commands.json'),
        link=identity(build/'CMakeFiles/WJH27PrincipalProbe.dir/link.txt')))

if __name__=='__main__':main()
