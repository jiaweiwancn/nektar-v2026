#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""New H26 native executable, isolated from shared installation/build state."""
import argparse,json,os,shutil,subprocess,time
from pathlib import Path
from run import THREADS
os.environ.update(THREADS)
from wj_h15_observer import identity,atomic


def build(root,version):
    source=root/f'source-v{version}';output=root/f'build-v{version}';runtime=root/f'runtime-v{version}'
    assert source.exists()and not output.exists()and not runtime.exists()
    cmake=['cmake','-S',str(source/'solvers/IncNavierStokesSolver/Tests/RANS/h26'),'-B',str(output),'-DCMAKE_BUILD_TYPE=Release','-DCMAKE_CXX_FLAGS=-ffp-contract=off -fno-fast-math','-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake','-DCMAKE_EXPORT_COMPILE_COMMANDS=ON']
    subprocess.run(cmake,check=True);command=['cmake','--build',str(output),'--parallel','4'];subprocess.run(command,check=True)
    runtime.mkdir();elf=runtime/'WJNativeFlowDriver.elf';shutil.copy2(output/'WJNativeFlowDriver',elf)
    old=Path('/home/jay/build/wj-h22-early-phase-replay-20261115/runtime');wrapper=runtime/'WJNativeFlowDriver'
    loader=old/'glibc/ld-linux-x86-64.so.2';libraries=str(old/'glibc')+':'+str(old/'geometry')+':/home/jay/.local/opt/nektar-v2026-mpi/lib64'
    wrapper.write_text('#!/bin/sh\nexec '+str(loader)+' --library-path '+libraries+' '+str(elf)+' "$@"\n');wrapper.chmod(0o755)
    atomic(root/f'science/runtime-v{version}.json',dict(status='PASS',source_manifest=identity(root/f'source-v{version}.json'),configure=cmake,build=command,compiler=subprocess.check_output(['c++','--version'],text=True),executable=identity(elf),wrapper=identity(wrapper),loader=identity(loader),libraries=libraries,compile_commands=identity(output/'compile_commands.json'),link=identity(output/'CMakeFiles/WJNativeFlowDriver.dir/link.txt')))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--version',type=int,required=True);a=p.parse_args();build(a.root,a.version)
