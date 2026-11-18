#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Build opt-in H27 native targets in an immutable source namespace."""
import argparse,json,os,shutil,subprocess
from pathlib import Path
from resources import THREADS,atomic,identity

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    p.add_argument('--version',type=int,required=True)
    p.add_argument('--targets',default='WJH27PrincipalProbe');a=p.parse_args()
    targets=a.targets.split(',');assert set(targets)<=set(('WJH27PrincipalProbe','WJH27MixedProbe','WJH27ScalarProbe','WJH27ReferenceProbe','WJH27AdmissibilityProbe','WJH27CoupledMMSProbe','WJH27ConstitutiveProbe','WJH27SteadyKrylovProbe','WJH27InitialFieldsProbe','WJH27PerturbationProbe','WJH27FlowTangentProbe','WJH27TraceGeometryProbe','WJH27CutoffPanelProbe','WJH27CutoffPanelsUnit','WJH27InnerKFixture','WJNativeFlowDriver','IncNavierStokesSolver'))
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
    command=['cmake','--build',str(build),'--target']+targets+['--parallel','4']
    subprocess.run(command,check=True)
    runtime.mkdir()
    old=Path('/home/jay/build/wj-h22-early-phase-replay-20261115/runtime')
    loader=old/'glibc/ld-linux-x86-64.so.2'
    libs=str(old/'glibc')+':'+str(old/'geometry')+':/home/jay/.local/opt/nektar-v2026-mpi/lib64'
    executables=[]
    for name in targets:
        binary=build/name
        if not binary.exists():binary=build/'h26-kernel'/name
        elf=runtime/(name+'.elf');shutil.copy2(binary,elf);wrapper=runtime/name
        wrapper.write_text('#!/bin/sh\nexec '+str(loader)+' --library-path '+libs+' '+str(elf)+' "$@"\n');wrapper.chmod(0o755)
        link=build/f'CMakeFiles/{name}.dir/link.txt'
        if not link.exists():link=build/f'h26-kernel/CMakeFiles/{name}.dir/link.txt'
        executables.append(dict(executable=identity(elf),wrapper=identity(wrapper),link=identity(link)))
    atomic(root/'science'/f'runtime-v{a.version}.json',dict(
        scope='opt-in H27 native targets; qualification scope belongs to each test',
        source=identity(root/f'source-v{a.version}.json'),
        configure=configure,build=command,executables=executables,
        loader=identity(loader),libraries=libs,
        compiler=subprocess.check_output(['c++','--version'],text=True),
        compile_commands=identity(build/'compile_commands.json')))

if __name__=='__main__':main()
