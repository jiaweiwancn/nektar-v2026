#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""New H26 native executable, isolated from shared installation/build state."""
import argparse,json,os,shutil,subprocess,time,shlex
from pathlib import Path
from run import THREADS
os.environ.update(THREADS)
from wj_h15_observer import identity,atomic


def build(root,version):
    source=root/f'source-v{version}';output=root/f'build-v{version}';runtime=root/f'runtime-v{version}'
    assert source.exists()and not output.exists()and not runtime.exists()
    cmake=['cmake','-S',str(source/'solvers/IncNavierStokesSolver/Tests/RANS/h26'),'-B',str(output),'-DCMAKE_BUILD_TYPE=Release','-DCMAKE_CXX_FLAGS=-ffp-contract=off -fno-fast-math','-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake','-DCMAKE_EXPORT_COMPILE_COMMANDS=ON']
    reuse=[]
    if(source/'solvers/IncNavierStokesSolver/Tests/RANS/h26/stock.cpp').exists():
        previous=Path('/home/jay/build/wj-h25-production-20261116/resumption-r1/build-memo-v1');oldsource=previous.parent/'source-v3';cmake+=['-DH26_STOCK_BUILD='+str(previous)]
        commands=json.loads((previous/'compile_commands.json').read_text())
        for obj in sorted((previous/'CMakeFiles/IncNavierStokesSolver.dir').rglob('*.o')):
            if obj.name.startswith('RANS')or obj.name=='IncNavierStokesSolver.cpp.o':continue
            dep=Path(str(obj)+'.d');assert dep.exists(),str(dep)
            pairs=[]
            for raw in shlex.split(dep.read_text().replace('\\\n',' ').split(':',1)[1]):
                p=Path(raw);new=source/p.relative_to(oldsource)if p.is_relative_to(oldsource)else p
                if not new.exists()and '/Tests/RANS/include/'in str(p):new=p
                assert p.exists()and new.exists()and identity(p)['sha256']==identity(new)['sha256'],(p,new)
                pairs.append(identity(new))
            matching=[c for c in commands if str(obj.relative_to(previous))in c['command']];assert len(matching)==1
            reuse.append(dict(object=identity(obj),dependency=identity(dep),command=matching[0],inputs=pairs))
    subprocess.run(cmake,check=True);command=['cmake','--build',str(output),'--parallel','4'];subprocess.run(command,check=True)
    runtime.mkdir();old=Path('/home/jay/build/wj-h22-early-phase-replay-20261115/runtime')
    loader=old/'glibc/ld-linux-x86-64.so.2';libraries=str(old/'glibc')+':'+str(old/'geometry')+':/home/jay/.local/opt/nektar-v2026-mpi/lib64'
    executables=[]
    for name in ('WJNativeFlowDriver','WJStockFlowDriver'):
        if not(output/name).exists():continue
        elf=runtime/(name+'.elf');shutil.copy2(output/name,elf);wrapper=runtime/name
        wrapper.write_text('#!/bin/sh\nexec '+str(loader)+' --library-path '+libraries+' '+str(elf)+' "$@"\n');wrapper.chmod(0o755);executables.append(dict(executable=identity(elf),wrapper=identity(wrapper),link=identity(output/f'CMakeFiles/{name}.dir/link.txt')))
    atomic(root/f'science/runtime-v{version}.json',dict(status='PASS',source_manifest=identity(root/f'source-v{version}.json'),configure=cmake,build=command,compiler=subprocess.check_output(['c++','--version'],text=True),executables=executables,qualified_stock_objects=reuse,loader=identity(loader),libraries=libraries,compile_commands=identity(output/'compile_commands.json')))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--version',type=int,required=True);a=p.parse_args();build(a.root,a.version)
