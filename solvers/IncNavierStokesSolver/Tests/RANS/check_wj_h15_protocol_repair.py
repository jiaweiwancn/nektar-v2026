#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Reproduce the publication race and test its repair without a CFD solve."""
import argparse
import os
from pathlib import Path
import shlex
import subprocess
import time
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,atomic

def run(root):
    out=root/'science/protocol-repair-checks-v2';out.mkdir()
    build=root/'build-protocol-v5';source=root/'source-protocol-v4'
    commands=load(root/'build/compile_commands.json')
    unit=next(v for v in commands if v['file'].endswith('/test_rans_wall_monitor.cpp'))
    argv=shlex.split(unit['command']);old_object=out/'old-scanner.o';old_exe=out/'old-scanner'
    argv[argv.index('-o')+1]=str(old_object);argv[argv.index(unit['file'])]=str(source/'solvers/IncNavierStokesSolver/Tests/RANS/test_rans_wall_monitor.cpp')
    assert '-ffp-contract=off'in argv and '-fno-fast-math'in argv
    subprocess.run(argv,cwd=unit['directory'],check=True)
    subprocess.run(['/usr/bin/c++',str(old_object),'-o',str(old_exe)],check=True)
    loader=root/'runtime/glibc/ld-linux-x86-64.so.2'
    libpath=str(root/'runtime/glibc')+':'+str(root/'runtime/geometry')+':/home/jay/.local/opt/nektar-v2026-mpi/lib64'
    records=[]
    for name,exe,expected in (('old-scanner',old_exe,1),('corrected-scanner',build/'RANSWallMonitorChecks',0)):
        log=out/(name+'.log');files=set();fixture=out/(name+'-records')
        command=[str(loader),'--library-path',libpath,str(exe),str(fixture)]
        with log.open('x')as f:
            p=subprocess.Popen(command,stdout=f,stderr=subprocess.STDOUT)
            while p.poll()is None:
                try:
                    files.update(str(Path(s.split(maxsplit=5)[5]).resolve())for s in Path(f'/proc/{p.pid}/maps').read_text().splitlines()if len(s.split(maxsplit=5))==6 and s.split(maxsplit=5)[5].startswith('/')and not s.endswith(' (deleted)'))
                except FileNotFoundError:pass
                time.sleep(.1)
        assert p.returncode==expected,(name,p.returncode)
        if expected:assert 'incorrect protocol decision: atomic-dual-name'in log.read_text()
        else:assert 'bounded catalog and protocol checks PASS'in log.read_text()
        libc=[v for v in files if Path(v).name=='libc.so.6'];assert libc==[str(root/'runtime/glibc/libc.so.6')]
        calls=(fixture/'catalog/calls.jsonl').read_text().splitlines()
        import json
        calls=[json.loads(v)for v in calls]
        assert [v['call']for v in calls]==list(range(1,1101))
        assert max(v['resident_definitions']for v in calls)<=16 and max(v['resident_key_bytes']for v in calls)<=32*1024**2
        records.append(dict(case=name,status='EXPECTED_OLD_FAILURE'if expected else'PASS',exit=p.returncode,command=command,ELF=identity(exe),log=identity(log),maps=[identity(p)for p in sorted(files)if Path(p).is_file()],calls=len(calls),resident_definitions=calls[-1]['resident_definitions'],evictions=calls[-1]['evictions']))
    cpp=load(build/'compile_commands.json');assert all('-ffp-contract=off'in v['command']and'-fno-fast-math'in v['command']for v in cpp)
    atomic(out/'result.json',dict(status='PASS',governing_attempts=0,CP0s=0,strict_views=0,records=records,original_header=identity(root/'source/solvers/IncNavierStokesSolver/EquationSystems/RANSWallMonitor.h'),corrected_header=identity(source/'solvers/IncNavierStokesSolver/EquationSystems/RANSWallMonitor.h'),
        repaired_build=dict(compile_commands=identity(build/'compile_commands.json'),compiled_inputs=[dict(source=identity(v['file']),command=v['command'])for v in cpp],binaries=[identity(build/name)for name in ('IncNavierStokesSolver','RANSWallResolutionProbe','RANSWallMonitorChecks')],governor_and_reader_execution='NOT_RUN after repair; original target reservation not repeated')))
    print('old atomic-publication regression FAIL as expected; corrected protocol/catalog PASS; zero CFD attempts',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);run(p.parse_args().root)
