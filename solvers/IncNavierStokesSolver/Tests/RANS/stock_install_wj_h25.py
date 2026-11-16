#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check the installed stock executable independently of the RANS capsule ABI."""
import argparse,os,re,subprocess,time
from pathlib import Path
import xml.etree.ElementTree as ET
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,atomic
from cls_epawf_m3_cases import write

def run(root):
    d=root/'science/S0-installed';d.mkdir();x=ET.parse(root/'science/S0-NoCollection/case.xml').getroot();write(x,d/'case.xml')
    elf=Path('/home/jay/.local/opt/nektar-v2026-mpi/bin/IncNavierStokesSolver')
    rt=Path('/home/jay/build/wj-h22-early-phase-replay-20261115/runtime')
    cmd=[str(rt/'glibc/ld-linux-x86-64.so.2'),'--library-path',str(rt/'glibc')+':'+str(rt/'geometry')+':/home/jay/.local/opt/nektar-v2026-mpi/lib64',str(elf),str(d/'case.xml')]
    start=time.monotonic()
    with(d/'native.log').open('x')as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
    text=(d/'native.log').read_text();samples=[dict(step=int(s),time=float(t),seconds=float(w))for s,t,w in re.findall(r'Steps:\s*(\d+)\s+Time:\s*([\d.eE+-]+)\s+CPU Time:\s*([\d.eE+-]+)s',text)]
    atomic(d/'result.json',dict(status='PASS'if r.returncode==0 and len(samples)==8 else'FAIL',exit=r.returncode,samples=samples,
        seconds=time.monotonic()-start,ELF=identity(elf),command=cmd,SVV=False,scope='installed stock capability/ABI control; compiled source qualification remains separately required'))
    assert r.returncode==0 and len(samples)==8

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();run(a.root)
