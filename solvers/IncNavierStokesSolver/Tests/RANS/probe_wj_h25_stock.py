#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded stock startup diagnosis, not an admissibility/production benchmark."""
import argparse,os,subprocess,time
from pathlib import Path
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import atomic,identity

def probe(root):
    folder=root/'science/stock-startup-diagnosis-v1';folder.mkdir()
    runtime=Path('/home/jay/build/wj-h22-early-phase-replay-20261115/runtime')
    old=runtime/'IncNavierStokesSolver.elf';loader=runtime/'glibc/ld-linux-x86-64.so.2'
    libs=str(runtime/'glibc')+':'+str(runtime/'geometry')+':/home/jay/.local/opt/nektar-v2026-mpi/lib64'
    case=root/'science/S0-v2/case.xml'
    command=['gdb','--batch','-ex','set pagination off','-ex','set confirm off','-ex','run','-ex','bt 24','-ex','info sharedlibrary',
        '--args',str(loader),'--library-path',libs,str(old),str(case)]
    began=time.monotonic()
    with(folder/'gdb.log').open('x')as log:
        try:r=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,timeout=25);code=r.returncode
        except subprocess.TimeoutExpired:code='TIMEOUT_DIAGNOSTIC_ONLY'
    atomic(folder/'result.json',dict(role='startup allocator failure diagnosis; no WJ trajectory',command=command,exit=code,
        wall_seconds=time.monotonic()-began,ELF=identity(old),loader=identity(loader),input=identity(case),log=identity(folder/'gdb.log')))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();probe(a.root)
