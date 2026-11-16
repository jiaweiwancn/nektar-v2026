#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Supported explicit opt-file avoids SessionReader's missing-opt autotuning."""
import argparse,os,re,subprocess,time
from pathlib import Path
import xml.etree.ElementTree as ET
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,atomic
from cls_epawf_m3_cases import write
from study_wj_h25 import executable

def run(root):
    d=root/'science/S0-explicit-opt';d.mkdir();case=d/'case.xml'
    write(ET.parse(root/'science/S0-NoCollection/case.xml').getroot(),case)
    opt=d/'case.opt'
    opt.write_text('<?xml version="1.0"?><NEKTAR><COLLECTIONS DEFAULT="NoCollection"/></NEKTAR>\n')
    command=[str(executable(root,'IncNavierStokesSolver')),'--use-opt-file',str(opt),str(case)]
    atomic(d/'execution-intent.json',dict(command=command,case=identity(case),opt=identity(opt),SVV=False,
        source_reason='SessionReader.cpp missing-opt path sets m_updateOptFile=true; explicit --use-opt-file preserves false; ExpList.cpp otherwise invokes SetWithTimings',
        scientific_changes=False,previous_failures=['stock-S0-700-v1','stock-S0-NoCollection-v2','stock-installed-capability-v1']))
    started=time.monotonic();maps=None
    with(d/'native.log').open('x')as log:
        p=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT)
        while p.poll()is None:
            try:
                files=sorted({s.split(maxsplit=5)[5]for s in Path(f'/proc/{p.pid}/maps').read_text().splitlines()if len(s.split(maxsplit=5))==6 and s.split(maxsplit=5)[5].startswith('/')and not s.endswith(' (deleted)')})
                if any('IncNavierStokesSolver.elf'in f for f in files)and any(f.endswith('/libc.so.6')for f in files):maps=files
            except FileNotFoundError:pass
            time.sleep(.025)
        code=p.returncode
    text=(d/'native.log').read_text();samples=[dict(step=int(s),time=float(t),seconds=float(w))for s,t,w in re.findall(r'Steps:\s*(\d+)\s+Time:\s*([\d.eE+-]+)\s+CPU Time:\s*([\d.eE+-]+)s',text)]
    atomic(d/'result.json',dict(status='PASS'if code==0 and len(samples)==8 else'FAIL',exit=code,seconds=time.monotonic()-started,samples=samples,
        maps=maps,source=load(root/'science/runtime-build-v1.json'),input=identity(case),opt=identity(opt),log=identity(d/'native.log'),SVV=False,
        interval=[.04,.06],dt=.0025,threads=THREADS,scope='S0 fair same700physicalmesh/Q3/globalQ2/noRANS comparison; unavoidable physics/BC differences remain disclosed'))
    assert code==0 and len(samples)==8,'explicit-opt S0 failed: preserve new evidence'

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();run(a.root)
