#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""S0 using supported NoCollection after diagnosed autotuner allocator failure."""
import argparse,json,os,re,subprocess,time
from pathlib import Path
import xml.etree.ElementTree as ET
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,verify,atomic
from cls_epawf_m3_cases import write
from study_wj_h25 import executable

def run(root):
    old=root/'science/S0-v2/case.xml';folder=root/'science/S0-NoCollection';folder.mkdir()
    x=ET.parse(old).getroot()
    assert x.find('COLLECTIONS')is None
    ET.SubElement(x,'COLLECTIONS',DEFAULT='NoCollection');case=folder/'case.xml';write(x,case)
    atomic(folder/'execution-diff.json',dict(original=identity(old),actual=identity(case),only_change='COLLECTIONS DEFAULT=NoCollection',
        cause='GDB abort in CollectionOptimisation::SetWithTimings before any advance',SVV=False,
        parity_scope='same stock equations, physical mesh, fields, dt and boundary Functions; supported execution backend'))
    exe=executable(root,'IncNavierStokesSolver');maps=None;began=time.monotonic()
    with(folder/'native.log').open('x')as log:
        p=subprocess.Popen([str(exe),str(case)],stdout=log,stderr=subprocess.STDOUT)
        while p.poll()is None:
            try:
                lines=Path(f'/proc/{p.pid}/maps').read_text().splitlines()
                files=sorted({s.split(maxsplit=5)[5]for s in lines if len(s.split(maxsplit=5))==6 and s.split(maxsplit=5)[5].startswith('/')and not s.endswith(' (deleted)')})
                if any('IncNavierStokesSolver.elf'in f for f in files)and any(f.endswith('/libc.so.6')for f in files):maps=files
            except FileNotFoundError:pass
            time.sleep(.025)
        code=p.returncode
    text=(folder/'native.log').read_text()
    samples=[dict(step=int(s),time=float(t),seconds=float(w))for s,t,w in re.findall(r'Steps:\s*(\d+)\s+Time:\s*([\d.eE+-]+)\s+CPU Time:\s*([\d.eE+-]+)s',text)]
    atomic(folder/'result.json',dict(status='PASS'if code==0 and len(samples)==8 else'FAIL',exit=code,wall_seconds=time.monotonic()-began,
        samples=samples,SVV=False,ELF=load(root/'science/runtime-build-v1.json')['executables'][0],maps=maps,
        private_libc=bool(maps)and [f for f in maps if f.endswith('/libc.so.6')]==['/home/jay/build/wj-h22-early-phase-replay-20261115/runtime/glibc/libc.so.6'],
        log=identity(folder/'native.log'),input=identity(case),scope='same700physicalmesh, Q3/globalQ2, dt.0025, same.04→.06 interval, one worker'))
    assert code==0 and len(samples)==8,'stock NoCollection failure; preserve input/output, no silent stabilization'
    assert all('nan'not in s.lower()for s in text.splitlines()if 'error'in s.lower())

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();run(a.root)
