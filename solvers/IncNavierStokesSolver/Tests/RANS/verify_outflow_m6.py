#!/usr/bin/env python3
"""Compare native and matrix-free total-traction outlets on complete fields."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET
sys.dont_write_bytecode=True
from generate_m6 import generate_abl
from generate_m4 import parameter


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if not sys.platform.startswith('linux'):raise RuntimeError('Use Linux for authoritative solver comparisons')
    build=args.build.resolve();output=args.output.resolve();output.mkdir(parents=True,exist_ok=True)
    solver=build/'IncNavierStokesSolver';audit=build/'RANSFieldAudit'
    env=os.environ.copy();env.update({key:'1' for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')})
    commands=[];results={}
    def execute(command,cwd,label):
        command=list(map(str,command));commands.append(dict(command=command,cwd=str(cwd)))
        result=subprocess.run(command,cwd=cwd,env=env,capture_output=True,text=True,timeout=1800)
        text=result.stdout+result.stderr;(cwd/(label+'.log')).write_text(text)
        if result.returncode:raise RuntimeError(str(cwd/(label+'.log')))
        return text
    for dim in (2,3):
        for model in ('standard','cls_fixed'):
            pair=[]
            for backend in ('static','matrixfree'):
                name=f'outlet_d{dim}_{model}_{backend}';directory=output/name;directory.mkdir(exist_ok=True)
                session=directory/(name+'.xml')
                generate_abl(session,model=model,dimension=dim,degree=3,layers=4,nx=1,steps=2,dt=.5,
                    backend=backend,boundary='outflow',monitor=True,tolerance=1e-8)
                tree=ET.parse(session);parameter(tree.getroot(),'RANSMaxIterations',10000);tree.write(session,encoding='utf-8',xml_declaration=True)
                print('Running '+name,flush=True)
                command=(['mpirun','--bind-to','none','-np','4'] if dim==3 else [])+[solver,session.name]
                execute(command,directory,'solver');pair.append((session,directory/(name+'.fld')))
            raw=execute([audit,pair[0][0],pair[0][1],pair[1][1]],output,f'compare_d{dim}_{model}')
            value=json.loads(raw[raw.index('{"time"'):])
            largest=max(f['difference_linf'] for f in value['fields'].values())
            if largest>=1e-6:raise RuntimeError('Outlet backend difference exceeds 1e-6: '+str(largest))
            results[f'd{dim}_{model}']=value
            print(f'd{dim} {model}: max difference {largest:.4g}',flush=True)
    report=dict(status='passed',maximum_workers=4,commands=commands,comparisons=results,
        binaries={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (solver,audit)})
    (output/'results.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
