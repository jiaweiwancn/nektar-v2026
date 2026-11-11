#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Isolated CT2 build/freeze. Invoke through ct_accepted_command.py."""
import argparse,os,subprocess,time,resource,shutil,json
from pathlib import Path
from rans_qualification import THREADS,write_json
os.environ.update(THREADS)
from rans_request_manifest import sha,runtime_identity
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--previous',type=Path,required=True);a=p.parse_args()
    r=a.root.resolve();s=r/'source/solvers/IncNavierStokesSolver/Tests/RANS'
    (r/'logs').mkdir(exist_ok=True);(r/'dependency').mkdir(exist_ok=False)
    dep=r/'dependency/libSpatialDomains.so.5.9.0';shutil.copy2(a.previous/'dependency'/dep.name,dep)
    os.environ.update(TMPDIR=str(r/'tmp'),LD_LIBRARY_PATH=str(dep.parent))
    commands=[['cmake','-S',str(s),'-B',str(r/'build'),'-DCMAKE_BUILD_TYPE=Release','-DRANS_BUILD_SOLVER=ON','-DRANS_PREBUILT_SPATIAL_LIBRARY='+str(dep),'-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake'],
      ['cmake','--build',str(r/'build'),'--parallel','4','--target','IncNavierStokesSolver','RANSPreservationProbe','RANSNewModelsProbe','RANSNativeDiagnosticChecks']]
    receipts=[]
    for n,c in enumerate(commands):
        t=time.monotonic();b=resource.getrusage(resource.RUSAGE_CHILDREN)
        with (r/'logs'/('build-'+str(n)+'.log')).open('x') as f:code=subprocess.run(c,stdout=f,stderr=subprocess.STDOUT).returncode
        e=resource.getrusage(resource.RUSAGE_CHILDREN)
        receipts.append(dict(command=c,exit_code=code,wall_seconds=time.monotonic()-t,process_cpu_seconds=e.ru_utime+e.ru_stime-b.ru_utime-b.ru_stime))
        write_json(r/'build-execution.json',receipts)
        if code:return code
    rt=r/'runtime';rt.mkdir()
    shutil.copytree(a.previous/'runtime-v2/lib',rt/'lib',symlinks=True)
    for name in ('IncNavierStokesSolver','RANSPreservationProbe','RANSNewModelsProbe','RANSNativeDiagnosticChecks'):shutil.copy2(r/'build'/name,rt/name)
    os.environ['LD_LIBRARY_PATH']=str(rt/'lib')
    write_json(r/'runtime-identity.json',dict(source_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=r/'source',text=True).strip(),
        source_status=subprocess.check_output(['git','status','--porcelain'],cwd=r/'source',text=True),runtimes={n:runtime_identity(rt/n) for n in ('IncNavierStokesSolver','RANSPreservationProbe','RANSNewModelsProbe','RANSNativeDiagnosticChecks')},
        reused_frozen_libraries={p.name:sha(p) for p in (rt/'lib').iterdir() if p.is_file()},dependency_sha256=sha(dep),threads=THREADS,compiler=subprocess.check_output(['c++','--version'],text=True),build_receipt_sha256=sha(r/'build-execution.json')))
    return 0
if __name__=='__main__':raise SystemExit(main())
