#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Isolated two-job CT4 native solver/probe build; no shared install mutation."""
import argparse,json,os,shutil,subprocess
from pathlib import Path
from rans_qualification import THREADS,write_json
os.environ.update(THREADS)
from ct_accepted_open import run,require
from rans_request_manifest import RequestManifest,runtime_identity,sha

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','ct2','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();source=a.root/'source';g=['git','-c','index.threads=1','-c','core.preloadIndex=false']
    m=RequestManifest(a.output,dict(schema='CT4NativeBuildV1',source_revision=subprocess.check_output(g+['rev-parse','HEAD'],cwd=source,text=True).strip(),
        source_status=subprocess.check_output(g+['status','--porcelain'],cwd=source,text=True),source_sha256=sha(Path(__file__)),compiler_jobs=2))
    os.environ['CMAKE_BUILD_PARALLEL_LEVEL']='2';result=dict(status='fail',executions=[])
    cmds=[['cmake','-S',str(Path(__file__).parent),'-B',str(a.output/'build'),'-DCMAKE_BUILD_TYPE=Release','-DRANS_BUILD_SOLVER=ON',
        '-DRANS_PREBUILT_SPATIAL_LIBRARY='+str(a.ct2/'dependency/libSpatialDomains.so.5.9.0'),'-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake'],
        ['cmake','--build',str(a.output/'build'),'--parallel','2','--target','IncNavierStokesSolver','RANSPreservationProbe','RANSNewModelsProbe','RANSNativeDiagnosticChecks']]
    for i,cmd in enumerate(cmds):
        e=run(cmd,a.output,'build-'+str(i)+'.log');result['executions'].append(e);write_json(a.output/'progress.json',result)
        require(e['exit_code']==0,'CT4_native_build_failure')
    rt=a.output/'runtime';rt.mkdir();shutil.copytree(a.ct2/'runtime/lib',rt/'lib',symlinks=True)
    for name in ('IncNavierStokesSolver','RANSPreservationProbe','RANSNewModelsProbe','RANSNativeDiagnosticChecks'):shutil.copy2(a.output/'build'/name,rt/name)
    os.environ['LD_LIBRARY_PATH']=str(rt/'lib')
    result.update(status='pass',runtimes={n:runtime_identity(rt/n) for n in ('IncNavierStokesSolver','RANSPreservationProbe','RANSNewModelsProbe','RANSNativeDiagnosticChecks')},
        libraries={p.name:sha(p) for p in (rt/'lib').iterdir() if p.is_file()},scope='new continuous-reference input provider only; original governing closures, wall kernel, ownership and PIMPLE method unchanged')
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);print('CT4 build pass');return 0
if __name__=='__main__':raise SystemExit(main())
