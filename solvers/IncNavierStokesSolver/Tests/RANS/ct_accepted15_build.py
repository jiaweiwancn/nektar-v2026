#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Isolated CT15 build, two compiler jobs; preserve shared/historical runtimes."""
import argparse,os,resource,shutil,subprocess
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
from ct_accepted_open import run,require
from ct_accepted_feasibility import write_json
from rans_request_manifest import RequestManifest,runtime_identity,sha

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','ct9','output'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--phase',choices=['capacity','runtime'],required=True)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3))
    source=a.root/'source';build=a.root/'build';old=a.ct9/'build-v8/runtime';rt=a.output/'runtime'
    targets=['CT13SparseOrderingProbe','CT9SparseAlgebraChecks'] if a.phase=='capacity' else ['IncNavierStokesSolver','RANSPreservationProbe','CT13SparseOrderingProbe','CT9SparseAlgebraChecks']
    m=RequestManifest(a.output,dict(schema='CT15PrivateBuildV1',phase=a.phase,source_revision=subprocess.check_output(['git','-c','index.threads=1','rev-parse','HEAD'],cwd=source,text=True).strip(),source_sha256=sha(Path(__file__)),governing_ancestor=runtime_identity(old/'IncNavierStokesSolver'),jobs=2,targets=targets))
    executions=[]
    for i,cmd in enumerate([['cmake','-S',str(Path(__file__).parent),'-B',str(build),'-DCMAKE_BUILD_TYPE=Release','-DRANS_BUILD_SOLVER=ON','-DRANS_ENABLE_SERIAL_SPARSE=ON','-DRANS_PREBUILT_SPATIAL_LIBRARY='+str(old/'lib/libSpatialDomains.so.5.9.0'),'-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake'],['cmake','--build',str(build),'--parallel','2','--target',*targets]]):
        e=run(cmd,a.output,str(i)+'.log');executions.append(e);write_json(a.output/'progress.json',dict(executions=executions));require(e['exit_code']==0,'CT15 private build failed:'+str(i))
    rt.mkdir();shutil.copytree(old/'lib',rt/'lib',symlinks=True)
    for name in targets:shutil.copy2(build/name,rt/name)
    shutil.copy2(old/'RANSNewModelsProbe',rt/'RANSNewModelsProbe');os.environ['LD_LIBRARY_PATH']=str(rt/'lib')
    result=dict(schema='CT15PrivateBuildV1',status='pass',executions=executions,runtime={n:runtime_identity(rt/n) for n in targets+['RANSNewModelsProbe']},compiler=subprocess.check_output(['c++','--version'],text=True),scotch_header_sha256=sha('/usr/include/scotch/scotch.h'),scotch_library=runtime_identity(Path('/usr/lib/x86_64-linux-gnu/libscotch.so')),flags=(build/'CMakeFiles/CT13SparseOrderingProbe.dir/flags.make').read_text())
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);return 0
if __name__=='__main__':raise SystemExit(main())
