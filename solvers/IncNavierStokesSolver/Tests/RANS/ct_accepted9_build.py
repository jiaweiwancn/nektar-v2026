#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Private, explicitly Eigen-enabled two-job CT9 runtime; no installation."""
import argparse,json,os,resource,shutil,subprocess
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
from ct_accepted_open import run,require
from ct_accepted_feasibility import write_json
from rans_request_manifest import RequestManifest,runtime_identity,sha

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','ct8','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));s=a.root/'source';build=a.root/'build';rt=a.output/'runtime';previous=a.ct8/'build-v3/runtime'
    g=['git','-c','core.preloadIndex=false','-c','index.threads=1','-c','core.filemode=false'];m=RequestManifest(a.output,dict(schema='CT9SparseBuildV1',revision=subprocess.check_output(g+['rev-parse','HEAD'],cwd=s,text=True).strip(),status=subprocess.check_output(g+['status','--porcelain'],cwd=s,text=True),jobs=2,source=sha(Path(__file__))))
    dep=a.root/'dependency';dep.mkdir(exist_ok=True);spatial=dep/'libSpatialDomains.so.5.9.0'
    if not spatial.exists():shutil.copy2(previous/'lib/libSpatialDomains.so.5.9.0',spatial)
    executions=[];targets=['IncNavierStokesSolver','RANSPreservationProbe','RANSNewModelsProbe','RANSCurvedAcceptedKernelProbe','CT9SparseAlgebraChecks']
    for i,cmd in enumerate([['cmake','-S',str(Path(__file__).parent),'-B',str(build),'-DCMAKE_BUILD_TYPE=Release','-DRANS_BUILD_SOLVER=ON','-DRANS_ENABLE_SERIAL_SPARSE=ON','-DRANS_PREBUILT_SPATIAL_LIBRARY='+str(spatial),'-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake'],['cmake','--build',str(build),'--parallel','2','--target',*targets]]):
        e=run(cmd,a.output,str(i)+'.log');executions.append(e);write_json(a.output/'progress.json',dict(executions=executions));require(e['exit_code']==0,'CT9 isolated build failed; preserve log:'+str(i))
    rt.mkdir();shutil.copytree(previous/'lib',rt/'lib',symlinks=True)
    for name in targets:shutil.copy2(build/name,rt/name)
    os.environ['LD_LIBRARY_PATH']=str(rt/'lib')
    result=dict(schema='CT9SparseBuildV1',status='pass',executions=executions,runtime={n:runtime_identity(rt/n) for n in targets},compiler=subprocess.check_output(['c++','--version'],text=True),flags=(build/'CMakeFiles/IncNavierStokesSolver.dir/flags.make').read_text(),Eigen_headers={str(f):sha(f) for f in (Path('/usr/include/eigen3/Eigen/SparseLU'),Path('/usr/include/eigen3/Eigen/SparseCholesky'),Path('/usr/include/eigen3/Eigen/src/Core/util/Macros.h'))},threads=THREADS)
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);return 0
if __name__=='__main__':raise SystemExit(main())
