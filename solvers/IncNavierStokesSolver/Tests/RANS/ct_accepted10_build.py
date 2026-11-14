#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Private read-only probe build; the sealed CT9 governing solver is not rebuilt."""
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
 a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));s=a.root/'source';old=a.ct9/'build-v8/runtime';build=a.root/'build';rt=a.output/'runtime'
 m=RequestManifest(a.output,dict(schema='CT10ProbeBuildV1',revision=subprocess.check_output(['git','-c','core.preloadIndex=false','-c','index.threads=1','rev-parse','HEAD'],cwd=s,text=True).strip(),jobs=2,governing=runtime_identity(old/'IncNavierStokesSolver'),source=sha(Path(__file__))))
 executions=[]
 for i,cmd in enumerate([['cmake','-S',str(Path(__file__).parent),'-B',str(build),'-DCMAKE_BUILD_TYPE=Release','-DRANS_BUILD_SOLVER=ON','-DRANS_ENABLE_SERIAL_SPARSE=ON','-DRANS_PREBUILT_SPATIAL_LIBRARY='+str(old/'lib/libSpatialDomains.so.5.9.0'),'-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake'],['cmake','--build',str(build),'--parallel','2','--target','RANSPreservationProbe']]):
  e=run(cmd,a.output,str(i)+'.log');executions.append(e);write_json(a.output/'progress.json',dict(executions=executions));require(e['exit_code']==0,'CT10 probe build failed:'+str(i))
 rt.mkdir();shutil.copytree(old/'lib',rt/'lib',symlinks=True);shutil.copy2(build/'RANSPreservationProbe',rt/'RANSPreservationProbe');os.environ['LD_LIBRARY_PATH']=str(rt/'lib')
 result=dict(schema='CT10ProbeBuildV1',status='pass',executions=executions,diagnostic=runtime_identity(rt/'RANSPreservationProbe'),governing=runtime_identity(old/'IncNavierStokesSolver'),compiler=subprocess.check_output(['c++','--version'],text=True),flags=(build/'CMakeFiles/RANSPreservationProbe.dir/flags.make').read_text(),scope='only diagnostic target rebuilt; no new governing runtime, numerical component or physical policy')
 write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);return 0
if __name__=='__main__':raise SystemExit(main())
