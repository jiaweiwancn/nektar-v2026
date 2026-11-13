#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Isolated CT6 two-job build and immutable runtime; use the owned launcher."""
import argparse,json,os,shutil,subprocess
from pathlib import Path
from rans_qualification import THREADS,write_json
os.environ.update(THREADS)
from ct_accepted_open import run,require
from rans_request_manifest import RequestManifest,sha,runtime_identity

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','previous','output'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--build',type=Path);a=p.parse_args();s=a.root/'source';build=a.build or a.root/'build'
    g=['git','-c','index.threads=1','-c','core.preloadIndex=false'];m=RequestManifest(a.output,dict(schema='CT6IsolatedBuildV1',revision=subprocess.check_output(g+['rev-parse','HEAD'],cwd=s,text=True).strip(),dirty=subprocess.check_output(g+['status','--porcelain'],cwd=s,text=True),compiler_jobs=2,source_sha256=sha(Path(__file__))))
    dep=a.root/'dependency';dep.mkdir(exist_ok=True);spatial=dep/'libSpatialDomains.so.5.9.0'
    if not spatial.exists():shutil.copy2(a.previous/'build-v3/runtime/lib/libSpatialDomains.so.5.9.0',spatial)
    targets=('IncNavierStokesSolver','RANSPreservationProbe','RANSNewModelsProbe','RANSCurvedAcceptedKernelProbe','CT6ManufacturedJetChecks')
    cmds=[['cmake','-S',str(Path(__file__).parent),'-B',str(build),'-DCMAKE_BUILD_TYPE=Release','-DRANS_BUILD_SOLVER=ON','-DRANS_PREBUILT_SPATIAL_LIBRARY='+str(spatial),'-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake'],['cmake','--build',str(build),'--parallel','2','--target',*targets]]
    out=dict(status='incomplete',executions=[])
    for i,cmd in enumerate(cmds):
        row=run(cmd,a.output,'build-'+str(i)+'.log');out['executions'].append(row);write_json(a.output/'progress.json',out);require(row['exit_code']==0,'CT6_isolated_build_failure')
    rt=a.output/'runtime';rt.mkdir();shutil.copytree(a.previous/'build-v3/runtime/lib',rt/'lib',symlinks=True)
    for name in targets:shutil.copy2(build/name,rt/name)
    os.environ['LD_LIBRARY_PATH']=str(rt/'lib')
    out.update(status='pass',runtimes={n:runtime_identity(rt/n) for n in targets},compiler=subprocess.check_output(['c++','--version'],text=True),flags=(build/'CMakeFiles/IncNavierStokesSolver.dir/flags.make').read_text(),threads=THREADS)
    write_json(a.output/'results.json',out);m.seal(a.output,['results.json']);print('CT6 isolated build pass');return 0
if __name__=='__main__':raise SystemExit(main())
