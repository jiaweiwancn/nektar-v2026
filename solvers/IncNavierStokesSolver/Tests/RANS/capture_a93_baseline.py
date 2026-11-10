#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Capture baseline before implementation; also rerun unchanged CLS regression."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
sys.dont_write_bytecode=True
from run_cls_epawf_picard import THREADS,sha
from continue_cls_epawf_pimple import freeze

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--build',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False);os.environ.update(THREADS)
    def git(*args):return subprocess.check_output(['git','-c','core.filemode=false','-C',str(a.source),*args],text=True).strip()
    record=dict(schema=1,revision=git('rev-parse','HEAD'),status=git('status','--porcelain'),threads=THREADS,ranks=1,commands=[],source_files={})
    for root in (a.source/'library/SolverUtils/Turbulence',a.source/'solvers/IncNavierStokesSolver'):
        for f in root.rglob('*'):
            if f.is_file() and (f.suffix in ('.cpp','.h') or f.name=='CMakeLists.txt'):record['source_files'][str(f.relative_to(a.source))]=sha(f)
    record['source_manifest_sha256']=hashlib.sha256(json.dumps(record['source_files'],sort_keys=True).encode()).hexdigest()
    record['cmake_cache_sha256']=sha(a.build/'CMakeCache.txt')
    (a.output/'CMakeCache.txt').write_bytes((a.build/'CMakeCache.txt').read_bytes())
    def run(command,name):
        with (a.output/(name+'.log')).open('w') as f:r=subprocess.run(list(map(str,command)),stdout=f,stderr=subprocess.STDOUT)
        record['commands'].append(dict(command=list(map(str,command)),exit_code=r.returncode,log=name+'.log'))
        (a.output/'baseline.json').write_text(json.dumps(record,indent=2)+'\n')
        if r.returncode:raise RuntimeError(name+' failed; evidence retained')
    tests=a.source/'solvers/IncNavierStokesSolver/Tests/RANS'
    run(['python3','-B',tests/'validate_cls_epawf_m1.py','--kernel-probe',a.build/'RANSKernelProbe','--constitutive-probe',a.build/'RANSConstitutiveProbe','--output',a.output/'kernel-oracle'],'existing-kernel-oracle')
    run([a.build/'RANSAcceptedAWFProbe','--output',a.output/'accepted-kernel.json'],'accepted-kernel')
    run([a.build/'RANSProfileBranchChecks'],'profile-branches')
    solver,identity=freeze(a.build/'IncNavierStokesSolver',a.output/'runtime')
    record['runtime']=identity
    run(['python3','-B',tests/'run_cls_epawf_pimple.py','--solver',solver,'--output',a.output/'cls-control','--steps','2','--dt','.000625','--diagnostic','--wall-policy','UniformEPFlat'],'native-cls-control')
    report=json.loads((a.output/'cls-control/case.evolution.json').read_text())
    record.update(accepted_steps=report['steps'],failure=report['failure'],initial_native_audit=report.get('state_audits',[])[:1],status='pass' if report['steps']==2 and not report['failure'] else 'fail')
    (a.output/'baseline.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(dict(status=record['status'],steps=report['steps'],failure=report['failure'])))
    return 0 if record['status']=='pass' else 1
if __name__=='__main__':sys.exit(main())
