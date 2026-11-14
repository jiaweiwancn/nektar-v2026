#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Seal compact CT9 publication after the numerical collector has exited.

The outer owned launcher closes after this script. Copy its closed receipt
separately; do not create a circular archive that claims to contain its own
future resource receipt. No fields, matrices or native histories enter the tar.
"""
import argparse,json,os,shutil,subprocess,tarfile
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
from rans_request_manifest import RequestManifest,sha
from ct_accepted_feasibility import write_json,require

def read(path):return json.loads(path.read_text())

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','qualified','collector','build','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();source=a.root/'source';request=read(a.qualified/'request-manifest.json');seal=read(a.qualified/'completed-artifacts.json')
    require(seal['request_sha256']==request['sha256'],'foreign collected request')
    for name,digest in seal['files'].items():require(sha(a.qualified/name)==digest,'changed compact numerical evidence:'+name)
    collector=read(a.collector/'execution.json');require(collector['exit_code']==0 and collector['resource']['project_status']==collector['resource']['aggregate_status']=='pass','collector not closed/resource compliant')
    g=['git','-c','core.preloadIndex=false','-c','index.threads=1','-c','core.filemode=false'];revision=subprocess.check_output(g+['rev-parse','HEAD'],cwd=source,text=True).strip();status=subprocess.check_output(g+['status','--porcelain'],cwd=source,text=True)
    require(not status,'publication source is dirty')
    m=RequestManifest(a.output,dict(schema='CT9PublicationV1',qualified_seal=sha(a.qualified/'completed-artifacts.json'),collector_execution=sha(a.collector/'execution.json'),build=sha(a.build/'results.json'),source_revision=revision,source=sha(Path(__file__)),scope='compact publication; numerical outcomes unchanged; outer launcher closure bound separately after exit'))
    closed=a.output/'collector-resource';closed.mkdir()
    for f in a.collector.iterdir():
        if f.is_file() and f.suffix in ('.json','.jsonl'):shutil.copy2(f,closed/f.name)
    for name in ('ctest-v4.log','ctest-v4.xml'):shutil.copy2(a.root/name,a.output/name)
    sdk=Path('/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake/Nektar++Config.cmake');require(sdk.is_file(),'actual SDK config missing')
    infrastructure=['library/LibUtilities/LinearAlgebra/'+n for n in ('SparseMatrix.hpp','SparseUtils.hpp','NekLinSysIter.h','NekLinSysIterGMRES.h','NekLinSysIterGMRES.cpp')]
    config_lines=[line for line in sdk.read_text().splitlines() if any(k in line.upper() for k in ('PETSC','NEKTAR_USE_MPI','EIGEN'))]
    cache=Path('/home/jay/build/nektar-v2026-mpi/CMakeCache.txt');flags=[line for line in cache.read_text().splitlines() if line.startswith(('CMAKE_BUILD_TYPE:','NEKTAR_USE_MPI:','NEKTAR_USE_PETSC:'))]
    require('NEKTAR_USE_PETSC:BOOL=OFF' in flags and 'NEKTAR_USE_MPI:BOOL=ON' in flags,'installed build-cache dependency flags changed')
    package=subprocess.run(['dpkg-query','-W','-f=${Version}','libeigen3-dev'],capture_output=True,text=True)
    dependency=dict(installed_sdk_config=dict(path=str(sdk),sha256=sha(sdk),selected_lines=config_lines),installed_build_cache=dict(path=str(cache),sha256=sha(cache),flags=flags),eigen_package=dict(exit_code=package.returncode,version=package.stdout.strip(),stderr=package.stderr.strip()),infrastructure={n:sha(source/n) for n in infrastructure},choice='private Eigen SparseLU/positive LDLT and installed Nektar GMRES callbacks; direct CSR adapter for combined regional/sampler links; no PETSc or shared installation change',audit_note='SDK package exports MPI but omits PETSc; its actual OFF status is verified in the installed build cache, alongside executable/loaded-library hashes')
    write_json(a.output/'dependency.json',dependency)
    baseline='0c28610e685aef95f7e6f21b9a4fdcd305dbf319';unchanged=['library/SolverUtils/Turbulence/RANSModels.cpp','library/SolverUtils/Turbulence/EquilibriumAWFProfiles.cpp','docs/literature/accepted_awf_appendix_A93_A95_E94_E95.tex','docs/literature/A93_A95_E94_E95_neutral_ABL_EP_AWF_combined.tex','docs/plans/cao-tamura-accepted-8-evidence/physical-reference-manifest.json']
    import hashlib
    hashes={}
    for name in unchanged:
        old=hashlib.sha256(subprocess.check_output(g+['show',baseline+':'+name],cwd=source)).hexdigest();require(old==sha(source/name),'unintended physical/specification change:'+name);hashes[name]=old
    scope=dict(baseline=baseline,publication_controller_revision=revision,compiled_revision=read(a.build/'request-manifest.json')['request']['revision'],unchanged=hashes,history='fresh sparse component-bound histories; old physical histories reused read-only only at the pinned dense identity; no signature edits',physical_continuation=False,decision=read(a.qualified/'decision.json'))
    write_json(a.output/'source-scope.json',scope);write_json(a.output/'files.json',{str(f.relative_to(a.output)):sha(f) for f in a.output.rglob('*') if f.is_file() and f.name not in ('request-manifest.json','completed-artifacts.json','files.json')})
    m.seal(a.output,[str(f.relative_to(a.output)) for f in a.output.rglob('*') if f.is_file() and f.name not in ('request-manifest.json','completed-artifacts.json')])
    archive=a.root/(a.output.name+'.tar.gz');require(not archive.exists(),'publication archive already exists')
    with tarfile.open(archive,'w:gz') as out:
        out.add(a.qualified,arcname=a.qualified.name);out.add(a.output,arcname=a.output.name)
    print('CT9 publication',archive,sha(archive),flush=True);return 0
if __name__=='__main__':raise SystemExit(main())
