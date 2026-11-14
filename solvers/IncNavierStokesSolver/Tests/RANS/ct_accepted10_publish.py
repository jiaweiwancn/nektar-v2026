#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Seal compact CT10 publication after the scientific collector has exited."""
import argparse,json,os,shutil,subprocess,tarfile
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
from rans_request_manifest import RequestManifest,sha
from ct_accepted_feasibility import write_json,require


def read(p):return json.loads(p.read_text())


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','qualified','collector','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();request=read(a.qualified/'request-manifest.json');seal=read(a.qualified/'completed-artifacts.json');require(request['sha256']==seal['request_sha256'],'foreign numerical collection')
    for name,digest in seal['files'].items():require(sha(a.qualified/name)==digest,'changed collected evidence:'+name)
    closed=read(a.collector/'execution.json');require(closed['exit_code']==0 and closed['resource']['aggregate_status']==closed['resource']['project_status']=='pass','collector resource receipt not closed and passing')
    source=a.root/'source';git=['git','-c','core.preloadIndex=false','-c','index.threads=1','-c','core.filemode=false'];revision=subprocess.check_output(git+['rev-parse','HEAD'],cwd=source,text=True).strip();require(not subprocess.check_output(git+['status','--porcelain'],cwd=source,text=True),'dirty publication source')
    m=RequestManifest(a.output,dict(schema='CT10PublicationV1',qualified_seal=sha(a.qualified/'completed-artifacts.json'),collector=sha(a.collector/'execution.json'),controller_revision=revision,source=sha(Path(__file__)),scope='compact evidence only; this publisher outer resource receipt is bound separately after it exits'))
    dest=a.output/'collector-resource';dest.mkdir()
    for f in a.collector.iterdir():
        if f.is_file() and f.suffix in ('.json','.jsonl'):shutil.copy2(f,dest/f.name)
    mount=subprocess.check_output(['findmnt','--json','-T',str(a.root)],text=True);sdk=Path('/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake/Nektar++Config.cmake');cache=Path('/home/jay/build/nektar-v2026-mpi/CMakeCache.txt');flags=[s for s in cache.read_text().splitlines() if s.startswith(('NEKTAR_USE_MPI:','NEKTAR_USE_PETSC:','CMAKE_BUILD_TYPE:'))]
    require('NEKTAR_USE_MPI:BOOL=ON' in flags and 'NEKTAR_USE_PETSC:BOOL=OFF' in flags,'runtime dependency build flags changed')
    write_json(a.output/'provenance.json',dict(source_revision=revision,source_root=str(source),experiment_root=str(a.root),mount=json.loads(mount),sdk_config=dict(path=str(sdk),sha256=sha(sdk)),installed_build_cache=dict(path=str(cache),sha256=sha(cache),flags=flags),compiler=subprocess.check_output(['c++','--version'],text=True),libraries='complete actual mapped-file catalogs and executable identities are in qualified/runtime and loaded-libraries',governing='sealed CT9 build-v8,unchanged; new builds are read-only probes only',historical='CT4–CT9 records,all failed physical qualifications and extraction corners unchanged'))
    files=[str(f.relative_to(a.output)) for f in a.output.rglob('*') if f.is_file() and f.name not in ('request-manifest.json','completed-artifacts.json')];m.seal(a.output,files)
    archive=a.root/(a.output.name+'.tar.gz');require(not archive.exists(),'immutable publication destination exists')
    with tarfile.open(archive,'w:gz') as tar:tar.add(a.qualified,arcname=a.qualified.name);tar.add(a.output,arcname=a.output.name)
    print('CT10 publication',archive,sha(archive),flush=True);return 0
if __name__=='__main__':raise SystemExit(main())
