#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Affected-target build; old objects reused only with exact dependency proof."""
import argparse,hashlib,os,re,shlex,shutil,subprocess,time
from pathlib import Path
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,verify,atomic

OLD=Path('/home/jay/build/wj-h19-target-timestep-20261115')

def build(root):
    start=time.monotonic();source=root/'source';build=root/'build'
    cache={}
    for line in (OLD/'build/CMakeCache.txt').read_text().splitlines():
        if line and not line.startswith(('#','//'))and '='in line:
            left,value=line.split('=',1);cache[left.split(':')[0]]=value
    keys=['CMAKE_BUILD_TYPE','CMAKE_CXX_COMPILER','CMAKE_CXX_FLAGS','CMAKE_CXX_FLAGS_RELEASE',
          'RANS_BUILD_SOLVER','RANS_PREBUILT_SPATIAL_LIBRARY']
    keys.extend(k for k in cache if 'nektar' in k.lower()and k.endswith('_DIR'))
    command=['cmake','-S',str(source/'solvers/IncNavierStokesSolver/Tests/RANS'),'-B',str(build),'-DCMAKE_EXPORT_COMPILE_COMMANDS=ON']
    command.extend('-D'+k+'='+cache[k]for k in keys if k in cache)
    assert cache['RANS_BUILD_SOLVER']=='ON'and '-ffp-contract=off'in cache['CMAKE_CXX_FLAGS']and '-fno-fast-math'in cache['CMAKE_CXX_FLAGS']
    subprocess.run(command,check=True)
    old=load(OLD/'build/compile_commands.json');new=load(build/'compile_commands.json');records=[];digests={}
    def digest(path):
        path=Path(path)
        if path not in digests:digests[path]=identity(path)['sha256']
        return digests[path]
    def key(v):
        args=shlex.split(v['command']);out=args[args.index('-o')+1]
        return out.split('/')[1],str(Path(v['file']).relative_to(OLD/'source'if str(v['file']).startswith(str(OLD))else source))
    previous={key(v):v for v in old};wanted={'IncNavierStokesSolver.dir','RANSWallResolutionProbe.dir','RANSKernels.dir'}
    for v in new:
        k=key(v)
        if k[0]not in wanted or k not in previous:continue
        oldv=previous[k];oa=shlex.split(oldv['command']);na=shlex.split(v['command'])
        normalized=[x.replace(str(OLD),str(root))for x in oa]
        if normalized!=na:continue
        oldobj=OLD/'build'/oa[oa.index('-o')+1];newobj=build/na[na.index('-o')+1];dep=Path(str(oldobj)+'.d')
        if not(oldobj.exists()and dep.exists()):continue
        deps=shlex.split(dep.read_text().replace('\\\n',' ').split(':',1)[1]);proof=[];compatible=True
        for p in deps:
            oldp=Path(p);newp=Path(p.replace(str(OLD/'source'),str(source)))
            if not oldp.is_file()or not newp.is_file()or digest(oldp)!=digest(newp):compatible=False;break
            if str(oldp).startswith(str(OLD/'source')):proof.append(dict(old=identity(oldp),actual=identity(newp)))
        if not compatible:continue
        newobj.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(oldobj,newobj)
        Path(str(newobj)+'.d').write_text(dep.read_text().replace(str(OLD),str(root)))
        records.append(dict(source=v['file'],object=identity(newobj),historical_object=identity(oldobj),historical_command=oldv['command'],
                            actual_equivalent_command=v['command'],source_dependencies=proof))
    atomic(root/'science/object-reuse-proof.json',dict(reused=records,configure_command=command,
          criterion='same normalized compiler flags and exact dependency bytes; changed-header objects are rebuilt',new_native_revision=load(root/'source-intent.json')['native_revision']))
    subprocess.run(['cmake','--build',str(build),'--target','IncNavierStokesSolver','RANSWallResolutionProbe','WJFlowVectorChecks','--parallel','8'],check=True)
    atomic(root/'science/build-result.json',dict(status='PASS',seconds=time.monotonic()-start,configure=command,
        binaries=[identity(build/n)for n in ('IncNavierStokesSolver','RANSWallResolutionProbe','WJFlowVectorChecks')],
        compiler=subprocess.check_output(['/usr/bin/c++','--version'],text=True),reused_objects=len(records)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();build(a.root)
