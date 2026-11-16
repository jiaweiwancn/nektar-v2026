#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Isolated H25 affected build, with exact old-object dependency qualification."""
import argparse,hashlib,os,re,shlex,shutil,subprocess,time
from pathlib import Path
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,verify,atomic

OLD=Path('/home/jay/build/wj-h22-early-phase-replay-20261115')

def build(root,variant='profile'):
    if variant=='profile':
        source=root/'source-v4';output=root/'build-profile-v2';previous=OLD/'build';previous_source=OLD/'source';tag='profile';manifest='source-manifest-v4.json'
    elif variant=='lookup':
        source=root/'source-v5';output=root/'build-lookup-v1';previous=root/'build-profile-v2';previous_source=root/'source-v4';tag='lookup';manifest='source-manifest-v5.json'
    elif variant=='seed':
        source=root/'resumption-r1/source-v1';output=root/'resumption-r1/build-seed-v1';previous=root/'build-scalar-v3';previous_source=root/'source-v8';tag='seed';manifest='resumption-r1/source-manifest-v1.json'
    else:
        assert variant=='scalar'
        source=root/'source-v8';output=root/'build-scalar-v3';previous=root/'build-scalar-v2';previous_source=root/'source-v8';tag='scalar';manifest='source-manifest-v8.json'
    started=time.monotonic();proof=(root/'resumption-r1/science' if variant=='seed' else root/'science')/(output.name+'-object-reuse.json')
    assert not output.exists()and not proof.exists(),'versioned build output/receipt already exists'
    cmake=(source/'solvers/IncNavierStokesSolver/Tests/RANS/CMakeLists.txt').read_text()
    for relative in re.findall(r'\$\{RANS_SOURCE_ROOT\}/([A-Za-z0-9_./-]+\.(?:cpp|h))',cmake):
        assert(source/relative).is_file(),'missing configured source fixture: '+relative
    cache={}
    for line in(OLD/'build/CMakeCache.txt').read_text().splitlines():
        if line and not line.startswith(('#','//'))and '='in line:
            left,value=line.split('=',1);cache[left.split(':')[0]]=value
    keys=['CMAKE_BUILD_TYPE','CMAKE_CXX_COMPILER','CMAKE_CXX_FLAGS','CMAKE_CXX_FLAGS_RELEASE','RANS_BUILD_SOLVER','RANS_PREBUILT_SPATIAL_LIBRARY']
    keys.extend(k for k in cache if k.endswith('_DIR')and 'nektar'in k.lower())
    assert cache['RANS_BUILD_SOLVER']=='ON'and cache['CMAKE_BUILD_TYPE']=='Release'
    assert all(v in cache['CMAKE_CXX_FLAGS']for v in ('-ffp-contract=off','-fno-fast-math'))
    config=['cmake','-S',str(source/'solvers/IncNavierStokesSolver/Tests/RANS'),'-B',str(output),'-DCMAKE_EXPORT_COMPILE_COMMANDS=ON']
    config+=['-D'+k+'='+cache[k]for k in keys if k in cache]
    subprocess.run(config,check=True)
    old=load(previous/'compile_commands.json');new=load(output/'compile_commands.json');records=[];digests={}
    def digest(p):
        p=Path(p)
        if p not in digests:digests[p]=identity(p)['sha256']
        return digests[p]
    for n in new:
        suffix=str(Path(n['file']).relative_to(source))
        matches=[o for o in old if o['file'].endswith('/'+suffix)]
        if not matches:continue
        o=matches[0];oa,na=shlex.split(o['command']),shlex.split(n['command'])
        flags=lambda a:[v for v in a if v.startswith(('-O','-f','-m','-std=','-D'))]
        if flags(oa)!=flags(na):continue
        oldobj=previous/oa[oa.index('-o')+1];newobj=output/na[na.index('-o')+1];dep=Path(str(oldobj)+'.d')
        if not oldobj.exists()or not dep.exists():continue
        paths=shlex.split(dep.read_text().replace('\\\n',' ').split(':',1)[1]);pairs=[];okay=True
        for path in paths:
            p=Path(path);newpath=Path(str(p).replace(str(previous_source),str(source),1))if p.is_relative_to(previous_source)else p
            # Explicit old private compatibility include roots remain compiler
            # inputs if that generated header is not in the new source tree.
            if not newpath.exists()and '/Tests/RANS/include/'in str(p):newpath=p
            if not p.exists()or not newpath.exists()or digest(p)!=digest(newpath):okay=False;break
            pairs.append(dict(old=str(p),actual=str(newpath),sha256=digest(newpath)))
        if not okay:continue
        newobj.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(oldobj,newobj)
        # CMake's dependency scanner needs the new source paths too.
        Path(str(newobj)+'.d').write_text(str(newobj)+': '+ ' '.join(shlex.quote(p['actual'])for p in pairs)+'\n')
        now=time.time()+1;os.utime(newobj,(now,now))
        records.append(dict(source=identity(n['file']),object=identity(newobj),old_object=identity(oldobj),dependencies=pairs))
    atomic(proof,dict(status='QUALIFIED_EXACT_DEPENDENCIES',objects=records,configure=config,cache=identity(OLD/'build/CMakeCache.txt')))
    names=('IncNavierStokesSolver','RANSWallResolutionProbe','WJH25CommonMassChecks','RANSWallTimingChecks','RANSWallMonitorChecks')
    if variant=='seed':names=('RANSWallResolutionProbe',)
    if variant in ('lookup','scalar'):names+=('WJH25ExactLookupChecks',)
    command=['cmake','--build',str(output),'--target',*names,'--parallel','8']
    subprocess.run(command,check=True)
    runtime=(root/'resumption-r1'/('runtime-'+tag+'-v1')) if variant=='seed' else root/('runtime-'+tag+'-v1');runtime.mkdir()
    oldrt=OLD/'runtime';executables=[]
    for name in names:
        elf=runtime/(name+'.elf');shutil.copy2(output/name,elf)
        wrapper=runtime/name
        text='#!/bin/sh\nexec '+shlex.quote(str(oldrt/'glibc/ld-linux-x86-64.so.2'))+' --library-path '+shlex.quote(str(oldrt/'glibc')+':'+str(oldrt/'geometry')+':/home/jay/.local/opt/nektar-v2026-mpi/lib64')+' '+shlex.quote(str(elf))+' "$@"\n'
        wrapper.write_text(text);wrapper.chmod(0o755);executables.append(dict(ELF=identity(elf),wrapper=identity(wrapper),build=identity(output/name)))
    objects=[]
    for c in new:
        a=shlex.split(c['command']);obj=output/a[a.index('-o')+1]
        if obj.exists():objects.append(dict(source=identity(c['file']),command=c['command'],object=identity(obj),dependency=identity(str(obj)+'.d')if Path(str(obj)+'.d').exists()else None))
    runtime_record=root/('resumption-r1/science/runtime-build-seed-v1.json') if variant=='seed' else root/('science/runtime-build-'+tag+'-v1.json')
    atomic(runtime_record,dict(status='PASS',source_manifest=identity(root/manifest),
        compiler=subprocess.check_output(['/usr/bin/c++','--version'],text=True),configure=config,build=command,seconds=time.monotonic()-started,
        compile_commands=identity(output/'compile_commands.json'),objects=objects,executables=executables,
        link_commands=[identity(p)for p in sorted(output.glob('CMakeFiles/*/link.txt'))if p.parent.name in ('IncNavierStokesSolver.dir','RANSWallResolutionProbe.dir','WJH25CommonMassChecks.dir')],
        private_runtime=identity(root/'science/runtime-dependencies-v1.json'),old_object_proof=identity(proof),purpose='H25R1 target initial-data producer; existing governing solver unchanged' if variant=='seed' else 'profiling-only baseline'if variant=='profile'else'opt-in exact containment/geometry lookup and scalar value rows; same governing mathematics; explicit H25 transient initial-data helper'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--variant',choices=('profile','lookup','scalar','seed'),default='profile');a=p.parse_args();build(a.root,a.variant)
