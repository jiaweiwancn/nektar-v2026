#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Capture read-only identities for an isolated WJ-1 source/build/run root."""
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys

root=Path(sys.argv[1]).resolve()
source,build=root/'source',root/'build'

def command(*args):
    return subprocess.check_output(args,text=True).strip()

def digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):
            h.update(block)
    return h.hexdigest()

entries=subprocess.check_output(['git','-C',str(source),'ls-tree','-rz','HEAD']).split(b'\0')
manifest=[]
changed=[]
for entry in entries:
    if not entry: continue
    header,path=entry.split(b'\t',1)
    mode,kind,oid=header.split()
    name=os.fsdecode(path)
    file=source/name
    if mode==b'120000' and file.is_symlink():
        content=os.fsencode(os.readlink(file))
    else:
        content=file.read_bytes()
    actual=hashlib.sha1(b'blob '+str(len(content)).encode()+b'\0'+content).hexdigest()
    if actual!=oid.decode(): changed.append(name)
    manifest.append(f'{hashlib.sha256(content).hexdigest()}  {name}\n')
manifest.sort()
(root/'source-content.sha256').write_text(''.join(manifest))
assert not changed, changed
executables=['WJKernelProbe','RANSKernelProbe','RANSConstitutiveProbe','RANSAcceptedAWFProbe',
    'RANSEquilibriumProfilesProbe','RANSProfileBranchChecks','RANSKBoundQPChecks',
    'RANSHistoryRejectionChecks','RANSStressTangentAudit']
identities={}
libraries={}
for name in executables:
    path=build/name
    links=command('ldd',str(path))
    identities[name]=dict(path=str(path),sha256=digest(path),ldd=links)
    for lib in re.findall(r'(/\S+)\s+\(0x[0-9a-fA-F]+\)',links):
        actual=Path(lib).resolve()
        libraries[str(actual)]=digest(actual)
inventory=json.loads((root/'test-inventory.json').read_text())
registered={t['name']:t for t in inventory['tests']}
required=['WJKernelContracts','WJIndependentOracle','RANSKernelUnitTests','RANSCLSEPAWFM1Oracle',
 'RANSCLSEPAWFM2Oracle','RANSAcceptedAWFKernel','RANSEquilibriumProfilesNative','RANSProfileBranchChecks',
 'RANSKBoundQPChecks','RANSHistoryRejectionChecks','RANSAcceptedAWFCases','RANSStressTangentOracle']
assert set(required)<=registered.keys()
flags={name:(build/'CMakeFiles'/f'{name}.dir'/'flags.make').read_text()
       for name in ('RANSKernels','WJKernelProbe')}
report=dict(revision=command('git','-C',str(source),'rev-parse','HEAD'),
    git_tree=command('git','-C',str(source),'rev-parse','HEAD^{tree}'),
    content_files=len(manifest),source_manifest_sha256=digest(root/'source-content.sha256'),
    changed_content=changed,working_tree_ignoring_ntfs_modes=command('git','-c','core.filemode=false','-C',str(source),'status','--porcelain'),
    platform=platform.platform(),hostname=platform.node(),compiler=command('/usr/bin/c++','--version'),
    cmake=command('cmake','--version'),python=sys.version,
    specs={name:digest(source/'docs/literature'/name) for name in
       ('wj_earsm_neutral_abl_ep_awf_note.tex','wj_earsm_budget_cutoff_awf_appendix.tex')},
    registered_tests=sorted(registered),required_tests=required,
    required_test_properties={name:registered[name].get('properties',[]) for name in required},
    flags=flags,executables=identities,linked_libraries_sha256=libraries,
    static_kernels_sha256=digest(build/'libRANSKernels.a'),
    threads={name:os.environ.get(name) for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')})
(root/'identities.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:report[k] for k in ('revision','content_files','source_manifest_sha256','changed_content','threads')},indent=2))
