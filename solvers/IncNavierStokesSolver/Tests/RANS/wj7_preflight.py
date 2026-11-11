#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only final runtime/test/import equivalence and pre-ledger I/O charge."""
import argparse
import math
from pathlib import Path
from wj_immutable import read,write_once,file_identity,file_hash,verify_artifacts,Ledger,digest
from run_cls_epawf_picard import require


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--parent',type=Path,required=True)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--study-name',default='study');a=p.parse_args();root=a.root.resolve();parent=a.parent.resolve()
    aliases=read(parent/'verified-runtime-aliases-v4.json')['aliases'];binaries={};libraries={}
    for path,r in aliases.items():
        require(Path(path).resolve()==Path(r['path']).resolve() and file_hash(path)==r['sha256'],'reused regression executable differs')
        binaries[path]=r['sha256']
        for name,sha in r['libraries'].items():
            require(file_hash(name)==sha,'reused regression dependency differs');libraries[name]=sha
    require(len(binaries)==35,'regression runtime alias coverage differs')
    eq=read(root/'native-equivalence.json')
    for name,sha in eq['native_inputs'].items():require(file_hash(a.source/name)==sha,'native source input differs')
    manifest=read(parent/'compact-v3/manifest.json')
    for name,wanted in manifest['effective_flags'].items():
        require((parent/'build-v2/CMakeFiles'/name/'flags.make').read_text()==wanted,'effective compiler flags differ')
    for r in eq['runtime'].values():
        require(file_hash(r['path'])==r['sha256'],'scientific runtime changed')
        for name,sha in r['libraries'].items():require(file_hash(name)==sha,'scientific dependency changed')
    imports={}
    for policy in ('EP','Budget'):
        old=read(parent/'compact-v3'/(policy+'.json'));new=read(root/a.study_name/policy/'imports.json')
        require(new['status']=='PASS' and set(new['endpoints'])=={'coarse','fine','baseline'},'six-endpoint import incomplete')
        design=read(root/a.study_name/policy/'design.json')
        for r in (design['runtime'][str(Path(manifest['runtime'][n]['path']).resolve())] for n in ('IncNavierStokesSolver','RANSWallResolutionProbe')):
            require(any(r==v for v in eq['runtime'].values()),'imported runtime does not match reviewed native producer')
        for tag,t in old['trajectories'].items():
            for c in (t['initial'],t['trajectory'][-1]):
                require(file_hash(Path(c['checkpoint'])/'State.bin')==c['checkpoint_sha256'] and file_hash(c['report'])==c['report_sha256'],
                        'original endpoint/CP0/report bytes changed during import')
        imports[policy]=file_identity(root/a.study_name/policy/'imports.json')
    # Earlier source worktree checkouts preceded command monitoring. Derive
    # their filesystem interval and add a fixed120s allowance for unwrapped
    # fetch/SSH/load inspection. This conservative charge is written once.
    checkout=[]
    for n in range(1,6):
        f=root/('source-v'+str(n))/'.git';gitdir=Path(f.read_text().split(': ',1)[1].strip());index=gitdir/'index'
        checkout.append(dict(source=str(f.parent),start_epoch=f.stat().st_mtime,end_epoch=index.stat().st_mtime,
            observed_interval_seconds=max(0,index.stat().st_mtime-f.stat().st_mtime)))
    seconds=math.ceil(sum(r['observed_interval_seconds'] for r in checkout)*1.25+120)
    charge=dict(schema='WJ7PreledgerSetupAllowanceV1',checkout_intervals=checkout,headroom=1.25,
        additional_unwrapped_fetch_SSH_inspection_allowance_seconds=120,charged_seconds_per_policy=seconds,
        scope='reconstructed conservative setup/I/O charge; monitored source-v6-and-later work is charged by its actual receipts')
    write_once(root/'setup-cost-allowance.json',charge)
    for policy in ('EP','Budget'):
        l=Ledger(root,policy,cap=800);label='unwrapped-initial-source-and-inspection'
        l.reserve(label,charge,0,seconds);l.finish(label,0,seconds,str(root/'setup-cost-allowance.json'))
    value=dict(schema='WJ7FinalPreflightV1',status='PASS',regression_executables=binaries,regression_libraries=libraries,
        reused_test_count=38,reused_regression_evidence=file_identity(parent/'compact-v3/regression-verification.json'),
        imports=imports,native_equivalence=file_identity(root/'native-equivalence.json'),setup_allowance=file_identity(root/'setup-cost-allowance.json'),
        source=str(a.source.resolve()),code=file_identity(__file__))
    write_once(root/'preflight.json',value);print(value['status'],seconds)


if __name__=='__main__':main()
