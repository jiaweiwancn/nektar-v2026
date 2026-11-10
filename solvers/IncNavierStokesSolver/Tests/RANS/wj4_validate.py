#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Registered WJ-3 controls plus WJ-4 checks, with package accounting."""
import argparse
import os
from pathlib import Path
import subprocess
import time
from wj_immutable import Ledger,read,write_once,artifacts,file_hash
from run_cls_epawf_picard import THREADS

CONTROL='WJKernelContracts|WJIndependentOracle|RANSKernelUnitTests|RANSCLSEPAWFM1Oracle|RANSCLSEPAWFM2Oracle|RANSAcceptedAWFKernel|RANSEquilibriumProfilesNative|RANSProfileBranchChecks|RANSKBoundQPChecks|RANSHistoryRejectionChecks|RANSAcceptedAWFCases|RANSStressTangentOracle|RANSReturnedFieldsChecks|RANSObservationKernels|RANSWallTimingChecks|RANSEquilibriumTimeNative|RANSReducedModelChecks|RANSStandardTransportSmoke|RANSCLSTransportSmoke|RANSAWFOriginalSmoke|RANSAWFEquilibriumSmoke|RANSStandardAWFEquilibriumSmoke|RANSEquilibriumWallEvolutionSmoke|RANSEquilibriumWallStageSmoke'
SUITES=dict(controls=(CONTROL,24,0),native=('WJNativeEP|WJNativeBudget|RANSWJSharedLegacyTransaction',3,100),
    readonly=('WJReadOnlyStageEP|WJReadOnlyStageBudget',2,10),metadata=('WJQualificationMetadata',1,0),
    composite=('WJCompositeObservationEP|WJCompositeObservationBudget',2,0))

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--label',default='qualification');p.add_argument('--suite',choices=SUITES,required=True);a=p.parse_args()
    root=a.root.resolve();build=root/'build';scripts=Path(__file__).resolve().parent
    env=os.environ.copy();env.update(THREADS);env['WJ4_PACKAGE_ROOT']=os.environ.get('WJ4_PACKAGE_ROOT',str(root));env['LD_LIBRARY_PATH']=str(build/'geometry-lib')+':'+env.get('LD_LIBRARY_PATH','')
    regex,count,reserve=SUITES[a.suite];regex='^('+regex+')$';label=a.label+'-'+a.suite
    inventory=subprocess.check_output(['ctest','--test-dir',str(build),'--show-only=json-v1','-R',regex],env=env,text=True)
    import json
    tests=json.loads(inventory)['tests']
    if len(tests)!=count:raise ValueError('required registration count mismatch: '+str([r['name'] for r in tests]))
    write_once(root/(label+'-inventory.json'),json.loads(inventory))
    allocation=root/(label+'-allocation.json')
    subprocess.run(['python3','-B',str(scripts/'wj4_resources.py'),'--reserve','10','--output',str(allocation)],env=env,check=True)
    ledgers={p:Ledger(Path(env['WJ4_PACKAGE_ROOT']),p) for p in ('EP','Budget')}
    if a.suite!='composite':
        for l in ledgers.values():l.reserve('ctest/'+label,dict(regex=regex,inventory_sha256=file_hash(root/(label+'-inventory.json'))),reserve,1200)
    start=time.time()
    rc=subprocess.run(['python3','-B',str(scripts/'wj3_receipts.py'),'--root',str(root),'--label','ctest-'+label,'--',
        'ctest','--test-dir',str(build),'-j','1','--output-on-failure','-R',regex],env=env).returncode
    elapsed=time.time()-start
    index=[]
    if a.suite in ('native','readonly'):
        scopes=[build/('wj2-native-'+p) for p in ('EP','Budget')]+[build/'wj2-legacy-transaction'] if a.suite=='native' else [build/'wj3-smoke']
        for scope in scopes:
            for file in scope.rglob('*.evolution.json'):
                d=read(file)
                if 'pimple_attempts' in d and file.stat().st_mtime>=start:
                    policy=('Budget' if 'Budget' in d['wall_policy'] else 'EP') if d.get('turbulence_model')=='WJHighRe' else 'SharedLegacy'
                    index.append(dict(path=str(file),sha256=file_hash(file),policy=policy,attempts=len(d['pimple_attempts']),failure=d.get('failure')))
    if a.suite!='composite':
        for policy,l in ledgers.items():l.finish('ctest/'+label,sum(r['attempts'] for r in index if r['policy'] in (policy,'SharedLegacy')),elapsed,str(root/('ctest-'+label+'.log')))
    write_once(root/(label+'-accounting.json'),dict(exit_code=rc,seconds=elapsed,attempt_index=index,
        wall_attribution='full shared/control CTest time and shared legacy scheduled attempts conservatively charged to each policy; composite sessions account their own policy wall time',
        ledgers={p:l.totals() for p,l in ledgers.items()}))
    log=build/'Testing/Temporary/LastTest.log'
    if log.exists():(root/(label+'-LastTest.log')).write_bytes(log.read_bytes())
    raise SystemExit(rc)

if __name__=='__main__':main()
