#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Native quadrature trajectories from byte-identical initial coefficients.

A fresh time-zero IMEX history is created by the native API with the TARGET
session signature. No evolved history/signature is patched or relabelled.
"""
import argparse
import copy
import os
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET
from rans_qualification import write_json,exact_vector_hash,compare_samples,coverage,THREADS
from qualify_rans_preservation import load,probe,replay_request,Q
from rans_preservation_oracle import Reference
from run_a93_a95_e94e95_matrix import case_xml
from run_cls_epawf_picard import set_info,write,sha
from continue_cls_epawf_pimple import header

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path);p.add_argument('--steps',type=int,default=2);a=p.parse_args();r=a.root.resolve();out=a.output or r/'frozen-boundary-regression';out.mkdir(parents=True,exist_ok=False)
    os.environ.update(THREADS,LD_LIBRARY_PATH=str(r/'runtime/lib'),TMPDIR=str(r/'tmp'))
    sampler=r/'runtime/RANSPreservationProbe';result=dict(schema='FrozenEffectiveInputsV1',checks=[],status='running')
    from rans_preservation_oracle import physical_grid
    for model in ('E94E95','A95','A93'):
        try:
            ref=Reference(model,'Natural','Accepted');grid,weights=physical_grid(ref.knots(),6,3);packets=[];reports=[]
            for qi in (0,2):
                folder=out/(model+'-Q'+str(qi));case_xml(folder,model,'Natural','Accepted',steps=a.steps,descriptor_probe=r/'runtime/RANSNewModelsProbe',radial_knots=ref.knots(),degree=3,quadrature=Q[qi],cache='ExactStepV1',checkpoint_interval=a.steps,boundary_policy='FrozenReferenceV1')
                if qi==0:
                    seed=probe(sampler,folder,dict(operation='neutral',grid=grid,export_basis=(model=='E94E95')))
                request=replay_request(seed,grid,state=seed['reference'],time=0.);request['write_initial_history']=str(folder/'initial-native-history')
                exported=probe(sampler,folder,request);receipt=header(Path(exported['initial_history']),True)
                assert receipt['steps']==0 and receipt['time']==0. and receipt['order']==1
                assert exact_vector_hash(receipt['coupled_state'])==exact_vector_hash(seed['reference'])
                root=ET.parse(folder/'case.xml').getroot();set_info(root,'RANSWallRestart',exported['initial_history']);write(root,folder/'case.xml')
                with (folder/'solver.log').open('w') as f:c=subprocess.run([str(r/'runtime/IncNavierStokesSolver'),'case.xml'],cwd=folder,stdout=f,stderr=subprocess.STDOUT,timeout=600)
                report=load(folder/'case.evolution.json');assert c.returncode==0 and not report['failure'] and report['steps']==a.steps,report.get('failure')
                ir=report['initialization_receipt'];assert exact_vector_hash(report['accepted_states'][0])==exact_vector_hash(seed['reference'])
                assert ir['raw_initial_sha256']==ir['effective_initial_sha256']==ir['fixed_reference_sha256']
                assert ir['explicit_boundary_sha256']==ir['implicit_boundary_sha256']
                assert ir['raw_effective_max_difference']==0.
                request=replay_request(seed,grid,state=report['state'],time=report['time']);request['fixed_boundary']=ir['live_fixed_boundary']
                packets.append(probe(sampler,folder,request)['observation']);reports.append(report)
                if model=='E94E95' and qi==2:
                    for kind in ('descriptor','identity','nonzero_time','different_reference','boundary'):
                        guard=replay_request(seed,[],state=seed['reference']);guard['write_initial_history']=str(folder/('must-not-exist-'+kind))
                        if kind in ('descriptor','identity'):guard[kind]='changed-'+guard[kind]
                        elif kind=='nonzero_time':guard['time']=.0016
                        elif kind=='boundary':guard['fixed_boundary']=copy.deepcopy(guard['fixed_boundary']);guard['fixed_boundary'][0][1]+=1e-3
                        else:guard['state']=list(guard['state']);guard['state'][-1]+=1e-3
                        write_json(folder/('guard-'+kind+'.json'),guard)
                        with (folder/('guard-'+kind+'.log')).open('w') as f:gc=subprocess.run([str(sampler),'case.xml','guard-'+kind+'.json','guard-'+kind+'-result.json'],cwd=folder,stdout=f,stderr=subprocess.STDOUT)
                        failure=load(folder/('guard-'+kind+'-result.json')).get('failure','');reason='descriptor mismatch' if kind=='descriptor' else 'identity mismatch' if kind=='identity' else 'external boundary identity mismatch' if kind=='boundary' else 'initial history export requires'
                        result['checks'].append(dict(name='guard-'+kind,status='pass' if gc.returncode!=0 and reason in failure and not Path(guard['write_initial_history']).exists() else 'fail',observed_failure=failure,intended_rejection=reason,exit_code=gc.returncode))
            keys=('raw_initial_sha256','effective_initial_sha256','fixed_reference_sha256','explicit_boundary_sha256','implicit_boundary_sha256')
            assert all(reports[0]['initialization_receipt'][k]==reports[1]['initialization_receipt'][k] for k in keys),'effective input mismatch'
            comparison=compare_samples(packets[0],packets[1],weights)
            assert max(v['weighted_RMS'] for v in comparison['fields'].values())<1e-10
            result['checks'].append(dict(name='Q2-Q0-'+model+'-Accepted',status='pass',comparison=comparison,initialization_receipts=[x['initialization_receipt'] for x in reports],steps=a.steps))
        except Exception as e:result['checks'].append(dict(name='Q2-Q0-'+model+'-Accepted',status='fail',failure=repr(e)))
        write_json(out/'results.json',result)
    expected={'Q2-Q0-'+m+'-Accepted' for m in ('E94E95','A95','A93')}|{'guard-'+k for k in ('descriptor','identity','nonzero_time','different_reference','boundary')}
    result['coverage']=coverage(result['checks'],expected,key=lambda x:x['name']);result['status']='pass' if result['coverage']['status']=='pass' and all(x['status']=='pass' for x in result['checks']) else 'fail';write_json(out/'results.json',result);print(result['status']);return 0 if result['status']=='pass' else 1
if __name__=='__main__':sys.exit(main())
