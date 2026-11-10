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
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--suffix',default='');a=p.parse_args();r=a.root;out=r/('frozen-initial-history-controls'+a.suffix);out.mkdir(exist_ok=False)
    os.environ.update(THREADS,LD_LIBRARY_PATH=str(r/'runtime/lib'),TMPDIR=str(r/'tmp'))
    sampler=r/'runtime-frozen-history/RANSPreservationProbe';grid=load(r/'preservation-spatial/evaluation-grid.json');weights=grid['weights'];result=dict(checks=[],status='running',execution_scope='fresh initialization artifact; no evolved checkpoint conversion')
    for model in ('E94E95','A95','A93'):
        name=model+'-Natural-Accepted-dt0'
        try:
            ref=Reference(model,'Natural','Accepted');lowdir=r/('preservation-temporal-Q0-'+model)/name;lo=load(lowdir/'state-packet.json');lor=load(lowdir/'case.evolution.json')
            folder=out/name;case_xml(folder,model,'Natural','Accepted',steps=16,descriptor_probe=r/'runtime/RANSNewModelsProbe',radial_knots=ref.knots(),degree=3,quadrature=Q[2],cache='ExactStepV1',checkpoint_interval=16)
            request=replay_request(lo,[],state=lo['reference'],time=0.);request['write_initial_history']=str(folder/'initial-native-history')
            initial=probe(sampler,folder,request);receipt=header(Path(initial['initial_history']),True)
            assert receipt['steps']==0 and receipt['time']==0. and receipt['order']==1
            assert exact_vector_hash(receipt['coupled_state'])==exact_vector_hash(lor['accepted_states'][0]),'export changed initial coefficients'
            differential=initial['observation']['offsets'][4]
            assert receipt['history'][0]==lo['reference'][:differential] and all(v==0. for v in receipt['history'][1]),'inconsistent fresh native history'
            root=ET.parse(folder/'case.xml').getroot();set_info(root,'RANSWallRestart',initial['initial_history']);write(root,folder/'case.xml')
            with (folder/'solver.log').open('w') as f:child=subprocess.run([str(r/'runtime/IncNavierStokesSolver'),'case.xml'],cwd=folder,stdout=f,stderr=subprocess.STDOUT,timeout=900)
            hi=load(folder/'case.evolution.json');assert child.returncode==0 and hi['steps']==16 and not hi['failure'] and abs(hi['time']-.0016)<1e-14, 'target native execution: '+str(hi.get('failure'))
            assert exact_vector_hash(hi['accepted_states'][0])==exact_vector_hash(lor['accepted_states'][0]),'native restart initialization changed coefficients'
            final=probe(sampler,folder,replay_request(lo,grid['grid'],state=hi['state'],time=hi['time']))
            comparison=compare_samples(dict(final['observation'],grid_weights=weights),dict(lo['observation'],grid_weights=weights),weights)
            result['checks'].append(dict(name='Q2-Q0-'+model+'-Accepted',status='pass' if max(v['weighted_RMS'] for v in comparison['fields'].values())<1e-10 else 'fail',
                comparison=comparison,initial_state_sha256=exact_vector_hash(receipt['coupled_state']),fresh_history_sha256=receipt['sha256'],native_initial_state_identical=True,
                source_initializer_quadrature=Q[0],target_equation_quadrature=Q[2],target_signature_verified_by_native_restart=True,xml_sha256=sha(folder/'case.xml')))
            # Negative replay/signature controls use an explicit preserved input.
            if model=='E94E95':
                for kind in ('descriptor','identity','nonzero_time','different_reference'):
                    guard=copy.deepcopy(request)
                    if kind in ('descriptor','identity'):guard[kind]='changed-'+guard[kind]
                    elif kind=='nonzero_time':guard['time']=.0016
                    else:
                        guard['state']=list(guard['state']) # break an intentional state/reference alias before perturbing
                        guard['state'][-1]+=1e-3
                    guard['write_initial_history']=str(folder/('must-not-exist-'+kind));write_json(folder/('guard-'+kind+'.json'),guard)
                    with (folder/('guard-'+kind+'.log')).open('w') as f:c=subprocess.run([str(sampler),'case.xml','guard-'+kind+'.json','guard-'+kind+'-result.json'],cwd=folder,stdout=f,stderr=subprocess.STDOUT)
                    failure=load(folder/('guard-'+kind+'-result.json')).get('failure','')
                    reason='descriptor mismatch' if kind=='descriptor' else 'identity mismatch' if kind=='identity' else 'initial history export requires'
                    result['checks'].append(dict(name='guard-'+kind,status='pass' if c.returncode!=0 and reason in failure and not Path(guard['write_initial_history']).exists() else 'fail',intended_rejection=reason,observed_failure=failure,exit_code=c.returncode))
        except Exception as error:result['checks'].append(dict(name='Q2-Q0-'+model+'-Accepted',status='fail',failure=repr(error)))
        write_json(out/'results.json',result)
    expected={'Q2-Q0-'+m+'-Accepted' for m in ('E94E95','A95','A93')}|{'guard-'+s for s in ('descriptor','identity','nonzero_time','different_reference')}
    result['coverage']=coverage(result['checks'],expected,key=lambda x:x['name']);result['status']='pass' if result['coverage']['status']=='pass' and all(c['status']=='pass' for c in result['checks']) else 'fail';write_json(out/'results.json',result)
    original=load(r/'uncertainty/results.json');combined=dict(original);combined['checks']=[c for c in original['checks'] if not c['name'].endswith('-Accepted') or c['name'].startswith('evaluation-')]
    combined['checks'].extend(c for c in result['checks'] if c['name'].startswith('Q2-Q0-'))
    combined['strict_initialization_controls']=result;combined['initial_quadrature_failure_receipt']=str(r/'uncertainty/results.json')
    wanted={'Q2-Q0-'+m+'-'+w for m in ('E94E95','A95','A93') for w in ('EP','Accepted')}|{'evaluation-GL10-GL14-'+w for w in ('EP','Accepted')}
    combined['coverage']=coverage(combined['checks'],wanted,key=lambda x:x['name']);combined['status']='pass' if result['status']=='pass' and combined['coverage']['status']=='pass' and all(x['status']=='pass' for x in combined['checks']) else 'fail'
    write_json(r/'uncertainty-final.json',combined);print(combined['status']);return 0 if combined['status']=='pass' else 1
if __name__=='__main__':sys.exit(main())
