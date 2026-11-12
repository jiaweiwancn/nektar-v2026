#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""CT4 input-provider negatives, exact histories and unchanged-route controls."""
import argparse,copy,json,os,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from ct_accepted4_run import make,construction
from ct_accepted4_reference import Column,parameters,write_json
from ct_accepted_open import require,run,probe
from compare_cls_epawf_pimple import same_history
from run_cls_epawf_picard import set_info,set_parameter,write
from rans_request_manifest import RequestManifest,runtime_identity,sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','ct2','output'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--runtime',type=Path);p.add_argument('--initial-root',type=Path)
    a=p.parse_args();r=a.root;rt=a.runtime or r/'build-v1/runtime';ir=a.initial_root or r/'initial-v1';os.environ['LD_LIBRARY_PATH']=str(rt/'lib')
    manifest=RequestManifest(a.output,dict(schema='CT4ProviderControlsV1',runtime=runtime_identity(rt/'IncNavierStokesSolver'),
        references_sha256=sha(r/'reference-v1/results.json'),source_sha256=sha(Path(__file__))))
    result=dict(implementation='fail',continuous_and_consumed_inputs={},negative_controls={},history={},executions=[])
    for model in ('A93','A95','E94E95'):
        for h in (5.5,6.,6.5):
            key=model+'-'+str(h)+'-P3';file=ir/key/'results.json'
            q=json.loads(file.read_text());rows=q['members']
            require(q['identical_pair_initial_coefficients'],'wrong_dt_pair_initial_state')
            for dt,row in rows.items():
                require(row['continuous_oracle']['status']=='pass' and row['constitutive']['status']=='pass','continuous_oracle_or_full_law_mismatch')
                z=json.loads((ir/key/dt/'initial.json').read_text())
                rec=row['initialization_receipt'];prec=z['initialization_receipt']
                require(rec['explicit_boundary_sha256']==rec['implicit_boundary_sha256']==prec['explicit_boundary_sha256']==prec['implicit_boundary_sha256'],
                    'driver_probe_explicit_implicit_live_boundary_mismatch')
                require(rec['effective_initial_sha256']==prec['effective_initial_sha256'],'driver_probe_effective_initial_mismatch')
            result['continuous_and_consumed_inputs'][key]=dict(status='pass',receipt_sha256=sha(file),
                continuous_oracle={dt:row['continuous_oracle'] for dt,row in rows.items()},native_slots={dt:row['slots'] for dt,row in rows.items()},
                projected_numerical_admission=q['status'])
            # These are reference-verification negatives, not rejection rules
            # for arbitrary nonequilibrium states in the governing equations.
            c=q['case'];pi=c['reference_column']['outer_equilibrium']['Pi'];signal=rows['0.002']['continuous_oracle']['nonconstant_inner_k_signal']
            require(signal>1e-3,'constant_inner_K_substitution_not_discriminated')
            result['negative_controls'][key+'-constant-inner-k']=dict(status='pass',expected='reference comparison fails',scaled_defect=signal)
            if model!='E94E95':
                require(abs(pi)/c['ustar']**2>1e-3,'missing_pressure_offset_not_discriminated')
                result['negative_controls'][key+'-missing-pressure-offset']=dict(status='pass',expected='full normal traction reference comparison fails',
                    defect=-pi,normalized=abs(pi)/c['ustar']**2)
    initial=ir/'A95-6.0-P3/0.002';c=json.loads((initial/'construction.json').read_text());zero=json.loads((initial/'case.evolution.json').read_text())
    for name,change,reason in [('wrong_model_equilibrium',lambda root:set_parameter(root,'RANSOpenEquilibriumQ',3.4012654174),'wrong model equilibrium data'),
                              ('unknown_provider',lambda root:set_info(root,'RANSOpenReference','NeutralCompositeV999'),'unknown open reference provider version')]:
        folder=a.output/name;make(folder,rt,c,.002,0,template=initial/'case.xml');root=ET.parse(folder/'case.xml').getroot();change(root);write(root,folder/'case.xml')
        packet,e=probe(rt,folder,dict(operation='functions',grid=[],grid_sides=[]),'negative')
        require(e['exit_code']!=0 and reason in packet.get('failure',''),'wrong_negative_reason:'+name)
        result['negative_controls'][name]=dict(status='pass',expected_reason=reason,observed=packet['failure'],execution=e)
    # Same method, real native histories. No signature rewriting or changed dt.
    reports={}
    for name,steps,parent in [('full',2,Path(zero['checkpoint'])),('first',1,Path(zero['checkpoint']))]:
        folder=a.output/name;make(folder,rt,c,.002,steps,parent,initial/'case.xml');e=run([rt/'IncNavierStokesSolver','case.xml'],folder,'solver.log')
        n=json.loads((folder/'case.evolution.json').read_text());require(e['exit_code']==0 and not n['failure'],'unchanged_method_history_control_failed');reports[name]=n;result['executions'].append(e)
    folder=a.output/'restart';make(folder,rt,c,.002,1,Path(reports['first']['checkpoint']),initial/'case.xml');e=run([rt/'IncNavierStokesSolver','case.xml'],folder,'solver.log')
    n=json.loads((folder/'case.evolution.json').read_text());require(e['exit_code']==0 and same_history(reports['full'],n),'neutral_provider_restart_not_equivalent')
    result['history']['same_method_restart']=dict(status='pass',state_sha256=exact_vector_hash(n['state']),execution=e)
    folder=a.output/'rollback';make(folder,rt,c,.002,1,Path(zero['checkpoint']),initial/'case.xml');root=ET.parse(folder/'case.xml').getroot();set_parameter(root,'RANSWallPimpleRejectStep',1);write(root,folder/'case.xml')
    e=run([rt/'IncNavierStokesSolver','case.xml'],folder,'solver.log');n=json.loads((folder/'case.evolution.json').read_text())
    require(e['exit_code']!=0 and same_history(zero,n) and n.get('rejected_state'),'neutral_provider_rollback_failed')
    result['history']['rejected_rollback']=dict(status='pass',execution=e,last_state_sha256=exact_vector_hash(n['state']),rejected_state_sha256=exact_vector_hash(n['rejected_state']))
    folder=a.output/'changed_physical_identity';make(folder,rt,c,.002,0,Path(zero['checkpoint']),initial/'case.xml');root=ET.parse(folder/'case.xml').getroot()
    us=c['ustar']*(1+1e-6);q=Column(c['model'],parameters(c['model'])['q']*us*us,c['H'],c['roughness_input'],c['nu'])
    set_parameter(root,'RANSWallFrictionVelocity',us);set_parameter(root,'RANSOpenMeanIntercept',q.t0*q.J0-us/q.p['kappa']*np.log(c['H']));write(root,folder/'case.xml')
    e=run([rt/'IncNavierStokesSolver','case.xml'],folder,'solver.log');text=(folder/'solver.log').read_text()
    require(e['exit_code']!=0 and 'checkpoint signature mismatch' in text,'changed_physical_reference_history_not_rejected')
    result['history']['changed_physical_identity']=dict(status='pass',expected_reason='checkpoint signature mismatch',execution=e)
    result['implementation']='pass';write_json(a.output/'results.json',result);manifest.seal(a.output,['results.json']);print('CT4 provider controls pass');return 0
if __name__=='__main__':raise SystemExit(main())
