#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only independent first-rejection witnesses, including effective traces."""
import argparse,json,os,shutil
from pathlib import Path
from rans_qualification import THREADS,write_json,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from ct_accepted_open import require,probe,row_value,run
from ct_accepted_pilot import replay
from ct_accepted3_metrics import point
from validate_a93_a95_e94e95 import oracle
from ct_accepted2_storage import artifact
from rans_request_manifest import RequestManifest,sha,runtime_identity
from ct_accepted3_run import make_case
from compare_cls_epawf_pimple import same_history


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('root','ct2','output'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--reproduce',action='store_true',help='one unchanged-method native attempt from each retained last accepted history')
    a=p.parse_args();rt=a.ct2/'runtime';os.environ['LD_LIBRARY_PATH']=str(rt/'lib')
    files=[]
    for f in a.root.glob('**/case.evolution.json'):
        if not any(case in f.parts for case in ('A95-XZminus-P3','A95-XZminus-P4','A93-Zminus-P3','E94E95-Zminus-P3')):continue
        r=json.loads(f.read_text())
        if r['failure']:files.append((f,r))
    m=RequestManifest(a.output,dict(schema='CT3IndependentFailureWitnessesV1',runtime=runtime_identity(rt/'RANSPreservationProbe'),
        reports={str(f):sha(f) for f,_ in files},source_sha256=sha(Path(__file__)),reproduce=a.reproduce))
    result=dict(diagnostic_execution='pass',failures=[])
    for f,r in files:
        if not r.get('physical_range_events'):continue
        c=json.loads((f.parent/'construction.json').read_text());model=c['model'];dt=str(r['dt'])
        case=next(x for x in f.parts if x in ('A95-XZminus-P3','A95-XZminus-P4','A93-Zminus-P3','E94E95-Zminus-P3'))
        initial=next(x for x in (a.root/case).rglob('initial.json') if x.parent.name==dt)
        ini=json.loads(initial.read_text());event=r['physical_range_events'][0];audit=event['audit'];cov=audit['accepted_awf_covariance']
        region=1 if cov[1]<cov[0] else 0;side='outer' if region else 'inner'
        coords=audit['accepted_awf_'+side+'_minimum_point'][:3];native_inputs=audit['accepted_awf_'+side+'_minimum_point']
        provenance=audit['native_'+side+'_minimum_geometry_element_side_block'];G=np.array(audit['native_'+side+'_covariance_witness_gradient_row_major']).reshape(3,3)
        native_R=np.array(audit['native_'+side+'_covariance_tensor_row_major']).reshape(3,3)
        if region:
            k,e=native_inputs[3:5];law=oracle(model,k,e,G);expected=law['R'];P=law['production'];Q=law['Q']
        else:
            k,nt=native_inputs[3],native_inputs[5];e=None;expected=np.eye(3)*2*k/3-nt*(G+G.T);P=-np.sum(expected*G);Q=np.zeros((3,3))
        eig,eigenvectors=np.linalg.eigh(expected);defect=float(np.max(abs(expected-native_R)))
        require(defect<2e-10 and eig[0]<-1e-10*c['K0'],'independent_expected_covariance_exit_not_reproduced')
        require(exact_vector_hash(event['state'])==exact_vector_hash(r['rejected_state']) and event['stopped_before_history_advance'],'rejected_candidate_or_history_event_mismatch')
        folder=a.output/(case+'-'+dt);folder.mkdir();shutil.copy2(f.parent/'case.xml',folder/'case.xml')
        req=replay(ini,r['state'],r['time'],[coords],export_test_basis=True)
        req.update(grid_sides=[int(provenance[2])],grid_test_geometries=[int(provenance[0])])
        packet,execution=probe(rt,folder,req,'witness-basis');require(execution['exit_code']==0,'last_accepted_witness_basis_failed')
        q=packet['test_basis'][0];off=packet['observation']['offsets'];stages={'last_accepted':r['state']}
        attempt=r['pimple_attempts'][-1]
        if 'flow_snapshot' in attempt:
            flow=json.loads(artifact(attempt['flow_snapshot'],attempt['flow_snapshot_sha256']).read_text())
            for name in ('input','momentum_predictor','corrected_flow'):stages[name]=flow[name]
        if 'k_snapshot' in attempt:
            scalar=json.loads(artifact(attempt['k_snapshot'],attempt['k_snapshot_sha256']).read_text())
            stages.update(post_epsilon_pre_k=scalar['pre_k'],post_k_raw=scalar['raw'])
        stages['rejected_candidate']=event['state'];observations={}
        for name,state in stages.items():
            if region:
                v=point(q,np.array(state),off,model)
                if 'R' in v:v['covariance_eigenvalues']=np.linalg.eigvalsh(v['R']).tolist()
            else:
                Gs=np.zeros((3,3))
                for b,axis in enumerate((0,2)):Gs[axis]=row_value(q['velocity_k'],np.array(state),off[b])[1]
                pi,dp=row_value(q['pressure'],np.array(state),off[4])
                v=dict(k=float(row_value(q['velocity_k'],np.array(state),off[2])[0]),raw_gradient=Gs.tolist(),raw_divergence=float(np.trace(Gs)),
                    Pi=float(pi),pressure_gradient=dp.tolist(),scope='inner raw field; authoritative native wall viscosity inputs retained separately')
            observations[name]=v
        row=dict(case=case,dt=r['dt'],last_accepted_time=r['time'],attempted_time=event['time'],failure=r['failure'],reasons=event['reasons'],
            coordinate=coords,geometry=int(provenance[0]),native_element=int(provenance[1]),side=int(provenance[2]),region=side,
            independent_native_audit_inputs=dict(k=k,epsilon=e,raw_gradient=G.tolist(),R=expected.tolist(),Q=Q.tolist(),production=float(P),
                divergence=float(np.trace(G)),isotropic_work=float(-2*k*np.trace(G)/3),covariance_eigenvalues=eig.tolist(),minimum_eigenvector=eigenvectors[:,0].tolist(),tensor_defect=defect),
            stages_at_same_physical_point=observations,actual_initial_state_sha256=exact_vector_hash(ini['observation']['state']),
            rejected_state_sha256=exact_vector_hash(event['state']),last_accepted_state_sha256=exact_vector_hash(r['state']),
            native_report_sha256=sha(f),native_event_audit=audit,execution=execution,
            stopping='native original covariance gate, no clipping or retry; diagnostic intermediate values do not add a rejection rule')
        if region:
            row['native_audit_epsilon_minus_raw_rejected_field']=float(e-observations['rejected_candidate']['epsilon'])
            row['epsilon_scope']='native audit may include the declared refreshed dependent trace; raw frozen-stage field retained separately, never repaired'
        if all(name in observations for name in ('input','momentum_predictor','corrected_flow')):
            gi,gp,gc=[np.array(observations[name]['raw_gradient']) for name in ('input','momentum_predictor','corrected_flow')]
            row['predictor_gradient_increment']=(gp-gi).tolist();row['pressure_mass_correction_gradient_increment']=(gc-gp).tolist()
        if a.reproduce:
            sub=folder/'single-attempt-reproducer';make_case(sub,f.parent/'case.xml',c,c['degree'],r['dt'],1,Path(r['checkpoint']))
            execution=run([rt/'IncNavierStokesSolver','case.xml'],sub,'solver.log');repeated=json.loads((sub/'case.evolution.json').read_text())
            require(execution['exit_code']!=0 and repeated['failure']==r['failure'] and same_history(r,repeated)
                and exact_vector_hash(repeated['rejected_state'])==exact_vector_hash(r['rejected_state']), 'first_rejection_or_complete_rollback_not_reproducible')
            row['minimal_native_reproducer']=dict(status='expected_physical_rejection_reproduced',unchanged_last_accepted_state_history=True,
                rejected_state_bitwise_equal=True,report_sha256=sha(sub/'case.evolution.json'),execution=execution)
        result['failures'].append(row);write_json(a.output/'progress.json',result)
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);print('independent covariance failures',len(result['failures']));return 0
if __name__=='__main__':raise SystemExit(main())
