#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact, honestly labelled numerical outcomes from immutable native reports."""
import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
from rans_qualification import dependencies,coverage,audited_extrema,native_check_coverage
dependencies("numpy","scipy")
import numpy as np
from validate_a93_a95_e94e95 import oracle,constants
from a93_a95_e94e95_note_wall_oracle import Kernel


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def maximum(rows,key,default=0.):
    values=[r[key] for r in rows if key in r]
    return max(values,default=default)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--matrix',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--qualification',type=Path);a=p.parse_args()
    m=json.loads(a.matrix.read_text());result=dict(schema=1,tested_source=m['source'],ranks=m['ranks'],threads=m['threads'],
        scope='bounded serial flat2D native one-sweep; native polynomial reference is not a roundoff equilibrium',
        matrix_sha256=digest(a.matrix),cases=[])
    for c in m['cases']:
        try:
            path=Path(c['path'])/'case.evolution.json';r=json.loads(path.read_text());audits=r.get('state_audits',[]);initial=audits[0] if audits else {};last=audits[-1] if audits else {};attempts=r.get('pimple_attempts',[])
            d=last.get('accepted_awf_scalar_domain',[])
            extrema=audited_extrema(r)
            completed=[x for x in attempts if x.get('scheduled_complete')]
            requested=c.get('options',{}).get('steps',16)
            qm,kap,sigma=constants(c['model'],c['kappa_policy'])
            wall_kernel=Kernel(qm*.16,1.5e-5,2.,.1,.2295 if c['wall']=='EP' else kap/math.sqrt(qm),2.55 if c['wall']=='EP' else kap*qm**1.5)
            reference_gradient=np.zeros((3,3));reference_gradient[0,2]=.4/(kap*2.)
            expected_pressure_jump=float(oracle(c['model'],qm*.4**2,.4**3/(kap*2.),reference_gradient)['Q'][2,2])
            actual_pressure_jump=initial['interface'][1]/initial['interface'][0]
            # Actual one-sided native tensors against the independent dimensional
            # oracle, including the spanwise normal stress in these 2D runs.
            witness_errors={};raw_divergence_witness={}
            for side in ('inner','outer'):
                tensor=np.asarray(last['native_'+side+'_covariance_tensor_row_major']).reshape(3,3)
                gradient=np.asarray(last['native_'+side+'_covariance_witness_gradient_row_major']).reshape(3,3)
                raw_divergence_witness[side]=float(np.trace(gradient))
                point=last['accepted_awf_'+side+'_minimum_point'];k=point[3]
                expected=(2*k*np.eye(3)/3-point[5]*(gradient+gradient.T) if side=='inner' else oracle(c['model'],k,point[4],gradient)['R'])
                witness_errors[side]=float(np.max(np.abs(tensor-expected))/max(k,1e-30))
            checks=dict(unanchored=r.get('eta')==0 and r.get('raw_reference_norms')==r.get('corrected_reference_norms'),
                exactly_one_sweep=len(completed)==requested and len(attempts)==requested and all(s.get('outer_sweeps')==1 and all(s.get(k)==1 for k in ('momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves')) for s in attempts if s.get('scheduled_complete')),
                no_bound_correction=all(s.get('bound_calls',0)==0 for s in attempts),
                same_state_sources=maximum(attempts,'k_same_state_source_storage_defect')<1e-9 and maximum(attempts,'epsilon_same_state_source_storage_defect')<1e-9,
                compatible_stage_constraints=maximum(attempts,'imposed_trace_defect')<1e-12 and all(max(s.get('corrected_continuity_gauge',[0]))<1e-10 for s in attempts),
                positive_covariance=all(x['value']/(qm*.16)>=-1e-10 for key,x in extrema['over_all_audited_states_and_stages'].items() if key in ('inner_covariance','outer_covariance')) and all(key in extrema['over_all_audited_states_and_stages'] for key in ('inner_covariance','outer_covariance')),
                independent_full_native_tensor_witness=max(witness_errors.values())<1e-11,
                independent_reference_pressure_offset=abs(actual_pressure_jump-expected_pressure_jump)<1e-11)
            result['cases'].append(dict(model=c['model'],kappa_policy=c['kappa_policy'],wall=c['wall'],dt=r.get('dt'),accepted_steps=r.get('steps',0),time=r.get('time'),
                independent_analytical_reference=dict(q=qm,kappa=kap,sigma_epsilon=sigma,friction_velocity=.4,physical_nu=1.5e-5,H=2.,roughness_input=.1,
                    K=qm*.16,epsilon_H=.4**3/(kap*2.),velocity_H=.16*wall_kernel.J0(2.),aerodynamic_z0=2.*math.exp(-kap*.4*wall_kernel.J0(2.))),
                bounded_integration_pass=c['status']=='pass' and all(checks.values()),checks=checks,failure=r.get('failure'),
                audited_extrema=extrema,execution_status=c.get('execution_status',c['status']),accuracy_status='bounded_checks_only',convergence_status='not_established_by_endpoint',
                initial_native_spatial_residual=r.get('raw_reference_norms'),
                initial_k_inner_outer_shared_interface=initial.get('native_k_inner_outer_shared_interface_dual'),
                final_stage_residual=last.get('coupled_norms'),final_k_inner_outer_shared_interface=last.get('native_k_inner_outer_shared_interface_dual'),
                drift_at_1_2_8_16=[s for s in r.get('step_audits',[]) if int(s[1]) in (1,2,8,16)],
                covariance_minima=last.get('accepted_awf_covariance'),minimum_inner_k=d[2] if d else None,minimum_outer_k=d[3] if d else None,minimum_outer_epsilon=d[4] if d else None,
                native_full_tensor_independent_relative_errors=witness_errors,
                native_raw_divergence_at_covariance_witness=raw_divergence_witness,
                native_inner_full_covariance=last.get('native_inner_covariance_tensor_row_major'),native_outer_full_covariance=last.get('native_outer_covariance_tensor_row_major'),
                reference_Pi_inner_minus_outer=dict(native=actual_pressure_jump,independent=expected_pressure_jump),
                bound_correction_budgets=dict(calls=sum(s.get('bound_calls',0) for s in attempts),
                    signed_integral=0. if checks['no_bound_correction'] else None,absolute_integral=0. if checks['no_bound_correction'] else None),
                minimum_budget_ratio=d[7] if d else None,maximum_budget_relative_defect=d[8] if d else None,
                raw_divergence_relative_reference_shear=last.get('outer_constitutive_range',[None]*8)[7],
                pointwise_interface_diagnostic=last.get('interface'),common_velocity_k_dependent_epsilon_trace=last.get('common_velocity_k_dependent_epsilon_trace_defects'),
                max_refreshed_epsilon_trace_defect=maximum(attempts,'refreshed_trace_defect'),
                report_sha256=digest(path),xml_sha256=c['xml_sha256'],report=str(path)))
        except Exception as error:
            result['cases'].append(dict(model=c.get('model'),kappa_policy=c.get('kappa_policy'),wall=c.get('wall'),bounded_integration_pass=False,failure=repr(error)))
    expected={(m,p,w) for m in ('A93','A95','E94E95') for p in ('Natural','Prescribed') for w in ('EP','Accepted')}
    result['coverage']=coverage(result['cases'],expected)
    result['bounded_matrix_pass']=result['coverage']['status']=='pass' and all(c['bounded_integration_pass'] for c in result['cases'])
    if a.qualification:
        q=json.loads(a.qualification.read_text());result['continuation_controls']=q['checks'];result['continuation_status']=q['status']
        wanted=q.get('expected_check_names')
        result['continuation_coverage']=native_check_coverage(q['checks'],wanted)
        if not wanted or result['continuation_coverage']['status']!='pass':result['continuation_status']='fail'
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(bounded_matrix_pass=result['bounded_matrix_pass'],configurations=len(result['cases']))))
    return 0 if result['bounded_matrix_pass'] and result.get('continuation_status','pass')=='pass' else 1
if __name__=='__main__':sys.exit(main())
