# SPDX-License-Identifier: MIT
"""A passing enclosing label must not admit incomplete/failed native proofs."""
import copy,subprocess,sys
from pathlib import Path
from ct_accepted9_gate import backend_controls,backend_trajectory,capacity

def capacity_fixture():
    result=dict(schema='CT9FlatCapacityV1',status='pass',cases={})
    for name,count in [('intermediate-p3',8864),('target-p3',33267),('target-p4',60907)]:
        identity='A95:OriginalHighRe-v1:SamePhysicalWallDisplacedAcceptedV1:CT8PhysicalBelow25mmEqualIds012V1:EigenLU-LDLT-IR3-F40M'
        q=dict(identity=identity,stored_slots=count,pressure_C_nnz=0,pressure_iterations=2,masses=[dict(rows=3,nnz=7) for _ in range(6)],scalars=[dict(block=b,rhs_maximum=1.,maximum_true_residual=0.) for b in (2,3)])
        q.update(dict.fromkeys(('flow_rows','flow_nnz','velocity_rows','pressure_rows','predictor_factor_nnz','consistent_mass_factor_nnz','pressure_preconditioner_factor_nnz'),1))
        q.update(dict.fromkeys(('predictor_rhs_maximum','pressure_rhs_maximum','consistent_mass_correction_rhs_maximum'),1.))
        q.update(dict.fromkeys(('predictor_residual','true_pressure_residual','consistent_mass_correction_true_residual','full_native_G_action_defect','true_continuity','splitting_identity'),0.))
        result['cases'][name]=dict(status='static_capacity_pass',zero_time=0.,startup_physical_failure='',native_signature=identity,history_sha256='history',zero_state_sha256='state',native_report_sha256='native',native_count=dict(stored_slots=count),sparse_backend=q,reference_agreement=dict(status='pass'),startup_execution=dict(exit_code=0),frozen_execution=dict(exit_code=0))
    return result

def trajectory_fixture():
    fields={f:dict(weighted_RMS=0.,scale=1.) for f in ('U','w','k','outer_epsilon','Pi')}
    equation=dict(physical_scale=1.,storage_backward_difference=0.,stage_source=0.,stage_physical_boundary_flux=0.,frozen_stage_physical_defect=0.)
    budget=dict(equations={n:copy.deepcopy(equation) for n in ('inner_k','outer_k','outer_epsilon','whole_k')},reconstruction={n:dict(scaled_defect=0.) for n in ('k','outer_epsilon')})
    metric=dict(time=.008,status='pass',drift=fields,mass=dict(outward_net=0.,physical_scale=1.),traction=dict(normalized_wall_traction_error=0.,normalized_interface_traction_RMS=0.))
    member=dict(time=.008,failure='',signature=':resolved=A95:OriginalHighRe-v1:SamePhysicalWallDisplacedAcceptedV1:CT8PhysicalBelow25mmEqualIds012V1',native_report_sha256='sparse',accepted_metrics=[metric],budgets=dict(budgets={n:copy.deepcopy(budget) for n in ('native','q48f16','q64f16','q64f32')},captured_k_matrix=dict(maximum_defect=0.,scale=1.)),initialization_receipt=dict(open_explicit_implicit_agreement=True,explicit_boundary_sha256='b',implicit_boundary_sha256='b'),domain=dict(status='pass'),constitutive=dict(status='pass'),corrections=0,extra_solves=0)
    result=dict(schema='CT9PhysicalPairV1',required_dt=[.002,.001],target=.008,same_actual_initial=True,members={},equivalence={},temporal={f:dict(measured_drift=0.,difference=v) for f,v in fields.items()},traction_temporal={f:dict(spread=0.,measured_drift=0.) for f in ('ground','interface','signed_mean')})
    result['temporal']['Pi']['difference']=dict(weighted_RMS=2e-6,scale=1.)
    for dt in (.002,.001):
        name=str(dt);result['members'][name]=dict(copy.deepcopy(member),dt=dt)
        times=[i*dt for i in range(round(.008/dt)+1)]
        result['equivalence'][name]=dict(status='pass',dense_report='dense',sparse_report='sparse',common_times=[dict(time=t,fields=copy.deepcopy(fields),traction=dict(profiles={f:dict(normalized_weighted_RMS=0.) for f in ('ground','interface')})) for t in times])
    return result

def fixture():
    keys={'raw_residual','constrained_residual','external_constraints','traction_load','diffusive_load','advection_load','predictor','pressure_increment','consistent_mass_correction','scalar_increment_2','scalar_increment_3'}|{'mass_'+str(i) for i in range(6)}|{'block_'+str(i)+'_'+key for i in (0,2,3) for key in ('matrix','residual')}
    q=dict(status='pass',defects=dict.fromkeys(keys,0.),artifacts=[dict(independent=dict(status='pass',defects={'all_row':0.},pressure_minimum_singular_ratio=.001,constant_pressure_action=1.,additional_pressure_null_modes=0)) for _ in range(2)])
    q.update({k:1. for k in ('negative_roughness_pressure_omission_signal','negative_sign_signal','negative_mass_replacement_signal','negative_stale_mask_signal')})
    rows={k:copy.deepcopy(q) for k in ('d0-p3','d0-p4','displaced-p3','displaced-p4')}
    rows['restart_rollback_backend_identity']=dict(status='pass',uninterrupted_state='same',continued_state='same',dense_sparse_state_maximum=0.,rebuilt_dense_CT8_history_bitwise_equal=True,negative=[dict(name=n,status='pass',reason=r) for n,r in [('backend_history','signature'),('stale_fit','signature'),('wrong_model','A95'),('unsupported_cache','caching')]])
    rows['absolute_pressure_level']=dict(status='pass',shift=.125,maximum_coefficient_defect=0.)
    return dict(schema='CT9BackendControlsV2',status='pass',checks=rows)

def main():
    good=fixture();assert backend_controls(good)['status']=='pass'
    if '--failed-required' in sys.argv:
        good['checks']['displaced-p4']['defects']['pressure_increment']=1.;return 0 if backend_controls(good)['status']=='pass' else 1
    for mutate in [lambda q:q['checks']['displaced-p4']['defects'].__setitem__('pressure_increment',1.),lambda q:q['checks']['d0-p3']['defects'].__setitem__('mass_0',float('nan')),lambda q:q['checks'].pop('absolute_pressure_level'),lambda q:q['checks']['displaced-p3']['artifacts'].pop(),lambda q:q['checks']['d0-p4'].__setitem__('negative_mass_replacement_signal',0.),lambda q:q.__setitem__('schema','ForeignModelControls'),lambda q:q['checks']['restart_rollback_backend_identity'].__setitem__('continued_state','foreign')]:
        bad=copy.deepcopy(good);mutate(bad);assert backend_controls(bad)['status']=='fail'
    good=trajectory_fixture();assert backend_trajectory(good)['status']=='pass'
    for mutate in [lambda q:q['equivalence'].pop('0.001'),lambda q:q['equivalence']['0.001']['common_times'].pop(),lambda q:q['equivalence']['0.002']['common_times'][1]['fields']['Pi'].__setitem__('weighted_RMS',1e-6),lambda q:q['equivalence']['0.001']['common_times'][1]['fields']['U'].__setitem__('weighted_RMS',float('nan')),lambda q:q['equivalence']['0.002'].__setitem__('sparse_report','foreign'),lambda q:q['temporal']['Pi']['difference'].__setitem__('weighted_RMS',0.)]:
        bad=copy.deepcopy(good);mutate(bad);assert backend_trajectory(bad)['status']=='fail'
    good=capacity_fixture();assert capacity(good)['status']=='pass'
    for mutate in [lambda q:q['cases'].pop('target-p4'),lambda q:q['cases']['target-p3']['sparse_backend'].__setitem__('true_pressure_residual',1e-6),lambda q:q['cases']['intermediate-p3']['sparse_backend']['scalars'].pop(),lambda q:q['cases']['target-p4']['sparse_backend'].__setitem__('pressure_rhs_maximum',float('nan')),lambda q:q['cases']['target-p3']['frozen_execution'].__setitem__('exit_code',1),lambda q:q['cases']['target-p3'].__setitem__('native_signature','old'),lambda q:q['cases']['target-p4'].__setitem__('history_sha256','')]:
        bad=copy.deepcopy(good);mutate(bad);assert capacity(bad)['status']=='fail'
    code=subprocess.run([sys.executable,'-B',str(Path(__file__).resolve()),'--failed-required']).returncode;assert code==1
    print('CT9 required-scope primitives and process-exit negative controls pass');return 0
if __name__=='__main__':raise SystemExit(main())
