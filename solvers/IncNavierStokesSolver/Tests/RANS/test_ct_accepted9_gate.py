# SPDX-License-Identifier: MIT
"""A passing enclosing label must not admit incomplete/failed native proofs."""
import copy,subprocess,sys
from pathlib import Path
from ct_accepted9_gate import backend_controls

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
    code=subprocess.run([sys.executable,'-B',str(Path(__file__).resolve()),'--failed-required']).returncode;assert code==1
    print('CT9 required-scope primitives and process-exit negative controls pass');return 0
if __name__=='__main__':raise SystemExit(main())
