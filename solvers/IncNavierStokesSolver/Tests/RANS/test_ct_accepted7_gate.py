#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
import copy,unittest
from ct_accepted7_gate import top,temporal,resources,history,TOP_CASES,TOP_EXTRA

class CT7Gates(unittest.TestCase):
    def setUp(self):
        row=dict(status='pass',defects={k:0. for k in ['all_row','boundary_advection','boundary_traction','boundary_diffusion',
            *('frozen_'+str(i) for i in (0,2,3)),*('residual_'+str(i) for i in (0,2,3))]},
            continuous_native_test_maximum=1e-14,additional_pressure_null_modes=0,pressure_minimum_singular_ratio=.01,constant_pressure_action=.1,
            zero_epsilon_flux_negative_signal=.001,reversed_shear_negative_signal=.1,
            initialization_receipt=dict(open_explicit_implicit_agreement=True,explicit_boundary_sha256='same',implicit_boundary_sha256='same'))
        checks={k:copy.deepcopy(row) for k in TOP_CASES}
        checks.update(pressure_offset=dict(status='pass',coefficient_defect=0.,free_residual_defect=1e-15),
            state_independent_epsilon_flux=dict(status='pass',datum=-.001,wrong_evolving_diffusivity_signal=.0001),
            displaced_original_wall=dict(status='pass',reason='unqualified displacement'),unsupported_cache=dict(status='pass',reason='caching unqualified'))
        self.r=dict(status='pass',checks=checks)
    def test_complete_top(self):self.assertEqual(top(self.r)['status'],'pass')
    def test_missing_model(self):
        del self.r['checks']['A93-p3'];self.assertIn('missing_or_foreign_top_control',top(self.r)['reasons'])
    def test_pass_label_operator_failure(self):
        self.r['checks']['A95-p4']['defects']['all_row']=3e-10;self.assertIn('A95-p4:operator_accuracy',top(self.r)['reasons'])
    def test_missing_mode_primitive(self):
        del self.r['checks']['A95-p3']['defects']['frozen_3'];self.assertIn('A95-p3:missing_operator_primitive',top(self.r)['reasons'])
    def test_extra_pressure_null_mode(self):
        self.r['checks']['E94E95-p3']['additional_pressure_null_modes']=1;self.assertIn('E94E95-p3:pressure_rank',top(self.r)['reasons'])
    def test_negative_must_be_discriminating(self):
        self.r['checks']['A95-p3']['zero_epsilon_flux_negative_signal']=0.;self.assertNotEqual(top(self.r)['status'],'pass')
    def test_boundary_identity(self):
        self.r['checks']['A95-p3']['initialization_receipt']['implicit_boundary_sha256']='different';self.assertIn('A95-p3:consumed_boundary_identity',top(self.r)['reasons'])
    def test_wrong_expected_rejection(self):
        self.r['checks']['displaced_original_wall']['reason']='some unrelated failure';self.assertIn('wrong_expected_rejection:displaced_original_wall',top(self.r)['reasons'])
    def test_temporal_primitive_overrides_label(self):
        r=dict(temporal_status='pass',orders={k:dict(D002_001=.0002,D001_0005=.0001) for k in ('U','w','k','outer_epsilon','Pi')})
        self.assertEqual(temporal(r)['status'],'pass');r['orders']['Pi']['D001_0005']=.0002;self.assertIn('Pi',temporal(r)['reasons'])
    def test_closed_resources_required(self):
        self.assertEqual(resources({})['status'],'fail')
        r=dict(maximum_project=16,maximum_aggregate=72,project_status='pass',aggregate_status='pass',control_helpers_counted=2)
        self.assertEqual(resources({'closed':r})['status'],'pass');r['maximum_project']=17;self.assertEqual(resources({'closed':r})['status'],'fail')
    def test_foreign_native_smoke_model(self):
        receipt=dict(checks={k:dict(status='pass') for k in ('restart','rollback','old_boundary_history','pressure_level_evolution')},models={});primitives={}
        for model in ('A93','A95','E94E95'):
            desc=model+':OriginalHighRe-v1';binding=dict(model=model,resolved_descriptor=desc,native_report_sha256=model+'-report',state_sha256=model+'-state',configuration_sha256=model+'-xml')
            receipt['models'][model]=dict(completion='pass',time=2e-6,model=model,resolved_descriptor=desc,native_signature='patch:resolved='+desc,
                report_sha256=binding['native_report_sha256'],state_sha256=binding['state_sha256'],configuration_sha256=binding['configuration_sha256'],
                mass=dict(outward_net=0.,scale=1.),captured_k_matrix=dict(maximum_defect=0.,scale=1.))
            eq={k:dict(physical_scale=1.,storage_backward_difference=0.,stage_source=0.,stage_physical_boundary_flux=0.,frozen_stage_physical_defect=0.) for k in ('inner_k','outer_k','outer_epsilon','whole_k')}
            rec={k:dict(scaled_defect=0.,normalization=1.,maximum_defect=0.) for k in ('k','outer_epsilon')}
            primitives[model]=dict(budget_binding=binding,budgets={k:dict(binding=binding,equations=eq,reconstruction=rec) for k in ('native','q48f16','q64f16','q64f32')})
        self.assertEqual(history(receipt,primitives)['status'],'pass')
        receipt['models']['A93']['model']='A95'
        self.assertIn('A93:foreign_native_model_configuration',history(receipt,primitives)['reasons'])
        receipt['models']['A93']['model']='A93';primitives['E94E95']['budget_binding']['model']='A95'
        self.assertIn('E94E95:budget_native_identity_mismatch',history(receipt,primitives)['reasons'])

if __name__=='__main__':unittest.main()
