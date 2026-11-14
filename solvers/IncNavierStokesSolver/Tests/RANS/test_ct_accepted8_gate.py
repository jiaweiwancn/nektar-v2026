# SPDX-License-Identifier: MIT
import copy,unittest
from ct_accepted8_gate import pair

class Checks(unittest.TestCase):
    def complete(self):
        fields={k:dict(scale=1.,weighted_RMS=1e-5) for k in ('U','w','k','outer_epsilon','Pi')}
        eq={k:dict(physical_scale=1.,storage_backward_difference=0.,stage_source=0.,stage_physical_boundary_flux=0.,frozen_stage_physical_defect=0.) for k in ('inner_k','outer_k','outer_epsilon','whole_k')}
        rec={k:dict(scaled_defect=0.,normalization=1.,maximum_defect=0.) for k in ('k','outer_epsilon')}
        budgets=dict(budgets={k:dict(equations=copy.deepcopy(eq),reconstruction=copy.deepcopy(rec)) for k in ('native','q48f16','q64f16','q64f32')},captured_k_matrix=dict(maximum_defect=0.,scale=1.))
        row=dict(time=.008,failure='',signature=':resolved=A95:OriginalHighRe-v1:SamePhysicalWallDisplacedAcceptedV1:CT8PhysicalBelow25mmEqualIds012V1',accepted_metrics=[dict(time=.008,status='pass',drift=fields,mass=dict(outward_net=0.,physical_scale=1.),traction=dict(normalized_wall_traction_error=0.,normalized_interface_traction_RMS=0.))],
            budgets=budgets,domain=dict(status='pass'),constitutive=dict(status='pass'),corrections=0,extra_solves=0,initialization_receipt=dict(open_explicit_implicit_agreement=True,explicit_boundary_sha256='same',implicit_boundary_sha256='same'))
        return dict(target=.008,members={k:copy.deepcopy(row) for k in ('0.002','0.001')},temporal={k:dict(measured_drift=1e-5,difference=dict(scale=1.,weighted_RMS=1e-8)) for k in fields},traction_temporal={k:dict(spread=1e-8,measured_drift=1e-5) for k in ('ground','interface','signed_mean')})
    def test_complete_primitives(self):self.assertEqual(pair(self.complete())['status'],'pass')
    def test_missing_epsilon_budget(self):
        q=self.complete();del q['members']['0.001']['budgets']['budgets']['native']['equations']['outer_epsilon'];self.assertIn('0.001:budget_primitives',pair(q)['reasons'])
    def test_false_captured_matrix_label(self):
        q=self.complete();q['members']['0.002']['budgets']['captured_k_matrix']['maximum_defect']=1e-4;self.assertIn('0.002:captured_frozen_matrix',pair(q)['reasons'])
    def test_boundary_mismatch(self):
        q=self.complete();q['members']['0.002']['initialization_receipt']['implicit_boundary_sha256']='other';self.assertIn('0.002:live_boundary_identity',pair(q)['reasons'])
    def test_empty_pass_label_rejected(self):
        self.assertEqual(pair(dict(status='pass',qualification=dict(status='pass')))['status'],'fail')
    def test_missing_companion_reason(self):
        r=pair(dict(target=.008,members={},temporal={}))
        self.assertIn('missing_timestep_member',r['reasons'])
    def test_foreign_model(self):
        r=pair(dict(target=.008,members={'0.002':dict(time=.008,failure='',signature='A93',accepted_metrics=[])},temporal={}))
        self.assertIn('0.002:foreign_identity',r['reasons'])
    def test_nonfinite_time(self):
        r=pair(dict(target=.008,members={'0.002':dict(time=float('nan'),signature='',accepted_metrics=[])},temporal={}))
        self.assertIn('0.002:incomplete_or_rejected',r['reasons'])
    def test_passing_label_cannot_hide_temporal_defect(self):
        r=pair(dict(target=.008,members={},temporal={'U':dict(status='pass',measured_drift=1e-5,difference=dict(weighted_RMS=1e-3,scale=1.))}))
        self.assertIn('temporal:U',r['reasons'])
if __name__=='__main__':unittest.main()
