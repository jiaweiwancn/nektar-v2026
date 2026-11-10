#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Adversarial pass-labelled evidence must not satisfy a configuration gate."""
import copy,unittest
from rans_duration_contract import aggregate, required_times, scalar_coverage, MEMBERS


class DurationContract(unittest.TestCase):
    def setUp(self):
        self.config=dict(model='A93',policy='Natural',wall='Accepted',audit_times=[.016,.1,.5])
        def row(t, norm):
            return dict(time=t,analytical_weak_balance=dict(status='pass'),spatial=dict(status='pass',
                native_full_mass_dual=dict(u=norm,w=norm,k=norm,epsilon_outer=1e-9),
                native_k_partition_mass_dual=dict(inner=norm)),common_test=dict(fields={n:dict(fine_restricted_common_mass_dual=norm) for n in ('u','w','k')}))
        self.good=dict(configuration=self.config,companions={m:dict(configuration=self.config,status='pass',scalar_coverage=dict(last_time=.5)) for m in MEMBERS},spatial_audits={m:dict(configuration=self.config,status='pass',checks=[row(t,n) for t in self.config['audit_times']]) for m,n in [('B5',1e-3),('B11',1e-4)]})

    def reject(self, receipt, reason):
        r=aggregate(self.config,receipt)
        self.assertEqual(r['status'],'fail')
        self.assertTrue(any(reason in e for e in r['failures']),r)

    def test_positive(self):self.assertEqual(aggregate(self.config,self.good)['status'],'pass')
    def test_missing_time_despite_pass(self):
        r=copy.deepcopy(self.good);r['spatial_audits']['B11']['checks'].pop(1);self.reject(r,'missing_required_time')
    def test_foreign_configuration(self):
        for key,value in [('model','E94E95'),('wall','EP')]:
            r=copy.deepcopy(self.good);r['configuration']=dict(r['configuration'],**{key:value});self.reject(r,'configuration_mismatch')
    def test_missing_each_companion(self):
        for member in MEMBERS:
            r=copy.deepcopy(self.good);del r['companions'][member];self.reject(r,'missing_companion')
    def test_short_scalar_observer(self):
        r=copy.deepcopy(self.good);r['companions']['B11-dt2']['scalar_coverage']['last_time']=.1;self.reject(r,'scalar_horizon_incomplete')
    def test_weak_and_reconstruction_do_not_replace_reduction(self):
        r=copy.deepcopy(self.good)
        for row in r['spatial_audits']['B11']['checks']:row['spatial']['native_full_mass_dual']['u']=1.1e-3
        self.reject(r,'evolved_reduction_failed')
    def test_both_short_times(self):
        with self.assertRaisesRegex(ValueError,'missing_required_time'):required_times([dict(time=.016)],(.0016,.016))
        required_times([dict(time=.0016+1e-17),dict(time=.016-1e-17)],(.0016,.016))
    def test_each_scalar_state_identity(self):
        from rans_qualification import exact_vector_hash
        report=dict(audit_times=[0.,.001,.002],accepted_states=[[0.],[1.],[2.]],time=.002,
                    pimple_attempts=[dict(alpha=1000.),dict(alpha=1000.)])
        observations=[dict(stage='accepted',time=t,state_sha256=exact_vector_hash(state))
                      for t,state in zip(report['audit_times'],report['accepted_states'])]
        observer=dict(status='pass',configuration=self.config,observations=observations)
        self.assertEqual(scalar_coverage(observer,[report],self.config,.002)['unique_accepted_states'],3)
        del observations[-1]['state_sha256']
        with self.assertRaisesRegex(ValueError,'scalar_state_identity_missing'):
            scalar_coverage(observer,[report],self.config,.002)


if __name__=='__main__':unittest.main()
