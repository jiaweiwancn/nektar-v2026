# SPDX-License-Identifier: MIT
import math,unittest
from study_wj_h25 import measure_matches
from seed_wj_h25 import checked_seed
from analyze_wj_h25_profile import category

class ProductionStrategy(unittest.TestCase):
    def test_derived_physical_measure_roundoff_is_not_state_payload_equality(self):
        self.assertTrue(measure_matches(25*(20-.0875-.15),494.0625))
        self.assertTrue(measure_matches(math.nextafter(497.8125,math.inf),497.8125))
        self.assertFalse(measure_matches(25*(20-.15),494.0625))
        self.assertFalse(measure_matches(494.0625+1e-7,494.0625))

    def test_invalid_or_wrong_measures_fail(self):
        for x in (math.nan,math.inf,-math.inf,3.75):self.assertFalse(measure_matches(x,497.8125))

    def seed_fixture(self):
        cp=dict(steps=0,time=0.,dt=.005,signature='own .005',state=[1.,-0.,2.],
            fixed_reference=[3.,-0.],fixed_defect=[4.,0.],eta=0.,
            native_history=[[1.,-0.],[0.,2.]],native_history_times=[0.,.005])
        seed=dict(cp,schema='WJH25NativeSeedV2',status='PASS',governing_solves=0,
            initial_velocity_applications=0,native_history_initializations=1,
            spatial_work=[0,0],perturbed=False,last_stage_alpha=0,last_stage_history=[])
        return seed,cp,dict(cp)

    def test_new_history_is_own_dt_and_original_spatial_bits_are_preserved(self):
        seed,cp,source=self.seed_fixture();source['dt']=.0025;source['signature']='old .0025'
        source['native_history_times']=[0.,.0025]
        checked_seed(seed,cp,source,False,.005)
        cp['native_history_times']=[0.,.0025]
        with self.assertRaises(AssertionError):checked_seed(seed,cp,source,False,.005)

    def test_signed_zero_spatial_mismatch_is_not_a_timestep_tolerance(self):
        seed,cp,source=self.seed_fixture();source['state']=[1.,0.,2.]
        with self.assertRaises(AssertionError):checked_seed(seed,cp,source,False,.005)

    def test_evolved_or_callback_role_cannot_be_relabelled_initial(self):
        for change in (dict(steps=1),dict(time=.005),dict(last_stage_alpha=200),dict(last_stage_history=[1.,2.,3.])):
            seed,cp,source=self.seed_fixture();seed.update(change)
            with self.assertRaises(AssertionError):checked_seed(seed,cp,source,False,.005)

    def test_new_transient_does_not_recompute_internal_reference(self):
        seed,cp,source=self.seed_fixture();cp['state']=[1.01,-0.,2.01];seed.update(state=list(cp['state']),
            perturbed=True,scalar_fraction=.01,velocity_scale=.001,spatial_work=[1,1])
        checked_seed(seed,cp,source,True,.005)
        source['fixed_reference']=[3.01,-0.]
        with self.assertRaises(AssertionError):checked_seed(seed,cp,source,True,.005)

    def test_history_initialization_and_governing_counts_are_distinct(self):
        for change in (dict(governing_solves=1),dict(initial_velocity_applications=1),dict(native_history_initializations=0),dict(spatial_work=[1,1])):
            seed,cp,source=self.seed_fixture();seed.update(change)
            with self.assertRaises(AssertionError):checked_seed(seed,cp,source,False,.005)

    def test_guards_cannot_inflate_equation_work_fraction(self):
        self.assertEqual(category('domain.parts.moving_ray_rows',['domain.compute','native.advance']),
            'guards_bounds_residual_verification')
        self.assertEqual(category('mapped.band_solve',['audit.WJ.total']), 'full_scientific_audit')
        self.assertEqual(category('pressure.build',['pimple.sweep']), 'equation_evaluation_assembly_solution_update')

if __name__=='__main__':unittest.main()
