#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Custom-normal physical rules and immutable requested-intent fixtures."""
import copy
from pathlib import Path
import tempfile
import unittest
import numpy as np
from numpy.testing import assert_allclose
from wj2_cases import local_normal_knots,validated_knots,normal_knots,case_xml
from wj4_measure import rule,boundary_rule,owned_rule
from wj_immutable import write_once,xml_identity

class Meshes(unittest.TestCase):
    def test_retained_metadata_integer_fingerprint(self):
        from run_wj5_native import read_execution_intent
        from wj_immutable import digest,read
        import json
        with tempfile.TemporaryDirectory() as t:
            path=Path(t)/'intent.json';intent=dict(steps=32,mesh=dict(level=1,nx=2),runtime=dict(bytes=1088),dt=3.2e-5)
            path.write_text(json.dumps(intent))
            self.assertEqual(digest(read_execution_intent(path)),digest(intent))
            self.assertNotEqual(digest(read(path)),digest(intent))
    def test_local_mesh_is_refinement_only(self):
        l1=normal_knots(1);m1=local_normal_knots('M1');m2=local_normal_knots('M2')
        self.assertEqual((len(l1),len(m1),len(m2)),(31,35,43))
        self.assertTrue(set(l1)<set(m1)<set(m2))
        self.assertEqual([v for v in m1 if v not in l1],[.085,.095,.125,.175])
        self.assertEqual([v for v in m2 if v not in m1],[.0825,.0875,.0925,.0975,.1125,.1375,.1625,.1875])
    def test_invalid_knots(self):
        for knots in ([0,.1,1,2,2],[0,1,2],[0,.1,float('nan'),1,2],[0,.1,1,1.5],[-1,.1,1,2]):
            with self.assertRaises(ValueError):validated_knots(knots,0.,2.,(.1,1.))
    def test_actual_composite_and_surface_measures(self):
        meshes=[local_normal_knots('M1'),local_normal_knots('M2')]
        for order,xorder in ((6,4),(8,6)):
            rows=rule(zmeshes=meshes,order=order,xorder=xorder,all_faces=True)
            q=np.array([r for r in rows if r[3]>0]);x,z,w=q[:,0],q[:,2],q[:,3]
            assert_allclose([sum(w),w@x**6,w@z**6],[2,2/7,128/7],rtol=3e-14,atol=3e-14)
            surface=boundary_rule(zmeshes=meshes,order=order,xorder=xorder)
            assert_allclose(sum(surface.values()),5.,rtol=0,atol=3e-14)
            self.assertTrue(all(r[3]==0 for r in rows if r[6]))
        physical=rule(zmeshes=meshes,all_faces=True)
        rows=owned_rule(physical,1,[0,.5,1],meshes[0])
        for z in meshes[0][1:-1]:
            corner=[r for r in rows if r[:3]==[.5,0.,z]]
            self.assertEqual(len({r[5] for r in corner}),4)
            self.assertEqual({(r[6],r[7]) for r in corner},{(-1,-1),(-1,1),(1,-1),(1,1)})
    def test_knots_and_comparison_role_are_immutable(self):
        with tempfile.TemporaryDirectory() as t:
            path=Path(t)/'intent.json';a=dict(level=1,zknots=local_normal_knots('M1'),role='horizontal',dt=3.2e-5,target=128)
            write_once(path,a);write_once(path,a)
            for change in (dict(zknots=local_normal_knots('M2')),dict(role='timestep'),dict(dt=1.6e-5),dict(target=32)):
                b=copy.deepcopy(a);b.update(change)
                with self.assertRaisesRegex(ValueError,'immutable intent'):write_once(path,b)
    def test_case_defaults_and_custom_identity(self):
        with tempfile.TemporaryDirectory() as t:
            base=case_xml(Path(t)/'a',level=1)
            same=case_xml(Path(t)/'b',level=1,zknots=normal_knots(1))
            custom=case_xml(Path(t)/'c',level=1,zknots=local_normal_knots('M1'))
            self.assertNotEqual(xml_identity(same)['semantic_sha256'],xml_identity(custom)['semantic_sha256'])
            # Checkpoint destinations are intentionally case-specific.
            import xml.etree.ElementTree as ET
            for xml in (base,same):
                root=ET.parse(xml).getroot();self.assertEqual(len(root.findall('./GEOMETRY/ELEMENT/Q')),30)

class ConstraintUncertainty(unittest.TestCase):
    def test_fixed_free_prescribed_and_pressure_rows(self):
        from wj5_measure import constraint_estimate
        fixed=dict(row=[[0,1.]],free_row=[],free_row_L2=0.,gross_row_L2=1.,absolute_accumulation=7.,prescribed_accumulation=7.)
        a=constraint_estimate(fixed,1e-2)
        self.assertGreater(a['total'],0.);self.assertLess(a['total'],1e-9)
        self.assertGreater(constraint_estimate(fixed,1e-2,prescribed_uncertainty=1e-6)['total'],1e-5)
        free=dict(fixed,free_row=[[4,.5]],free_row_L2=.5)
        self.assertGreaterEqual(constraint_estimate(free,1e-2)['total'],.005)
        self.assertGreaterEqual(constraint_estimate(free,1e-8,free_state_bound=.201)['total'],.201)
        pressure=dict(fixed,free_row=[[8,1.]],free_row_L2=1.)
        self.assertGreaterEqual(constraint_estimate(pressure,1e-2)['total'],.01)
        for key in ('free_row_L2','absolute_accumulation'):
            invalid=dict(fixed);invalid[key]=float('nan')
            with self.assertRaisesRegex(ValueError,'invalid constraint'):constraint_estimate(invalid,1e-2)
    def test_dependent_epsilon_uses_its_k_input_uncertainty(self):
        from wj5_measure import boundary_uncertainty
        from wj3_measure import FIELDS
        fixed=dict(block=0,row=[[0,1.]],free_row=[],free_row_L2=0.,gross_row_L2=1.,absolute_accumulation=1.,prescribed_accumulation=1.)
        epsilon=dict(fixed,block=3,row=[[30,1.]],free_row=[[20,3.],[30,4.]],free_row_L2=5.)
        pressure=dict(fixed,block=4,row=[[100,1.]],free_row=[[100,1.]],free_row_L2=1.)
        m=dict(observations=[dict(position=[0.,0.,1.],geometry=0,inner=False,constraint_rows=[fixed,fixed,fixed,epsilon,pressure])],
            dependent_trace_action=[dict(row=300,free_K_action=[[20,3.]])],arithmetic_evaluation_sensitivity=[0.]*6,
            passive_diagnostics=dict(pressure_regions=[dict(offset=100,gauge_mean_row=[.5]),dict(offset=200,gauge_mean_row=[.5])]))
        packet=dict(systems={k:dict(physical_metric_min_eigenvalue=1.) for k in ('momentum','k','epsilon')},
            estimated_comparison_floor={f:1e-9 for f in FIELDS+('velocity',)},estimated_sampled_comparison_floor={f:1e-9 for f in FIELDS+('velocity',)})
        packet['systems']['pressure']=dict(physical_mass_operator_bound=1.)
        packet['estimated_comparison_floor']['k']=1e-3
        first=boundary_uncertainty(m,packet,{})['fields']['outer_epsilon']['sampled_Linf']
        packet['estimated_comparison_floor']['k']=2e-3
        second=boundary_uncertainty(m,packet,{})['fields']['outer_epsilon']['sampled_Linf']
        self.assertGreaterEqual(first,.003+4e-9)
        self.assertAlmostEqual(second-first,.003,places=14)
        self.assertGreater(boundary_uncertainty(m,packet,{})['fields']['reduced_pressure']['sampled_Linf'],0.)

    def test_every_boundary_row_map_is_propagated_and_hashed(self):
        from wj5_measure import boundary_envelope
        from wj3_measure import FIELDS
        from wj_immutable import file_hash
        from unittest.mock import patch
        import json
        packet=dict(safety_factor=100,systems={k:dict(physical_metric_min_eigenvalue=1.) for k in ('momentum','k','epsilon')},
            estimated_comparison_floor={f:1e-9 for f in FIELDS+('velocity',)},estimated_sampled_comparison_floor={f:1e-8 for f in FIELDS+('velocity',)})
        packet['systems']['pressure']=dict(physical_mass_operator_bound=1.)
        with tempfile.TemporaryDirectory() as t:
            records=[]
            for step in (0,32,64):
                path=Path(t)/(str(step)+'.json');path.write_text(json.dumps(dict(measurements=dict(marker=step))))
                c=dict(path=str(path),sha256=file_hash(path),steps=step,state_sha256=str(step),rule_settings={})
                if step==32:c['arithmetic']=packet
                records.append(c)
            seen=[]
            def estimate(m,a,settings):
                seen.append(m['marker']);return dict(fields={f:dict(L2=1e-10,sampled_Linf=1e-10) for f in FIELDS+('velocity',)},maximum_boundary_free_rows={})
            with patch('wj5_measure.boundary_uncertainty',side_effect=estimate):
                _,_,proof=boundary_envelope(dict(mesh=dict(initial=records[0],trajectory=records[1:])))
                self.assertEqual(seen,[0,32,64]);self.assertEqual(len(proof['coverage']),3)
                Path(records[-1]['path']).write_text('{}')
                with self.assertRaisesRegex(ValueError,'boundary observation bytes changed'):boundary_envelope(dict(mesh=dict(initial=records[0],trajectory=records[1:])))

    def test_complete_scope_estimate_is_required(self):
        from wj4_qualification import sensitivity_sequence,FIELDS,SCOPES
        design=dict(amplitudes={f:1. for f in FIELDS},budgets=dict(L2=1e-4,sampled_Linf=1e-3),sensitivity_fraction=.1)
        floor={f:dict(L2=1e-8,sampled_Linf=1.59e-5) for f in FIELDS}
        pair={scope:dict(fields={f:dict(measure=1.,normalized_L2=0.,normalized_sampled_Linf=0.) for f in FIELDS}) for scope in SCOPES}
        # Original global boundary estimate stays unresolved. Only a separately
        # justified complete estimate can change the new classification.
        self.assertEqual(sensitivity_sequence({'0':pair},floor,design,['0'])['status'],'UNRESOLVED')
        tight={f:dict(L2=1e-9,sampled_Linf=1e-9) for f in FIELDS}
        result=sensitivity_sequence({'0':pair},floor,design,['0'],scope_uncertainty={'constrained_boundary':tight})
        self.assertEqual(result['status'],'PASS')
        del tight['k']
        self.assertEqual(sensitivity_sequence({'0':pair},floor,design,['0'],scope_uncertainty={'constrained_boundary':tight})['status'],'FAIL')

if __name__=='__main__':unittest.main()
