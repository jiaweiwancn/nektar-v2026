#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
import json
import sys
import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch
import numpy as np
from numpy.testing import assert_allclose
import run_wj_h7 as journal
import run_wj_h8 as admission
from study_wj_h7 import clone,compatible
from study_wj_h6 import configuration
from wj_h5_cases import case_xml
from wj_h8_analysis import PressureProjection,ratio,strict_context
from wj_h7_analysis import Rows,principal_batch
from wj_h5_neutral import NeutralReference
from wj_h6_projection import project_field
from test_wj_h6_projection import fixture
from wj_h4_points import vector_hash,aligned

class LocalTimestepTests(unittest.TestCase):
    def test_child_halt_survives_successful_parent_and_refuses_native(self):
        # A real child publishes a halt and another reservation while execute()
        # still owns the pre-spawn journal. This reproduced the original bug.
        with tempfile.TemporaryDirectory()as d:
            root=Path(d)/'wj-h7-neutral-transport-20261113';root.mkdir()
            script=Path(d)/'child.py'
            script.write_text('from pathlib import Path\nfrom run_wj_h7 import record_coverage_halt,update_journal\np=Path('+repr(str(root/'accounting.json'))+')\nrecord_coverage_halt(p,dict(status="OBSTRUCTED",role="accepted",trajectory="6/fine",time=.03,step=6,witness={"marker":"first"}))\nupdate_journal(p,lambda b:b["entries"].update(child_reserved=dict(status="RESERVED",reserved_seconds=7)))\n')
            import os
            with patch.dict(os.environ,PYTHONPATH=str(Path(journal.__file__).parent)+os.pathsep+os.environ.get('PYTHONPATH','')),patch.object(journal,'scientific_workers',return_value=dict(workers=0)):
                self.assertTrue(journal.execute(root,'child-halt',[sys.executable,str(script)],seconds=10))
            b=journal.load(root/'accounting.json');self.assertEqual(b['coverage_halt']['witness']['marker'],'first');self.assertEqual(b['entries']['child_reserved']['reserved_seconds'],7)
            self.assertEqual(b['entries']['child-halt']['status'],'PASS')
            for cp0 in (True,False):
                with self.assertRaisesRegex(AssertionError,'coverage obstruction'):journal.execute(root,'refused'+str(cp0),['IncNavierStokesSolver','unused.xml'],cp0=cp0,governing=not cp0,seconds=1)

    def test_first_witness_and_earliest_physical_time_are_distinct(self):
        with tempfile.TemporaryDirectory()as d:
            p=Path(d)/'book.json';journal.save(p,dict(entries={'parent':dict(status='RESERVED',reserved_seconds=3)}))
            first=dict(status='UNRESOLVED',role='accepted',trajectory='6/fine',time=.03,step=6,witness={'id':1})
            journal.record_coverage_halt(p,first);journal.record_coverage_halt(p,dict(first,time=.01,step=2,witness={'id':2}))
            journal.complete_operation(p,'parent',dict(status='PASS',seconds=.1));b=journal.load(p)
            self.assertEqual(b['coverage_halt'],first);self.assertEqual(b['coverage_earliest']['6/fine']['time'],.01)
            with self.assertRaises(ValueError):journal.record_coverage_halt(p,dict(first,role='rejected'))

    def test_mixed_statuses_and_rejected_candidate_are_not_passed(self):
        for status in ('OBSTRUCTED','UNRESOLVED'):
            r=[dict(role='accepted',status='PASS'),dict(role='accepted',status=status),dict(role='rejected',status='QUALITY_LOST')]
            self.assertEqual(journal.coverage_quality(r),dict(coverage=status,quality=status))
        self.assertEqual(journal.coverage_quality([dict(role='rejected',status='PASS')])['coverage'],'NOT_RUN')
        self.assertEqual(journal.coverage_quality([dict(role='accepted',status='QUALITY_LOST')])['quality'],'QUALITY_LOST')

    def test_fixed_case_cp0_and_sixteen_attempt_design(self):
        with tempfile.TemporaryDirectory()as d:
            p=Path(d);src=case_xml(p/'src',configuration(6.),bisect_once=True)
            x=clone(src,p/'fine.xml',.0025,0);self.assertEqual(admission.admit_xml(x,True)[1],0)
            x=clone(src,p/'run.xml',.0025,16,p/'CP');self.assertEqual(admission.admit_xml(x)[1],16)
            for dt,n in ((.005,16),(.0025,17),(.0025,8)):
                clone(src,x,dt,n,p/'CP')
                with self.assertRaises(AssertionError):admission.admit_xml(x)
            wrong=case_xml(p/'wrong',configuration(5.5),bisect_once=True);clone(wrong,x,.0025,16,p/'CP')
            with self.assertRaises(AssertionError):admission.admit_xml(x)

    def test_reservation_view_cp0_and_attempt_ceilings(self):
        with tempfile.TemporaryDirectory()as d:
            root=Path(d)/'wj-h8-local-timestep-20261113';root.mkdir()
            src=case_xml(Path(d)/'src',configuration(6.),bisect_once=True);x=clone(src,Path(d)/'cp0.xml',.0025,0)
            book=dict(seconds_ceiling=1200,setup_seconds=60,closure_reserve=180,entries={},cp0_reservations=1,native_reservations=[],strict_view_reservations=7)
            journal.save(root/'accounting.json',book)
            with self.assertRaisesRegex(AssertionError,'one fresh'):admission.execute(root,'cp0',['IncNavierStokesSolver',str(x)],cp0=True,seconds=1)
            with self.assertRaisesRegex(AssertionError,'seven strict'):admission.execute(root,'view',['true'],strict_view=True,seconds=1)
            book['entries']['unknown']=dict(status='RESERVED',reserved_seconds=1,units=16);journal.save(root/'accounting.json',book);clone(src,x,.0025,16,root/'CP')
            with self.assertRaisesRegex(AssertionError,'16 governing'):admission.execute(root,'run',['IncNavierStokesSolver',str(x)],governing=True,seconds=1)

    def test_cross_dt_state_and_reference_are_not_relabelled(self):
        a=dict(dt=.005,state=[-0.,1.],fixed_reference=[-0.,2.],signature='five',steps=0,time=0)
        with self.assertRaisesRegex(ValueError,'timestep'):compatible(a,a,.0025)
        compatible(a,dict(a,dt=.0025,signature='quarter'),.0025,same_dt=False)
        for key in ('state','fixed_reference'):
            with self.assertRaisesRegex(ValueError,'state|reference'):compatible(a,dict(a,**{key:[0.,1.]}),.005)

    def test_checkpoint_names_do_not_define_absolute_steps(self):
        from study_wj_h8 import checkpoint_index
        expected={'/run/CP.0/State.bin':{'sha256':'first'},'/run/CP.1/State.bin':{'sha256':'second'}}
        with patch('study_wj_h8.checked'),patch('study_wj_h8.checkpoint_payload',side_effect=lambda p:dict(steps=1 if p.name=='CP.0'else 2)):
            index=checkpoint_index(Path('/run/case.evolution.json'),expected)
        self.assertEqual(index[1][0],Path('/run/CP.0'));self.assertEqual(index[2][0],Path('/run/CP.1'))

    def test_genuine_intermediate_callback_and_rejected_ineligibility(self):
        s=[[0.],[1.],[2.]];hist=[[0.],[.8],[1.8]];packets=[dict(state_sha256=vector_hash(v),stage_history_sha256=vector_hash(h),callback_history=h,alpha_mass=0 if i==0 else 400,time=i*.0025)for i,(v,h)in enumerate(zip(s,hist))]
        r=dict(dt=.0025,steps=2,accepted_states=s,wj_audit_packets=packets,physical_range_events=[dict(state=[3.],history=[2.8],alpha=400,time=.0075,wj_packet={})])
        c=strict_context(r,1);self.assertEqual(c['history'],hist[1]);self.assertEqual(c['alpha'],400)
        with self.assertRaises(ValueError):strict_context(r,3)
        with self.assertRaisesRegex(ValueError,'binding'):strict_context(dict(r,accepted_states=[[0.],[7.],[2.]]),1)

    def test_common_integer_times_and_zero_safe_three_level_ratio(self):
        self.assertEqual(aligned(3,.01,.03),aligned(12,.0025,.03+2e-17));self.assertNotEqual(aligned(11,.0025,.0275),aligned(12,.0025,.03))
        self.assertEqual(ratio(.002,.0008,4e-11)['rho'],.4)
        self.assertEqual(ratio(0,0,4e-11)['status'],'UNRESOLVED');self.assertEqual(ratio(.001,.002,4e-11)['status'],'DOES_NOT_CONTRACT')
        with self.assertRaises(ValueError):ratio(float('nan'),0,1e-11)

    def test_actual_pressure_spaces_polynomial_and_complement_identity(self):
        r=fixture();r['geometry_signature']='fixture:CompositeTunnelV1:0.2:2';c=NeutralReference(.5,.0056762703212168,d0=.2,top=2)
        rows=Rows(r,[[0,0,.7],[8,0,.7]],c);projection=PressureProjection(rows)
        x=np.array([p[0]for p,o in rows.rows]);data=projection.measure(dict(divergence=1+.3*x))
        for v in data['regions'].values():
            self.assertLess(v['continuous_Q1_complement_L2'],4e-14);self.assertLess(v['broken_Q1_complement_L2'],4e-14);self.assertLess(v['identity_roundoff'],1e-12)
        data=projection.measure(dict(divergence=x*x))
        for block,v in data['regions'].items():
            height=.5 if block=='4'else 1.3
            self.assertAlmostEqual(v['broken_Q1_complement_L2']**2,8**5*height/180,places=10)
            self.assertLess(v['identity_roundoff'],2e-11);self.assertLess(v['broken_identity_roundoff'],2e-11)
        self.assertEqual(data['governing_solves'],0);self.assertEqual(data['mass_solves'],5)

    def test_early_noncontraction_is_not_hidden_by_terminal_pass(self):
        from publish_wj_h8 import contraction_status
        primary=[dict(time=t,geometry=17,rho=r)for t,r in zip((.01,.02,.03,.04),(1.86,.82,.31,.05))]
        envelope=[dict(time=v['time'],rho=v['rho'])for v in primary]
        self.assertEqual(contraction_status(primary,envelope)['status'],'UNRESOLVED')
        self.assertEqual(contraction_status(primary[1:],envelope)['status'],'NOT_RUN')
        primary[0]['rho']=.4;envelope[0]['rho']=.4
        self.assertEqual(contraction_status(primary,envelope)['status'],'PASS')

    def test_direction_detail_does_not_change_principal_action(self):
        G=np.array([[[.02,0,.7],[0,.01,0],[.01,0,-.04]]]);k=np.array([.009]);e=np.array([.002])
        a=principal_batch(k,e,G);b=principal_batch(k,e,G,details=True)
        for u,v in zip(a,b[:6]):assert_allclose(u,v,rtol=0,atol=2e-14)
        self.assertAlmostEqual(np.linalg.norm(b[6]['transverse_xi'][0]),1.)

if __name__=='__main__':unittest.main()
