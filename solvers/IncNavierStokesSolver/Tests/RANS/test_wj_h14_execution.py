# SPDX-License-Identifier: MIT
import copy
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET
from cls_epawf_m3_cases import set_info,set_parameter,write
from wj_h14_execution import protected,validate_pair,clone,target_preflight,QUALIFICATION,file_id,native_command,OwnerIndex

SOURCE=Path(__file__).resolve().parents[4]
PARENT=Path(os.environ.get('WJ_H14_PARENT',str(SOURCE/'docs/plans/wj-h13-linux/publication-v1/paired-case')))

class FixedHeightExecution(unittest.TestCase):
    def test_immutable_basis_cache_keeps_actual_bits_and_signed_zero_keys(self):
        import numpy as np
        from validate_cls_epawf_m3 import basis
        r=dict(fields=[dict(mesh=[dict(geometry=0,lower=[0,0,0],upper=[1,0,1])])]);cache=OwnerIndex(r)
        for p in(1,3):
            for x in(-1.,-.375,-0.,0.,.375,1.):
                expected=basis(p,[x]);a=cache.basis(p,[x]);b=cache.basis(p,[x])
                self.assertIs(a,b)
                for old,new in zip(expected,a):self.assertEqual(old.tobytes(),new.tobytes());self.assertFalse(new.flags.writeable)
        self.assertEqual(len(cache.basis_values),12)
    def test_index_preserves_every_physical_owner_and_roundoff_role(self):
        from wj_h13_domain import topology
        from wj_h4_points import owners,physical
        t=topology(PARENT/'flat/case.xml');r=dict(geometry_signature='Plane:CompositeTunnelV1:0.0875:20',parameters=[6.423173803526447e-5,.15],
            fields=[dict(mesh=[dict(e,geometry=e['id'])for e in t['elements']])])
        index=OwnerIndex(r)
        for x,z in [(-10.,0.),(0.,.15),(1.25,t['chart_normal_knots'][17]),(15.,20.),(0.,t['chart_normal_knots'][18])]:
            p=physical(r,[x,0,z]);self.assertEqual(index(p),owners(r,p))
            for a in(0,2):
                for s in(-1,1):
                    q=p.copy();q[a]+=s*1e-6
                    if -10<=q[0]<=15 and .0875<=q[2]<=20:self.assertEqual(index(q),owners(r,q))
        with self.assertRaises(ValueError):index([-10.001,0,.2])
    def test_native_command_is_bound_to_exact_inputs(self):
        with tempfile.TemporaryDirectory()as d:
            f=Path(d)/'case.xml';f.write_text('immutable execution intent')
            cmd=['governor',str(f)];a=dict(status='PASS',commands={'target_cp0':cmd},command_inputs={'target_cp0':[file_id(f)]})
            native_command(a,'target_cp0',cmd)
            with self.assertRaises(ValueError):native_command(a,'target_cp0',['other',str(f)])
            f.write_text('changed intent')
            with self.assertRaises(ValueError):native_command(a,'target_cp0',cmd)
    def test_closed_or_wrong_native_role_is_refused(self):
        for status,phase in [('NOT_RUN','target_cp0'),('PASS','governing')]:
            with self.assertRaises(ValueError):native_command(dict(status=status),phase,[])
    def test_stored_reference_never_refits_or_inverts_roughness(self):
        from study_wj_h14 import stored_reference
        with patch('wj_h5_neutral.required_roughness',side_effect=AssertionError('H14 must not calibrate')):
            c=stored_reference(PARENT/'paired-domain.json')
        self.assertEqual(c.H,.15);self.assertEqual(c.h,.10854503542535461)
        self.assertEqual(c.K,.008809015145515192)
        self.assertAlmostEqual(c.w['rc'],.026419898309769336,places=14)
    def test_execution_only_diff_and_exact_physical_functions(self):
        with tempfile.TemporaryDirectory()as d:
            a=PARENT/'flat/case.xml';b=clone(a,Path(d)/'case',dt=.0025)
            self.assertEqual(validate_pair(a,b)['status'],'PASS')
            self.assertEqual([ET.tostring(v)for v in ET.parse(a).findall('./CONDITIONS/FUNCTION')],
                             [ET.tostring(v)for v in ET.parse(b).findall('./CONDITIONS/FUNCTION')])
            self.assertEqual(protected(ET.parse(a).getroot()),protected(ET.parse(b).getroot()))
    def test_protected_mesh_and_time_cannot_be_relabelled(self):
        a=ET.parse(PARENT/'flat/case.xml').getroot()
        for key,value in [('Time',.0025),('Kinvis',1e-4)]:
            b=copy.deepcopy(a);set_parameter(b,key,value)
            with self.assertRaises(ValueError):validate_pair(a,b)
        b=copy.deepcopy(a);b.find('./GEOMETRY/VERTEX/V').text='-9.9 0 0'
        with self.assertRaises(ValueError):validate_pair(a,b)
    def test_no_height_generator_and_preserved_observation_roles(self):
        m=json.loads((PARENT/'paired-domain.json').read_text())
        self.assertEqual(m['primary_reference']['H_w'],.15)
        self.assertEqual(m['primary_reference']['physical_z_H_mm'],9.5)
        self.assertEqual(m['experimental']['physical_7mm_region'],'inner AWF')
        self.assertEqual(m['experimental']['intersection'],'EMPTY')
    def test_incomplete_backend_never_opens_target(self):
        a=PARENT/'flat/case.xml';m=json.loads((PARENT/'paired-domain.json').read_text())
        with tempfile.TemporaryDirectory()as d:
            b=clone(a,Path(d)/'case',dt=.0025)
            for missing in QUALIFICATION:
                q={k:'PASS'for k in QUALIFICATION};q[missing]='NOT_RUN'
                gate=target_preflight(a,b,m,q,{'status':'PASS'},{'status':'PASS','projected_total_seconds':3000})
                self.assertFalse(gate['native_launch_allowed']);self.assertIn(missing,gate['missing_gates'])
    def test_finite_pilot_and_closure_resources(self):
        a=ET.parse(PARENT/'flat/case.xml').getroot();b=copy.deepcopy(a);set_parameter(b,'NumSteps',5)
        with self.assertRaises(ValueError):validate_pair(a,b)
        with tempfile.TemporaryDirectory()as d:
            path=clone(PARENT/'flat/case.xml',Path(d)/'case',dt=.0025)
            q={k:'PASS'for k in QUALIFICATION};m=json.loads((PARENT/'paired-domain.json').read_text())
            gate=target_preflight(PARENT/'flat/case.xml',path,m,q,{'status':'PASS'},{'status':'PASS','projected_total_seconds':7201})
            self.assertFalse(gate['native_launch_allowed']);self.assertIn('resources',gate['missing_gates'])

if __name__=='__main__':unittest.main()
