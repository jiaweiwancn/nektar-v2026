# SPDX-License-Identifier: MIT
import copy
import json
import os
from pathlib import Path
import tempfile
import unittest
import xml.etree.ElementTree as ET
from cls_epawf_m3_cases import set_info,set_parameter,write
from wj_h14_execution import protected,validate_pair,clone,target_preflight,QUALIFICATION

SOURCE=Path(__file__).resolve().parents[4]
PARENT=Path(os.environ.get('WJ_H14_PARENT',str(SOURCE/'docs/plans/wj-h13-linux/publication-v1/paired-case')))

class FixedHeightExecution(unittest.TestCase):
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
