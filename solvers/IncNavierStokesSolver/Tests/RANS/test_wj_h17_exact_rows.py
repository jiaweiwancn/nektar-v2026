# SPDX-License-Identifier: MIT
import copy
import json
from pathlib import Path
import tempfile
import unittest
import xml.etree.ElementTree as ET
from cls_epawf_m3_cases import set_parameter,set_info
from wj_h17_observer import equal,exact,native_load,phase_parity,independent_parity
from study_wj_h17 import protected,CAP_BYTES,CAP_ENTRIES,POLICIES
from run_wj_h17 import halt,reserve_strict,small_selection

class ExactRows(unittest.TestCase):
    def test_coordinator_command_matches_observer_parser(self):
        from study_wj_h17 import observer_command
        from wj_h17_observer import argument_parser
        command=observer_command(Path('/run'),Path('/run/intent'),Path('/run/compat'))
        args=argument_parser().parse_args(command[3:])
        self.assertEqual((args.root,args.intent,args.compat,args.behavior),(Path('/run'),Path('/run/intent'),Path('/run/compat'),'PASS'))

    def test_exact_small_commands_cannot_admit_CP0_or_replacement_work(self):
        import sys
        import run_wj_h17
        root=Path('/declared/H17');cmd=[sys.executable,'-B',str(Path(run_wj_h17.__file__).with_name('study_wj_h17.py')),'--root',str(root),'--mode','small','--case','A']
        self.assertEqual(small_selection(root,cmd,2),'A')
        with self.assertRaises(AssertionError):small_selection(root,cmd,1)
        for bad in (cmd[:-1]+['C'],cmd[:-3]+['setup'],['governor','NumSteps=0']):
            with self.assertRaisesRegex(AssertionError,'undeclared'):small_selection(root,bad,2)

    def test_bitwise_stream_comparison_retains_signed_zero_ULP_order(self):
        import math
        for a,b in (([0.],[-0.]),([1.],[math.nextafter(1.,2.)]),([17.,18.],[18.,17.]),({'lower':1.},{'upper':1.})):
            with self.assertRaisesRegex(AssertionError,'bitwise'):equal(a,b,'fixture')
        equal({'ids':[17,18],'g':[-0.,1.]},{'g':[-0.,1.],'ids':[17,18]},'ordered keys')
        with tempfile.TemporaryDirectory()as d:
            p=Path(d)/'native.json';p.write_text('{"v":[-0,0,1.0000000000000002]}')
            self.assertNotEqual(exact(native_load(p)['v'][0]),exact(native_load(p)['v'][1]))

    def test_execution_policy_excludes_only_execution_fields(self):
        x=ET.fromstring('<NEKTAR><CONDITIONS><PARAMETERS><P>Time=0</P><P>TimeStep=.0025</P><P>NumSteps=4</P><P>RANSWallHeight=.15</P></PARAMETERS><SOLVERINFO><I PROPERTY="RANSWallEvaluationReuse" VALUE="ExactStepV1"/></SOLVERINFO><FUNCTION NAME="Initial"><E VAR="u" VALUE="frozen"/></FUNCTION></CONDITIONS><GEOMETRY DIM="2"/></NEKTAR>')
        z=copy.deepcopy(x)
        for policy in POLICIES:
            set_info(z,'RANSWallImmutableRowReuse',policy);set_info(z,'RANSWallImmutableRowEntries',CAP_ENTRIES);set_info(z,'RANSWallImmutableRowBytes',CAP_BYTES);set_parameter(z,'NumSteps',0)
            self.assertEqual(protected(x),protected(z))
        for k,v in (('Time',.01),('TimeStep',.005),('RANSWallHeight',.2)):
            bad=copy.deepcopy(z);set_parameter(bad,k,v);self.assertNotEqual(protected(x),protected(bad))
        self.assertLessEqual(CAP_BYTES,512*1024**2);self.assertLessEqual(CAP_ENTRIES,524288)

    def test_native_replay_requires_callback_reference_and_actual_history(self):
        r=dict(schema='WJCompletedStepV1',time=.0025,steps=1,dt=.0025,signature='original',last_stage_alpha=400.,last_stage_time=.0025,last_stage_history=[-0.,1.],state=[2.],fixed_reference=[3.],fixed_defect=[4.],native_history=[[2.]],native_history_times=[.0025],stage_audits=[[.0025,0.,400.]],step_audits=[[.0025]],explicit_audits=[[0.,1.]],first_retained_explicit_call=1,actual_explicit_calls=1,native_hard_audit={'admissible':[1.]},accepted_states=[[2.]],wj_audit_packets=[dict(state_sha256='state',stage_history_sha256='callback',time=.0025,alpha_mass=400.,principal_witnesses={'upper':[17.,1.]})])
        phase_parity(r,copy.deepcopy(r))
        for key,value in (('last_stage_history',[0.,1.]),('native_history',[[2.0000000000000004]]),('last_stage_alpha',0.),('fixed_defect',[0.])):
            bad=copy.deepcopy(r);bad[key]=value
            with self.assertRaises(AssertionError):phase_parity(bad,r)
        bad=copy.deepcopy(r);bad['wj_audit_packets'][0]['principal_witnesses']={'upper':[18.,1.]}
        with self.assertRaises(AssertionError):phase_parity(bad,r)

    def test_live_inventory_parity_never_substitutes_final_union(self):
        inv=dict(rule='frozen',points=[[1.,0.,.15]],owned=[{'geometry':17,'side':'upper'},{'geometry':18,'side':'lower'}],physical_rows=3,quadrature_rows=1)
        result=dict(status='PASS',summary={'m':.003},all_regional_Q1_moments={'inner':[-0.,1e-16],'outer':[2e-17]})
        independent_parity(result,inv,copy.deepcopy(result),copy.deepcopy(inv))
        for update in ({'points':inv['points']+[[2.,0.,.15]]},{'owned':list(reversed(inv['owned']))},{'quadrature_rows':0}):
            bad=dict(inv,**update)
            with self.assertRaises(AssertionError):independent_parity(result,bad,result,inv)
        bad=copy.deepcopy(result);bad['all_regional_Q1_moments']['inner'][0]=0.
        with self.assertRaises(AssertionError):independent_parity(bad,inv,result,inv)

    def test_halt_and_final_readback_are_monotone_separate_roles(self):
        with tempfile.TemporaryDirectory()as d:
            p=Path(d);(p/'accounting.json').write_text(json.dumps(dict(closed=False,full_audit_reserved=12,small_charged_attempts=3)))
            halt(p,'qualification_halt',dict(first='owner17'));halt(p,'qualification_halt',dict(first='later'))
            reserve_strict(p,'target')
            with self.assertRaises(AssertionError):reserve_strict(p,'target')
            b=json.loads((p/'accounting.json').read_text());self.assertEqual(b['qualification_halt']['first'],'owner17');self.assertEqual(b['full_audit_reserved'],12);self.assertEqual(b['small_charged_attempts'],3)

if __name__=='__main__':unittest.main()
