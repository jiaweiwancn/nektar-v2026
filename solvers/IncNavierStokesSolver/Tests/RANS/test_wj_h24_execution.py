# SPDX-License-Identifier: MIT
import copy,json,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import xml.etree.ElementTree as ET
from wj_h15_observer import identity
from wj_h23_continuation import ADDITIONAL,ORDER,active_execution_bindings
from run_wj_h24 import command_kind,reserve_native,frozen,NATIVE
from study_wj_h22 import protected_case

def put(path,value):
    path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value));return identity(path)

class H24Execution(unittest.TestCase):
    def test_only_unchanged_three_branches_and_own_terminal_reads(self):
        root=Path('/home/jay/build/wj-h24-flat-continuation-20261115')
        base=[sys.executable,'-B',str(Path(__file__).with_name('study_wj_h24.py')),'--root',str(root),'--mode']
        for k in NATIVE:self.assertEqual(command_kind(root,base+[k]),k)
        for k in ('cp0','bootstrap','build','probe','H22_replay','M@unregistered'):
            with self.assertRaises(AssertionError):command_kind(root,base+[k])
        b=dict(entries={},native_reservations=[],reserved_attempts=0)
        for k in NATIVE:reserve_native(b,k)
        self.assertEqual(b['reserved_attempts'],28)
        with self.assertRaises(AssertionError):reserve_native(b,'M')

    def test_frozen_files_do_not_substitute_for_actual_import_check(self):
        with tempfile.TemporaryDirectory()as t:
            root=Path(t);pin=put(root/'compatibility.json',{})
            put(root/'execution-config.json',dict(files=[pin],compatibility=pin))
            with patch('wj_h23_continuation.active_execution_bindings',side_effect=AssertionError('stale actual path'))as guard:
                with self.assertRaisesRegex(AssertionError,'stale actual path'):frozen(root)
                guard.assert_called_once_with(pin)

    def test_same_import_bytes_at_other_paths_need_a_separate_record(self):
        with tempfile.TemporaryDirectory()as t:
            root=Path(t);historical=put(root/'historical.json',{})
            pins={k:put(root/(k+'.json'),{})for k in ('active_adapter','new_native_identity','runtime_dependencies')}
            actual=[dict(active=dict(path='/final/capsule/numerics.py',sha256='same',bytes=1))]
            wrong=[dict(active=dict(path='/prior/capsule/numerics.py',sha256='same',bytes=1))]
            old=put(root/'old.json',dict(status='PASS',actual_imports_verified=True,historical_manifest=historical,active_sources=wrong,**pins))
            new=put(root/'new.json',dict(status='PASS',actual_imports_verified=True,historical_manifest=historical,active_sources=actual,**pins))
            with patch('wj_h23_continuation.active_source_bindings',return_value=actual):
                with self.assertRaises(AssertionError):active_execution_bindings(old)
                self.assertEqual(active_execution_bindings(new)['status'],'PASS')
            self.assertEqual(identity(root/'old.json'),old)

    def test_unchanged_scientific_xml_still_rejects_dt_model_and_functions(self):
        x=ET.fromstring('<NEKTAR><CONDITIONS><PARAMETERS><P>Time = 0</P><P>TimeStep = .0025</P><P>NumSteps = 8</P><P>RANSWallH = .15</P></PARAMETERS><SOLVERINFO><I PROPERTY="RANSWallModel" VALUE="Original"/></SOLVERINFO><FUNCTION NAME="InitialConditions"><E VAR="u" VALUE="log(y)"/></FUNCTION></CONDITIONS></NEKTAR>')
        original=protected_case(x)
        for text in ('TimeStep = .00125','RANSWallH = .2'):
            bad=copy.deepcopy(x);key=text.split('=')[0].strip()
            next(n for n in bad.findall('./CONDITIONS/PARAMETERS/P')if n.text.split('=')[0].strip()==key).text=text
            self.assertNotEqual(protected_case(bad),original)
        bad=copy.deepcopy(x);bad.find('./CONDITIONS/FUNCTION/E').set('VALUE','0');self.assertNotEqual(protected_case(bad),original)
        bad=copy.deepcopy(x);bad.find('./CONDITIONS/SOLVERINFO/I').set('VALUE','EP');self.assertNotEqual(protected_case(bad),original)

    def test_binding_hook_is_called_by_actual_coordinator_before_native_Popen(self):
        # Inspect the actual callback seam, not a separate permissive launcher.
        import inspect,study_wj_h23,study_wj_h24
        source=inspect.getsource(study_wj_h23.target)
        self.assertLess(source.index('if prelaunch is not None:prelaunch()'),source.index("native=subprocess.Popen"))
        with patch('study_wj_h24.frozen'),patch('study_wj_h24.check_bindings')as guard,patch('study_wj_h23.target')as target:
            study_wj_h24.target(Path('/fixture'),'M')
            target.call_args.kwargs['prelaunch']()
            self.assertEqual(guard.call_count,2)

if __name__=='__main__':unittest.main()
