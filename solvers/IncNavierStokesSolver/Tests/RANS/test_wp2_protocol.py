# SPDX-License-Identifier: MIT
"""Physical comparison conventions and fail-fast WP2 gate contracts."""
import copy, unittest,tempfile,json,sys
from pathlib import Path
from unittest.mock import patch
from wp2_qualification import common_differences,quarter_budget_differences,physical_gates


class PhysicalComparison(unittest.TestCase):
    def test_adjacent_dt_signature_exclusion_is_narrow_and_checked(self):
        from analyze_wp2_global_recovery_dt import comparison_signature
        left=dict(signature='same|TIMESTEP=0.01|roughness=0.002',dt=.01)
        right=dict(signature='same|TIMESTEP=0.0050000000000000001|roughness=0.002',dt=.005)
        self.assertEqual(comparison_signature(left),comparison_signature(right))
        wrong=dict(right,signature=right['signature'].replace('roughness=0.002','roughness=0.003'))
        self.assertNotEqual(comparison_signature(left),comparison_signature(wrong))
        with self.assertRaisesRegex(ValueError,'signature timestep'):
            comparison_signature(dict(right,dt=.01))
    def test_window_envelope_detects_internal_excursion(self):
        import numpy as np
        from audit_wp2_window_envelope import range_checks
        sample=self.samples();v=np.asarray(sample['values']);R=np.zeros((2,3,3))
        lo=dict(values=v.copy(),R=R.copy(),T=R.copy());hi=copy.deepcopy(lo)
        lo['values'][:,0]-=.001;hi['values'][:,0]+=.001
        self.assertTrue(quarter_budget_differences(sample,sample,.1,.02,.0065)['passed'])
        check=range_checks(lo,hi,dict(values=v),sample['columns'],np.asarray(sample['weights']),np.asarray(sample['regions']),.02)
        self.assertFalse(check['passed'])
        self.assertFalse(check['principal']['gates']['0_U']['passed'])
    def test_long_recovery_time_bound_preserves_actual_additions(self):
        from wp2_qualification import physical_time_tolerance
        import math
        for dt,n in [(.02,60),(.01,120),(.01,450),(.005,900)]:
            actual=0.
            for _ in range(n):actual+=dt
            tolerance=physical_time_tolerance(n,dt,1.1942387693600194)
            self.assertLessEqual(abs(actual-n*dt),tolerance)
            self.assertGreater(abs(actual+1e-7-n*dt),tolerance)
        self.assertEqual(physical_time_tolerance(120,.01,1.1942387693600194),64*math.ulp(1.1942387693600194))
    def samples(self):
        r=[1.,0.,.02,.01,0.,-.013,.002,.005,.003,.001,0.,2.,0.,0.]
        return dict(columns='U,W,k,epsilon,pi,mechanical_pressure,nut,total_shear,total_normal_stress,production,U_x,U_z,W_x,W_z',
                    values=[r[:],r[:]],reference_values=[r[:],r[:]],weights=[.2,.8],regions=[0,1])
    def test_absolute_regional_pressure_and_complete_outer_epsilon(self):
        a=self.samples();b=copy.deepcopy(a);b['values'][0][3]+=100
        for row in b['values']:row[4]+=.003
        d=common_differences(a,b,.1,.02,.0065)
        self.assertEqual(d['epsilon']['dimensional_L2'],0.)
        self.assertAlmostEqual(d['pi']['dimensional_L2'],.003)
        self.assertFalse(quarter_budget_differences(a,b,.1,.02,.0065)['passed'])
    def test_nonzero_small_difference_and_full_gradient(self):
        a=self.samples();b=copy.deepcopy(a);b['values'][1][0]+=1e-5
        self.assertTrue(quarter_budget_differences(a,b,.1,.02,.0065)['passed'])
        b['values'][0][12]=.1
        self.assertFalse(quarter_budget_differences(a,b,.1,.02,.0065)['passed'])
    def test_recovery_interface_accuracy_is_recorded_until_final_qualification(self):
        import wp2_qualification as q
        d=dict(wp2_per_state_monitoring=[dict(time=.05,physical_checks={})],accepted_states=[[]],
            state_audits=[dict(accepted_awf_hard_admissible=[1],accepted_awf_covariance=[.005,.006],interface=[0,0,.00015,.00015,0,0])],pimple_attempts=[])
        with patch.object(q,'physical_gates',side_effect=lambda *args,**kwargs:dict(passed=True,gates={})):
            with self.assertRaisesRegex(ValueError,'native interface budget'):q.verify_report(d,.1,.02)
            result=q.verify_report(d,.1,.02,allow_recovery=True)
            self.assertFalse(result[-1]['passed']);self.assertAlmostEqual(result[-1]['gates']['interface_normal_full_RMS']['value'],.015)
            d['state_audits'][0]['accepted_awf_covariance'][0]=-.001
            with self.assertRaisesRegex(ValueError,'covariance margin'):q.verify_report(d,.1,.02,allow_recovery=True)

    def test_changed_grid_rejected(self):
        a=self.samples();b=copy.deepcopy(a);b['weights'][0]*=2
        with self.assertRaises(Exception):common_differences(a,b,.1,.02,.0065)


class ContinuationOwnership(unittest.TestCase):
    def test_confirmation_readback_requires_checked_terminal_snapshot(self):
        from complete_wp2_confirmation_mesh_readback import qualified_terminal
        from run_cls_epawf_picard import sha
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            self.assertIsNone(qualified_terminal(root,dict(status='physical')))
            with self.assertRaisesRegex(ValueError,'supervisor stopped'):
                qualified_terminal(root,dict(status='terminal_stop'))
            path=root/'confirmation-recovery-terminal-progress.json'
            path.write_text(json.dumps(dict(status='completed_bounded_trajectory',accepted_steps=120,actual_time=1.2,Treq=1.194)))
            s=dict(status='assessment',terminal_progress_sha256=sha(path),jobs=dict(physical=dict(exit_code=0)))
            self.assertEqual(qualified_terminal(root,s),path)
            s['jobs']['physical']['exit_code']=1
            with self.assertRaisesRegex(ValueError,'physical job did not succeed'):qualified_terminal(root,s)
            s['jobs']['physical']['exit_code']=0;path.write_text('{}')
            with self.assertRaisesRegex(ValueError,'snapshot changed'):qualified_terminal(root,s)
    def test_cumulative_resource_extension_preserves_origin_and_acceptance(self):
        import hashlib
        from wp2_budget import wall_limit
        from run_cls_epawf_picard import sha
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);(root/'origin.json').write_text('{}')
            self.assertEqual(wall_limit(root),43200)
            p=root/'source/docs/plans/cls-cao-tamura-budget-wp2-evidence/policy.json';p.parent.mkdir(parents=True)
            policy=dict(acceptance=dict(covariance=1e-4));p.write_text(json.dumps(policy))
            amendment=dict(version='WP2CumulativeResourceAmendmentV1',original_wall_seconds=43200,cumulative_wall_seconds=46800,
                origin_sha256=sha(root/'origin.json'),unchanged_acceptance_sha256=hashlib.sha256(json.dumps(policy['acceptance'],sort_keys=True,separators=(',',':')).encode()).hexdigest())
            (root/'resource-extension-13h.json').write_text(json.dumps(amendment));self.assertEqual(wall_limit(root),46800)
            extension=dict(version='WP2CumulativeResourceAmendmentV2',original_wall_seconds=43200,
                previous_cumulative_wall_seconds=46800,cumulative_wall_seconds=57600,additional_physical_attempt_limit=180,
                predecessor_sha256=sha(root/'resource-extension-13h.json'),origin_sha256=sha(root/'origin.json'),
                unchanged_acceptance_sha256=amendment['unchanged_acceptance_sha256'])
            (root/'resource-extension-16h.json').write_text(json.dumps(extension));self.assertEqual(wall_limit(root),57600)
            extension['additional_physical_attempt_limit']=181
            (root/'resource-extension-16h.json').write_text(json.dumps(extension))
            with self.assertRaisesRegex(ValueError,'unsupported WP2 global-control amendment'):wall_limit(root)
            extension['additional_physical_attempt_limit']=180
            (root/'resource-extension-16h.json').write_text(json.dumps(extension))
            confirmation=dict(version='WP2CumulativeResourceAmendmentV3',original_wall_seconds=43200,
                previous_cumulative_wall_seconds=57600,cumulative_wall_seconds=86400,additional_physical_attempt_limit=120,
                reserved_attempts_before=1646,predecessor_sha256=sha(root/'resource-extension-16h.json'),
                origin_sha256=sha(root/'origin.json'),unchanged_acceptance_sha256=amendment['unchanged_acceptance_sha256'])
            (root/'resource-extension-24h.json').write_text(json.dumps(confirmation));self.assertEqual(wall_limit(root),86400)
            confirmation['additional_physical_attempt_limit']=121
            (root/'resource-extension-24h.json').write_text(json.dumps(confirmation))
            with self.assertRaisesRegex(ValueError,'unsupported WP2 confirmation-recovery amendment'):wall_limit(root)
            confirmation['additional_physical_attempt_limit']=120
            (root/'resource-extension-24h.json').write_text(json.dumps(confirmation))
            from wp2_budget import pressure_scope
            scoped=dict(version='WP2PressureControlScopeV1',cumulative_wall_seconds=86400,reserved_attempts_before=1766,
                additional_physical_attempt_limit=742,spatial_case_limit=40,physical_attempt_limit_unchanged=6000,
                predecessor_sha256=sha(root/'resource-extension-24h.json'),origin_sha256=sha(root/'origin.json'),
                unchanged_acceptance_sha256=amendment['unchanged_acceptance_sha256'])
            (root/'pressure-control-scope.json').write_text(json.dumps(scoped))
            final=dict(version='WP2CumulativeResourceAmendmentV4',original_wall_seconds=43200,previous_cumulative_wall_seconds=86400,
                cumulative_wall_seconds=100800,reserved_attempts_before=1766,additional_physical_attempt_limit=1102,
                predecessor_sha256=sha(root/'resource-extension-24h.json'),origin_sha256=sha(root/'origin.json'),
                unchanged_acceptance_sha256=amendment['unchanged_acceptance_sha256'],pressure_scope_sha256=sha(root/'pressure-control-scope.json'))
            (root/'resource-extension-pressure-28h.json').write_text(json.dumps(final))
            self.assertEqual(wall_limit(root),100800);self.assertEqual(pressure_scope(root)['additional_physical_attempt_limit'],1102)
            final['additional_physical_attempt_limit']=1103
            (root/'resource-extension-pressure-28h.json').write_text(json.dumps(final))
            with self.assertRaisesRegex(ValueError,'unsupported WP2 pressure resource amendment'):wall_limit(root)
            final['additional_physical_attempt_limit']=1102
            (root/'resource-extension-pressure-28h.json').write_text(json.dumps(final))
            (root/'origin.json').write_text('{"reset":true}')
            with self.assertRaisesRegex(ValueError,'clock origin changed'):wall_limit(root)
    def test_pressure_gate_positive_bookkeeping_and_requested_transaction_dt(self):
        import qualify_wp2_pressure_control as driver
        from continue_cls_epawf_pimple import atomic
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);solver=root/'solver';solver.write_text('fixture identity')
            atomic(root/'origin.json',dict(started_epoch=__import__('time').time()))
            for batch in ['pressure-cold-primary','pressure-cold-confirmation']:
                folder=root/batch;folder.mkdir();atomic(folder/'progress.json',dict(status='completed_native_batch'))
                for model in ['CLS','standard']:
                    out=folder/(model+'-native');out.mkdir()
                    atomic(out/'result.json',dict(boundary=dict(passed=True),spatial=dict(analytical_free_block_norms=[0]*4)))
            requested=[]
            class Native:
                def __init__(self,*args):
                    self.output=args[4];self.output.mkdir();self.result=dict(results={})
                def save(self):atomic(self.output/'progress.json',self.result)
                def native(self,label,n,dt):
                    folder=self.output/label;folder.mkdir();report=dict(accepted_states=[[]],pimple_attempts=[dict(native_full_G_action_defect=0)],wp2_common_samples={})
                    atomic(folder/'report.json',report);(folder/'timing.jsonl').write_text(json.dumps(dict(name='native.step',wall_seconds=3.)))
                    self.result['results'][label]=dict(report=str(folder/'report.json'),usage='10 9 1 1000',seconds=10.)
                    return report
                def fixture(self,dt):requested.append(dt);self.save()
            class Geometry:
                def __init__(self,*args):pass
                def check_native(self,*args):return dict(passed=True)
                def tuple(self,*args):return {}
                def volume_audit(self,*args,**kwargs):return dict(passed=True)
            args=['pressure','--root',str(root),'--solver',str(solver),'--probe','/unused','--h1','/unused','--preflight','/unused','--model','standard']
            with patch.object(sys,'argv',args),patch.object(driver,'Run',Native),patch.object(driver,'Rows',Geometry),patch.object(driver,'charge'),patch.object(driver,'verify_report'),patch.object(driver,'physical_gates',return_value=dict(passed=True)),patch.object(driver,'quarter_budget_differences',return_value=dict(passed=True)):
                driver.main()
            d=json.loads((root/'pressure-initial-gates-standard.json').read_text())
            self.assertEqual(d['status'],'passed_initial_pressure_gates');self.assertEqual(requested,[.02])
            self.assertEqual([e['name'] for e in d['events'] if e['event_type']=='trial_passed'],['coarse','fine'])

    def test_readback_identity_resolves_actual_overlay_dependency(self):
        import wp2_qualification
        from analyze_wp2_approach import source_identity
        from run_cls_epawf_picard import sha
        with tempfile.TemporaryDirectory() as temp:
            entry=Path(temp)/'isolated.py';entry.write_text('# isolated analysis entry')
            self.assertFalse(entry.with_name('wp2_qualification.py').exists())
            identity=source_identity(['isolated.py','wp2_qualification.py'],entry)
            self.assertEqual(identity['isolated.py'],sha(entry))
            self.assertEqual(identity['wp2_qualification.py'],sha(Path(wp2_qualification.__file__)))

    @unittest.skipUnless(sys.platform.startswith('linux'),'Linux process-group ownership fixture')
    def test_native_wrapper_timeout_leaves_no_executing_child(self):
        import os,subprocess,time
        from wp2_native import execute_owned
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);pidfile=root/'child.pid'
            script="import subprocess,sys,time;from pathlib import Path;p=subprocess.Popen([sys.executable,'-c','import time;time.sleep(20)']);Path("+repr(str(pidfile))+").write_text(str(p.pid));time.sleep(20)"
            with (root/'log').open('w') as log:
                with self.assertRaises(subprocess.TimeoutExpired):
                    execute_owned([sys.executable,'-c',script],root,dict(os.environ),log,.5)
            self.assertTrue(pidfile.exists());pid=int(pidfile.read_text());status=Path('/proc')/str(pid)/'status'
            for _ in range(20):
                if not status.exists() or any(line.startswith('State:') and 'Z' in line for line in status.read_text().splitlines()):break
                time.sleep(.03)
            else:self.fail('native-like child survived the wrapper timeout')
            with (root/'success').open('w') as log:
                self.assertEqual(execute_owned([sys.executable,'-c','pass'],root,dict(os.environ),log,5).returncode,0)

    def test_fixture_honors_requested_production_timestep(self):
        from wp2_native import Run
        x=Run.__new__(Run);calls=[]
        # Stop at the first actual invocation: this is a dataflow test, not a
        # mocked numerical restart claim.
        def native(*args,**kwargs):
            calls.append((args,kwargs));raise RuntimeError('captured native dt')
        x.native=native;x.result={};x.save=lambda:None
        with self.assertRaisesRegex(RuntimeError,'captured native dt'):x.fixture(.02)
        self.assertEqual(calls[0][1]['dt'],.02)

    def test_new_pressure_assessment_does_not_mix_legacy_prefixes(self):
        from analyze_wp2_approach import paths
        root=Path('/immutable');progress=dict(case_definition=dict(pressure_space='inner-p2-v1'),observations=[dict(report='/immutable/new/report.json')])
        self.assertEqual(paths(root,'standard','coarse',progress),[Path('/immutable/new/report.json')])
        progress['prefix_observations']=[dict(report='/immutable/accepted-prefix/report.json')]
        self.assertEqual(paths(root,'standard','fine',progress),[Path('/immutable/accepted-prefix/report.json'),Path('/immutable/new/report.json')])

    def test_execution_index_includes_nested_legacy_and_deduplicates(self):
        from wp2_execution_index import executions
        from run_cls_epawf_picard import sha
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);run=root/'nested'/'run';run.mkdir(parents=True)
            report=run/'case.evolution.json';report.write_text(json.dumps(dict(implicit_callback_count=2,step_audits=[{}])))
            (run/'case.xml').write_text('<NEKTAR><CONDITIONS><PARAMETERS><P>NumSteps = 2</P></PARAMETERS></CONDITIONS></NEKTAR>')
            entry=dict(report=str(report),report_sha256=sha(report))
            (run/'progress.json').write_text(json.dumps(dict(results=dict(old=entry))))
            mirror=root/'mirror';mirror.mkdir()
            (mirror/'progress.json').write_text(json.dumps(dict(results=dict(reused=entry))))
            found=executions(root)
            self.assertEqual(len(found),1)
            self.assertEqual(found[str(report)]['observed_physical_attempts'],2)
            self.assertEqual(found[str(report)]['accepted_step_executions'],1)
            self.assertEqual(found[str(report)]['requested_attempt_work'],2)
            entry['observed_physical_attempts']=1
            (run/'progress.json').write_text(json.dumps(dict(results=dict(wrong=entry))))
            with self.assertRaisesRegex(ValueError,'execution count differs'):executions(root)
    def test_retained_state_audit_does_not_use_terminal_endpoint_cache(self):
        import wp2_common_quadrature as quadrature
        from wp2_all_state import Rows
        import numpy as np
        captured=[]
        def stop_at_quadrature(report,order):
            captured.append(report['physical_samples']['wall_K'])
            raise RuntimeError('quadrature input captured')
        rows=Rows.__new__(Rows)
        state_K=np.array([.017,.018])
        report=dict(wp2_common_samples={},physical_samples=dict(wall_K=[99.,99.]))
        with patch.object(quadrature,'resample',stop_at_quadrature):
            with self.assertRaisesRegex(RuntimeError,'quadrature input captured'):
                rows.volume_audit(report,dict(values=np.zeros((2,14)),endpoint_K=state_K))
        self.assertTrue(np.array_equal(captured[0],state_K))
        self.assertEqual(report['physical_samples']['wall_K'],[99.,99.])
    def test_terminal_recheck_rejects_unfinished_transaction_before_probe(self):
        import run_wp2_final_checks as driver
        with tempfile.TemporaryDirectory() as temp:
            r=Path(temp);p=r/'transit-transaction-CLS';p.mkdir()
            (p/'progress.json').write_text(json.dumps(dict(status='running')))
            args=['wp2','--root',str(r),'--solver','/unused','--probe','/unused','--h1','/unused',
                  '--preflight','/unused','--terminal-only','--record','terminal-recheck.json']
            with patch.object(sys,'argv',args),patch.object(driver,'run_probe') as native:
                with self.assertRaisesRegex(ValueError,'existing transaction incomplete'):driver.main()
            native.assert_not_called()
            self.assertEqual(json.loads((r/'terminal-recheck.json').read_text())['status'],'terminal_stop')
    def test_completion_requires_observed_physical_horizon_and_pass(self):
        from run_wp2_completion_waves import completed
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'progress.json'
            p.write_text(json.dumps(dict(status='terminal_stop',failure='native rejection')))
            with self.assertRaises(Exception):completed(p)
            p.write_text(json.dumps(dict(status='running')));self.assertFalse(completed(p))
            d=dict(status='completed_bounded_trajectory',actual_time=.1,final_physical_budgets_passed=True)
            p.write_text(json.dumps(d))
            with self.assertRaises(Exception):completed(p)
            d.update(actual_time=1.2,final_physical_budgets_passed=False);p.write_text(json.dumps(d))
            with self.assertRaises(Exception):completed(p)
            d['final_physical_budgets_passed']=True;p.write_text(json.dumps(d));self.assertTrue(completed(p))
    def test_changed_runtime_and_case_do_not_resume(self):
        from wp2_native import Run
        from run_cls_epawf_picard import sha
        with tempfile.TemporaryDirectory() as temp:
            r=Path(temp);s=r/'solver';s.write_text('immutable fixture');pre=r/'pre.json';pre.write_text('{}');out=r/'run';out.mkdir()
            d=dict(status='completed_bounded_trajectory',model='standard',solver_sha256=sha(s),case_definition=dict(level=4,order=16,projection='RadialEnergyV1'))
            (out/'progress.json').write_text(json.dumps(d))
            Run(r,s,r,pre,out,'standard',4,16,'RadialEnergyV1',True)
            with self.assertRaises(Exception):Run(r,s,r,pre,out,'standard',5,16,'RadialEnergyV1',True)
            unchanged=sha(out/'progress.json')
            with self.assertRaises(Exception):Run(r,s,r,pre,out,'standard',4,16,'RadialEnergyV1',True,'inner-p2-v1')
            self.assertEqual(sha(out/'progress.json'),unchanged)
            s.write_text('different runtime')
            with self.assertRaises(Exception):Run(r,s,r,pre,out,'standard',4,16,'RadialEnergyV1',True)
    def test_unexpected_failure_stops_all_later_blocks(self):
        import run_wp2_preservation as driver
        made=[]
        class Failed:
            def __init__(self,*args):self.result={};self.calls=0;made.append(self)
            def save(self):pass
            def native(self,*args,**kwargs):self.calls+=1;raise RuntimeError('unexpected native rejection')
        args=['wp2','--root','/unused','--solver','/unused','--h1','/unused','--preflight','/unused','--output','/unused','--model','standard','--mode','preserve','--dt','.02']
        with patch.object(driver,'Run',Failed),patch.object(sys,'argv',args):
            with self.assertRaisesRegex(RuntimeError,'unexpected native rejection'):driver.main()
        self.assertEqual(made[0].calls,1)
        self.assertEqual(made[0].result['status'],'terminal_stop')


if __name__=='__main__':unittest.main()
