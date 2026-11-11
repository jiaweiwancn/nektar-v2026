#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
import tempfile
import json
from pathlib import Path
import unittest
from unittest.mock import patch
from cls_awf_native_metrics import interface_metrics,gl4_differences,difference_localization
from watch_cls_awf_qualification_windows import eligible_endpoints,publish_stop,TERMINAL
from run_cls_epawf_picard import load
import watch_cls_awf_qualification_windows as observer


class CompleteMetrics(unittest.TestCase):
    def progress(self,status):
        return dict(status=status,dt=.000625,segments=[dict(end_step=n,failure=None,exit_code=0) for n in (7200,8000,8800,9600)])

    def test_terminal_cap_assesses_final_windows(self):
        for status in TERMINAL+('running',):
            self.assertEqual(eligible_endpoints(self.progress(status),[]),[9600])
        self.assertEqual(eligible_endpoints(self.progress('targets_completed'),[9600]),[])

    def test_terminal_main_runs_final_assessment_before_exit(self):
        for terminal in ('targets_completed','numerically_rejected'):
            with self.subTest(terminal=terminal),tempfile.TemporaryDirectory() as d:
                root=Path(d);campaign=root/'campaign';campaign.mkdir();(campaign/'progress.json').write_text('{}')
                policy=root/'policy.json';policy.write_text(json.dumps(dict(primary_physical_window=.5,required_consecutive_windows=3)))
                probe=root/'probe';probe.write_text('immutable probe fixture')
                progress=self.progress(terminal);progress.update(requested_targets=[9600],accepted_steps=9600,accepted_time=6.,checkpoint='saved')
                commands=[]
                def run(cmd,**kw):
                    commands.append(cmd);out=Path(cmd[cmd.index('--output')+1]);out.mkdir(parents=True)
                    q=dict(cases=dict(Refined05=dict(selected_steps=[7200,8000,8800,9600]))) if len(commands)==1 else dict(passing_windows=3)
                    (out/'summary.json').write_text(json.dumps(q))
                    return type('Result',(),dict(returncode=0))()
                argv=['observer','--campaign',str(campaign),'--policy',str(policy),'--probe',str(probe),'--output',str(root/'result')]
                with patch('sys.argv',argv),patch.object(observer,'load_progress',return_value=progress),patch.object(observer.subprocess,'run',side_effect=run),patch.object(observer,'terminal_readback',return_value=dict(status='verified')) as final:observer.main()
                q=load(root/'result/progress.json')
                self.assertEqual(q['status'],terminal);self.assertEqual(q['assessed_steps'],[9600])
                self.assertTrue(q['final_assessment_complete']);self.assertEqual(len(commands),2)
                final.assert_called_once();self.assertEqual(q['final_retained_state']['status'],'verified')
                self.assertFalse((campaign/'stop-after-block.json').exists())

    def test_rejected_attempt_not_accepted_endpoint(self):
        q=self.progress('numerically_rejected');q['segments'][-1]['failure']='covariance'
        self.assertEqual(eligible_endpoints(q,[]),[])

    def test_existing_request_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'stop.json';old={'owner':'other','reason':'retain me'}
            self.assertEqual(publish_stop(p,old)['status'],'published');before=p.read_bytes()
            r=publish_stop(p,{'owner':'new'})
            self.assertEqual(r['status'],'existing_stop_preserved');self.assertEqual(p.read_bytes(),before)
            self.assertEqual(load(p),old)

    def test_interface_index_names(self):
        d=interface_metrics(list(range(8)))
        self.assertEqual(d['normal_traction_jump_RMS'],2)
        self.assertEqual(d['full_vector_traction_jump_RMS'],3)
        self.assertEqual(d['maximum_absolute_sampled_normal_jump'],4)
        self.assertEqual(d['k_flux_jump_RMS'],5)

    def test_independent_norm_keeps_relative_pressure_offset(self):
        import math
        from numpy.polynomial.legendre import leggauss
        x,w=leggauss(8);z,v=leggauss(16);knots=[0.,.0025,.005,.02,.04,.5,1.,2.,4.];xknots=[-4.,-1.,2.]
        rows={k:[] for k in ('u','v','k','epsilon','pressure')}
        for a,b in zip(xknots,xknots[1:]):
            for t in x:
                xx=a+(t+1)*(b-a)/2
                for c,d in zip(knots,knots[1:]):
                    for r in z:
                        zz=c+(r+1)*(d-c)/2
                        rows['u'].append(xx**3+2*zz**2);rows['v'].append(xx*zz);rows['k'].append(zz**3)
                        rows['pressure'].append(1. if zz<.5 else -1.)
                        if zz>.5:rows['epsilon'].append(zz**3)
        base=dict(H=.5,common_reporting_knots=knots,common_samples=rows)
        zero=dict(H=.5,common_reporting_knots=knots,common_samples={k:[0.]*len(v) for k,v in rows.items()})
        result=gl4_differences(base,zero,xknots)
        self.assertAlmostEqual(result['dimensional']['pressure'],math.sqrt(10.5),places=12)
        self.assertAlmostEqual(result['dimensional']['k'],math.sqrt(6*4**7/7),places=10)
        self.assertAlmostEqual(result['physical_areas']['epsilon'],21.,places=12)
        local=difference_localization(base,zero,xknots)
        self.assertAlmostEqual(local['strong_divergence_difference_L2_inner_outer'][0],math.sqrt(.5*1564.8),places=10)
        self.assertAlmostEqual(local['strong_divergence_difference_L2_inner_outer'][1],math.sqrt(3.5*1564.8),places=10)


if __name__=='__main__':unittest.main()
