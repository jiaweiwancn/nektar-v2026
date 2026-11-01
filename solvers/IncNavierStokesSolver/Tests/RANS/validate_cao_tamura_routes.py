#!/usr/bin/env python3
"""Independent guards for native composite-route feasibility evidence."""
import csv
import math
from pathlib import Path
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

sys.dont_write_bytecode=True
from cao_tamura_experiment import config,generate,inlet_values,reference,wall_epsilon
from cao_tamura_routes import coverage,envelope,route_metadata
from generate_m4 import info


class RouteGuards(unittest.TestCase):
    def test_explicit_routes_and_declared_actual_mismatch(self):
        with tempfile.TemporaryDirectory() as d:
            for kind in ('EPTransition','UniformEP'):
                p=Path(d)/kind;generate(dict(model='cls_fixed',wall_route=kind),p)
                self.assertEqual(route_metadata(p/'case.xml')['kind'],kind)
            case=Path(d)/'UniformEP/case.xml';tree=ET.parse(case)
            info(tree.getroot(),'RANSExperimentalRoute','EPTransition');tree.write(case)
            with self.assertRaises(ValueError):route_metadata(case)
        with self.assertRaises(ValueError):config(dict(model='cls_original',wall_route='EPTransition'))

    def test_spatial_epsilon_trace_and_common_outer_inputs(self):
        for surface in ('S','R'):
            cs=[config(dict(surface=surface,model=m,wall_route='EPTransition',
                inlet_delta=.21,k_profile='MeasuredRms',epsilon_profile='ShearLength',inlet_wall_collar=True,
                inlet_reference_q=10/3,inlet_reference_kappa=math.sqrt(.48*1.3/(10/3))))
                for m in ('standard_fixed','cls_fixed')]
            rs=[reference(c) for c in cs]
            for x in (-1.,-.28,-.26,-.24,-.22,-.2,0.,1.2):
                for c,r in zip(cs,rs):
                    ep=r['wall_epsilon_near'];original=r['k_near']**1.5/(2.55*c['hm'])
                    expected=original+envelope(x,-.28,-.2)*(ep-original)
                    self.assertAlmostEqual(wall_epsilon(c,r,x),expected,places=13)
                    self.assertAlmostEqual(inlet_values(c['hm'],c,r,x)['epsilon'],expected,places=13)
                for z in (.022,.03,.04,.1,.21,.3,.8):
                    a,b=[inlet_values(z,c,r,x) for c,r in zip(cs,rs)]
                    for key in a:self.assertAlmostEqual(a[key],b[key],places=13)

    def test_zero_downstream_is_not_fallback_but_lost_upstream_is(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);generate(dict(model='cls_fixed'),d);case=d/'case.xml'
            (d/'case.rans.csv').write_text('time\n0.01\n0.02\n')
            rows=[]
            for t in (.01,.02):
                for x in (-.6,-.24,0.):
                    w=envelope(x,-.28,-.2)
                    rows.append(dict(phase='1',time=t,x=x,mask=w,requested_chi=w,chi=w,reason=0 if w else 1))
            def save():
                with (d/'case.terrain-points-r0.csv').open('w',newline='') as f:
                    writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
            save();result=coverage(case)
            self.assertEqual(result['fallback_native_traces'],0)
            self.assertEqual(result['accepted_times'],2)
            rows[0].update(chi=0,reason=2);save()
            with self.assertRaises(ValueError):coverage(case)
            rows[0].update(chi=1,reason=0);rows[:]=[r for r in rows if r['x']!=-.24];save()
            with self.assertRaises(ValueError):coverage(case)
            # Even complete regions at one time cannot cover a missing frame.
            rows[:]=[dict(phase='1',time=.01,x=x,mask=envelope(x,-.28,-.2),
                requested_chi=envelope(x,-.28,-.2),chi=envelope(x,-.28,-.2),reason=0)
                for x in (-.6,-.24,0.)]
            save()
            with self.assertRaises(ValueError):coverage(case)


if __name__=='__main__':unittest.main()
