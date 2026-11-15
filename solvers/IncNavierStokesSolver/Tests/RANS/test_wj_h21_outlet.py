# SPDX-License-Identifier: MIT
import copy
import tempfile
from pathlib import Path
import unittest
import numpy as np
from wj_h21_outlet import (ROLES, POINT, OWNERS, TIMES, DT, scale, checked_owner,
                           linear, endpoint_window, comparison, complete_structure, geometry_compatible)
from wj_h12_phase import local, decompositions, cancellation, DIRECTION
from wj_h4_points import chart, owner_point
from wj_h15_observer import identity, verify
from run_wj_h21 import passive_intent


def fixture():
    r=dict(parameters=[.051*.005/3.97,.15,.1085450354253546],
           geometry_signature='fixture:CompositeTunnelV1:0.0875:20',offsets=[i*32 for i in range(6)])
    faces=[.24678758944336157,POINT[2],.2753627683300847]
    mesh=[]
    for j,(lo,hi)in enumerate(zip(faces[:-1],faces[1:])):
        mesh.append(dict(geometry=OWNERS[j][0],lower=chart(r,[12.5,0,lo]).tolist(),upper=chart(r,[15,0,hi]).tolist(),
            degree=3,ids=list(range(j*16,(j+1)*16)),sign=[-1 if i%3==0 else 1 for i in range(16)]))
    r['fields']=[dict(mesh=copy.deepcopy(mesh))for _ in range(6)]
    s=np.zeros(192)
    for block in range(6):
        for j,(lo,hi)in enumerate(zip(faces[:-1],faces[1:])):
            c=np.zeros((4,4));c[0,:2]=[block+25+3*lo,block+30+3*lo]
            c[1,:2]=[block+25+3*hi,block+30+3*hi]
            el=mesh[j];s[np.array(el['ids'])+r['offsets'][block]]=c.ravel()*el['sign']
    return r,s,[dict(physical=POINT.copy(),owner=owner_point(r,e,POINT))for e in mesh]


class OutletDiagnosis(unittest.TestCase):
    def test_exact_24_endpoint_roles_and_no_stage_substitution(self):
        self.assertEqual(len(ROLES),24)
        self.assertEqual(sorted({(t,z,n)for t,z,n,g,s in ROLES}),
                         sorted((t,z,round(t/dt))for t in TIMES for z,dt in DT.items()))
        self.assertEqual({(g,s)for t,z,n,g,s in ROLES},set(OWNERS))
        self.assertEqual({t for t,z,n,g,s in ROLES},{.005,.01,.015,.02})
        self.assertEqual(len({r for r in ROLES}),24)

    def test_native_signs_and_physical_jacobian_at_both_traces(self):
        r,s,owners=fixture()
        for j,o in enumerate(owners):
            v=linear(r,s,o,j)
            np.testing.assert_allclose(v['U'],[30+3*POINT[2],0,31+3*POINT[2]],atol=3e-13)
            np.testing.assert_allclose(v['G'],[[2,0,3],[0,0,0],[2,0,3]],atol=1e-12)
            self.assertAlmostEqual(v['pi'],35+3*POINT[2])

    def test_bounded_structure_requires_pinned_original_parameters(self):
        r,s,o=fixture();r.update(map_parameters=[0.]*12,quadrature_orders=[12,16])
        bounded={k:r[k]for k in ('fields','offsets','geometry_signature')}
        result=complete_structure(bounded,r)
        self.assertEqual(result['parameters'],r['parameters'])
        bad=copy.deepcopy(bounded);bad['fields'][0]['mesh'][0]['sign'][0]*=-1
        with self.assertRaises(AssertionError):complete_structure(bad,r)
        bad=copy.deepcopy(bounded);bad['parameters']=[0.,.15,.1]
        with self.assertRaises(AssertionError):complete_structure(bad,r)

    def test_dt_bearing_structure_identity_does_not_relax_spatial_equality(self):
        r,s,o=fixture();a=dict(fields=r['fields'],offsets=r['offsets'],geometry_signature=r['geometry_signature'],
                             physical_signature_sha256='own-dt',intent='one-run')
        b=copy.deepcopy(a);b.update(physical_signature_sha256='another-dt',intent='other-run')
        geometry_compatible(a,b)
        b['fields'][0]['mesh'][0]['degree']=2
        with self.assertRaises(AssertionError):geometry_compatible(a,b)

    def test_wrong_owner_snap_and_coordinate_bits_rejected(self):
        r,s,o=fixture()
        for mutation in ('geometry','normal_side','roundoff_snap','physical'):
            bad=copy.deepcopy(o[0])
            if mutation=='geometry':bad['owner']['geometry']=265
            elif mutation=='normal_side':bad['owner']['normal_side']='lower'
            elif mutation=='roundoff_snap':bad['owner']['roundoff_snap'][2]=1e-8
            else:bad['physical'][2]=float(np.nextafter(POINT[2],np.inf))
            with self.assertRaises((ValueError,AssertionError)):checked_owner(r,bad,0)
        bad=copy.deepcopy(o[0]);bad['physical'][1]=-0.
        with self.assertRaises(AssertionError):checked_owner(r,bad,0)

    def test_real_closure_tensor_normalization_and_gap_identities(self):
        G=np.array([[.003,0,.79],[0,0,0],[.001,0,-.002]])
        before=local(dict(G=G,k=.0088,epsilon=.0021))
        after=local(dict(G=G*1.002,k=.00881,epsilon=.00212))
        v=decompositions(before,after)
        self.assertEqual(v['status'],'PASS')
        self.assertGreater(np.linalg.norm(v['delta_B']['gradient_term']),0)
        self.assertGreater(np.linalg.norm(v['delta_B']['timescale_term']),0)
        self.assertLess(v['normalized_identity']['maximum_error'],1e-13)
        self.assertLess(v['optimizer_identity']['maximum_error'],1e-13)
        self.assertGreater(abs(v['normalization_term']),0)
        self.assertGreater(abs(v['dimensional_term']),0)
        for r in (before,after):
            self.assertAlmostEqual(r['a'],r['d']/r['nu_eff'])
            self.assertGreaterEqual(r['g'],-2e-11)
            self.assertLess(r['implicit_derivative_error'],2e-11)
        self.assertAlmostEqual(np.linalg.norm(DIRECTION),1.)
        self.assertAlmostEqual(np.trace(DIRECTION),0.)

    def test_terminal_nonzero_incoming_and_signed_cancellation(self):
        c0=dict(mhat=.004-3e-6);f0=dict(mhat=.004)
        c1=dict(mhat=.004+4e-6);f1=dict(mhat=.004+3.99e-6)
        v=endpoint_window(c0,c1,f0,f1,'mhat')
        self.assertAlmostEqual(float(v['incoming']),-3e-6)
        self.assertAlmostEqual(float(v['net_update']),3.01e-6)
        self.assertAlmostEqual(float(v['endpoint']),1e-8)
        self.assertGreater(float(v['net_cancellation']['cancellation']),.998)
        self.assertLess(v['identity']['maximum_error'],1e-17)
        # A componentwise sign reversal must not be replaced by a tensor norm.
        a={k:dict(G=np.array([x,-x]))for k,x in [('c0',2.),('f0',1.),('c1',1.1),('f1',1.2)]}
        w=endpoint_window(a['c0'],a['c1'],a['f0'],a['f1'],'G')
        np.testing.assert_allclose(w['endpoint'],[-.1,.1])
        self.assertTrue(np.all(w['net_cancellation']['cancellation']>.9))

    def test_floor_classifies_without_zeroing_signed_values(self):
        v=endpoint_window(dict(mhat=1e-12),dict(mhat=-2e-12),dict(mhat=0.),dict(mhat=0.),'mhat')
        self.assertEqual(float(v['endpoint']),-2e-12)
        self.assertFalse(v['net_cancellation']['resolved'])
        self.assertTrue(np.isnan(v['net_cancellation']['cancellation']))
        self.assertEqual(comparison(0.,1e-16,scale('pi'))['absolute_error'],1e-16)
        with self.assertRaises(AssertionError):comparison(0.,1e-6,1.)

    def test_hash_binding_rejects_changed_saved_input(self):
        with tempfile.TemporaryDirectory()as d:
            p=Path(d)/'packet.json';p.write_text('original');v=identity(p);verify(v)
            p.write_text('changed')
            with self.assertRaises(AssertionError):verify(v)

    def test_zero_native_mass_and_arbitrary_bridge_scope(self):
        for native in ('IncNavierStokesSolver','RANSWallResolutionProbe','WJKernelProbe','cmake','ninja','mpirun'):
            with self.assertRaises(AssertionError):passive_intent([native])
        for flag in ('--governing','--strict','--cp0','--restart','--build','--mass-solve'):
            with self.assertRaises(AssertionError):passive_intent(['python3','wj_h21_outlet.py',flag])
        with self.assertRaises(AssertionError):passive_intent(['python3','-c','print(1)'])
        self.assertIsNone(passive_intent(['python3','-B','wj_h21_outlet.py','--mode','science']))


if __name__=='__main__':unittest.main()
