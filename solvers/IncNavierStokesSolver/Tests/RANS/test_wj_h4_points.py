# SPDX-License-Identifier: MIT
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from numpy.testing import assert_allclose
from wj_h4_points import mapping,physical,owners,owner_point,field,evaluate,scan,pair_rows,aligned,digest
from validate_cls_epawf_m3 import basis


def fixture(zknots=(0.,.5,1.,2.)):
    H=.5;d0=.0875;top=20.;J=(top-d0-H)/(top-H);mesh=[]
    for lo,hi in zip(zknots[:-1],zknots[1:]):
        for x0,x1 in ((0.,1.),(1.,2.)):
            mesh.append(dict(geometry=len(mesh),lower=[x0,0,lo],upper=[x1,0,hi]))
    report=dict(parameters=[6.423173803526447e-5,H,.0372],geometry_signature=f'fixture:CompositeTunnelV1:{d0}:{top}',fields=[],offsets=[])
    state=[]
    def mean(x,r):
        u=2*r if r<=H else 2*H+2*J*(r-H) if r<=1 else 2*H+2*J*(1-H)+3.5*J*(r-1)
        return u+.02*abs(x-1)
    for b in range(6):
        report['offsets'].append(len(state));fm=[];local=0
        for el in mesh:
            inner=el['upper'][2]<=H
            if b==3 and inner or b==4 and not inner or b==5 and inner:continue
            p=1 if b>=4 else 3;n=(p+1)**2;g=np.linspace(-1,1,p+1);px=basis(p,g)[0];V=np.kron(px,px)
            xx,rr=np.meshgrid((el['lower'][0]+el['upper'][0])/2+(el['upper'][0]-el['lower'][0])/2*g,
                (el['lower'][2]+el['upper'][2])/2+(el['upper'][2]-el['lower'][2])/2*g)
            f=np.array([mean(x,r) if b==0 else 0. if b==1 or b>=4 else 1. for x,r in zip(xx.ravel(),rr.ravel())])
            coeff=np.linalg.solve(V,f);fm.append(dict(el,degree=p,ids=list(range(local,local+n)),sign=[1]*n))
            local+=n;state.extend(coeff)
        report['fields'].append(dict(mesh=fm))
    report['offsets'].append(len(state));return report,state


class OwnedPhysicalPoints(unittest.TestCase):
    def test_faces_intersection_and_adverse_high_owner(self):
        r,s=fixture();p=physical(r,[.5,0,1.]);o=owners(r,p);self.assertEqual(len(o),2)
        vals=[evaluate(r,s,p,v) for v in o]
        self.assertGreater(vals[0]['principal']['in_plane'],0)
        self.assertLess(vals[1]['principal']['in_plane'],0)
        self.assertEqual(len(owners(r,physical(r,[1.,0,1.25]))),2)
        self.assertEqual(len(owners(r,physical(r,[1.,0,1.]))),4)
        assert_allclose([v['G'][0][2] for v in vals],[2.,3.5],atol=2e-14,rtol=0)

    def test_roles_inward_offsets_and_no_extrapolation(self):
        r,s=fixture();p=physical(r,[.5,0,.5]);rows=[evaluate(r,s,p,o)for o in owners(r,p)]
        self.assertEqual({v['owner']['region']for v in rows},{'inner','outer'})
        self.assertIsNone(next(v for v in rows if v['owner']['region']=='inner')['epsilon'])
        self.assertIsNone(next(v for v in rows if v['owner']['region']=='inner')['principal'])
        q=physical(r,[.5,0,1.]);lower=q.copy();lower[2]-=1e-6;upper=q.copy();upper[2]+=1e-6
        self.assertEqual(len(owners(r,lower)),1);self.assertEqual(len(owners(r,upper)),1)
        el=next(e for e in r['fields'][0]['mesh']if e['geometry']==owners(r,lower)[0]['geometry'])
        with self.assertRaises(ValueError):owner_point(r,el,upper)
        forged=dict(owners(r,lower)[0],chart=[.5,0,1.1])
        with self.assertRaises(ValueError):field(r,s,0,forged)

    def test_nonmatching_mesh_pairs_by_physical_coordinate(self):
        a,sa=fixture();b,sb=fixture((0.,.5,2.));p=physical(a,[.5,0,1.])
        inv=dict(points=[dict(id=0,physical=p,tags=['fixture'])],sha256='same-common-inventory')
        x,y=scan(a,sa,inv),scan(b,sb,inv)
        self.assertEqual(x['status'],'FAIL');self.assertEqual(len(pair_rows(x,y,2,.01,.02,4,.005,.020000000000000004)),2)
        with self.assertRaises(ValueError):pair_rows(x,y,2,.01,.02,3,.005,.015)
        with self.assertRaises(ValueError):aligned(4,.005,.0200001)

    def test_chart_scaling_is_not_physical_height(self):
        r,s=fixture();p=physical(r,[.5,0,1.25]);o=owners(r,p)[0];u,g=field(r,s,0,o)
        h=1e-6;pl=p.copy();pu=p.copy();pl[2]-=h;pu[2]+=h
        el=next(e for e in r['fields'][0]['mesh']if e['geometry']==o['geometry'])
        ul=field(r,s,0,owner_point(r,el,pl))[0];uu=field(r,s,0,owner_point(r,el,pu))[0]
        assert_allclose(g[2],(uu-ul)/(2*h),atol=3e-9,rtol=0)
        self.assertNotEqual(p[2],o['chart'][2]);assert_allclose(g[2],3.5,atol=1e-13,rtol=0)

    def test_optional_physical_knots_preserve_defaults(self):
        from wj_h2_cases import case_xml
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);a=case_xml(d/'a',steps=0);b=case_xml(d/'b',steps=0,extra_physical_z=None)
            self.assertEqual(a.read_text().replace(str(d/'a'),'<case>'),b.read_text().replace(str(d/'b'),'<case>'))
            knots=[.8818157002137355,1.1947077409162703]
            p=case_xml(d/'extra',steps=0,extra_physical_z=knots)
            design=json.loads(p.with_name('case-design.json').read_text());H=design['completion']['H_w']
            zs=[.0875+H+(rho-H)*(20-.0875-H)/(20-H) for rho in design['chart_z'] if rho>=H]
            for z in knots:self.assertLess(min(abs(z-v)for v in zs),3e-14)
            self.assertEqual(len(design['chart_z']),26)
            for bad in ([float('nan')],[.2],[20.],knots[::-1],[.8,.9,1.]):
                with self.assertRaises(ValueError):case_xml(d/'invalid',steps=0,extra_physical_z=bad)

    def test_partial_failure_preserves_original_reason_and_counts(self):
        from wj_h4_study import attempt_counts,COUNTERS
        r=dict(failure='original principal/input exit',pimple_attempts=[dict(zip(COUNTERS,[1,1,0,0,0]))])
        v=attempt_counts(r);self.assertEqual(v['native_failure'],r['failure']);self.assertEqual(v['partial'],[[1,1,0,0,0]])
        self.assertTrue(v['no_extra_solves']);self.assertEqual(v['complete'],0)
        r['failure']='';self.assertFalse(attempt_counts(r)['success_counts_valid'])
        r['pimple_attempts'][0][COUNTERS[0]]=2;self.assertFalse(attempt_counts(r)['no_extra_solves'])

if __name__=='__main__':unittest.main()
