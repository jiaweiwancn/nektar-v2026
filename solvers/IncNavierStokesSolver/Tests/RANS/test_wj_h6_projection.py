#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from numpy.testing import assert_allclose
from wj_h5_neutral import NeutralReference,mean_fit,SCALE_MM
from wj_h5_cases import mesh_knots,case_xml
from wj_h6_projection import project_field,constrained_segment,Functions,expression,reconstruct
from validate_cls_epawf_m3 import basis
from wj_h3_analysis import values


def fixture():
    H=.5;z=[0.,H,1.,2.];fields=[];offset=[0]
    for b in range(6):
        p=1 if b>=4 else 3;m=p+1;mesh=[];last=None;n=0
        for j,(lo,hi)in enumerate(zip(z[:-1],z[1:])):
            if (b==3 and hi<=H)or(b==4 and lo>=H)or(b==5 and hi<=H):continue
            ids=[]
            for ir in range(m):
                for ix in range(m):
                    if last is not None and ir==0:ids.append(last[m+ix])
                    else:ids.append(n);n+=1
            mesh.append(dict(geometry=j,degree=p,lower=[0,0,lo],upper=[8,0,hi],ids=ids,sign=[1]*(m*m)));last=ids
        fields.append(dict(size=n,mesh=mesh));offset.append(offset[-1]+n)
    s=[0.]*offset[-1]
    return dict(parameters=[.000064,.5,.1,.5,.5,.0088],ustar=.051,geometry_signature='fixture:CompositeTunnelV1:0.0875:20',
        fields=fields,offsets=offset,state=s,steps=0,time=0,implicit_callback_count=0,pimple_attempts=[],failure='',accepted_states=[s])


class ProjectionChecks(unittest.TestCase):
    def test_prescribed_midpoints_counts_and_functions(self):
        import xml.etree.ElementTree as ET
        with tempfile.TemporaryDirectory()as folder:
            for H in (5.,5.5,6.):
                c=NeutralReference(H/SCALE_MM,mean_fit()['z0_e']);a,b,z=mesh_knots(c);aa,bb,zz=mesh_knots(c,bisect_once=True)
                self.assertEqual((len(a)-1,len(b)-1),(8,17));self.assertEqual((len(aa)-1,len(bb)-1),(16,34))
                self.assertEqual(38*16+50*34+20,2328);self.assertEqual(zz[::2],z);self.assertEqual(zz.count(c.H),1)
                J=(c.top-c.d0-c.H)/(c.top-c.H)
                physical=lambda q:c.d0+q if q<=c.H else c.d0+c.H+(q-c.H)*J
                assert_allclose([physical(q)for q in zz[1::2]],[(physical(x)+physical(y))/2 for x,y in zip(z[:-1],z[1:])],atol=4e-15,rtol=0)
                old=case_xml(Path(folder)/str(H),c);new=case_xml(Path(folder)/(str(H)+'-bisect'),c,bisect_once=True)
                xml=lambda p:ET.parse(p).findall('./CONDITIONS/FUNCTION')
                self.assertEqual([ET.tostring(v)for v in xml(old)],[ET.tostring(v)for v in xml(new)])
                self.assertEqual(json.loads(new.with_name('case-design.json').read_text())['predicted_state_size'],2328)
                self.assertIn('NumSteps = 0',new.read_text());self.assertNotIn('RANSWallRestart',new.read_text())
        with self.assertRaises(ValueError):mesh_knots(c,bisect_once=2)

    def test_physical_mass_constants_and_representable_polynomial(self):
        r=fixture();H=.5;J=(20-.0875-H)/(20-H)
        r['fields'][0]['mesh'][0]['sign'][10]=-1  # an owned interior orientation
        physical_z=lambda rho:.0875+np.where(rho<=H,rho,H+(rho-H)*J)
        fn=lambda x,rho:1+.2*x+.3*physical_z(rho)**2-.01*x*physical_z(rho)
        c,meta=project_field(r,0,fn);self.assertLess(meta['backward_residual'],1e-15)
        state=np.zeros(len(r['state']));state[:len(c)]=c
        for el in r['fields'][0]['mesh']:
            x=np.array([0.,.7,4.,8.]);rho=np.linspace(el['lower'][2],el['upper'][2],4)
            v,g=values(r,state,0,el['geometry'],x,rho)
            assert_allclose(v,fn(x,rho),rtol=2e-14,atol=3e-14)
            assert_allclose(g[:,2],.6*physical_z(rho)-.01*x,rtol=2e-13,atol=5e-14)
        one,_=project_field(r,0,lambda x,z:np.full_like(x,3.));s=np.zeros_like(state);s[:len(one)]=one
        for el in r['fields'][0]['mesh']:
            assert_allclose(values(r,s,0,el['geometry'],[2.],[sum([el['lower'][2],el['upper'][2]])/2])[0],[3.],atol=8e-15,rtol=0)

    def test_endpoint_constrained_trace_and_expression_safety(self):
        known=np.array([.7,1.2,.3,-.4]);c=constrained_segment(lambda x:basis(3,x)[0]@known)
        assert_allclose(c,known,atol=2e-15,rtol=0)
        assert_allclose(expression('(y<1)*(y^2)+(y>=1)*2')([0,0],[.5,2]),[.25,2.],atol=0,rtol=0)
        with self.assertRaises(ValueError):expression('__import__("os")')

    def test_lifting_and_actual_K_trace_are_separate(self):
        r=fixture();H=r['parameters'][1];d=.0875;K=r['parameters'][5];A=(K/(81/(2*np.sqrt(143))))**1.5/.381917863546733
        import xml.etree.ElementTree as ET
        root=ET.Element('NEKTAR');co=ET.SubElement(root,'CONDITIONS')
        for family in ('InitialConditions','BoundaryValues'):
            for side in ('Inner','Outer'):
                f=ET.SubElement(co,'FUNCTION',NAME=family+side)
                for key,value in dict(u=f'(y-{d})^2',v='0',k=f'{K}*(1+.2*(1-(y-{d})/{H})^4)'if side=='Inner'else str(K),
                    epsilon=f'{A}/(y-{d})',p='-.002'if side=='Inner'else '0').items():ET.SubElement(f,'E',VAR=key,VALUE=value)
        with tempfile.TemporaryDirectory()as folder:
            p=Path(folder)/'functions.xml';ET.ElementTree(root).write(p)
            stages,meta=reconstruct(r,p)
            b,c,dv=[stages[k]for k in ('physical_L2','exterior_lift','dependent_epsilon_lift')]
            self.assertGreater(np.max(abs(c-b)),1e-10)
            a,z=r['offsets'][3:5];self.assertGreater(np.max(abs(dv[a:z]-c[a:z])),1e-10)
            assert_allclose(dv[:a],c[:a],atol=0,rtol=0);assert_allclose(dv[z:],c[z:],atol=0,rtol=0)
            for x in (0.,4.,8.):
                k=values(r,dv,2,1,[x],[H])[0][0];eps=values(r,dv,3,1,[x],[H])[0][0]
                self.assertAlmostEqual(eps,(k/(81/(2*np.sqrt(143))))**1.5/(.381917863546733*H),places=12)
            self.assertNotAlmostEqual(values(r,dv,2,1,[8.],[H])[0][0],K,places=10)
            r.update(state=dv.tolist(),accepted_states=[dv.tolist()]);_,meta=reconstruct(r,p);self.assertEqual(meta['status'],'PASS')
            r.update(failure='native startup exit',checkpoint='',accepted_states=[],physical_range_events=[dict(state=dv.tolist(),alpha=0,time=0)])
            _,meta=reconstruct(r,p);self.assertIn('uncommitted',meta['native_lineage'])
            r['checkpoint']='invented accepted checkpoint'
            with self.assertRaises(AssertionError):reconstruct(r,p)

    def test_executed_unresolved_is_not_not_run(self):
        from study_wj_h6 import weak_status
        self.assertEqual(weak_status(True,'UNRESOLVED'),'UNRESOLVED')
        self.assertEqual(weak_status(False,None),'NOT_RUN')
        self.assertEqual(weak_status(True,'PASS'),'PASS')
        with self.assertRaises(ValueError):weak_status(True,None)

    def test_two_distinct_reading_intervals(self):
        from wj_h5_neutral import KAPPA
        f=mean_fit();lo=[];hi=[]
        for row in f['points']:
            lo.append((row['r_mm']-row['reading_height_mm'])*np.exp(-KAPPA*(row['value']+row['reading_value'])/f['u_tau']))
            hi.append((row['r_mm']+row['reading_height_mm'])*np.exp(-KAPPA*(row['value']-row['reading_value'])/f['u_tau']))
        assert_allclose([max(lo),min(hi)],[.2261849783,.2309394239],atol=6e-11,rtol=0)
        assert_allclose(f['reading_bounds_mm']['combined'],[.2171963971,.2370377766],atol=6e-11,rtol=0)
        self.assertLess(max(lo),f['z0_e_mm']);self.assertGreater(min(hi),f['z0_e_mm'])


if __name__=='__main__':unittest.main()
