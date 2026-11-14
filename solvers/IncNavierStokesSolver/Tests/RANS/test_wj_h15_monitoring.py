# SPDX-License-Identifier: MIT
import copy
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET
import numpy as np
from wj_h15_observer import PersistentRows,atomic,bits,identity,verify,pkey
from wj_h7_analysis import Rows
from wj_h4_points import physical,owners,aligned,vector_hash
from study_wj_h14 import stored_reference
from study_wj_h15 import xml_case

def tensor_fixture():
    """Coefficient/owner fixture at fixed6mm; no PDE run or reference seed."""
    c=stored_reference();mesh=[]
    for a,b in ((0.,c.H),(c.H,c.top)):
        for x,y in ((0.,4.),(4.,8.)):mesh.append(dict(geometry=len(mesh),lower=[x,0.,a],upper=[y,0.,b]))
    r=dict(fields=[],offsets=[],parameters=[c.nu,c.H,c.h],quadrature_orders=[12,16],geometry_signature=f'Plane:CompositeTunnelV1:{c.d0}:{c.top}',map_parameters=[0.]*14,signature='synthetic-fixed6mm-no-PDE')
    state=[]
    for b in range(6):
        r['offsets'].append(len(state));fm=[];local=0;p=1 if b>=4 else 3
        for e in mesh:
            inner=e['upper'][2]<=c.H
            if b==3and inner or b==4and not inner or b==5and inner:continue
            n=(p+1)**2;v=np.zeros((p+1,p+1));sign=np.ones(n);sign[2::5]=-1
            for j,rho in enumerate((e['lower'][2],e['upper'][2])):
                value=.01*rho if b==0 else c.K if b==2 else .001 if b==3 else 0.
                v[j,:2]=value
            state.extend(v.ravel()*sign);fm.append(dict(e,degree=p,ids=list(range(local,local+n)),sign=sign.tolist()));local+=n
        r['fields'].append(dict(size=local,mesh=fm))
    r['offsets'].append(len(state));return r,state,c

class PersistentMonitoring(unittest.TestCase):
    def test_atomic_pin_and_signed_zero(self):
        with tempfile.TemporaryDirectory()as d:
            p=Path(d)/'receipt.json';atomic(p,dict(value=-0.));self.assertEqual(verify(identity(p)),p)
            with self.assertRaises(FileExistsError):atomic(p,dict(value=1.))
            self.assertNotEqual(bits(-0.),bits(0.))
            p.write_text('changed')
            with self.assertRaises(AssertionError):verify(dict(path=str(p),bytes=1,sha256='old'))
    def test_integer_time_and_identity_are_not_nearby_step_aliases(self):
        aligned(3,.0025,.0075+2*np.finfo(float).eps)
        with self.assertRaises(ValueError):aligned(4,.0025,.0075)
        self.assertNotEqual(vector_hash([0.,1.]),vector_hash([-0.,1.]))
    def test_cache_rows_match_original_full_and_both_H_owners(self):
        r,s,c=tensor_fixture();points=[physical(r,[x,0,z])for x,z in ((0.,0.),(4.,c.H),(8.,.3))]
        old=Rows(r,points,c);cache=PersistentRows(r,points,s,dict(qh0=vector_hash(s)))
        new,initial,_=cache.view(points)
        self.assertEqual(new.rows,old.rows);self.assertEqual(new.weights.tobytes(),old.weights.tobytes())
        for b in range(6):
            self.assertEqual(new.maps[b][0].tobytes(),old.maps[b][0].tobytes())
            for a,z in zip(new.maps[b][1],old.maps[b][1]):self.assertEqual(a.tobytes(),z.tobytes())
        self.assertEqual(len(owners(r,points[1])),4)
        before=old.evaluate(s,principal=False);after=new.evaluate(s,principal=False)
        for k in ('U','G','k','epsilon','R','K','rc'):self.assertEqual(before[k].tobytes(),after[k].tobytes());self.assertEqual(before[k].tobytes(),initial[k].tobytes())
    def test_variable_eviction_does_not_remove_required_points(self):
        r,s,c=tensor_fixture();base=[physical(r,[0.,0,.1]),physical(r,[8.,0,.3])];cache=PersistentRows(r,base,s,{},variable_limit=2)
        points=base+[physical(r,[x,0,.6])for x in (1.,2.,3.,5.,6.)]
        view,_,work=cache.view(points);self.assertEqual(work['variable_points'],2);self.assertGreater(work['evictions'],0)
        self.assertEqual({pkey(p)for p,o in view.rows[:view.point_count]},{pkey(p)for p in points})
        again,_,_=cache.view(points)
        for b in range(6):
            for a,z in zip(view.maps[b][1],again.maps[b][1]):self.assertEqual(a.tobytes(),z.tobytes())
    def test_state_dependencies_and_failure_are_not_cached(self):
        r,s,c=tensor_fixture();cache=PersistentRows(r,[physical(r,[4.,0,c.H])],s,{})
        rows,_,_=cache.view([physical(r,[4.,0,c.H])]);a=rows.evaluate(s,principal=False);changed=np.array(s)
        changed[r['offsets'][2]:r['offsets'][3]]*=1.01;b=rows.evaluate(changed,principal=False)
        self.assertGreater(np.max(abs(a['K']-b['K'])),0);self.assertGreater(np.max(abs(a['rc']-b['rc'])),0)
        changed[r['offsets'][2]:r['offsets'][3]]=-abs(changed[r['offsets'][2]:r['offsets'][3]])
        with self.assertRaises(ValueError):rows.evaluate(changed)
        self.assertEqual(a['R'].tobytes(),rows.evaluate(s,principal=False)['R'].tobytes())
    def test_execution_diff_keeps_time_mesh_functions_and_height(self):
        default=Path(__file__).resolve().parents[4]/'docs/plans/wj-h13-linux/publication-v1/paired-case'
        parent=Path(os.environ.get('WJ_H14_PARENT',str(default)))/'flat/case.xml'
        with tempfile.TemporaryDirectory()as d:
            p=xml_case(parent,Path(d)/'case',Path(d)/'old-CP0',4)
            a=ET.parse(parent);b=ET.parse(p)
            for path in ('./GEOMETRY','./EXPANSIONS'):
                self.assertEqual(ET.tostring(a.find(path)),ET.tostring(b.find(path)))
            self.assertEqual([ET.tostring(v)for v in a.findall('./CONDITIONS/FUNCTION')],[ET.tostring(v)for v in b.findall('./CONDITIONS/FUNCTION')])
            params=lambda t:{p.text.split('=')[0].strip():p.text.split('=')[1].strip()for p in t.findall('./CONDITIONS/PARAMETERS/P')}
            for k in ('Time','TimeStep','RANSWallHeight','RANSWallRoughness'):self.assertEqual(params(a)[k],params(b)[k])
            self.assertEqual(params(b)['NumSteps'],'4')
    def test_small_gate_scope_and_no_target_reinitialization(self):
        from run_wj_h15 import halt_admission,CAPS
        b=dict(closed=False,coverage_halt={'first':'retained'})
        with self.assertRaises(AssertionError):halt_admission(b,'target')
        halt_admission(b,'closure');self.assertEqual(b['coverage_halt'],{'first':'retained'})
        self.assertEqual(CAPS['target'],1800);self.assertNotIn('target_cp0',CAPS)

if __name__=='__main__':unittest.main()
