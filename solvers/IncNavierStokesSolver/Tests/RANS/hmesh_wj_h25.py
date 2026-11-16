#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H25 parent-knot representation screen; no claim of evolved spatial accuracy."""
import argparse,collections,hashlib,os,subprocess
from pathlib import Path
import xml.etree.ElementTree as ET
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
import numpy as np
from wj_h15_observer import identity,verify,atomic,pkey,BoundedIndex,add_witnesses,assess
from wj_h17_observer import native_load,load_compatible
from wj_h4_points import physical
from wj_h7_analysis import Rows
from wj_h19_analysis import compare
from cls_epawf_m4_cases import tensor_mesh
from cls_epawf_m3_cases import set_info,set_parameter,write

def directory(root,label):
    assert label.replace('-','').isalnum()
    return root/'resumption-r1/science'/label

def initialize(root,label,successor=None):
    inputs=load(root/'inputs/baseline.json');original=verify(inputs['own_cp0']['xml'])
    x=ET.parse(original).getroot();vertices=[list(map(float,v.text.split()))for v in x.findall('./GEOMETRY/VERTEX/V')]
    xs=sorted({v[0]for v in vertices});fine=sorted({v[1]for v in vertices});parent=fine[::2]
    assert len(xs)==15 and len(fine)==51 and len(parent)==26 and parent[8]==.15
    restored=[]
    if successor:
        proposal=load(successor);assert proposal['status']=='SUCCESSOR_PROPOSAL_NOT_ACCEPTANCE'
        verify(proposal['source']);restored=proposal['restore_pairs'];assert set(restored)<=set(range(25))
        parent=sorted(set(parent+[fine[2*i+1] for i in restored]))
        assert (len(parent)-1)*14==proposal['cells']
    inner=parent.index(.15);outer=len(parent)-1-inner
    states=3*43*(3*(inner+outer)+1)+43*(3*outer+1)+15*(inner+outer+2)
    d=directory(root,label);d.mkdir();tensor_mesh(x,parent,3,2,8.,1.,points=[8,8],xknots=xs)
    for e in x.findall('./EXPANSIONS/E'):e.set('FIELDS',e.get('FIELDS').replace('eps_tilde','epsilon'))
    for e in x.findall('./CONDITIONS/VARIABLES/V'):
        if e.text=='eps_tilde':e.text='epsilon'
    design=dict(mesh=f'14x({inner}+{outer}) stored parent intervals and diagnosed restored children',reference_xml=identity(original),reference_knots=fine,parent_knots=parent,xknots=xs,
        states=states,cells=14*(inner+outer),restored_pairs=restored,successor=identity(successor) if successor else None,
        governing_quadrature=[12,16],fresh_representation=True,not_a_H24_restart=True)
    atomic(d/'design.json',design)
    for k,v in dict(RANSWallNativeLookup='TensorExactV1',RANSWallScalarEvaluation='Off',RANSWallMappedResourceDesign=identity(d/'design.json')['sha256'],
        RANSWallRestart='None',RANSWallCheckpoint=str(d/'CP'),RANSWallReport=str(d/'native.json'),RANSWallRetention='Off',RANSWallMonitoringPolicy='Off',
        RANSWallPimpleSnapshots='Off',RANSWallPimpleFlowVectors='Off',RANSWallPimpleArithmetic='Off',RANSWallPimpleVerify='Off').items():set_info(x,k,v)
    for k,v in dict(NumSteps=0,TimeStep=.0025,Time=0,IO_CheckSteps=1,RANSWallDiagnostics=1).items():set_parameter(x,k,v)
    write(x,d/'case.xml');atomic(d/'execution-intent.json',dict(xml=identity(d/'case.xml'),design=identity(d/'design.json'),native=identity(root/'science/runtime-build-scalar-v1.json'),
        governing_attempts=0,initialization_role='new h representation from unchanged physical Functions; native projection/lift/FlatInletTrace and history counted, never a relabelled H24 field'))
    with(d/'native.log').open('x')as f:code=subprocess.run([str(root/'runtime-scalar-v1/IncNavierStokesSolver'),str(d/'case.xml')],stdout=f,stderr=subprocess.STDOUT).returncode
    atomic(d/'launch-result.json',dict(exit=code));assert code==0,'h350 native initial admission failed; preserve original report'
    r=native_load(d/'case.evolution.json');assert int(r['steps'])==0 and r['time']==0 and r['failure']==''and len(r['state'])==states
    atomic(d/'initial-source.json',dict(report=identity(d/'case.evolution.json'),checkpoint=identity(d/'CP/State.bin'),history_role='genuine native initial history, alpha0'))

def analyze(root,label,candidate_report=None):
    from wj_h16_observer import compatibility
    d=directory(root,label);candidate_report=candidate_report or d/'case.evolution.json';r=native_load(candidate_report)
    prior=load(root/'science/cache-compatible-v2.json')
    bound=compatibility(root,verify(prior['historical_manifest']),d/'fine-cache-compatibility-base.json')
    intent=load(d/'execution-intent.json');candidate_runtime=intent.get('runtime',intent.get('native'))
    assert candidate_runtime is not None;verify(candidate_runtime)
    bound.update(new_native_identity=candidate_runtime,
        permitted_difference='fine immutable reference rows only; coarse rows rebuilt independently on its new field/native map; actual numerical helper paths verified')
    atomic(d/'fine-cache-compatibility.json',bound)
    cache,_=load_compatible(d/'fine-cache-compatibility.json');c=cache.config
    fixed=load(verify(load(root/'inputs/baseline.json')['support']))['points'];points=collections.OrderedDict((pkey(p),p)for p in fixed)
    lobatto=[-1.,-np.sqrt(3/7),0.,np.sqrt(3/7),1.]
    for e in r['fields'][0]['mesh']:
        for a in lobatto:
            for b in lobatto:
                p=physical(r,[(e['lower'][0]+e['upper'][0])/2+a*(e['upper'][0]-e['lower'][0])/2,0.,
                    (e['lower'][2]+e['upper'][2])/2+b*(e['upper'][2]-e['lower'][2])/2])
                points.setdefault(pkey(p),list(p))
    for packet in r['wj_audit_packets']:add_witnesses(points,packet,cache.report,c)
    xs=sorted({p[0]for p in points.values()if abs(p[2]-(c.d0+c.H))<1e-14})
    index=BoundedIndex(r)
    for report,state,idx in ((r,r['state'],index),(cache.report,cache.qh0,cache.index)):
        probes=Rows(report,[[x,0.,c.d0+c.H]for x in xs],c,owner_lookup=idx,geometry_cache=idx,include_volume=False)
        values=probes.evaluate(state,principal=False)
        for i,(p,owner)in enumerate(probes.rows):
            for delta in (-1e-6,0.,1e-6):
                q=[p[0],0.,c.d0+float(values['rc'][i])+delta];points.setdefault(pkey(q),q)
    fine_rows,initial,stats=cache.view(list(points.values()));fine=fine_rows.evaluate(cache.qh0)
    selected=[]
    for p,o in fine_rows.rows:
        owners=[v for v in index(p)if v['region']==o['region']and all(min(v['upper'][a],o['upper'][a])>max(v['lower'][a],o['lower'][a])for a in (0,2))]
        assert len(owners)==1,'ambiguous/missing common one-sided spatial owner'
        selected.append(owners[0])
    it=iter(selected)
    rows=Rows(r,[p for p,o in fine_rows.rows],c,owner_lookup=lambda p:[next(it)],geometry_cache=index,include_volume=False)
    rows.weights=fine_rows.weights.copy();rows.point_count=fine_rows.point_count
    assert np.array_equal(rows.outer,fine_rows.outer)and len(rows.rows)==len(fine_rows.rows)
    coarse_physical,coarse=assess(rows,r['state'],initial,0);fine_physical,_=assess(fine_rows,cache.qh0,initial,0)
    result,actions=compare(rows,coarse,fine,initial)
    wall=[i for i,(p,o)in enumerate(rows.rows[:rows.point_count])if o['chart'][2]==0]
    wall_error=float(np.max(abs(coarse['T'][wall,0,2]-fine['T'][wall,0,2]))/c.u_tau**2)
    status='PASS_INITIAL_REPRESENTATION_ONLY'if result['status']=='PASS'and wall_error<=1e-3 and coarse_physical['summary']['status']=='PASS'and fine_physical['summary']['status']=='PASS'else'REJECTED'
    payload={name+'_'+k:v for name,data in (('coarse',coarse),('fine',fine),('fine_qh0',initial))for k,v in data.items()}
    payload.update(physical=np.asarray([p for p,o in rows.rows]),weights=rows.weights,
        coarse_geometry=np.asarray([o['geometry']for p,o in rows.rows]),fine_geometry=np.asarray([o['geometry']for p,o in fine_rows.rows]))
    np.savez_compressed(d/'common-signed-fields.npz',**payload)
    atomic(d/'comparison.json',dict(status=status,numerical_screen=result,coarse_physical=coarse_physical,fine_reassessed_physical=fine_physical,
        wall_stress_scaled_max=wall_error,wall_limit=1e-3,source=identity(d/'initial-source.json') if (d/'initial-source.json').exists() else identity(candidate_report),signed_fields=identity(d/'common-signed-fields.npz'),
        common_inventory=dict(fixed_points=len(points),owned_rows=rows.point_count,total_rows=len(rows.rows),fine_cache=stats,
            rule='fine fixed support + all coarse5x5 and native witnesses/neighborhoods + both current cutoffs; all one-sided traces paired by positive interior overlap'),
        distinction='initial target mesh/quadrature representation; physical discrepancy measured, no evolved discretization qualification or H24 state migration claimed'))
    print(status,result['status'],'wall',wall_error,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('initialize','analyze'),required=True)
    p.add_argument('--label',default='h350-r1-v1')
    p.add_argument('--candidate-report',type=Path)
    p.add_argument('--successor',type=Path)
    a=p.parse_args()
    if a.mode=='initialize':initialize(a.root,a.label,a.successor)
    else:analyze(a.root,a.label,a.candidate_report)
