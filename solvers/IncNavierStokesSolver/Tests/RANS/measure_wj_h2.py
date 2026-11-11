#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Versioned read-only H2 observations and incident reconstruction; no SVD."""
import argparse
import json
import math
import os
from pathlib import Path
import numpy as np
from numpy.testing import assert_allclose
from run_wj_h2 import execute,load,save,sha
from study_wj_h2 import observation_rule,comparison
from wj_h2_profiles import Completion,QE,KAPPA,D0,TOP,US,constitutive
from wj3_measure import checkpoint_payload,replay_witnesses
from validate_cls_epawf_m3 import basis
from validate_wj1 import Probe
from wj_h2_weak import independent_weak
from test_wj_h2_native import independent_boundary


def checkpoint_for(record,step):
    report=load(record['report']);choices=[]
    for p in Path(record['xml']).parent.glob('CP*'):
        if not (p/'State.bin').is_file():continue
        v=checkpoint_payload(p)
        if v['steps']==step:
            assert v['signature']==report['signature']
            first=report['steps']+1-len(report['accepted_states']);index=step-first
            assert index>=0
            assert np.asarray(v['state'],dtype='<f8').tobytes()==np.asarray(report['accepted_states'][index],dtype='<f8').tobytes()
            choices.append((str(p),v))
    assert choices,'no retained native checkpoint at required absolute step'
    return sorted(choices)[-1]


def field(report,state,block,q,geometry=None):
    H=report['parameters'][1]
    elements=report['fields'][block]['mesh']
    candidates=[e for e in elements if (geometry is None or e['geometry']==geometry) and
        all(e['lower'][a]-1e-12<=q[a]<=e['upper'][a]+1e-12 for a in (0,2))]
    assert candidates,'native geometry absent from independent reconstruction'
    e=min(candidates,key=lambda a:a['geometry']);p=e['degree'];phi=[];grad=[]
    for axis in (0,2):
        lo,hi=e['lower'][axis],e['upper'][axis];x=2*(q[axis]-lo)/(hi-lo)-1
        v,d=basis(p,[x]);phi.append(v[0]);grad.append(d[0]*2/(hi-lo))
    P=np.outer(phi[1],phi[0]).ravel();dx=np.outer(phi[1],grad[0]).ravel();dr=np.outer(grad[1],phi[0]).ravel()
    values=np.asarray(state)[np.array(e['ids'],int)+int(report['offsets'][block])]*e['sign']
    jac=1 if e['upper'][2]<=H else (TOP-D0-H)/(TOP-H)
    return float(P@values),np.array([dx@values,0,dr@values/jac]),e


def incident(report,kernel):
    event=report['physical_range_events'][-1];packet=event['wj_packet'];state=event['state']
    independent=replay_witnesses(packet);compiled={};reconstruction={}
    for name,w in packet['principal_witnesses'].items():
        q=w['chart'];geometry=w['geometry'];G=np.zeros((3,3));U=[]
        for b,axis in ((0,0),(1,2)):
            u,g,_=field(report,state,b,q,geometry);G[axis]=g;U.append(u)
        k,_,_=field(report,state,2,q,geometry);e,_,_=field(report,state,3,q,geometry)
        assert_allclose(G.ravel(),w['raw_G_row_major'],rtol=2e-11,atol=2e-11)
        assert_allclose([k,e],[w['k'],w['epsilon']],rtol=2e-12,atol=2e-13)
        xi=np.array(w['xi']);a=np.array(w['polarization']);H=np.outer(a,xi)
        c=kernel.closure(k,e,G);d=kernel.closure(k,e,G,(0,0,H))
        value=-np.sum(H*d['R'])/c['scalar'][1]
        assert_allclose(value,w['minimum_over_nu_eff'],rtol=2e-9,atol=2e-11)
        assert_allclose(np.linalg.eigvalsh(c['R']),w['covariance_eigenvalues'],rtol=2e-11,atol=2e-12)
        compiled[name]=dict(value=float(value),N=float(c['scalar'][12]),nu_eff=float(c['scalar'][1]),nu_d=float(c['scalar'][2]),
            covariance_eigenvalues=np.linalg.eigvalsh(c['R']).tolist(),Praw=float(c['scalar'][3]),Pdev=float(c['scalar'][4]))
        reconstruction[name]=dict(geometry=geometry,element=w['element'],chart=q,physical=w['position'],
            U_xz=U,G=G.tolist(),k=k,epsilon=e,scope='independent Jacobi tensor basis and affine composite inverse on declared native element; both normal-face owners remain distinct')
    return dict(time=event['time'],alpha=event['alpha'],reason=event['reasons'],state_sha256=packet['state_sha256'],
        callback_history_sha256=packet['stage_history_sha256'],independent_replay=independent,compiled_replay=compiled,
        independent_native_coefficient_reconstruction=reconstruction,witnesses=packet['principal_witnesses'],
        accepted_predecessor_time=report['time'],accepted_predecessor_step=report['steps'],
        correction_activations=sum(v['bound_activations'] for v in report['pimple_attempts']))


class Measurements:
    def __init__(self,root,observer,version='observations-v1'):
        self.root=Path(root).resolve();self.out=self.root/version;self.out.mkdir(exist_ok=False)
        self.observer=Path(observer).resolve();self.records={};self.env=dict(LD_LIBRARY_PATH=os.environ['LD_LIBRARY_PATH'])
        save(self.out/'intent.json',dict(observer=str(self.observer),observer_sha256=sha(self.observer),
            code_sha256=sha(__file__),role='immutable accepted-state observations, no continuation',normal_orders=[8,10],horizontal_orders=[4,6],
            quadrature_absolute_norm_floor=1e-10,original_study_sha256=sha(self.root/'study-v1/progress.json')))

    def observe(self,label,record,step,points,weights=None,weak=False,assembly=False):
        folder=self.out/label;folder.mkdir();cp,payload=checkpoint_for(record,step);report=load(record['report'])
        assert sha(record['report'])==record['report_sha256']
        packet=dict(operation='open_checks',inner_stress_model='Linear',native_report=record['report'],
            native_report_sha256=record['report_sha256'],native_checkpoint=cp,checkpoint_sha256=sha(Path(cp)/'State.bin'),
            allow_failed_producer=bool(report['failure']),physical_observations=points,weak=weak,assembly=assembly,reuse=False)
        save(folder/'request.json',packet)
        if weights is not None:save(folder/'weights.json',weights)
        assert execute(self.root,'measure-'+label,[self.observer,record['xml'],folder/'request.json',folder/'native.json'],folder,0,600,extra_env=self.env)
        result=load(folder/'native.json');assert result['selected_absolute_step']==step and result['initialization_work']==[0,0]
        source=dict(report_sha256=record['report_sha256'],checkpoint_sha256=packet['checkpoint_sha256'],
            callback_history_sha256=result['callback_history_sha256'],checkpoint=cp,step=step,time=payload['time'],native=str(folder/'native.json'),native_sha256=sha(folder/'native.json'))
        if weak:
            selected=dict(report,open_boundary=result['open_boundary'])
            source['independent_weak']=independent_weak(result['weak_checks'],selected)
            source['independent_boundary']=independent_boundary(selected)
        self.records[label]=source;save(self.out/'index.json',self.records)
        return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--observer',type=Path,required=True)
    p.add_argument('--kernel',type=Path,required=True);p.add_argument('--output-version',default='observations-v1')
    a=p.parse_args();m=Measurements(a.root,a.observer,a.output_version)
    study=load(a.root/'study-v1/progress.json');records=study['runs']
    c=Completion();designs=[load(Path(records[k]['xml']).with_name('case-design.json')) for k in ('base-pilot-on','refined-pilot')]
    result=dict(status='OBSERVATIONS_COMPLETE',original_study='OBSTRUCTED',unreached_target=8.,comparisons={},profiles={})
    kernel=Probe(str(a.kernel.resolve()))
    result['first_incident']=incident(load(records['refined-pilot']['report']),kernel)
    half=load(a.root/'diagnosis-v1/result.json');hrec=dict(xml=half['xml'],report=half['report'],report_sha256=half['report_sha256'])
    hd=load(half['report'])
    if hd['failure']:result['half_dt_incident']=incident(hd,kernel)
    result['half_dt']=half;kernel.close()
    K=c.K;E=(K/QE)**1.5/(KAPPA*c.H)
    low,weights=observation_rule(designs,8,4);high,wh=observation_rule(designs,10,6)
    # Matched accepted common time .04 across both meshes and timesteps.
    selected=[('base04',records['base-pilot-on'],2),('fine04',records['fine-pilot'],4),('refined04',records['refined-pilot'],2),
              ('base08',records['base-pilot-on'],4),('fine08',records['fine-pilot'],8)]
    payloads={}
    for label,record,step in selected:
        for order,points,w in ((8,low,weights),(10,high,wh)):
            payloads[(label,order)]=m.observe(f'{label}-q{order}',record,step,points,w,
                weak=label=='base08' and order==8,assembly=label=='base08' and order==8)
    for order,w in ((8,weights),(10,wh)):
        for label,x,y in [('time04','base04','fine04'),('normal04','base04','refined04'),('time08','base08','fine08')]:
            result['comparisons'][f'{label}-q{order}']=comparison(payloads[(x,order)],payloads[(y,order)],w,K,E)
    # Polynomial field norm quadrature, independently increased in x and z.
    quadrature={}
    for label in ('time04','normal04','time08'):
        l=result['comparisons'][f'{label}-q8']['fields'];h=result['comparisons'][f'{label}-q10']['fields']
        errors={f:abs(l[f]['normalized_L2']-h[f]['normalized_L2']) for f in ('U','k','epsilon','pi')}
        assert all(v<1e-10 for v in errors.values())
        quadrature[label]=dict(status='PASS',normalized_L2_changes=errors,floor=1e-10)
    result['quadrature']=quadrature
    heights=sorted(set([D0,D0+c.H,TOP]+[float(r['y_over_H']) for r in c.rows]+np.linspace(D0,TOP,161).tolist()))
    points=[[x,0,z] for x in (0.,2.,4.,6.,8.) for z in heights]
    for label,record,step in [('startup',records['base-cp0'],0),('base04',records['base-pilot-on'],2),('base08',records['base-pilot-on'],4),('refined04',records['refined-pilot'],2)]:
        r=m.observe('profiles-'+label,record,step,points)
        result['profiles'][label]=m.records['profiles-'+label]
        result['profiles'][label]['spatial_norms']=r['physical_mass_dual_norms']
    save(m.out/'summary.json',result);print(json.dumps(dict(status='COMPLETE',output=str(m.out),original_study='OBSTRUCTED')))
if __name__=='__main__':main()
