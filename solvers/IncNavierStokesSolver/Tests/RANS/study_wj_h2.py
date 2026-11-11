#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded physical precursor blocks, common-point measurements and admission.

No checkpoint/signature editing, hidden retries, projection or diagnostic solve.
The finite target and comparison criteria originate in wj-h2-contract.md.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import time
import numpy as np
from numpy.testing import assert_allclose
from scipy.special import roots_legendre
from wj_h2_profiles import Completion,NU,US,QE,KAPPA,D0,TOP
from wj_h2_cases import case_xml
from run_wj_h2 import execute,load,save,sha,THREADS
from compare_cls_epawf_pimple import same_history


def accepted(report):
    if report.get('failure'):return False
    if not report.get('wj_audit_packets'):return False
    for r in report.get('pimple_attempts',[]):
        for name in ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls'):
            if r.get(name)!=1:return False
        if r.get('bound_activations')!=0:return False
        if r.get('k_snapshot'):
            s=load(r['k_snapshot'])
            if np.asarray(s['raw'],dtype='<f8').tobytes()!=np.asarray(s['corrected'],dtype='<f8').tobytes():return False
    if report['state_audits'][-1]['wj_hard_admissible']!=[1]:return False
    return report['open_boundary']['gauge_unknowns']==0 and report['open_boundary']['boundary_flux_correction']==0


def observation_rule(designs,order=8,xorder=4):
    c=Completion();H=c.H
    knots={D0,D0+H,TOP}
    for d in designs:
        for rho in d['chart_z']:
            r=rho if rho<=H else H+(rho-H)*(TOP-D0-H)/(TOP-H)
            knots.add(D0+r)
        knots.add(D0+d['reference_cutoff']);knots.add(D0+c.h)
    z,wz=roots_legendre(order);x,wx=roots_legendre(xorder)
    observations=[];weights=[]
    for lo,hi in zip(sorted(knots)[:-1],sorted(knots)[1:]):
        for xi,wi in zip(x,wx):
            for zi,wj in zip(z,wz):
                observations.append([4*(1+float(xi)),0,(lo+hi)/2+(hi-lo)/2*float(zi)])
                weights.append(4*float(wi)*(hi-lo)/2*float(wj))
    assert_allclose(sum(weights),8*(TOP-D0),rtol=2e-15)
    return observations,weights


def comparison(a,b,weights,K,epsilonH):
    aa=a['open_boundary']['stations'];bb=b['open_boundary']['stations'];assert len(aa)==len(bb)==len(weights)
    assert abs(a['open_boundary']['time']-b['open_boundary']['time'])<=1e-12,'comparison physical times differ'
    result={}
    for field,amplitude in (('U',1.),('k',K),('epsilon',epsilonH),('pi',US**2)):
        values=[];ww=[]
        for x,y,w in zip(aa,bb,weights):
            assert_allclose(x['physical'],y['physical'],rtol=0,atol=3e-13)
            if field=='epsilon' and x['inner']:continue
            d=np.asarray(x[field])-np.asarray(y[field]);values.append(float(np.sum(d*d)));ww.append(w)
        l2=math.sqrt(sum(w*d for w,d in zip(ww,values))/sum(ww))/amplitude
        maximum=math.sqrt(max(values))/amplitude
        result[field]=dict(normalized_L2=l2,normalized_sampled_maximum=maximum,
            limit_L2=.02 if field!='pi' else None,limit_maximum=.05 if field!='pi' else None,
            status='PASS' if field=='pi' or (l2<=.02 and maximum<=.05) else 'UNRESOLVED')
    for key,amplitude,limit in (('K',K,.02),('cutoff',1.,.02)):
        a0=a['open_boundary']['interface'];b0=b['open_boundary']['interface'];assert len(a0)==len(b0)
        delta=max(abs(x[key]-y[key])/(amplitude if key=='K' else max(abs(x[key]),abs(y[key]),1e-15)) for x,y in zip(a0,b0))
        result[key]=dict(normalized_sampled_maximum=delta,limit=limit,status='PASS' if delta<=limit else 'UNRESOLVED')
    # Signed wall stress, measured on identical physical ground quadrature.
    a0=[v for v in a['open_boundary']['boundary_rows'] if v['role']=='ground']
    b0=[v for v in b['open_boundary']['boundary_rows'] if v['role']=='ground']
    assert len(a0)==len(b0)
    delta=max(abs(x['actual_traction'][0]-y['actual_traction'][0])/US**2 for x,y in zip(a0,b0))
    result['wall_stress']=dict(normalized_sampled_maximum=delta,limit=.05,status='PASS' if delta<=.05 else 'UNRESOLVED')
    return dict(status='PASS' if all(v['status']=='PASS' for v in result.values()) else 'UNRESOLVED',fields=result,
        pressure_convention='fixed physical traction level; no gauge realignment',points=len(weights),measure=sum(weights))


class Study:
    def __init__(self,root,runtime,observer):
        self.root=Path(root).resolve();self.runtime=Path(runtime).resolve();self.observer=Path(observer).resolve()
        self.folder=self.root/'study-v1';self.folder.mkdir(exist_ok=True)
        self.status=self.folder/'progress.json'
        self.summary=load(self.status) if self.status.exists() else dict(status='STARTED',runs={},observations={},comparisons={},target=8.)
        self.env=dict(LD_LIBRARY_PATH=os.environ['LD_LIBRARY_PATH'])
        identity=dict(schema='WJH2StudyV1',solver=dict(path=str(self.runtime),sha256=sha(self.runtime)),
            observer=dict(path=str(self.observer),sha256=sha(self.observer)),driver_sha256=sha(__file__),
            target=8.,dts=[.02,.01],refinement='normal only, identical measured completion/H/h_s',
            boundary='OpenTunnelV1/SlipV1/stop-backflow',maximum_block_steps=32,
            maximum_observation_interval=.5,unit_ceiling=2200,seconds_ceiling=21600)
        intent=self.folder/'intent.json'
        if intent.exists():assert load(intent)==identity,'changed study intent; new output root required'
        else:save(intent,identity)

    def persist(self):save(self.status,self.summary)

    def native(self,key,steps,dt=.02,refined=False,restart=None,reuse='ExactStepV1',length=.09):
        if key in self.summary['runs']:raise ValueError('immutable run key exists; inspect retained result before explicit continuation')
        folder=self.folder/key;xml=case_xml(folder,'physical',steps,dt,refined,reuse,length=length,restart=restart)
        start=time.monotonic()
        okay=execute(self.root,'study-'+key,[self.runtime,xml],folder,steps,900,extra_env=self.env)
        rp=xml.with_suffix('.evolution.json');report=load(rp) if rp.exists() else None
        passed=okay and report is not None and accepted(report)
        record=dict(status='PASS' if passed else 'OBSTRUCTED',xml=str(xml),xml_sha256=sha(xml),report=str(rp),
            report_sha256=sha(rp) if rp.exists() else None,requested_steps=steps,dt=dt,refined=refined,
            reuse=reuse,seconds=time.monotonic()-start,checkpoint=report.get('checkpoint') if report else None,
            steps=report.get('steps') if report else None,time=report.get('time') if report else None,
            failure=report.get('failure') if report else 'missing native report')
        self.summary['runs'][key]=record;self.persist()
        return report,record

    def observe(self,key,run,observations,weights=None,weak=False,assembly=False):
        if key in self.summary['observations']:raise ValueError('immutable observation exists')
        report=load(run['report']);assert sha(run['report'])==run['report_sha256'] and accepted(report)
        folder=self.folder/key;folder.mkdir();packet=dict(operation='open_checks',inner_stress_model='Linear',
            native_report=run['report'],native_report_sha256=run['report_sha256'],native_checkpoint=report['checkpoint'],
            checkpoint_sha256=sha(Path(report['checkpoint'])/'State.bin'),physical_observations=observations,
            assembly=assembly,reuse=False,weak=weak)
        save(folder/'request.json',packet)
        if weights is not None:save(folder/'weights.json',weights)
        okay=execute(self.root,'observe-'+key,[self.observer,run['xml'],folder/'request.json',folder/'result.json'],folder,0,600,extra_env=self.env)
        result=load(folder/'result.json') if okay else None
        r=dict(status='PASS' if okay else 'FAIL',path=str(folder/'result.json'),sha256=sha(folder/'result.json') if okay else None,
            producer_report_sha256=run['report_sha256'],points=len(observations),weak=weak,assembly=assembly)
        self.summary['observations'][key]=r;self.persist();return result

    def pilot(self):
        # Each dt has its own immutable native CP0; no changed-dt checkpoint reuse.
        cp0,z=self.native('base-cp0',0,reuse='Off')
        if z['status']!='PASS':self.summary.update(status='OBSTRUCTED',first_limiting_gate='represented physical CP0');self.persist();return
        base,b=self.native('base-pilot-off',4,restart=cp0['checkpoint'],reuse='Off')
        if b['status']!='PASS':self.summary.update(status='OBSTRUCTED',first_limiting_gate='coarse physical pilot');self.persist();return
        on,o=self.native('base-pilot-on',4,restart=cp0['checkpoint'])
        if o['status']!='PASS':self.summary.update(status='OBSTRUCTED',first_limiting_gate='exact-cache pilot');self.persist();return
        same_history(base,on)
        self.summary['physical_off_on_exact']=True
        fcp0,fz=self.native('fine-cp0',0,dt=.01)
        if fz['status']!='PASS':self.summary.update(status='OBSTRUCTED',first_limiting_gate='fine represented CP0');self.persist();return
        assert np.asarray(cp0['state'],dtype='<f8').tobytes()==np.asarray(fcp0['state'],dtype='<f8').tobytes()
        fine,f=self.native('fine-pilot',8,dt=.01,restart=fcp0['checkpoint'])
        if f['status']!='PASS':self.summary.update(status='OBSTRUCTED',first_limiting_gate='fine physical pilot');self.persist();return
        rcp0,rz=self.native('refined-cp0',0,refined=True)
        if rz['status']!='PASS':self.summary.update(status='OBSTRUCTED',first_limiting_gate='refined represented CP0');self.persist();return
        refined,r=self.native('refined-pilot',4,refined=True,restart=rcp0['checkpoint'])
        if r['status']!='PASS':self.summary.update(status='OBSTRUCTED',first_limiting_gate='refined physical pilot');self.persist();return
        designs=[load(Path(b['xml']).with_name('case-design.json')),load(Path(r['xml']).with_name('case-design.json'))]
        points,weights=observation_rule(designs)
        obs=[]
        for label,record in [('coarse',o),('fine',f),('refined',r)]:
            v=self.observe('pilot-'+label,record,points,weights,weak=label=='coarse',assembly=label=='coarse');assert v is not None;obs.append(v)
        K=designs[0]['K'];eps=(K/QE)**1.5/(KAPPA*.5)
        self.summary['comparisons']['pilot_time']=comparison(obs[0],obs[1],weights,K,eps)
        self.summary['comparisons']['pilot_normal']=comparison(obs[0],obs[2],weights,K,eps)
        save(self.folder/'pilot-comparisons.json',self.summary['comparisons'])
        rates={name:self.summary['runs'][key]['seconds']/n for name,key,n in [('base','base-pilot-on',4),('fine','fine-pilot',8),('refined','refined-pilot',4)]}
        book=load(self.root/'accounting.json');spent=sum(e.get('seconds',e['reserved_seconds']) for e in book['entries'].values())+book['prior_profile_setup_seconds']
        remaining=396*rates['base']+792*rates['fine']+396*rates['refined']+400*rates['base']+90*15
        projection=dict(measured_seconds_per_step_including_startup=rates,spent_seconds=spent,
            remaining_raw_seconds=remaining,headroom=1.25,closure_seconds=600,
            projected_total_seconds=spent+1.25*remaining+600,ceiling_seconds=21600,
            projected_units=2004+60,status='PASS' if spent+1.25*remaining+600<=21600 else 'RESOURCE_LIMIT')
        save(self.folder/'admission-v1.json',projection)
        statuses=[v['status'] for v in self.summary['comparisons'].values()]
        self.summary.update(status='PASS' if all(s=='PASS' for s in statuses) and projection['status']=='PASS' else 'UNRESOLVED' if any(s!='PASS' for s in statuses) else 'RESOURCE_LIMIT',
            first_limiting_gate=None if all(s=='PASS' for s in statuses) and projection['status']=='PASS' else 'pilot numerical comparisons' if any(s!='PASS' for s in statuses) else 'remaining cost projection')
        self.persist()

    def extend(self):
        assert self.summary.get('physical_off_on_exact') and self.summary['status']=='PASS','full continuation requires completed pilot admission'
        assert load(self.folder/'admission-v1.json')['status']=='PASS'
        last={name:self.summary['runs'][key] for name,key in [('base','base-pilot-on'),('fine','fine-pilot'),('refined','refined-pilot')]}
        designs=[load(Path(v['xml']).with_name('case-design.json')) for v in last.values()]
        points,weights=observation_rule(designs);K=designs[0]['K'];eps=(K/QE)**1.5/(KAPPA*.5)
        for target in [i*.5 for i in range(1,17)]:
            current={}
            for name in ('base','fine','refined'):
                dt=.01 if name=='fine' else .02;goal=round(target/dt);record=last[name]
                while record['steps']<goal:
                    count=min(25,goal-record['steps']);key=f'{name}-to-{record["steps"]+count:04d}'
                    report,record=self.native(key,count,dt,name=='refined',record['checkpoint'])
                    if record['status']!='PASS':
                        self.summary.update(status='OBSTRUCTED',first_limiting_gate=key);self.persist();return
                last[name]=record
                current[name]=self.observe(f'{name}-measure-{goal:04d}',record,points,weights,weak=target==8 and name=='base')
                if current[name] is None:self.summary.update(status='OBSTRUCTED',first_limiting_gate='passive measurement');self.persist();return
            for key,pair in [('time',('base','fine')),('normal',('base','refined'))]:
                c=comparison(current[pair[0]],current[pair[1]],weights,K,eps)
                self.summary['comparisons'][f'{target:g}-{key}']=c;self.persist()
                if c['status']!='PASS':self.summary.update(status='UNRESOLVED',first_limiting_gate=f'{target:g}-{key}');self.persist();return
        self.summary.update(status='PASS',completed_time=8.,development='not assessed from completion; station/residual criteria are separate');self.persist()


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--solver',type=Path,required=True)
    p.add_argument('--observer',type=Path,required=True);p.add_argument('--phase',choices=('pilot','extend'),required=True);a=p.parse_args()
    s=Study(a.root,a.solver,a.observer);getattr(s,a.phase)();print(json.dumps(s.summary,indent=2))
if __name__=='__main__':main()
