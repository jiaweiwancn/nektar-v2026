#!/usr/bin/env python3
"""Source-pinned Linux verification of the prescribed upstream transition.

Small manufactured runs verify the new coupling, not long-time CLS stability.
Raw meshes, states, logs and full-history checkpoints stay in the run directory.
"""
import argparse
from datetime import datetime,timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
import uuid
import xml.etree.ElementTree as ET
sys.dont_write_bytecode=True
THREADS={k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')}
os.environ.update(THREADS)
import numpy as np
from validate_cls_epawf_m4 import Runner as NativeRunner
from validate_cls_epawf_m5 import spaces
from cls_epawf_m6_cases import transition_session,arc_offset
from cls_epawf_m3_cases import write,set_info,set_parameter

checks=0
errors={}
def require(value,message):
    global checks
    checks+=1
    if not value:raise RuntimeError(message)
def near(a,b,name,scale=1.,tolerance=1e-10):
    a,b=np.asarray(a),np.asarray(b)
    require(a.shape==b.shape and np.all(np.isfinite(a)) and np.all(np.isfinite(b)),'invalid '+name)
    error=float(np.max(np.abs(a-b),initial=0.))/max(scale,float(np.max(np.abs(b),initial=0.)))
    errors[name]=max(errors.get(name,0.),error)
    require(error<=tolerance,f'{name}: {error} > {tolerance}')
    return error
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def state_difference(a,b,name,tolerance=1e-8):
    require(a['dimension']==b['dimension'] and a['offsets']==b['offsets'],'incompatible physical comparison spaces')
    off=a['offsets'];d=a['dimension'];values=[]
    for block,field in enumerate(spaces(a)):
        delta=np.asarray(a['state'][off[block]:off[block+1]])-b['state'][off[block]:off[block+1]]
        scale=.4 if block<d else a['parameters'][5] if block==d else .4**3/(a['parameters'][7]*2) if block==d+1 else .16
        value=math.sqrt(max(0.,float(delta@field.mass()@delta)))/math.sqrt(field.volume)/scale
        require(value<tolerance,f'{name} physical block {block}: {value}');values.append(value)
    near(a['state'][-1],b['state'][-1],name+'_gauge',tolerance=tolerance)
    errors[name]=max(values);return values

def expression_audit(run,family):
    name='expressions-'+family;directory=run.output/name;directory.mkdir()
    xml=transition_session(run.template,directory/'case.xml',family=family,steps=0)
    run.command(name,'RANSWallExpressionProbe',[xml,directory/'samples.json'])
    raw=json.loads((directory/'samples.json').read_text());qstar=3.1906750816242875;kappa=.4405539051866587
    K=qstar*.16*(1+.02*math.sin(.3));kd=qstar*.16*.02*math.cos(.3);nu=1.5e-5
    A=K**1.5/2.55;delta=(K/qstar)**1.5/kappa-A;offset=arc_offset() if family=='feet' else 0.
    xx=.1*math.sqrt(K)/(70*nu);m=max(.5-.4*xx**.7,1-.79*xx**(-.28));zv=10.7*nu*(1-xx**m)/math.sqrt(K)
    maximum=0.
    for sample in raw['samples']:
        x,r,_=sample['point'];s=min(1,max(0,(x-offset+3)/.75))
        chi=1-10*s**3+15*s**4-6*s**5;first=-30*s*s*(1-s)**2/.75;second=-60*s*(1-s)*(1-2*s)/.75**2
        B=A+delta*chi;nt=.2295*math.sqrt(K)*max(0,r-zv)
        p0=(K/qstar)**2*nt/(nu+nt)**2;ea=A/max(r,5.1*nu/math.sqrt(K))
        for entry in sample['values']:
            key=entry['name']
            if key.endswith('/k'):
                expected=kd+(1-chi)*ea+chi*p0 if key.startswith('ScalarSourceInner') else kd+B/r if key.startswith('ScalarSourceOuter') else K
            elif key.startswith('InitialConditionsOuter'):expected=B/r
            else:
                E=B/r;cmu=.3*(1-math.exp(-.36));trBr=0.
                if family=='feet' and abs(x)<1:
                    fx=-.02*math.pi/2*math.sin(math.pi*x);fxx=-.02*math.pi**2/2*math.cos(math.pi*x)
                    curvature=-fxx/(1+fx*fx)**1.5;trBr=curvature/(1+r*curvature)
                logBxx=delta*second/B-(delta*first/B)**2
                diffusion=cmu*K*K/1.3*(1/r**2-trBr/r+logBxx)
                ratio=K**1.5/(2.5*E*r);yap=.83*E*E/K*max(ratio-1,0)*ratio**2
                expected=1.5*kd/K*E-diffusion+1.92*E*E/K-yap
            error=near(entry['value'],expected,'native_prescribed_expression_'+key,tolerance=1e-10)
            maximum=max(maximum,error)
    return {'family':family,'samples':len(raw['samples']),'maximum_scaled_error':maximum,
            'meaning':'actual installed Nektar interpreter values versus independent strong source formulas'}

class Runner(NativeRunner):
    def probe(self,name,**kwargs):
        directory=self.output/name;directory.mkdir()
        xml=transition_session(self.template,directory/'case.xml',probe=True,**kwargs)
        self.command(name,'RANSCurvedWallProbe',[xml,directory/'probe.json'],kwargs.pop('ranks',1) if 'ranks' in kwargs else 1)
        return json.loads((directory/'probe.json').read_text())

def manufactured_error(raw,stationary=False):
    d=raw['dimension'];f=spaces(raw);off=raw['offsets'];time=raw['time']
    q,kappa=raw['parameters'][6:8];K0=q*.16;K=K0 if stationary else K0*(1+.02*math.sin(time))
    A=K**1.5/2.55;delta=(K/q)**1.5/kappa-A
    offset=arc_offset() if int(raw['map_parameters'][0])==4 else 0.
    square=np.zeros(d+4);volume=np.zeros(d+4)
    for b,field in enumerate(f):
        for e in field.elements:
            points,w=field.quadrature(e,24);ids,p,g=field.rows(e,points)
            actual=p@np.asarray(raw['state'])[off[b]+ids]
            if b<d or b>=d+2:exact=np.zeros_like(actual)
            elif b==d:exact=np.full_like(actual,K)
            else:
                s=np.clip((points[:,0]-offset+3)/.75,0,1)
                chi=1-10*s**3+15*s**4-6*s**5
                exact=(A+delta*chi)/points[:,2]
            square[b]+=np.dot(w,(actual-exact)**2);volume[b]+=np.sum(w)
    rms=np.sqrt(square/volume)
    values={'velocity':float(np.linalg.norm(rms[:d])),'k':float(rms[d]),'epsilon':float(rms[d+1]),
            'pressure':float(np.linalg.norm(rms[d+2:]))}
    require(raw['eta']==0,'transition/terrain anchor is nonzero')
    require(raw['upstream_transition'],'missing actual production transition')
    require(all(np.isfinite(list(values.values()))),'nonfinite manufactured error')
    trace_value=trace_gradient=measure=omission=0.
    z,w=np.polynomial.legendre.leggauss(32)
    for e in f[d+1].elements:
        if abs(e['lower'][2]-2)>1e-11:continue
        half=(e['upper'][0]-e['lower'][0])/2
        for xx,ww in zip(e['lower'][0]+half*(z+1),half*w):
            point=[xx,.37 if d==3 else 0.,2.]
            value,grad=f[d+1].evaluate(raw['state'],off[d+1],point,1)
            s=float(np.clip((xx-offset+3)/.75,0,1));chi=1-10*s**3+15*s**4-6*s**5
            derivative=-30*s*s*(1-s)**2/.75
            trace_value+=ww*(value-(A+delta*chi)/2)**2
            # The nonzero analytic tangential derivative occurs only on the
            # exactly flat fetch, where the physical and chart x axes agree.
            if xx-offset<=-2.:
                trace_gradient+=ww*(grad[0]-delta*derivative/2)**2
                omission+=ww*(delta*derivative/2)**2
            measure+=ww
    values.update(trace_value_rms=math.sqrt(trace_value/measure),trace_gradient_rms=math.sqrt(trace_gradient/measure),
                  omitted_grad_chi_rms=math.sqrt(omission/measure))
    require(values['omitted_grad_chi_rms']>1e-5,'manufactured trace has no detectable grad(chi) contribution')
    for sample in raw['transition_samples']:
        ell=sample['input'][0];gates=sample['gates']
        require(0<=gates[0]<=1,'analytic chi left [0,1]')
        if ell>=-2.25:near(gates[0],0,'no_downstream_reactivation')
        require(gates[2]==1 and gates[3]==1,'MMS left its declared state/admissibility unit plateau')
    return values

def stage_audit(raw):
    require(raw['eta']==0 and raw['upstream_transition'],'invalid terrain/transition stage policy')
    require(bool(raw['stage_audits']),'no native stage was solved')
    maximum=max(max(row[6:]) for row in raw['stage_audits'])
    require(maximum<1e-9,'unconverged transition coupled stage')
    near(np.sum(raw['regional_stage_loads'],axis=0),raw['final_stage_loads'],'regional_transition_stage_sum',tolerance=1e-10)
    require(abs(raw['state'][-1])<1e-9,'gauge multiplier would create a mass source')
    return {'maximum_scaled_residual':maximum,'unknowns':len(raw['state']),
            'boundary_surface_rms_correction':raw['boundary_flux_correction'],
            'energy_identity_error':abs(raw['energy_budget'][7])}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--build',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--mode',choices=['all','probes','smoke'],default='all')
    args=p.parse_args();require(platform.system()=='Linux','native numerical qualification belongs on Linux')
    if args.output.exists():args.output=args.output/('replay-'+uuid.uuid4().hex[:12])
    run=Runner(args.build,args.output);sections={};source=Path(__file__).resolve().parents[4]
    summary={'schema':1,'milestone':'M6','started_utc':datetime.now(timezone.utc).isoformat(),
             'source_root':str(source),'thread_environment':THREADS,'status':'running','sections':sections}
    def save():
        summary.update(checks=checks,errors=errors,commands=run.commands)
        (run.output/'validation.json').write_text(json.dumps(summary,indent=2)+'\n')
    try:
        if args.mode in ('all','probes'):
            from cls_epawf_m6_operator_oracle import verify
            sections['operators']=[]
            for family,dimension in [('flat',2),('feet',2),('flat',3)]:
                raw=run.probe(f'operator-{family}-d{dimension}',family=family,dimension=dimension,quadrature=8)
                sections['operators'].append(verify(raw,near,require,directional=dimension==2));save()
        if args.mode in ('all','smoke'):
            sections['prescribed_expressions']=[expression_audit(run,family) for family in ('flat','feet')];save()
            sections['manufactured']=[];cases={}
            for family,degrees in [('flat',(3,4,5)),('feet',(3,4))]:
                records=[]
                for degree in degrees:
                    name=f'mms-{family}-p{degree}'
                    raw=run.evolution(name,generator=transition_session,family=family,degree=degree,steps=2,dt=.0025,quadrature=12)
                    value={'family':family,'degree':degree,'refine':1,'error':manufactured_error(raw),'stage':stage_audit(raw)}
                    records.append(value);sections['manufactured'].append(value);cases[(family,degree)]=raw;save()
                require(records[-1]['error']['epsilon']<records[0]['error']['epsilon']*.8,'outer epsilon p refinement failed: '+family)
                require(records[-1]['error']['trace_gradient_rms']<records[0]['error']['trace_gradient_rms']*.8,'trace gradient p refinement failed: '+family)
            refined=run.evolution('mms-flat-h2',generator=transition_session,family='flat',degree=3,refine=2,steps=2,dt=.0025)
            er=manufactured_error(refined);stage_audit(refined)
            require(er['epsilon']<sections['manufactured'][0]['error']['epsilon']*.8,'fixed-physical-position h refinement failed')
            sections['h_refinement']={'refine':2,'error':er};save()
            # Identical prescribed MMS and mesh: compare actual full/Schur
            # stages, collective ownership, and every native history slot.
            base=cases[('feet',3)]
            condensed=run.evolution('condensed-feet',generator=transition_session,family='feet',degree=3,steps=2,dt=.0025,condensed=True)
            stage_audit(condensed);state_difference(base,condensed,'M6_full_condensed')
            sections['parallel_restart']=[]
            for ranks in (2,4):
                other=run.evolution('mpi-feet-'+str(ranks),generator=transition_session,family='feet',degree=3,steps=2,dt=.0025,ranks=ranks,ownership_shift=1)
                stage_audit(other);state_difference(base,other,'M6_MPI_'+str(ranks))
                near(other['native_history'],base['native_history'],'M6_MPI_history_'+str(ranks),tolerance=1e-8)
                sections['parallel_restart'].append({'ranks':ranks,'stage':stage_audit(other)});save()
            checkpoint=run.output/'restart-first'/'state'
            first=run.evolution('restart-first',generator=transition_session,family='feet',degree=3,steps=1,dt=.0025,checkpoint=str(checkpoint))
            resumed=run.evolution('restart-last',generator=transition_session,family='feet',degree=3,steps=1,dt=.0025,
                                  restart=first['checkpoint'],ranks=4,ownership_shift=2)
            stage_audit(resumed);state_difference(base,resumed,'M6_full_history_restart')
            near(resumed['native_history'],base['native_history'],'M6_complete_restart_history',tolerance=1e-8)
            near(resumed['native_history_times'],base['native_history_times'],'M6_complete_restart_history_times',tolerance=1e-12)
            sections['parallel_restart'].append({'restart_ranks':4,'all_history_slots':len(resumed['native_history'])});save()
            # Actual native startup IMEX1 on the same newly coupled operator.
            one=run.evolution('imex1-transition',generator=transition_session,family='flat',degree=3,steps=2,order=1,dt=.0025)
            sections['imex1']={'error':manufactured_error(one),'stage':stage_audit(one)};save()
            control=run.evolution('standard-kepsilon-transition',generator=transition_session,family='flat',degree=3,
                                  steps=2,dt=.0025,model='StandardKEpsilon')
            sections['standard_control']={'error':manufactured_error(control),'stage':stage_audit(control),
                                           'Q':control['parameters'][6],'kappa':control['parameters'][7]};save()
            # Fixed mesh/model/mask and common final physical time. Compare
            # actual complete coefficients to the finest native IMEX2 map;
            # this cancels spatial MMS error rather than fitting its floor.
            sections['temporal']=[];temporal=[]
            for dt in (.02,.01,.005,.0025):
                raw=run.evolution('temporal-'+str(dt),generator=transition_session,family='flat',degree=3,
                                  steps=round(.04/dt),dt=dt,order=2)
                stage_audit(raw);temporal.append(raw)
            fine=temporal[-1];fields=spaces(fine);off=fine['offsets'];temporal_errors=[]
            for raw in temporal[:-1]:
                square=0.
                for block in (raw['dimension'],raw['dimension']+1):
                    field=fields[block];delta=np.asarray(raw['state'][off[block]:off[block+1]])-fine['state'][off[block]:off[block+1]]
                    scale=fine['parameters'][5] if block==raw['dimension'] else .4**3/(fine['parameters'][7]*2)
                    square+=float(delta@field.mass()@delta)/field.volume/scale**2
                error=math.sqrt(max(0.,square));temporal_errors.append(error)
                sections['temporal'].append({'dt':raw['dt'],'scaled_scalar_mass_difference_to_dt_0025':error})
            require(temporal_errors[0]>temporal_errors[1]>temporal_errors[2]>0,'native transition temporal differences did not decrease')
            require(temporal_errors[0]/temporal_errors[1]>2.8 and temporal_errors[1]/temporal_errors[2]>2.8,'native transition IMEX2 did not show its expected temporal order')
            save()
            sections['rejections']=[]
            invalid=[('anchor',None,'RANSWallEta',1),('sensor-at-H',None,'RANSWallSensorDistanceA',2),
                     ('zero-buffer',None,'RANSWallBaselineStart',-2.25),('crossed-plateau',None,'RANSWallTransitionStart',-3.1),
                     ('crossed-baseline',None,'RANSWallBaselineStart',-1.9),('zero-weight',None,'RANSWallSensorWeightPressure',0),
                     ('unsafe-Re',None,'RANSWallSensorReynoldsOff',20),('outside-reference',None,'RANSWallSensorReynoldsOn',1e9),
                     ('backward',None,'RANSWallDirectionX',-1),('no-flat-fetch','RANSWallMap','Gaussian',None),
                     ('hill-support',None,'RANSWallTransitionEnd',-.5)]
            for label,info,key,value in invalid:
                name='reject-'+label;directory=run.output/name;directory.mkdir()
                xml=transition_session(run.template,directory/'case.xml',family='feet',steps=0)
                root=ET.parse(xml).getroot()
                if info:set_info(root,info,key)
                else:set_parameter(root,key,value)
                write(root,xml);run.command(name,'IncNavierStokesSolver',[xml],expect_failure=True,timeout=120)
                sections['rejections'].append(label);save()
        summary.update(status='passed',finished_utc=datetime.now(timezone.utc).isoformat());save()
        print(json.dumps({'status':summary['status'],'checks':checks,'commands':len(run.commands)}),flush=True)
    except Exception as error:
        summary.update(status='failed',error=str(error),finished_utc=datetime.now(timezone.utc).isoformat());save();raise

if __name__=='__main__':main()
