#!/usr/bin/env python3
"""M8-R2 exact affine and complete pressure-active native action verification.

No native time integration is performed. Independent Jacobi rows, dimensional
AWF primitives and explicit scalar differentiation give the answer; the native
transition/profile/action exports are comparisons, never oracle inputs.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import traceback
import xml.etree.ElementTree as ET

sys.dont_write_bytecode=True
THREADS={name:'1' for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
    'BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')}
os.environ.update(THREADS)
import numpy as np
from cls_epawf_m3_cases import set_info,set_parameter,write
from cls_epawf_m4_cases import tensor_mesh
from cls_epawf_m6_cases import transition_session
from cls_epawf_m8_cases import apply_variant
from cls_epawf_m6_operator_oracle import _diagnostic,_segment_projection
from cls_epawf_m7_energy_oracle import _fields,_horizontal,_profile
from cls_epawf_m8_native_oracle import normal_rule
from cls_epawf_m8_r1_nonlinear import Audit
from validate_cls_epawf_m2 import Profiles

R1_X=[-4.,-3.,-2.,-1.,1.,2.]
ACTIVE_X=[-4.,-3.9,-3.1,-3.,-2.,-1.,1.,2.]
ACTIVE_R=[0.,.1,.97,1.,1.03,1.97,2.,2.1,2.5,3.5,4.]
LEVELS=[2e-4,5e-5,1e-5,2e-6,5e-7]


def stamp():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def session(template,destination,active):
    transition_session(template,destination,family='flat',degree=3,dimension=2,steps=0,probe=True)
    root=ET.parse(destination).getroot()
    tensor_mesh(root,ACTIVE_R if active else [0.,2.,4.],3,2,6.,1.,points=8,
                xknots=ACTIVE_X if active else R1_X)
    apply_variant(root,'Regularized',.4,0.,.2,.1,False,1.)
    set_info(root,'RANSWallInnerStress','SaturatedR2' if active else 'Linear')
    set_info(root,'RANSWallNonlinearPolicy','Legacy')
    set_info(root,'RANSWallReference','Neutral')
    for name,value in {'RANSWallTransitionStart':-3.,'RANSWallTransitionEnd':-2.25,
        'RANSWallBaselineStart':-2.,'RANSWallSensorDistanceA':2.5,'RANSWallSensorDistanceB':3.5,
        'RANSWallSensorOn':.02,'RANSWallSensorOff':.1,
        'RANSWallTangentialOrder':8,'RANSWallNormalOrder':16,'RANSWallEta':0.,'RANSWallChi':1.,
        'RANSWallPressureProbeActive':int(active),'RANSWallPressureProbeDose':.06,
        'RANSWallPressureProbeDepth':.99,'RANSWallPressureProbeWidth':.03}.items():
        set_parameter(root,name,str(value))
    conditions=root.find('CONDITIONS')
    for node in list(conditions.findall('FUNCTION')):
        conditions.remove(node)
    return write(root,destination)


class PressureOracle:
    def __init__(self,raw,state=None):
        self.raw=raw;self.fields=_fields(raw);self.state=np.asarray(raw['state'] if state is None else state)
        self.d,self.off=raw['dimension'],raw['offsets'];self.H=raw['parameters'][1]
        self.audit=Audit(raw,Path(__file__).resolve().parents[4]);self.audit.gauge_variable_scale=raw['ustar']/self.H
        self.rays={};self.row_cache={}
        self.view=dict(raw,transition=dict(controls=raw['transition_controls'],weights=raw['transition_weights'],
                                          reference=raw['transition_reference']))

    def rows(self,block,points):
        points=np.atleast_2d(points);key=(block,points.tobytes())
        if key in self.row_cache:
            return self.row_cache[key]
        field=self.fields[block];groups={}
        for i,point in enumerate(points):
            groups.setdefault(field.element(point)['geometry'],[]).append(i)
        output=[]
        for geometry,indices in groups.items():
            ids,phi,gradient=field.rows(field.by_id[geometry],points[indices])
            output.append((np.asarray(indices,int),ids,phi,gradient))
        self.row_cache[key]=output
        return output

    def values(self,block,points,coefficients=None):
        points=np.atleast_2d(points);values=np.empty(len(points));gradient=np.empty((len(points),3))
        state=self.state if coefficients is None else np.asarray(coefficients)
        for indices,ids,phi,derivative in self.rows(block,points):
            c=state[self.off[block]+ids]
            values[indices]=phi@c;gradient[indices]=np.einsum('qaj,j->qa',derivative,c)
        return values,gradient

    def ray(self,chart):
        key=(float(chart[0]),float(chart[1]))
        if key in self.rays:
            return self.rays[key]
        diagnostic=_diagnostic(self.view,self.state,self.fields,chart,zero_body=True)
        profile,metric=_profile(self.raw,self.fields,self.state,chart)
        r,dr=normal_rule(profile,self.raw,self.fields,metric)
        points=np.column_stack([np.full_like(r,chart[0]),np.zeros_like(r),r])
        k=self.values(self.d,points)[0];pre=np.asarray(profile.fields(r),float)
        E,D=pre[2],pre[2]-pre[3];weights=dr*np.asarray(profile.jacobian(r),float)
        theta,delta,smooth=self.raw['robust_parameters'][2:5]
        t=np.minimum(1.,k/(theta*float(profile.k)))
        if np.any(t<0):
            raise ValueError('Pressure oracle received negative normalization k')
        limiter=10*t**3-15*t**4+6*t**5 if smooth else t
        B=float(weights@D);b=float(weights@(D*limiter))/B;m=max(b,delta);lam=b/m
        if profile.chi>0:
            corrected=float(profile.A)/np.maximum(r,float(profile.rc))
            original=float(profile.k)**1.5/(2.55*np.maximum(r,float(profile.ze)))
            dE=corrected-original
            # M2 exports the actual reference production as the second field.
            dD=pre[1]-original
            Bchi=float(weights@dD);Lchi=float(weights@(dD*limiter))
            dbchi=(Lchi-b*Bchi)/B
        else:
            dE=dD=np.zeros_like(r);dbchi=0.
        answer=dict(chart=np.asarray(chart),diagnostic=diagnostic,profile=profile,points=points,r=r,
            weights=weights,k=k,E=E,D=D,limiter=limiter,b=b,m=m,lam=lam,dE=dE,dD=dD,dbchi=dbchi)
        self.rays[key]=answer
        return answer

    def dchi(self,ray,direction):
        chart=ray['chart'].copy();chart[2]=self.H
        gradient=self.fields[self.d+3].evaluate(direction,self.off[self.d+3],chart,1)[1]
        ds=self.H*np.array([gradient[0],0.,0.])/(float(ray['profile'].k)/self.raw['parameters'][6]+
             (self.raw['parameters'][0]/self.H)**2)
        value=ray['diagnostic'];d2=value['gates'][-1];controls=self.raw['transition_controls']
        on,off=controls[5:7]
        if d2<=on*on or d2>=off*off or not self.raw['chi']:
            return 0.
        t=(d2-on*on)/(off*off-on*on)
        derivative=-30*t*t*(1-t)*(1-t)/(off*off-on*on)
        dd2=2*self.raw['transition_weights'][5]*float(value['sensors'][9:12]@ds)
        return float(value['gates'][1]*value['gates'][3]*derivative*dd2)

    def affine(self,direction):
        direction=np.asarray(direction);result=np.zeros(len(direction),np.longdouble)
        nu,H,h,rK,rm,K0,Q,kappa=self.raw['parameters']
        seed=Profiles(dict(geometry=[nu,H,h,rK,rm],metric=[0.,0.],k=K0,q=Q,kappa=kappa,chi=0.))
        metric=dict(height=H,scale=1.)
        r,dr=normal_rule(seed,self.raw,self.fields,metric)
        for chart,horizontal in _horizontal(self.fields,self.d,H,self.raw['quadrature_orders'][0]):
            points=np.column_stack([np.full_like(r,chart[0]),np.zeros_like(r),r])
            pressure=self.values(self.d+2,points,direction)[0]
            gradient=self.fields[self.d+3].evaluate(direction,self.off[self.d+3],chart,1)[1]
            force=np.array([gradient[0],0.,0.])
            for block,axis in enumerate((0,2)):
                for indices,ids,phi,derivative in self.rows(block,points):
                    weights=horizontal*dr[indices]
                    load=-(derivative[:,axis].T@(weights*pressure[indices]))
                    load-=phi.T@(weights*(r[indices]<=h)*force[axis])
                    result[self.off[block]+ids]+=load.astype(np.longdouble)
        for element in self.fields[self.d+3].elements:
            points,weights=self.fields[self.d+3].quadrature(element,self.raw['quadrature_orders'][0])
            pressure=self.values(self.d+3,points,direction)[0]
            for block,axis in enumerate((0,2)):
                for indices,ids,phi,derivative in self.rows(block,points):
                    del phi
                    result[self.off[block]+ids]-=(derivative[:,axis].T@(weights[indices]*pressure[indices])).astype(np.longdouble)
        for block in (self.d+2,self.d+3):
            field=self.fields[block];ids=np.arange(field.size)+self.off[block]
            mean=field.mass()@field.constant()/self.fields[0].volume
            result[ids]+=direction[-1]*mean
            result[-1]+=mean@direction[ids]
        for row,value in self.raw['constraints']:
            del value
            result[int(row)]=0.
        return np.asarray(result,float)

    def source_and_trace(self,direction):
        result=np.zeros(len(self.state));trace={};baseline_ids=set();records=[]
        for chart,horizontal in _horizontal(self.fields,self.d,self.H,self.raw['quadrature_orders'][0]):
            ray=self.ray(chart);dc=self.dchi(ray,direction)
            db=ray['dbchi']*dc;delta=self.raw['robust_parameters'][3]
            dm=db if ray['b']>delta else 0. if ray['b']<delta else max(db,0.)
            dD=dc*ray['dD']*ray['limiter']/ray['m']-ray['D']*ray['limiter']*dm/ray['m']**2
            for indices,ids,phi,gradient in self.rows(self.d,ray['points']):
                del gradient
                result[self.off[self.d]+ids]+=phi.T@(horizontal*ray['weights'][indices]*dD[indices])
            records.append(dict(chart=np.asarray(chart).tolist(),chi=float(ray['profile'].chi),dchi=dc,b=ray['b'],lambda_value=ray['lam']))
        field=self.fields[self.d+1]
        for element in field.elements:
            if abs(element['lower'][2]-self.H)>1e-11:
                continue
            nodes,transform=_segment_projection(element['degree'])
            xs=element['lower'][0]+(nodes+1)*(element['upper'][0]-element['lower'][0])/2
            values=[]
            for x in xs:
                ray=self.ray(np.array([x,0.,self.H]));dc=self.dchi(ray,direction)
                db=ray['dbchi']*dc;delta=self.raw['robust_parameters'][3]
                dm=db if ray['b']>delta else 0. if ray['b']<delta else max(db,0.)
                dl=db/ray['m']-ray['b']*dm/ray['m']**2
                profile=ray['profile'];preH=float(np.asarray(profile.fields(np.array([self.H])),float)[2,0])
                if dc:
                    original=float(profile.k)**1.5/(2.55*max(self.H,float(profile.ze)))
                    corrected=float(profile.A)/max(self.H,float(profile.rc))
                    value=dl*preH+ray['lam']*dc*(corrected-original)
                else:
                    value=0.
                values.append(value)
            coefficients=transform@values;n=element['degree']+1
            for index,value in zip(element['ids'][:n],np.asarray(element['sign'][:n])*coefficients):
                trace[int(index)]=float(value)
            if element['lower'][0]>=self.raw['transition_controls'][4]-1e-11:
                baseline_ids.update(map(int,element['ids'][:n]))
        for row,value in self.raw['constraints']:
            del value
            result[int(row)]=0.
        for index,value in trace.items():
            result[self.off[self.d+1]+index]=0. if index in baseline_ids else -value
        return result,records


def verify(raw):
    oracle=PressureOracle(raw);audit=oracle.audit;receipts=[]
    active=raw['fixture']=='pressure_sensor_and_capped_trace_active'
    samples=raw['samples']
    if active:
        audit.require(any(0<s['gates'][0]<1 and 0<s['gates'][2]<1 for s in samples),'Pressure transition interior was not exercised')
        audit.require(any(s['b']<.1 and 0<s['lambda']<1 for s in samples),'Active pressure fixture did not exercise capped attenuation')
    for sample in samples:
        ray=oracle.ray(np.asarray(sample['chart']))
        gates=ray['diagnostic']['gates']
        audit.near(sample['pressure_sensor'],ray['diagnostic']['sensors'][9:12],'independent_pressure_sensor')
        audit.near(sample['gates'],gates[[0,1,2,3,6]],'independent_pressure_transition_gates')
        audit.near([sample['K'],sample['b'],sample['m'],sample['lambda']],
                   [float(ray['profile'].k),ray['b'],ray['m'],ray['lam']],'independent_pressure_ray_normalization')
    constraints={int(row):value for row,value in raw['constraints']}
    for case in raw['directions']:
        direction=np.asarray(case['direction']);affine=oracle.affine(direction)
        nonlinear,records=oracle.source_and_trace(direction);complete=affine+nonlinear
        audit.near(case['affine_action'],affine,'independent_exact_affine_pressure_action')
        audit.near(case['complete_action'],complete,'independent_complete_pressure_action')
        error=audit.norm_record(np.asarray(case['complete_action'])-complete,constraints)
        audit.require(error['maximum']<=1e-7,'Complete pressure action fails fixed physical-norm1e-7 gate')
        if active and case['kind']=='outer_pressure':
            audit.require(np.max(np.abs(nonlinear))>1e-6,'Pressure/source/trace derivative is inactive')
            audit.require(np.max(np.abs(nonlinear[raw['offsets'][raw['dimension']+1]:raw['offsets'][raw['dimension']+2]]))>1e-6,
                          'Nonlinear projected epsilon trace derivative is inactive')
        levels=[]
        scale=max(1.,max(audit.norms(complete,constraints)))
        for level in case['levels']:
            action=(np.asarray(level['plus_residual'])-np.asarray(level['minus_residual']))/(2*level['step'])
            difference=audit.norm_record(action-complete,constraints)
            levels.append(dict(step=level['step'],full_residual_difference_error=difference,
                relative_full_residual_difference_error=difference['maximum']/scale,
                legacy_difference_fixed_threshold_met=difference['maximum']/scale<=1e-7,
                affine_only_omission_error=audit.norm_record(nonlinear,constraints)))
        receipts.append(dict(kind=case['kind'],independent_action_error=error,
            nonlinear_source_and_trace_norm=audit.norm_record(nonlinear,constraints),pressure_records=records,levels=levels))
    return dict(success=True,checks=audit.checks,errors=audit.errors,active_fixture=active,directions=receipts,
        scope='Exact affine and complete analytic pressure actions; full-residual finite-difference cancellation/branch errors retained without tolerance relaxation')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build',required=True,type=Path)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--r1-trace',type=Path)
    args=parser.parse_args()
    if platform.system()!='Linux':
        raise SystemExit('Scientific pressure verification must run on Linux')
    source=Path(__file__).resolve().parents[4]
    output=args.output or args.build/('rans-pressure-r2-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%f'))
    output.mkdir(parents=True,exist_ok=False)
    candidates=[args.build/'RANSPressureActionProbe',args.build/'library/UnitTests/SolverUtils/RANSPressureActionProbe']
    program=next((p for p in candidates if p.is_file()),None)
    report=dict(schema=1,milestone='M8-R2 C pressure actions',started_utc=stamp(),success=False,
        source_commit=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip(),
        source_clean_start=not bool(subprocess.check_output(['git','-C',str(source),'status','--porcelain'],text=True).strip()),
        driver_sha256=sha(__file__),threads=THREADS,commands=[],cases={},failures=[],native_time_trajectories=0)
    preregistration=dict(scope='operator-only; no time integration or physical qualification',
        selector_quiet='Linear',selector_active='SaturatedR2',nonlinear_policy='Legacy',
        component_tolerance=1e-10,physical_action_tolerance=1e-7,steps=LEVELS,
        pressure_dose=.06,depletion=.99,width=.03,inner_knots=ACTIVE_R,
        active_x=ACTIVE_X,model=dict(beta=.4,GammaScale=1.,theta=.2,delta=.1),
        active_requirements=['interior pressure-driven chi','b<delta with0<lambda<1',
                             'nonzero nonlinear k-source and projected epsilon-trace derivative'])
    (output/'preregistration.json').write_text(json.dumps(preregistration,indent=2,sort_keys=True)+'\n')
    report['preregistration_sha256']=sha(output/'preregistration.json')
    try:
        if not report['source_clean_start']:
            raise ValueError('Pressure validation requires the pinned clean source')
        if program is None:
            raise ValueError('Missing native RANSPressureActionProbe')
        report['program']=dict(path=str(program),sha256=sha(program))
        template=Path(__file__).with_name('EquilibriumWallEvolution.xml')
        for active,name in [(False,'affine'),(True,'pressure-active')]:
            directory=output/name;directory.mkdir();xml=session(template,directory/'case.xml',active)
            command=[str(program),str(xml),str(directory/'probe.json')]
            started=time.monotonic()
            with (directory/'run.log').open('w') as log:
                result=subprocess.run(command,cwd=directory,env=os.environ.copy(),stdout=log,stderr=subprocess.STDOUT,timeout=900)
            report['commands'].append(dict(command=command,cwd=str(directory),exit_code=result.returncode,
                seconds=time.monotonic()-started,ranks=1,threads=THREADS,xml_sha256=sha(xml),log_sha256=sha(directory/'run.log')))
            if result.returncode:
                raise ValueError('Native pressure probe failed: '+name)
            raw=json.loads((directory/'probe.json').read_text());report['cases'][name]=verify(raw)
            report['commands'][-1]['report_sha256']=sha(directory/'probe.json')
        if args.r1_trace:
            old=json.loads(args.r1_trace.read_text());old_raw=json.loads(args.r1_trace.with_name('case.evolution.json').read_text())
            old_raw=dict(old_raw,state=old['state'],constraints=old['constraints'])
            quiet=json.loads((output/'affine/probe.json').read_text())
            if quiet['offsets']!=old['offsets'] or quiet['fields']!=old_raw['fields']:
                raise ValueError('Retained R1 pressure localization has different native maps')
            oracle=PressureOracle(old_raw,np.asarray(old['state']))
            columns={int(item['column']):np.asarray(item['action']) for item in quiet['affine_columns']}
            case=next(item for item in old['directional_checks'] if item['kind']=='pressure_and_gauge')
            direction=np.asarray(case['direction']);exact=sum(direction[i]*value for i,value in columns.items())
            nonlinear,_=oracle.source_and_trace(direction);complete=exact+nonlinear
            constraints={int(row):value for row,value in old['constraints']};errors=[]
            independent=oracle.affine(direction)+nonlinear
            oracle.audit.near(complete,independent,'retained_R1_exact_native_affine_plus_independent_nonlinear')
            direct_error=oracle.audit.norm_record(complete-independent,constraints)
            oracle.audit.require(direct_error['maximum']<=1e-7,'Retained R1 exact action fails fixed physical-norm threshold')
            scale=max(1.,max(oracle.audit.norms(complete,constraints)))
            for level in case['levels']:
                action=(np.asarray(level['plus']['residual'])-np.asarray(level['minus']['residual']))/(2*level['step'])
                difference=oracle.audit.norm_record(action-complete,constraints)
                errors.append(dict(step=level['step'],full_residual_difference_error=difference,
                    relative_full_residual_difference_error=difference['maximum']/scale,
                    cached_error=oracle.audit.norm_record(np.asarray(case['cached_action'])-complete,constraints),
                    old_fresh_error=oracle.audit.norm_record(np.asarray(case['fresh_action'])-complete,constraints)))
            report['retained_R1_localization']=dict(trace_path=str(args.r1_trace),trace_sha256=sha(args.r1_trace),
                failing_state_is_not_physical_warmstart=True,nonlinear_pressure_action=oracle.audit.norm_record(nonlinear,constraints),
                exact_action_error=direct_error,levels=errors,
                interpretation='Exact native affine columns plus independently reconstructed nonlinear pressure source/trace; no inactive-chi assumption')
        report['success']=True
    except Exception as error:
        report['failures'].append(dict(type=type(error).__name__,message=str(error),traceback=traceback.format_exc()))
    report['finished_utc']=stamp()
    report['source_commit_end']=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()
    report['source_clean_end']=not bool(subprocess.check_output(['git','-C',str(source),'status','--porcelain'],text=True).strip())
    if not report['source_clean_end'] or report['source_commit_end']!=report['source_commit']:
        report['success']=False;report['failures'].append(dict(type='ProvenanceError',message='Source changed during pressure verification'))
    report['imported_source_files']={str(Path(module.__file__).resolve()):sha(module.__file__) for module in list(sys.modules.values())
        if getattr(module,'__file__',None) and str(Path(module.__file__).resolve()).startswith(str(source))}
    (output/'validation.json').write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print(json.dumps(dict(success=report['success'],output=str(output/'validation.json'))))
    return 0 if report['success'] else 1


if __name__=='__main__':
    raise SystemExit(main())
