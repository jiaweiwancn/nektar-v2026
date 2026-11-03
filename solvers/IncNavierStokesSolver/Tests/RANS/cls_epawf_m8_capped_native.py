#!/usr/bin/env python3
"""Native severely capped source/epsilon-trace mechanism, not robust flow.

An exactly represented positive Q3 k field retains K(rK)=K0 while most of
each middle ray is depleted. The physical capped source must attenuate its
budget and dependent epsilon trace. Independent ray quadrature and native
trace-basis reconstruction check that coupling, including free-local-k
directions that leave the sampler and every exterior k value unchanged.
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
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
THREADS = {name:'1' for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
                                'BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')}
os.environ.update(THREADS)
import numpy as np
from cls_epawf_m3_cases import set_info, set_parameter, write
from cls_epawf_m4_cases import tensor_mesh, _evolution
from cls_epawf_m8_cases import apply_variant
from cls_epawf_m8_native_oracle import verify as verify_energy, ray_source
from cls_epawf_m7_energy_oracle import Audit, _fields
from cls_epawf_m6_operator_oracle import _segment_projection
from validate_cls_epawf_m7_transition_flow import failure_observation
import validate_cls_epawf_m7 as native_audit

DEPTH, WIDTH, DT, TOLERANCE = .99, .03, 1e-4, 1e-10
DIRECTIONS = (2e-4,1e-4)
X_KNOTS, R_KNOTS = (0.,.25,.75,1.), (0.,.1,.97,1.,1.03,2.,4.)


def stamp():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1048576),b''):
            value.update(block)
    return value.hexdigest()


def git(source,*arguments):
    return subprocess.check_output(['git','-c','gc.auto=0','-C',str(source),*arguments],text=True).strip()


def session(template,destination,steps=0,direction=0.):
    root = ET.parse(template).getroot()
    tensor_mesh(root,list(R_KNOTS),3,2,1.,1.,xknots=list(X_KNOTS))
    _evolution(root,steps,2,DT,0.,False,None,None)
    for name,value in {'RANSWallReference':'Functions','RANSWallTransition':'Fixed',
        'RANSWallMap':'Flat','RANSWallFunctionCoordinates':'Chart','TurbulenceModel':'CLSHighRe'}.items():
        set_info(root,name,value)
    for name,value in {'RANSWallHeight':2.,'RANSWallRoughness':.1,'RANSWallSampleHeight':1.,
        'RANSWallTailStart':1.,'Kinvis':1.5e-5,'RANSWallFrictionVelocity':.4,
        'RANSWallDirectionX':1.,'RANSWallDirectionY':0.,'RANSWallEta':0.,'RANSWallChi':1.,
        'RANSWallTangentialOrder':8,'RANSWallNormalOrder':16,'RANSWallMaxIterations':80,
        'RANSWallNonlinearTolerance':TOLERANCE,'RANSWallExpressionDimension':2,
        'CapDepth':DEPTH,'CapWidth':WIDTH,'CapDirection':direction}.items():
        set_parameter(root,name,repr(value))
    apply_variant(root,'CubicScaled',.4,0.,.2,.1,False,0.)
    # Dynamic zero/one operands avoid the independently documented installed
    # expression-interpreter constant-binary-operand bug; they change no values.
    left = '(min((1+0*t),max((0*t),(4*x))))'
    right = '(min((1+0*t),max((0*t),(4*(1-x)))))'
    ramp = f'((3*{left}^2-2*{left}^3)*(3*{right}^2-2*{right}^3))'
    radial = '(min((1+0*t),abs((y-1)/CapWidth)))'
    bump = f'(1-3*{radial}^2+2*{radial}^3)'
    base = f'(K0*(1-CapDepth*{ramp}*(1-{bump})))'
    local_direction = f'(K0*{ramp}*y*(y-1)*(y-2))'
    matching = f'(K0*(1-CapDepth*{ramp}))'
    outer_k = f'({matching}+(K0-{matching})*(y-2)/2)'
    initial_inner = dict(u='0',v='0',k=f'({base}+CapDirection*{local_direction})',p='0')
    initial_outer = dict(u='0',v='0',k=outer_k,eps_tilde='E0*RANSWallHeight/y',p='0')
    functions = {'InitialConditionsInner':initial_inner,'InitialConditionsOuter':initial_outer,
        'BoundaryValuesInner':dict(u='0',v='0',k=base),
        'BoundaryValuesOuter':dict(u='0',v='0',k=outer_k,eps_tilde='E0*RANSWallHeight/y')}
    conditions = root.find('CONDITIONS')
    for node in list(conditions.findall('FUNCTION')):
        conditions.remove(node)
    for name,values in functions.items():
        node = ET.SubElement(conditions,'FUNCTION',{'NAME':name})
        for variable,value in values.items():
            ET.SubElement(node,'E',{'VAR':variable,'VALUE':value})
    # No dummy perturbation, body force, scalar source, anchor or limiter of the
    # represented field: production Functions+Capped startup must impose trace.
    return write(root,destination)


def ramp(x):
    left,right = np.minimum(1.,np.maximum(0.,4*np.asarray(x))),np.minimum(1.,np.maximum(0.,4*(1-np.asarray(x))))
    return (3*left**2-2*left**3)*(3*right**2-2*right**3)


def exact_k(points,K0,direction=0.):
    points = np.asarray(points); x,r = points[:,0],points[:,2]
    t = np.minimum(1.,np.abs((r-1)/WIDTH)); f = 1-3*t*t+2*t*t*t
    inside = K0*(1-DEPTH*ramp(x)*(1-f))+direction*K0*ramp(x)*r*(r-1)*(r-2)
    lower = K0*(1-DEPTH*ramp(x)); outside = lower+(K0-lower)*(r-2)/2
    return np.where(r<=2,inside,outside)


def projected_trace(raw,state,fields,analytic_direction=False):
    field,H,d,offsets = fields[raw['dimension']+1],raw['parameters'][1],raw['dimension'],raw['offsets']
    expected,records = {},[]
    for element in field.elements:
        if abs(element['lower'][2]-H)>1e-12:
            continue
        nodes,transform = _segment_projection(element['degree'])
        xs = element['lower'][0]+(nodes+1)*(element['upper'][0]-element['lower'][0])/2
        values = []
        for x in xs:
            chart = np.array([x,0.,H]); ray = ray_source(raw,fields,state,chart)
            profile = ray['profile']; pre_h = float(np.asarray(profile.fields(np.array([H])),float)[2,0])
            value = ray['lam']*pre_h
            if analytic_direction:
                delta_k = raw['parameters'][5]*float(ramp(x))*ray['r']*(ray['r']-1)*(ray['r']-2)
                ratio = ray['local_k']/float(profile.k); theta = raw['robust_parameters'][2]
                slope = np.where(ratio<theta,delta_k/(theta*float(profile.k)),0.)
                tied = np.abs(ratio-theta)<1e-13
                slope[tied] = .5*delta_k[tied]/(theta*float(profile.k))
                db = float(ray['measure']@(ray['pre_d']*slope))/float(ray['measure']@ray['pre_d'])
                delta = raw['robust_parameters'][3]
                dlambda = db/delta if ray['b']<delta-1e-13 else 0. if ray['b']>delta+1e-13 else .5*db/delta
                value = dlambda*pre_h
            values.append(value)
            actual = fields[d+1].evaluate(state,offsets[d+1],chart,side=1)[0]
            records.append(dict(x=float(x),K=float(profile.k),b=ray['b'],m=ray['m'],lambda_value=ray['lam'],
                pre_trace=pre_h,analytic_target=ray['lam']*pre_h,actual_projected_trace=actual,
                direction_target=value if analytic_direction else None))
        coefficients = transform@np.asarray(values)
        count = element['degree']+1
        for index,value in zip(element['ids'][:count],np.asarray(element['sign'][:count])*coefficients):
            index,value = int(index),float(value)
            if index in expected and abs(expected[index]-value)>1e-10:
                raise RuntimeError('Independent trace transform has inconsistent shared endpoints')
            expected[index] = value
    return expected,records


def verify_state(raw,direction=0.,initial=False):
    audit = Audit()
    audit.require(not raw.get('failure') and raw.get('diagnostics_schema')==1,'Missing successful native cap report')
    audit.require(raw['model_variant']=='CubicScaled' and raw['eta']==0 and raw['chi']==1 and
        not raw['upstream_transition'] and not raw['mapped'],'Wrong isolated native capped-source mode')
    audit.near(raw['robust_parameters'],[.4,0.,.2,.1,0.,0.,1.],'fixed_declared_model')
    audit.require(not any('|'+name+'/' in raw.get('signature','') for name in
        ('BodyForceInner','BodyForceOuter','ScalarSourceInner','ScalarSourceOuter','InitialPerturbationInner','InitialPerturbationOuter')),
        'Fixture must not hide body/scalar sources or a startup workaround')
    fields,state,offsets,d = _fields(raw),np.asarray(raw['state']),raw['offsets'],raw['dimension']
    source = raw['state_audits'][-1]['robust_source']
    audit.require(source[5]>0 and 0<source[6]<1,'Severe capped regime did not attenuate the actual native ray budget')
    audit.near(source[3],0.,'native_zero_integral_redistribution')
    trace,records = projected_trace(raw,state,fields)
    ids = np.array(sorted(trace),int); target = np.array([trace[i] for i in ids])
    audit.near(state[offsets[d+1]+ids],target,'actual_dependent_trace_coefficients')
    audit.require(min(r['actual_projected_trace'] for r in records)>0,'Actual native matching epsilon trace is not positive')
    middle = [r for r in records if .25+1e-10<r['x']<.75-1e-10]
    audit.require(bool(middle) and all(0<r['b']<.1 and r['m']==.1 and 0<r['lambda_value']<1 for r in middle),
        'Middle matching faces did not use b<delta and lambda=b/delta')
    if initial:
        square = 0.
        for element in fields[d].elements:
            points,weights = fields[d].quadrature(element,6); ids_k,phi,_ = fields[d].rows(element,points)
            error = phi@state[offsets[d]+ids_k]-exact_k(points,raw['parameters'][5],direction)
            square += float(weights@(error*error))
        projection = math.sqrt(max(0.,square))/(raw['parameters'][5]*math.sqrt(fields[d].volume))
        audit.require(projection<1e-10,'Initial positive tensor-cubic k field was not represented exactly')
        audit.near([r['K'] for r in records],np.full(len(records),raw['parameters'][5]),'fixed_actual_K_sampler')
        for block in (0,1,d+2,d+3):
            audit.near(state[offsets[block]:offsets[block+1]],np.zeros(fields[block].size),'zero_initial_velocity_and_regional_pressure')
        audit.require(all(raw['state_audits'][-1][name][3]>0 for name in
            ('source_volume_inner','source_volume_outer','interface_inner','interface_outer')),
            'The zero-velocity positive initial fixture must have positive sampled covariance')
    else:
        projection = None
    native = native_audit.audits(raw,'severe-native-cap')
    energy = verify_energy(raw,full_volume=True)
    audit.require(energy['direct_volume_components_reconstructed'],'Missing independent full physical source/energy reconstruction')
    audit.require(raw['native_time_map'] and len(raw['native_history'])==2*raw['order'],'Missing native full-history state')
    audit.near(raw['native_history'][0],state[:offsets[d+2]],'complete_differential_native_history')
    return dict(success=True,checks=audit.checks+energy['checks'],errors=audit.errors,
        initial_k_projection_error=projection,trace_records=records,trace_ids=ids.tolist(),
        native_audit=native,independent_energy=energy,capped_volume=source[5],minimum_lambda=source[6],
        minimum_matching_epsilon=min(r['actual_projected_trace'] for r in records),
        pointwise_trace_projection_error=max(abs(r['actual_projected_trace']-r['analytic_target']) for r in records),
        scope='Actual projected dependent trace and discrete source budgets; pointwise analytic trace equality outside the representable plateau is not asserted')


def verify_directions(base,reports):
    audit = Audit(); fields = _fields(base); d,offsets = base['dimension'],base['offsets']; state = np.asarray(base['state'])
    expected,records = projected_trace(base,state,fields,True)
    ids = np.array(sorted(expected),int); target = np.array([expected[i] for i in ids])
    audit.require(np.max(np.abs(target))>1e-5,'The free-local-k trace direction is not a nonzero sensitivity control')
    outputs = []
    for step in DIRECTIONS:
        plus,minus = reports[(step,1)],reports[(step,-1)]
        audit.require(plus['fields']==base['fields']==minus['fields'],'Directional fixtures changed native spaces')
        yp,ym = np.asarray(plus['state']),np.asarray(minus['state'])
        actual = (yp[offsets[d+1]+ids]-ym[offsets[d+1]+ids])/(2*step)
        audit.near(actual,target,'native_free_k_to_trace_direction_'+str(step))
        for x in (0.,.125,.25,.5,.75,.875,1.):
            kp = fields[d].evaluate(yp,offsets[d],[x,0.,1.])[0]
            km = fields[d].evaluate(ym,offsets[d],[x,0.,1.])[0]
            audit.near([kp,km],[base['parameters'][5]]*2,'sampled_K_remains_fixed')
            for r in (0.,2.,4.):
                values=[fields[d].evaluate(y,offsets[d],[x,0.,r])[0] for y in (yp,ym,state)]
                audit.near(values[:2],[values[2]]*2,'fixed_exterior_and_matching_k')
        for x in (0.,1.):
            for r in R_KNOTS:
                values=[fields[d].evaluate(y,offsets[d],[x,0.,r])[0] for y in (yp,ym,state)]
                audit.near(values[:2],[values[2]]*2,'fixed_sidewall_k')
        other=np.setdiff1d(np.arange(fields[d+1].size),ids)
        audit.near(yp[offsets[d+1]+other],ym[offsets[d+1]+other],'epsilon_changes_only_in_matching_trace_coefficients')
        outputs.append(dict(step=step,actual_trace_direction=actual.tolist(),independent_direction=target.tolist(),
            maximum_error=float(np.max(np.abs(actual-target))),omitted_local_k_response_error=float(np.max(np.abs(actual)))))
    return dict(success=True,checks=audit.checks,errors=audit.errors,records=outputs,
        sampled_K_fixed=True,exterior_k_fixed=True,nonzero_response_to_free_local_k=True,
        analytic_middle_direction=[r for r in records if .25<r['x']<.75],
        scope='Native constrained trace finite differences with independently integrated full local-k normalization; not a separately exported Newton matrix')


def preregistration():
    return dict(schema=1,milestone='M8',purpose='exercise actual native b<delta,lambda<1 and its coupled epsilon trace',
        degree=3,xknots=list(X_KNOTS),rknots=list(R_KNOTS),H=2.,Z=4.,rK=1.,rm=1.,width=WIDTH,depletion=DEPTH,
        conditioning_choice='Before first execution, depth .99 leaves positive .01K; theta/delta/model are unchanged',
        model=dict(variant='CubicScaled',beta=.4,gamma=0.,theta=.2,delta=.1,smoothing='Linear',chi=1,eta=0),
        initial='U=0,regional pi=0,k=K0*(1-.99*A(x)*(1-F(r))),K(rK)=K0; outer k linearly recovers K0, epsilon seed E0*H/r',
        physical_initialization='Normal Functions+Capped startup Lift before SetReference; no dummy perturbation',
        source='No body/scalar/manufactured forcing; ordinary fixed essential boundary functions',
        directions=dict(steps=list(DIRECTIONS),formula='h*K0*A(x)*r*(r-1)*(r-2), inner only',
            invariants='zero at r=0,H,rK and x=0,1; sampler K and exterior k unchanged'),
        evolution=dict(steps=1,dt=DT,order=2,startup='native first-order startup',nonlinear_tolerance=TOLERANCE,maximum_iterations=80),
        acceptance=dict(isolated_identity=1e-10,coupled_residual=1e-9,minimum_lambda='strictly between0and1',
            source_budget='zero total redistribution with actual native normalization',trace='independent constrained-GLL trace coefficients'),
        failure_policy='Preserve native stop/rejected state; no field clipping, model retuning or automatic smaller-step retry',
        interpretation='A source/trace mechanism fixture; high-order discrete positivity and covariance remain independent, no robust-flow qualification')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--timeout',type=float,default=1800.)
    args = parser.parse_args()
    if platform.system()!='Linux':
        raise RuntimeError('Native capped-source validation belongs on Linux')
    source,build,output = Path(__file__).resolve().parents[4],args.build.resolve(),args.output.resolve()
    output.mkdir(parents=True,exist_ok=False);(output/'preregistration.json').write_text(json.dumps(preregistration(),indent=2)+'\n')
    report = dict(schema=1,milestone='M8',success=False,mechanism_verified=False,evolution_verified=False,
        candidate_qualified=False,started_utc=stamp(),source_commit=git(source,'rev-parse','HEAD'),
        source_clean_start=not git(source,'status','--porcelain'),script_sha256=sha(__file__),
        preregistration_sha256=sha(output/'preregistration.json'),commands=[],cases=[],failures=[],
        environment=dict(THREADS,TMPDIR=os.environ.get('TMPDIR'),cpu_affinity=sorted(os.sched_getaffinity(0))))
    solver = build/'IncNavierStokesSolver';report['program'] = dict(path=str(solver),sha256=sha(solver))
    environment = dict(os.environ,**THREADS,PYTHONDONTWRITEBYTECODE='1')
    environment['LD_LIBRARY_PATH'] = str(build/'geometry-lib')+':'+environment.get('LD_LIBRARY_PATH','')
    report['environment']['LD_LIBRARY_PATH'] = environment['LD_LIBRARY_PATH']
    def save():
        (output/'validation.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    def case(name,steps=0,direction=0.):
        directory=output/name;directory.mkdir();xml=session(Path(__file__).with_name('EquilibriumWallEvolution.xml'),directory/'case.xml',steps,direction)
        argv=[str(solver),str(xml)];log=directory/'run.log';started=time.monotonic()
        record=dict(name=name,steps=steps,direction=direction,xml=str(xml),xml_sha256=sha(xml),status='running')
        report['cases'].append(record);save()
        try:
            with log.open('w') as stream:code=subprocess.run(argv,cwd=directory,env=environment,stdout=stream,stderr=subprocess.STDOUT,timeout=args.timeout).returncode
        except subprocess.TimeoutExpired:
            code=124
        record.update(exit_code=code,seconds=time.monotonic()-started)
        report['commands'].append(dict(argv=argv,cwd=str(directory),ranks=1,threads=1,exit_code=code,seconds=record['seconds'],
            timeout=args.timeout,log=str(log),log_sha256=sha(log)))
        path=directory/'case.evolution.json';raw=json.loads(path.read_text()) if path.exists() else {}
        if path.exists():record.update(native_report=str(path),native_report_sha256=sha(path))
        if code or raw.get('failure') or not raw:
            record.update(status='stopped',observation=failure_observation(raw,code));save();return None,raw
        record['audit']=verify_state(raw,direction,initial=steps==0);record['status']='verified';save()
        print(name+': native cap mechanism verified',flush=True)
        return record,raw
    start=time.monotonic();save()
    try:
        if not report['source_clean_start']:
            raise RuntimeError('Capped native validation source must start pinned and clean')
        base_record,base=case('base-zero')
        if base_record is None:
            raise RuntimeError('Initial severe-cap fixture failed; no time-map or trace derivative is qualified')
        directions={}
        for step in DIRECTIONS:
            for sign in (1,-1):
                record,raw=case(('plus-' if sign>0 else 'minus-')+format(step,'.0e'),direction=sign*step)
                if record is None:
                    raise RuntimeError('A preregistered positive free-k directional fixture failed')
                directions[(step,sign)]=raw
        report['trace_direction']=verify_directions(base,directions);report['mechanism_verified']=True;save()
        record,raw=case('one-step',steps=1)
        if record is not None:
            audit=Audit();audit.require(len(raw['stage_audits'])==1 and len(raw['accepted_states'])==2,'One actual native stage was not completed')
            audit.require(max(raw['stage_audits'][0][4:5]+raw['stage_audits'][0][6:])<=1e-9,'Actual coupled stage did not converge')
            audit.near(raw['stage_audits'][0][0:3],[DT,DT,1/DT],'actual_native_startup_epoch')
            report['evolution_verified']=True;report['evolution_checks']=audit.checks
        else:
            report['failures'].append(dict(gate='native_time_map',interpretation='Stopped trajectory is not a positivity proof or qualified time map',observation=report['cases'][-1]['observation']))
        report['success']=report['mechanism_verified'] and report['evolution_verified']
    except Exception as error:
        report['failures'].append(dict(gate='mechanism_or_provenance',error=str(error)))
    try:
        report.update(source_commit_end=git(source,'rev-parse','HEAD'),source_clean_end=not git(source,'status','--porcelain'))
        report['program']['sha256_end']=sha(solver)
        if not report['source_clean_end'] or report['source_commit_end']!=report['source_commit'] or report['program']['sha256_end']!=report['program']['sha256']:
            raise RuntimeError('Capped native source or solver image changed')
    except Exception as error:
        report['success']=False;report['failures'].append(dict(gate='provenance',error=str(error)))
    report.update(status='verified' if report['success'] else 'mechanism_verified_time_map_failed' if report['mechanism_verified'] else 'failed',
        finished_utc=stamp(),seconds=time.monotonic()-start);save()
    print(json.dumps({key:report[key] for key in ('success','mechanism_verified','evolution_verified','candidate_qualified','status','failures')}),flush=True)
    return 0 if report['success'] else 1


if __name__=='__main__':
    raise SystemExit(main())
