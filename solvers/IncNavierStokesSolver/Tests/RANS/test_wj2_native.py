#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compiled WJ native serial-flat startup, independent loads and transactions.

Adapted from the transaction structure of test_cls_accepted_awf_native.py;
all WJ model, wall, audit and acceptance choices are explicit here.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import struct
import subprocess
import time
import uuid
import xml.etree.ElementTree as ET
import numpy as np
from numpy.testing import assert_allclose
from scipy.integrate import quad
from wj2_cases import case_xml
from cls_epawf_m3_cases import set_info,set_parameter,write
from run_cls_epawf_picard import THREADS,load,state_hash,sha,require,write_json
from compare_cls_epawf_pimple import same_history
from validate_wj1 import independent,CSTAR,QE,KAPPA


def stable(value):
    if isinstance(value,dict):return {k:stable(v) for k,v in sorted(value.items())}
    if isinstance(value,list):return [stable(v) for v in value]
    if isinstance(value,(float,int)):return struct.pack('<d',value).hex()
    return value


def digest_object(value):
    return hashlib.sha256(json.dumps(stable(value),sort_keys=True).encode()).hexdigest()


def wall(K,nu=1.5e-5,H=1.,h=.1):
    t=K/QE;A=(K/QE)**1.5/KAPPA
    cl=KAPPA*QE**1.5;a=KAPPA*math.sqrt(K/QE)
    x=h*math.sqrt(K)/(70*nu)
    b=10.7 if h==0 else 10.7*(1-x**max(.5-.4*x**.7,1-.79*x**(-.28)))
    zv=nu*b/math.sqrt(K)
    start=min(H,max(0.,zv));D=nu+a*max(start-zv,0)
    J0=start/nu+math.log1p(a*(H-start)/D)/a
    J2=start/nu**2+(H-start)/(D*(D+a*(H-start)))
    budget=t*t*(J0-nu*J2)
    c=H*math.exp(1-budget/A)
    return a,zv,t,A,c


def resting_functions(tree):
    """Physical direct-epsilon functions; no legacy eps_tilde facade."""
    set_info(tree,'RANSWallReference','Functions')
    K=QE*.25;epsilon=.125/KAPPA
    for family in ('InitialConditions','BoundaryValues'):
        for side in ('Inner','Outer'):
            f=ET.SubElement(tree.find('CONDITIONS'),'FUNCTION',{'NAME':family+side})
            values=dict(u=0.,v=0.,k=K,p=0.)
            if side=='Outer':values['epsilon']=epsilon
            for variable,value in values.items():ET.SubElement(f,'E',{'VAR':variable,'VALUE':format(value,'.17g')})


def independent_loads(checks,policy,offsets):
    n=int(offsets[-1])+1
    names=('inner_k_production','inner_k_destruction','outer_k_production','outer_k_destruction',
           'k_diffusion','epsilon_diffusion','epsilon_production','epsilon_destruction')
    expected={name:np.zeros(n) for name in names}
    source_count=0
    for row in checks['independent_quadrature']:
        k,e,G,w=row['k'],row['epsilon'],np.array(row['G']).reshape(3,3),row['weight']
        ids=np.array(row['ids'],dtype=int)+int(offsets[2]);phi=np.array(row['phi'])
        dphi=np.array(row['derivative_xz'])
        if row['inner']:
            a,zv,t,A,c=wall(row['sampled_K'])
            z=row['coordinate'][2];nut=a*max(z-zv,0);D=1.5e-5+nut
            strain=G+G.T
            P=.5*nut*np.sum(strain*strain)-2*k*np.trace(G)/3
            destruction=t*t*nut/D**2 if policy=='EP' else A/max(z,c)
            expected['inner_k_production'][ids]-=w*phi*P
            expected['inner_k_destruction'][ids]+=w*phi*destruction
            diffusivity=D
        else:
            c=independent(k,e,G);P=c['scalar'][3];diffusivity=CSTAR*k*k/e
            expected['outer_k_production'][ids]-=w*phi*P
            expected['outer_k_destruction'][ids]+=w*phi*e
            eids=np.array(row['epsilon_ids'],dtype=int)+int(offsets[3])
            ephi=np.array(row['epsilon_phi']);edphi=np.array(row['epsilon_derivative_xz'])
            expected['epsilon_production'][eids]-=w*ephi*1.44*e/k*P
            expected['epsilon_destruction'][eids]+=w*ephi*1.82*e*e/k
            ge=np.array(row['grad_epsilon'])[[0,2]]
            expected['epsilon_diffusion'][eids]+=w*diffusivity/1.3*np.sum(edphi*ge[:,None],axis=0)
        gk=np.array(row['grad_k'])[[0,2]]
        expected['k_diffusion'][ids]+=w*diffusivity*np.sum(dphi*gk[:,None],axis=0)
        source_count+=1
    errors={}
    for name in names:
        actual=np.array(checks['terms'][name])
        assert_allclose(actual,expected[name],rtol=1e-10,atol=2e-10,err_msg=name)
        errors[name]=float(np.max(abs(actual-expected[name])))
    for row in checks['extracted_constitutive']:
        k,e=row['k'],row['epsilon'];G=np.array(row['G']).reshape(3,3);H=np.array(row['direction_G']).reshape(3,3)
        require(np.max(abs(G[1]))==0 and np.max(abs(G[:,1]))==0 and np.max(abs(H[1]))==0 and
            np.max(abs(H[:,1]))==0,'native fixture lost the physical x-z embedding')
        c=independent(k,e,G)
        assert_allclose(np.array(row['R']).reshape(3,3),c['R'],rtol=3e-12,atol=3e-12)
        require(abs(row['divergence'])>.01 and abs(row['nu_eff']-row['nu_d'])>1e-5,'off-reference native action fixture is ineffective')
        assert_allclose([row['Praw'],row['Pdev']],c['scalar'][[3,4]],rtol=3e-12,atol=3e-12)
        h=1e-30
        d=independent(k+.13j*h,e-.04j*h,G.astype(complex)+1j*h*H,c['scalar'][12])
        assert_allclose(np.array(row['direction_R']).reshape(3,3),d['R'].imag/h,rtol=3e-10,atol=5e-11)
        assert_allclose(row['direction_sources_diffusion'],d['scalar'][[5,8,9,10]].imag/h,rtol=3e-10,atol=5e-11)
    assert_allclose(checks['endpoint_K_cutoff_plus_minus'],
        [wall(k)[4] for k in checks['endpoint_K_plus_minus']],rtol=3e-12,atol=3e-14)
    return dict(quadrature_points=source_count,load_max_absolute_errors=errors,
                extracted_states=len(checks['extracted_constitutive']))


def reference_comparison(checks,policy):
    """Continuum conditions and represented errors are separate receipts."""
    K=QE*.25;nu=1.5e-5;H=1.;a,zv,t,A,c=wall(K);d=nu-a*zv
    def J(z):return math.log1p(a*z/d)/a
    def Q(z):
        E=A*(z/c if z<c else 1+math.log(z/c))
        return E-t*t*(J(z)-nu*z/(d*(d+a*z)))
    require(zv<=0 and d>=2*nu and t*t*(J(H)-nu*H/(d*(d+a*H)))/A>=1,
        'continuum rough-reference proposition hypotheses failed')
    assert_allclose([Q(0),Q(H)],[0.,0.],atol=1e-12,rtol=0)
    cache={}
    def kref(z):
        if policy=='EP':return K
        if z not in cache:
            cuts=sorted(set([z,H]+[p for p in (c,.1) if z<p<H]))
            cache[z]=K-sum(quad(lambda s:Q(s)/(d+a*s),lo,hi,epsabs=1e-12,epsrel=3e-12)[0]
                for lo,hi in zip(cuts[:-1],cuts[1:]))
        return cache[z]
    errors=dict(velocity=0.,k=0.,outer_epsilon=0.,inner_pi=0.,outer_pi=0.,mechanical_pressure=0.)
    jump=-K/(81/20);pii=jump/2;pio=-jump/2
    for row in checks['independent_quadrature']:
        z=row['coordinate'][2];inner=row['inner'];kr=kref(z) if inner else K
        u=t*J(z) if inner else t*J(H)+.5/KAPPA*math.log(z/H)
        errors['velocity']=max(errors['velocity'],float(np.max(abs(np.array(row['velocity'])-[u,0,0]))))
        errors['k']=max(errors['k'],abs(row['k']-kr))
        key='inner_pi' if inner else 'outer_pi';pi=pii if inner else pio
        errors[key]=max(errors[key],abs(row['reduced_pressure']-pi))
        errors['mechanical_pressure']=max(errors['mechanical_pressure'],abs(row['mechanical_pressure_over_rho']-(pi-2*kr/3)))
        if not inner:errors['outer_epsilon']=max(errors['outer_epsilon'],abs(row['epsilon']-.125/(KAPPA*z)))
    interface=dict(shared_k_defect=0.,epsilon_trace_defect=0.,mechanical_jump_defect=0.,
                   maximum_traction_jump=0.,maximum_k_flux_jump=0.)
    for f in checks['interface_flux']:
        assert_allclose(f['opposing_normals'],[1,-1],atol=0,rtol=0)
        ki,ko=f['inner_outer_k'];ei,et=f['epsilon_trace_target'];pi,po=f['inner_outer_mechanical_pressure']
        interface['shared_k_defect']=max(interface['shared_k_defect'],abs(ki-ko))
        interface['epsilon_trace_defect']=max(interface['epsilon_trace_defect'],abs(ei-et))
        interface['mechanical_jump_defect']=max(interface['mechanical_jump_defect'],abs(pi-po-jump))
        interface['maximum_traction_jump']=max(interface['maximum_traction_jump'],
            float(np.max(abs(np.array(f['inner_outward_traction'])+f['outer_outward_traction']))))
        interface['maximum_k_flux_jump']=max(interface['maximum_k_flux_jump'],abs(sum(f['opposing_k_diffusive_flux'])))
    require(interface['shared_k_defect']<1e-12 and interface['epsilon_trace_defect']<1e-10,'shared value/one epsilon trace differs')
    return dict(continuum=dict(K=K,zv=zv,d=d,twice_nu=2*nu,Lambda=t*t*(J(H)-nu*H/(d*(d+a*H)))/A,
        rough_proposition_hypotheses=zv<=0 and d>=2*nu and t*t*(J(H)-nu*H/(d*(d+a*H)))/A>=1,
        Q0_endpoints=[Q(0),Q(H)],ground_reference_k=kref(0),pressure_jump=jump,
        inner_covariance_lower_bound=(2*QE/3-1)*t),represented_max_absolute_errors=errors,interface=interface)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--solver',type=Path,required=True);p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--policy',choices=('EP','Budget'),required=True)
    a=p.parse_args();require(platform.system()=='Linux','compiled Linux verification required')
    os.environ.update(THREADS)
    a.solver=a.solver.resolve();a.probe=a.probe.resolve();a.output=a.output.resolve()
    a.output.mkdir(parents=True,exist_ok=True)
    root=a.output/('run-'+uuid.uuid4().hex[:12]);root.mkdir()
    summary=dict(status='started',policy=a.policy,dt=1e-6,dt_ustar_over_H=5e-7,ranks=1,threads=THREADS,
        solver_sha256=sha(a.solver),probe_sha256=sha(a.probe),receipts=[],scope='native short-step transaction verification')
    def save():write_json(root/'summary.json',summary)
    def execute(command,folder):
        started=time.time()
        with (folder/'solver.log').open('w') as stream:
            result=subprocess.run(list(map(str,command)),cwd=folder,stdout=stream,stderr=subprocess.STDOUT,timeout=1800,env=os.environ.copy())
        receipt=dict(command=list(map(str,command)),exit_code=result.returncode,seconds=time.time()-started,
                     log=str(folder/'solver.log'),log_sha256=sha(folder/'solver.log'))
        write_json(folder/'receipt.json',receipt);summary['receipts'].append(receipt);save()
        return result.returncode
    def native(name,steps=0,cp=None,mutate=None,**kwargs):
        folder=root/name;xml=case_xml(folder,a.policy,steps,cp,**kwargs)
        if mutate:
            tree=ET.parse(xml).getroot();mutate(tree);write(tree,xml)
        rc=execute([a.solver,'case.xml'],folder)
        report=load(folder/'case.evolution.json') if (folder/'case.evolution.json').exists() else None
        return report,rc,folder
    def valid(report,rc):
        require(report is not None and rc==0 and not report.get('failure'),'native run failed: '+str(report.get('failure') if report else rc))
        require(report['turbulence_model']=='WJHighRe' and report['epsilon_semantics']=='epsilon','wrong native model/epsilon semantics')
        require(report['eta']==0 and report['order']==1 and report['audit_revision']=='NativeElementSideV2','wrong native configuration/provenance')
        require(report['initial_k_bound_check']['activations']==0,'startup projection changed represented reference')
        def domain_packets(value):
            if isinstance(value,dict):
                if 'dependent_trace_evaluation' in value:
                    evaluated=value['dependent_trace_evaluation']=='evaluated'
                    require((value['minimum_projected_trace_epsilon'] is not None)==evaluated,
                        'unevaluated dependent trace leaked a numeric sentinel')
                    if not evaluated:require('full declared' not in value['scope'],'partial domain labelled full')
                for child in value.values():domain_packets(child)
            elif isinstance(value,list):
                for child in value:domain_packets(child)
        domain_packets(report)
        for audit in report['state_audits']:
            require(audit['wj_hard_admissible']==[1.],'WJ hard audit failed')
            require(audit['wj_cls_max_beta04_R1_R2_cap_applicable']==[0.]*5,'CLS qualification leaked into WJ')
            witness=audit['wj_separate_analytic_adverse_witnesses']
            assert_allclose(witness,[4.130462147295278,-.00633378103539413,.27877891195383653,1.0357609543534572,
                3,-2/15,-2/15,34/15,48/5],rtol=3e-12,atol=3e-12)
        for attempt in report['pimple_attempts']:
            require(attempt['scheduled_complete'],'incomplete scheduled attempt')
            require(attempt['scalar_order']=='k_then_epsilon' and attempt['refreshed_trace_defect']<1e-11,
                'WJ sweep did not retain its final endpoint-K epsilon trace')
            for key in ('momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves','bound_calls'):
                require(attempt[key]==1,'scheduled solve count differs: '+key)
            require(attempt['additional_governing_equation_solves']==0 and attempt['bound_activations']==0,'hidden solve or material baseline projection')
            require(attempt['epsilon_frozen_linear_residual']<1e-8 and attempt['k_frozen_linear_residual']<1e-8,'frozen scalar solve failed')
    def frozen_readback(name,report,folder):
        out=root/name;out.mkdir();cp=Path(report['checkpoint'])
        payload=dict(operation='wj_checks',inner_stress_model='Linear',native_checkpoint=str(cp),checkpoint_sha256=sha(cp/'State.bin'),
            native_signature=report['signature'],state=report['state'],history=report['state'],time=report['time'],alpha=0.,
            native_differential_state=report['state'][:int(report['offsets'][4])],
            native_history=report['native_history'],native_history_times=report['native_history_times'])
        write_json(out/'input.json',payload)
        require(execute([a.probe,folder/'case.xml','input.json','result.json'],out)==0,'immutable native checks failed')
        result=load(out/'result.json')
        require(result['replay']['governing_solves']==0 and result['replay']['initialization_lift_calls_iterations']==[0.,0.],
                'read-only diagnostic solved or lifted state')
        require(result['source_balance']['normalization'] is None,'uncapped source readback normalized the wall')
        checks=result['checks'];require(checks['governing_solves']==0 and checks['state_unchanged'],'diagnostic mutated state')
        independent_result=independent_loads(checks,a.policy,report['offsets'])
        return dict(result=str(out/'result.json'),result_sha256=sha(out/'result.json'),independent=independent_result,
                    reference_comparison=reference_comparison(checks,a.policy),
                    max_same_state_defect=checks['max_same_state_defect'],pressure_modes=checks['pressure_modes'])
    try:
        zero,rc,initial=native('initial');valid(zero,rc)
        cp0=Path(zero['checkpoint']);immutable0={str(f.relative_to(cp0)):sha(f) for f in cp0.rglob('*') if f.is_file()}
        summary['zero_step']=dict(raw_reference_norms=zero['raw_reference_norms'],audit=zero['state_audits'][-1],
            reference_state_sha256=state_hash(zero['state']),checkpoint=str(cp0),checkpoint_hashes=immutable0)
        save()
        summary['native_checks']=frozen_readback('native-checks',zero,initial);save()
        # One bounded normal-mesh refinement and a separate source quadrature
        # check preserve physical H, zm, h and endpoint sampling.
        refinement,rc,_=native('normal-refinement',refined=True);valid(refinement,rc)
        quadrature,rc,_=native('quadrature-32',normal=32);valid(quadrature,rc)
        summary['representation_check']=dict(coarse=zero['raw_reference_norms'],refined=refinement['raw_reference_norms'],
            quadrature32=quadrature['raw_reference_norms'],coarse_unknowns=len(zero['state']),refined_unknowns=len(refinement['state']))
        functions,rc,_=native('direct-epsilon-functions',mutate=resting_functions);valid(functions,rc)
        assert_allclose(functions['state_audits'][-1]['wj_Praw_Pdev_ranges'],[0.]*4,atol=1e-13,rtol=0)
        assert_allclose(functions['state_audits'][-1]['wj_root_minimum_scaled_residual'][0],1.8,atol=1e-14,rtol=0)
        functions_back,rc,_=native('functions-zero-restart',0,Path(functions['checkpoint']),mutate=resting_functions)
        valid(functions_back,rc)
        require(same_history(functions,functions_back) and functions_back['initialization_lift_calls_iterations']==[0.,0.],
            'direct epsilon function checkpoint changed its state/history or lifted on readback')
        summary['direct_epsilon_function_semantics']=dict(exact_checkpoint_readback=True,
            raw_zero_velocity_norms=functions['raw_reference_norms'],epsilon_name=functions['epsilon_semantics'])
        whole,rc,_=native('whole2',2,cp0);valid(whole,rc)
        first,rc,_=native('split1',1,cp0);valid(first,rc);cp1=Path(first['checkpoint']);immutable1=sha(cp1/'State.bin')
        second,rc,second_folder=native('restart1',1,cp1);valid(second,rc)
        require(same_history(whole,second),'same-runtime restart primary/history mismatch')
        require(digest_object(whole['state_audits'][-1])==digest_object(second['state_audits'][-1]),'restart accepted audit packet differs')
        require(digest_object(whole['physical_samples'])==digest_object(second['physical_samples']),'restart physical samples differ')
        require(state_hash(whole['last_stage_history'])==state_hash(second['last_stage_history']) and
            whole['last_stage_time']==second['last_stage_time'],'restart last implicit history/time differs')
        summary['restart_exact']=True
        readback,rc,_=native('accepted-zero-step',0,cp1);valid(readback,rc)
        require(readback['initialization_lift_calls_iterations']==[0.,0.],'accepted WJ readback lifted a new reference/state')
        require(same_history(first,readback),'accepted zero-step primary/history mismatch')
        require(digest_object(first['physical_samples'])==digest_object(readback['physical_samples']),'zero-step readback physical samples differ')
        summary['accepted_readonly']=frozen_readback('accepted-readonly',first,root/'split1');save()
        rejected,rc,_=native('reject-after-hard-audit',1,cp1,hard_reject=2)
        require(rejected and 'deliberate' in rejected.get('failure',''),'scheduled hard rejection not exercised')
        require(same_history(readback,rejected),'rejection rollback primary/history mismatch')
        require(digest_object(readback['state_audits'])==digest_object(rejected['state_audits']),'rollback changed accepted audit packets')
        require(digest_object(readback['physical_samples'])==digest_object(rejected['physical_samples']),'rollback changed accepted physical samples')
        require(rejected['rejected_state'] and state_hash(rejected['rejected_state'])!=state_hash(rejected['state']),'failed candidate not retained separately')
        attempt=rejected['pimple_attempts'][-1]
        for key in ('momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves','bound_calls'):
            require(attempt[key]==1,'end-of-step rejection lost scheduled count '+key)
        summary['rollback_exact']=True
        def eps_semantics(t):
            for e in t.findall('./EXPANSIONS/E'):e.set('FIELDS',e.get('FIELDS').replace('epsilon','eps_tilde'))
            for v in t.findall('./CONDITIONS/VARIABLES/V'):
                if v.text=='epsilon':v.text='eps_tilde'
        def other_policy(t):
            other='Budget' if a.policy=='EP' else 'EP'
            set_info(t,'RANSWallFormulation','WJCalibrated'+('EPAWF' if other=='EP' else 'BudgetAWF'))
            set_info(t,'RANSWallTransition','UniformEPFlat' if other=='EP' else 'Fixed')
            set_parameter(t,'RANSWallChi',1 if other=='EP' else 0);set_parameter(t,'RANSWallTailStart',.5 if other=='EP' else 1)
        failures={}
        for name,mutate,expected,options in [
            ('positive-wall-nu',lambda t:set_parameter(t,'Kinvis',1.5001e-5),'signature mismatch',{}),
            ('valid-other-wall-policy',other_policy,'signature mismatch',{}),
            ('valid-refined-mesh',None,'signature mismatch',dict(refined=True)),
            ('changed-dt',None,'signature mismatch',dict(dt=2e-6)),
            ('unsupported-Cd',lambda t:set_parameter(t,'RANSWJCd',.12),'only the reviewed default coefficients',{}),
            ('gradient-policy',lambda t:set_info(t,'RANSWJGradientPolicy','Incompressible'),'requires explicit DevStrainRawWorkV1',{}),
            ('epsilon-semantics',eps_semantics,'WJ evolution variables must',{}),
            ('legacy-model',lambda t:set_info(t,'TurbulenceModel','CLSHighRe'),'legacy evolution variables must',{}),
            ('sampling',lambda t:set_info(t,'RANSWallKSampling','Point'),'WJ requires InterfaceTrace sampling',{}),
            ('physical-cap',lambda t:set_info(t,'RANSWallDestruction','Capped'),'WJ requires Original outer, Linear full-strain inner and uncapped destruction',{}),
            ('transition',lambda t:set_info(t,'RANSWallTransition','Upstream'),'WJ transitions are unsupported',{}),
            ('curved',lambda t:set_info(t,'RANSWallMap','Cylinder'),'WJ native supports serial flat 2D only',{}),
            ('evaluation-reuse',lambda t:set_info(t,'RANSWallEvaluationReuse','ExactStepV1'),'WJ requires explicit WJNativeV1 MassProjection and evaluation reuse Off',{}),
            ('extra-verification-solves',lambda t:set_info(t,'RANSWallPimpleVerify','On'),'WJ production route forbids RANSWallPimpleVerify=On',{})]:
            report,rc,folder=native('incompatible-'+name,0,cp1,mutate=mutate,**options)
            text=(folder/'solver.log').read_text()
            require(rc!=0 or (report and report.get('failure')),'incompatible input accepted '+name)
            require(expected in text,'wrong incompatible-input failure reason: '+name+' expected '+expected)
            failures[name]=dict(exit_code=rc,expected_reason=expected,signature_mismatch='signature mismatch' in text,log_sha256=sha(folder/'solver.log'))
        summary['incompatible']=failures
        if a.policy=='Budget':
            report,rc,folder=native('smooth-negative-reference',mutate=lambda t:set_parameter(t,'RANSWallRoughness',0))
            text=(folder/'solver.log').read_text();require(rc!=0 and 'reference_nonpositive' in text and '-0.330127' in text,
                'smooth negative reference was hidden or misclassified')
            summary['smooth_reference_nonpositive']=dict(log=str(folder/'solver.log'),sha256=sha(folder/'solver.log'))
        for name,nu in [('cutoff-failure',1.),('invalid-nu',0.)]:
            report,rc,folder=native(name,mutate=lambda t,nu=nu:set_parameter(t,'Kinvis',nu))
            expected='Calibrated AWF cutoff/regime inadmissible' if nu else 'EP-AWF requires positive input'
            require(rc!=0 and expected in (folder/'solver.log').read_text(),'wrong cutoff/scalar failure reason: '+name)
            summary[name]=dict(expected_reason=expected,log=str(folder/'solver.log'),sha256=sha(folder/'solver.log'))
        require(immutable0=={str(f.relative_to(cp0)):sha(f) for f in cp0.rglob('*') if f.is_file()} and sha(cp1/'State.bin')==immutable1,
                'saved checkpoint bytes changed')
        summary.update(status='PASS',reference_unchanged=True,source_revision=subprocess.check_output(
            ['git','-C',str(Path(__file__).resolve().parents[4]),'rev-parse','HEAD'],text=True).strip())
        summary['native_json_sha256']={str(f.relative_to(root)):sha(f) for f in root.rglob('*.json') if f.name!='summary.json'}
        save();print(root/'summary.json')
    except Exception as error:
        summary.update(status='FAIL',failure=str(error));save();raise


if __name__=='__main__':main()
