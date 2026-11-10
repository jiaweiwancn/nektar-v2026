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
        c=independent(k,e,G)
        assert_allclose(np.array(row['R']).reshape(3,3),c['R'],rtol=3e-12,atol=3e-12)
        require(abs(row['divergence'])>.01 and abs(row['nu_eff']-row['nu_d'])>1e-5,'off-reference native action fixture is ineffective')
        assert_allclose([row['Praw'],row['Pdev']],c['scalar'][[3,4]],rtol=3e-12,atol=3e-12)
        h=1e-30
        d=independent(k+.13j*h,e-.04j*h,G.astype(complex)+1j*h*H,c['scalar'][12])
        assert_allclose(np.array(row['direction_R']).reshape(3,3),d['R'].imag/h,rtol=3e-10,atol=5e-11)
        assert_allclose(row['direction_sources_diffusion'],d['scalar'][[5,8,9,10]].imag/h,rtol=3e-10,atol=5e-11)
    return dict(quadrature_points=source_count,load_max_absolute_errors=errors,
                extracted_states=len(checks['extracted_constitutive']))


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
        for audit in report['state_audits']:
            require(audit['wj_hard_admissible']==[1.],'WJ hard audit failed')
            require(audit['wj_cls_max_beta04_R1_R2_cap_applicable']==[0.]*5,'CLS qualification leaked into WJ')
        for attempt in report['pimple_attempts']:
            require(attempt['scheduled_complete'],'incomplete scheduled attempt')
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
        checks=result['checks'];require(checks['governing_solves']==0 and checks['state_unchanged'],'diagnostic mutated state')
        independent_result=independent_loads(checks,a.policy,report['offsets'])
        return dict(result=str(out/'result.json'),result_sha256=sha(out/'result.json'),independent=independent_result,
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
        for name,mutate,signature,options in [
            ('positive-wall-nu',lambda t:set_parameter(t,'Kinvis',1.5001e-5),True,{}),
            ('valid-other-wall-policy',other_policy,True,{}),
            ('valid-refined-mesh',None,False,dict(refined=True)),
            ('unsupported-Cd',lambda t:set_parameter(t,'RANSWJCd',.12),False,{}),
            ('gradient-policy',lambda t:set_info(t,'RANSWJGradientPolicy','Incompressible'),False,{}),
            ('epsilon-semantics',eps_semantics,False,{}),
            ('legacy-model',lambda t:set_info(t,'TurbulenceModel','CLSHighRe'),False,{}),
            ('sampling',lambda t:set_info(t,'RANSWallKSampling','Point'),False,{}),
            ('physical-cap',lambda t:set_info(t,'RANSWallDestruction','Capped'),False,{}),
            ('transition',lambda t:set_info(t,'RANSWallTransition','Upstream'),False,{}),
            ('curved',lambda t:set_info(t,'RANSWallMap','Cylinder'),False,{}),
            ('evaluation-reuse',lambda t:set_info(t,'RANSWallEvaluationReuse','ExactStepV1'),False,{})]:
            report,rc,folder=native('incompatible-'+name,0,cp1,mutate=mutate,**options)
            text=(folder/'solver.log').read_text()
            require(rc!=0 or (report and report.get('failure')),'incompatible input accepted '+name)
            if signature:require('signature mismatch' in text,'valid changed configuration did not reach signature validation: '+name)
            failures[name]=dict(exit_code=rc,signature_mismatch='signature mismatch' in text,log_sha256=sha(folder/'solver.log'))
        summary['incompatible']=failures
        if a.policy=='Budget':
            report,rc,folder=native('smooth-negative-reference',mutate=lambda t:set_parameter(t,'RANSWallRoughness',0))
            text=(folder/'solver.log').read_text();require(rc!=0 and 'reference_nonpositive' in text and '-0.330127' in text,
                'smooth negative reference was hidden or misclassified')
            summary['smooth_reference_nonpositive']=dict(log=str(folder/'solver.log'),sha256=sha(folder/'solver.log'))
        for name,nu in [('cutoff-failure',1.),('invalid-nu',0.)]:
            report,rc,folder=native(name,mutate=lambda t,nu=nu:set_parameter(t,'Kinvis',nu))
            require(rc!=0,'invalid cutoff/scalar case accepted')
            summary[name]=dict(log=str(folder/'solver.log'),sha256=sha(folder/'solver.log'))
        require(immutable0=={str(f.relative_to(cp0)):sha(f) for f in cp0.rglob('*') if f.is_file()} and sha(cp1/'State.bin')==immutable1,
                'saved checkpoint bytes changed')
        summary.update(status='PASS',reference_unchanged=True,source_revision=subprocess.check_output(
            ['git','-C',str(Path(__file__).resolve().parents[4]),'rev-parse','HEAD'],text=True).strip())
        summary['native_json_sha256']={str(f.relative_to(root)):sha(f) for f in root.rglob('*.json') if f.name!='summary.json'}
        save();print(root/'summary.json')
    except Exception as error:
        summary.update(status='FAIL',failure=str(error));save();raise


if __name__=='__main__':main()
