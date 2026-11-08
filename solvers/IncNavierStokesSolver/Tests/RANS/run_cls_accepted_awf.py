#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Fresh flat native cases for the accepted CLS-calibrated budget AWF.

This family changes the wall formulation, not merely chi. Legacy Upstream and
Uniform fixtures and checkpoints are never used as evolving initial states.
"""
import argparse
import math
import os
from pathlib import Path
import platform
import subprocess
import time
import xml.etree.ElementTree as ET

from cls_epawf_m4_cases import tensor_mesh
from cls_epawf_m8_r1_cases import X
from run_cls_epawf_picard import (ROOT,ORIGINAL,THREADS,load,write_json,require,
    state_hash,sha,set_info,set_parameter,write,session,run)
from run_cls_epawf_pimple import BOUND_WORK,identity
from continue_cls_epawf_pimple import freeze,header,CampaignLock,atomic

FAMILY='AcceptedAWFV1'
FORMULATION='CLSCalibratedBudgetAWF'
INVENTORY='AcceptedBudgetV1'
SAMPLER='InterfaceTrace'
FUNCTIONAL_VERSION='InterfaceTrace-v1:common-CG-matching-face:physical-local-H'
CASES={'Accepted-AWF-2G':(2.,'graded'),'Accepted-AWF-05G':(.5,'graded'),
       'Accepted-AWF-05O':(.5,'one'),'Accepted-AWF-2O':(2.,'one')}
COMPANIONS=('faithful','regularized','saturated','regularized-saturated')
REFINEMENTS=('none','cutoff-v1')
KN=[0.,.01,.1,.25,.5,1.,2.]+[2+j/3 for j in range(1,7)]
WORK=dict(add=2000,total=20000,separations=2000,seconds=900)
SPEC=ROOT/'docs/literature/accepted_awf_appendix.tex'
SPEC_SHA256='9998548d33bd2ca856bb6cc85551c82b86d3b499fd54ae8dc9ba26453450fdbf'


def model_choices(companion='faithful'):
    require(companion in COMPANIONS,'unknown explicitly named companion')
    regularized=companion in ('regularized','regularized-saturated')
    saturated=companion in ('saturated','regularized-saturated')
    return dict(outer_variant='Regularized' if regularized else 'Original',
        inner_stress='SaturatedR2' if saturated else 'Linear',
        cubic_beta=.4 if regularized else 1.,regularization_gamma=0.,
        regularization_scale=1. if regularized else 0.)


def case_xml(case,folder,steps,dt=.000625,restart=None,normal=16,hard_reject=0,
             companion='faithful',companion_reason='',perturbation=0.,verify=False,inner_refinement='none'):
    require(case in CASES and int(steps)==steps and steps>=0 and dt>0,'invalid accepted-AWF case/step controls')
    require(math.isfinite(perturbation) and abs(perturbation)<=.05,'small admissible perturbation amplitude required')
    require(companion=='faithful' or companion_reason.strip(),'a companion must identify the preserved faithful-model blocker')
    require(sha(SPEC)==SPEC_SHA256,'approved accepted-AWF appendix differs; version the changed specification explicitly')
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=False)
    xml=session(ORIGINAL,folder/'case.xml','outer6')
    root=ET.parse(xml).getroot();H,partition=CASES[case];choices=model_choices(companion)
    knots=KN if partition=='graded' else [z for z in KN if z==0 or z>=H]
    require(inner_refinement in REFINEMENTS,'unknown focused inner refinement')
    require(inner_refinement=='none' or partition=='graded','cutoff-v1 is a separately named graded follow-up')
    added=[z for z in (.0025,.005,.02,.04) if z<H] if inner_refinement=='cutoff-v1' else []
    if added:knots=sorted(set(knots+added))
    nx=len(X)-1;ni=sum(z<=H for z in knots[1:]);no=len(knots)-1-ni
    sizes=[(3*nx+1)*(3*(ni+no)+1)]*3+[(3*nx+1)*(3*no+1),(nx+1)*(ni+1),(nx+1)*(no+1)]
    require(sum(sizes)+1<=3000,'declared case exceeds the serial dense pilot size limit')
    tensor_mesh(root,knots,3,2,6.,1.,points=[8,8],xknots=X)
    root.find('./CONDITIONS/TIMEINTEGRATIONSCHEME/ORDER').text='1'
    for key,value in {
        'RANSWallFormulation':FORMULATION,'RANSWallTransition':'Fixed',
        'RANSWallTemporalCoupling':'PimpleSingleSweep','RANSWallNonlinearPolicy':'Legacy',
        'RANSWallPhysicalRangeGuard':'AcceptedAWFDiagnosticV1','RANSWallVariant':choices['outer_variant'],
        'RANSWallInnerStress':choices['inner_stress'],'RANSWallDestruction':'Original',
        'RANSWallKSampling':SAMPLER,'RANSWallReference':'Neutral','RANSWallMap':'Flat',
        'RANSWallComparison':FAMILY,'RANSWallComparisonCase':case,
        'RANSWallPimpleKBound':'MassProjection','RANSWallPimpleKInventory':INVENTORY,
        'RANSWallPimpleKBoundSolver':'RetainedPool','RANSWallPimpleSnapshots':'On',
        'RANSWallPimpleVerify':'On' if verify else 'Off','RANSWallSolveDiagnostics':'Off',
        'RANSWallRestart':str(Path(restart).resolve()) if restart else 'None',
        'RANSWallCheckpoint':str(folder.resolve()/'checkpoint'),
        'RANSWallPicardReduced':'Off','RANSWallPicardGroundProjection':'Off',
        'RANSWallPicardKNewton':'Off','RANSWallPicardEquilibrateK':'Off',
        'RANSWallPicardContinuation':'Off'}.items():set_info(root,key,value)
    for key,value in {
        'NumSteps':int(steps),'TimeStep':dt,'Time':0.,'RANSWallHeight':H,
        'RANSWallSampleHeight':H,'RANSWallTailStart':H,'RANSWallRoughness':.1,'Kinvis':1.5e-5,
        'RANSWallFrictionVelocity':.4,'RANSWallEta':0.,'RANSWallChi':0.,
        'RANSWallPicardFlowRelaxation':1.,'RANSWallPicardKRelaxation':1.,
        'RANSWallPicardEpsilonRelaxation':1.,'RANSWallPicardPseudoAlpha':0.,
        'RANSWallCubicBeta':choices['cubic_beta'],'RANSWallRegularizationGamma':choices['regularization_gamma'],
        'RANSWallRegularizationScale':choices['regularization_scale'],
        'RANSWallPimpleKFloorRelative':1e-8,'RANSWallNormalOrder':normal,'RANSWallTangentialOrder':8,
        'RANSWallAcceptedKAmplitude':perturbation,'RANSWallPimpleHardAuditRejectStep':hard_reject,
        'RANSWallPimpleRejectStep':0,'RANSWallDiagnostics':1,'IO_CheckSteps':8,
        'Slope':'Kappa*sqrt(K0/Qstar)',
        'E0':'(K0/Qstar)^1.5/(Kappa*RANSWallHeight)',
        'YapScale':'K0^1.5/(2.5*E0)'}.items():set_parameter(root,key,value)
    for key,value in WORK.items():set_parameter(root,BOUND_WORK[key],value)
    # These archived expressions are identically zero, and the native solver
    # deliberately preserves its zero-perturbation initialization fast path.
    # The separately selected accepted-k amplitude is implemented natively.
    if H!=2:
        for e in root.findall('./CONDITIONS/FUNCTION[@NAME="InitialPerturbationOuter"]/E'):
            if e.get('VAR')=='eps_tilde':e.set('VALUE',e.get('VALUE').replace('(y-(2.0))/(2.0)','(y-RANSWallHeight)/(4-RANSWallHeight)'))
    write(root,xml)
    definition=dict(family=FAMILY,case=case,formulation=FORMULATION,formulation_version='accepted-appendix-v1',
        case_identity=case+(':'+inner_refinement if inner_refinement!='none' else ''),inner_refinement=inner_refinement,
        focused_refinement_knots=added,
        appendix_sha256=SPEC_SHA256,wall_policy=FORMULATION,companion=companion,companion_reason=companion_reason,
        model_choices=choices,H=H,z_K=H,z_m=H,h=.1,Z=4.,molecular_nu=1.5e-5,friction_velocity=.4,
        K_sampling=SAMPLER,K_functional_version=FUNCTIONAL_VERSION,source='P_h-E_C(K); prescribed redistribution exactly zero',
        scalar_source_split='frozen K/viscosity/cutoff/E_C; signed source; no local reaction rate',
        outer_equations=choices['outer_variant']+' high-Re CLS/Yap; no molecular outer diffusion',inner_viscosity='nu+nu_tA',
        wall_coefficients='derived from selected CLS q_star/kappa: c_mu=q_star^-2; c_length=kappa*q_star^1.5',
        chi=None,capped_normalization=None,reference='native projection of nonconstant diffusion-balanced k_w,0',
        temporal_coupling='PimpleSingleSweep',native_order=1,dt=dt,steps=int(steps),relaxation=1.,eta=0.,prescribed_forcing=0.,
        radial_knots=knots,tangential_knots=X,velocity_scalar_degree=3,pressure_degree=1,
        tensor_unknown_estimate=sum(sizes)+1,tensor_field_size_estimates=sizes,
        native_points=[8,8],source_normal_order=normal,source_tangential_order=8,
        inventory=INVENTORY,bound_solver='RetainedPool',k_floor_relative=1e-8,work=WORK,
        perturbation=dict(amplitude=perturbation,scale='K0',
            shape='sin(pi*(x-xmin)/(xmax-xmin))^2*(1-3*(z/H)^2+2*(z/H)^3) on inner layer; zero outside',
            endpoint_and_external_compatible=True,maximum_relative_amplitude=abs(perturbation)),
        xml_sha256=sha(xml),fresh_native_initialization=restart is None,
        restart=str(Path(restart).resolve()) if restart else None,zero_prescribed_perturbation_functions_retained=True)
    write_json(folder/'definition.json',definition)
    return xml


def check_report(report,restarting=False,definition=None):
    """Numerical completion is distinct from source qualification and rejection.

    Failed attempts may complete all governing solves and still fail their
    projection/audit. Keep those outcomes rather than asserting acceptance.
    """
    require(report['temporal_coupling']=='PimpleSingleSweep' and report['order']==1,'wrong native temporal scheme')
    require(report['wall_policy']==FORMULATION and report['eta']==0.,'wrong wall formulation or anchoring')
    require(report['k_bound_mode']=='MassProjection','missing declared final correction')
    if definition:
        expected=definition['model_choices'];p=report['parameters'];r=report['robust_parameters']
        require(report['inner_stress_model']==expected['inner_stress'],'unexpected inner stress')
        require(report['model_variant']==expected['outer_variant'],'unexpected outer constitutive variant')
        require(r[0]==expected['cubic_beta'] and r[1]==expected['regularization_gamma'] and
            r[5]==expected['regularization_scale'] and r[6]==0.,'changed outer regularization or capped destruction')
        require(p[0]==definition['molecular_nu'] and p[1]==p[3]==p[4]==definition['H'] and
            p[2]==definition['h'],'changed physical geometry or endpoint/tail semantics')
    require(report['inner_stress_model'] in ('Linear','SaturatedR2'),'unreported inner stress model')
    attempts=report['pimple_attempts'];require(report['implicit_callback_count']==len(attempts),'implicit-callback/attempt mismatch')
    accepted_times={x[0] for x in report.get('stage_audits',[])}
    for t in attempts:
        t['numerically_accepted']=t['time'] in accepted_times
        require(t.get('additional_governing_equation_solves',0)==0,'hidden governing-equation solve')
        for name in ('outer_sweeps','momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves','bound_calls'):
            require(t.get(name,0) in (0,1),'extra scheduled pass: '+name)
        if not t['numerically_accepted']:continue
        require(t['scheduled_complete'],'accepted step lacks the scheduled pass')
        for name in ('outer_sweeps','momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves','bound_calls'):
            require(t[name]==1,'accepted step missing scheduled work: '+name)
        require(t['mass_inverse_applications']==t['pressure_gauge_size'],'unexpected pressure correction work')
        require(t['additional_correction_or_restoration_solves']==t['bound_kkt_solves'],'unreported algebraic correction work')
        b=t['k_bound'];require(b['work_controls'].get('finite_qp_solver')=='RetainedPool','wrong finite-QP solver')
        if b['activations']:
            require(b['success'] and b['finite_pool_complete'] and b['primal_floor_defect']==0 and b['ground_floor_defect']==0,
                'accepted correction lacks complete finite/physical feasibility')
        require(max(t['corrected_continuity_gauge'])<1e-10 and t['boundary_defect']<1e-13 and
            t['momentum_splitting_identity_defect']<1e-9,'flow constraints or splitting identity failed')
        require(t['imposed_trace_defect']<1e-12 and t['trace_mass_column_max']>0 and t['retained_trace_column_max']>0,
            'single frozen epsilon trace target/columns failed')
        require(t['epsilon_frozen_linear_residual']<1e-9 and t['k_frozen_linear_residual']<1e-9,'frozen scalar solve failed')
    return attempts


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--solver',type=Path)
    p.add_argument('--probe',type=Path);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--cases',choices=tuple(CASES),nargs='+',default=['Accepted-AWF-2G'])
    p.add_argument('--companion',choices=COMPANIONS,default='faithful');p.add_argument('--companion-reason',default='')
    p.add_argument('--steps',type=int,default=32);p.add_argument('--dt',type=float,default=.000625)
    p.add_argument('--normal-order',type=int,default=16);p.add_argument('--perturbation',type=float,default=0.)
    p.add_argument('--inner-refinement',choices=REFINEMENTS,default='none')
    p.add_argument('--timeout',type=float,default=14400);p.add_argument('--prepare-only',action='store_true')
    p.add_argument('--verify',action='store_true');a=p.parse_args();os.environ.update(THREADS)
    require(a.steps>=0 and a.dt>0 and a.timeout>0,'invalid fixed-step controls')
    a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False)
    summary=dict(schema=1,family=FAMILY,formulation=FORMULATION,specification_sha256=SPEC_SHA256,
        source_revision=subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),
        source_status=subprocess.check_output(['git','-C',str(ROOT),'status','--porcelain'],text=True).strip(),
        ranks_per_runtime=1,threads=THREADS,dt=a.dt,requested_steps=a.steps,status='prepared',cases={})
    if a.prepare_only:
        for case in a.cases:
            xml=case_xml(case,a.output/case,a.steps,a.dt,normal=a.normal_order,companion=a.companion,
                companion_reason=a.companion_reason,perturbation=a.perturbation,verify=a.verify,inner_refinement=a.inner_refinement)
            summary['cases'][case]=load(xml.parent/'definition.json')
        write_json(a.output/'manifest.json',summary);return
    require(platform.system()=='Linux' and a.solver,'authoritative native execution requires Linux and --solver')
    require(a.output.is_relative_to(Path('/media/jay/新加卷/nektar-v2026-runs')),'heavy run storage')
    lock=CampaignLock(a.output);(a.output/'tmp').mkdir();os.environ['TMPDIR']=str(a.output/'tmp')
    solver,solver_id=freeze(a.solver.resolve(),a.output/'runtime');probe=None
    if a.probe:probe,_=freeze(a.probe.resolve(),a.output/'probe-runtime')
    os.environ['LD_LIBRARY_PATH']=str(a.output/'runtime/lib');summary.update(solver=solver_id,status='running')
    atomic(a.output/'manifest.json',summary)
    def native(executable,xml,folder,extra=(),timeout=None):
        started=time.monotonic();command=[str(executable),str(xml),*map(str,extra)]
        with (folder/'solver.log').open('w') as stream:
            child=subprocess.Popen(command,cwd=folder,env=os.environ.copy(),stdout=stream,
                stderr=subprocess.STDOUT,pass_fds=(lock.stream.fileno(),))
            atomic(folder/'process.json',dict(pid=child.pid,command=command))
            try:code=child.wait(timeout=timeout or a.timeout)
            except subprocess.TimeoutExpired:child.terminate();child.wait(timeout=30);code=124
        receipt=dict(command=command,exit_code=code,seconds=time.monotonic()-started,log_sha256=sha(folder/'solver.log'))
        write_json(folder/'receipt.json',receipt);return receipt
    for case in a.cases:
        base=a.output/case;base.mkdir();initial=base/'initial'
        xml=case_xml(case,initial,0,a.dt,normal=a.normal_order,companion=a.companion,
            companion_reason=a.companion_reason,perturbation=a.perturbation,inner_refinement=a.inner_refinement)
        receipt=native(solver,xml,initial);reportfile=initial/'case.evolution.json'
        result=dict(case=case,definition=load(initial/'definition.json'),initial_receipt=receipt,status='initialization_rejected')
        summary['cases'][case]=result
        if not reportfile.exists():
            result.update(status='native_exit_without_report',failure='inspect preserved initial solver.log')
            atomic(base/'result.json',result);atomic(a.output/'manifest.json',summary);continue
        report=load(reportfile);result.update(initial_failure=report.get('failure'),
            initial_native_report_sha256=sha(reportfile),initialization=report.get('accepted_awf_final'),
            initialization_complete=report.get('initialization_complete'),
            initial_candidate_sha256=state_hash(report['state']) if report.get('state') else None,
            initial_candidate_history_hashes=[state_hash(x) for x in report.get('native_history',[])])
        if receipt['exit_code'] or report.get('initialization_complete') is False or report.get('failure'):
            result['failure']=report.get('failure');atomic(base/'result.json',result);atomic(a.output/'manifest.json',summary);continue
        check_report(report,False,result['definition'])
        require('accepted_awf_final' in report,'native model-aware initialization readback missing')
        require(not report['pimple_attempts'] and report['implicit_callback_count']==0,'initialization performed a physical step')
        require(report['initial_k_bound_check']['activations']==0,'native h0 was changed by startup repair')
        cp=Path(report['checkpoint']);h0=state_hash(report['state'])
        result.update(actual_h0_sha256=h0,initial_native_history_hashes=[state_hash(x) for x in report['native_history']],
            initial_native_history_times=report['native_history_times'],native_h0_checkpoint=header(cp),
            initialization=report.get('accepted_awf_final'),initialization_work=report['initialization_lift_calls_iterations'])
        if not a.steps:
            result.update(status='initialized',steps=0,time=0.,accepted_checkpoint=header(cp))
            atomic(base/'result.json',result);atomic(a.output/'manifest.json',summary);continue
        trajectory=base/'trajectory';txml=case_xml(case,trajectory,a.steps,a.dt,restart=cp,
            normal=a.normal_order,companion=a.companion,companion_reason=a.companion_reason,
            perturbation=a.perturbation,verify=a.verify,inner_refinement=a.inner_refinement)
        rc=native(solver,txml,trajectory);nativefile=trajectory/'case.evolution.json'
        result['trajectory_receipt']=rc
        if not nativefile.exists():
            result.update(status='infrastructure_interrupted',failure='native exit without durable final report; inspect checkpoint/process evidence')
            atomic(base/'result.json',result);atomic(a.output/'manifest.json',summary);continue
        d=load(nativefile);attempts=check_report(d,True,result['definition'])
        require(state_hash(d['accepted_states'][0])==h0,'compatible restart changed actual native h0')
        result.update(status='infrastructure_interrupted' if rc['exit_code']==124 else
            ('numerical_trajectory_completed' if rc['exit_code']==0 and not d['failure'] else 'numerically_rejected'),
            steps=d['steps'],time=d['time'],failure=d['failure'],final=d.get('accepted_awf_final'),
            accepted_checkpoint=header(Path(d['checkpoint'])),implicit_callbacks=d['implicit_callback_count'],
            correction_integral=sum(t.get('k_bound',{}).get('integral_total',0.) for t in attempts if t['numerically_accepted']),
            snapshots={})
        if rc['exit_code']==0:require(d['steps']==a.steps,'successful native trajectory missed the fixed horizon')
        if probe:
            for step in (8,16,32,80,160):
                if step>d['steps']:continue
                folder=base/f'observe-{step:06d}';folder.mkdir();state=d['accepted_states'][step]
                write_json(folder/'input.json',dict(operation='accepted_awf',state=state,history=state,
                    alpha=0.,time=step*a.dt,require_history_h0=False,
                    inner_stress_model=result['definition']['model_choices']['inner_stress']))
                read=native(probe,txml,folder,('input.json','result.json'),1800)
                result['snapshots'][str(step)]=dict(time=step*a.dt,exit_code=read['exit_code'],
                    path=str(folder/'result.json'),sha256=sha(folder/'result.json') if (folder/'result.json').exists() else None,
                    interpretation='read-only optional query; cannot replace the native trajectory outcome')
        atomic(base/'result.json',result);atomic(a.output/'manifest.json',summary)
    summary['status']='suite_recorded';atomic(a.output/'manifest.json',summary);print(a.output/'manifest.json')


if __name__=='__main__':main()
