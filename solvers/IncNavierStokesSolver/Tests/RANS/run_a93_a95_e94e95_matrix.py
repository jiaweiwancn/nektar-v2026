#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded serial native matrix; keeps every rejected case and raw report."""
import argparse
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import xml.etree.ElementTree as ET
sys.dont_write_bytecode=True
from rans_qualification import dependencies,coverage,audited_extrema
dependencies('numpy','scipy')
from run_cls_epawf_picard import ROOT,ORIGINAL,THREADS,session,set_info,set_parameter,write,sha,state_hash
from cls_epawf_m4_cases import tensor_mesh
from run_cls_epawf_pimple import identity
from validate_a93_a95_e94e95 import Probe


def case_xml(folder,model,policy,wall,steps=16,dt=.0001,restart=None,cache='Off',reject=0,refine=False,perturb=0.,descriptor_probe=None,
             radial_knots=None,degree=3,quadrature=None,checkpoint_interval=1):
    folder.mkdir(parents=True,exist_ok=False)
    xml=session(ORIGINAL,folder/'case.xml','outer6');root=ET.parse(xml).getroot()
    knots=[0.,.0025,.005,.01,.02,.04,.1,.25,.5,1.,2.,2.25,2.5,3.,3.5,4.]
    if refine:knots=sorted(set(knots+[(a+b)/2 for a,b in zip(knots[:-1],knots[1:])]))
    if radial_knots is not None:knots=radial_knots
    tensor_mesh(root,knots,degree,2,2.,1.,points=[max(8,degree+3)]*2,xknots=[-1.,1.])
    if descriptor_probe is None:raise ValueError('case construction requires the resolved descriptor probe')
    provider=Probe(descriptor_probe)
    resolved=provider.closure(model,[[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]],policy=policy)
    provider.close()
    if not resolved['ok']:raise ValueError(resolved)
    q,kap,sigma=resolved['q'],resolved['kappa'],resolved['sigma']
    root.find('./CONDITIONS/TIMEINTEGRATIONSCHEME/ORDER').text='1'
    infos={'TurbulenceModel':model,'RANSKappaPolicy':policy,'RANSWallFormulation':'EPAWF' if wall=='EP' else 'AcceptedBudgetAWF',
        'RANSWallTransition':'UniformEPFlat' if wall=='EP' else 'Fixed','RANSWallTemporalCoupling':'PimpleSingleSweep',
        'RANSWallNonlinearPolicy':'Legacy','RANSWallVariant':'Original','RANSWallInnerStress':'Linear','RANSWallDestruction':'Original',
        'RANSWallPhysicalRangeGuard':'NewModelDiagnosticV1','RANSWallReference':'Neutral','RANSWallMap':'Flat',
        'RANSWallKSampling':'Point' if wall=='EP' else 'InterfaceTrace','RANSWallComparison':'Off','RANSWallComparisonCase':'None',
        'RANSWallPimpleKBound':'Off','RANSWallPimpleKInventory':'LegacyV1','RANSWallPimpleSnapshots':'Off',
        'RANSWallEvaluationReuse':cache,'RANSWallPimpleVerify':'Off','RANSWallSolveDiagnostics':'Off',
        'RANSWallRestart':str(restart.resolve()) if restart else 'None','RANSWallCheckpoint':str(folder.resolve()/'checkpoint'),
        'RANSWallPicardGroundProjection':'Off','RANSWallPicardKNewton':'Off','RANSWallPicardReduced':'Off',
        'RANSWallPicardEquilibrateK':'Off','RANSWallPicardContinuation':'Off'}
    for key,value in infos.items():set_info(root,key,value)
    params={'NumSteps':steps,'TimeStep':dt,'Time':0.,'RANSWallHeight':2.,'RANSWallSampleHeight':2.,'RANSWallTailStart':2.,
        'RANSWallRoughness':.1,'Kinvis':1.5e-5,'RANSWallFrictionVelocity':.4,'RANSWallEta':0.,'RANSWallChi':1. if wall=='EP' else 0.,
        'RANSKappa':.41,'Qstar':q,'Kappa':kap,'K0':q*.16,'RANSWallNormalOrder':48 if refine else 32,'RANSWallTangentialOrder':8,
        'RANSWallDiagnostics':1,'IO_CheckSteps':checkpoint_interval,'RANSWallPimpleRejectStep':reject,'RANSWallPimpleHardAuditRejectStep':0,
        'RANSWallPicardFlowRelaxation':1.,'RANSWallPicardKRelaxation':1.,'RANSWallPicardEpsilonRelaxation':1.,
        'RANSWallPicardPseudoAlpha':0.,'RANSWallCubicBeta':1.,'RANSWallRegularizationGamma':0.,'RANSWallRegularizationScale':0.,
        'RANSWallAcceptedKAmplitude':perturb if wall=='Accepted' else 0.}
    for key,value in params.items():set_parameter(root,key,value)
    if quadrature:
        for key,value in zip(('RANSWallNormalOrder','RANSWallTangentialOrder','RANSWallMassOrder','RANSWallVolumeOrder','RANSWallInterfaceOrder'),quadrature,strict=True):set_parameter(root,key,value)
    # Initial perturbation functions of the CLS archive are not part of this
    # reference definition. A separately declared native inner-k perturbation
    # uses fixed external and interface traces when requested.
    for function in list(root.findall('./CONDITIONS/FUNCTION')):
        if function.get('NAME','').startswith('InitialPerturbation'):root.find('CONDITIONS').remove(function)
    if perturb and wall=='EP':
        function=ET.SubElement(root.find('CONDITIONS'),'FUNCTION',{'NAME':'InitialPerturbationInner'})
        ET.SubElement(function,'E',{'VAR':'k','VALUE':str(perturb)+'*K0*sin(3.141592653589793*(x+1)/2)^2*(1-3*(y/RANSWallHeight)^2+2*(y/RANSWallHeight)^3)'})
    write(root,xml)
    return xml


def run_case(solver,folder,model,policy,wall,**kwargs):
    xml=case_xml(folder,model,policy,wall,**kwargs)
    with (folder/'solver.log').open('w') as f:
        try:completed=subprocess.run([str(solver),'case.xml'],cwd=folder,env=dict(os.environ,**THREADS),stdout=f,stderr=subprocess.STDOUT,timeout=1800);code=completed.returncode
        except subprocess.TimeoutExpired:code=124
    row=dict(model=model,kappa_policy=policy,wall=wall,options={k:str(v) if isinstance(v,Path) else v for k,v in kwargs.items()},exit_code=code,
        xml_sha256=sha(xml),path=str(folder),status='fail',report=None)
    report=folder/'case.evolution.json'
    if report.exists():
        r=json.loads(report.read_text());row['report_sha256']=sha(report)
        row.update(accepted_steps=r.get('steps',0),time=r.get('time',0.),failure=r.get('failure','initialization incomplete'),initial_residual=r.get('raw_reference_norms'),
            final_stage_norms=r.get('stage_audits',[])[-1:] ,initial_state_sha256=state_hash((r.get('accepted_states') or [r.get('initial',[])])[0]),
            state_sha256=state_hash(r.get('state',[])),history=r.get('native_history'),physics_signature=r.get('signature'),
            state_audits=r.get('state_audits'),step_diagnostics=[{k:v for k,v in s.items() if k in ('time','refreshed_original_block_norms','refreshed_trace_defect',
                'corrected_continuity_gauge','momentum_splitting_identity_defect','epsilon_same_state_source_storage_defect','k_same_state_source_storage_defect','bound_calls','k_bound')} for s in r.get('pimple_attempts',[])])
        requested=kwargs.get('steps',16);start_steps=0;start_time=0.
        if kwargs.get('restart'):
            from continue_cls_epawf_pimple import header
            receipt=header(kwargs['restart']);start_steps=receipt['steps'];start_time=receipt['time']
        end=start_time+requested*kwargs.get('dt',.0001)
        complete=(code==0 and not r.get('failure') and r.get('steps')==start_steps+requested and
            abs(r.get('time',-1.)-end)<1e-13 and len(r.get('stage_audits',[]))==requested)
        if complete:row['status']='pass'
        row['execution_status']='completed' if complete else 'failed_or_incomplete'
        row['numerical_accuracy_status']='not_assessed_by_trajectory_runner';row['convergence_status']='not_assessed'
        row['audited_extrema']=audited_extrema(r)
        # Keep compact reports; full states and every native audit stay on Linux.
        row['drift_history']=r.get('step_audits',[])
        row['milestones']=[s for s in r.get('step_audits',[]) if int(s[1]) in (1,2,8,16)]
        audits=r.get('state_audits',[])
        row['initial_diagnostics']=audits[0] if audits else None
        row['final_diagnostics']=audits[-1] if audits else None
        row.pop('history',None);row.pop('state_audits',None)
    (folder/'result.json').write_text(json.dumps(row,indent=2)+'\n');return row


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--solver',type=Path,required=True);p.add_argument('--descriptor-probe',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--steps',type=int,default=16);p.add_argument('--models',nargs='+',default=['E94E95','A95','A93']);p.add_argument('--policies',nargs='+',default=['Natural','Prescribed']);p.add_argument('--walls',nargs='+',default=['EP','Accepted'])
    a=p.parse_args();assert platform.system()=='Linux';os.environ.update(THREADS);a.output.mkdir(parents=True,exist_ok=False)
    summary=dict(schema=1,source=subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),source_status=subprocess.check_output(['git','-c','core.filemode=false','-C',str(ROOT),'status','--porcelain'],text=True),runtime=identity(a.solver),threads=THREADS,ranks=1,cases=[])
    for model in a.models:
        for policy in a.policies:
            for wall in a.walls:
                folder=a.output/(model+'-'+policy+'-'+wall)
                row=run_case(a.solver,folder,model,policy,wall,steps=a.steps,descriptor_probe=a.descriptor_probe)
                summary['cases'].append(row);(a.output/'matrix.json').write_text(json.dumps(summary,indent=2)+'\n')
                print(json.dumps({k:row.get(k) for k in ('model','kappa_policy','wall','status','accepted_steps','failure')}),flush=True)
    expected={(m,p,w) for m in a.models for p in a.policies for w in a.walls}
    summary['coverage']=coverage(summary['cases'],expected)
    summary['status']='pass' if summary['coverage']['status']=='pass' and all(r['status']=='pass' for r in summary['cases']) else 'fail'
    summary['accuracy_status']='not_assessed';summary['convergence_status']='not_assessed'
    (a.output/'matrix.json').write_text(json.dumps(summary,indent=2)+'\n')
    return 0 if summary['status']=='pass' else 1
if __name__=='__main__':sys.exit(main())
