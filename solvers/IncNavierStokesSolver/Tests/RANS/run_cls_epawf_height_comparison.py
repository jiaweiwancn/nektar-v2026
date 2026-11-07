#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Declared H/tail/K/partition experiment. Native h0 and independent fixed dt."""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from cls_epawf_m4_cases import tensor_mesh
from cls_epawf_m8_r1_cases import X
from run_cls_epawf_picard import ROOT,THREADS,load,write_json,require,state_hash,sha,set_info,set_parameter,write
from run_cls_epawf_pimple import prepare,check_report,identity
from continue_cls_epawf_pimple import freeze,header,CampaignLock,atomic,process_identity

KN= [0.,.01,.1,.25,.5,1.,2.]+[2+j/3 for j in range(1,7)]
CASES=dict(A=(2.,'graded','Point',1.),B=(2.,'graded','Point',2.),C=(2.,'graded','ColumnAverage',2.),D=(2.,'one','ColumnAverage',2.),
           E=(1.,'graded','ColumnAverage',1.),F=(1.,'one','ColumnAverage',1.),G=(.5,'graded','ColumnAverage',.5),H=(.5,'one','ColumnAverage',.5))
UPSTREAM_CASES={'U05-graded':CASES['G'],'U05-one':CASES['H'],'U20-graded':CASES['C'],'U20-one':CASES['D']}
UPSTREAM_BRIDGES={'U20-point':(2.,'graded','Point',2.)}
WORK=dict(add=2000,total=20000,separations=2000,seconds=900)


def case_xml(case,folder,steps,dt=.000625,restart=None,normal=16,hard_reject=0,outer_refined=False,
             inventory='ConstitutiveV2',wall_policy=None):
    upstream=case in UPSTREAM_CASES or case in UPSTREAM_BRIDGES
    H,partition,sampler,tail=({**UPSTREAM_CASES,**UPSTREAM_BRIDGES} if upstream else CASES)[case]
    wall_policy=wall_policy or ('Upstream' if upstream else 'UniformEPFlat')
    if upstream:
        require(wall_policy=='Upstream' and inventory=='ConstitutiveV3','new Upstream family requires branch-complete inventory')
        require(tail==H,'every new Upstream physical case requires z_m=H')
    xml=prepare(folder,dt,steps,True,restart=restart,checkpoint=True,k_bound='MassProjection',floor_relative=1e-8,
        snapshots=True,inventory=inventory,wall_policy=wall_policy,bound_work=WORK,hard_reject_step=hard_reject)
    root=ET.parse(xml).getroot()
    knots=KN if partition=='graded' else [z for z in KN if z==0 or z>=H]
    require(not outer_refined or case in ('G','H','U05-graded','U05-one'),'focused outer refinement is declared only for the H=.5 pair')
    if outer_refined:knots=sorted(set(knots+[.75,1.5]))
    if partition=='one' or outer_refined:tensor_mesh(root,knots,3,2,6.,1.,points=[8,8],xknots=X)
    for key,value in dict(RANSWallHeight=H,RANSWallSampleHeight=1. if sampler=='Point' else H/2,RANSWallTailStart=tail,
        RANSWallNormalOrder=normal,IO_CheckSteps=8).items():set_parameter(root,key,value)
    for key,value in dict(RANSWallKSampling=sampler,RANSWallComparison='HeightV1',RANSWallComparisonCase=case).items():set_info(root,key,value)
    # Preserve the actual zero-perturbation initialization path. Update only
    # its geometric helper in the smaller-H outer function, still identically zero.
    if H!=2:
        for e in root.findall('./CONDITIONS/FUNCTION[@NAME="InitialPerturbationOuter"]/E'):
            if e.get('VAR')=='eps_tilde':e.set('VALUE',e.get('VALUE').replace('(y-(2.0))/(2.0)','(y-RANSWallHeight)/(4-RANSWallHeight)'))
    write(root,xml)
    family=('UpstreamPointBridgeV1' if case in UPSTREAM_BRIDGES else
        'UpstreamHeightV1'+('OuterHalves' if outer_refined else '') if upstream else
        'HeightSamplerOuterHalvesV1' if outer_refined else 'HeightSamplerV1')
    write_json(folder/'definition.json',dict(family=family,
        case=case,wall_policy=wall_policy,inventory=inventory,H=H,h=.1,Z=4,z_m=tail,
        temporal_coupling='PimpleSingleSweep',native_order=1,relaxation=1,eta=0,prescribed_forcing=0,k_floor_relative=1e-8,
        K_sampling=sampler,z_K=1. if sampler=='Point' else None,placeholder=H/2 if sampler!='Point' else None,
        K_functional_version='ColumnAverage-v1:physical-dz/H:fixed-native-panel-Gauss8:global-CG-row' if sampler=='ColumnAverage' else 'Point',
        transition=dict(start=-3,end=-2.25,baseline_start=-2,sensor_heights=[2.5,3.5],chi_enable=1),
        radial_knots=knots,tangential_knots=X,velocity_scalar_degree=3,pressure_degree=1,source_normal_order=normal,
        work=WORK,dt=dt,steps=steps,xml_sha256=sha(xml),zero_perturbation_functions_retained=True))
    return xml


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--solver',type=Path,required=True);p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--cases',nargs='+',choices=list(CASES)+list(UPSTREAM_CASES)+list(UPSTREAM_BRIDGES),default=list(CASES))
    p.add_argument('--inventory',choices=('ConstitutiveV2','ConstitutiveV3'),default='ConstitutiveV2')
    p.add_argument('--workers',type=int,default=4);p.add_argument('--steps',type=int,default=32)
    p.add_argument('--dt',type=float,default=.000625);p.add_argument('--normal-order',type=int,default=16)
    p.add_argument('--outer-refined',action='store_true',help='G/H only: add outer knots .75 and 1.5; retain all knots above 2')
    p.add_argument('--resume',action='store_true');a=p.parse_args();require(sys.platform=='linux','Linux authoritative')
    require(1<=a.workers<=4 and a.steps>=0,'bounded suite controls');a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=a.resume)
    os.environ.update(THREADS);(a.output/'tmp').mkdir(exist_ok=True);os.environ['TMPDIR']=str(a.output/'tmp')
    lock=CampaignLock(a.output)
    if not a.resume:
        solver,_=freeze(a.solver.resolve(),a.output/'runtime');probe,_=freeze(a.probe.resolve(),a.output/'probe-runtime')
    else:solver=a.output/'runtime'/a.solver.name;probe=a.output/'probe-runtime'/a.probe.name
    # Both frozen closures use the same dependencies; loader path is explicit.
    os.environ['LD_LIBRARY_PATH']=str(a.output/'runtime/lib')
    summary=dict(source=subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),
        source_status=subprocess.check_output(['git','-C',str(ROOT),'status','--porcelain'],text=True).strip(),
        solver=identity(solver),probe=identity(probe),status='running',workers=a.workers,ranks_per_case=1,threads=THREADS,cases={})
    atomic(a.output/'progress.json',summary)
    def native(exe,xml,folder,extra=(),timeout=7200):
        before=time.monotonic();command=[str(exe),str(xml),*map(str,extra)]
        with (folder/'solver.log').open('w') as log:
            child=subprocess.Popen(command,cwd=folder,stdout=log,stderr=subprocess.STDOUT,env=os.environ.copy(),pass_fds=(lock.stream.fileno(),))
            atomic(folder/'process.json',dict(identity=process_identity(child.pid),command=command))
            try:code=child.wait(timeout=timeout)
            except subprocess.TimeoutExpired:child.terminate();child.wait(timeout=30);code=124
        receipt=dict(command=command,exit_code=code,seconds=time.monotonic()-before,log_sha256=sha(folder/'solver.log'))
        write_json(folder/'receipt.json',receipt);return receipt
    def run_case(case):
        root=a.output/case
        if a.resume and (root/'result.json').exists():return load(root/'result.json')
        require(not root.exists(),'incomplete case requires actual process/checkpoint reconciliation, not an unchanged retry: '+case)
        root.mkdir();init=root/'initial';xml=case_xml(case,init,0,a.dt,normal=a.normal_order,outer_refined=a.outer_refined,inventory=a.inventory)
        receipt=native(solver,xml,init);report=load(init/'case.evolution.json')
        result=dict(case=case,definition=load(init/'definition.json'),initial_receipt=receipt,status='initialization_rejected',initial_failure=report.get('failure'))
        if report.get('initialization_complete') is False or receipt['exit_code']:
            write_json(root/'result.json',result);return result
        check_report(report,True);require(report['initial_k_bound_check']['activations']==0,'initial state was repaired')
        result.update(actual_h0_sha256=state_hash(report['state']),native_h0_checkpoint=header(Path(report['checkpoint'])),
            initialization=report['height_comparison_final'],initial_native_history_hashes=[state_hash(v) for v in report['native_history']])
        checks=root/'sampling';checks.mkdir();write_json(checks/'input.json',dict(operation='sampling_checks',state=report['state'],history=report['state'],time=0.,alpha=0.,require_history_h0=False,inner_stress_model='SaturatedR2'))
        test=native(probe,xml,checks,('input.json','result.json'),1800);require(test['exit_code']==0,'sampling verification failed: '+case)
        result['sampling_checks']=load(checks/'result.json')
        trajectory=root/'trajectory';xml=case_xml(case,trajectory,a.steps,a.dt,restart=Path(report['checkpoint']),normal=a.normal_order,outer_refined=a.outer_refined,inventory=a.inventory)
        receipt=native(solver,xml,trajectory,timeout=14400);e=load(trajectory/'case.evolution.json')
        check_report(e,True);result.update(status='completed' if not e['failure'] and e['steps']==a.steps else 'numerically_rejected',
            failure=e['failure'],steps=e['steps'],time=e['time'],trajectory_receipt=receipt,final=e.get('height_comparison_final'),
            accepted_checkpoint=header(Path(e['checkpoint'])),qualification_events=[{k:v for k,v in t.items() if k not in ('state','history','audit')} for t in e['physical_range_events']],
            observations={})
        require(state_hash(e['accepted_states'][0])==result['actual_h0_sha256'],'restart modified actual h0')
        for step in (8,16,32):
            if step>e['steps']:continue
            folder=root/f'probe-{step:04d}';folder.mkdir();state=e['accepted_states'][step]
            write_json(folder/'input.json',dict(operation='height_comparison',state=state,history=state,time=step*a.dt,alpha=0.,require_history_h0=False,inner_stress_model='SaturatedR2'))
            rc=native(probe,xml,folder,('input.json','result.json'),1800)
            # Optional readback cannot replace the actual native outcome.
            result['observations'][str(step)]=dict(exit_code=rc['exit_code'],path=str(folder/'result.json'),
                sha256=sha(folder/'result.json') if (folder/'result.json').exists() else None)
        result['attempts']=[{k:v for k,v in row.items() if k not in ('input_domain','k_domain')} for row in e['pimple_attempts']]
        write_json(root/'result.json',result);return result
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        jobs={pool.submit(run_case,c):c for c in a.cases}
        for f in as_completed(jobs):
            name=jobs[f]
            try:summary['cases'][name]=f.result()
            except BaseException as error:summary['cases'][name]=dict(status='infrastructure_or_verification_failure',failure=repr(error))
            atomic(a.output/'progress.json',summary)
    summary['status']='suite_completed';atomic(a.output/'progress.json',summary);print(a.output/'progress.json')

if __name__=='__main__':main()
