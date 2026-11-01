#!/usr/bin/env python3
"""Native CV3 inlet/startup verification, separate from calibration outcomes."""
import argparse
import json
import math
import os
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode=True
from cao_tamura_experiment import (ROOT, MODELS, THREADS, allocation, environment,
    execute, generate, prepare, records, run, sha, write_json)
from generate_m4 import info, parameter
from validate_terrain import compare
from cao_tamura_cv3_budgets import flat_momentum_budget
from terrain_continuity import audit as strip_continuity


def require(value, message):
    if not value:raise RuntimeError(message)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--cpus',required=True)
    args=parser.parse_args();require(sys.platform.startswith('linux'),'Linux is authoritative')
    build=environment(args.build);allocation(4,args.cpus)
    output=args.output.resolve();output.mkdir(parents=True,exist_ok=False)
    os.sched_setaffinity(0,set(map(int,args.cpus.split(','))))
    first=args.cpus.split(',')[0];results={}
    common=dict(geometry='floor',degree=3,pressure_degree=2,layers=6,outer_blend_cells=3,
        x_columns=[-1.,0.,1.2],inlet_delta=.21,velocity_scaling='incident_ReH',
        inlet_reference_q=10/3,inlet_reference_kappa=math.sqrt(.48*1.3/(10/3)),
        k_profile='MeasuredRms',epsilon_profile='ShearLength',inlet_wall_collar=True,
        dt=.005,checkpoint=1,steps=2)
    preparations={}
    for surface in ('S','R'):
        for model in ('standard','standard_original','standard_fixed','cls_original','cls_fixed'):
            folder=output/(surface+'-'+model)
            m=generate({**common,'surface':surface,'model':model,'ks':0 if surface=='S' else .004},folder)
            receipt=prepare(folder/'case.xml',build)
            preparations[surface+'-'+model]={k:receipt[k] for k in ('maximum_relative_input_error','minimum_projected','raw_interface_normal_max')}
            rows=records(folder/'initial.fld.csv')
            epsilon=[r for r in rows if r['field'] in ('eps','eps_tilde') and abs(float(r['distance'])-.011)<1e-10]
            require(epsilon and max(abs(float(r['input'])-m['reference']['wall_epsilon_near']) for r in epsilon)<5e-10,
                    'Native epsilon initializer differs from the policy trace')
    results['native_inlet_all_primary_policies']=preparations
    print('native inlet expressions PASS',flush=True)

    def solve(name,steps,duration,pressure=2):
        folder=output/name
        # Keep the requested duration in both generated XML and its sidecar.
        uinf=1.260222109045355
        m=generate({**common,'model':'standard_original','steps':steps,'advection':'Implicit',
                    'pressure_degree':pressure,'startup_flow_times':duration*uinf/2.2},folder)
        tree=ET.parse(folder/'case.xml');root=tree.getroot()
        info(root,'GlobalSysSoln','IterativeMultiLevelStaticCond');parameter(root,'RANSStageAudit',1)
        parameter(root,'RANSStartupDuration',duration)
        parameter(root,'RANSGradDiv',.01)
        tree.write(folder/'case.xml',encoding='utf-8',xml_declaration=True)
        m['case_sha256']=sha(folder/'case.xml');write_json(folder/'manifest.json',m)
        prep=prepare(folder/'case.xml',build);receipt=run(folder/'case.xml',build,1,first,600)
        residual=records(folder/'case.rans-residual.csv')
        free=max(float(r['free_coefficient_l2']) for r in residual if r['phase']=='1')
        require(free<1e-7 and max(receipt['stages'][k] for k in ('maximum_weak_divergence','maximum_projected_normal'))<1e-7,
                'Native startup stage residual/continuity failed')
        return folder,dict(run=receipt,preparation=prep,maximum_free_residual=free)

    whole,results['startup_whole']=solve('startup-whole',8,.01)
    prefix,results['startup_prefix']=solve('startup-prefix',4,.01)
    require(results['startup_whole']['run']['stages']['final_startup_weight']==0,'Startup did not end exactly')
    weights=records(whole/'case.rans.csv')
    require(all(float(r['startup_weight'])==0 for r in weights if float(r['time'])>=.01),'A startup load survived its duration')
    resumed=output/'startup-resumed';resumed.mkdir()
    helper=ROOT/'solvers/IncNavierStokesSolver/Examples/NeutralABL/restart.py'
    execute([sys.executable,helper,prefix/'case.xml',prefix/'case.fld',resumed/'case.xml','--steps','4'],resumed,'prepare-restart')
    m=json.loads((prefix/'manifest.json').read_text());m['case_sha256']=sha(resumed/'case.xml');write_json(resumed/'manifest.json',m)
    receipt=run(resumed/'case.xml',build,4,args.cpus,600)
    paths=[]
    for folder in (whole,resumed):
        path=folder/'full.csv'
        execute(allocation(1,first)+[build/'RANSWallFieldAudit',folder/'case.xml',folder/'case.fld',path],folder,'full')
        paths.append(path)
    errors=compare(*paths);require(max(errors.values())<1e-6,'Post-startup 1-to-4-rank history differs')
    results['post_startup_MPI_history']=dict(run=receipt,scaled_errors=errors)
    results['strip_momentum_inventory']={folder.name:flat_momentum_budget(folder/'case.xml') for folder in (whole,resumed)}
    rejected=output/'mid-startup-rejected';rejected.mkdir()
    execute([sys.executable,helper,prefix/'case.xml',prefix/'case_1.chk',rejected/'case.xml','--steps','1'],rejected,'prepare-restart')
    try:execute(allocation(1,first)+[build/'IncNavierStokesSolver','case.xml'],rejected,'reject')
    except RuntimeError:
        require('startup weak load is not checkpointed' in (rejected/'reject.log').read_text(),'Wrong mid-startup rejection')
    else:raise RuntimeError('Unsupported mid-startup restart was accepted')
    results['mid_startup_restart_rejected']=True
    for degree in (1,2):
        folder,result=solve('hold-pressure-'+str(degree),2,10.,degree)
        residual=records(folder/'case.rans-residual.csv')
        initial=max(float(r['free_coefficient_l2']) for r in residual if r['phase']=='0' and r['field']=='continuity')
        require(initial<1e-10,'Parallel inlet initial continuity audit has a boundary orientation error')
        result['initial_continuity_residual']=initial;results['hold-pressure-'+str(degree)]=result
    # This real precursor previously stopped at t=1.5 because its linear
    # Picard image was signed. Verify a finite, positive physical trajectory
    # through that point; this short regression is not a developed precursor.
    folder=output/'newton-signed-image'
    m=generate({**common,'surface':'R','model':'cls_original','ks':.008,
        'steps':35,'checkpoint':35,'dt':.05,'advection':'Implicit'},folder)
    tree=ET.parse(folder/'case.xml');root=tree.getroot()
    info(root,'GlobalSysSoln','DirectMultiLevelStaticCond')
    parameter(root,'RANSStageAudit',1);parameter(root,'RANSGradDiv',.01)
    tree.write(folder/'case.xml',encoding='utf-8',xml_declaration=True)
    m['case_sha256']=sha(folder/'case.xml');write_json(folder/'manifest.json',m)
    prepare(folder/'case.xml',build);receipt=run(folder/'case.xml',build,1,first,900)
    stages=receipt['stages']
    require(stages['maximum_nonpositive_linear_images']>0,'Signed-image regression did not exercise its branch')
    require(min(stages['minimum_k'],stages['minimum_epsilon'])>0 and stages['maximum_guard_activations']==0,
        'A signed linear image reached a physical scalar state')
    require(max(stages[k] for k in ('maximum_weak_divergence','maximum_projected_normal',
        'maximum_normal_integral_defect','maximum_resolved_mass_defect'))<1e-7,'Signed-image continuity/mass failed')
    residual=records(folder/'case.rans-residual.csv')
    free=max(float(r['free_coefficient_l2']) for r in residual if r['phase']=='1')
    require(free<1e-7,'Signed-image full weak stage residual failed')
    results['newton_signed_image']=dict(run=receipt,maximum_free_residual=free,
        strip_continuity=strip_continuity(folder,'case'),strip_momentum=flat_momentum_budget(folder/'case.xml'),
        limitation='Short accepted physical evolution through a formerly rejected linear image, not a developed calibration')
    receipt=dict(status='PASS',scope='native inlet/startup/continuity/history verification; not developed calibration',
        source_parent=__import__('subprocess').check_output(['git','-C',ROOT,'rev-parse','HEAD'],text=True).strip(),
        binary_sha256={name:sha(build/name) for name in ('IncNavierStokesSolver','RANSWallFieldAudit','RANSObservationAudit')},
        cpus=args.cpus,maximum_ranks=4,threads=THREADS,results=results)
    write_json(output/'evidence.json',receipt)
    print('native CV3 inlet/startup/history checks PASS',flush=True)


if __name__=='__main__':main()
