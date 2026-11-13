#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Conditional H14 target admission; never generate a scientific mesh/profile."""
import argparse
import math
from pathlib import Path
import subprocess
import sys
import time
import numpy as np
from run_wj_h4 import load,save,sha
from wj_h14_execution import clone,file_id,target_preflight,validate_pair
from study_wj_h14 import H13,verify,stored_reference,assessment

def prepare(root):
    gate=load(root/'science/qualification-gates.json');assert gate['status']=='PASS'
    assert load(root/'science/small/unsupported/result.json')['status']=='PASS'
    assert load(root/'science/runtime-qualified-v4.json')['status']=='PASS'
    manifest=load(H13/'publication-v1/manifest.json');parent=H13/'publication-v1/paired-case'
    selected=[]
    for p in(parent/'flat/case.xml',parent/'paired-domain.json',parent/'preservation-protocol.json',parent/'mapped-backend-prerequisite.json'):
        rec=next(v for v in manifest['artifacts']if v['path']==str(p));verify(rec);selected.append(rec)
    m=load(parent/'paired-domain.json');c=stored_reference();dt=.0025
    # Cell-based mesh scales and profile-variation scales are different. IMEX
    # diffusion is implicit; neither number below certifies temporal adequacy.
    from wj_h13_domain import topology
    t=topology(parent/'flat/case.xml');z=t['chart_normal_knots'];dx=min(np.diff(t['xknots']));J=(c.top-c.d0-c.H)/(c.top-c.H)
    local=[]
    for a,b in zip(z[:-1],z[1:]):
        r=(a+b)/2;inner=b<=c.H;physical_r=r if inner else c.H+(r-c.H)*J;D=c.D(r)if inner else c.a*physical_r
        dz=(b-a)*(1 if inner else J);length=D/c.a if inner else physical_r
        local.append(dict(chart_interval=[a,b],physical_normal_width=dz,region='inner'if inner else'outer',diffusivity=D,
            dt_D_over_cell_width_squared=dt*D/dz**2,profile_diffusion_time=length*length/D,dt_over_profile_diffusion_time=dt*D/length**2))
    clock=dict(dt=dt,maximum_steps=4,dt_source='existing qualified6mm small timestep; selected for a capacity/admissibility pilot only',
        minimum_horizontal_cell=dx,max_reference_U=c.U_outer(c.top-c.d0),
        cell_advective_CFL=dt*c.U_outer(c.top-c.d0)/dx,Q3_degree_squared_advective_indicator=9*dt*c.U_outer(c.top-c.d0)/dx,
        max_cell_diffusion_indicator=max(v['dt_D_over_cell_width_squared']for v in local),
        max_profile_diffusion_indicator=max(v['dt_over_profile_diffusion_time']for v in local),normal_intervals=local,
        limitation='implicit diffusion has no asserted explicit stability bound; high-frequency cell response and temporal preservation remain unqualified; no timestep search')
    book=load(root/'accounting.json');charged=60+sum(v.get('seconds',v['reserved_seconds'])for v in book['entries'].values())
    forecast=dict(status='PASS',charged_before_admission=charged,headroom=1.25,untouched_closure=300,
        remaining_nominal_seconds=dict(target_cp0=600,target_pilot=600,passive_strict_and_final_analysis=700,transport_and_verification=150),
        projected_total_seconds=charged+1.25*(600+600+700+150)+300,
        measured_small_cp0_seconds=book['entries']['small-mapped-CP0-v1']['seconds'],
        measured_small_two_step_seconds=book['entries']['small-mapped-two-step-v1']['seconds'],
        target_rate='unknown until actual CP0/pilot; finite phase ceilings retained, no small-to-target cost claim')
    if forecast['projected_total_seconds']>7200:forecast['status']='RESOURCE_LIMIT'
    design=dict(schema='WJH14Fixed6mmResourceDesignV1',height_mm=6,target_state=24688,selected_H13=selected,
        source_allocation=file_id(root/'science/target-source-allocation-v1.json'),qualified_gates=file_id(root/'science/qualification-gates.json'),
        unsupported=file_id(root/'science/small/unsupported/result.json'),runtime=file_id(root/'science/runtime-qualified-v4.json'),
        limits=dict(state=24688,pressure=1024,factor_payload_bytes=512*1024**2,native_RSS_bytes=8*1024**3,science_bytes=2*1024**3,pilot_seconds=600),
        forecast=forecast,pilot_rule='one-step genuine checkpoint continuations, at most4, after each independent hard/quality/moment assessment; before each next step require25% cost headroom within the same600s reservation; no retry',
        timestep=clock,scientific_data_unchanged=True)
    path=root/'science/target-resource-design-v1.json';save(path,design)
    original=parent/'flat/case.xml';xml=clone(original,root/'science/target/cp0',dt=dt,design=sha(path))
    command=[str(root/'runtime-v4/IncNavierStokesSolver'),str(xml)]
    script=Path(__file__).resolve();pilot=[sys.executable,'-B',str(script),'--root',str(root),'--mode','pilot']
    admission=target_preflight(original,xml,m,gate['gates'],{'status':'PASS'},forecast)
    admission.update(commands=dict(target_cp0=command,pilot=pilot),command_inputs=dict(target_cp0=[file_id(xml),file_id(root/'runtime-v4/IncNavierStokesSolver'),file_id(root/'runtime-v4/IncNavierStokesSolver.elf')],
        pilot=[file_id(script),file_id(Path(assessment.__code__.co_filename)),file_id(xml)]),resource_design=file_id(path),maximum_pilot_attempts=4)
    save(root/'target-admission.json',admission);print('target admission',admission['status'],'projected charged seconds',forecast['projected_total_seconds'],flush=True)

def pilot(root):
    start=time.monotonic();admission=load(root/'target-admission.json');assert admission['status']=='PASS'
    cp0=root/'science/target/cp0';a=load(cp0/'assessment/result.json');assert a['status']=='PASS'
    original=H13/'publication-v1/paired-case/flat/case.xml';design=admission['resource_design']['sha256'];dt=.0025
    prior=load(cp0/'case.evolution.json');assert prior['steps']==0 and prior['dt']==dt
    records=[];stop=None;cp0cost=load(root/'accounting.json')['entries']['target-full-domain-CP0-v1']['seconds']
    estimate=max(30.,2*cp0cost);assessment_estimate=load(root/'accounting.json')['entries']['target-CP0-owned-assessment-v1']['seconds']
    for step in range(1,5):
        elapsed=time.monotonic()-start;required=1.25*(estimate+assessment_estimate)
        projection=dict(next_absolute_step=step,elapsed=elapsed,remaining=600-elapsed,required_with25percent_headroom=required,
            native_estimate_seconds=estimate,assessment_estimate_seconds=assessment_estimate,prior_records=len(records))
        save(root/f'science/target/pilot-projection-{step}.json',projection)
        if elapsed+required+5>600:stop='RESOURCE_LIMIT';break
        folder=root/f'science/target/pilot{step}';xml=clone(original,folder,dt=dt,steps=1,restart=prior['checkpoint'],design=design)
        launch=time.monotonic()
        with(folder/'native.log').open('x')as out:
            code=subprocess.run([str(root/'runtime-v4/IncNavierStokesSolver'),str(xml)],stdout=out,stderr=subprocess.STDOUT).returncode
        native_seconds=time.monotonic()-launch
        if not(folder/'case.evolution.json').exists():
            records.append(dict(step=step,native_exit=code,status='UNRESOLVED',native_seconds=native_seconds));stop='UNRESOLVED';break
        r=load(folder/'case.evolution.json');rec=dict(step=step,native_exit=code,native_seconds=native_seconds,report=file_id(folder/'case.evolution.json'),native_failure=r['failure'],attempts=len(r['pimple_attempts']))
        records.append(rec)
        if code or r['failure']:rec['status']='OBSTRUCTED';stop='OBSTRUCTED';break
        assert r['steps']==step and r['dt']==dt and r['initial_velocity_calls']==0 and r['initialization_lift_calls_iterations']==[0,0]
        from wj_h9_local import same,phase_bindings
        assert same(r['accepted_states'][0],prior['state'])
        for k in('fixed_reference','fixed_defect'):assert same(r[k],prior[k])
        p=r['pimple_attempts'][0];assert all(p[k]==1 for k in('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls'))
        assert p['bound_activations']==p['additional_correction_or_restoration_solves']==0
        snap=load(p['k_snapshot']);assert sha(p['k_snapshot'])==p['k_snapshot_sha256'];_,binding=phase_bindings(r,snap,step);assert binding['P_C_bitwise'];rec['phase_binding']=binding
        scan=time.monotonic();assessment(root,f'target/pilot{step}');assessment_seconds=time.monotonic()-scan
        result=load(folder/'assessment/result.json');rec.update(status=result['status'],assessment_seconds=assessment_seconds,assessment=file_id(folder/'assessment/result.json'))
        if rec['status']!='PASS':stop=rec['status'];break
        prior=r;estimate=max(native_seconds,estimate*.5);assessment_estimate=max(assessment_seconds,assessment_estimate*.5)
    complete=sum(v.get('status')=='PASS'for v in records)
    save(root/'science/target/pilot-result.json',dict(status=stop if stop and stop!='RESOURCE_LIMIT'else'PASS'if complete else'RESOURCE_LIMIT',
        completed_steps=complete,actual_attempts=sum(v.get('attempts',1)for v in records),reserved_maximum_attempts=4,records=records,
        remaining_steps_status='RESOURCE_LIMIT'if stop=='RESOURCE_LIMIT'else'NOT_RUN'if complete<4 else'COMPLETE',stop_reason=stop,
        elapsed_seconds=time.monotonic()-start,scope='capacity/admissibility pilot only; no duration preservation or physical experimental qualification'))
    print('full-domain pilot',complete,'completed steps; stop',stop,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('prepare','pilot'),required=True);a=p.parse_args();globals()[a.mode](a.root)
