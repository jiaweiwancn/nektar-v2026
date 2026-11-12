#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H3 versioned passive observations, actual callback contexts and comparisons."""
import argparse
import hashlib
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from numpy.testing import assert_allclose
from run_wj_h3 import execute,load,save,sha
from wj_h3_analysis import actual,identity,point,principal,quadrature,state_context,pressure_complement,interface_view
from study_wj_h3 import locations,native
from wj_h2_profiles import Completion,D0,TOP
from wj_h2_cases import continuous_outer,reduced_normal_offset
from measure_wj_h2 import checkpoint_for
from wj_h2_weak import independent_weak


def view(report,step=None,rejected=False):
    state,history,alpha,t,packet=state_context(report,step,rejected)
    digest=hashlib.sha256(np.asarray(state,dtype='<f8').tobytes()).hexdigest()
    assert digest==packet['state_sha256']
    rows=locations(report,state)
    # Keep the native minimizing locations even when outside the fixed profile
    # grid. Their witness is a finite native inventory, not a global statement.
    for name,w in packet['principal_witnesses'].items():
        x,_,rho=w['chart'];v=point(report,state,w['geometry'],x,rho)
        assert_allclose(v['G'],np.array(w['raw_G_row_major']).reshape(3,3),rtol=3e-10,atol=2e-12)
        v['principal']=principal(v['k'],v['epsilon'],np.array(v['G']));v['native_witness']=name;rows.append(v)
    xs=sorted({0.,4.,6.,7.5,7.9,7.99,8.}|{float(e[a][0]) for e in report['fields'][0]['mesh'] for a in ('lower','upper')})
    interface=[v for x in xs for v in interface_view(report,state,x)]
    return dict(time=t,role='rejected candidate' if rejected else 'accepted state',state_sha256=digest,
        callback_history_sha256=packet['stage_history_sha256'],alpha_mass=alpha,packet=packet,rows=rows,interface=interface,
        minimum_principal=min(v['principal']['in_plane'] for v in rows),maximum_sampled_divergence=max(abs(v['divergence']) for v in rows),
        pressure_test_projection=pressure_complement(report,state))


def analyze(root):
    out=root/'analysis-v2';out.mkdir(exist_ok=False);summary=dict(schema='WJH3ShortControlAnalysisV2',cases={})
    # Old/new comparisons use native CP0 coefficients, not a prescribed profile.
    h=actual('/media/jay/新加卷/nektar-v2026-runs/wj-h2-20261111')
    records={}
    for key in ('refined-cp0','refined-pilot'):
        rp=h/'study-v1'/key/'case.evolution.json';records['H2-'+key]=rp
    records['H2-half']=h/'diagnosis-v1/half-dt/case.evolution.json'
    for f in sorted((root/'cases').glob('*/result.json')):records[f.parent.name]=f.with_name('case.evolution.json')
    for key,rp in records.items():
        producer=identity(rp);r=load(rp);first=r['steps']+1-len(r['accepted_states']);states=[]
        for step in range(first,r['steps']+1):
            # H3 all common times are retained; H2 old refined CP0 and .02/.04,
            # half-dt .04/.06/.08 suffice for the predetermined causal contrast.
            if key=='H2-half' and step not in (0,4,6,8):continue
            v=view(r,step);v['absolute_step']=step;states.append(v)
        if r['failure'] and r['physical_range_events']:states.append(view(r,rejected=True))
        result=dict(producer=producer,xml=identity(rp.with_name('case.xml')),states=states,status='OBSTRUCTED' if r['failure'] else 'PASS',
            failure=r['failure'],accepted_steps=r['steps'],accepted_time=r['time'],state_size=len(r['state']),
            open_boundary=r['open_boundary'],raw_reference_norms=r['raw_reference_norms'],
            regional_stage_loads=r['regional_stage_loads'],final_stage_constraints=r['final_stage_constraints'])
        if key.startswith('B-') and not key.endswith('-cp0'):
            # Decisive accepted AND rejected endpoints retain a fixed q_star;
            # a later integration failure is reported, never adaptively rerun.
            result['terminal_quadratures']={}
            for rejected in ([False,True] if r['failure'] else [False]):
                s,hist,alpha,t,p=state_context(r,rejected=rejected)
                result['terminal_quadratures']['rejected' if rejected else 'accepted']=quadrature(r,s,hist,alpha)
        save(out/(key+'.json'),result)
        assert sha(rp)==producer['sha256'],'passive analysis changed producing bytes'
        summary['cases'][key]=dict(path=str(out/(key+'.json')),sha256=sha(out/(key+'.json')),status=result['status'],
            accepted_steps=r['steps'],accepted_time=r['time'],state_size=len(r['state']),failure=r['failure'],
            principal_by_time=[dict(time=v['time'],role=v['role'],minimum=v['minimum_principal']) for v in states])
    save(out/'summary.json',summary);print('H3 passive all-state analysis complete',flush=True)


def observe(root):
    out=root/'native-observations-v1';out.mkdir(exist_ok=False);index={}
    observer=actual('/media/jay/新加卷/nektar-v2026-runs/wj-h2-20261111/observer-v6/RANSWallResolutionProbe')
    selected=['A-cp0','balanced-q8-cp0','A-dt02','A-dt01','B-zero-x1-cp0','B-frozen-x1-cp0',
        'B-zero-x1','B-frozen-x1','B-zero-x3','B-frozen-x3']
    for key in selected:
        rp=root/'cases'/key/'case.evolution.json'
        if not rp.exists():continue
        r=load(rp);folder=out/key;folder.mkdir()
        cp,payload=checkpoint_for(dict(xml=str(rp.with_name('case.xml')),report=str(rp)),r['steps'])
        request=dict(operation='open_checks',inner_stress_model='Linear',native_report=str(rp),native_report_sha256=sha(rp),
            native_checkpoint=cp,checkpoint_sha256=sha(Path(cp)/'State.bin'),allow_failed_producer=bool(r['failure']),
            physical_observations=[[x,0,z] for x in (0.,4.,7.5,7.99,8.) for z in (.5875,.6720111978440669,1.,2.,6.25,12.5,20.)],
            weak=key.startswith('B-'),weak_outer_order=32,assembly=False,reuse=False)
        save(folder/'request.json',request)
        assert execute(root,'passive-native-'+key,[observer,rp.with_name('case.xml'),folder/'request.json',folder/'native.json'],folder,0,120)
        p=load(folder/'native.json');assert p['initialization_work']==[0,0]
        v=dict(producer=identity(rp),checkpoint=identity(Path(cp)/'State.bin'),native=identity(folder/'native.json'),
            callback_history_sha256=p['callback_history_sha256'],selected_absolute_step=p['selected_absolute_step'],
            stage_norms=p['physical_mass_dual_norms'])
        if request['weak']:
            assert p['weak_checks']['governing_solves']==0
            v['independent_weak']=independent_weak(p['weak_checks'],dict(r,open_boundary=p['open_boundary']),enforce=False)
        c=Completion();C=reduced_normal_offset(c.K,continuous_outer(c,D0+c.H)[1][2,2]);errors=[]
        for row in p['open_boundary']['boundary_rows']:
            if row['role']!='outlet':continue
            if 'frozen' not in key:expected=np.zeros(3)
            elif row['inner']:expected=np.array([-C,0,c.UH/c.w['J0']])
            else:
                R=continuous_outer(c,row['physical'][2])[1]
                expected=np.array([R[2,2]-R[0,0]-C,0,-R[2,0]])
            errors.append(float(np.max(abs(expected-row['applied_traction']))))
        v['independent_applied_outlet_traction_max_error']=max(errors,default=0.)
        assert v['independent_applied_outlet_traction_max_error']<2e-11
        index[key]=v;save(out/'index.json',index)
    print('Native read-only diagnostics complete',flush=True)


def count_and_signature(root):
    q=load(root/'quadrature-selection.json')['q_star']
    r,v=native(root,'native-size-x2',dt=.01,seed='VerticalBalanceV1',nx=2,q=q)
    assert len(r['state'])==1846
    # A genuine supported seed change must reject the old signature before any
    # governing attempt. It is not a corrected-initial-state restart operation.
    from wj_h2_cases import case_xml
    source=root/'cases/A-cp0/CP';folder=root/'seed-signature-rejection'
    xml=case_xml(folder,steps=0,refined=True,pressure_seed='VerticalBalanceV1',restart=source,reuse='ExactStepV1')
    runtime=actual('/media/jay/新加卷/nektar-v2026-runs/wj-h2-20261111/runtime-v4/IncNavierStokesSolver')
    original=identity(source/'State.bin')
    assert execute(root,'seed-signature-rejection',[runtime,xml],folder,0,120,expected='signature')
    assert sha(source/'State.bin')==original['sha256']
    save(folder/'result.json',dict(status='PASS',original_checkpoint=original,reason='changed supported pressure-seed Function identity'))


def startup_pair(root):
    # One additional zero-step field isolates the two seed choices on identical
    # base mesh/q8 inputs. No trajectory, new timestep or fitted profile is added.
    save(root/'matched-startup-intent.json',dict(role='passive initializer comparison only',native_attempts=0,
        seed='SignCorrectedV1',comparison='balanced-q8-cp0',normal_mesh='base23',q=8,dt=.01,
        preserves_original_U_k_epsilon_roughness=True))
    r,v=native(root,'sign-only-base-cp0',dt=.01)
    folder=root/'matched-startup-observation';folder.mkdir()
    rp=root/'cases/sign-only-base-cp0/case.evolution.json';cp=Path(r['checkpoint'])
    request=dict(operation='open_checks',inner_stress_model='Linear',native_report=str(rp),native_report_sha256=sha(rp),
        native_checkpoint=str(cp),checkpoint_sha256=sha(cp/'State.bin'),physical_observations=[[4.,0,.5875]],
        weak=False,assembly=False,reuse=False)
    save(folder/'request.json',request)
    observer=actual('/media/jay/新加卷/nektar-v2026-runs/wj-h2-20261111/observer-v6/RANSWallResolutionProbe')
    assert execute(root,'matched-startup-observation',[observer,rp.with_name('case.xml'),folder/'request.json',folder/'native.json'],folder,0,120)
    p=load(folder/'native.json');assert p['initialization_work']==[0,0]
    old=load(root/'cases/balanced-q8-cp0/case.evolution.json')
    for block in range(4):
        lo,hi=r['offsets'][block:block+2]
        assert np.asarray(r['state'][int(lo):int(hi)],dtype='<f8').tobytes()==np.asarray(old['state'][int(lo):int(hi)],dtype='<f8').tobytes()
    save(folder/'result.json',dict(status='PASS',producer=identity(rp),native=identity(folder/'native.json'),
        minimal_offset_spatial=p['physical_mass_dual_norms'],
        balanced_spatial=load(root/'native-observations-v1/index.json')['balanced-q8-cp0']['stage_norms'],
        differential_initial_coefficients_identical=True,scope='same mesh/order/physical completion; pressure initialization only'))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('analyze','observe','count-and-signature','startup-pair'),required=True)
    a=p.parse_args();globals()[a.mode.replace('-','_')](a.root)
