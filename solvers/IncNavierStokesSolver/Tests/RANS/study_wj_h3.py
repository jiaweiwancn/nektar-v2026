#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Fresh bounded H3 experiments and passive H2 diagnosis, never H2 continuation."""
import argparse
import hashlib
import os
from pathlib import Path
import time
import xml.etree.ElementTree as ET
import numpy as np
from numpy.testing import assert_allclose
from run_wj_h3 import execute,load,save,sha
from wj_h2_cases import case_xml,continuous_outer
from wj_h2_profiles import Completion,D0,TOP,constitutive
from wj_h3_analysis import actual,identity,point,principal,quadrature,state_context,pressure_complement,algebra,values
from measure_wj_h2 import checkpoint_for


def locations(r,state):
    H=r['parameters'][1];J=(TOP-D0-H)/(TOP-H);out=[]
    # Include all owners of the failing normal face and the first outer cells.
    heights=[.5875,.6720111978440669,.7565223956881337,.8,1.,2.,6.25,12.5,20.]
    xs=sorted(set([0.,4.,6.,7.5,7.9,7.99,7.999999,8.]+[v[a][0] for v in r['fields'][3]['mesh'] for a in ('lower','upper')]))
    for z in heights:
        rho=H+(z-D0-H)/J
        for el in r['fields'][3]['mesh']:
            if not el['lower'][2]-1e-12<=rho<=el['upper'][2]+1e-12:continue
            for x in xs:
                if not el['lower'][0]-1e-12<=x<=el['upper'][0]+1e-12:continue
                v=point(r,state,el['geometry'],x,rho);v['principal']=principal(v['k'],v['epsilon'],np.array(v['G']))
                out.append(v)
    return out


def retained(root):
    h=actual('/media/jay/新加卷/nektar-v2026-runs/wj-h2-20261111');folder=root/'retained-v1';folder.mkdir(exist_ok=False)
    progress=load(h/'study-v1/progress.json');half=load(h/'diagnosis-v1/result.json')
    selections=[('base-cp0',progress['runs']['base-cp0'],0,False),('refined-cp0',progress['runs']['refined-cp0'],0,False),
        ('refined04',progress['runs']['refined-pilot'],2,False),('rejected06',progress['runs']['refined-pilot'],None,True),
        ('half04',half,4,False),('half06',half,6,False),('half08',half,8,False)]
    result=dict(schema='WJH3RetainedDiagnosisV1',original_h2_status='OBSTRUCTED',h2_long_continuation='DENIED_NOT_RUN',states={})
    for key,record,step,rejected in selections:
        rp=actual(record['report']);assert sha(rp)==record['report_sha256'];r=load(rp)
        state,history,alpha,t,packet=state_context(r,step,rejected)
        binding=dict(producer=identity(rp),role='rejected candidate' if rejected else 'accepted checkpoint',time=t,
            state_sha256=hashlib.sha256(np.asarray(state,dtype='<f8').tobytes()).hexdigest(),packet=packet)
        assert binding['state_sha256']==packet['state_sha256']
        if not rejected:
            cp,payload=checkpoint_for(dict(xml=str(actual(record['xml'])),report=str(rp)),step)
            binding['checkpoint']=identity(Path(cp)/'State.bin')
            assert np.asarray(payload['state'],dtype='<f8').tobytes()==np.asarray(state,dtype='<f8').tobytes()
        rows=locations(r,state);q=quadrature(r,state,history,alpha)
        if key=='half08':
            # Previously hash-bound native samples are a separate validation of
            # reconstruction and the SAME-node algebra, not a high-order oracle.
            idx=load(h/'observations-v3/index.json')['refined-fine08-nativeweak']
            fp=actual(idx['native']);assert sha(fp)==idx['native_sha256'];old=load(fp)['weak_checks']
            for sample in [v for v in old['samples'] if not v['inner']][::97]:
                x,_,rho=sample['chart'];v=point(r,state,sample['geometry'],x,rho)
                assert_allclose(v['G'],sample['G'],rtol=2e-10,atol=2e-12)
                assert_allclose([v['k'],v['epsilon'],v['pi']],[sample['k'],sample['epsilon'],sample['pi']],rtol=2e-12,atol=2e-13)
            assert_allclose(q['views'][8]['components']['outer_momentum_stress'],old['native_terms']['outer_momentum_stress'],rtol=2e-11,atol=2e-13)
            binding['native_sample_validation']=identity(fp)
        record=dict(binding=binding,rows=rows,quadrature=q,pressure_test_projection=pressure_complement(r,state),
            interface=r['open_boundary']['interface'] if not rejected and step==r['steps'] else None,
            note='terminal interface rows only belong to terminal accepted state; earlier/rejected interfaces reconstructed separately')
        save(folder/(key+'.json'),record)
        result['states'][key]=dict(path=str(folder/(key+'.json')),sha256=sha(folder/(key+'.json')),time=t,
            minimum_sampled_principal=min(x['principal']['in_plane'] for x in rows),quadrature=q['comparisons'],
            q_star=q['smallest_adequate'],reference_status=q['reference_status'])
    # Numerical independent algebra must also agree with compiled kernel at the
    # actual native failing state, not merely at the prescribed profile.
    from validate_wj1 import Probe
    kernel=actual('/media/jay/新加卷/nektar-v2026-runs/wj-h2-20261111/build-v1/WJKernelProbe');probe=Probe(str(kernel))
    event=load(actual(progress['runs']['refined-pilot']['report']))['physical_range_events'][-1]
    w=event['wj_packet']['principal_witnesses']['in_plane'];G=np.array(w['raw_G_row_major']).reshape(3,3);k=w['k'];e=w['epsilon']
    c=Completion();z=w['position'][2];Gp,Rp=continuous_outer(c,z);kp,ep=c.k(z),c.epsilon(z)
    comparisons={}
    for name,kk,ee,gg in [('failed',k,e,G),('failed_G_profile_scalars',kp,ep,G),('profile_G_failed_scalars',k,e,Gp),
        ('continuous_profile',kp,ep,Gp),('trace_removed_diagnostic',k,e,G-np.trace(G)*np.eye(3)/3)]:
        p=principal(kk,ee,gg);hdir=np.outer(p['polarization'],p['xi']);v=probe.closure(kk,ee,gg)
        d=probe.closure(kk,ee,gg,(0,0,hdir));actual_form=-np.sum(hdir*d['R'])/v['scalar'][1]
        assert_allclose(actual_form,p['in_plane'],rtol=1e-9,atol=2e-11)
        assert_allclose(algebra(np.array([kk]),np.array([ee]),np.array([gg]))[0][0],v['R'],rtol=3e-12,atol=3e-14)
        comparisons[name]=dict(principal=p,compiled_form=float(actual_form),G=gg.tolist(),tau=float(kk/ee))
    probe.close();result['independent_compiled_replay']=comparisons;result['kernel']=identity(kernel)
    save(folder/'summary.json',result);print('Retained passive analysis complete',flush=True)


def native(root,key,steps=0,dt=.02,refined=False,seed='SignCorrectedV1',outlet='ZeroTractionV1',nx=1,q=8,restart=None):
    folder=root/'cases'/key
    xml=case_xml(folder,'physical',steps,dt,refined,'ExactStepV1',restart=restart,pressure_seed=seed,outlet=outlet,
        xknots=np.linspace(0,8,nx+1).tolist(),tangential_order=q)
    runtime=actual('/media/jay/新加卷/nektar-v2026-runs/wj-h2-20261111/runtime-v4/IncNavierStokesSolver')
    intent=dict(xml=identity(xml),runtime=identity(runtime),seed=seed,outlet=outlet,nx=nx,q=q,dt=dt,steps=steps,
        restart=identity(Path(restart)/'State.bin') if restart else None,source_revision=os.environ.get('WJ_H3_SOURCE','recorded git/source hashes in receipt'))
    save(folder/'intent.json',intent)
    okay=execute(root,key,[runtime,xml],folder,steps,300)
    rp=xml.with_suffix('.evolution.json')
    if not rp.is_file():raise RuntimeError('Missing native report; stop, inspect resource/implementation obstruction')
    r=load(rp);attempts=r['pimple_attempts']
    for a in attempts:
        assert [a[k] for k in ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls')]==[1]*5
        assert a['bound_activations']==0,'material correction is a hard stop'
    record=dict(intent=intent,report=identity(rp),status='PASS' if okay and not r['failure'] else 'OBSTRUCTED',
        failure=r['failure'],accepted_step=r['steps'],time=r['time'],attempts=len(attempts),state_size=len(r['state']),
        corrections=sum(a['bound_activations'] for a in attempts),solve_counts=[sum(a[k] for a in attempts) for k in
            ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls')],
        interface=r['open_boundary']['interface'],initialization_work=r['initialization_lift_calls_iterations'])
    save(folder/'result.json',record);print(key,record['status'],record['time'],record['failure'],flush=True)
    return r,record


def setup(root):
    # These are fresh corrected initializations, not imported H2 histories.
    a,ar=native(root,'A-cp0',refined=True)
    old=actual('/media/jay/新加卷/nektar-v2026-runs/wj-h2-20261111/study-v1/refined-cp0/case.xml')
    get=lambda p:{(f.get('NAME'),e.get('VAR')):e.get('VALUE') for f in ET.parse(p).getroot().findall('./CONDITIONS/FUNCTION') for e in f}
    before=get(old);after=get(root/'cases/A-cp0/case.xml')
    for name,value in before.items():
        if name[1]!='p':assert value==after[name],('non-pressure change in sign-only control',name)
    b,br=native(root,'balanced-q8-cp0',dt=.01,seed='VerticalBalanceV1')
    save(root/'setup-result.json',dict(A=ar,balanced=br,sign_only_all_nonpressure_Function_bytes_equal=True,
        note='Balanced q8 CP0 is a passive integration field; B evolution gets fresh q_star CP0/signature'))


def select_q(root):
    retained=load(root/'retained-v1/summary.json');paths=[v['path'] for v in retained['states'].values()]
    qstar=max(load(p)['quadrature']['smallest_adequate'] or 999 for p in paths)
    results={}
    for key in ('A-cp0','balanced-q8-cp0'):
        r=load(root/'cases'/key/'case.evolution.json');state,h,a,t,p=state_context(r)
        q=quadrature(r,state,h,a);save(root/'cases'/key/'passive-quadrature.json',q)
        results[key]=dict(q_star=q['smallest_adequate'],reference_status=q['reference_status'],comparisons=q['comparisons'])
        qstar=max(qstar,q['smallest_adequate'] or 999)
    refs=all(load(p)['quadrature']['reference_status']=='PASS' for p in paths) and all(v['reference_status']=='PASS' for v in results.values())
    result=dict(status='PASS' if qstar<=16 and refs else 'UNRESOLVED',q_star=qstar if qstar<=16 else None,
        corrected_initial_fields=results,consumers=['outer volume x and normal quadrature','inner tangential ray sampling','pressure/source/diagnostic rows using those rays'],
        unchanged_consumers=['inner normal quadrature16','exterior quadrature max(16,q_star)'],
        role='passive integration view; no historical checkpoint numerical signature reinterpreted')
    save(root/'quadrature-selection.json',result);print(result['status'],'q_star',result['q_star'],flush=True)


def controls(root):
    a=load(root/'cases/A-cp0/case.evolution.json')
    r,record=native(root,'A-dt02',4,refined=True,restart=a['checkpoint'])
    # The old refined .02 run failed at .06 whereas .01 reached .08. A changed
    # outcome or repeated domain exit needs precisely the permitted dt companion.
    f,fr=native(root,'A-dt01-cp0',dt=.01,refined=True)
    rr,rec=native(root,'A-dt01',8,dt=.01,refined=True,restart=f['checkpoint'])
    q=load(root/'quadrature-selection.json')
    if q['status']!='PASS':save(root/'B-NOT_RUN.json',dict(status='NOT_RUN',reason='unqualified outer load quadrature'));return
    for name,nx,outlet in (('B-zero-x1',1,'ZeroTractionV1'),('B-frozen-x1',1,'FrozenSeedTractionV1'),
            ('B-zero-x3',3,'ZeroTractionV1'),('B-frozen-x3',3,'FrozenSeedTractionV1')):
        cp,cr=native(root,name+'-cp0',dt=.01,seed='VerticalBalanceV1',outlet=outlet,nx=nx,q=q['q_star'])
        if cr['status']!='PASS':continue
        native(root,name,8,dt=.01,seed='VerticalBalanceV1',outlet=outlet,nx=nx,q=q['q_star'],restart=cp['checkpoint'])


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('retained','setup','select-q','controls'),required=True)
    a=p.parse_args();assert a.root.name=='wj-h3-20261112';globals()[a.mode.replace('-','_')](a.root)
if __name__=='__main__':main()
