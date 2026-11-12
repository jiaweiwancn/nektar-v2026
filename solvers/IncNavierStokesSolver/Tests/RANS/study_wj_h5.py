#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H5 continuum calibration and immutable CP0 observations; never time advance."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import numpy as np
from numpy.testing import assert_allclose
from run_wj_h5 import load,save,sha
from wj_h5_neutral import *
from wj_h5_cases import case_xml
from wj_h5_oracle import check_reference
from wj_h4_points import owners,evaluate,vector_hash
from wj_h3_analysis import outer_loads
from wj_h2_weak import independent_weak
from validate_wj1 import Probe


def identity(p):
    p=Path(p);return dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p))


def configurations():
    fit=mean_fit()
    return [NeutralReference(H/SCALE_MM,fit['z0_e'])for H in (5.,5.5,6.)]


def verify(root):
    source=Path(__file__).resolve().parents[4]
    # These are the already-qualified compilation and loaded-ELF inventories,
    # not a claim that the current Python revision built the native executables.
    sealed=source/'docs/plans/wj-h2-linux';inputs=load(sealed/'native-source-inputs.json')
    checks=[]
    for item in inputs:
        path=source/item['path'];actual=sha(path)
        assert actual==item['sha256'],f'changed compiled source input: {path}'
        checks.append(dict(path=item['path'],sha256=actual))
    old=load(sealed/'loaded-runtime.json');elf=[];other=[]
    for v in old['loaded_files']:
        path=Path(v['path'].replace('/media/jay/新加卷/','/media/jay/新加卷1/',1))
        if not path.is_file():other.append(dict(path=str(path),status='not a required ELF input'));continue
        with path.open('rb')as f:magic=f.read(4)
        if magic!=b'\x7fELF':other.append(dict(path=str(path),recorded_sha256=v['sha256'],actual_sha256=sha(path)));continue
        actual=identity(path);assert actual['sha256']==v['sha256'],f'changed qualified ELF: {path}'
        elf.append(actual)
    runs=root.parent;h2=runs/'wj-h2-20261111'
    executables={
        'governing':(h2/'runtime-v4/IncNavierStokesSolver','9ae37a1aba62755b10ad99c9484c498cef9d927a75fb950f00344f4ffc5d06fe'),
        'observer':(h2/'observer-v6/RANSWallResolutionProbe','6f35b837f0dd13eb155e0fe34f0dc670c609ee6cbe739609a5de288dd49ed06d'),
        'kernel':(h2/'build-v1/WJKernelProbe','3d4edf79a2718a281175be1a575776a33dd5f988c48d5a1b594a9b5f2dd7844c')}
    exe={}
    for role,(path,expected)in executables.items():
        exe[role]=identity(path);assert exe[role]['sha256']==expected
    specs=[]
    for name,expected in (('wj_earsm_neutral_abl_ep_awf_note.tex','1959e96ea0b61f9deb9ff9ad1f5b26488a932a8451c2de6d3b6f7b7c9dbd9124'),
        ('wj_earsm_budget_cutoff_awf_appendix.tex','0ab52e0d813b8fb344973f75ca6d8da97ff389ef9f93476500ca9d1df1fbc47c')):
        v=identity(source/'docs/literature'/name);assert v['sha256']==expected;specs.append(v)
    save(root/'runtime-equivalence.json',dict(status='PASS',compiled_inputs=checks,qualified_ELFs=elf,informational_non_ELF=other,
        executables=exe,specifications=specs,generating_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=source,text=True).strip(),
        governing_producer='c344dc2209f030e426e2ee5501d5e5b72a6d2c46',observer_producer='e62b5f4d265f068dd5af3f090bfb0e385c3ae53e',
        reused_arithmetic_policy='GNU11.4 C++17 Release -O3 -DNDEBUG -ffp-contract=off -fno-fast-math; RANS_BUILD_SOLVER=ON',
        source_inventory=identity(sealed/'native-source-inputs.json'),runtime_inventory=identity(sealed/'loaded-runtime.json'),
        loader=os.environ.get('LD_LIBRARY_PATH'),native_rebuild=False))
    print(f'PASS: {len(checks)} compiled source inputs, {len(elf)} qualified ELF identities; separate Python/native producers',flush=True)


def write_csv(path,rows):
    with Path(path).open('w',newline='')as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def reference(root,kernel):
    out=root/'reference';out.mkdir(exist_ok=False);fit=mean_fit();cfg=configurations()
    probe=Probe(str(kernel));oracles={}
    try:
        for c in cfg:
            assert all(c.summary()['positivity_extension'].values())
            oracles[str(c.H*SCALE_MM)]=check_reference(c,probe)
    finally:probe.close()
    fits=dict(primary=fit,omit_lowest=mean_fit(omit_lowest=True),through_50mm=mean_fit(window_mm=50))
    rows=[];obs=[]
    for c in cfg:
        s=c.summary();s['printed_z0_required_h_mm']=SCALE_MM*required_roughness(c.H,.2/SCALE_MM)['h']
        rows.append(s);obs+=observation_regions(c)
    save(out/'calibration.json',dict(fits=fits,feasibility=feasibility(fit),candidates=rows,oracles=oracles,
        source_data=identity(DATA),digitization=identity(DIGITIZATION),governing_timesteps=0,
        old_family='PAUSED; H2-H4 evidence and denials unchanged'))
    keys=('H_w_mm','physical_z_H_mm','h_s_mm','printed_z0_required_h_mm','rc_mm','Lambda','d_over_nu',
        'k_ground_over_K','reference_covariance_min_over_utau2','H_plus','outer_molecular_turbulent_ratio',
        'inner_molecular_fraction','h_over_H','physical_plane_over_cube_height','gradient_ratio')
    write_csv(out/'candidates.csv',[{k:s[k]for k in keys}for s in rows])
    write_csv(out/'observations.csv',obs)
    write_csv(out/'fit-residuals.csv',[dict(fit=name,**r)for name,v in fits.items()for r in v['points']])
    print('PASS: fixed mean fit, feasibility, three independent continuum/compiled kernel references',flush=True)


def generate(root):
    intent=[]
    for c in configurations():
        xml=case_xml(root/'cases'/f'H{SCALE_MM*c.H:g}mm',c)
        intent.append(dict(height_mm=SCALE_MM*c.H,xml=identity(xml),design=identity(xml.with_name('case-design.json')),
            governing_timesteps=0,reference='Functions',quadrature=[12,16],boundary='compatible complete regional outlet + reference top shear/epsilon flux'))
    save(root/'initialization-intent.json',intent)
    print('Generated exactly three zero-step inputs; no continuation option',flush=True)


def compare_loads(a,b):
    errors={n:np.asarray(a['components'][n])-b['components'][n]for n in a['components']}
    conditions={n:bool(np.all(np.isfinite(v))and np.all(abs(v)<=2e-9+2e-8*abs(np.asarray(b['components'][n]))))for n,v in errors.items()}
    regional=np.asarray(a['regional'])-b['regional']
    conditions['regional_sum']=bool(np.all(np.isfinite(regional))and np.all(abs(regional)<=4e-9+2e-8*abs(np.asarray(b['regional']))))
    return dict(status='PASS'if all(conditions.values())else 'UNRESOLVED',conditions=conditions,
        component_max_errors={n:float(np.max(abs(v)))for n,v in errors.items()},regional_max_error=float(np.max(abs(regional))),
        reference='same physical represented fields, independent outer element quadrature; fixed inner/exterior loads')


def assess(root,height,observer):
    c=next(c for c in configurations()if abs(c.H*SCALE_MM-height)<1e-12)
    folder=root/'cases'/f'H{height:g}mm';rp=folder/'case.evolution.json';r=load(rp)
    assert r['steps']==0 and r['time']==0 and not r['pimple_attempts'] and r['implicit_callback_count']==0
    assert r['last_stage_alpha']==0
    assert len(r['state'])<=3000 and r['parameters'][:3]==[c.nu,c.H,c.h]
    failed=bool(r['failure']);cp=Path(r['checkpoint'])if r['checkpoint']else None
    original={str(p):sha(p)for p in (rp,folder/'case.xml')}
    if cp:original[str(cp/'State.bin')]=sha(cp/'State.bin')
    if failed:
        event=r['physical_range_events'][-1];packet=event['wj_packet']
        assert not r['accepted_states'] and not r['wj_audit_packets'] and cp is None
        assert event['alpha']==0 and event['time']==0 and vector_hash(event['state'])==vector_hash(r['state'])
    else:
        assert len(r['wj_audit_packets'])==1;packet=r['wj_audit_packets'][0]
    out=root/'observations'/f'H{height:g}mm';out.mkdir(parents=True,exist_ok=False)
    zs=sorted(set([c.d0,c.d0+c.H,c.d0+c.w['rc'],c.d0+c.h,7/SCALE_MM,250/SCALE_MM,c.top]+[
        row['z_mm']/SCALE_MM for row in measurements()]))
    pts=[[x,0,z]for x in (0.,4.,8.)for z in zs]
    for w in packet['principal_witnesses'].values():
        p=list(w['position']);pts.append(p)
        for axis in (0,2):
            for sign in (-1,1):
                q=p.copy();q[axis]+=sign*1e-6
                if 0<=q[0]<=8 and c.d0<=q[2]<=c.top:pts.append(q)
    pts=[list(v)for v in sorted(set(tuple(v)for v in pts))]
    p=None
    if not failed:
        request=dict(operation='open_checks',inner_stress_model='Linear',native_report=str(rp),native_report_sha256=sha(rp),
            native_checkpoint=str(cp),checkpoint_sha256=sha(cp/'State.bin'),physical_observations=pts,
            weak=True,weak_outer_order=32,assembly=False,reuse=False)
        save(out/'request.json',request)
        subprocess.run([str(observer),str(folder/'case.xml'),str(out/'request.json'),str(out/'native.json')],check=True)
        p=load(out/'native.json')
        assert p['initialization_work']==[0,0] and p['governing_solves']==0 and p['selected_absolute_step']==0
        assert not p['actual_completed_stage_residual'] and not p['physical_mass_dual_norms']['completed_stage']
        assert p['state_sha256']==vector_hash(r['state'])
    assert all(sha(k)==v for k,v in original.items()),'read-only observer changed producing bytes'
    weak=independent_weak(p['weak_checks'],dict(r,open_boundary=p['open_boundary']),enforce=False)if p else dict(
        status='NOT_RUN',reason='startup rejected before native checkpoint/accepted lineage; strict immutable reader cannot admit it')
    # q12 is fixed, with no adaptive governing order or historical rerun.
    views={q:outer_loads(r,r['state'],r['state'],0.,q)for q in (12,24,32)}
    quadrature={f'{q}-32':compare_loads(views[q],views[32])for q in (12,24)}
    save(out/'outer-loads.json',views)
    rows=[];interface=[]
    for physical_point in pts:
        for o in owners(r,physical_point):
            v=evaluate(r,r['state'],physical_point,o);dist=physical_point[2]-c.d0
            inside=o['region']=='inner'
            expected=dict(U=[c.U_inner(dist)if inside else c.U_outer(dist),0,0],
                k=c.k_inner(dist)if inside else c.K,epsilon=None if inside else c.epsilon(dist),
                pi=-c.K/NE if inside else 0.,R=(c.covariance_inner(dist)if inside else c.Rc).tolist())
            v['continuum']=expected
            v['errors']={k:float(np.max(abs(np.asarray(v[k])-expected[k])))for k in ('U','k','pi','R')}
            if not inside:v['errors']['epsilon']=abs(v['epsilon']-expected['epsilon'])
            rows.append(v)
        if abs(physical_point[2]-(c.d0+c.H))<1e-14:
            both=rows[-2:];assert {v['owner']['region']for v in both}=={'inner','outer'}
            i=next(v for v in both if v['owner']['region']=='inner');o=next(v for v in both if v['owner']['region']=='outer')
            interface.append(dict(physical=physical_point,shared_U_max=float(np.max(abs(np.asarray(i['U'])-o['U']))),
                shared_k=i['k']-o['k'],endpoint_K=o['sampled_K'],cutoff=i['cutoff'],
                traction_jump=(np.asarray(i['T'])[:,2]-np.asarray(o['T'])[:,2]).tolist(),
                k_flux_jump=-i['nu_d']*i['grad_k'][2]+o['nu_d']*o['grad_k'][2],
                epsilon_trace=o['epsilon'],epsilon_target=(o['sampled_K']/QE)**1.5/(KAPPA*c.H),
                pressure_jump=i['pi']-o['pi'],reference_pressure_jump=-c.K/NE,
                inner_normal_gradient=i['G'][0][2],outer_normal_gradient=o['G'][0][2]))
    save(out/'physical-rows.json',rows)
    outer=[v for v in rows if v['owner']['region']=='outer'];inner=[v for v in rows if v['owner']['region']=='inner']
    minimum=lambda field,region:min(float(v[field])for v in region)
    principal_min=min(min(v['principal']['in_plane'],v['principal']['transverse'])for v in outer)
    cov_min=min(min(v['covariance'])for v in rows)
    kernel=Probe(str(root.parent/'wj-h2-20261111/build-v1/WJKernelProbe'));witnesses={}
    try:
        for name,w in packet['principal_witnesses'].items():
            v=next(v for v in rows if v['owner']['geometry']==w['geometry']and
                np.max(abs(np.asarray(v['physical'])-w['position']))<1e-14)
            G=np.asarray(v['G']);H=np.asarray(w['H_row_major']).reshape(3,3)
            assert_allclose(G.ravel(),w['raw_G_row_major'],rtol=2e-11,atol=2e-12)
            value=kernel.closure(v['k'],v['epsilon'],G);action=kernel.closure(v['k'],v['epsilon'],G,(0.,0.,H))
            replay=-np.sum(H*action['R'])/value['scalar'][1]
            assert_allclose(replay,w['minimum_over_nu_eff'],rtol=2e-9,atol=2e-11)
            rphys=w['position'][2]-c.d0;Gc=np.zeros((3,3));Gc[0,2]=c.u_tau/(KAPPA*rphys)
            from wj_h3_analysis import principal
            witnesses[name]=dict(native=w,compiled_replay=float(replay),independent_principal=v['principal'],
                independent_coefficient_row=v,continuum_G=Gc.tolist(),continuum_epsilon=c.epsilon(rphys),
                continuum_principal=principal(c.K,c.epsilon(rphys),Gc),
                note='same physical point and declared owner; no replacement of rejected coefficients or admitted checkpoint')
    finally:kernel.close()
    represent=dict(min_k=minimum('k',rows),min_epsilon=minimum('epsilon',outer),min_covariance=cov_min,
        min_planar_principal_over_nu_eff=principal_min,max_root_residual=max(abs(v['root_residual'])for v in outer),
        max_raw_divergence=max(abs(v['raw_divergence'])for v in rows),
        max_errors={key:max(v['errors'].get(key,0)for v in rows)for key in ('U','k','epsilon','pi','R')},
        scope='declared finite physical sample set, all legitimate owners; startup projection, no evolved-flow claim')
    gates=dict(native_initialization=not r['failure'],positive_scalars=represent['min_k']>0 and represent['min_epsilon']>0,
        covariance=cov_min>0,principal=principal_min>0,root=represent['max_root_residual']<1e-12,
        independent_weak=weak['status']=='PASS'if p else None,quadrature=all(v['status']=='PASS'for v in quadrature.values()),
        exact_shared_traces=all(v['shared_U_max']<2e-13 and abs(v['shared_k'])<2e-13 and abs(v['epsilon_trace']-v['epsilon_target'])<2e-13 for v in interface))
    result=dict(status='OBSTRUCTED'if failed else 'PASS'if all(gates.values())else 'UNRESOLVED',gates=gates,producer=identity(rp),xml=identity(folder/'case.xml'),
        checkpoint=identity(cp/'State.bin')if cp else None,native=identity(out/'native.json')if p else None,source_sha256=sha(__file__),
        state_size=len(r['state']),normal_cells=len(r['fields'][0]['mesh']),rank=1,threads=1,steps=0,
        governing_solves=0,completed_stage='NOT_APPLICABLE at CP0',native_initialization_work=r['initialization_lift_calls_iterations'],
        passive_initialization_work=p['initialization_work']if p else 'NOT_RUN:no accepted checkpoint',
        spatial_norms=p['physical_mass_dual_norms']if p else dict(native_startup_coupled_norms=event['audit']['coupled_norms'],
            completed_stage='NOT_APPLICABLE',scope='native startup alpha_mass=0, eta=0; reader-scaled packet unavailable'),
        native_audit=packet,principal_witness_replays=witnesses,represented=represent,interface=interface,independent_weak=weak,
        quadrature=quadrature,open_boundary=p['open_boundary']if p else r['open_boundary'],scalar_domain=p['scalar_domain']if p else None,
        literal_reference=c.summary(),original_bytes_unchanged=original)
    save(out/'result.json',result)
    print(json.dumps(dict(height_mm=height,status=result['status'],gates=gates,state_size=len(r['state']),represented=represent)),flush=True)


def plot(root):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out=root/'figures';out.mkdir(exist_ok=False);cal=load(root/'reference/calibration.json');cfg=configurations();fit=cal['fits']['primary']
    fig,ax=plt.subplots(1,2,figsize=(9,3.7),constrained_layout=True)
    rows=[r for r in measurements()if r['panel']=='a'];z=np.array([r['z_mm']for r in rows]);u=[r['value']for r in rows]
    ax[0].plot(u,z,'o',mfc='none',label='Fig.3(a) digitized rough symbols')
    zz=np.geomspace(8,800,400);ax[0].plot([cfg[1].U_outer(v/SCALE_MM-cfg[1].d0)for v in zz],zz,label='fixed-slope neutral log')
    ax[0].axhline(25,color='.4',ls=':',label='primary fit upper height');ax[0].set_ylim(5,800);ax[0].set_yscale('log')
    ax[0].set(xlabel=r'$U/U_\infty$',ylabel='physical height z [mm]',title='Mean fit and finite-depth mismatch');ax[0].legend(fontsize=7)
    for c in cfg:
        rr=np.linspace(0,c.H,300);ax[1].plot([c.k_inner(r)/c.K for r in rr],rr*SCALE_MM,label=f'Hw={c.H*SCALE_MM:g}mm')
    ax[1].set(xlabel=r'$k_w/K$',ylabel='effective-wall distance r [mm]',title='Nonconstant budget reference');ax[1].legend()
    for a in ax:a.grid(alpha=.2)
    fig.savefig(out/'mean-and-budget-profile.png',dpi=170);fig.savefig(out/'mean-and-budget-profile.pdf');plt.close(fig)
    fig,ax=plt.subplots(figsize=(7.8,3.5),constrained_layout=True)
    hs=np.linspace(3.1,6.3,200);required=[SCALE_MM*required_roughness(H/SCALE_MM,fit['z0_e'])['h']for H in hs]
    ax.plot(hs,required,label='required sand-grain h_s (primary mean fit)');ax.plot(hs,hs,'--',label='h_s=H_w boundary')
    ax.axvspan(3.1,3.5,color='tab:red',alpha=.12,label='all 7mm observations outside: H_w<3.5mm')
    ax.axvline(cal['feasibility']['fig3_conservative_ceiling_mm'],ls=':',color='.4',label='conservative Fig.3 ceiling')
    for c in cfg:ax.plot(c.H*SCALE_MM,c.h*SCALE_MM,'o',color='tab:blue')
    ax.set(xlabel='matching height H_w [mm from effective wall]',ylabel='sand-grain roughness h_s [mm]',title='Separate observation scopes: no all-7mm-outside solution')
    ax.grid(alpha=.2);ax.legend(fontsize=8);fig.savefig(out/'height-feasibility.png',dpi=170);fig.savefig(out/'height-feasibility.pdf');plt.close(fig)
    stress=[];c=cfg[1]
    for row in measurements():
        if row['panel']=='a':continue
        panel=row['panel'];model=c.u_tau**2 if panel=='b'else math.sqrt(c.Rc[0,0]if panel=='c'else c.Rc[2,2])/c.u_tau
        stress.append(dict(row,model=model,model_minus_digitized=model-row['value'],
            normalization='-R_xz/U_infinity^2'if panel=='b'else 'sigma_x/u_tau'if panel=='c'else 'sigma_vertical/u_tau',
            scope='exact neutral outer stress; Fig.3 samples above all three candidate interfaces'))
    write_csv(out/'neutral-stress-discrepancy.csv',stress)
    fig,ax=plt.subplots(1,2,figsize=(8.8,3.6),constrained_layout=True)
    for panel,label,color in (('c','longitudinal','tab:blue'),('d','vertical','tab:orange')):
        rows=[v for v in stress if v['panel']==panel]
        ax[0].plot([v['value']for v in rows],[v['z_mm']for v in rows],'o',mfc='none',color=color,label=f'Fig.3 {label}')
        ax[0].axvline(rows[0]['model'],color=color,ls='--',label=f'neutral {label}')
    ax[0].set(xlabel=r'$\sigma/u_\tau$',ylabel='physical height z [mm]',title='Normal-stress bias remains visible');ax[0].legend(fontsize=7)
    rows=[v for v in stress if v['panel']=='b'];ax[1].plot([v['value']for v in rows],[v['z_mm']for v in rows],'o',mfc='none',label='Fig.3 shear')
    ax[1].axvline(c.u_tau**2,ls='--',label='neutral constant shear');ax[1].set(xlabel=r'$-R_{xz}/U_\infty^2$',ylabel='physical height z [mm]',title='Idealized equilibrium is not a finite layer')
    ax[1].legend(fontsize=8)
    for a in ax:a.grid(alpha=.2)
    fig.savefig(out/'neutral-stress-discrepancy.png',dpi=170);fig.savefig(out/'neutral-stress-discrepancy.pdf');plt.close(fig)
    save(out/'seven-mm-model-regions.json',[dict(H_w_mm=c.H*SCALE_MM,physical_z_mm=7.,r_mm=7.-c.d0*SCALE_MM,
        region='inner AWF',R=c.covariance_inner(7/SCALE_MM-c.d0).tolist(),
        longitudinal_rms_over_utau=math.sqrt(c.covariance_inner(7/SCALE_MM-c.d0)[0,0])/c.u_tau,
        scope='reference prediction only, no Fig.12 rms value fitted or digitized')for c in cfg])
    print('Compact calibration/profile/stress plots written; Fig.12 7mm target retained inside AWF',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('verify','reference','generate','assess','plot'),required=True)
    p.add_argument('--kernel',type=Path);p.add_argument('--observer',type=Path);p.add_argument('--height-mm',type=float)
    a=p.parse_args()
    if a.mode=='verify':verify(a.root)
    elif a.mode=='reference':reference(a.root,a.kernel)
    elif a.mode=='generate':generate(a.root)
    elif a.mode=='assess':assess(a.root,a.height_mm,a.observer)
    else:plot(a.root)
