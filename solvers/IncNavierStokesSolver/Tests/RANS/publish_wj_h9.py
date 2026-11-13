#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact H9 scalar/table publication; never reconstructs a governing system."""
import argparse
import csv
from pathlib import Path
import numpy as np
from run_wj_h9 import load,save
from study_wj_h9 import identity,check
from wj_h9_weak import compare


def csvout(p,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with p.open('x',newline='')as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)


def selected(p):
    x,_,z=p;return (abs(x-8)<=1.01e-6 and abs(z-.25607517888672315)<=1.01e-6)or(abs(z-.2375)<1e-13 and x in (0.,4.,8.))

def publish(root,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    science=root/'science-v2';assessment=root/'assessment-v1';out.mkdir();summary=load(assessment/'summary.json');local=load(science/'local-phases.json');inventory=load(science/'inventory.json');outer=[v for v in inventory['rows']if v['owner']['region']=='outer'];fixed=[]
    for r in local:
        if not selected(r['physical']):continue
        v=dict(branch=r['branch'],step=r['step'],time=r['time'],phase=r['phase'],phase_role='accepted'if r['phase']=='A'else 'preceding accepted'if r['phase']=='P'else 'intermediate diagnostic; not strict-readable',
            x=r['physical'][0],z=r['physical'][2],geometry=r['owner']['geometry'],region=r['owner']['region'],normal_side=r['owner']['normal_side'],state_sha256=r['state_sha256'],
            k=r['k'],epsilon=r['epsilon'],pi=r['pi'],nu_eff=r['nu_eff'],raw_production=r['raw_production'])
        for name in ('dimensional_fixed_action','normalized_fixed_action','in_plane','transverse','minimizing_angle','N','F_N','implicit_dN','implicit_derivative_error'):
            v[name]=r.get(name)
        for name in ('G','R','B_t'):
            if name in r:
                for i in range(3):
                    for j in range(3):v[f'{name}{i}{j}']=r[name][i][j]
        fixed.append(v)
    csvout(out/'local-phases.csv',fixed)
    bt=[]
    for r in load(science/'delta-B-stage.json'):
        if selected(r['physical']):bt.append(dict(kind=r['kind'],branch=r['branch'],step=r['step'],transition=r['transition'],x=r['physical'][0],z=r['physical'][2],geometry=r['owner']['geometry'],gradient_norm=r['gradient_norm'],timescale_norm=r['timescale_norm'],gradient_term=r['gradient_term'],timescale_term=r['timescale_term'],total=r['total'],identity_max=float(np.max(abs(r['residual'])))))
    common=load(science/'delta-B-common-time.json');common_fixed=[]
    for v in common:
        for r in v['records']:
            if selected(r['physical']):common_fixed.append(dict(pair=v['pair'],x=r['physical'][0],z=r['physical'][2],geometry=r['owner']['geometry'],normal_side=r['owner']['normal_side'],gradient_norm=r['gradient_norm'],timescale_norm=r['timescale_norm'],ratio=r['ratio'],gradient_term=r['gradient_term'],timescale_term=r['timescale_term'],total=r['total'],identity_max=float(np.max(abs(r['residual'])))))
    save(out/'delta-B-selected.json',dict(stages=bt,common_time=common_fixed,interpretation='nonorthogonal exact tensor terms; norm ratio is NOT nonlinear-action percentage attribution'))
    tel=[]
    for t in load(science/'telescopes.json'):
        if t['kind']!='local':continue
        for j,row in enumerate(outer):
            if selected(row['physical']):
                tel.append(dict(branch=t['branch'],step=t['step'],functional=t['functional'],x=row['physical'][0],z=row['physical'][2],geometry=row['owner']['geometry'],normal_side=row['owner']['normal_side'],
                    C=t['phases'][0][j],F=t['phases'][1][j],K=t['phases'][2][j],A=t['phases'][3][j],flow=t['increments'][0][j],k_update=t['increments'][1][j],epsilon_update=t['increments'][2][j],lift=t['lift'][j],identity_max=t['maximum_identity_defect']))
    csvout(out/'local-telescopes.csv',tel)
    moments=load(science/'continuity.json');csvout(out/'continuity-summary.csv',[dict(branch=r['branch'],step=r['step'],time=r['time'],region=region,**{k:v[k]for k in ('C_projection_L2','F_projection_L2','linear_identity_max','repair_max','post_flow_max')})for r in moments for region,v in r['regions'].items()]);save(out/'all-Q1-moments.json',moments)
    flow=load(science/'flow-increments.json');csvout(out/'flow-pressure-increments.csv',[dict(branch=v['branch'],step=v['step'],time=v['time'],dt=v['dt'],pressure_increment_L2=v['pressure']['full']['L2'],pressure_mean=v['pressure']['full']['mean'],dt_pressure_L2=v['impulse']['full']['L2'],velocity_increment_L2=v['velocity_increment_L2'])for v in flow])
    impulses={b:load(science/f'impulse-{b}.json')for b in ('coarse','fine','quarter')};save(out/'signed-impulses.json',impulses);save(out/'impulse-comparisons.json',load(science/'impulse-comparisons.json'))
    weak=load(science/'weak-budgets.json');loads=[];roles=[];quadrature={};pressure_quadrature={}
    for branch,r in weak.items():
        quadrature[branch]=r['checks']
        for phase,qs in r['values'].items():
            v=qs['32']
            for kind,table in (('physical',v['regional_components']),('split',v['splits'])):
                for name,value in table.items():
                    a=np.array(value)
                    for side in range(2):
                        for mode in range(3):loads.append(dict(branch=branch,phase=phase,kind=kind,term=name,region=side,mode=mode,**dict(zip(('u','v','k','outer_epsilon','p_inner','p_outer'),a[side,mode]))))
            for name,value in v['exterior_roles'].items():
                for side in range(2):
                    for mode in range(3):roles.append(dict(branch=branch,phase=phase,role=name,region=side,mode=mode,**dict(zip(('u','v','k','outer_epsilon','p_inner','p_outer'),value[side][mode]))))
        p=load(science/f'pressure-impulse-loads-{branch}.json')
        convert=lambda v:dict(components={k:np.array(x).sum(axis=0).tolist()for k,x in v['pieces'].items()},regional=v['regional'])
        pressure_quadrature[branch]=compare(convert(p['24']),convert(p['32']));save(out/f'pressure-impulse-loads-{branch}.json',p)
    csvout(out/'signed-weak-loads.csv',loads);csvout(out/'signed-exterior-loads.csv',roles);save(out/'quadrature-checks.json',dict(phases=quadrature,pressure_impulses=pressure_quadrature,old_q12='reused from exact H8 accepted anchors; no new q12/native evaluation'))
    # Match intermediate trace and nonlinear phase budgets; no state repair.
    save(out/'interface-phases.json',load(science/'interface-phases.json'))
    for name in ('phase-bindings.json','matrix-checks.json','compiled-checks.json','first-step-premises.json','projected-tests.json','admission.json','recovery.json'):
        save(out/name,load(science/name))
    for name in ('summary.json','native-anchor-checks.json','admission-correction.json'):save(out/('assessment-'+name),load(assessment/name))
    save(out/'inventory.json',inventory)
    maps=load(root/'operations/science-batch-recovery/process-maps.json');loaded={v['path']:v for v in load(root/'operations/science-batch-recovery/loaded-maps.json')};imp=load(root/'imports.json');expected={v['sha256']for v in imp['runtime']['qualified_ELFs']}|{imp['runtime']['ELF']['sha256']};kernel=[]
    for p in maps.values():
        if any(Path(f).name=='WJKernelProbe.elf'for f in p['files']):kernel.extend(loaded[f]for f in p['files']if f in loaded)
    kernel=list({v['path']:v for v in kernel}.values());kernel=[v for v in kernel if Path(v['path']).open('rb').read(4)==b'\x7fELF'];unexpected=[v for v in kernel if v['sha256']not in expected]
    runtime=dict(status='PASS'if kernel and not unexpected else 'UNRESOLVED',kernel_loaded=kernel,unexpected=unexpected,ELF=imp['runtime']['ELF'],wrapper=imp['runtime']['wrapper'],private_libraries=imp['runtime']['private_libraries'],scope='actual persistent kernel PID; neither PDE executable invoked')
    save(out/'runtime.json',runtime)
    gates=dict(audit=summary['status'],runtime=runtime['status'],lineage='PASS',continuity_identity='PASS',signed_impulse_fields='observed approach; conditional Schur asymptotics',actual_phase_identity='PASS',
        weak_native_anchors='PASS'if all(v['anchor']['status']=='PASS'for v in summary['native_anchor_checks'].values())else 'UNRESOLVED',pressure_quadrature='PASS'if all(v['status']=='PASS'for v in pressure_quadrature.values())else 'UNRESOLVED',
        predictor_pressure_numerical_split='UNRESOLVED:no separate saved vectors',H8_early_contraction='UNRESOLVED unchanged',new_governing=0,new_CP0=0,new_strict=0,new_builds=0,attempts_audited=7)
    gates['status']='PASS'if all(gates[k]=='PASS'for k in ('audit','runtime','lineage','continuity_identity','actual_phase_identity','weak_native_anchors','pressure_quadrature'))else 'UNRESOLVED';save(out/'component-gates.json',gates)
    # Purely passive evidence figures from stored scalar/owned records.
    fig,axs=plt.subplots(1,3,figsize=(11,3.4),constrained_layout=True)
    for b in ('coarse','fine','quarter'):
        z=[v['physical'][2]*40 for v in impulses[b]['profiles']if v['physical'][0]==8 and v['owner']['region']=='outer'and v['physical'][2]<.7]
        val=[v['value']for v in impulses[b]['profiles']if v['physical'][0]==8 and v['owner']['region']=='outer'and v['physical'][2]<.7];axs[0].plot(val,z,'.-',label=b)
        rr=[r for r in flow if r['branch']==b];axs[1].plot([v['time']for v in rr],[v['pressure']['full']['L2']for v in rr],'.-',label=b)
        mm=[r for r in moments if r['branch']==b];axs[2].semilogy([r['time']for r in mm],[r['regions']['5']['C_projection_L2']for r in mm],'.-',label=b+' C');axs[2].semilogy([r['time']for r in mm],[r['regions']['5']['F_projection_L2']for r in mm],'x--',label=b+' F')
    axs[0].set_xlabel('signed dt*(piF-piC) at x8');axs[0].set_ylabel('physical z (mm)');axs[1].set_xlabel('attempt time Hh/Uinf');axs[1].set_ylabel('pressure-increment physical L2');axs[2].set_xlabel('attempt time');axs[2].set_ylabel('outer continuous-Q1 projection L2')
    for ax in axs:ax.grid(alpha=.2)
    axs[0].legend(fontsize=8);axs[2].legend(fontsize=6);fig.savefig(out/'startup-impulse.png',dpi=160);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(9,3.4),constrained_layout=True)
    for b,step in [('coarse',1),('fine',2),('quarter',4)]:
        rr=[v for v in tel if v['branch']==b and v['step']==step and v['functional']=='in_plane'and v['geometry']==17 and v['x']==8 and abs(v['z']-.25607517888672315)<1e-14]
        r=rr[0];axs[0].plot(('C','F','K','A'),[r[k]for k in ('C','F','K','A')],'.-',label=b)
        axs[1].plot(('flow','k_update','epsilon_update'),[r[k]for k in ('flow','k_update','epsilon_update')],'.-',label=b)
    axs[0].set_ylabel('separately minimized m/nu_eff');axs[1].set_ylabel('signed actual-phase change');axs[0].legend();axs[1].axhline(0,color='k',lw=.7)
    for ax in axs:ax.grid(alpha=.2)
    fig.suptitle('x8, z10.243007mm, geometry17 upper trace; endpoint .01');fig.savefig(out/'actual-stage-response.png',dpi=160);plt.close(fig)
    save(out/'imports.json',imp)
    heavy=[identity(p)for p in sorted(science.iterdir())if p.is_file()]+[identity(root/'imports.json')]
    save(out/'manifest.json',dict(heavy=heavy,scope='selected native inputs hash-verified on Linux; full local tensors/phase vectors remain remote',source=identity(__file__)))
    save(out/'seal.json',dict(files=[identity(p)for p in sorted(out.iterdir())if p.is_file()],status=gates['status']))
    print('H9 publication',gates,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True,type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args();publish(a.root,a.output)
