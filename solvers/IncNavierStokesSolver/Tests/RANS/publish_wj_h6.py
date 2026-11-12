#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact H6 scientific tables/plots, retaining heavy reconstructed fields."""
import argparse
from pathlib import Path
import shutil
import subprocess
import numpy as np
from run_wj_h6 import load,save,sha,THREADS
from study_wj_h5 import identity,write_csv
from study_wj_h6 import weak_status,record,configuration
from wj_h5_neutral import KAPPA


def publish(root):
    out=root/'publication';out.mkdir(exist_ok=False);staged=[];quality=[];models=[];heavy=[];checks={};witnesses={};interfaces={};residuals={};comparisons={}
    for height in (5.,5.5,6.):
        key=f'{height:g}';c=configuration(height);r,rp,xml=record(root,height)
        n_i=sum(e['upper'][2]<=c.H for e in r['fields'][0]['mesh']);n_o=len(r['fields'][0]['mesh'])-n_i
        assert (n_i,n_o,len(r['state']))==(16,34,2328)
        assert not r['failure'] and r['steps']==0 and r['time']==0 and not r['pimple_attempts'] and r['implicit_callback_count']==0
        observation=load(root/f'observations/H{key}mm/result.json');comparison=load(root/f'comparisons/H{key}mm/result.json')
        comparisons[key]=comparison
        weak=weak_status(observation['native']is not None,observation['independent_weak']['status'])
        stages={}
        for role in ('coarse','refined'):
            folder=root/f'replays-v2/{role}/H{key}mm';v=load(folder/'result.json');m=load(folder/'reconstruction.json');assert v['status']=='PASS'
            stages[role]=dict(coefficient_max_scaled=max(x['scaled']for x in m['coefficient_errors']),
                field_max_scaled=v['field_replay_max'],physical_G_max=v['physical_G_replay_max'],
                backward_residual_max=max(x['backward_residual']for x in m['physical_linear_systems']),
                XML_errors=m['XML_continuum_max_errors'],lineage=m['native_lineage'],
                reconstruction=identity(folder/'reconstruction.json'),result=identity(folder/'result.json'))
            heavy.extend([identity(folder/'stage-rows.json'),identity(folder/'reconstructed-stages.npz')])
            for w in v['witnesses']:
                for label,row in w['stages'].items():
                    native=w['native'];rphys=native['position'][2]-c.d0;g=c.u_tau/(c.summary()['kappa']*rphys)
                    staged.append(dict(H_w_mm=height,mesh=role,witness_kind=w['kind'],stage=label,
                        x=native['position'][0],z_mm=native['position'][2]*40,geometry=native['geometry'],
                        Gxz=row['G'][0][2],Gxx=row['G'][0][0],Gzz=row['G'][2][2],
                        k=row['k'],epsilon=row['epsilon'],tau=row['tau'],
                        shear_error=row['G'][0][2]/g-1,epsilon_error=row['epsilon']/c.epsilon(rphys)-1,
                        m_in_plane=row['principal']['in_plane'],m_transverse=row['principal']['transverse'],
                        min_covariance=min(row['covariance']),raw_work=row['Praw'],N=row['N'],root_residual=row['root_residual']))
            witnesses[f'{role}-{key}']=v['witnesses']
        A=bool(all(max(v['XML_errors'].values())<2e-11 for v in stages.values()))
        B=bool(all(v['passed']is True for v in observation['native_audit']['conditions'].values()))
        D='PASS'if weak=='PASS'and all(v['status']=='PASS'for v in observation['quadrature'].values())and observation['gates']['exact_shared_traces']else 'UNRESOLVED'
        checks[key]=dict(reference_A='PASS'if A else 'FAIL',native_B='PASS'if B else 'OBSTRUCTED',
            representation_C=comparison['status'],weak_interface_D=D,executed_weak=weak,
            reconstruction=stages,common_inventory=comparison['inventory'],
            status='PASS'if A and B and comparison['status']=='PASS'and D=='PASS'else 'OBSTRUCTED'if not B or comparison['status']=='OBSTRUCTED'else 'UNRESOLVED')
        for mesh,v in comparison['summaries'].items():
            quality.append(dict(H_w_mm=height,mesh=mesh,**v['errors'],min_in_plane=v['min_in_plane'],
                min_transverse=v['min_transverse'],min_covariance=v['min_covariance'],
                exact_margin=v['exact_min'],margin_fraction=v['retained_fraction'],outer_owner_rows=v['outer_owner_rows']))
        interfaces[key]=observation['interface'];residuals[key]=dict(H6=observation['spatial_norms'],
            H5=load(root.parent/'wj-h5-equilibrium-inlet-20261112/observations/H5mm/result.json')['spatial_norms']if height==5 else None,
            H5_scope='comparable accepted-reader mass-dual packet'if height==5 else 'NOT_AVAILABLE: rejected H5 startup coupled_norms use different normalization')
        models.append(dict(H_w_mm=height,Ni=n_i,No=n_o,state_size=len(r['state']),
            native_min_in_plane=observation['native_audit']['principal_witnesses']['in_plane']['minimum_over_nu_eff'],
            native_min_transverse=observation['native_audit']['principal_witnesses']['transverse']['minimum_over_nu_eff'],
            status=checks[key]['status'],weak=weak,governing_steps=0,governing_solves=0,
            max_traction_component_jump=max(max(abs(z)for z in v['traction_jump'])for v in observation['interface']),
            max_k_flux_jump=max(abs(v['k_flux_jump'])for v in observation['interface']),
            max_epsilon_trace_defect=max(abs(v['epsilon_trace']-v['epsilon_target'])for v in observation['interface'])))
        for path in (rp,xml,xml.with_name('case-design.json'),root/f'observations/H{key}mm/result.json',root/f'observations/H{key}mm/native.json',root/f'comparisons/H{key}mm/result.json',root/f'comparisons/H{key}mm/inventory.json'):
            heavy.append(identity(path))
        heavy.extend([identity(Path(r['checkpoint'])/'State.bin'),*comparison['rows'].values()])
    assert len(models)==3
    write_csv(out/'stage-witnesses.csv',staged);write_csv(out/'quality.csv',quality);write_csv(out/'native-startups.csv',models)
    overall='PASS'if all(v['status']=='PASS'for v in checks.values())else 'OBSTRUCTED'if any(v['status']=='OBSTRUCTED'for v in checks.values())else 'UNRESOLVED'
    save(out/'component-gates.json',dict(status=overall,cases=checks,
        governing_timesteps=0,extra_initializations=0,scope='single prescribed bisection and finite physical/owner inventory, no temporal or global stability claim'))
    save(out/'principal-witnesses.json',witnesses);save(out/'interface-traces.json',interfaces);save(out/'spatial-residuals.json',residuals)
    save(out/'comparisons.json',comparisons)
    shutil.copyfile(root/'imports.json',out/'imports.json');shutil.copyfile(root/'run-intent.json',out/'run-intent.json')
    save(out/'weak-quadrature.json',{key:{k:load(root/f'observations/H{key}mm/result.json')[k]for k in ('independent_weak','quadrature','scalar_domain','native_audit')}for key in checks})
    plot(root,out,staged,quality)
    receipts={};logs=out/'logs';logs.mkdir()
    for path in sorted((root/'operations').glob('*/receipt.json')):
        receipts[path.parent.name]=dict(load(path),receipt=identity(path));shutil.copyfile(path.with_name('command.log'),logs/(path.parent.name+'.log'))
    save(out/'commands.json',receipts)
    shutil.copyfile(root/'accounting.json',out/'accounting-reserved.json');shutil.copyfile(root/'setup.json',out/'setup.json')
    shutil.copyfile(root/'setup.log',logs/'failed-source-fetch.log')
    eq=load(root/'runtime-equivalence.json');mapped=[]
    for label in ('CP0-H5mm','CP0-H5.5mm','CP0-H6mm','H6-observe-5','H6-observe-5.5','H6-observe-6'):
        path=root/'operations'/label/'loaded-maps.json';v=load(path)
        spatial=[x for x in v if Path(x['path']).name=='libSpatialDomains.so.5.9.0'];assert len(spatial)==1
        assert spatial[0]['sha256']=='0c23e0de180a233bd88f70b265ae9865aeac456e6482f8ff0fa1b4be84d69dab'
        mapped.append(dict(operation=label,SpatialDomains=spatial[0],maps=identity(path)))
    save(out/'runtime.json',dict(governing_producer=eq['governing_producer'],observer_producer=eq['observer_producer'],
        executables=eq['executables'],verified_native_inputs=len(eq['compiled_inputs']),verified_ELFs=len(eq['qualified_ELFs']),
        mappings=mapped,threads=THREADS,arithmetic_policy=eq['reused_arithmetic_policy'],native_rebuild=False,
        original_equivalence=identity(root/'runtime-equivalence.json')))
    revision=lambda v:subprocess.check_output(['git','rev-parse',v],cwd=Path(__file__).parent,text=True).strip()
    save(out/'manifest.json',dict(root=str(root),baseline=revision('34937bb'),science_source=revision('14bf437'),publisher_source=revision('HEAD'),
        heavy=heavy,specifications=eq['specifications'],H5_parent='sealed H5 artifacts imported read-only; rejected candidates have no accepted checkpoint',
        physical_reference='WJNeutralBudgetReferenceV1 unchanged',publication_accounting='reserved snapshot, final closure receipt is separate'))
    save(out/'seal.json',dict(files=[dict(name=str(p.relative_to(out)),**identity(p))for p in sorted(out.rglob('*'))if p.is_file()],
        source='native startups + independent reconstructed stages and fixed comparison inventory'))
    print(models,flush=True)


def plot(root,out,staged,quality):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(11,3.5),constrained_layout=True)
    names=('physical_L2','exterior_lift','dependent_epsilon_lift','native')
    for ax,H in zip(axes,(5.,5.5,6.)):
        rows=[v for v in staged if v['H_w_mm']==H and v['mesh']=='coarse'and v['witness_kind']=='in_plane']
        ax.bar(np.arange(4),[next(v for v in rows if v['stage']==k)['m_in_plane']for k in names],color=['.6','tab:orange','tab:green','tab:blue'])
        exact=next(v['exact_margin']for v in quality if v['H_w_mm']==H and v['mesh']=='H6')
        ax.axhline(exact,color='k',ls=':',label='exact neutral')
        ax.axhline(0,color='.3',lw=.6);ax.set_xticks(np.arange(4),['L2','exterior','epsilon','native'],rotation=20)
        ax.set(title=f'H_w={H:g}mm, H5 limiting point',ylabel=r'$m/\nu_{eff}$');ax.grid(axis='y',alpha=.2)
    axes[0].legend(fontsize=8);fig.savefig(out/'stage-attribution.png',dpi=170);fig.savefig(out/'stage-attribution.pdf');plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(11,6.3),constrained_layout=True)
    for j,H in enumerate((5.,5.5,6.)):
        c=configuration(H);exact=next(v['exact_margin']for v in quality if v['H_w_mm']==H and v['mesh']=='H6')
        for label,color in (('H5','tab:orange'),('H6','tab:blue')):
            v=load(root/f'comparisons/H{H:g}mm/{label}-rows.json')
            rows=sorted((x for x in v if x['physical'][0]==8 and x['physical'][2]*40<=30),key=lambda x:x['physical'][2])
            z=[x['physical'][2]*40 for x in rows]
            axes[0,j].plot(z,[100*(x['G'][0][2]/(c.u_tau/KAPPA/(x['physical'][2]-c.d0))-1)for x in rows],'.-',ms=2,lw=.7,color=color,label=label)
            axes[1,j].plot(z,[x['principal']['in_plane']for x in rows],'.-',ms=2,lw=.7,color=color,label=label)
        axes[0,j].set(title=f'H_w={H:g}mm, x=8',ylabel='Gxz error [%]')
        axes[1,j].set(xlabel='physical height z [mm]',ylabel=r'$m/\nu_{eff}$')
        axes[1,j].axhline(.5*exact,color='k',ls=':',label='half-reference gate')
        for a in axes[:,j]:a.grid(alpha=.2)
    axes[0,0].legend();axes[1,0].legend(fontsize=8)
    fig.savefig(out/'normal-refinement.png',dpi=170);fig.savefig(out/'normal-refinement.pdf');plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();publish(a.root)
