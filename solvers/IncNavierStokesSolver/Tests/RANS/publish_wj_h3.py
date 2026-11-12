#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact H3 tables/figures and immutable heavy-data index; no evolution."""
import argparse
import csv
import hashlib
import math
from pathlib import Path
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_wj_h3 import load,save,sha
from wj_h3_analysis import actual,identity


def csvfile(path,rows):
    with path.open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();root=a.root
    out=root/'publication-v1';out.mkdir(exist_ok=False)
    summary=load(root/'analysis-v2/summary.json');cases={k:load(v['path']) for k,v in summary['cases'].items()}
    series=[];interfaces=[];local=[];witnesses={};runs=[];heavy={}
    def index(path):
        path=Path(path)
        if path.is_file():heavy[str(path)]=identity(path)
    for key,d in cases.items():
        index(d['producer']['path']);index(d['xml']['path']);index(root/'analysis-v2'/(key+'.json'))
        for s in d['states']:
            w=s['packet']['principal_witnesses']['in_plane'];G=np.array(w['raw_G_row_major']).reshape(3,3)
            series.append(dict(case=key,time=s['time'],role=s['role'],absolute_step=s.get('absolute_step','candidate'),
                m_over_nu_eff=w['minimum_over_nu_eff'],x=w['position'][0],z=w['position'][2],geometry=w['geometry'],
                Gxx=G[0,0],Gxz=G[0,2],Gzx=G[2,0],Gzz=G[2,2],divergence=float(np.trace(G)),
                k=w['k'],epsilon=w['epsilon'],tau=w['k']/w['epsilon'],
                covariance_min=min(w['covariance_eigenvalues']),nu_eff=w['nu_eff'],nu_d=w['nu_d'],root=w['N'],
                state_sha256=s['state_sha256'],callback_history_sha256=s['callback_history_sha256']))
            for v in s['interface']:
                interfaces.append(dict(case=key,time=s['time'],role=s['role'],x=v['x'],inner_geometry=v['geometry_inner_outer'][0],
                    outer_geometry=v['geometry_inner_outer'][1],K=v['K'],cutoff=v['cutoff'],pi_inner=v['pi_inner'],pi_outer=v['pi_outer'],
                    normal_traction_defect=v['traction_defect'][2],tangent_traction_defect=v['traction_defect'][0],
                    k_flux_defect=v['k_flux_defect'],mass_flux_defect=v['mass_flux_defect'],epsilon_defect=v['epsilon_defect'],
                    roughness_force_x=v['roughness_force_x'],inner_production=v['inner_production']))
            for v in s['rows']:
                if abs(v['physical'][2]-.6720111978440669)>1e-9:continue
                local.append(dict(case=key,time=s['time'],role=s['role'],x=v['physical'][0],z=v['physical'][2],geometry=v['geometry'],
                    Gxz=v['G'][0][2],Gxx=v['G'][0][0],Gzx=v['G'][2][0],Gzz=v['G'][2][2],tau=v['k']/v['epsilon'],
                    principal=v['principal']['in_plane'],divergence=v['divergence'],Praw=v['Praw'],Pdev=v['Pdev']))
            if s['role']=='rejected candidate':witnesses[key]=s
    csvfile(out/'principal-history.csv',series);csvfile(out/'interface-history.csv',interfaces);csvfile(out/'failing-face.csv',local)
    save(out/'adverse-witnesses.json',witnesses)
    for f in sorted((root/'cases').glob('*/result.json')):
        d=load(f);r=load(f.with_name('case.evolution.json'));index(f);index(f.with_name('case-design.json'));index(f.with_name('intent.json'))
        index(f.with_name('case.evolution.json'));index(f.with_name('case.xml'))
        snapshots=[]
        for attempt in r['pimple_attempts']:
            path=actual(attempt['k_snapshot']);snapshot=load(path)
            assert np.asarray(snapshot['raw'],dtype='<f8').tobytes()==np.asarray(snapshot['corrected'],dtype='<f8').tobytes()
            index(path);snapshots.append(identity(path))
        for cp in f.parent.glob('CP*/State.bin'):index(cp)
        runs.append(dict(case=f.parent.name,status=d['status'],state_size=d['state_size'],dt=d['intent']['dt'],
            seed=d['intent']['seed'],outlet=d['intent']['outlet'],Nx=d['intent']['nx'],q=d['intent']['q'],
            accepted_steps=d['accepted_step'],accepted_time=d['time'],attempts=d['attempts'],solve_counts=d['solve_counts'],
            corrections=d['corrections'],raw_corrected_k_exact=True,snapshots=snapshots,producer=d['report'],failure=d['failure']))
    attempts=sum(r['attempts'] for r in runs);accepted=sum(r['accepted_steps'] for r in runs)
    assert attempts==43 and accepted==42 and all(r['corrections']==0 for r in runs)
    save(out/'runs.json',runs)
    obs=load(root/'native-observations-v1/index.json');save(out/'native-observations.json',obs)
    for v in obs.values():
        index(v['native']['path']);index(v['checkpoint']['path'])
    budgets={}
    for key,v in obs.items():
        if 'independent_weak' not in v:continue
        p=load(v['native']['path']);rp=load(v['producer']['path'])
        from wj_h2_profiles import wall,NU,QE
        pref=destruction=production=0.
        for row in p['weak_checks']['samples']:
            if not row['inner']:continue
            K=row['K'];r=row['chart'][2];w=wall(K,rp['parameters'][1],rp['parameters'][2]);nt=w['a']*max(r-w['zv'],0)
            M=np.array(row['G'])+np.array(row['G']).T
            pref+=row['weight']*(K/QE)**2*nt/(NU+nt)**2
            destruction+=row['weight']*w['A']/max(r,w['rc'])
            production+=row['weight']*(.5*nt*np.sum(M*M)-2*row['k']*np.trace(row['G'])/3)
        ground=[r for r in p['open_boundary']['boundary_rows'] if r['role']=='ground']
        budgets[key]=dict(actual_inner_production=production,current_K_counterfactual=pref,destruction=destruction,
            reference_budget_difference=destruction-pref,
            signed_mean_wall_stress=-sum(r['weight']*r['actual_traction'][0] for r in ground)/sum(r['weight'] for r in ground),
            wall_stress_scope='positive fluid-on-wall x shear; opposite of outward fluid traction',
            mass_flux=p['open_boundary']['mass_flux'],boundary_load=p['open_boundary']['boundary_load'],
            independent_components=v['independent_weak']['independent_components'],
            regional_independent=v['independent_weak']['regional_independent'],
            regional_native=v['independent_weak']['regional_native'],gate=v['independent_weak']['status'])
    save(out/'weak-budgets.json',budgets)
    retained=load(root/'retained-v1/summary.json');save(out/'retained-summary.json',retained)
    selection=load(root/'quadrature-selection.json');save(out/'quadrature-selection.json',selection)
    qtable=[]
    for key,v in retained['states'].items():
        d=load(v['path']);index(v['path'])
        for pair,c in d['quadrature']['comparisons'].items():qtable.append(dict(case=key,pair=pair,status=c['status'],
            outer_stress_error=c['component_errors']['outer_momentum_stress'],regional_error=c['regional_max'],worst_element=c['worst_element'][1]))
    for key,d in cases.items():
        for role,q in d.get('terminal_quadratures',{}).items():
            for pair,c in q['comparisons'].items():qtable.append(dict(case=key+'-'+role,pair=pair,status=c['status'],
                outer_stress_error=c['component_errors']['outer_momentum_stress'],regional_error=c['regional_max'],worst_element=c['worst_element'][1]))
    csvfile(out/'quadrature.csv',qtable)
    # Same-time matched points include both face owners. Bulk norms cannot stand
    # in for the local gradient and principal comparisons below.
    contrasts=[]
    pairs=[('outlet-x1','B-zero-x1','B-frozen-x1'),('outlet-x3','B-zero-x3','B-frozen-x3'),
        ('horizontal-zero','B-zero-x1','B-zero-x3'),('horizontal-frozen','B-frozen-x1','B-frozen-x3')]
    for label,ka,kb in pairs:
        owner={}
        for key in (ka,kb):
            rr=load(cases[key]['producer']['path'])
            owner[key]={e['geometry']:e for e in rr['fields'][3]['mesh']}
        for sa,sb in zip(cases[ka]['states'],cases[kb]['states']):
            assert sa['time']==sb['time']
            # Exclude extra horizontal-face ownership only from this common
            # point pairing, while retaining all owners in the full observations.
            def get(s,key):
                result={}
                for v in s['rows']:
                    if 'native_witness' in v:continue
                    e=owner[key][v['geometry']];q=(v['physical'][0],v['physical'][2],e['lower'][2],e['upper'][2])
                    result.setdefault(q,[]).append(v)
                return result
            aa,bb=get(sa,ka),get(sb,kb);shared=sorted(aa.keys()&bb.keys())
            for x,z,lo,hi in shared:
                for va in aa[(x,z,lo,hi)]:
                    for vb in bb[(x,z,lo,hi)]:
                        ea,eb=owner[ka][va['geometry']],owner[kb][vb['geometry']]
                        if len(aa[(x,z,lo,hi)])==len(bb[(x,z,lo,hi)])==2 and ea['lower'][0]!=eb['lower'][0]:continue
                        contrasts.append(dict(pair=label,time=sa['time'],x=x,z=z,delta_Gxz=vb['G'][0][2]-va['G'][0][2],
                            geometry_a=va['geometry'],geometry_b=vb['geometry'],normal_owner_lo=lo,normal_owner_hi=hi,
                            delta_principal=vb['principal']['in_plane']-va['principal']['in_plane'],
                            principal_a=va['principal']['in_plane'],principal_b=vb['principal']['in_plane'],
                            delta_divergence=vb['divergence']-va['divergence'],delta_tau=vb['k']/vb['epsilon']-va['k']/va['epsilon'],
                            delta_pi=vb['pi']-va['pi'],delta_Ux=vb['U'][0]-va['U'][0],delta_Uz=vb['U'][2]-va['U'][2],
                            delta_k=vb['k']-va['k'],delta_epsilon=vb['epsilon']-va['epsilon']))
    csvfile(out/'paired-local.csv',contrasts)
    # Scientific figures use native finite-inventory minima and independently
    # reconstructed owned-face samples, with explicit dimensional scopes.
    plt.rcParams.update({'font.size':9,'figure.dpi':150})
    fig,ax=plt.subplots(figsize=(7,3.8))
    for key in ('H2-refined-pilot','A-dt02','A-dt01','B-zero-x1','B-frozen-x1','B-zero-x3','B-frozen-x3'):
        rows=[v for v in series if v['case']==key];ax.plot([v['time'] for v in rows],[v['m_over_nu_eff'] for v in rows],'-o',ms=3,label=key)
        for v in rows:
            if v['role']=='rejected candidate':ax.plot(v['time'],v['m_over_nu_eff'],'rx',ms=9)
    ax.axhline(0,color='k',lw=.8);ax.set(xlabel='T = t U∞/Hh',ylabel='native sampled principal / νeff',title='Fresh H3 short controls; rejected candidates are crosses')
    ax.legend(ncol=2,fontsize=7);fig.tight_layout();fig.savefig(out/'principal-history.png');plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(8,3.4))
    for key,t,role in (('H2-refined-pilot',0.,'accepted state'),('H2-refined-pilot',.04,'accepted state'),('H2-refined-pilot',.06,'rejected candidate')):
        for geom in (17,18):
            rr=sorted([v for v in local if v['case']==key and v['time']==t and v['role']==role and v['geometry']==geom],key=lambda v:v['x'])
            if not rr:continue
            label=f'T={t:g}, owner {geom}'
            axes[0].plot([v['x'] for v in rr],[v['Gxz'] for v in rr],'-o',ms=2,label=label)
            axes[1].plot([v['x'] for v in rr],[v['principal'] for v in rr],'-o',ms=2,label=label)
    axes[0].axhline(.24047947338930892,color='grey',ls=':',label='continuous prescribed shear')
    axes[1].axhline(0,color='k',lw=.8)
    for ax in axes:ax.set_xlabel('x/Hh');ax.legend(fontsize=6)
    axes[0].set_ylabel('Gxz Hh/U∞');axes[1].set_ylabel('principal / νeff')
    fig.suptitle('Retained H2 normal-face traces at z/Hh=.6720112; separate owners')
    fig.tight_layout();fig.savefig(out/'retained-face.png');plt.close(fig)
    fig,ax=plt.subplots(figsize=(6,3.4))
    for key in ('refined-cp0','rejected06','half08'):
        rr=[v for v in qtable if v['case']==key and v['pair'] in ('8-32','12-32','16-32','24-32')]
        ax.semilogy([int(v['pair'].split('-')[0]) for v in rr],[v['outer_stress_error'] for v in rr],'-o',label=key)
    ax.axhline(2e-9,color='k',ls=':',label='absolute component floor (relative term also tested)')
    ax.set(xlabel='passive integration order',ylabel='max outer-stress load difference to q32',title='Original fields; fixed regions, modes and physical Jacobians')
    ax.legend(fontsize=7);fig.tight_layout();fig.savefig(out/'quadrature.png');plt.close(fig)
    startup=load(root/'matched-startup-observation/result.json');save(out/'matched-startup.json',startup)
    h2=actual('/media/jay/新加卷/nektar-v2026-runs/wj-h2-20261111')
    entry=load(h2/'observations-v3/index.json')['profiles-startup'];oldpath=actual(entry['native'])
    assert sha(oldpath)==entry['native_sha256'];index(oldpath)
    oldnorm=load(oldpath)['physical_mass_dual_norms']['uncancelled_spatial'][1]
    vals=[oldnorm,startup['minimal_offset_spatial']['uncancelled_spatial'][1],startup['balanced_spatial']['uncancelled_spatial'][1]]
    fig,ax=plt.subplots(figsize=(6,3.2));ax.bar(['H2 wrong offset','H3 corrected offset','H3 vertical balance'],vals)
    ax.set_yscale('log');ax.set(ylabel='vertical spatial residual, consistent-mass dual',title='Identical base mesh/q8 and differential initial coefficients')
    fig.tight_layout();fig.savefig(out/'startup-balance.png');plt.close(fig)
    # Retain exact command receipts, including failed setup attempts and the
    # intended native hard exit. No receipt is rewritten as a successful run.
    commands=[]
    for f in sorted((root/'operations').glob('*/receipt.json')):
        d=load(f);commands.append(d);index(f);index(f.with_name('request.json'));index(f.with_name('command.log'))
    save(out/'commands.json',commands)
    for pth in (root/'receipts').glob('*'):index(pth)
    for pth in (root/'matched-startup-observation').glob('*'):index(pth)
    index(root/'seed-signature-rejection/result.json');index(root/'short-control-admission-v1.json')
    source=Path(__file__).resolve().parents[4]
    index(source/'docs/plans/wj-h2-linux/continuation-denial.json')
    save(out/'preedit-verification.json',load(root/'receipts/preedit-verification.json'))
    inventory=load(root/'receipts/focused-test-inventory.json')
    save(out/'focused-tests.json',{k:inventory[k] for k in ('registrations','selected','selected_count')})
    save(out/'heavy-index.json',heavy)
    gates=dict(schema='WJH3ComponentGatesV1',initializer_sign='PASS',old_sign_fixture='EXPECTED_FAIL_REPRODUCED',
        balanced_continuous_normal_balance='PASS',represented_balance='MEASURED_PROJECTION_ERROR',
        completion_label_only='PASS',retained_principal_replay='PASS',quadrature='PASS',q_star=12,
        sign_only_coarse_control='OBSTRUCTED',sign_only_half_dt_control='PASS_SHORT_TRANSIENT',
        predetermined_B={k:'PASS_SHORT_TRANSIENT' for k in ('B-zero-x1','B-frozen-x1','B-zero-x3','B-frozen-x3')},
        evolved_B_loads='PASS' if all(v.get('independent_weak',{}).get('status','PASS')=='PASS' for v in obs.values()) else 'UNRESOLVED',
        attempts=attempts,accepted_steps=accepted,rejected_attempts=attempts-accepted,
        solve_counts=[sum(v['solve_counts'][i] for v in runs) for i in range(5)],correction_activations=0,
        developed_precursor='OBSTRUCTED_NOT_QUALIFIED',H2_T8='DENIED_NOT_RUN_UNCHANGED',new_maximum_time=.08,
        new_flow_through_fraction=.01,scope='finite predetermined diagnostic design; no global stability or developed calibration claim')
    save(out/'component-gates.json',gates)
    manifest=dict(schema='WJH3CompactPublicationV1',root=str(root),publisher_source_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=source,text=True).strip(),
        publisher_sha256=sha(__file__),native_source='c344dc2209f030e426e2ee5501d5e5b72a6d2c46',
        native_runtime=identity(actual('/media/jay/新加卷/nektar-v2026-runs/wj-h2-20261111/runtime-v4/IncNavierStokesSolver')),
        new_input_generator_source='c1f96524',native_observer_source='e62b5f4d265f068dd5af3f090bfb0e385c3ae53e',
        compact={f.name:identity(f) for f in out.iterdir() if f.is_file()},
        historical_mount_alias=dict(old='/media/jay/新加卷',actual='/media/jay/新加卷1',uuid='5AB8A669B8A642FD',writes_to_originals=False))
    save(out/'manifest.json',manifest);print(gates,flush=True)
if __name__=='__main__':main()
