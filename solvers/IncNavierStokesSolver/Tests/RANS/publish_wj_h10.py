#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H10 final common-owner audit and compact publication; no native invocation."""
import argparse
import csv
from pathlib import Path
import numpy as np
from run_wj_h10 import load,save
from study_wj_h9 import identity,check,serial
from study_wj_h10 import selected_analysis
from wj_h10_analysis import MATRIX,QUALITY
from wj_h7_analysis import Rows,classify
from wj_h8_analysis import extend_witnesses,FACE
from study_wj_h6 import configuration
from wj_h4_points import vector_hash
from wj_h9_local import MomentAudit


def table(path,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('x',newline='')as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)


def common(root):
    out=root/'science/final-common';out.mkdir();reports={};inv=load(root/'science/imports/inventory.json')
    for branch,dt in MATRIX:
        for phase in ('cp0','run'):
            path=root/f'science/cases/{branch}/{phase}/case.evolution.json';r=load(path);assert not r['failure'];reports[branch+'-'+phase]=(r,path)
            inv=extend_witnesses(r,inv,configuration(6.))
    save(out/'inventory.json',inv);base=reports['coarse-cp0'][0];rows=Rows(base,[v['physical']for v in inv['points']],configuration(6.));mom=MomentAudit(rows);cache={};results={}
    for name,(r,path)in reports.items():
        assert r['fields']==base['fields']and r['parameters']==base['parameters'];key=vector_hash(r['state'])
        if key not in cache:
            d=rows.evaluate(r['state']);summary=rows.summary(d,d);ii=np.flatnonzero(rows.outer);j=ii[np.argmin(d['margin'][ii])]
            summary['minimum_including_physical_volume_nodes']=float(d['margin'][j]);summary['volume_inclusive_witness']=rows.witness(d,j)
            summary['volume_inclusive_status']=classify(float(d['margin'][j]),float(np.min(d['transverse'][ii])),float(np.min(d['covariance'])),float(min(np.min(d['k']),np.min(d['epsilon'][ii]))))
            mm=mom.moments(r['state']);summary['moments']={k:dict(values=v['moments'].tolist(),projection_L2=v['projection_L2'])for k,v in mm.items()};cache[key]=summary
        results[name]=dict(producer=identity(path),state_sha256=key,**cache[key])
    ok=all(v['status']=='PASS'and v['volume_inclusive_status']=='PASS'for v in results.values())
    save(out/'result.json',serial(dict(status='PASS'if ok else 'OBSTRUCTED',states=results,unique_states=len(cache),same_owner_rule='union of frozen H8 inventory and all actual H10 native witnesses/neighborhoods; all owners, no adaptive point search',mass_factorizations=mom.projection.factorizations,mass_solves=mom.solves,additional_native=0)))
    if not ok:
        from study_wj_h10 import halt
        halt(root,'OBSTRUCTED','final common inventory quality/coverage failure',witness=serial(results))
    print('H10 final common inventory',ok,'unique states',len(cache),flush=True)


def publish(root):
    out=root/'publication-v1';out.mkdir();book=load(root/'accounting.json');common_result=load(root/'science/final-common/result.json');timeline=[];impulses=[];counts=[];moments={};arithmetic_files=[];weak=[];local=[];norms=[];interface=[];records={};runtime=[]
    for branch,dt in MATRIX:
        for phase in ('cp0','run'):
            v=selected_analysis(root,branch,phase);check(v['producer']);check(v['checkpoint']);r=load(v['producer']['path']);records[branch+'-'+phase]=v
            s=common_result['states'][branch+'-'+phase];strict=root/f'science/views/{branch}-{phase}/result.json';view=load(strict);assert view['status']=='PASS'
            timeline.append(dict(branch=branch,dt=dt,phase=phase,time=r['time'],steps=r['steps'],native='PASS'if not r['failure']else 'OBSTRUCTED',expanded=s['status'],volume_inclusive=s['volume_inclusive_status'],minimum=s['min_in_plane'],minimum_with_volume=s['minimum_including_physical_volume_nodes'],transverse=s['min_transverse'],covariance=s['min_covariance'],raw_divergence=s['raw_divergence_max'],initializer_calls=r['initial_velocity_calls'],initializer_expected=1 if phase=='cp0'else 0,Project_Lift=str(r['initialization_lift_calls_iterations']),state_sha256=vector_hash(r['state']),signature=r['signature']))
            moments[branch+'-'+phase]=dict(native=view['native_Q1_moments'],independent=v['summary']['moments'],final_common=s['moments'])
            norms.append(dict(branch=branch,phase=phase,native_mass_dual=view['norms'],reference_error=v.get('reference_error',v['summary']['reference_error']),initial_drift=v['summary']['initial_drift'],velocity_initialization_change=v.get('velocity_change'),errors=v['summary']['errors']))
            interface.append(dict(branch=branch,phase=phase,**v['summary']['interface_budget']))
            for z in load(Path(v['producer']['path']).parents[3]/'missing.json')if False else load(Path(v['arithmetic']['path']).parent/'local.json')if phase=='run'else load(root/f'science/analysis/{branch}-cp0/local.json'):
                pt=z['physical'];o=z['owner']
                if pt[0]==8. and(abs(pt[2]-FACE)<=1.01e-6 or abs(pt[2]-.2375)<=1.01e-6):local.append(dict(branch=branch,phase=phase,time=r['time'],**z))
            weak.append(dict(branch=branch,phase=phase,result=identity(strict),signed_weak=view['signed_weak'],anchor=view['anchor'],quadrature=view['quadrature'],weak_evaluations=view['weak_evaluations']))
            if phase=='run':
                ar=load(v['arithmetic']['path']);arithmetic_files.append(v['arithmetic']);old=load(root/f'science/imports/{branch}.impulse.json')
                impulses.append(dict(branch=branch,dt=dt,time=r['time'],new_pressure_L2=ar['pressure_physical']['full']['L2'],new_j_L2=ar['impulse_physical']['full']['L2'],old_pressure_L2=old['pressure_increment']['full']['L2'],old_j_L2=old['physical']['full']['L2'],new_pressure_mean=ar['pressure_physical']['full']['mean'],new_j_mean=ar['impulse_physical']['full']['mean'],predictor_L2=ar['predictor_physical_L2'],pressure_correction_L2=ar['pressure_correction_physical_L2'],total_du_L2=ar['total_physical_L2'],predictor_per_dt_L2=ar['predictor_per_dt_physical_L2']))
                save(out/(branch+'-arithmetic.json'),ar)
                a=r['pimple_attempts'][0];counts.append(dict(branch=branch,**{k:a[k]for k in ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls','bound_activations')}))
                for phase_name,vec in [('C',load(a['k_snapshot'])['committed']),('F',load(a['k_snapshot'])['post_flow']),('K',load(a['k_snapshot'])['raw']),('A',r['state'])]:
                    assert len(vec)==2328  # retain aliases/trace roles, never lift here
            label=branch+'-'+('step'if phase=='run'else 'cp0');samples=load(root/'operations'/label/'native-exec-samples.json');assert samples and all(x['status']=='PASS'for x in samples);runtime.append(dict(branch=branch,phase=phase,actual_exec_samples=identity(root/'operations'/label/'native-exec-samples.json')))
    table(out/'timeline.csv',timeline);table(out/'pressure-impulses.csv',impulses);table(out/'solve-counts.csv',counts);save(out/'all-Q1-moments.json',moments);save(out/'norms-errors.json',norms);save(out/'interface-budgets.json',interface);save(out/'local-owned-states.json',serial(local));save(out/'weak-checks.json',weak);save(out/'common-inventory-result.json',common_result)
    loads=[];roles=[]
    for v in weak:
        data=load(v['signed_weak']['path'])['32']
        for kind,d in [('physical',data['regional_components']),('split',data['splits']),('exterior_role',data['exterior_roles'])]:
            for name,rows in d.items():
                for region in range(2):
                    for mode in range(3):loads.append(dict(branch=v['branch'],phase=v['phase'],kind=kind,term=name,region=region,mode=mode,**dict(zip(('u','v','k','outer_epsilon','p_inner','p_outer'),rows[region][mode]))))
    table(out/'signed-weak-components.csv',loads)
    gates=dict(implementation_build='PASS',native_fixture='PASS',startup_compatibility='PASS'if all(v['expanded']=='PASS'for v in timeline if v['phase']=='cp0')else 'OBSTRUCTED',native_hard_admission='PASS'if all(v['native']=='PASS'for v in timeline)else 'OBSTRUCTED',principal_quality=common_result['status'],restart_history='PASS',strict_weak='PASS'if all(v['anchor']['status']==v['quadrature']['status']=='PASS'for v in weak)else 'UNRESOLVED',actual_predictor_correction='PASS',coverage='finite inventory only',H8_common_time_contraction='UNRESOLVED unchanged',new_CP0=book['cp0_reservations'],new_governing=book['governing_attempts'],new_strict=book['strict_readbacks'],closure_transport='PENDING')
    assert(gates['new_CP0'],gates['new_governing'],gates['new_strict'])==(3,3,6)and len(timeline)==6
    gates['status']='PASS'if not book.get('h10_halt')and all(gates[k]=='PASS'for k in ('implementation_build','startup_compatibility','native_hard_admission','principal_quality','strict_weak','actual_predictor_correction'))else 'OBSTRUCTED';save(out/'component-gates.json',gates)
    save(out/'state-bindings.json',{k:dict(producer=v['producer'],checkpoint=v['checkpoint'],invariants=v['invariants'],arithmetic=v.get('arithmetic'))for k,v in records.items()})
    save(out/'runtime.json',dict(build=load(root/'science/runtime-build.json'),qualification=identity(root/'runtime-qualified.json'),actual_cases=runtime))
    save(out/'inputs.json',load(root/'science/input-intent.json'))
    # All figures use retained scalar/profile bytes, no recomputation of a state.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    x=np.array([v['dt']for v in impulses]);fig,axs=plt.subplots(1,3,figsize=(11.5,3.5),constrained_layout=True)
    for name,label in [('pressure_L2','pressure increment'),('j_L2','dt * pressure increment')]:
        ax=axs[0 if name=='pressure_L2'else 1];ax.loglog(x,[v['old_'+name]for v in impulses],'o--',label='H9 incompatible CP0');ax.loglog(x,[v['new_'+name]for v in impulses],'o-',label='H10 compatible CP0');ax.set(xlabel='dt = first-step time (Hh/Uinf)',ylabel=label+' physical L2');ax.legend(fontsize=7)
    for k,label in [('predictor_L2','predictor'),('pressure_correction_L2','pressure correction (added)'),('total_du_L2','total velocity increment')]:axs[2].loglog(x,[v[k]for v in impulses],'o-',label=label)
    axs[2].set(xlabel='dt = first-step time',ylabel='velocity physical L2');axs[2].legend(fontsize=7)
    for ax in axs:ax.grid(alpha=.2)
    fig.suptitle('Different first-step times; no regauging or temporal-order fit');fig.savefig(out/'first-step-impulses.png',dpi=160);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(10,3.6),constrained_layout=True)
    for branch,dt in MATRIX:
        ar=load(out/(branch+'-arithmetic.json'));p=sorted((v for v in ar['pressure_profiles']if v['physical'][0]==8 and v['owner']['region']=='outer'and v['physical'][2]<.7),key=lambda v:v['physical'][2]);axs[0].plot([v['j_dt']for v in p],[40*v['physical'][2]for v in p],'.-',label=f'dt={dt:g}')
        l=[v for v in local if v['branch']==branch and v['phase']=='run'and v['owner']['geometry']==17 and abs(v['physical'][2]-FACE)<1e-14];axs[1].plot(dt,l[0]['normalized_fixed_action'],'o',label=f'dt={dt:g}')
    axs[0].set(xlabel='signed first-step j at x8',ylabel='physical z (mm)');axs[1].axhline(QUALITY,ls='--',color='k',label='half-reference quality');axs[1].set(xlabel='first-step time (Hh/Uinf)',ylabel='fixed-direction action / nu_eff, face17 upper')
    for ax in axs:ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.savefig(out/'owned-profiles-and-action.png',dpi=160);plt.close(fig)
    heavy=[]
    for p in sorted((root/'science').rglob('*')):
        if p.is_file():heavy.append(identity(p))
    save(out/'manifest.json',dict(heavy=heavy,ledger=identity(root/'accounting.json'),source=identity(__file__),scope='full histories, fields, matrices and native view rows stay remote; compact scalars have exact bridges'))
    save(out/'seal.json',dict(files=[identity(p)for p in sorted(out.iterdir())if p.is_file()]))
    print('H10 publication',gates,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('common','publish'),required=True);a=p.parse_args();globals()[a.mode](a.root)
