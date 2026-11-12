#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Small H4 scientific tables/figures and byte/producer bindings; no evolution."""
import argparse
import csv
import gzip
import hashlib
from pathlib import Path
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_wj_h4 import load,save,sha
from wj_h3_analysis import identity,state_context,interface_view
from wj_h4_results import samples
from wj_h4_points import mapping,physical
from wj_h2_profiles import wall,NU,QE


def csvfile(path,rows):
    if not rows:return
    with path.open('x',newline='')as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def publish(root):
    out=root/'publication-v1';out.mkdir();heavy={}
    def index(p):
        p=Path(p)
        if str(p)not in heavy:heavy[str(p)]=identity(p)
        return heavy[str(p)]
    comparisons=load(root/'comparisons-v1/summary.json');minima=load(root/'comparisons-v1/minima.json')
    inv,groups=samples(root);save(out/'inventory.json',inv)
    save(out/'comparison-summary.json',comparisons)
    save(out/'initial-physical-equality.json',load(root/'comparisons-v1/initial-physical-equality.json'))
    csvfile(out/'principal-history.csv',minima)
    runs=[];meshes={};quadrature=[];interfaces=[];stages=[];adverse=[]
    for entry in load(root/'catalog.json'):
        index(entry['path']);v=load(entry['path']);index(v['producer']['path'])
        if v['failures']:adverse.append(dict(case=entry['case'],time=entry['time'],native_status=entry['native_status'],
            expanded_status=entry['expanded_status'],producer=v['producer'],state_sha256=v['state_sha256'],
            callback_history_sha256=v['callback_history_sha256'],inventory_sha256=v['inventory_sha256'],failures=v['failures']))
    save(out/'adverse-witnesses.json',adverse)
    for name,g in groups.items():
        for e in g.values():
            r=load(e['producer']['path']);state,h,alpha,t,packet=state_context(r,e['step'],e['role']=='rejected')
            for x in (0.,8/3,3.,4.,16/3,8.):
                for v in interface_view(r,state,x):
                    interfaces.append(dict(case=name,time=t,step=e['step'],role=e['role'],x=x,
                        inner_geometry=v['geometry_inner_outer'][0],outer_geometry=v['geometry_inner_outer'][1],
                        sampled_K=v['K'],cutoff=v['cutoff'],mass_defect=v['mass_flux_defect'],
                        normal_traction=v['normal_traction_defect'],shear_traction=v['traction_defect'][0],
                        k_flux=v['k_flux_defect'],epsilon_trace=v['epsilon_defect'],pi_inner=v['pi_inner'],pi_outer=v['pi_outer'],
                        mechanical_jump=v['mechanical_jump'],roughness_force_x=v['roughness_force_x'],inner_production_H=v['inner_production']))
            stages.append(dict(case=name,step=e['step'],time=t,role=e['role'],alpha_mass=alpha,
                state_sha256=e['state_sha256'],callback_history_sha256=e['callback_history_sha256'],producer=e['producer'],
                packet_sha256=hashlib.sha256(__import__('json').dumps(packet,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
                startup_inapplicable_completed_stage=e['step']==0))
    csvfile(out/'interface-history.csv',interfaces);save(out/'stage-bindings.json',stages)
    for f in sorted((root/'cases').glob('*/trajectory.json')):
        d=load(f);name=f.parent.name;calls=[]
        for v in d['calls']:
            rp=Path(v['producer']['path']);r=load(rp);index(rp);index(v['xml']['path']);index(v['log']['path'])
            for attempt in r['pimple_attempts']:
                snap=Path(attempt['k_snapshot']);q=load(snap)
                assert np.asarray(q['raw'],dtype='<f8').tobytes()==np.asarray(q['corrected'],dtype='<f8').tobytes()
                index(snap)
            for cp in rp.parent.glob('CP*/State.bin'):index(cp)
            calls.append(dict(producer=v['producer'],xml=v['xml'],seconds=v['seconds'],accounting=v['accounting'],raw_corrected_equal=v['raw_corrected_equal']))
        cp0=f.parent/'cp0/case.evolution.json';index(cp0);ir=load(cp0);design=load(cp0.with_name('case-design.json'))
        index(cp0.with_name('case.xml'));index(Path(ir['checkpoint'])/'State.bin')
        H,d0,top,J=mapping(ir)
        meshes[name]=dict(configuration=d['configuration'],actual_states=len(ir['state']),chart_knots=design['chart_z'],
            physical_knots=[physical(ir,[0,0,z])[2]for z in design['chart_z']],xknots=design['x'],H=H,d0=d0,top=top,J_outer=J,
            generator=identity(cp0.with_name('case-design.json')))
        runs.append(dict(case=name,**{k:d[k]for k in ('configuration','accepted_steps','accepted_time','requested_target','stop','accounting','status')},
            calls=calls,cp0=index(cp0),execution_path=d['execution_path']))
        for role,q in d['outer_checks'].items():
            for pair,v in q['comparisons'].items():quadrature.append(dict(case=name,role=role,time=q['time'],pair=pair,status=v['status'],
                max_component=max(v['errors'].values()),outer_stress=v['errors']['outer_momentum_stress'],regional=v['regional_error']))
    save(out/'runs.json',runs);save(out/'meshes.json',meshes);csvfile(out/'quadrature.csv',quadrature)
    native=load(root/'native-observations-v1/index.json');budgets={};passive={}
    for name,v in native.items():
        index(v['native']['path']);index(v['checkpoint']['path']);index(v['producer']['path']);p=load(v['native']['path']);r=load(v['producer']['path'])
        pref=dest=prod=0.
        for row in p['weak_checks']['samples']:
            if not row['inner']:continue
            K=row['K'];w=wall(K,r['parameters'][1],r['parameters'][2]);nt=w['a']*max(row['chart'][2]-w['zv'],0)
            M=np.array(row['G'])+np.array(row['G']).T
            pref+=row['weight']*(K/QE)**2*nt/(NU+nt)**2
            dest+=row['weight']*w['A']/max(row['chart'][2],w['rc'])
            prod+=row['weight']*(.5*nt*np.sum(M*M)-2*row['k']*np.trace(row['G'])/3)
        ground=[z for z in p['open_boundary']['boundary_rows']if z['role']=='ground']
        budgets[name]=dict(actual_inner_production=prod,counterfactual= pref,destruction=dest,counterfactual_difference=dest-pref,
            mean_wall_stress=-sum(z['weight']*z['actual_traction'][0]for z in ground)/sum(z['weight']for z in ground),
            mass_flux=p['open_boundary']['mass_flux'],exterior_load=p['open_boundary']['boundary_load'],
            native_final_stage_loads=r['final_stage_loads'],native_regional_stage_loads=r['regional_stage_loads'],
            independent_weak=v['independent_weak'])
        passive[name]=dict(producer=v['producer'],checkpoint=v['checkpoint'],native=v['native'],
            callback_history_sha256=v['callback_history_sha256'],physical_mass_dual_norms=v['stage_norms'],
            pressure_projection=v['offline_pressure_projection'],initialization_work=p['initialization_work'],governing_solves=p['weak_checks']['governing_solves'])
    save(out/'weak-budgets.json',budgets);save(out/'passive-readbacks.json',passive)
    # Export a compact owner-labelled subset; full union CSVs remain immutable
    # heavy data. Never omit a point from the scientific assessment itself.
    with gzip.open(out/'paired-declared.csv.gz','wt',newline='')as ff:
        writer=None
        for label,v in comparisons.items():
            if 'full'not in v:continue
            index(v['full']['path']);index(v['limiting']['path'])
            with gzip.open(v['full']['path'],'rt',newline='')as f:
                for row in csv.DictReader(f):
                    x,z=float(row['x']),float(row['z'])
                    if not any(abs(x-t)<1e-12 for t in (3.,4.,7.999999,8.)):continue
                    if not any(abs(z-t)<2e-12 for t in (.5875,.6720111978440669,.8818157002137355,.9080394908593612,1.1947077409162705)):continue
                    row=dict(comparison=label,**row)
                    if writer is None:writer=csv.DictWriter(ff,fieldnames=list(row));writer.writeheader()
                    writer.writerow(row)
    compact=out/'paired-declared.csv.gz';decoded=gzip.decompress(compact.read_bytes())
    save(out/'compressed-identities.json',dict(file=identity(compact),decoded_bytes=len(decoded),decoded_sha256=hashlib.sha256(decoded).hexdigest()))
    # Scientific figures: all minima stay separate from same-position changes.
    plt.rcParams.update({'font.size':8,'figure.dpi':160})
    fig,axes=plt.subplots(1,2,figsize=(9,3.5))
    for name in ('base-x1','base-x1-half','refined-x1','refined-x1-half','base-x3','base-x3-half','local-x3','local-x3-half'):
        rr=sorted((v for v in minima if v['case']==name and v['role']=='accepted'),key=lambda v:v['time'])
        if not rr:continue
        for ax in axes:ax.plot([v['time']for v in rr],[v['margin']for v in rr],'-o',ms=2,label=name)
    axes[0].set_ylim(bottom=0);axes[1].set_ylim(-.0005,.002)
    for ax in axes:ax.axhline(0,color='k',lw=.8);ax.set(xlabel='T=t U∞/Hh',ylabel='expanded finite-inventory principal / νeff')
    axes[0].legend(fontsize=6,ncol=2);axes[1].set_title('Detail near zero; own minima, different locations')
    fig.tight_layout();fig.savefig(out/'principal-history.png');plt.close(fig)
    with gzip.open(compact,'rt')as f:declared=list(csv.DictReader(f))
    fig,axes=plt.subplots(1,2,figsize=(9,3.5))
    for label in ('time-base-x1','time-refined-x1','time-base-x3','normal-dt01','normal-dt005','local-dt005'):
        rr=sorted((v for v in declared if v['comparison']==label and abs(float(v['x'])-3)<1e-12 and abs(float(v['z'])-1.1947077409162705)<2e-12 and v['normal_side_a']!='lower' and v['normal_side_b']!='lower'),key=lambda v:float(v['intended_time']))
        if not rr:continue
        axes[0].plot([float(v['intended_time'])for v in rr],[float(v['delta_principal'])for v in rr],'-o',ms=2,label=label)
        axes[1].plot([float(v['intended_time'])for v in rr],[float(v['delta_Gxz'])for v in rr],'-o',ms=2,label=label)
    for ax in axes:ax.set_xlabel('T=t U∞/Hh');ax.axhline(0,color='k',lw=.6)
    axes[0].set_ylabel('signed Δ(principal / νeff)');axes[1].set_ylabel('signed ΔGxz Hh/U∞')
    axes[0].legend(fontsize=6);fig.suptitle('Same point x=3, z=1.19470774; below-face/interior owners')
    fig.tight_layout();fig.savefig(out/'same-point-sensitivity.png');plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(9,3.5))
    for name in ('base-x1-half','refined-x1-half','base-x3-half','local-x3-half'):
        rr=[v for v in interfaces if v['case']==name and abs(v['x']-8)<1e-12]
        axes[0].plot([v['time']for v in rr],[v['normal_traction']for v in rr],'-o',ms=2,label=name)
        axes[1].plot([v['time']for v in rr],[v['shear_traction']for v in rr],'-o',ms=2,label=name)
    axes[0].set_ylabel('H normal traction defect / U∞²');axes[1].set_ylabel('H shear traction defect / U∞²')
    for ax in axes:ax.set_xlabel('T=t U∞/Hh');ax.axhline(0,color='k',lw=.6)
    axes[0].legend(fontsize=6);fig.suptitle('Outlet/interface corner x=8; pointwise defects, no added strong constraint')
    fig.tight_layout();fig.savefig(out/'interface-defects.png');plt.close(fig)
    counts=[sum(v['accounting']['counts'][i]for v in runs)for i in range(5)]
    required=[v for k,v in comparisons.items()if k!='retained-A-time']
    core_names={'base-x1-half','refined-x1','refined-x1-half','base-x3-half'}
    core={v['case']:v['status']for v in runs if v['case']in core_names}
    expected_readbacks={v['case']for v in runs}|{'base-x1-half-cp0','local-x3-cp0'}
    gates=dict(schema='WJH4ComponentGatesV1',core='PASS'if set(core)==core_names and all(v=='PASS_SHORT_CONTROL'for v in core.values())else 'INCOMPLETE',
        fresh_cases={v['case']:v['status']for v in runs},expanded_native_accepted='PASS'if all(v['expanded_status']=='PASS'for v in minima if v['role']=='accepted')else 'OBSTRUCTED',
        common_points='PASS'if all(v['status']=='COMPLETE_DATA'for v in required)else 'UNRESOLVED',
        evolved_quadrature='PASS'if all(v['status']=='PASS'for v in quadrature)else 'UNRESOLVED',
        independent_weak='PASS'if set(native)==expected_readbacks and all(v['independent_weak']['status']=='PASS'for v in native.values())else 'UNRESOLVED',
        attempts=sum(v['accounting']['attempts']for v in runs),solve_bound_counts=counts,
        numerical_corrections=sum(v['accounting']['corrections']for v in runs),
        execution_path_check='NOT_RUN:no inventory gap',max_time=max(v['accepted_time']for v in runs),
        H2_T8='DENIED_NOT_RUN_UNCHANGED',later_T04='NOT_RUN_OUT_OF_SCOPE',developed_precursor='NOT_QUALIFIED',
        scope='finite physical union, all owners and local angular action; sensitivity, not asymptotic order or global stability')
    gates['matched_local_admissibility']='OBSTRUCTED'if any(v['status']=='OBSTRUCTED'for v in runs)else 'PASS_FINITE_INVENTORY'if all(v['status']=='PASS_SHORT_CONTROL'for v in runs)and len(runs)==6 else 'INCOMPLETE'
    save(out/'component-gates.json',gates)
    commands=[]
    for p in sorted((root/'operations').glob('*/receipt.json')):
        commands.append(load(p));index(p);index(p.with_name('request.json'));index(p.with_name('command.log'))
    save(out/'commands.json',commands)
    for p in (root/'receipts').glob('*.json'):index(p)
    for p in root.glob('admission-*.json'):index(p)
    for p in root.glob('optional-activation.json'):index(p)
    for p in (root/'cases').glob('*/completed-cp0-recovery.json'):index(p)
    for p in (root/'scan-extensions').glob('*/completed-passive-recovery.json'):index(p)
    index(root/'focused-config/Testing/Temporary/LastTest.log')
    save(out/'heavy-index.json',heavy)
    source=Path(__file__).resolve().parents[4]
    manifest=dict(schema='WJH4PublicationV1',root=str(root),publisher_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=source,text=True).strip(),
        publisher=identity(__file__),native_source='c344dc2209f030e426e2ee5501d5e5b72a6d2c46',observer_native_source='e62b5f4d265f068dd5af3f090bfb0e385c3ae53e',
        producing_python=['a9d9023','bd72c24'],runtime_equivalence=load(root/'receipts/runtime-equivalence.json'),
        compact={p.name:identity(p)for p in out.iterdir()if p.is_file()},historical_paths_translated_read_only=True)
    save(out/'manifest.json',manifest);print(gates,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();publish(a.root)
