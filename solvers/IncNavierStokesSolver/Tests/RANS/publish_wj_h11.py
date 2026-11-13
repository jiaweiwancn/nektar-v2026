#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H11 compact scalar publication; full rows/matrices/checkpoints stay remote."""
import argparse
import csv
from pathlib import Path
import re
import numpy as np
from run_wj_h11 import H10,load,save
from study_wj_h9 import identity,check
from wj_h11_analysis import MATRIX,VIEWS
from study_wj_h11 import parent,audit


def table(path,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('x',newline='')as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)


def publish(root):
    out=root/'publication-v1';out.mkdir();book=load(root/'accounting.json');save(out/'accounting-at-publication.json',book)
    common=load(root/'science/common-v1/result.json');records=load(root/'science/common-v1/states.json');timeline=[];norms=[];interfaces=[];pressure=[];counts=[];views=[];loads=[];weak=[];sources=[]
    save(out/'common-union.json',common);save(out/'principal-pairs.json',load(root/'science/common-v1/principal-pairs.json'));table(out/'field-pairs.csv',load(root/'science/common-v1/field-pairs.json'))
    for v in records:
        s=v['summary'];w=s['witness'];p=w['physical'];o=w['owner'];timeline.append(dict(branch=v['branch'],step=v['step'],recorded_time=v['time'],comparison_time=v['comparison_time'],alpha=v['alpha'],
            native='ACCEPTED',origin='H10 imported'if v['step']<=1 else 'H11 new',expanded=v['status'],minimum=s['min_in_plane'],minimum_volume=s['minimum_volume_inclusive'],transverse=s['min_transverse'],covariance=s['min_covariance'],min_scalar=s['min_scalar'],strong_divergence=s['raw_divergence_max'],
            x=p[0],physical_z=p[2],physical_z_mm=40*p[2],geometry=o['geometry'],normal_side=o['normal_side'],x_side=o['x_side'],producer_sha256=v['producer']['sha256'],state_sha256=v['state_sha256'],history_sha256=v['history_sha256'],compiled_action=v['compiled_replay']))
        norms.append(dict(branch=v['branch'],step=v['step'],time=v['time'],errors=s['errors'],reference_error=s['reference_error'],compatible_CP0_drift=s['initial_drift']))
        interfaces.append(dict(branch=v['branch'],step=v['step'],time=v['time'],**s['interface_budget']))
        pressure.append(dict(branch=v['branch'],step=v['step'],time=v['time'],regions={k:{kk:vv for kk,vv in x.items()if kk!='elements'}for k,x in v['projection']['regions'].items()}))
    table(out/'timeline.csv',timeline);save(out/'errors-drift.json',norms);save(out/'interface-history.json',interfaces);save(out/'Q1-continuity.json',pressure)
    fixed=[]
    for v in load(root/'science/common-v1/fixed-owner-history.json'):
        d=v['record'];p=d['physical'];o=d['owner'];G=np.array(d['G']);R=np.array(d['R'])
        fixed.append(dict(branch=v['branch'],step=v['step'],time=v['time'],x=p[0],physical_z=p[2],physical_z_mm=40*p[2],geometry=o['geometry'],normal_side=o['normal_side'],x_side=o['x_side'],
            U=d['U'][0],W=d['U'][2],k=d['k'],epsilon=d['epsilon'],k_over_epsilon=d['k']/d['epsilon'],pi=d['pi'],mechanical_pressure=d['mechanical_pressure'],K=d['K'],rc=d['rc'],root=d['root'],nu_eff=d['nu_eff'],nu_d=d['nu_d'],margin=d['margin'],transverse=d['transverse'],minimizing_angle=d['angle'],raw_production=d['production'],
            **{f'G{i}{j}':G[i,j]for i in range(3)for j in range(3)},**{f'R{i}{j}':R[i,j]for i in range(3)for j in range(3)}))
    table(out/'fixed-owner-history.csv',fixed)
    for branch,dt,n,target in MATRIX:
        case=parent(root,branch);a=audit(root,branch);check(a['producer']);r=load(case['report']);old=load(case['parent']['producer']['path'])
        assert r['geometry_signature']==old['geometry_signature']and r['map_parameters']==old['map_parameters']
        counts.append(dict(branch=branch,reserved=n,attempted=a['attempted'],accepted=a['accepted'],qualified=a['qualified'],absolute_target=target,actual_absolute_step=r['steps'],time=r['time'],failure=r['failure'],initializer_calls=r['initial_velocity_calls'],initialization_lift=r['initialization_lift_calls_iterations'],
            **{k:sum(x[k]for x in r['pimple_attempts'])for k in('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls','bound_activations','additional_governing_equation_solves','additional_correction_or_restoration_solves')}))
        sources.append(dict(branch=branch,xml=case['xml'],parent=case['parent'],producer=a['producer'],checkpoint_index=a['index'],phases=a['phases'],every_step=a['states']))
        for step in VIEWS[branch]:
            folder=root/f'science/views/{branch}-{step}';v=load(folder/'result.json');assert v['status']=='PASS';views.append(dict(branch=branch,step=step,**v))
    table(out/'attempts-solves.csv',counts);save(out/'source-state-bindings.json',sources)
    old_views=load(root/'source/docs/plans/wj-h10-linux/publication-v1/weak-checks.json');old=next(v for v in old_views if v['branch']=='coarse'and v['phase']=='run');check(old['result']);reused=load(old['result']['path']);views.insert(0,dict(branch='coarse',step=1,reused_H10=True,**reused))
    for v in views:
        cx=v['context'];weak.append(dict(branch=v['branch'],step=v['step'],reused_H10=v.get('reused_H10',False),context=cx,native=v['native'],signed_weak=v['signed_weak'],anchor=v['anchor'],quadrature=v['quadrature'],norms=v['norms'],open_boundary=v['open_boundary']))
        data=load(v['signed_weak']['path'])['32']
        for kind,components in [('physical',data['regional_components']),('split',data['splits']),('exterior_role',data['exterior_roles'])]:
            for name,rows in components.items():
                for region in range(2):
                    for mode in range(3):loads.append(dict(branch=v['branch'],step=v['step'],time=cx['time'],kind=kind,term=name,region=region,mode=mode,**dict(zip(('u','w','k','outer_epsilon','p_inner','p_outer'),rows[region][mode]))))
    save(out/'strict-weak-stage.json',weak);table(out/'signed-weak-components.csv',loads)
    ar=load(root/'science/arithmetic-v1/result.json');save(out/'arithmetic-bindings.json',ar);imp=[]
    for branch,dt,n,target in MATRIX:
        packet=load(ar['branches'][branch]['packet']['path']);save(out/(branch+'-arithmetic.json'),packet)
        imp.append(dict(branch=branch,dt=dt,absolute_step=packet['absolute_step'],time=packet['time'],alpha=packet['alpha'],pressure_L2=packet['pressure_physical']['full']['L2'],j_L2=packet['impulse_physical']['full']['L2'],pressure_mean=packet['pressure_physical']['full']['mean'],j_mean=packet['impulse_physical']['full']['mean'],predictor_L2=packet['predictor_physical_L2'],correction_L2=packet['pressure_correction_physical_L2'],total_du_L2=packet['total_physical_L2']))
    table(out/'common-time-increments.csv',imp)
    tests=set()
    for label in ('focused-tests','focused-h11-final','focused-h11-callback'):
        text=(root/'operations'/label/'command.log').read_text();assert '\nOK\n'in text;found=re.findall(r'^(test_\w+) \(([^)]+)\)',text,re.M);assert len(found)==int(re.search(r'Ran (\d+) tests',text).group(1));tests.update(found)
    assert len(tests)==75
    native=all(not c['failure']and c['accepted']==c['reserved']and c['qualified']==c['reserved']for c in counts)
    strict=len(views)==12 and all(v['anchor']['status']==v['quadrature']['status']=='PASS'for v in views)
    validation=native and strict and common['status']=='PASS'and common['complete']and ar['status']=='PASS'and len(records)==31 and not book.get('h11_halt')and not book.get('coverage_halt')
    gates=dict(verification='PASS'if validation else 'OBSTRUCTED',native_completion='PASS'if native else 'OBSTRUCTED',hard_quality_coverage=common['status'],strict_weak='PASS'if strict else 'UNRESOLVED',selected_arithmetic=ar['status'],
        temporal_common_support=common['temporal'],status='PASS'if validation and common['temporal']=='CONTRACTS'else 'UNRESOLVED'if validation else 'OBSTRUCTED',
        new_governing_reserved=book['reserved_attempts'],new_governing_attempted=sum(c['attempted']for c in counts),new_governing_accepted=sum(c['accepted']for c in counts),native_restart_blocks=book['governing_blocks'],new_CP0=0,new_builds=0,new_strict=book['strict_readbacks'],unique_passing_methods=len(tests),closure_transport='PENDING',
        scope='bounded compatible-initialization control toT=.04; no temporal order/global stability/physical inlet claim')
    assert(gates['new_governing_reserved'],gates['new_governing_attempted'],gates['new_governing_accepted'],gates['native_restart_blocks'],gates['new_strict'])==(25,25,25,3,11)
    save(out/'component-gates.json',gates);save(out/'methods.json',dict(methods=[list(v)for v in sorted(tests)],unique=len(tests),CTest_execution='none; direct Python invocation only'))
    save(out/'runtime.json',dict(qualified_H10=load(root/'science/input-intent.json')['runtime'],source_import=identity(root/'science/import-compatibility-v1.json'),original_compiled_inputs=identity(H10/'closure/compiled-inputs.json'),source_capsule=load(root/'source-archive.json')))
    # Figures use scalar results only, without new field evaluation.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    pairs=load(out/'principal-pairs.json');fig,axs=plt.subplots(1,2,figsize=(10,3.7),constrained_layout=True)
    for geometry,side,label in [(17,'upper','10.243mm:17 upper'),(18,'lower','10.243mm:18 lower'),(18,'upper','12.779mm:18 upper'),(19,'lower','12.779mm:19 lower')]:
        p=[v for v in pairs if v['time']>0 and v.get('owner',{}).get('geometry')==geometry and v['owner']['normal_side']==side];axs[0].plot([v['time']for v in p],[v.get('rho',np.nan)for v in p],'o-',label=label)
    p=[v for v in pairs if v['time']>0 and 'scope'in v];axs[1].plot([v['time']for v in p],[v['D1']for v in p],'o-',label='D1:.01/.005');axs[1].plot([v['time']for v in p],[v['D2']for v in p],'o-',label='D2:.005/.0025')
    axs[0].axhline(1,color='k',ls='--');axs[0].set(xlabel='common time Hh/Uinf',ylabel='same-owner principal D2/D1');axs[1].set(xlabel='common time Hh/Uinf',ylabel='maximum same-owner outer principal difference')
    for ax in axs:ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.savefig(out/'same-owner-contraction.png',dpi=160);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(10,3.7),constrained_layout=True)
    for branch,dt,n,target in MATRIX:
        t=[v for v in timeline if v['branch']==branch];axs[0].plot([v['recorded_time']for v in t],[v['minimum']for v in t],'o-',label=f'dt={dt:g}')
        q=[v for v in norms if v['branch']==branch];axs[1].plot([v['time']for v in q],[v['reference_error']['full']['U']['L2']for v in q],'-',label=f'reference dt={dt:g}');axs[1].plot([v['time']for v in q],[v['compatible_CP0_drift']['full']['U']['L2']for v in q],'--',label=f'drift dt={dt:g}')
    axs[0].axhline(.002068763923454,color='k',ls='--',label='half-reference quality');axs[0].set(xlabel='actual time Hh/Uinf',ylabel='minimum normalized principal action');axs[1].set(xlabel='actual time Hh/Uinf',ylabel='velocity full physical L2')
    for ax in axs:ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.savefig(out/'margin-error-drift.png',dpi=160);plt.close(fig)
    heavy=[identity(p)for p in sorted((root/'science').rglob('*'))if p.is_file()];save(out/'manifest.json',dict(heavy=heavy,accounting_snapshot=identity(out/'accounting-at-publication.json'),source=identity(__file__),scope='remote heavy hashes; compact transport independently checked on Mac'))
    save(out/'seal.json',dict(files=[identity(p)for p in sorted(out.iterdir())if p.is_file()],status=gates['status']));print('H11 publication',gates,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();publish(a.root)
