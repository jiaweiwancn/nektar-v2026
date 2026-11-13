#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact signed H12 publication and verification; no new field evaluation."""
import argparse
import csv
import gzip
import hashlib
import re
from pathlib import Path
import numpy as np
from run_wj_h4 import load,save,sha
from study_wj_h12 import serial,file_id,check,H11
from wj_h12_phase import identity,scale,comparison_scale
from wj_h11_analysis import ratio

def table(path,rows):
    keys=list(dict.fromkeys(k for v in rows for k in v))
    with Path(path).open('x',newline='')as f:
        w=csv.DictWriter(f,keys);w.writeheader();w.writerows(rows)

def methods(root):
    ids=set()
    for p in(root/'operations').glob('*/command.log'):
        for name,cls in re.findall(r'^(test_\w+) \(([^)]+)\) \.\.\. ok$',p.read_text(),re.M):
            # Python's unittest includes the method in parentheses on some
            # versions; normalize without confusing invocation/unique counts.
            ids.add(cls if cls.endswith('.'+name)else cls+'.'+name)
    return sorted(ids)

def publish(root):
    r=root/'science/results-v2';out=root/'publication-v1';out.mkdir(exist_ok=False)
    result=load(r/'result.json');assert result['status']=='PASS'
    roles=load(r/'local-phase-roles.json');cross=load(r/'cross-window-accounting.json');bindings=load(r/'phase-bindings.json')
    transitions=load(r/'phase-decompositions.json');imp=load(root/'science/input-intent.json')
    # This immutable snapshot deliberately retains the active publication
    # reservation. It is not relabelled as final closure accounting.
    save(out/'accounting-at-publication.json',load(root/'accounting.json'))
    save(out/'source-input-bindings.json',imp)
    save(out/'phase-bindings.json',bindings);save(out/'endpoint-replay.json',load(r/'endpoint-replay.json'))
    save(out/'compiled-checks.json',load(r/'compiled-checks.json'));save(out/'kernel-native-exec.json',load(r/'kernel-native-exec.json'))
    save(out/'predictor-coverage.json',load(r/'predictor-coverage.json'))
    save(out/'operations-summary.json',[dict(operation=p.parent.name,**load(p))for p in sorted((root/'operations').glob('*/receipt.json'))])
    save(out/'passive-repair.json',load(root/'science/passive-repair-v3.json'))
    save(out/'retained-pilot.json',load(root/'science/results-v1/retained-pilot.json'))
    save(out/'methods.json',dict(unique_passing_methods=len(methods(root)),methods=methods(root),
        excluded='H9 all-Q1/mass and whole-volume tests; no CTest/configure/build invocation'))
    compression=[]
    for name in ('local-phase-roles.json','phase-decompositions.json','window-cumulative.json','cross-window-accounting.json','endpoint-decompositions.json'):
        original=r/name;target=out/(name+'.gz');raw=original.read_bytes();target.write_bytes(gzip.compress(raw,mtime=0))
        compression.append(dict(transport=file_id(target),decoded=file_id(original),encoding='gzip; original JSON bytes, including signed zeros'))
    save(out/'compressed-identities.json',compression)
    rows=[]
    for v in roles:
        x=v['value'];o=v['owner'];row={k:v[k]for k in ('branch','step','window','phase','time','attempt_time','state_absolute_step','accepted_role','state_sha256','alias_key','owner_index')}
        row.update(geometry=o['geometry'],normal_side=o['normal_side'],x_side=o['x_side'],x=v['physical'][0],physical_z=v['physical'][2],physical_z_mm=40*v['physical'][2],constitutive_status=x['constitutive_status'],physical_status=x.get('physical_status'))
        for k in ('k','epsilon','tau','pi','mechanical_pressure','trace_G','N','F_N','denominator_Q','nu_eff','nu_d','Praw','d','a','g','mhat','transverse','minimizing_angle','root_residual','implicit_dN','implicit_derivative_error'):
            row[k]=x.get(k)
        if x['constitutive_status']=='DEFINED':row['covariance_min']=min(x['covariance'])
        for k in ('G','B','R'):
            for i,z in enumerate(np.asarray(x[k]).ravel()):row[k+str(i//3)+str(i%3)]=z
        for k,z in zip(('U','V','W'),x['U']):row[k]=z
        rows.append(row)
    table(out/'signed-local-phases.csv',rows)
    rows=[]
    for v in cross:
        if v['status']!='PASS':continue
        shape=np.shape(v['endpoint'])
        for index in np.ndindex(shape):
            val=lambda a:float(np.asarray(a)[index]);p=v['contributions'];can=v['cancellation']
            row=dict(window=v['window'],pair=v['pair'],owner_index=v['owner_index'],observable=v['observable'],component=','.join(map(str,index))if index else 'scalar',
                     incoming=val(v['incoming']),lift=val(p['lift']),flow=val(p['flow']),k=val(p['k']),epsilon=val(p['epsilon']),endpoint=val(v['endpoint']),
                     telescope_error=val(v['identity']['residual']),roundoff_bound=val(v['identity']['bound']),S=val(can['S']),cancellation_floor=can['arithmetic_floor'])
            cv=np.asarray(can['cancellation'],dtype=object)[index];row['cancellation']=cv;row['ratio_status']='RESOLVED'if np.asarray(can['resolved'])[index]else 'UNRESOLVED'
            rows.append(row)
    table(out/'signed-window-contributions.csv',rows)
    # Reproduce actual published CSV values as well as the hash-bound H11 JSON.
    old=H11/'publication-v1/fixed-owner-history.csv';frozen=root/'source/docs/plans/wj-h11-linux/publication-v1/fixed-owner-history.csv'
    assert sha(old)==sha(frozen)
    with old.open()as f:published=list(csv.DictReader(f))
    accepted=[v for v in roles if v['phase']=='A'and abs(v['time']-([.01,.04][v['window']]))<1e-15]
    verification=[]
    for v in accepted:
        candidates=[z for z in published if z['branch']==v['branch']and int(z['step'])==v['step']and int(z['geometry'])==v['owner']['geometry']and z['normal_side']==v['owner']['normal_side']]
        assert len(candidates)==1;z=candidates[0];error=abs(float(z['margin'])-v['value']['mhat']);assert error<2e-11
        verification.append(dict(branch=v['branch'],step=v['step'],owner_index=v['owner_index'],margin_error=error))
    assert len(verification)==24
    save(out/'published-csv-replay.json',dict(status='PASS',parent=file_id(old),records=verification))
    pair_rows=[]
    for w in range(2):
        for i in range(4):
            a={v['observable']:v for v in cross if v['window']==w and v['owner_index']==i and v['pair']=='coarse-fine'}
            b={v['observable']:v for v in cross if v['window']==w and v['owner_index']==i and v['pair']=='fine-quarter'}
            for name in ('a','mhat','g','d','nu_eff','tau'):
                v=ratio(abs(a[name]['endpoint']),abs(b[name]['endpoint']),4e-11*comparison_scale(name))
                pair_rows.append(dict(window=w,time=[.01,.04][w],owner_index=i,observable=name,coarse_minus_fine=a[name]['endpoint'],fine_minus_quarter=b[name]['endpoint'],**v))
    save(out/'endpoint-pairs.json',pair_rows)
    # Accumulate the exact dimensional/normalization/gap identities, with
    # every saved attempt and the late incoming difference retained.
    normalized=[]
    for w in range(2):
        branch={}
        for b in ('coarse','fine','quarter'):
            steps={v['step']for v in roles if v['window']==w and v['branch']==b}
            for i in range(4):
                q=[v['decomposition']for v in transitions if v['branch']==b and v['step']in steps and v['owner_index']==i]
                branch[b,i]={k:sum(z[k]for z in q)for k in ('dimensional_term','normalization_term','delta_g','delta_mhat')}
        for a,b in (('coarse','fine'),('fine','quarter')):
            for i in range(4):
                prior=next(v for v in cross if v['window']==w and v['pair']==a+'-'+b and v['owner_index']==i and v['observable']=='mhat')
                d={k:branch[a,i][k]-branch[b,i][k]for k in branch[a,i]};start=prior['incoming'];end=prior['endpoint']
                error=end-start-d['dimensional_term']-d['normalization_term']+d['delta_g']
                check_identity=identity(error,[end,start,d['dimensional_term'],d['normalization_term'],d['delta_g']])
                normalized.append(dict(window=w,pair=a+'-'+b,owner_index=i,incoming=start,endpoint=end,
                    dimensional_term=d['dimensional_term'],normalization_term=d['normalization_term'],negative_gap_term=-d['delta_g'],identity=check_identity))
    save(out/'window-normalization-gap.json',serial(normalized))
    gates=dict(diagnostic_completion='PASS',input_roles='PASS',endpoint_reproduction='PASS',signed_algebra='PASS',
        compiled_replay='PASS',bounded_predictor_coverage='PASS: 3/7 early, 0/7 late, unequal incoming durations',
        H11_historical_temporal='UNRESOLVED',new_governing_attempts=0,new_CP0=0,new_strict_readbacks=0,new_builds=0,
        diagnostic_mass_factorizations=0,diagnostic_mass_solves=0,saved_governing_matrix_matvecs=0,
        saved_governing_matrix_factorizations=0,saved_governing_matrix_solves=0,method_count=len(methods(root)),
        attempts=14,raw_phase_owner_roles=224,predecessor_owner_roles=56,kernel_requests=result['kernel_requests'],
        kernel_unique_inputs=result['kernel_unique_inputs'],unique_state_evaluations=result['unique_physical_state_evaluations'],
        physical_warnings={phase:sum(v['value']['physical_status']!='SAMPLED_POSITIVE'for v in roles if v['phase']==phase)for phase in ('P','C','F','K','A')},
        retained_failed_passive_operation='sparse-science-v2; endpoint schema KeyError; all phase calculations retained and reused',
        interpretation='telescopes locate response in the saved sequence, not a unique discretization cause or temporal-order proof')
    save(out/'component-gates.json',gates)
    # Two signed figures on the retained physical owner, never magnitudes only.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(9,5.8),sharey=True)
    colors=['.5','.8','#4477aa','#228833','#ee7733','black']
    for w in range(2):
        for j,pair in enumerate(('coarse-fine','fine-quarter')):
            v=next(z for z in cross if z['window']==w and z['pair']==pair and z['owner_index']==0 and z['observable']=='mhat')
            values=[v['incoming'],*[v['contributions'][k]for k in ('lift','flow','k','epsilon')],v['endpoint']]
            ax=axes[w,j];ax.bar(range(6),np.array(values)*1e6,color=colors);ax.axhline(0,color='.4',lw=.6)
            ax.set_xticks(range(6),['incoming','lift','flow','k','epsilon','end'],rotation=25);ax.set_title(pair+'; '+('[0,.01]'if w==0 else '[.03,.04]'))
            ax.set_ylabel('signed normalized principal difference × 10⁶')
    fig.suptitle('Geometry 17 upper, x=8, physical z=10.243007 mm\ncoarser minus finer; incoming history stays separate')
    fig.tight_layout();fig.savefig(out/'signed-window-contributions.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(9,3.8),sharey=True)
    for w,ax in enumerate(axes):
        packet=load(r/'endpoint-decompositions.json');items=[next(v['decomposition']for v in packet if v['window']==w and v['owner_index']==0 and v['pair']==pair)for pair in ('coarse-fine','fine-quarter')]
        x=np.arange(4);width=.36
        for j,(v,pair)in enumerate(zip(items,('coarse-fine','fine-quarter'))):
            ax.bar(x+(j-.5)*width,np.array([v['dimensional_term'],v['normalization_term'],-v['delta_g'],v['delta_mhat']])*1e6,width,label=pair)
        ax.set_xticks(x,['dimensional\naction term','normalization','negative gap','minimized\ndifference']);ax.axhline(0,color='.4',lw=.6)
        ax.set_title('same-owner endpoints at t='+('.01'if w==0 else '.04'));ax.legend(fontsize=8)
    axes[0].set_ylabel('signed normalized difference × 10⁶')
    fig.suptitle('Fixed direction versus minimized action: geometry 17 upper\nphysical z=10.243007 mm; dimensional term is scaled by mean(1/ν_eff)')
    fig.tight_layout();fig.savefig(out/'fixed-normalization-gap.png',dpi=170);plt.close(fig)
    files=[file_id(p)for p in sorted(out.iterdir())if p.is_file()]
    save(out/'manifest.json',dict(artifacts=files,heavy_results=[file_id(p)for p in sorted(r.iterdir())if p.is_file()],
        source=file_id(__file__),parent_intent=file_id(root/'science/input-intent.json'),scope='published compact bytes versus in-place original coefficients; publication snapshot is not final accounting'))
    save(out/'seal.json',dict(manifest=file_id(out/'manifest.json'),component_gates=file_id(out/'component-gates.json')))
    print('H12 publication PASS',gates,flush=True)
    for v in cross:
        if v['owner_index']==0 and v['observable']in ('mhat','a','G','B','tau'):
            print('WINDOW',v['window'],v['pair'],v['observable'],'incoming',v['incoming'],'contributions',v['contributions'],'end',v['endpoint'],flush=True)

def verify(root):
    out=root/'publication-v1';r=root/'science/results-v2';result=load(r/'result.json');manifest=load(out/'manifest.json')
    for v in manifest['artifacts']+manifest['heavy_results']:check(v)
    decoded=[]
    for v in load(out/'compressed-identities.json'):
        raw=gzip.decompress(Path(v['transport']['path']).read_bytes());assert len(raw)==v['decoded']['bytes']and hashlib.sha256(raw).hexdigest()==v['decoded']['sha256']
        decoded.append(v)
    plan=load(root/'science/input-intent.json')
    for rr in plan['reports'].values():
        for v in rr.values():check(v)
    for rec in plan['attempts']:check(rec['packet'])
    for v in plan['arithmetic'].values():check(v)
    assert len(plan['attempts'])==14 and len(load(r/'local-phase-roles.json'))==280
    assert result['new_governing_attempts']==result['new_CP0']==result['new_strict_readbacks']==result['new_builds']==0
    assert result['kernel_requests']<=48 and len(methods(root))==24
    assert all(v['P_C_bitwise']and v['raw_corrected_bitwise']for v in load(r/'phase-bindings.json'))
    closure=root/'closure-v1';closure.mkdir(exist_ok=False)
    save(closure/'verification.json',dict(status='PASS',verified_publication_manifest=file_id(out/'manifest.json'),
        verified_parent_intent=file_id(root/'science/input-intent.json'),decoded_identities=decoded,
        exact_saved_attempts=14,raw_C_F_K_A_owner_roles=224,predecessor_roles=56,unique_passing_methods=24,
        new_governing_attempts=0,new_CP0=0,new_strict_readbacks=0,new_builds=0,kernel_requests=result['kernel_requests'],
        H11_temporal='UNRESOLVED',scope='fresh Linux selected-heavy hash verification; independent Mac compact transport verification follows'))
    print('H12 compact/decoded/selected-heavy verification PASS',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('publish','verify'),required=True);a=p.parse_args()
    globals()[a.mode](a.root)
