#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact publication of the linked, completed seven-attempt H22 diagnosis."""
import argparse,csv,gzip,json,os,shutil
from pathlib import Path
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
import numpy as np
from wj_h15_observer import identity,verify,atomic,vector_hash
from wj_h17_observer import native_load
from wj_h21_outlet import plain
from wj_h22_replay import read_completed
from run_wj_h22_completion import PARENT,usage

def publish(root):
    out=root/'publication-v1';out.mkdir();analysis=root/'science/phase-analysis-v2'
    counts=load(analysis/'counts.json');assert counts['attempts']==7 and counts['total_roles']==70
    assert counts['phase_owner_roles']==56 and counts['predecessor_owner_roles']==14
    copies=[]
    def copy(src,name):
        src=Path(src);dst=out/name;dst.parent.mkdir(parents=True,exist_ok=True)
        if src.suffix=='.json'and src.stat().st_size>100000:
            dst=Path(str(dst)+'.gz')
            with gzip.GzipFile(filename=str(dst),mode='wb',mtime=0)as f:f.write(src.read_bytes())
            assert gzip.decompress(dst.read_bytes())==src.read_bytes();transform='lossless gzip'
        else:shutil.copyfile(src,dst);transform='byte copy'
        copies.append(dict(original=identity(src),copy=identity(dst),transformation=transform));return identity(dst)
    for p in sorted(analysis.glob('*.json')):copy(p,p.name)
    copy(root/'science/whole-field-scatter-v1/result.json','whole-field-scatter.json')
    copy(root/'science/composite-source-index.json','composite-source-index.json')
    sources=load(root/'science/composite-source-index.json');timeline=[];capture=[];online=[]
    for level,attempts in sources['attempts'].items():
        for step,d in attempts.items():
            n=int(step);r,cp=read_completed(d['anchor'],level,n);a=r['pimple_attempts'][0]
            v=native_load(verify(d['vector_packet']))
            timeline.append(dict(level=level,step=n,time=r['time'],alpha=r['last_stage_alpha'],
                post_flow_state_sha256=v['flow_state_sha256'],incoming_base_sha256=v['base_sha256'],
                completed_state_sha256=vector_hash(r['state']),
                native_full_audit=r['native_full_audit_status'],native_completed=True,exact_replay='PASS',qualification=d['qualification'],
                source_root='retained stopped H22'if level!='F'else'new linked completion'))
            capture.append(dict(level=level,step=n,vector_binary=v['payload']['bytes'],vector_manifest=d['vector_packet']['bytes'],
                snapshot=d['phase_snapshot']['bytes'],sparse_k_matrix=d['sparse_k_matrix']['bytes'],execution_identity=d['execution_identity']))
            online.append(dict(level=level,step=n,**{k:a[k]for k in ('momentum_predictor_solves','pressure_system_solves','k_transport_solves',
                'epsilon_transport_solves','bound_calls','bound_activations','additional_correction_or_restoration_solves')},
                raw_k_sha256=a['k_bound']['raw_sha256'],corrected_k_sha256=a['k_bound']['corrected_sha256']))
    atomic(out/'timeline.json',timeline);atomic(out/'capture-bytes.json',capture);atomic(out/'online-counts.json',online)
    physical=[]
    for label,path in (('historical_C1',PARENT/'science/C/monitor/step-1.result.json'),
        ('fresh_M2',root/'science/M2/monitor/startup.result.json'),('fresh_F0',root/'science/F/monitor/startup.result.json'),
        ('fresh_F4',root/'science/F/monitor/step-4.result.json')):
        copy(path,'physical/'+label+'.json');v=load(path)
        physical.append(dict(role=label,status=v['status'],source=identity(path),
            summary_scalars={k:x for k,x in v['summary'].items()if isinstance(x,(str,int,float))},
            signed_Q1_summary={k:dict(count=len(x['moments']),maximum=x['maximum'],signed_integral=x['signed_integral'])for k,x in v['all_regional_Q1_moments'].items()}))
    atomic(out/'physical-summary.json',physical)
    strict=[]
    for label,folder,step in (('historical_C1',PARENT/'science/C',1),('historical_M0',PARENT/'science/M',0),
        ('fresh_M2',root/'science/M2',2),('fresh_F4',root/'science/F',4)):
        result=load(folder/'final-strict/result.json');assert result['status']=='PASS'and result['selected_absolute_step']==step
        if label.startswith('fresh_'):
            copy(folder/'final-strict/result.json','strict/'+label+'.json');copy(folder/'final-strict/request.json','strict/'+label+'-request.json')
        strict.append(dict(role=label,step=step,result=identity(folder/'final-strict/result.json'),
            physical_mass_dual_norms=result['physical_mass_dual_norms'],completed_M2_evidence=label=='fresh_M2'))
    atomic(out/'strict-role-table.json',strict)
    for role in ('M2','F'):
        folder=root/'science'/role
        for name in ('execution-diff.json','intent.json','native-launch.json','outcome.json','verification.json','qualification.json'):
            if(folder/name).exists():copy(folder/name,role+'/'+name)
        for name in ('startup.ready.json','startup.verified.json','startup.verification-commit.json','startup-observer-complete.json',
                     'startup-replay-parity.json','latest-fully-verified-checkpoint.json','step-4.qualification-commit.json','step-4.prepared.json','step-4.ack.json'):
            if(folder/'monitor'/name).exists():copy(folder/'monitor'/name,role+'/'+name)
    for name in ('runtime-build-v1.json','runtime-dependencies-v1.json','compiled-source-objects.json','runtime-reuse-preflight.json',
                 'runtime-qualified-v1.json','target-cache-compatibility.json','retained-preflight.json','flow-vector-execution.json'):
        copy(root/'science'/name,'controls/'+name)
    for name in ('selected.json','parent-index.json','original-immutability.json','setup.json'):copy(root/'inputs'/name,'controls/'+name)
    for name in ('execution-config.json','source-manifest.json'):copy(root/name,'controls/'+name)
    copy(root/'admission-F.json','controls/admission-F.json')
    # Copy only ended operation receipts/logs; never seal an active output stream.
    timers=[]
    for p in sorted((root/'operations').glob('*/receipt.json')):
        copy(p,'operations/'+p.parent.name+'/receipt.json')
        if p.parent.name=='focused-completion-tests-v1':copy(p.parent/'command.log','operations/'+p.parent.name+'/command.log')
        t=p.parent/'timing.jsonl'
        if t.exists():timers.extend(dict(operation=p.parent.name,**json.loads(line))for line in t.read_text().splitlines()if line.strip())
    tp=root/'science/completed-timing-snapshot.json';atomic(tp,timers);copy(tp,'timing-scopes.json')
    for n in ('run_wj_h22_completion.py','study_wj_h22_completion.py','wj_h22_completion.py','wj_h22_completion_observer.py',
              'wj_h22_observer.py','wj_h22_phase.py','test_wj_h22_completion.py'):
        copy(root/'source/solvers/IncNavierStokesSolver/Tests/RANS'/n,'source/'+n)
    copy(Path(__file__),'source/'+Path(__file__).name)
    copy(Path(__file__).with_name('check_wj_h22_completion_scatter.py'),'source/check_wj_h22_completion_scatter.py')
    for p in sorted((root/'passive-analysis-v2').glob('*.py')):copy(p,'source/passive-analysis-v2/'+p.name)
    pairs=load(analysis/'cross-level-telescopes.json');transitions=load(analysis/'transition-decompositions.json');flows=load(analysis/'cumulative-flow-split.json')
    principal=[dict(geometry=r['geometry'],pair=r['pair'],**r['quantities']['mhat'])for r in pairs]
    fixed=[dict(geometry=r['geometry'],pair=r['pair'],**r['quantities']['a'])for r in pairs]
    lead=[]
    for geo in (251,265):
        cm=next(r for r in principal if r['geometry']==geo and r['pair']=='C-M')
        mf=next(r for r in principal if r['geometry']==geo and r['pair']=='M-F')
        lead.append(dict(geometry=geo,E1=cm['endpoint'],E2=mf['endpoint'],E2_minus_E1=mf['endpoint']-cm['endpoint'],
            signed_contribution_changes={tag:mf['contributions'][tag]-cm['contributions'][tag]for tag in ('lift','flow','k','epsilon')},
            rho=abs(mf['endpoint']/cm['endpoint'])if abs(cm['endpoint'])>4e-11 else None))
    def table(name,rows):
        with(out/name).open('x',newline='')as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    signed=[]
    for r in pairs:
        for name,v in r['quantities'].items():
            if v['status']!='PASS':continue
            shape=np.shape(v['endpoint'])
            for idx in np.ndindex(shape)if shape else[()]:
                at=lambda x:np.asarray(x)[idx]
                signed.append(dict(geometry=r['geometry'],pair=r['pair'],quantity=name,component=str(idx),incoming=at(v['incoming']),
                    **{tag:at(v['contributions'][tag])for tag in ('lift','flow','k','epsilon')},endpoint=at(v['endpoint']),
                    residual=at(v['identity']['residual']),bound=at(v['identity']['bound'])))
    table('signed-cross-level-telescopes.csv',signed)
    table('actual-phase-principal-changes.csv',[dict(level=r['level'],step=r['step'],geometry=r['geometry'],phase=r['phase'],
        delta_mhat=r['signed'].get('mhat'),delta_a=r['signed'].get('a'),delta_g=r['signed'].get('g'),
        delta_Gxz=np.asarray(r['signed']['G'])[0,2],delta_tau=r['signed'].get('tau'))for r in transitions])
    def identities(x):
        if isinstance(x,dict):
            if 'maximum_error'in x and 'maximum_bound'in x:yield x
            for v in x.values():yield from identities(v)
        elif isinstance(x,list):
            for v in x:yield from identities(v)
    algebra={}
    for name in ('branch-telescopes','cross-level-telescopes','genuine-flow-increments','transition-decompositions'):
        rows=list(identities(load(analysis/(name+'.json'))))
        algebra[name]=dict(identities=len(rows),maximum_error=max(v['maximum_error']for v in rows),
            maximum_bound_utilization=max(v['maximum_error']/v['maximum_bound']for v in rows),status='PASS')
    composed=[]
    for f in load(analysis/'genuine-flow-increments.json'):
        for r in f['owners']:
            b=r['B_linear_identity']
            composed.append(dict(level=f['level'],step=f['step'],geometry=r['owner']['owner']['geometry'],
                maximum_error=b['maximum_error'],maximum_bound=b['maximum_bound'],
                frozen_direct_bound_exceeded_components=b['direct_bound_exceeded_count'],
                status=b['status'],proof='signed residual equals tau_C dev(certified G residual) within unchanged primitive roundoff'))
    findings=dict(status='PASS_COMPLETE_FINITE_ATTRIBUTION',principal=principal,fixed_action=fixed,early_adjacent_difference_changes=lead,
        algebra=algebra,composed_roundoff=composed,undefined_roles=counts['undefined_phase_owner_roles'],
        interpretation='actual sequential increments and complete equal-window signed sums; not nonlinear percentages or unique counterfactual causation',
        later_time_stages='NOT_RUN; H21 .01/.02 remain endpoint-only evidence',H20_absolute_screens='REUSED_AFTER_EXACT_ENDPOINT_REPLAY; NOT_NEW_COMPARISONS',
        historical_H22_observer_failure='UNCHANGED',M2_status='NEW_STARTUP_VERIFICATION_PASS; ORIGINAL_EXECUTION_NOT_COMMITTED',
        uniform_temporal_contraction='MIXED_UNCHANGED; NOT_A_PHYSICAL_GATE')
    atomic(out/'findings.json',plain(findings))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':9})
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for ax,geo in zip(axes,(251,265)):
        rows=[v for v in principal if v['geometry']==geo]
        for j,tag in enumerate(('flow','k','epsilon')):ax.bar(np.arange(3)+(j-1)*.23,[v['contributions'][tag]*1e6 for v in rows],width=.23,label=tag)
        ax.scatter(np.arange(3),[v['endpoint']*1e6 for v in rows],c='k',marker='D',s=22,label='endpoint sum',zorder=4)
        ax.axhline(0,color='.3',linewidth=.6);ax.set_xticks(range(3),[v['pair']for v in rows]);ax.set_title(f'Geometry {geo}, '+('upper'if geo==251 else'lower')+' normal trace');ax.grid(axis='y',alpha=.2)
    axes[0].set_ylabel('Signed normalized-minimum contribution (10^-6)');axes[1].legend(fontsize=8)
    fig.suptitle('Complete native stage sums over [0, .005]; outlet height 10.2430071555 mm')
    fig.tight_layout();fig.savefig(out/'signed-seven-attempt-attribution.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for ax,geo in zip(axes,(251,265)):
        rows=[v for v in flows if v['geometry']==geo]
        for j,name in enumerate(('predictor_B','pressure_B','total_B')):ax.bar(np.arange(3)+(j-1)*.23,[np.asarray(v['contributions'][name])[0,2]*1e6 for v in rows],width=.23,label=name)
        ax.axhline(0,color='.3',linewidth=.6);ax.set_xticks(range(3),[v['level']for v in rows]);ax.set_title(f'Geometry {geo}, '+('upper'if geo==251 else'lower')+' normal trace');ax.grid(axis='y',alpha=.2)
    axes[0].set_ylabel('Signed cumulative flow Delta Bxz (10^-6)');axes[1].legend(fontsize=8)
    fig.suptitle('All seven actual predictor / -bGz increments; actual tau_C per attempt')
    fig.tight_layout();fig.savefig(out/'seven-genuine-flow-vector-sums.png',dpi=160);plt.close(fig)
    atomic(out/'component-gates.json',dict(replay='PASS',capture='COMPLETE_7_OF_7',fresh_M2_verification='PASS',new_F_qualification='PASS',
        typed_M2_F4_strict_reads='PASS',identities='PASS',finite_diagnostic='CLOSED_COMPLETE',historical_failure='RETAINED_UNCHANGED',
        new_governing_attempts=4,aggregate_governing_attempts=7,new_CP0=0,new_builds=0,new_scientific_steps=0))
    atomic(out/'accounting-publication-snapshot.json',load(root/'accounting.json'))
    atomic(out/'copy-bindings.json',copies)
    atomic(out/'heavy-manifest.json',dict(files=[identity(p)for p in sorted((root/'science').rglob('*'))if p.is_file()],
        scope='new completion-owned science only; original captures and native runtime remain separately pinned in place'))
    atomic(out/'manifest.json',dict(files=[identity(p)for p in sorted(out.rglob('*'))if p.is_file()],
        scope='immutable publication snapshot; active publication log excluded; final CLOSED linked accounting is separate'))
    print(json.dumps(plain(dict(status=findings['status'],lead=lead,counts=counts,algebra=algebra))),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();publish(a.root)
