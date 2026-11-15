#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact publication of immutable H22 replay/phase evidence, not a solver."""
import argparse,csv,gzip,json,os,shutil
from pathlib import Path
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
import numpy as np
from wj_h15_observer import identity,verify,atomic
from wj_h21_outlet import plain

def publish(root):
    out=root/'publication-v1';out.mkdir()
    derived=root/'science/partial-phase-analysis-v1'
    book=load(root/'accounting.json')
    assert book['target_levels']==['C','M']and book['target_halt']and book['integrity_halt']
    assert not(root/'science/F/native-launch.json').exists()
    copies=[]
    def copy(path,name=None):
        path=Path(path);name=name or path.name;dest=out/name;dest.parent.mkdir(parents=True,exist_ok=True)
        if path.suffix=='.json'and path.stat().st_size>100000:
            dest=Path(str(dest)+'.gz')
            with gzip.GzipFile(filename=str(dest),mode='wb',mtime=0)as f:f.write(path.read_bytes())
            assert gzip.decompress(dest.read_bytes())==path.read_bytes()
        else:shutil.copyfile(path,dest)
        copies.append(dict(original=identity(path),copy=identity(dest),transformation='lossless gzip'if dest.suffix=='.gz'else'byte copy'))
    for p in sorted(derived.glob('*.json')):copy(p)
    roles=load(derived/'roles.json');pairs=load(derived/'cross-level-telescopes.json')
    branches=load(derived/'branch-telescopes.json');flows=load(derived/'cumulative-flow-split.json')
    transitions=load(derived/'transition-decompositions.json');attempts=load(derived/'attempt-bindings.json')
    def table(name,rows):
        with(out/name).open('x',newline='')as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    scalar=[]
    for r in pairs:
        for name,v in r['quantities'].items():
            if v['status']!='PASS':continue
            shape=np.shape(v['endpoint'])
            for idx in np.ndindex(shape)if shape else[()]:
                at=lambda x:np.asarray(x)[idx]
                scalar.append(dict(geometry=r['geometry'],pair=r['pair'],quantity=name,component=str(idx),incoming=at(v['incoming']),
                    **{tag:at(v['contributions'][tag])for tag in ('lift','flow','k','epsilon')},endpoint=at(v['endpoint']),
                    identity_residual=at(v['identity']['residual']),identity_bound=at(v['identity']['bound'])))
    table('signed-cross-level-telescopes.csv',scalar)
    table('actual-phase-principal-changes.csv',[dict(level=r['level'],step=r['step'],geometry=r['geometry'],phase=r['phase'],
        delta_mhat=r['signed'].get('mhat'),delta_a=r['signed'].get('a'),delta_g=r['signed'].get('g'),
        delta_Gxz=np.asarray(r['signed']['G'])[0,2],delta_tau=r['signed'].get('tau'),
        dimensional_term=r['values'].get('dimensional_term'),normalization_term=r['values'].get('normalization_term'))for r in transitions])
    principal=[dict(geometry=r['geometry'],pair=r['pair'],**r['quantities']['mhat'])for r in pairs]
    findings=dict(status='PARTIAL_AFTER_OBSERVER_ROLE_FAILURE',principal=principal,
         early_adjacent_difference_changes='NOT_RUN: no newly captured F phases; old F endpoints cannot replace them',
         undefined_roles=[dict(level=r['level'],step=r['step'],phase=r['phase'],geometry=r['owner']['owner']['geometry'],reason=r['values'].get('reason'))
                          for r in roles if r['values']['constitutive_status']!='DEFINED'],
         provisional_warnings=[dict(level=r['level'],step=r['step'],phase=r['phase'],geometry=r['owner']['owner']['geometry'])
                          for r in roles if r['values'].get('physical_status')=='WARNING'],
         causal_scope='three actual C/M attempts and actual predictor/-bGz increments; no predictor-only constitutive states, no F or later-window stages',
         H20_absolute_screens='C-M historical screen reused after terminal field parity; not a new screen or an F phase audit',
         temporal_contraction='H20/H21 mixed finding unchanged; not an H22 physical gate')
    atomic(out/'findings.json',plain(findings))
    timeline=[];physical=[];heavy=[];timing=[];capture=[];online=[]
    from wj_h17_observer import native_load
    from wj_h22_replay import read_completed,recovery,END
    for z in ('C','M'):
        folder=root/'science'/z
        if(folder/'qualification.json').exists():copy(folder/'qualification.json',z+'/qualification.json')
        copy(folder/'outcome.json',z+'/outcome.json');copy(folder/'execution-diff.json',z+'/execution-diff.json')
        copy(folder/'intent.json',z+'/intent.json');copy(folder/'final-strict/result.json',z+'/final-strict-result.json')
        copy(folder/'final-strict/request.json',z+'/final-strict-request.json')
        for name in ('startup',*(['step-1']if z=='C'else [])):
            copy(folder/'monitor'/(name+'.result.json'),z+'/'+name+'-physical.json')
            r=load(folder/'monitor'/(name+'.result.json'));physical.append(dict(level=z,role=name,
                 status=r['status'],summary=r['summary'],all_regional_Q1_moments=r['all_regional_Q1_moments']))
        for n in range(1,END[z]+1):
            p=folder/'monitor'/f'step-{n}.completion.json';c=load(p)
            r,_=read_completed(recovery(folder/'monitor',n),z,n)
            commit=(folder/'monitor'/f'step-{n}.qualification-commit.json').exists()
            a=r['pimple_attempts'][0]
            online.append(dict(level=z,step=n,**{k:a[k]for k in ('momentum_predictor_solves','pressure_system_solves',
                'k_transport_solves','epsilon_transport_solves','bound_calls','bound_activations','additional_correction_or_restoration_solves')},
                raw_k_sha256=a['k_bound']['raw_sha256'],corrected_k_sha256=a['k_bound']['corrected_sha256']))
            timeline.append(dict(level=z,step=n,time_bits=c['binding']['time_bits'],alpha_bits=c['binding']['alpha_bits'],
                coupled_state_sha256=c['binding']['state_sha256'],history_sha256=c['binding']['callback_history_sha256'],
                native_completed=True,exact_replay='PASS_PASSIVE_RETAINED_BYTES'if z=='M'else'PASS_BEFORE_NATIVE_COMMIT',
                native_full_audit=r['native_full_audit_status'],independent_audit='PASS'if commit else'NOT_RUN'if n==END[z]else'NOT_SCHEDULED',
                native_qualification_commit=commit))
            copy(p,z+'/step-'+str(n)+'.completion.json')
            if n==END[z]:atomic(out/z/'terminal-native-hard-map.json',r['native_hard_audit'])
        for p in (folder/'monitor').glob('*.qualification-commit.json'):copy(p,z+'/'+p.name)
        for p in (folder/'monitor').glob('startup.verification-commit.json'):copy(p,z+'/'+p.name)
        for name in ('latest-native-completed.json','latest-durable-checkpoint.json','latest-fully-verified-checkpoint.json',
                     'first-replay-capture-failure.json','step-2.ack.json'):
            p=folder/'monitor'/name
            if p.exists():copy(p,z+'/'+name)
    for n in range(1,5):timeline.append(dict(level='F',step=n,native_completed=False,exact_replay='NOT_RUN',
        native_full_audit='NOT_RUN',independent_audit='NOT_RUN',native_qualification_commit=False,reason='M observer integrity halt'))
    for p in sorted((root/'science/passive-role-fix-v1').glob('*.json')):copy(p,'passive-correction/'+p.name)
    for a in attempts:
        p=load(verify(a['vectors']['packet']));capture.append(dict(level=a['level'],step=a['step'],vector_packet=p['payload']['bytes'],
             vector_manifest=a['vectors']['packet']['bytes'],phase_json=a['phase_snapshot']['bytes'],sparse_k_matrix=Path(a['matrix']['path']).stat().st_size))
    for p in sorted((root/'operations').glob('*/receipt.json')):
        if p.parent.name=='publication-v1':continue
        copy(p,'operations/'+p.parent.name+'/receipt.json')
        log=p.parent/'command.log'
        if log.exists()and p.parent.name in ('focused-tests-v1','prepared-session-test-v1','cursor-schema-test-v1',
            'vector-native-component-v1','audit-role-regressions-v1','retained-role-and-phase-proof-v1'):
            copy(log,'operations/'+p.parent.name+'/command.log')
        tf=p.parent/'timing.jsonl'
        if tf.exists():
            # Preserve scope definitions and nesting; do not add nested costs.
            for line in tf.read_text().splitlines():
                if line.strip():timing.append(dict(operation=p.parent.name,**json.loads(line)))
    for name in ('runtime-build-v1.json','runtime-dependencies-v1.json','compiled-source-objects.json','flow-vector-execution.json',
                 'initial-and-oracle-lineage.json','build-result.json','object-reuse-proof.json','command-admission-refusal-v1.json',
                 'prelaunch-cursor-correction.json','runtime-qualified-v1.json','target-cache-compatibility.json'):
        copy(root/'science'/name,'controls/'+name)
    for name in ('execution-config.json','execution-config-v2.json','source-intent.json'):copy(root/name,'controls/'+name)
    copy(root/'inputs/selected.json','controls/selected-inputs.json')
    for folder in ('prelaunch-v1','passive-code-v3'):
        for p in sorted((root/folder).glob('*')):
            if p.is_file():copy(p,'source-versions/'+folder+'/'+p.name)
    for name in ('wj_h22_observer.py','study_wj_h22.py','wj_h22_replay.py','wj_h22_phase.py','run_wj_h22.py'):
        p=root/'source/solvers/IncNavierStokesSolver/Tests/RANS'/name
        copy(p,'source-versions/executed/'+name)
    copy(Path(__file__),'source-versions/publication/'+Path(__file__).name)
    for p in sorted((root/'science').rglob('*')):
        if p.is_file()and not p.is_symlink():heavy.append(identity(p))
    atomic(out/'timeline.json',timeline);atomic(out/'online-counts.json',online);atomic(out/'physical-gates-and-signed-moments.json',physical)
    atomic(out/'capture-bytes.json',capture);atomic(out/'timing-scopes.json',timing)
    atomic(out/'heavy-manifest.json',dict(files=heavy,scope='owned new H22 science only; original H18/H19/H20 inputs remain separately pinned'))
    atomic(out/'copy-bindings.json',copies)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':9})
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for ax,geo in zip(axes,(251,265)):
        rr=[v for v in principal if v['geometry']==geo]
        for j,tag in enumerate(('flow','k','epsilon')):
            ax.bar(np.arange(len(rr))+(j-1)*.23,[r['contributions'][tag]*1e6 for r in rr],width=.23,label=tag)
        ax.scatter(np.arange(len(rr)),[r['endpoint']*1e6 for r in rr],c='k',marker='D',s=22,label='endpoint sum',zorder=4)
        ax.axhline(0,color='0.3',linewidth=.6);ax.set_xticks(range(len(rr)),[r['pair']for r in rr]);ax.set_title(f'Geometry {geo}, '+('upper'if geo==251 else'lower')+' normal trace');ax.grid(axis='y',alpha=.2)
    axes[0].set_ylabel('Signed normalized-minimum contribution (10^-6)');axes[1].legend(fontsize=8)
    fig.suptitle('Retained C/M stage sums over [0, .005]; outlet height 10.2430071555 mm')
    fig.tight_layout();fig.savefig(out/'signed-native-stage-attribution.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for ax,geo in zip(axes,(251,265)):
        rr=[v for v in flows if v['geometry']==geo]
        for j,name in enumerate(('predictor_B','pressure_B','total_B')):
            ax.bar(np.arange(len(rr))+(j-1)*.23,[np.array(r['contributions'][name])[0,2]*1e6 for r in rr],width=.23,label=name)
        ax.axhline(0,color='0.3',linewidth=.6);ax.set_xticks(range(len(rr)),[r['level']for r in rr]);ax.set_title(f'Geometry {geo}, '+('upper'if geo==251 else'lower')+' normal trace');ax.grid(axis='y',alpha=.2)
    axes[0].set_ylabel('Signed cumulative flow Delta Bxz (10^-6)');axes[1].legend(fontsize=8)
    fig.suptitle('Three captured C/M predictor / -bGz increments; actual tau_C at each attempt')
    fig.tight_layout();fig.savefig(out/'genuine-predictor-pressure-split.png',dpi=160);plt.close(fig)
    atomic(out/'accounting-publication-snapshot.json',load(root/'accounting.json'))
    atomic(out/'component-gates.json',dict(overall='OBSTRUCTED_OBSERVER_AUDIT_ROLE',replay_fidelity='PASS_FOR_THREE_RETAINED_ATTEMPTS',
        capture_completeness='PARTIAL_3_OF_7',phase_and_flow_identities='PASS_FOR_C_M',
        physical_qualification='C1_COMMITTED; M2_NATIVE_FULL_PASS_BUT_INDEPENDENT_NOT_RUN_UNQUALIFIED; F_NOT_RUN',
        passive_role_correction='PASS_PASSIVE_ONLY; CORRECTED_OBSERVER_NATIVE_QUALIFICATION_NOT_RUN',
        diagnostic_question='INCOMPLETE; NO_F_STAGE_ATTRIBUTION',uniform_temporal_contraction='NOT_REQUIRED_MIXED_HISTORICAL_FINDING',
        maximum_new_scientific_time=0,followup_executed=False))
    atomic(out/'manifest.json',dict(files=[identity(p)for p in sorted(out.rglob('*'))if p.is_file()],
        scope='H22 compact derived/copy evidence; publication accounting is OPEN, final closure is separate'))
    print(plain(findings),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();publish(a.root)
