#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Seal compact WJ-4 evidence without rewriting command or prior WJ receipts."""
import argparse
import csv
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
from wj_immutable import file_hash,artifacts,read,write_once,digest,verify_artifacts


def csv_file(path,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)

def native_accounting(summary,policy):
    """Summarize every original native attempt, including deliberate rejection.

    Full per-step telemetry stays in immutable producer reports. The compact
    bridge binds their bytes and records partial counts without filling them in.
    """
    from run_cls_epawf_picard import state_hash
    counts=('momentum_predictor_solves','pressure_system_solves','k_transport_solves',
            'epsilon_transport_solves','bound_calls','additional_governing_equation_solves',
            'additional_correction_or_restoration_solves')
    rows=[];lineage={}
    for name,entry in summary.get('execution',{}).items():
        path=Path(entry['report']);receipt=Path(entry['receipt']);r=read(receipt)
        verify_artifacts(r['artifacts'])
        if not path.exists():
            rows.append(dict(policy=policy,case=name,report=str(path),status=entry['status'],attempts=r['actual_attempts'],failure=r.get('failure')));continue
        if file_hash(path)!=entry['report_sha256']:raise ValueError('producer report changed: '+str(path))
        d=read(path);attempts=d['pimple_attempts'];bounds=[a['k_bound'] for a in attempts if 'k_bound' in a]
        if r['actual_attempts']!=len(attempts):raise ValueError('native attempt count differs from receipt')
        complete=[a for a in attempts if a.get('scheduled_complete')]
        if any(any(a.get(k)!=1 for k in counts[:5]) or any(a.get(k)!=0 for k in counts[5:]) for a in complete):
            raise ValueError('scheduled one-sweep count violated: '+str(path))
        rows.append(dict(policy=policy,case=name,status=entry['status'],report=str(path),report_sha256=file_hash(path),
            receipt_sha256=file_hash(receipt),xml_sha256=file_hash(path.parent/'case.xml'),steps=d['steps'],time=d['time'],dt=d['dt'],
            unknowns=len(d['state']),attempts=len(attempts),fully_scheduled_attempts=len(complete),
            partial_attempts=len(attempts)-len(complete),failure=d.get('failure'),
            **{k:sum(a.get(k,0) for a in attempts) for k in counts},
            bound_activations=sum(a.get('bound_activations',0) for a in attempts),
            maximum_mass_correction=max((b['mass_norm_correction'] for b in bounds),default=0),
            minimum_arithmetic_reserve=min((b['arithmetic_reserve'] for b in bounds),default=None),
            raw_corrected_k_hashes_equal=all(b['raw_sha256']==b['corrected_sha256'] for b in bounds),
            maximum_added_k_storage_residual=max((abs(a.get('k_added_storage_residual',0)) for a in attempts),default=0)))
        stage_present=bool(d['last_stage_history'])
        accepted_callback=stage_present and d['implicit_callback_count']>0 and d['last_stage_alpha']>0 and abs(d['last_stage_time']-d['time'])<=2e-15
        callback_context=('failed_attempt_not_accepted' if d.get('failure') and stage_present else
                          'completed_accepted' if accepted_callback else 'inapplicable_startup_or_frozen_readback')
        lineage[policy+'/'+name]=dict(report=str(path),report_sha256=file_hash(path),state_sha256=state_hash(d['state']),
            callback_context=callback_context,usable_as_accepted_stage=callback_context=='completed_accepted',
            native_history_sha256=digest(d['native_history']),native_history_times=d['native_history_times'],
            last_stage_history_sha256=state_hash(d['last_stage_history']) if d['last_stage_history'] else None,
            last_stage_time=d['last_stage_time'],last_stage_alpha=d['last_stage_alpha'],
            accepted_audits_sha256=digest(d['wj_audit_packets']),physical_samples_sha256=digest(d['physical_samples']),
            checkpoint=artifacts([Path(d['checkpoint'])/'State.bin']),signature_sha256=digest(d['signature']),
            failed_candidate_sha256=state_hash(d['rejected_state']) if d.get('rejected_state') else None)
    return rows,lineage

def load_tables(c,context,loads,weak):
    for name,error in c.get('independent_loads',{}).get('load_max_absolute_errors',{}).items():
        loads.append(dict(**context,load=name,max_absolute_error=error,
            quadrature_points=c['independent_loads']['quadrature_points'],state_sha256=c['state_sha256']))
    for mode,r in c.get('independent_weak_modes',{}).get('modes',{}).items():
        for region,fields in zip(('inner','outer'),r['regional_loads']):
            for field,terms in zip(('u','v','k','outer_epsilon'),fields):
                for term,value in zip(('storage','convection','stress_or_diffusion','negative_source','roughness'),terms):
                    weak.append(dict(**context,mode=mode,region=region,field=field,term=term,physical_load=value,max_independent_absolute_error=r['max_absolute_error']))
        for term,values in (('signed_volume_sum',r['signed_sum']),('prescribed_exterior_reaction',r['prescribed_exterior_reaction'])):
            for field,value in zip(('u','v','k','outer_epsilon'),values):
                weak.append(dict(**context,mode=mode,region='inner+outer',field=field,term=term,physical_load=value,max_independent_absolute_error=r['max_absolute_error']))

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--observer-root',type=Path);p.add_argument('--study',type=Path);p.add_argument('--figures',type=Path);a=p.parse_args()
    root=a.root.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    meta=root/'publication-lineage'/out.name;meta.mkdir(parents=True,exist_ok=False)
    figures=a.figures.resolve() if a.figures else root/'figures'
    observer=a.observer_root.resolve() if a.observer_root else root;study=a.study.resolve() if a.study else root/'study'
    source=observer/'source';build=observer/'build';git=lambda *v:subprocess.check_output(['git','-c','core.filemode=false','-C',str(source),*v],text=True).strip()
    tracked=[source/f for f in subprocess.check_output(['git','-C',str(source),'ls-files'],text=True).splitlines() if (source/f).is_file()]
    content={str(f.relative_to(source)):file_hash(f) for f in tracked};write_once(meta/'source-content.json',content)
    runtime={}
    for name in ('IncNavierStokesSolver','RANSWallResolutionProbe','WJKernelProbe'):
        from wj_immutable import runtime as identity
        runtime[name]=identity(build/name)
    index={}
    for scope in (root/'study',study,observer/'build/wj4-observations',observer/'build/wj3-smoke',observer/'build/wj2-native-EP',observer/'build/wj2-native-Budget',root/'build/wj4-observations',root/'build/wj3-smoke',root/'build/wj2-native-EP',root/'build/wj2-native-Budget'):
        if scope.exists():
            for f in scope.rglob('*'):
                if f.is_file() and (f.suffix in ('.json','.xml','.bin','.log') or f.name=='State.bin'):index[str(f)]=dict(bytes=f.stat().st_size,sha256=file_hash(f))
    write_once(meta/'artifact-index.json',index)
    inventory=read(observer/'inventory.json');tests={}
    for r in read(root/'command-receipts.json')+(read(observer/'command-receipts.json') if observer!=root else []):
        if not r['label'].startswith('ctest-'):continue
        matches=re.findall(r'(\d+)% tests passed, (\d+) tests failed out of (\d+)',Path(r['log']).read_text())
        tests[str(Path(r['log']).parent)+'/'+r['label']]=dict(exit_code=r['exit_code'],result=matches[-1] if matches else None,log=r['log'],sha256=r['log_sha256'],seconds=r['seconds'])
    manifest=dict(schema='WJ4LinuxManifestV1',root=str(root),source_revision=git('rev-parse','HEAD'),source_tree=git('rev-parse','HEAD^{tree}'),source_status=git('status','--short'),
        source_content=dict(path=str(meta/'source-content.json'),sha256=file_hash(meta/'source-content.json')),
        specification_sha256={name:content['docs/literature/'+name] for name in ('wj_earsm_neutral_abl_ep_awf_note.tex','wj_earsm_budget_cutoff_awf_appendix.tex')},
        runtime=runtime,arithmetic_policy='Release C++17 -ffp-contract=off -fno-fast-math; all numerical library threads=1',
        effective_flags={name:(build/'CMakeFiles'/ (name+'.dir')/'flags.make').read_text() for name in ('IncNavierStokesSolver','RANSWallResolutionProbe','RANSKernels')},
        inventory=dict(count=len(inventory['tests']),path=str(observer/'inventory.json'),sha256=file_hash(observer/'inventory.json')),
        artifact_index=dict(count=len(index),path=str(meta/'artifact-index.json'),sha256=file_hash(meta/'artifact-index.json')),
        command_receipts=dict(path=str(root/'command-receipts.json'),sha256=file_hash(root/'command-receipts.json')),tests=tests,
        ledger={policy:dict(path=str(root/('ledger-'+policy+'.json')),sha256=file_hash(root/('ledger-'+policy+'.json'))) for policy in ('EP','Budget') if (root/('ledger-'+policy+'.json')).exists()},
        allocation_receipts=artifacts(list(root.glob('*allocation*.json'))+list(root.glob('processes*.txt'))+list(observer.glob('*allocation*.json'))),
        configuration=artifacts([root/'configure.log',build/'CMakeCache.txt']),
        producer_runtime=read(root/'study/EP/design.json')['runtime'],
        accepted_readback_runtime=read(study/'EP/design.json')['runtime'],
        observer_root=str(observer),study_root=str(study),producer_root=str(root),
        observer_commands=artifacts([observer/'command-receipts.json']),source_equivalence=artifacts(list(root.glob('*equivalence*.json'))+list(observer.glob('*equivalence*.json'))),finite_scope='reference-Dirichlet strip [0,1]x[0,2]; serial flat2D p3; finite observation and angular inventories; no global or long-time theorem')
    manifest['numerical_environment']=dict(python=sys.version,numpy=__import__('numpy').__version__,scipy=__import__('scipy').__version__,
        sympy=__import__('sympy').__version__,compiler=subprocess.check_output(['c++','--version'],text=True),platform=subprocess.check_output(['uname','-a'],text=True),
        allocation=dict(initial_build_workers=1,later_build_workers=4,policy_concurrency=2,test_parallel=1,ranks=1,numerical_threads=1,total_worker_cap=72))
    manifest['independent_proofs']=artifacts([root/'regression-verification.json',root/'shared-legacy-accounting-supplement.json',root/'legacy-control-accounting-supplement.json',
        *root.glob('independent-gradient-*/receipt.json'),*root.glob('independent-gradient-*/result.json'),
        *root.glob('all-step-transactions-*.json'),root/'parent-readback-proof/EP/proof.json',root/'arithmetic-cache-proof/EP/proof.json'])
    for context in (manifest['producer_runtime'],manifest['accepted_readback_runtime']):
        for value in context.values():
            if file_hash(value['path'])!=value['sha256']:raise ValueError('recorded runtime binary changed')
            for path,h in value['libraries'].items():
                if file_hash(path)!=h:raise ValueError('recorded dependency changed')
    write_once(meta/'manifest.json',manifest);write_once(out/'manifest.json',manifest)
    rows=[];interfaces=[];windows=[];differences=[];residuals=[];horizontal=[];audits=[];witnesses={};policies={};runs=[];stages={}
    loads=[];weak=[];regional=[];trace=[];arithmetic=[];hypotheses={}
    for policy in ('EP','Budget'):
        path=study/policy/'summary.json'
        if not path.exists():continue
        s=read(path);policies[policy]={k:s.get(k) for k in ('qualification','components','longer_interpretation','stationarity','ledger_totals','cost_design','mesh_inventory','arithmetic_envelope')}
        # Cost scouts are inherited, not re-executed during observer replay.
        # Account the complete original producer inventory, including scouts.
        accounting,packets=native_accounting(read(root/'study'/policy/'summary.json'),policy);runs.extend(accounting);stages.update(packets)
        policies[policy]['read_only_reprocessing']=s.get('read_only_reprocessing')
        policies[policy]['all_step_transaction_verification']=read(root/('all-step-transactions-'+policy+'.json'))
        policies[policy]['publication_snapshot_ledger_totals']=__import__('wj_immutable').Ledger(root,policy).totals()
        policies[policy]['heavy_summary']=dict(path=str(path),sha256=file_hash(path),bytes=path.stat().st_size)
        for group in ('horizontal','longer'):
            for key,trajectory in s.get(group,{}).items():
                for c in [trajectory['initial']]+trajectory.get('trajectory',[]):
                    context=dict(policy=policy,group=group,trajectory=key,steps=c['steps'],time=c['time'])
                    load_tables(c,context,loads,weak)
                    if c.get('continuum_hypotheses'):hypotheses[policy+'/'+group+'/'+key]=c['continuum_hypotheses']
                    for system,r in c.get('arithmetic',{}).get('systems',{}).items():arithmetic.append(dict(**context,system=system,**r))
                    for row in c['errors']['interface_trace_moments']['moments']:
                        for i,field in enumerate(('traction_x','traction_y','traction_z','k_flux')):trace.append(dict(**context,mode=row['mode'],field=field,
                            native_Gauss8=row['native_Gauss8'][i],independent_composite=row['independent_composite'][i],max_absolute_difference=row['max_absolute_difference']))
                    for kind in ('initial_representation','continuum_error','initial_drift'):
                        for scope in ('full','central','constrained_boundary'):
                            for field,r in c['errors'][kind][scope]['fields'].items():
                                prov=r['maximum_provenance'];rows.append(dict(policy=policy,group=group,trajectory=key,steps=c['steps'],time=c['time'],kind=kind,scope=scope,field=field,
                                    L2=r['L2'],sampled_Linf=r['sampled_Linf'],measure=r['measure'],normalized_L2=r['normalized_L2'],normalized_sampled_Linf=r['normalized_sampled_Linf'],maximum_position=json.dumps(prov['position']),
                                    geometry=prov['geometry'],element=prov['element'],side=prov['side'],requested_x_side=prov['requested_x_side']))
                    for j,field in enumerate(('u','v','k','outer_epsilon','reduced_pressure','mechanical_pressure')):
                        span=max(c['errors']['horizontal_samples'],key=lambda r:r['sampled_horizontal_spread'][j])
                        horizontal.append(dict(policy=policy,group=group,trajectory=key,steps=c['steps'],time=c['time'],field=field,
                            maximum_horizontal_spread=span['sampled_horizontal_spread'][j],height=span['z'],side=span['side'],mean_at_height=span['sampled_horizontal_mean'][j],
                            mean_rule=span['mean_rule'],vertical_velocity_max=max(r['vertical_velocity_max'] for r in c['errors']['horizontal_samples'])))
                    raw=read(Path(c['path']))['measurements'];packet=raw['audit']
                    for region,values in zip(('inner','outer'),raw['regional_spatial_mass_dual_norms']):
                        for block,value in zip(('u','v','k','entire_outer_epsilon','inner_continuity','outer_continuity','gauge','boundary_trace'),values):
                            regional.append(dict(**context,region=region,block=block,norm=value,scope='native region-restricted operator with same constraints; separately normed pieces are not additive'))
                    for condition,r in packet['conditions'].items():audits.append(dict(policy=policy,group=group,trajectory=key,steps=c['steps'],time=c['time'],
                        condition=condition,passed=r['passed'],margins=json.dumps(r['margins']),state_sha256=c['state_sha256']))
                    witness_key=policy+'/'+group+'/'+key+'/'+str(c['steps'])
                    witnesses[witness_key]=dict(state_sha256=c['state_sha256'],stage_history_sha256=c['history_sha256'],
                        principal=packet['principal_witnesses'],gradient_replay=c.get('gradient_replay'),source=dict(path=c['path'],sha256=c['sha256']))
                    interfaces.append(dict(policy=policy,group=group,trajectory=key,steps=c['steps'],time=c['time'],**{k:v for k,v in c['errors']['interface'].items() if not isinstance(v,list)},
                        sampled_K_min=c['errors']['interface']['sampled_K'][0],sampled_K_max=c['errors']['interface']['sampled_K'][1],
                        cutoff_min=c['errors']['interface']['cutoff'][0],cutoff_max=c['errors']['interface']['cutoff'][1],divergence=c['errors']['strong_divergence'],central_divergence=c['errors']['divergence_by_scope']['central']['normalized_maximum'],
                        boundary_divergence=c['errors']['divergence_by_scope']['constrained_boundary']['normalized_maximum'],accuracy_status=c['errors']['accuracy_gate']['status']))
                    for i,name in enumerate(('u','v','k','entire_outer_epsilon','inner_continuity','outer_continuity','gauge','boundary_trace')):
                        residuals.append(dict(policy=policy,group=group,trajectory=key,steps=c['steps'],time=c['time'],block=name,spatial=c['spatial_residual'][i],stage=c['stage_defect'][i] if c['stage_defect'] else None))
                for name in ('independent_obstruction_loads',):
                    c=trajectory.get(name)
                    if c:load_tables(c,dict(policy=policy,group=group,trajectory=key,steps=c['steps'],time=c['time']),loads,weak)
        for key,c in s.get('horizontal_comparisons',{}).items():
            for scope in ('full','central','constrained_boundary'):
                for field,r in c[scope]['fields'].items():differences.append(dict(policy=policy,kind='horizontal',index=key,scope=scope,field=field,**{k:v for k,v in r.items() if k!='maximum_provenance'}))
        for key,c in s.get('timestep_comparisons',{}).items():
            for scope in ('full','central','constrained_boundary'):
                for field,r in c['comparison'][scope]['fields'].items():differences.append(dict(policy=policy,kind='timestep',index=key,time=c['time'],fine_time=c['fine_time'],scope=scope,field=field,**{k:v for k,v in r.items() if k!='maximum_provenance'}))
        for r in s.get('stationarity',{}).get('windows',[]):windows.append(dict(policy=policy,trajectory=r['trajectory'],steps=json.dumps(r['steps']),times=json.dumps(r['times']),status=r['status'],**r['V']))
        write_once(out/(policy+'.json'),policies[policy])
    for name,data in (('fields',rows),('interfaces',interfaces),('windows',windows),('differences',differences),('residuals',residuals),('horizontal-spread',horizontal),('audit-margins',audits),
        ('independent-loads',loads),('weak-loads',weak),('regional-residuals',regional),('trace-moments',trace),('arithmetic',arithmetic)):csv_file(out/(name+'.csv'),data)
    write_once(out/'startup-hypotheses.json',hypotheses)
    csv_file(out/'run-accounting.csv',runs);write_once(out/'accepted-stage-lineage.json',stages)
    write_once(out/'audit-witnesses.json',witnesses)
    for path in [root/'regression-verification.json',root/'shared-legacy-accounting-supplement.json',root/'legacy-control-accounting-supplement.json',root/'command-receipts.json',
                 root/'ledger-EP.json',root/'ledger-Budget.json',*root.glob('all-step-transactions-*.json')]:
        if path.exists():shutil.copyfile(path,out/path.name)
    if observer!=root:shutil.copyfile(observer/'command-receipts.json',out/'observer-command-receipts.json')
    receipt=dict(schema='WJ4PublicationV1',manifest_sha256=file_hash(meta/'manifest.json'),policies={p:v['qualification'] for p,v in policies.items()},
        tables=artifacts(list(out.glob('*.csv'))),figures=artifacts(list(figures.glob('*'))),
        evidence_bridge=artifacts([meta/'manifest.json',meta/'artifact-index.json',meta/'source-content.json',root/'command-receipts.json']),
        measurement_only_requalification=artifacts(list((root/'requalification').glob('*.json'))),
        original_wj3_receipts=artifacts([source/'docs/plans/wj3-linux'/name for name in ('publication.json','receipt.json','manifest.json','EP.json','Budget.json')]))
    receipt['compact_identities']=artifacts(list(out.glob('*.json')))
    write_once(out/'receipt.json',receipt)
    if figures.exists():shutil.copytree(figures,out/'figures')
    for file in out.rglob('*'):
        if file.is_file():file.chmod(0o644)
    print(out/'receipt.json')

if __name__=='__main__':main()
