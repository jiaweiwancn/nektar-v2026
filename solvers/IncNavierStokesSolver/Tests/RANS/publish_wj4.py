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
from wj_immutable import file_hash,artifacts,read,write_once


def csv_file(path,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    root=a.root.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    source=root/'source';build=root/'build';git=lambda *v:subprocess.check_output(['git','-c','core.filemode=false','-C',str(source),*v],text=True).strip()
    tracked=[source/f for f in subprocess.check_output(['git','-C',str(source),'ls-files'],text=True).splitlines() if (source/f).is_file()]
    content={str(f.relative_to(source)):file_hash(f) for f in tracked};write_once(root/'wj4-source-content.json',content)
    runtime={}
    for name in ('IncNavierStokesSolver','RANSWallResolutionProbe','WJKernelProbe'):
        from wj_immutable import runtime as identity
        runtime[name]=identity(build/name)
    index={}
    for scope in (root/'study',root/'build/wj4-observations',root/'build/wj3-smoke',root/'build/wj2-native-EP',root/'build/wj2-native-Budget'):
        if scope.exists():
            for f in scope.rglob('*'):
                if f.is_file() and (f.suffix in ('.json','.xml','.bin','.log') or f.name=='State.bin'):index[str(f)]=dict(bytes=f.stat().st_size,sha256=file_hash(f))
    write_once(root/'wj4-artifact-index.json',index)
    inventory=read(root/'inventory.json');tests={}
    for r in read(root/'command-receipts.json'):
        if not r['label'].startswith('ctest-'):continue
        matches=re.findall(r'(\d+)% tests passed, (\d+) tests failed out of (\d+)',Path(r['log']).read_text())
        tests[r['label']]=dict(exit_code=r['exit_code'],result=matches[-1] if matches else None,log=r['log'],sha256=r['log_sha256'],seconds=r['seconds'])
    manifest=dict(schema='WJ4LinuxManifestV1',root=str(root),source_revision=git('rev-parse','HEAD'),source_tree=git('rev-parse','HEAD^{tree}'),source_status=git('status','--short'),
        source_content=dict(path=str(root/'wj4-source-content.json'),sha256=file_hash(root/'wj4-source-content.json')),
        specification_sha256={name:content['docs/literature/'+name] for name in ('wj_earsm_neutral_abl_ep_awf_note.tex','wj_earsm_budget_cutoff_awf_appendix.tex')},
        runtime=runtime,arithmetic_policy='Release C++17 -ffp-contract=off -fno-fast-math; all numerical library threads=1',
        effective_flags={name:(build/'CMakeFiles'/ (name+'.dir')/'flags.make').read_text() for name in ('IncNavierStokesSolver','RANSWallResolutionProbe','RANSKernels')},
        inventory=dict(count=len(inventory['tests']),path=str(root/'inventory.json'),sha256=file_hash(root/'inventory.json')),
        artifact_index=dict(count=len(index),path=str(root/'wj4-artifact-index.json'),sha256=file_hash(root/'wj4-artifact-index.json')),
        command_receipts=dict(path=str(root/'command-receipts.json'),sha256=file_hash(root/'command-receipts.json')),tests=tests,
        ledger={policy:dict(path=str(root/('ledger-'+policy+'.json')),sha256=file_hash(root/('ledger-'+policy+'.json'))) for policy in ('EP','Budget') if (root/('ledger-'+policy+'.json')).exists()},
        allocation_receipts=artifacts(list(root.glob('*allocation*.json'))+list(root.glob('processes*.txt'))),
        configuration=artifacts([root/'configure.log',build/'CMakeCache.txt']),finite_scope='reference-Dirichlet strip [0,1]x[0,2]; serial flat2D p3; finite observation and angular inventories; no global or long-time theorem')
    write_once(root/'wj4-manifest.json',manifest);write_once(out/'manifest.json',manifest)
    rows=[];interfaces=[];windows=[];differences=[];residuals=[];horizontal=[];audits=[];witnesses={};policies={}
    for policy in ('EP','Budget'):
        path=root/'study'/policy/'summary.json'
        if not path.exists():continue
        s=read(path);policies[policy]={k:s.get(k) for k in ('qualification','components','longer_interpretation','stationarity','ledger_totals','cost_design','mesh_inventory','arithmetic_envelope')}
        policies[policy]['heavy_summary']=dict(path=str(path),sha256=file_hash(path),bytes=path.stat().st_size)
        for group in ('horizontal','longer'):
            for key,trajectory in s.get(group,{}).items():
                for c in [trajectory['initial']]+trajectory.get('trajectory',[]):
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
        for key,c in s.get('horizontal_comparisons',{}).items():
            for scope in ('full','central','constrained_boundary'):
                for field,r in c[scope]['fields'].items():differences.append(dict(policy=policy,kind='horizontal',index=key,scope=scope,field=field,**{k:v for k,v in r.items() if k!='maximum_provenance'}))
        for key,c in s.get('timestep_comparisons',{}).items():
            for scope in ('full','central','constrained_boundary'):
                for field,r in c['comparison'][scope]['fields'].items():differences.append(dict(policy=policy,kind='timestep',index=key,time=c['time'],fine_time=c['fine_time'],scope=scope,field=field,**{k:v for k,v in r.items() if k!='maximum_provenance'}))
        for r in s.get('stationarity',{}).get('windows',[]):windows.append(dict(policy=policy,trajectory=r['trajectory'],steps=json.dumps(r['steps']),times=json.dumps(r['times']),status=r['status'],**r['V']))
        write_once(out/(policy+'.json'),policies[policy])
    for name,data in (('fields',rows),('interfaces',interfaces),('windows',windows),('differences',differences),('residuals',residuals),('horizontal-spread',horizontal),('audit-margins',audits)):csv_file(out/(name+'.csv'),data)
    write_once(out/'audit-witnesses.json',witnesses)
    receipt=dict(schema='WJ4PublicationV1',manifest_sha256=file_hash(root/'wj4-manifest.json'),policies={p:v['qualification'] for p,v in policies.items()},
        tables=artifacts(list(out.glob('*.csv'))),figures=artifacts(list((root/'figures').glob('*'))),
        evidence_bridge=artifacts([root/'wj4-manifest.json',root/'wj4-artifact-index.json',root/'wj4-source-content.json',root/'command-receipts.json']),
        measurement_only_requalification=artifacts(list((root/'requalification').glob('*.json'))),
        original_wj3_receipts=artifacts([source/'docs/plans/wj3-linux'/name for name in ('publication.json','receipt.json','manifest.json','EP.json','Budget.json')]))
    write_once(out/'receipt.json',receipt)
    if (root/'figures').exists():shutil.copytree(root/'figures',out/'figures')
    for file in out.rglob('*'):
        if file.is_file():file.chmod(0o644)
    print(out/'receipt.json')

if __name__=='__main__':main()
