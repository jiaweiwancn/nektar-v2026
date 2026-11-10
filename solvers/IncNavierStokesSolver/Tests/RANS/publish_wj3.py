#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Publish compact WJ-3 receipts/tables; retain heavy native data at immutable paths."""
import argparse,csv,json,math,re,shutil
from pathlib import Path
from run_cls_epawf_picard import load,sha,write_json,require


def table(path,rows):
    if not rows:return
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    root=a.root.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
    manifest=load(root/'manifest.json');summaries={policy:load(root/'study-final'/policy/'summary.json') for policy in ('EP','Budget')}
    spatial=[];temporal=[];drift=[];residuals=[];identities=[];witnesses={}
    labels=('u','v','k','entire_outer_epsilon','inner_continuity','outer_continuity','gauge','boundary_trace')
    for policy,s in summaries.items():
        compact={k:v for k,v in s.items() if k!='receipts'}
        compact['heavy_summary']=dict(path=str(root/'study-final'/policy/'summary.json'),sha256=sha(root/'study-final'/policy/'summary.json'))
        write_json(out/(policy+'.json'),compact)
        for row in s['spatial']:
            e=row['comparison']['field_errors'];r=dict(policy=policy,level=row['level'],unknowns=row['unknowns'])
            r.update(dict(zip(('spatial_'+k for k in labels),row['raw_reference_norms'])))
            r['first_step_pi_drift_L2']=row['first_step_pressure_adjustment']['fields']['reduced_pressure']['L2']
            r.update({k+'_initial_L2':e['initial_representation']['fields'][k]['L2'] for k in ('u','k','outer_epsilon')});spatial.append(r)
        for field,row in s.get('temporal',{}).get('rates',{}).items():
            temporal.append(dict(policy=policy,field=field,T=s['frozen_design']['T'],T_ustar_over_H=s['frozen_design']['T_ustar_over_H'],**row))
        for row in s.get('main',[]):
            if 'comparison' not in row:continue
            e=row['comparison']['field_errors']
            for kind in ('initial_representation','initial_drift','continuum_error'):
                for field,value in e[kind]['fields'].items():
                    drift.append(dict(policy=policy,steps=row['steps'],dt=row['dt'],T=s['frozen_design']['T'],kind=kind,field=field,**value))
        folder=root/'study-final'/policy
        for file in sorted(folder.rglob('identity.json')):
            d=load(file);identities.append(dict(policy=policy,case=str(file.parent),**d))
        for file in sorted(folder.glob('*/comparison.json')):
            d=load(file);raw=load(Path(d['path']))['measurements']
            row=dict(policy=policy,measurement=str(file.parent),time=raw['time'],state_sha256=raw['state_sha256'],
                stage_history_sha256=raw['stage_history_sha256'],stage_alpha=raw['completed_stage']['alpha_mass'] if raw['completed_stage'] else None)
            row.update(dict(zip(('spatial_'+k for k in labels),raw['spatial_mass_dual_norms'])))
            row.update(dict(zip(('stage_'+k for k in labels),raw['completed_stage']['mass_dual_norms'] if raw['completed_stage'] else [None]*8)));residuals.append(row)
        witnesses[policy]={}
        for case in ('spatial-measure-L2','main-128-measure','cross-128-measure'):
            file=folder/case/'result.json'
            if file.exists():
                data=load(file)['measurements']['audit'];data.pop('callback_history',None)
                witnesses[policy][case]=dict(path=str(file),sha256=sha(file),packet=data)
    table(out/'spatial.csv',spatial);table(out/'temporal.csv',temporal);table(out/'drift.csv',drift);table(out/'residuals.csv',residuals);table(out/'case-identities.csv',identities)
    write_json(out/'principal-witnesses.json',witnesses)
    shutil.copy2(root/'inventory.log',out/'test-inventory.json')
    shutil.copy2(root/'command-receipts.json',out/'command-receipts.json')
    # Shared library hashes are deduplicated while every executable retains its
    # exact linked set. No library is copied/replaced or relinked by publication.
    libraries={};runtime={}
    for name,value in manifest['runtime'].items():
        libraries.update(value['resolved_libraries']);runtime[name]={k:v for k,v in value.items() if k not in ('ldd','resolved_libraries')}
        runtime[name]['resolved_library_paths']=list(value['resolved_libraries'])
    compact_manifest={k:v for k,v in manifest.items() if k not in ('runtime','command_receipts')}
    compact_manifest.update(runtime=runtime,resolved_libraries=libraries,heavy_manifest=dict(path=str(root/'manifest.json'),sha256=sha(root/'manifest.json')))
    write_json(out/'manifest.json',compact_manifest)
    attempts=dict(EP=0,Budget=0);attempt_index=[]
    for scope in (root/'study',root/'study-final',root/'build/wj3-smoke',root/'build/wj2-native-EP',root/'build/wj2-native-Budget'):
        if not scope.exists():continue
        for file in sorted(scope.rglob('case.evolution.json')):
            d=load(file)
            if d.get('turbulence_model')!='WJHighRe':continue
            policy='Budget' if 'Budget' in d['wall_policy'] else 'EP';n=len(d['pimple_attempts']);attempts[policy]+=n
            attempt_index.append(dict(path=str(file),sha256=sha(file),policy=policy,attempts=n,
                accepted_steps=d['steps'],failure=d.get('failure','')))
    require(all(n<=500 for n in attempts.values()),'total package attempt cap exceeded')
    write_json(root/'package-attempt-index.json',attempt_index)
    tests={}
    smoke_label='ctest-readonly-publication' if any(r['label']=='ctest-readonly-publication' for r in manifest['command_receipts']) else 'ctest-readonly-final'
    for row in manifest['command_receipts']:
        if row['label'] not in ('ctest-regressions','ctest-native-transactions',smoke_label):continue
        text=Path(row['log']).read_text();matches=re.findall(r'(\d+)% tests passed, (\d+) tests failed out of (\d+)',text)
        require(row['exit_code']==0 and matches,'selected CTest failure/missing result')
        percent,failed,total=map(int,matches[-1]);require(percent==100 and failed==0,'CTest result is not PASS')
        tests[row['label']]=dict(selected=total,seconds=row['seconds'],log=row['log'],sha256=row['log_sha256'])
    require(sum(x['selected'] for x in tests.values())==29,'expected exactly29 final selected CTests')
    extras={}
    for name in ('run-env.sh','processes-start.txt','processes-before-four-workers.txt','processes-before-main.txt',
        'processes-during-main-and-tests.txt','processes-before-cross-mesh.txt','cross-mesh-allocation.json',
        'parallel-main-returns.json','parallel-main-resumed-returns.json','native-source-equivalence.json','inventory-verification.json','superseded-build-stop.json','build-initial.log','configure-initial.log'):
        file=root/name
        if file.exists():extras[name]=dict(path=str(file),sha256=sha(file))
    receipt=dict(schema='WJ3PublicationV1',tested_source=manifest['source_revision'],tested_tree=manifest['source_tree'],
        specification_sha256=manifest['specification_sha256'],source_content_sha256=manifest['source_content_sha256'],
        manifest_sha256=sha(root/'manifest.json'),artifact_index=manifest['artifacts'],
        tests=tests,qualification_gates={p:{k:s.get(k) for k in ('spatial_gate','scout_gate','main_gate')} for p,s in summaries.items()},
        scheduled_package_attempts=attempts,attempt_cap_per_policy=500,
        package_attempt_index=dict(path=str(root/'package-attempt-index.json'),sha256=sha(root/'package-attempt-index.json')),
        summary_sha256={k:sha(root/'study-final'/k/'summary.json') for k in summaries},extra_receipts=extras,
        table_sha256={p.name:sha(p) for p in sorted(out.glob('*.csv'))},
        finite_scope='serial flat2D p3, one horizontal element, finite spatial/angle inventories, short T; no global positivity/realizability/stability or long-time/3D/MPI/terrain claim')
    write_json(out/'receipt.json',receipt);print(out/'receipt.json')

if __name__=='__main__':main()
