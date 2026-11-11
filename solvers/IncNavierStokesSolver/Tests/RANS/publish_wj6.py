#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Versioned compact WJ-6 native qualification tables and immutable lineage."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import time

from wj_immutable import (Ledger,read,write_once,artifacts,file_identity,file_hash,
                          verify_artifacts,runtime,digest)
from wj_qualification import component,aggregate
from run_cls_epawf_picard import require


def table(path,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('x',newline='') as stream:
        writer=csv.DictWriter(stream,keys);writer.writeheader();writer.writerows(rows)


def leaves(gate,prefix=''):
    for name,value in gate.get('components',{}).items():
        key=prefix+'/'+str(name) if prefix else str(name)
        if 'components' in value:yield from leaves(value,key)
        else:yield key,value


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('root','study','output','build','checks','producer-source','driver-source'):
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();r=a.root.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    start=time.time();ledgers={p:Ledger(r,p,cap=800) for p in ('EP','Budget')}
    label='publication/'+str(out)
    for l in ledgers.values():l.reserve(label,dict(code=artifacts([Path(__file__)])),0,300)
    rows={name:[] for name in ('mesh-runs','fields','residuals','interfaces','divergence','pressure-projection',
        'near-h-loads','independent-loads','sensitivity','accuracy-gates','resources','run-accounting',
        'dense-allocations','arithmetic-estimates')}
    diagnoses={};audits={};stages={};neighborhoods={};summaries={};frozen={};index={}
    try:
        for policy in ledgers:
            summary=a.study/policy/'summary.json';raw=summary.read_bytes()
            snapshot=r/('final-study-'+policy+'-'+file_hash(summary)[:12]+'.json')
            if snapshot.exists():require(snapshot.read_bytes()==raw,'immutable study snapshot differs')
            else:
                with snapshot.open('xb') as stream:stream.write(raw)
            s=json.loads(raw);summaries[policy]=s;frozen[policy]=file_identity(snapshot)
            write_once(out/(policy+'.json'),s)
            for tag,t in s['trajectories'].items():
                for c in [t['initial']]+t['trajectory']:
                    verify_artifacts({c['qualified_record']['path']:c['qualified_record']})
                    for key in ('path','report'):
                        expected=c['sha256'] if key=='path' else c['report_sha256']
                        require(file_hash(c[key])==expected,'native observer/producer changed')
                    require(file_hash(Path(c['checkpoint'])/'State.bin')==c['checkpoint_sha256'],'checkpoint changed')
                    prefix=dict(policy=policy,tag=tag,steps=c['steps'],time=c['time'],dt=c['dt'],unknowns=c['unknowns'],
                        state_sha256=c['state_sha256'],history_sha256=c['history_sha256'],
                        observation_order=c['rule_settings']['order'],horizontal_order=c['rule_settings']['xorder'])
                    rows['mesh-runs'].append(dict(prefix,zknots=json.dumps(c['zknots']),xknots=json.dumps([0.,1.] if tag=='baseline' else [0.,.5,1.]),
                        checkpoint=c['checkpoint'],checkpoint_sha256=c['checkpoint_sha256'],producer=c['report'],producer_sha256=c['report_sha256'],
                        observation=c['path'],observation_sha256=c['sha256'],qualified_record=json.dumps(c['qualified_record'])))
                    for kind in ('continuum_error','initial_drift','initial_representation'):
                        for scope,data in c['errors'][kind].items():
                            if not isinstance(data,dict) or 'fields' not in data:continue
                            for field,v in data['fields'].items():rows['fields'].append(dict(prefix,kind=kind,scope=scope,field=field,
                                **{k:json.dumps(value) if isinstance(value,(dict,list)) else value for k,value in v.items()}))
                    for i,value in enumerate(c['spatial_residual']):
                        rows['residuals'].append(dict(prefix,block=('u','w','k','entire_outer_epsilon','inner_continuity','outer_continuity','gauge','boundary_trace')[i],
                            uncancelled_spatial=value,actual_completed_stage=c['stage_defect'][i] if c['stage_defect'] else None))
                    interface=c['errors']['interface']
                    rows['interfaces'].append(dict(prefix,**{k:v for k,v in interface.items() if not isinstance(v,list)},
                        sampled_K_min=interface['sampled_K'][0],sampled_K_max=interface['sampled_K'][1],
                        cutoff_min=interface['cutoff'][0],cutoff_max=interface['cutoff'][1],
                        strong_divergence=c['errors']['strong_divergence']))
                    for key,v in leaves(c['errors']['accuracy_gate']):
                        rows['accuracy-gates'].append(dict(prefix,quantity=key,**v,margin=v.get('limit',0)-v.get('value',0)))
                    for kind in ('independent_loads','independent_near_h','independent_weak_modes'):
                        if kind in c:rows['independent-loads'].append(dict(prefix,kind=kind,packet=json.dumps(c[kind])))
                    if c.get('arithmetic'):
                        from wj6_qualification import individual_arithmetic
                        individual=individual_arithmetic(c['arithmetic'])
                        for field,value in individual['estimated_comparison_floor'].items():
                            rows['arithmetic-estimates'].append(dict(prefix,field=field,individual_physical_L2=value,
                                individual_sampled_Linf=individual['estimated_sampled_comparison_floor'][field],
                                native_stored_pair_L2=c['arithmetic']['estimated_comparison_floor'][field],
                                safety_factor=100,scope='individual envelope; paired gate adds both trajectories; finite frozen-system estimate'))
                    if c.get('diagnosis'):
                        d=c['diagnosis'];name=policy+'/'+tag+'/'+str(int(c['steps']));diagnoses[name]=d
                        maximum=d['strong_components']['full']['divergence']['maximum']
                        owner=next(e for e in d['element_modes'] if e['geometry']==maximum['geometry'])
                        adjacent=[e for e in d['element_modes'] if e['region']==owner['region'] and
                            e['lower'][0]<=maximum['position'][0]<=e['upper'][0] and
                            e['lower'][2]<=owner['upper'][2]+1e-13 and e['upper'][2]>=owner['lower'][2]-1e-13]
                        neighborhoods[name]=dict(state_sha256=c['state_sha256'],producer=file_identity(c['report']),maximum=maximum,
                            owner=owner,adjacent_elements=adjacent,scope='passive native element/modal view; continuous Q1 complement is not just high element degrees; near-h hat remains [.08,.2]')
                        for scope,values in d['strong_components'].items():
                            for field,v in values.items():rows['divergence'].append(dict(prefix,scope=scope,field=field,L2=v['L2'],sampled_Linf=v['sampled_Linf'],maximum=json.dumps(v['maximum'])))
                        for q in d['pressure_projection']:rows['pressure-projection'].append(dict(prefix,**{k:v for k,v in q.items() if not isinstance(v,list)}))
                        for field,values in zip(('u','w'),d['near_h_signed_loads']):
                            for term,value in zip(d['near_h_terms'],values):rows['near-h-loads'].append(dict(prefix,field=field,term=term,value=value))
                    m=read(c['path'])['measurements'];audits[policy+'/'+tag+'/'+str(int(c['steps']))]=m['audit']
            for phase in ('pilot','primary','extension'):
                for quantity,v in leaves(s['components'].get(phase,{})):
                    if 'uncertainty' in v:
                        rows['sensitivity'].append(dict(policy=policy,phase=phase,quantity=quantity,status=v['status'],value=v['value'],uncertainty=v['uncertainty'],
                            limit=v['limit'],margin=v['limit']-v['value']-v['uncertainty'],fraction=(v['value']+v['uncertainty'])/v['limit']))
            supplement=s.get('order10_row_supplement',{}).get('result')
            if supplement:
                verify_artifacts({supplement['path']:supplement});q=read(supplement['path'])
                if q.get('records'):
                    verify_artifacts({q['records']['path']:q['records']})
                    for key,c in read(q['records']['path']).items():
                        prefix=dict(policy=policy,tag=key.split('/')[0],steps=c['steps'],time=c['time'],dt=c['dt'],unknowns=c['unknowns'],
                            state_sha256=c['state_sha256'],history_sha256=c['history_sha256'],observation_order=10,horizontal_order=6)
                        for kind in ('continuum_error','initial_drift','initial_representation'):
                            for scope,data in c['errors'][kind].items():
                                if not isinstance(data,dict) or 'fields' not in data:continue
                                for field,v in data['fields'].items():rows['fields'].append(dict(prefix,kind=kind,scope=scope,field=field,
                                    **{k:json.dumps(value) if isinstance(value,(dict,list)) else value for k,value in v.items()}))
                    for quantity,v in leaves(q['gate']):
                        if 'uncertainty' in v:rows['sensitivity'].append(dict(policy=policy,phase='all_orders',quantity=quantity,status=v['status'],
                            value=v['value'],uncertainty=v['uncertainty'],limit=v['limit'],margin=v['limit']-v['value']-v['uncertainty'],
                            fraction=(v['value']+v['uncertainty'])/v['limit']))
            for name,e in s['execution'].items():
                rows['resources'].append(dict(policy=policy,case=name,exit_code=e['exit_code'],attempts=e['actual_attempts'],**{k:v for k,v in e['resource'].items() if k!='exit_code' and not isinstance(v,(list,dict))}))
                identity=e.get('report')
                if not identity:continue
                verify_artifacts({identity['path']:identity});d=read(identity['path']);attempts=d['pimple_attempts']
                row=dict(policy=policy,case=name,report=identity['path'],report_sha256=identity['sha256'],absolute_steps=d['steps'],time=d['time'],
                    scheduled_attempts=len(attempts),fully_scheduled=sum(bool(v['scheduled_complete']) for v in attempts),failure=d.get('failure'))
                for key in ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls',
                    'bound_activations','additional_governing_equation_solves','additional_correction_or_restoration_solves'):
                    row[key]=sum(v.get(key,0) for v in attempts)
                equal=True
                for v in attempts:
                    snap=read(v['k_snapshot'])
                    from run_cls_epawf_picard import state_hash
                    equal=equal and state_hash(snap['raw'])==state_hash(snap['corrected'])
                    for field in ('momentum','pressure','k','epsilon'):
                        key=field+'_arithmetic'
                        if key not in v:continue
                        packet=v[key];verify_artifacts({packet['path']:packet});data=read(packet['path'])
                        for kind in ('matrix','mass_or_correction_map'):
                            f=data[kind];identity=file_identity(f['path']);shape=list(map(int,f['shape']))
                            require(identity['sha256']==f['sha256'] and identity['bytes']==8*shape[0]*shape[1],
                                'native dense allocation packet shape/bytes/hash differs')
                            rows['dense-allocations'].append(dict(policy=policy,case=name,time=v['time'],field=field,kind=kind,
                                rows=shape[0],columns=shape[1],bytes=identity['bytes'],sha256=identity['sha256'],path=f['path'],
                                scope='actual exported binary64 payload; selected allocation inventory, not sum of simultaneous live allocations or a full nonlinear Jacobian'))
                row['raw_corrected_k_equal']=equal;rows['run-accounting'].append(row)
                stages[identity['path']]=dict(report=identity,checkpoint=d.get('checkpoint'),last_stage_history=d.get('last_stage_history'),
                    last_stage_time=d.get('last_stage_time'),last_stage_alpha=d.get('last_stage_alpha'),accepted_packets=d.get('wj_audit_packets'),
                    scope='actual producing callbacks; startup alpha0 is not substituted for a completed stage')
            # Only immutable completed case/observer artifacts enter this index.
            # Mutable progress summaries and live publisher/ledger files do not.
            for f in (a.study/policy).rglob('*'):
                if f.is_file() and f.name!='summary.json' and f.suffix!='.lock':index[str(f.resolve())]=file_identity(f)
            passive_receipts=list((a.study/policy).rglob('receipt.json'))
            for root in r.glob('row-qualified-v*'):passive_receipts.extend((root/policy).rglob('receipt.json'))
            for f in passive_receipts:
                if f.parent.name in s['execution']:continue
                d=read(f)
                if 'resource' not in d:continue
                rows['resources'].append(dict(policy=policy,case=str(f.parent.relative_to(r)),
                    exit_code=d['exit_code'],attempts=d['actual_attempts'],receipt=json.dumps(file_identity(f)),
                    **{k:v for k,v in d['resource'].items() if k!='exit_code' and not isinstance(v,(list,dict))}))
        for name,data in rows.items():table(out/(name+'.csv'),data)
        write_once(out/'diagnoses.json',diagnoses);write_once(out/'audit-witnesses.json',audits)
        write_once(out/'maximum-neighborhoods.json',neighborhoods);write_once(out/'accepted-stage-lineage.json',stages)
        for scope in (a.build,a.checks):
            for f in scope.rglob('*'):
                if f.is_file() and 'CMakeFiles' not in f.parts and f.suffix in ('.xml','.json','.bin','.log') and 'Testing' not in f.parts:
                    index[str(f.resolve())]=file_identity(f)
        for scope in [r/'resource-projections-v1',*r.glob('row-qualified-v*')]:
            for f in scope.rglob('*'):
                if f.is_file() and f.suffix!='.lock':index[str(f.resolve())]=file_identity(f)
        index_path=r/('artifact-index-'+out.name+'.json');write_once(index_path,index)
        git=lambda folder,*args:subprocess.check_output(['git','-c','core.filemode=false','-C',str(folder),*args],text=True).strip()
        producer=a.producer_source.resolve();driver=a.driver_source.resolve()
        names=subprocess.check_output(['git','-C',str(producer),'ls-files','-z']).decode().split('\0')
        contents={n:file_hash(producer/n) for n in names if n and (producer/n).is_file()}
        write_once(r/'producer-source-content.json',contents)
        native={n:h for n,h in contents.items() if n.endswith(('.cpp','.h','.hpp')) or n.endswith('CMakeLists.txt')}
        require(all(file_hash(driver/n)==h for n,h in native.items()),'driver/native C++ inputs differ')
        prior=meta(r/'prior-seals.json');verify_artifacts(prior['wj4']);verify_artifacts(prior['wj5'])
        observer_archives={}
        for f in r.glob('observer-*-source.json'):
            packet=meta(f);archive=packet['archive']
            require(file_hash(archive['path'])==archive['sha256'] and Path(archive['path']).stat().st_size==archive['bytes'],'observer source archive changed')
            for name,h in packet['members'].items():require(file_hash(f.with_suffix('')/name)==h,'observer code member changed')
            observer_archives[f.name]=file_identity(f)
        manifest=dict(schema='WJ6LinuxManifestV1',root=str(r),native_source_revision=git(producer,'rev-parse','HEAD'),native_tree=git(producer,'rev-parse','HEAD^{tree}'),
            producer_source_path=str(producer),driver_source_path=str(driver),
            driver_source_revision=git(driver,'rev-parse','HEAD'),driver_tree=git(driver,'rev-parse','HEAD^{tree}'),
            producer_status=git(producer,'status','--short'),driver_status=git(driver,'status','--short'),
            native_input_equivalence=dict(equal=True,files=len(native),content_sha256=digest(native)),
            producer_source_content=file_identity(r/'producer-source-content.json'),
            specifications={n:contents['docs/literature/'+n] for n in ('wj_earsm_neutral_abl_ep_awf_note.tex','wj_earsm_budget_cutoff_awf_appendix.tex')},
            runtime={n:runtime(a.build/n) for n in ('IncNavierStokesSolver','RANSWallResolutionProbe','WJKernelProbe','RANSReturnedFieldsProbe')},
            effective_flags={n:(a.build/'CMakeFiles'/n/'flags.make').read_text() for n in ('IncNavierStokesSolver.dir','RANSWallResolutionProbe.dir','RANSKernels.dir')},
            cmake_cache=file_identity(a.build/'CMakeCache.txt'),test_inventory=file_identity(r/'inventory-checks-v4.json'),
            compiler=subprocess.check_output(['c++','--version'],text=True),platform=subprocess.check_output(['uname','-a'],text=True),
            dependencies=dict(python=__import__('sys').version,numpy=__import__('numpy').__version__,scipy=__import__('scipy').__version__,
                environment=file_identity('/home/jay/.local/opt/nektar-v2026-mpi/env.sh')),
            summaries=frozen,artifact_index=file_identity(index_path),prior_seals=prior,
            arithmetic_policy='Release -ffp-contract=off -fno-fast-math; serial rank; numerical threads1; <=72 all workers including MPI helpers',
            capacity_control='SolverInfo requested/effective4331, default3000; resource-design SHA binds finite process ceilings, not physical identity',
            publication_code=artifacts([Path(__file__)]),observer_source_archives=observer_archives,
            phase_components={p:s['components'] for p,s in summaries.items()})
        write_once(out/'manifest.json',manifest)
        write_once(out/'receipt.json',dict(schema='WJ6CompactReceiptV1',manifest=file_identity(out/'manifest.json'),files=artifacts([f for f in out.iterdir() if f.is_file() and f.suffix!='.lock']),
            elapsed_seconds=time.time()-start,governing_attempts=0,immutable_numeric_snapshots=frozen))
    finally:
        for l in ledgers.values():l.finish(label,0,time.time()-start,str(out/'receipt.json'))
    print(out/'receipt.json')


def meta(path):return json.loads(Path(path).read_text())


if __name__=='__main__':main()
