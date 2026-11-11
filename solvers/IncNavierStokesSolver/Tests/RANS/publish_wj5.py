#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Versioned measurement-only WJ-5 qualification and immutable compact index."""
import argparse
import copy
import csv
import json
import math
import os
from pathlib import Path
import subprocess
import time
from wj_immutable import (read,write_once,file_hash,file_identity,artifacts,verify_artifacts,Ledger,runtime,digest)
from wj4_qualification import sensitivity_sequence,arithmetic_floors
from wj_qualification import aggregate,component
from wj5_measure import boundary_uncertainty,boundary_envelope
from wj4_measure import comparison
from run_cls_epawf_picard import THREADS,require
from publish_wj4 import native_accounting


def records(s):
    result=dict(s.get('retained',{}))
    result.update(s.get('retained_fine',{}))
    t=s.get('adverse_dt_control',{})
    for k in ('initial','terminal'):
        if k in t:result['control/'+k]=t[k]
    for mesh,t in s.get('screens',{}).items():
        for k in ('initial','terminal'):result['screen/'+mesh+'/'+k]=t[k]
    for tag,t in s.get('selected',{}).items():
        result['selected/'+tag+'/0']=t['initial']
        for c in t['trajectory']:result['selected/'+tag+'/'+str(int(c['steps']))]=c
    return result


def requalify(s,design):
    updated=copy.deepcopy(s);rec=records(updated);cache={};lineage=[]
    for key,c in rec.items():
        if 'arithmetic' not in c:continue
        if c['path'] not in cache:
            require(file_hash(c['path'])==c['sha256'],'native observation bytes changed')
            m=read(c['path'])['measurements'];cache[c['path']]=boundary_uncertainty(m,c['arithmetic'],c['rule_settings'])
        c['boundary_uncertainty']=cache[c['path']]
        lineage.append(dict(record=key,observation=file_identity(c['path']),state_sha256=c['state_sha256'],history_sha256=c['history_sha256']))
    # JSON does not preserve aliases between the screen and selected records;
    # rec contains both actual dictionaries, and each is re-estimated above.
    traces=updated.get('selected',{})
    for key,a,b in (('selected_horizontal','Nx1','Nx2'),('selected_timestep','Nx2','Nx2half')):
        if a not in traces or b not in traces:continue
        pair=updated.get(key+'_comparisons',{});ta,tb=traces[a],traces[b];floors=arithmetic_floors({a:ta,b:tb})
        boundary={f:dict(L2=0.,sampled_Linf=0.) for f in floors}
        for t in (ta,tb):
            for c in t['trajectory']:
                if 'boundary_uncertainty' not in c:continue
                for f,row in c['boundary_uncertainty']['fields'].items():
                    for norm in ('L2','sampled_Linf'):boundary[f][norm]=max(boundary[f][norm],row[norm])
        old=sensitivity_sequence(pair,floors,design,[str(i) for i in range(5)])
        floors,boundary,coverage=boundary_envelope({a:ta,b:tb})
        new=sensitivity_sequence(pair,floors,design,[str(i) for i in range(5)],scope_uncertainty={'constrained_boundary':boundary})
        target={}
        for tag,t in ((a,ta),(b,tb)):
            expected=[64,128,192,256] if tag=='Nx2half' else [32,64,96,128]
            rows=t.get('trajectory',[])
            ok=[c['steps'] for c in rows]==expected and not t.get('hard_obstruction') and not t.get('accuracy_obstruction') and t['initial']['errors']['accuracy_gate']['status']=='PASS' and all(c['errors']['accuracy_gate']['status']=='PASS' and abs(c['time']-c['steps']*t['dt'])<=2e-15 for c in rows)
            target[tag]=component('PASS' if ok else 'OBSTRUCTED','all frozen cadence targets and accuracy gates required')
        previous=updated['components'].get(key,component('NOT_RUN'))
        arithmetic=previous.get('components',{}).get('arithmetic',component('NOT_RUN'))
        updated['components'][key]=aggregate(dict(target=aggregate(target),arithmetic=arithmetic,sensitivity=new))
        updated[key+'_uncertainty']=dict(global_floors=floors,constraint_aware_boundary=boundary,boundary_coverage=coverage,original_global_gate=old)
    retained=updated.get('retained',{})
    if all('retained-Nx%d-32'%nx in retained for nx in (1,2)):
        tr={str(nx):dict(initial=retained['retained-Nx%d-0'%nx],trajectory=[retained['retained-Nx%d-32'%nx]]) for nx in (1,2)}
        a,b=(retained['retained-Nx%d-32'%nx] for nx in (1,2));ma=read(a['path'])['measurements'];mb=read(b['path'])['measurements']
        pairs={'Tx':comparison(ma,mb,updated['policy'],a['rule_settings'])};floors,boundary,coverage=boundary_envelope(tr)
        updated['retained_boundary_requalification']=dict(comparison=pairs,original_algorithm=sensitivity_sequence(pairs,floors,design,['Tx']),
            constraint_aware=sensitivity_sequence(pairs,floors,design,['Tx'],scope_uncertainty={'constrained_boundary':boundary}),
            coverage=coverage,measured_divergence_unchanged=[a['errors']['strong_divergence'],b['errors']['strong_divergence']],
            original_WJ4='original immutable publication remains UNRESOLVED for constrained velocity boundary uncertainty and OBSTRUCTED for measured divergence')
    if updated.get('adverse_dt_control'):
        from wj3_measure import checkpoint_payload
        from run_cls_epawf_picard import state_hash
        a=checkpoint_payload(updated['retained']['retained-Nx2-0']['checkpoint']);b=checkpoint_payload(updated['adverse_dt_control']['initial']['checkpoint'])
        updated['adverse_dt_control']['fresh_startup']=dict(exact_primary_gauge=state_hash(a['state'])==state_hash(b['state']),
            exact_physical_histories=[state_hash(v) for v in a['native_history']]==[state_hash(v) for v in b['native_history']],
            native_dt_labels_differ=a['native_history_times']!=b['native_history_times'],
            old_dt=a['dt'],new_dt=b['dt'],old_signature=a['signature'],new_signature=b['signature'])
    updated['components']['measurement_requalification']=component('PASS','versioned constraint-aware estimates applied to every required available observation; original native files unchanged')
    updated['qualification']=aggregate(updated['components'])
    return updated,dict(schema='WJ5ConstraintRequalificationV1',lineage=lineage,code=artifacts([Path(__file__),*[Path(__file__).with_name(n) for n in ('wj5_measure.py','wj4_measure.py','wj4_qualification.py','wj3_measure.py','publish_wj4.py','wj_immutable.py')]]),
        governing_attempts=0,change='dependent epsilon trace uses separate k and epsilon coefficient estimates through their own free rows; pressure uses its own full-gauge row')


def table(path,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('x',newline='') as f:
        w=csv.DictWriter(f,keys);w.writeheader();w.writerows(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--study',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    os.environ.update(THREADS);os.environ['WJ_PACKAGE_ATTEMPT_CAP']='800';r=a.root.resolve();output=a.output.resolve();output.mkdir(parents=True,exist_ok=False)
    start=time.time();ledgers={p:Ledger(r,p,cap=800) for p in ('EP','Budget')};label='publication/'+str(output)
    inputs=artifacts([a.study/p/name for p in ('EP','Budget') for name in ('summary.json','design.json')])
    for l in ledgers.values():l.reserve(label,dict(inputs=inputs,code=artifacts([Path(__file__)])),0,240)
    summaries={};rows={k:[] for k in ('fields','residuals','divergence','pressure-projection','near-h-loads','interfaces','sensitivity','diagnostic-sensitivity','mesh-runs','independent-loads','face-traces','run-accounting')};diagnoses={};histories={};audits={}
    try:
        for policy in ('EP','Budget'):
            s=read(a.study/policy/'summary.json');design=read(a.study/policy/'design.json');s,lineage=requalify(s,design);summaries[policy]=s
            write_once(output/(policy+'-requalification.json'),lineage)
            runs,stages=native_accounting(s,policy);rows['run-accounting'].extend(runs);histories.update(stages)
            for name,c in records(s).items():
                prefix=dict(policy=policy,case=name,mesh=c['mesh'],nx=c['nx'],unknowns=c['unknowns'],steps=c['steps'],time=c['time'],state_sha256=c['state_sha256'],history_sha256=c['history_sha256'])
                rows['mesh-runs'].append(dict(prefix,report=c['report'],report_sha256=c['report_sha256'],checkpoint=c['checkpoint'],checkpoint_sha256=c['checkpoint_sha256'],observation=c['path'],observation_sha256=c['sha256'],zknots=json.dumps(c['zknots'])))
                for kind in ('continuum_error','initial_drift','initial_representation'):
                    for scope,t in c['errors'][kind].items():
                        if not isinstance(t,dict) or 'fields' not in t:continue
                        for field,v in t['fields'].items():rows['fields'].append(dict(prefix,kind=kind,scope=scope,field=field,**{k:v[k] for k in ('L2','sampled_Linf','measure','normalized_L2','normalized_sampled_Linf')}))
                for j,v in enumerate(c['spatial_residual']):rows['residuals'].append(dict(prefix,block=('u','w','k','entire_outer_epsilon','inner_continuity','outer_continuity','gauge','boundary_trace')[j],spatial=v,actual_stage=c['stage_defect'][j] if c['stage_defect'] else None))
                rows['interfaces'].append(dict(prefix,**{k:v for k,v in c['errors']['interface'].items() if not isinstance(v,list)},strong_divergence=c['errors']['strong_divergence']))
                if 'diagnosis' in c:
                    d=c['diagnosis'];diagnoses[policy+'/'+name]=d
                    for scope,t in d['strong_components'].items():
                        for field,v in t.items():rows['divergence'].append(dict(prefix,scope=scope,field=field,L2=v['L2'],sampled_Linf=v['sampled_Linf'],witness=json.dumps(v['maximum'])))
                    for t in d['pressure_projection']:rows['pressure-projection'].append(dict(prefix,**{k:v for k,v in t.items() if not isinstance(v,list)}))
                    for field,load in zip(('u','w'),d['near_h_signed_loads']):
                        for term,v in zip(('storage','conservative_advection','viscous_stress','reduced_pressure','roughness_force'),load):rows['near-h-loads'].append(dict(prefix,field=field,term=term,value=v))
                for kind in ('independent_loads','independent_near_h','independent_weak_modes'):
                    if kind in c:rows['independent-loads'].append(dict(prefix,kind=kind,packet=json.dumps(c[kind])))
                # All owners remain in immutable heavy packets. Compact face
                # traces retain every actual normal face at five shared rays.
                m=read(c['path'])['measurements']
                audits[policy+'/'+name]=m['audit']
                for obs in m['observations']:
                    x,_,z=obs['position']
                    if obs['weight'] or x not in (0.,.25,.5,.75,1.) or z not in c['zknots'] or not any(obs['face_sides_xz']):continue
                    rows['face-traces'].append(dict(prefix,x=x,z=z,geometry=obs['geometry'],element=obs['element'],H_side=obs['side'],x_side=obs['requested_x_side'],normal_side=obs['requested_normal_side'],local_faces=json.dumps(obs['face_sides_xz']),Gxx=obs['raw_G'][0],Gzz=obs['raw_G'][8],divergence=obs['raw_G'][0]+obs['raw_G'][8],u=obs['value'][0],w=obs['value'][1],pi=obs['value'][4]))
            for group in ('selected_horizontal','selected_timestep'):
                gate=s['components'].get(group,{})
                for timekey,t in gate.get('components',{}).get('sensitivity',{}).get('components',{}).items():
                    for key,v in t.get('components',{}).items():
                        if 'value' in v:rows['sensitivity'].append(dict(policy=policy,group=group,cadence_index=timekey,quantity=key,status=v['status'],value=v['value'],uncertainty=v['uncertainty'],limit=v['limit'],margin=v['limit']-v['value']-v['uncertainty'],fraction=(v['value']+v['uncertainty'])/v['limit']))
            for estimator,gate in s.get('retained_boundary_requalification',{}).items():
                if not isinstance(gate,dict) or 'components' not in gate:continue
                for t in gate['components'].values():
                    for key,v in t.get('components',{}).items():
                        if 'value' in v:rows['diagnostic-sensitivity'].append(dict(policy=policy,group='retained_hx_Tx',estimator=estimator,quantity=key,status=v['status'],value=v['value'],uncertainty=v['uncertainty'],limit=v['limit'],margin=v['limit']-v['value']-v['uncertainty']))
            write_once(output/(policy+'.json'),dict(qualification=s['qualification'],components=s['components'],retained=s.get('retained'),retained_fine=s.get('retained_fine'),retained_boundary_requalification=s.get('retained_boundary_requalification'),screens=s.get('screens'),adverse_dt_control=s.get('adverse_dt_control'),
                selected=s.get('selected'),comparisons={k:v for k,v in s.items() if k.endswith('_comparisons') or k.endswith('_uncertainty')},
                frozen_design=design,heavy_summary=file_identity(a.study/policy/'summary.json'),stationarity=s['stationarity']))
        for key,value in rows.items():table(output/(key+'.csv'),value)
        write_once(output/'diagnoses.json',diagnoses)
        write_once(output/'audit-witnesses.json',audits)
        write_once(output/'accepted-stage-lineage.json',histories)
        # Preserve exact heavy input/lineage bytes. Every source/run artifact is
        # indexed once; build object files and repository copies are not fields.
        index={}
        for scope in [a.study,r/'study',r/'build/wj2-native-EP',r/'build/wj2-native-Budget',r/'build/wj2-legacy-transaction',r/'build/wj3-smoke',r/'build/wj4-observations',r/'build/wj5-observations']:
            if scope.exists():
                for f in scope.rglob('*'):
                    if f.is_file():index[str(f.resolve())]=file_identity(f)
        for f in r.iterdir():
            if f.is_file() and f.suffix in ('.json','.log','.sh'):index[str(f.resolve())]=file_identity(f)
        for f in (r/'build').rglob('*'):
            if f.is_file() and 'CMakeFiles' not in f.parts and f.suffix in ('.json','.xml','.bin','.fld','.opt','.log'):
                index[str(f.resolve())]=file_identity(f)
        index_path=r/('wj5-artifact-index-'+output.name+'.json');write_once(index_path,index)
        source=r/'source';git=lambda folder,*args:subprocess.check_output(['git','-c','core.filemode=false','-C',str(folder),*args],text=True).strip()
        content={name:file_hash(source/name) for name in subprocess.check_output(['git','-C',str(source),'ls-files','-z']).decode().split('\0') if name and (source/name).is_file()}
        write_once(r/'wj5-source-content.json',content)
        metadata_source=Path(__file__).resolve().parents[4]
        native_names=[name for name in content if name.endswith(('.cpp','.h','.hpp')) or name.endswith('CMakeLists.txt')]
        equivalent=all((metadata_source/name).is_file() and file_hash(metadata_source/name)==content[name] for name in native_names)
        require(equivalent,'publication native inputs differ from tested source')
        manifest=dict(schema='WJ5LinuxManifestV1',root=str(r),study=str(a.study.resolve()),producer_source_revision=git(source,'rev-parse','HEAD'),producer_tree=git(source,'rev-parse','HEAD^{tree}'),producer_status=git(source,'status','--short'),
            publication_source_revision=git(metadata_source,'rev-parse','HEAD'),publication_tree=git(metadata_source,'rev-parse','HEAD^{tree}'),
            native_build_source_revision=git(source,'rev-parse','da9518a'),
            native_source_equivalence=dict(equal=equivalent,files=len(native_names),content_sha256=digest({n:content[n] for n in native_names})),source_content=file_identity(r/'wj5-source-content.json'),
            specification_sha256={n:content['docs/literature/'+n] for n in ('wj_earsm_neutral_abl_ep_awf_note.tex','wj_earsm_budget_cutoff_awf_appendix.tex')},
            runtime={name:runtime(r/'build'/name) for name in ('IncNavierStokesSolver','RANSWallResolutionProbe','WJKernelProbe','RANSReturnedFieldsProbe')},
            effective_flags={name:(r/'build/CMakeFiles'/name/'flags.make').read_text() for name in ('IncNavierStokesSolver.dir','RANSWallResolutionProbe.dir','RANSKernels.dir')},
            cmake_cache=file_identity(r/'build/CMakeCache.txt'),inventory=file_identity(r/'inventory.log'),command_receipts=file_identity(r/'command-receipts.json'),
            compiler=subprocess.check_output(['c++','--version'],text=True),platform=subprocess.check_output(['uname','-a'],text=True),
            dependencies=dict(python=__import__('sys').version,numpy=__import__('numpy').__version__,scipy=__import__('scipy').__version__,sympy=__import__('sympy').__version__,environment=file_identity('/home/jay/.local/opt/nektar-v2026-mpi/env.sh')),
            arithmetic_policy='Release; -ffp-contract=off -fno-fast-math; serial; numerical threads1; <=72 total server workers',threads=THREADS,
            ledger_snapshot={p:dict(identity=file_identity(l.file),totals=l.totals()) for p,l in ledgers.items()},heavy_artifact_index=file_identity(index_path),
            observation_contract='actual compared x/z meshes and physical breakpoints; own subdomain measures; one full-domain pressure gauge; zero-weight owner duplicates; actual accepted callback only',
            immutable_inputs=inputs,configuration={p:read(a.study/p/'design.json') for p in ('EP','Budget')},selection=read(a.study/'selection.json') if (a.study/'selection.json').exists() else None)
        write_once(output/'manifest.json',manifest)
        verify_artifacts(inputs)
        write_once(output/'receipt.json',dict(schema='WJ5CompactReceiptV1',manifest=file_identity(output/'manifest.json'),files=artifacts([f for f in output.iterdir() if f.is_file()]),
            governing_attempts=0,measured_seconds=time.time()-start,original_wj4_receipts_unchanged=True))
    finally:
        for l in ledgers.values():l.finish(label,0,time.time()-start,str(output/'receipt.json'))
    print(output/'receipt.json')

if __name__=='__main__':main()
