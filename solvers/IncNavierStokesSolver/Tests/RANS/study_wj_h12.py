#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Exactly fourteen immutable attempts, four owners and one sparse science batch."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import time
import numpy as np
from run_wj_h4 import load,save,sha
from run_wj_h7 import update_journal
from run_wj_h10 import runtime_snapshot
from wj_h12_phase import (ATTEMPTS,BRANCHES,WINDOWS,OWNER_KEYS,OBSERVABLES,PHASES,SparseRows,checked_context,
    bind,local,decompositions,window,paired_window,identity,scale,comparison_scale,kernel_input,sparse_flow,predictor_eligible,DIRECTION,role_context)
from wj_h4_points import vector_hash
from wj_h9_local import same
from wj_h3_analysis import principal
from validate_wj1 import Probe

H10=Path('/home/jay/build/wj-h10-flat-inlet-20261113')
H11=Path('/home/jay/build/wj-h11-compatible-common-time-20261113')
KERNEL_SHA='3d4edf79a2718a281175be1a575776a33dd5f988c48d5a1b594a9b5f2dd7844c'

def serial(v):
    if isinstance(v,np.ndarray):return serial(v.tolist())
    if isinstance(v,np.generic):return serial(v.item())
    if isinstance(v,float)and not np.isfinite(v):return None
    if isinstance(v,dict):return {k:serial(x)for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [serial(x)for x in v]
    return v

def file_id(p):
    p=Path(p);return dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p))

def check(v):
    actual=file_id(v['path'])
    if actual!=v:raise ValueError('immutable file identity differs: '+v['path'])
    return actual

def manifest_item(m,p):
    rows=[v for v in m['heavy']if v['path']==str(p)]
    if len(rows)!=1:raise ValueError('required input absent/duplicated in sealed manifest: '+str(p))
    return check(rows[0])

def stage(root):
    out=root/'science';out.mkdir(exist_ok=False);verified=[];manifests={}
    frozen=root/'source/docs/plans'
    for label,base in (('H10',H10),('H11',H11)):
        local_manifest=frozen/('wj-h10-linux'if label=='H10'else 'wj-h11-linux')/'publication-v1/manifest.json'
        remote_manifest=base/'publication-v1/manifest.json'
        assert sha(local_manifest)==sha(remote_manifest),'published parent manifest differs'
        manifests[label]=load(remote_manifest);verified.append(file_id(remote_manifest))
        for name in ('publication-v1/seal.json',):
            local_path=local_manifest.parent.parent/name;assert sha(local_path)==sha(base/name)
            verified.append(file_id(base/name))
    bindings_path=H11/'publication-v1/source-state-bindings.json'
    assert sha(bindings_path)==sha(frozen/'wj-h11-linux/publication-v1/source-state-bindings.json')
    binding=load(bindings_path);verified.append(file_id(bindings_path));reports={};attempts=[]
    for b,dt in BRANCHES:
        rec=next(v for v in binding if v['branch']==b)
        reports[b]={}
        for label,producer in (('H10',rec['parent']['producer']),('H11',rec['producer'])):
            manifest_item(manifests[label],producer['path']);check(producer);reports[b][label]=producer
        old=load(rec['parent']['producer']['path']);new=load(rec['producer']['path'])
        assert old['dt']==new['dt']==dt and not old['failure']and not new['failure']
        assert len(old['state'])==len(new['state'])==2328
        for key in ('fixed_reference','fixed_defect'):assert same(old[key],new[key])
        assert same(checked_context(old,1)['state'],checked_context(new,1)['state'])
        # Never use the H11 restored step-1 alpha=0 as its completed callback.
        assert checked_context(new,1)['alpha']==0 and checked_context(old,1)['alpha']==1/dt
        for w,bb,d,step in ATTEMPTS:
            if bb!=b:continue
            label='H10'if step==1 else 'H11';r=old if step==1 else new
            from wj_h11_analysis import selected
            cx,a,index=selected(r,step);checked_context(r,step-1)
            packet=manifest_item(manifests[label],a['k_snapshot'])
            assert packet['sha256']==a['k_snapshot_sha256']
            attempts.append(dict(window=w,branch=b,dt=dt,step=step,producer_label=label,
                                 producer=reports[b][label],packet=packet,local_attempt=index))
    # Original owner metadata, not coordinates regenerated from rounded mm.
    local_id=manifest_item(manifests['H11'],H11/'science/common-v1/fixed-owner-history.json')
    records=load(local_id['path']);owner_records=[]
    for geo,side,z in OWNER_KEYS:
        candidates=[v['record']for v in records if v['branch']=='coarse'and v['step']==1 and
                    v['record']['owner']['geometry']==geo and v['record']['owner']['normal_side']==side]
        assert len(candidates)==1 and candidates[0]['physical']==[8.,0.,z]
        owner_records.append({k:candidates[0][k]for k in ('physical','owner')})
    anchors=[]
    for v in records:
        if v['step']in (round(.01/dict(BRANCHES)[v['branch']]),round(.04/dict(BRANCHES)[v['branch']])):
            anchors.append(v)
    assert len(anchors)==24
    arithmetic={}
    ar=load(H11/'publication-v1/arithmetic-bindings.json')
    for b,dt in BRANCHES:
        p=H11/f'science/arithmetic-v1/{b}.json';v=manifest_item(manifests['H11'],p)
        arithmetic[b]=v
    build=load(H10/'science/runtime-build.json');observed=load(H10/'closure/runtime-observed.json')
    compiled=file_id(H10/'closure/compiled-inputs.json')
    assert compiled['sha256']==sha(frozen/'wj-h10-linux/closure/compiled-inputs.json')
    kernel=file_id(H10/'runtime/WJKernelProbe.elf');wrapper=file_id(H10/'runtime/WJKernelProbe')
    assert kernel['sha256']==KERNEL_SHA
    assert build['kernel']['sha256']in (wrapper['sha256'],kernel['sha256'])
    runtime_files=[]
    for v in observed['libraries']:
        if Path(v['path']).name in ('IncNavierStokesSolver.elf','RANSWallResolutionProbe.elf'):continue
        runtime_files.append(check(v))
    if kernel not in runtime_files:runtime_files.append(kernel)
    # Do not inspect/rebuild unrelated dependency roots. Recorded H10 compiler
    # and source identities are retained as compiled provenance, not new HEAD.
    runtime=dict(ELF=kernel,wrapper=wrapper,compiled_inputs=compiled,build=file_id(H10/'science/runtime-build.json'),
                 original_source=load(compiled['path']),libraries=runtime_files,observed=file_id(H10/'closure/runtime-observed.json'))
    runtime['selected_compiled_sources']=[]
    for v in build['compiled_inputs']:
        if 'RANSKernels.dir'in v['command']or v['path'].endswith('/WJKernelProbe.cpp'):
            assert sha(v['path'])==v['sha256']and '-ffp-contract=off'in v['command']and '-fno-fast-math'in v['command']
            runtime['selected_compiled_sources'].append(file_id(v['path']))
    assert any(v['path'].endswith('/WJEARSM.cpp')for v in runtime['selected_compiled_sources'])
    plan=dict(schema='WJH12SavedPhaseIntentV1',attempts=sorted(attempts,key=lambda v:(v['window'],[x[0]for x in BRANCHES].index(v['branch']),v['step'])),
              reports=reports,owners=owner_records,arithmetic=arithmetic,endpoint_records=local_id,
              accepted_endpoint_anchors=anchors,runtime=runtime,verified_metadata=verified,
              parent_manifests={k:file_id(base/'publication-v1/manifest.json')for k,base in (('H10',H10),('H11',H11))},
              aliases=[],scope=dict(attempts=14,raw_phase_owner_roles=224,new_native=0,new_strict=0,new_builds=0,
                                   pressure_mass_factorizations=0,pressure_mass_solves=0,kernel_requests_max=48))
    save(out/'input-intent.json',serial(plan))
    update_journal(root/'accounting.json',lambda b:b.setdefault('intent_sha256',sha(out/'input-intent.json')))
    print('H12 selected imports PASS: six producers, fourteen snapshots, four owners, three sealed arithmetic sets',flush=True)

def maps_for_kernel(probe,runtime):
    allowed={v['path']:v for v in runtime['libraries']}
    limit=time.monotonic()+5
    while time.monotonic()<limit:
        paths=set()
        for line in Path(f'/proc/{probe.p.pid}/maps').read_text().splitlines():
            f=line.split(maxsplit=5)
            if len(f)==6 and f[5].startswith('/'):paths.add(str(Path(f[5]).resolve()))
        sample=runtime_snapshot(H10,sorted(paths))
        if sample:
            if sample['status']!='PASS':raise ValueError('mixed native-exec private libc')
            for p in paths:
                if '.so'in Path(p).name or p.endswith('.elf'):
                    if p not in allowed:raise ValueError('unqualified native runtime mapping: '+p)
                    check(allowed[p])
            return dict(sample,pid=probe.p.pid,ELF=runtime['ELF'],wrapper=runtime['wrapper'],
                        mapped_identities=[file_id(p)for p in sorted(paths)if Path(p).is_file()])
        time.sleep(.01)
    raise ValueError('no qualified native-exec snapshot')

def science(root):
    start=time.monotonic();intent=root/'science/input-intent.json';assert sha(intent)==load(root/'accounting.json')['intent_sha256']
    imp=load(intent);out=root/'science/results-v1';out.mkdir(exist_ok=False)
    reports={b:{label:load(check(v)['path'])for label,v in rr.items()}for b,rr in imp['reports'].items()}
    caches={b:SparseRows(reports[b]['H10'],imp['owners'])for b,dt in BRANCHES}
    # Native mesh/maps must match before sharing immutable basis data.
    for b in reports:
        assert [f['mesh']for f in reports[b]['H10']['fields']]==[f['mesh']for f in reports[b]['H11']['fields']]
        assert reports[b]['H10']['offsets']==reports[b]['H11']['offsets']
        assert reports[b]['H10']['geometry_signature']==reports[b]['H11']['geometry_signature']
    local_cache={};role_records=[];bindings=[];transition=[];phase_values={};state_arrays={};source_by_attempt={};pilot=None
    def evaluate_state(state):
        key=vector_hash(state)
        if key not in local_cache:
            local_cache[key]=[local(v)for v in caches['coarse'].linear(state)]
        return key,local_cache[key]
    for number,rec in enumerate(imp['attempts']):
        t0=time.monotonic();r=reports[rec['branch']][rec['producer_label']];packet=load(check(rec['packet'])['path'])
        states,binding=bind(r,packet,rec['step'],rec['producer'],imp['reports'][rec['branch']][rec['producer_label']])
        binding.update(branch=rec['branch'],window=rec['window'],snapshot=rec['packet'],source=rec['producer'])
        bindings.append(binding);key=(rec['branch'],rec['step']);source_by_attempt[key]=rec;state_arrays[key]=states
        phase_values[key]={}
        for phase in PHASES:
            h,values=evaluate_state(states[phase]);phase_values[key][phase]=values
            for i,v in enumerate(values):
                role_records.append(dict(branch=rec['branch'],step=rec['step'],window=rec['window'],phase=phase,owner_index=i,
                    **role_context(binding,phase),state_sha256=h,alias_key=h+':'+str(i),physical=imp['owners'][i]['physical'],owner=imp['owners'][i]['owner'],value=v))
        for i in range(4):
            for pa,pb in zip(PHASES[:-1],PHASES[1:]):
                transition.append(dict(branch=rec['branch'],step=rec['step'],owner_index=i,before=pa,after=pb,
                                       decomposition=decompositions(phase_values[key][pa][i],phase_values[key][pb][i])))
        if number==0:
            error=caches['coarse'].reference_check(states['A']);assert error<2e-11
            elapsed=time.monotonic()-start;cost=time.monotonic()-t0
            remaining=1.25*(13*cost+25.) # 25 s for bounded kernel/aggregation/publication preparation.
            book=load(root/'accounting.json');prior=sum(v.get('seconds',v['reserved_seconds'])for v in book['entries'].values()if v['phase']=='operation')
            admitted=elapsed+remaining<=360 and 30+prior+elapsed+remaining+90<=600
            pilot=dict(status='PASS'if admitted else 'RESOURCE_LIMIT',attempt=rec,elapsed_setup_and_pilot=elapsed,
                       first_attempt_seconds=cost,remaining_projection_with_25_percent=remaining,untouched_closure_seconds=90,
                       resource_total_projection=30+prior+elapsed+remaining+90,scalar_reference_check_error=error)
            save(out/'retained-pilot.json',serial(pilot));save(out/'pilot-roles.json',serial(role_records));save(out/'pilot-binding.json',serial(binding))
            if not admitted:raise RuntimeError('RESOURCE_LIMIT measured H12 remaining work')
    assert len(bindings)==14 and len([v for v in role_records if v['phase']!='P'])==224
    save(out/'phase-bindings.json',serial(bindings));save(out/'local-phase-roles.json',serial(role_records))
    save(out/'phase-decompositions.json',serial(transition))
    # Accepted endpoints are independently compared to all 24 published H11 owners.
    replay=[]
    for old in imp['accepted_endpoint_anchors']:
        b,s=old['branch'],old['step'];rec=old['record'];i=next(i for i,o in enumerate(imp['owners'])if o['owner']==rec['owner']and o['physical']==rec['physical'])
        new=phase_values[b,s]['A'][i]
        if new['constitutive_status']!='DEFINED':
            replay.append(dict(branch=b,step=s,time=old['time'],owner_index=i,status='NOT_EVALUABLE',reason=new['reason']));continue
        errors={name:float(np.max(abs(np.asarray(new[newkey])-np.asarray(rec[oldkey]))))for name,newkey,oldkey in
                (('U','U','U'),('G','G','G'),('k','k','k'),('epsilon','epsilon','epsilon'),('R','R','R'),('pi','pi','pi'))}
        errors['mhat']=abs(new['mhat']-rec['principal']['in_plane'])
        independent=principal(new['k'],new['epsilon'],new['G']);errors['independent_mhat']=abs(independent['in_plane']-new['mhat'])
        floors={name:(2e-11 if name in ('mhat','independent_mhat')else 4e-11*comparison_scale(name))for name in errors}
        assert all(v<floors[name]for name,v in errors.items()),'H11 endpoint replay differs'
        replay.append(dict(branch=b,step=s,time=old['time'],owner_index=i,errors=errors,mhat=new['mhat'],a=new['a'],g=new['g']))
    save(out/'endpoint-replay.json',serial(dict(status='PASS',records=replay)))
    # Full equal-window cumulative contributions retain the incoming .03 term.
    cumul=[];cross=[];pairs=[]
    for w in range(2):
        ws={}
        for b,dt in BRANCHES:
            steps=[s for ww,bb,d,s in ATTEMPTS if (ww,bb)==(w,b)]
            for i in range(4):
                rows=[{p:phase_values[b,s][p][i]for p in PHASES}for s in steps]
                for a,z in zip(rows[:-1],rows[1:]):
                    for name in OBSERVABLES:
                        if name in a['A']and name in z['P']:assert np.array_equal(a['A'][name],z['P'][name])
                for name in OBSERVABLES:
                    v=window(rows,name);ws[b,i,name]=v;cumul.append(dict(window=w,branch=b,owner_index=i,observable=name,**v))
        for a,b in (('coarse','fine'),('fine','quarter')):
            for i in range(4):
                for name in OBSERVABLES:
                    v=paired_window(ws[a,i,name],ws[b,i,name],name)
                    cross.append(dict(window=w,pair=a+'-'+b,owner_index=i,observable=name,**v))
                sa=phase_values[a,round(WINDOWS[w][1]/dict(BRANCHES)[a])]['A'][i]
                sb=phase_values[b,round(WINDOWS[w][1]/dict(BRANCHES)[b])]['A'][i]
                pairs.append(dict(window=w,pair=a+'-'+b,owner_index=i,convention='coarser minus finer',decomposition=decompositions(sb,sa)))
    save(out/'window-cumulative.json',serial(cumul));save(out/'cross-window-accounting.json',serial(cross));save(out/'endpoint-decompositions.json',serial(pairs))
    predictors=[]
    for b,dt in BRANCHES:
        s=round(.01/dt);rec=source_by_attempt[b,s];packet=load(check(imp['arithmetic'][b])['path'])
        assert packet['producer']==rec['producer']and packet['absolute_step']==s and packet['alpha']==1/dt
        v=sparse_flow(caches[b],state_arrays[b,s],packet,b,s)
        predictors.append(dict(branch=b,step=s,producer=rec['producer'],arithmetic=imp['arithmetic'][b],**v))
    save(out/'predictor-coverage.json',serial(dict(status='PASS',covered_attempts=3,early_attempts=7,late_covered=0,
        limitation='only attempts ending .01, with unequal incoming durations; no late predictor/correction split',records=predictors)))
    # Verify closure identity and actual loader before any limited requests.
    runtime=imp['runtime'];check(runtime['ELF']);check(runtime['wrapper']);assert runtime['ELF']['sha256']==KERNEL_SHA
    probe=Probe(runtime['wrapper']['path']);compiled=[];aliases=[];requested={};skip=[]
    try:
        maps=maps_for_kernel(probe,runtime);save(out/'kernel-native-exec.json',serial(maps))
        for w,(lo,hi)in enumerate(WINDOWS):
            for b,dt in BRANCHES:
                s=round(hi/dt)
                for phase in ('C','F','K','A'):
                    v=phase_values[b,s][phase][0];role=dict(branch=b,step=s,window=w,phase=phase,owner_index=0)
                    if v['constitutive_status']!='DEFINED':skip.append(dict(role,reason=v['reason']));continue
                    h=kernel_input(v);aliases.append(dict(role,input_sha256=h))
                    if h in requested:continue
                    value=probe.closure(v['k'],v['epsilon'],v['G']);derivative=probe.closure(v['k'],v['epsilon'],v['G'],(0.,0.,DIRECTION))
                    d=-float(np.sum(DIRECTION*derivative['R']));normalized=d/value['scalar'][1]
                    errs=dict(dimensional=abs(d-v['d']),implicit_dN=abs(derivative['scalar'][12]-v['implicit_dN']),normalized=abs(normalized-v['a']),
                              R=float(np.max(abs(value['R']-v['R']))),N=abs(value['scalar'][12]-v['N']))
                    assert errs['dimensional']<2e-13 and errs['implicit_dN']<2e-11 and errs['normalized']<2e-11
                    requested[h]=len(compiled);compiled.append(dict(input_sha256=h,first_role=role,errors=errs,value={k:x.tolist()for k,x in value.items()},
                        derivative={k:x.tolist()for k,x in derivative.items()}))
        assert probe.calls<=48
    finally:probe.close()
    save(out/'compiled-checks.json',serial(dict(status='PASS'if not skip else 'PARTIAL',requests=probe.calls,raw_roles=24,aliases=aliases,skipped=skip,records=compiled)))
    status='PASS'if all(v['value']['constitutive_status']=='DEFINED'for v in role_records)and not skip else 'PARTIAL'
    result=dict(status=status,input_intent=file_id(intent),source=file_id(__file__),attempts=14,raw_phase_owner_roles=224,predecessor_roles=56,
                unique_physical_state_evaluations=len(local_cache),sparse_owner_evaluations=4*len(local_cache),kernel_requests=probe.calls,
                kernel_raw_roles=24,kernel_unique_inputs=len(compiled),saved_velocity_vector_scatters=9,
                diagnostic_mass_factorizations=0,diagnostic_mass_solves=0,governing_matrix_factorizations=0,governing_matrix_solves=0,
                new_governing_attempts=0,new_CP0=0,new_strict_readbacks=0,new_builds=0,seconds=time.monotonic()-start,
                H11_historical_temporal='UNRESOLVED',scope='sparse actual saved phases only; mathematical telescopes are not unique causal attribution')
    save(out/'result.json',serial(result));print('H12 science',status,result,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('stage','science'),required=True);a=p.parse_args()
    try:globals()[a.mode](a.root)
    except Exception as err:
        update_journal(a.root/'accounting.json',lambda b:b.setdefault('h12_halt',dict(status='UNRESOLVED',operation=a.mode,reason=repr(err))))
        raise
