#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Verify H11 semantic producer roles and close its existing independent ledger."""
import argparse
import json
from pathlib import Path
import re
from run_wj_h11 import H10,load,save,size
from study_wj_h9 import identity,check
from study_wj_h11 import parent,audit
from wj_h11_analysis import MATRIX,VIEWS,selected
from wj_h9_local import same,phase_bindings
from wj_h4_points import vector_hash
from wj3_measure import checkpoint_payload


def verify(root):
    out=root/'closure';out.mkdir();pub=root/'publication-v1';verified=[]
    for directory in (pub,root/'local-comparisons-v1'):
        for v in load(directory/'seal.json')['files']:verified.append(check(v))
    for v in load(pub/'manifest.json')['heavy']:check(v)
    gates=load(pub/'component-gates.json');assert gates['verification']=='PASS'and gates['status']in('PASS','UNRESOLVED')
    assert(gates['new_governing_reserved'],gates['new_governing_attempted'],gates['new_governing_accepted'],gates['new_strict'],gates['new_CP0'],gates['new_builds'])==(25,25,25,11,0,0)
    bindings=[];count=0;runtime=[];libraries={};methods=set()
    for branch,dt,n,target in MATRIX:
        case=parent(root,branch);r=load(case['report']);old=load(case['parent']['producer']['path']);assert r['steps']==target and not r['failure']
        assert r['initialization_lift_calls_iterations']==[0,0]and r['initial_velocity_calls']==0
        assert same(r['accepted_states'][0],old['state'])and r['signature']==old['signature']
        for key in ('fixed_reference','fixed_defect'):assert same(r[key],old[key])
        starting=checkpoint_payload(Path(case['parent']['checkpoint']['path']).parent,include_reference=True)
        assert same(starting['native_history'],old['native_history'])and same(starting['native_history_times'],old['native_history_times'])
        phases=[]
        for step in range(2,target+1):
            c,a,i=selected(r,step);assert a['scheduled_complete']and a['governing_solves_complete']and not a['failure']
            assert all(v['passed']for v in c['packet']['conditions'].values())and not c['packet']['failure_categories']
            p=load(a['k_snapshot']);assert __import__('run_wj_h11').sha(a['k_snapshot'])==a['k_snapshot_sha256'];states,b=phase_bindings(r,p,step)
            assert b['P_C_bitwise']and b['post_flow_pre_k_bitwise']and b['raw_corrected_bitwise']and b['flow_A_U_pi_bitwise']
            if step==2:assert same(p['history'],old['state']),'first actual physical callback history differs from imported state'
            phases.append(b);count+=1
        index=load(audit(root,branch)['index']['path'])
        for step in VIEWS[branch]:
            folder=root/f'science/views/{branch}-{step}';view=load(folder/'result.json');assert view['status']=='PASS';cx=view['context'];check(cx['report']);check(cx['checkpoint']);check(view['native'])
            c,a,i=selected(r,step);assert cx['alpha']==1/dt and cx['state_sha256']==vector_hash(c['state'])and cx['history_sha256']==vector_hash(c['history'])
            p=Path(view['native']['path'])
            with p.open('rb')as f:
                prefix=f.read(8192).decode();f.seek(max(0,p.stat().st_size-4096));tail=f.read().decode()
            marker=',"open_boundary":';assert marker in prefix;header=json.loads(prefix.partition(marker)[0]+'}')
            assert header['operation']=='open_checks'and header['selected_absolute_step']==step and header['producer_absolute_step']==target and not header['producer_failure']
            assert header['report_sha256']==cx['report']['sha256']and header['checkpoint_sha256']==cx['checkpoint']['sha256']
            assert header['state_sha256']==cx['state_sha256']and header['callback_history_sha256']==cx['history_sha256']
            cp=next(v for v in index if v['checkpoint']==cx['checkpoint']);assert header['checkpoint_native_history_sha256']==cp['native_history']and same(header['checkpoint_native_history_times'],cp['native_history_times'])
            assert '"initial_velocity_calls":0'in tail and '"governing_solves":0'in tail and '"initial_velocity_mode":"FlatInletTraceV1"'in tail
            assert re.search(r'"initialization_work":\[0,0\]',tail)
            assert all(max(map(abs,v))<=2e-15 for v in view['native_Q1_moments'].values())
            assert view['anchor']['status']==view['quadrature']['status']=='PASS';bindings.append(dict(branch=branch,step=step,context=cx,native_header=header,strict_native=check(view['native']),zero_work_verified=True))
        for label in (branch+'-block',*[f'{branch}-{s}-strict'for s in VIEWS[branch]]):
            p=root/'operations'/label;samples=load(p/'native-exec-samples.json');assert samples and all(v['status']=='PASS'for v in samples)
            observed={v['path']:v for v in load(p/'loaded-maps.json')};files={name for v in samples for name in v['files']if '.so'in Path(name).name or name.endswith('.elf')}
            for name in files:
                if name not in libraries:libraries[name]=check(observed[name])
            runtime.append(dict(operation=label,sample_count=len(samples),first=samples[0],last=samples[-1],full_samples=identity(p/'native-exec-samples.json')))
        save(out/(branch+'-phase-bindings.json'),phases)
    assert count==25 and len(bindings)==11
    for label in ('focused-tests','focused-h11-final'):
        text=(root/'operations'/label/'command.log').read_text();assert '\nOK\n'in text;found=re.findall(r'^(test_\w+) \(([^)]+)\)',text,re.M);assert len(found)==int(re.search(r'Ran (\d+) tests',text).group(1));methods.update(found)
    assert len(methods)==74
    arithmetic=load(root/'science/arithmetic-v1/result.json');assert arithmetic['status']=='PASS'and arithmetic['new_attempt_level_arithmetic_sets']==2
    for b in ('coarse','fine','quarter'):
        packet=load(arithmetic['branches'][b]['packet']['path']);assert packet['time']==.01 and all(v['status']=='PASS'for v in packet['checks'].values())
    common=load(root/'science/common-v1/result.json');assert common['roles']==31 and common['complete']and common['status']=='PASS'
    book=load(root/'accounting.json');assert book['reserved_attempts']==25 and book['governing_blocks']==3 and book['strict_readbacks']==11 and not book.get('h11_halt')and not book.get('coverage_halt')
    mass_factors=common['diagnostic_mass_factorizations'];mass_solves=common['diagnostic_mass_solves'];queries=common['kernel_queries']
    for branch,dt,n,target in MATRIX:
        a=audit(root,branch);mass_factors+=a['diagnostic_mass_factorizations'];mass_solves+=a['diagnostic_mass_solves'];queries+=a['kernel_queries']
    for branch in ('fine','quarter'):
        a=load(root/f'science/arithmetic-v1/{branch}.json');mass_factors+=a['diagnostic_mass_factorizations'];mass_solves+=a['diagnostic_mass_solves']
    save(out/'runtime-observed.json',dict(status='PASS',actual_exec=runtime,libraries=list(libraries.values()),qualified_H10_source=identity(H10/'closure/compiled-inputs.json'),scope='actual native-exec snapshots; pre-exec PID-union host mappings are separate'))
    save(out/'verification.json',dict(status='PASS',scientific_temporal=common['temporal'],original_component_gates=identity(pub/'component-gates.json'),verified_compact=verified,native_bindings=bindings,
        unique_passing_methods=74,methods=[list(v)for v in sorted(methods)],new_governing_attempts=25,accepted=25,new_CP0=0,new_builds=0,new_strict=11,
        passive=dict(diagnostic_mass_factorizations=mass_factors,diagnostic_mass_solves=mass_solves,kernel_queries=queries,saved_governing_matrix_factorizations=0,saved_governing_matrix_solves=0),
        scope='Linux heavy-byte verification and independent native-header role checks; Mac transport is separately verified',source=identity(__file__)))
    print('H11 verification PASS; temporal',common['temporal'],'25 attempts/11 views/74 methods',flush=True)


def finalize(root):
    out=root/'closure';verification=load(out/'verification.json');assert verification['status']=='PASS';book=load(root/'accounting.json');assert not book.get('closed')
    assert all(v['status']!='RESERVED'and 'seconds'in v for v in book['entries'].values())
    actual=sum(v['seconds']for v in book['entries'].values());charge=60+actual+240;peak=max(v['peak_tree_rss_bytes']for v in book['entries'].values())
    copied=sum(size(root/p)for p in ('source','runtime'));science=size(root/'science');assert actual<=1200 and charge<=1500 and peak<=4*1024**3 and copied<=1024**3 and science<=3*1024**3
    gates=load(root/'publication-v1/component-gates.json');book.update(closed=True,status=gates['status'],resource_status='PASS',operation_seconds=actual,charged_seconds=charge,setup_charge_seconds=60,closure_charge_seconds=240,peak_process_tree_RSS=peak,copied_source_runtime_bytes=copied,science_bytes=science,
        closure_scope='full240s charged for parent inspections/source/figure QA/publication transport; no reserve reclaimed or historical ledger reopened')
    save(root/'accounting.json',book);save(out/'accounting-final.json',book)
    operations=[]
    for label in book['entries']:
        p=root/'operations'/label;operations.append(dict(label=label,request=identity(p/'request.json'),receipt=identity(p/'receipt.json'),log=identity(p/'command.log'),maps=identity(p/'process-maps.json')))
    save(out/'bridge.json',dict(status='PASS',qualification_status=gates['status'],verification=identity(out/'verification.json'),final_accounting=identity(out/'accounting-final.json'),publication_seal=identity(root/'publication-v1/seal.json'),
        local_tables_seal=identity(root/'local-comparisons-v1/seal.json'),original_H10_accounting=check(identity(H10/'closure/accounting-final.json')),input_intent=identity(root/'science/input-intent.json'),
        source_archive=load(root/'source-archive.json'),source_versions=[identity(p)for directory in sorted(root.glob('observer-v*'))for p in sorted(directory.glob('*.py'))],operations=operations,
        limits=dict(total_seconds=1500,operations_seconds=1200,RSS=4*1024**3,science=3*1024**3,source_runtime=1024**3,workers=72),scope='completed fixed H11 matrix; native decisions, finite coverage and UNRESOLVED temporal result remain separate',source=identity(__file__)))
    save(out/'component-gates-final.json',dict(gates,closure_transport='REMOTE_VERIFIED',resource_accounting='PASS',closure_verification=identity(out/'verification.json'),accounting=identity(out/'accounting-final.json')))
    save(out/'seal.json',dict(files=[identity(p)for p in sorted(out.iterdir())if p.is_file()]))
    print(dict(status=gates['status'],resource='PASS',charged_seconds=charge,peak_RSS=peak,science_bytes=science,copied_source_runtime=copied),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--finalize',action='store_true');a=p.parse_args();(finalize if a.finalize else verify)(a.root)
