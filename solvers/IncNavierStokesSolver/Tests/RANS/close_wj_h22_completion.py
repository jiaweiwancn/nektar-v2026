#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Close the new linked ledger without changing the stopped H22 allocation."""
import argparse,json,tarfile
from pathlib import Path
from run_wj_h4 import load
from run_wj_h7 import update_journal
from run_wj_h15 import size
from wj_h15_observer import identity,verify,atomic
from run_wj_h22_completion import usage,PARENT,ROOT

def close(root):
    book=load(root/'accounting.json');assert not book['closed']
    assert book['new_native_reservations']==['verify_M2','strict_M2','F','strict_F']and book['new_attempt_capacity']==4
    assert not any(book.get(k)for k in ('native_halt','coverage_halt','integrity_halt','runtime_halt'))
    assert not any(v['status']=='RESERVED'for v in book['entries'].values())
    for role,name in (('M2','verification.json'),('F','qualification.json')):assert load(root/'science'/role/name)['status']=='PASS'
    for role,step,alpha in (('M2',2,400),('F',4,800)):
        v=load(root/'science'/role/'final-strict/result.json');assert v['status']=='PASS'and v['selected_absolute_step']==step and v['genuine_alpha']==alpha
    original=load(root/'inputs/original-immutability.json')
    for p in original['files']:verify(p)
    for p in original['absent']:assert not Path(p).exists()
    parent=load(verify(book['parent_accounting']));assert parent['closed']and parent['total_charged_seconds']==book['carry']['total']
    assert parent['target_halt']and parent['integrity_halt']
    counts=load(root/'science/phase-analysis-v2/counts.json');assert counts['attempts']==7 and counts['total_roles']==70
    assert load(root/'science/whole-field-scatter-v1/result.json')['status']=='PASS'
    pub=root/'publication-v1'
    for v in load(pub/'manifest.json')['files']:verify(v)
    for v in load(pub/'heavy-manifest.json')['files']:verify(v)
    for v in load(root/'execution-config.json')['files']:verify(v)
    assert load(root/'operations/completed-publication-v1/receipt.json')['status']=='PASS'
    for a in load(pub/'online-counts.json'):
        assert all(a[k]==1 for k in ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls'))
        assert a['bound_activations']==a['additional_correction_or_restoration_solves']==0 and a['raw_k_sha256']==a['corrected_k_sha256']
    u=usage(book);assert u['total']<=4500 and u['native']<=1800 and u['analysis']<=600 and u['build']<=900
    new_output=size(root/'science')+size(pub);retained_output=202392692
    new_source=sum(size(p)for p in root.iterdir()if p.is_dir()and p.name.startswith(('source','publication-code')))
    assert new_output+retained_output<=2*1024**3 and new_source+156584056<=8*1024**3
    peak=max(v.get('peak_tree_rss_bytes',0)for v in book['entries'].values());assert peak<=8*1024**3
    out=root/'closure-v1';out.mkdir()
    execution=dict(new_governing_advances=4,aggregate_governing_advances=7,aggregate_native_completions=7,
        new_governing_blocks=1,aggregate_governing_blocks=3,new_zero_advance_verification_processes=1,
        new_strict_reads=2,aggregate_strict_reads=4,strict_roles=['old C1','old M0 (not M2 evidence)','new verified M2','new F4'],
        aggregate_native_full_audits=7,aggregate_independent_assessments=6,new_independent_assessments=3,
        new_vector_packets=4,aggregate_vector_packets=7,new_phase_owner_roles=32,new_predecessor_roles=8,
        aggregate_phase_owner_roles=56,aggregate_predecessor_roles=14,new_CP0=0,new_spatial_or_history_productions=0,
        new_native_builds=0,aggregate_native_builds=1,new_scientific_steps=0,kernel_requests=0,saved_matrix_operations=0,
        additional_Python_mass_solves=0,new_focused_Python_methods=21,new_whole_scatter_synthetic_checks=1,
        original_M2_execution_qualification='NOT_COMMITTED_UNCHANGED',subsequent_M2_startup_verification='PASS')
    atomic(out/'execution-counts.json',execution)
    for p in sorted((root/'operations').glob('*/receipt.json')):
        dest=out/'operation-receipts'/p.parent.name/'receipt.json';dest.parent.mkdir(parents=True,exist_ok=True);atomic(dest,load(p))
    atomic(out/'original-evidence-immutability.json',dict(status='PASS',all_original_pins_unchanged=True,files=original['files'],
        original_absent_paths_still_absent=original['absent'],old_failure_and_closed_ledger_unchanged=True,
        new_startup_commit=identity(root/'science/M2/monitor/startup.verification-commit.json')))
    atomic(out/'verification.json',dict(status='PASS_CLOSED_FINITE_H22_DIAGNOSIS',aggregate_usage=u,aggregate_ceiling_seconds=4500,
        child_charge_seconds=u['total']-book['carry']['total'],carried_charges=book['carry'],new_setup_fully_charged=60,new_closure_fully_charged=300,
        new_peak_scientific_RSS_bytes=peak,retained_plus_new_science_publication_bytes=retained_output+new_output,
        retained_plus_new_source_runtime_bytes=156584056+new_source,counts=execution,parent_accounting=book['parent_accounting'],
        publication=identity(pub/'manifest.json'),runtime=identity(root/'science/runtime-qualified-v1.json'),closure_source=identity(__file__)))
    def seal(b):
        assert not b['closed'];b['entries']['closure']=dict(status='PASS',phase='closure',reserved_seconds=300.,seconds=300.,fully_charged=True,
            scope='selective verification, manifests/archive and bounded compact transport; active logs excluded')
        b.update(closed=True,aggregate_total_charged_seconds=u['total'],child_total_charged_seconds=u['total']-b['carry']['total'],
            aggregate_native_charged_seconds=u['native'],aggregate_analysis_charged_seconds=u['analysis'],closed_status='PASS_FINITE_DIAGNOSIS',final_verification=identity(out/'verification.json'))
    update_journal(root/'accounting.json',seal)
    atomic(out/'accounting-final.json',load(root/'accounting.json'))
    atomic(out/'manifest.json',dict(files=[identity(p)for p in sorted(out.rglob('*'))if p.is_file()],scope='separate CLOSED linked accounting; publication snapshot remains OPEN'))
    atomic(root/'authoritative-index-v1.json',dict(publication=identity(pub/'manifest.json'),closure=identity(out/'manifest.json'),
        accounting=identity(out/'accounting-final.json'),verification=identity(out/'verification.json'),parent_index=book['parent_index']))
    # Essentials first; large exact copied payloads remain indexed for optional transport.
    essential=[root/'authoritative-index-v1.json',*[p for p in out.rglob('*')if p.is_file()]]
    excluded=('attempt-bindings','retained-preflight','strict/fresh_M2.json','strict/fresh_F4.json','compiled-source-objects','runtime-build-v1')
    essential.extend(p for p in pub.rglob('*')if p.is_file()and not any(x in str(p.relative_to(pub))for x in excluded))
    archive=root/'h22-completion-essential-v1.tar.gz'
    with tarfile.open(archive,'w:gz')as t:
        for p in sorted(set(essential)):t.add(p,arcname=str(p.relative_to(root)))
    all_files=[root/'authoritative-index-v1.json',*[p for folder in (pub,out)for p in folder.rglob('*')if p.is_file()]]
    full=root/'h22-completion-compact-v1.tar.gz'
    with tarfile.open(full,'w:gz')as t:
        for p in sorted(all_files):t.add(p,arcname=str(p.relative_to(root)))
    atomic(root/'transport-intent-v1.json',dict(essential=identity(archive),compact=identity(full),essential_files=[identity(p)for p in sorted(set(essential))],
        remote_optional=[identity(p)for p in sorted(set(all_files)-set(essential))],closure_transport_fully_reserved=True))
    print(json.dumps(dict(index=identity(root/'authoritative-index-v1.json'),essential=identity(archive),compact=identity(full),aggregate_usage=u,counts=execution)),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();assert a.root.name==ROOT;close(a.root)
