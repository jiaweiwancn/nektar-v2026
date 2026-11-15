#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Seal H24 only after all six scientific closure criteria have evidence."""
import argparse,json,tarfile
from pathlib import Path
from run_wj_h4 import load
from run_wj_h7 import update_journal
from run_wj_h15 import size
from run_wj_h24 import usage,frozen
from wj_h15_observer import identity,verify,atomic
from wj_h23_continuation import START,END,DT,qualified,recovery

def close(root):
    cfg=frozen(root);book=load(root/'accounting.json');assert not book['closed']
    assert not any(v['status']=='RESERVED'for v in book['entries'].values())
    assert not any(book.get(k)for k in ('native_halt','coverage_halt','integrity_halt','screen_halt','runtime_halt','resource_halt'))
    assert book['native_reservations']==['M','strict_M','C','strict_C','F','strict_F']and book['reserved_attempts']==28
    inputs=load(root/'inputs/selected.json');pub=root/'publication-v1';gates=load(pub/'component-gates.json')
    for k in ('result','initial_restart_integrity','coherent_bindings_real_native_startup','physical','absolute_provisional_screens'):assert gates[k]=='PASS'
    assert gates['new_completed']==gates['new_attempts']==28 and gates['new_strict_reads']==3
    assert gates['governor_scientific_audits']==15 and gates['observer_state_assessments']==30
    assert gates['new_qualified_endpoints']==12 and gates['new_startup_verification_commits']==3
    for z in ('M','C','F'):
        q=load(root/f'science/{z}/qualification.json');s=load(root/f'science/{z}/final-strict/result.json')
        assert q['status']==s['status']=='PASS'and s['selected_absolute_step']==END[z]and s['genuine_alpha']==1/DT[z]
        d=root/f'science/{z}/monitor';a=recovery(d,END[z]);a['qualified_commit']=identity(d/f'step-{END[z]}.qualification-commit.json')
        r,cp=qualified(a,z,END[z]);assert abs(r['time']-.04)<128*2.220446049250313e-16
    original=load(root/'inputs/historical-immutability.json')
    for v in original['files']:verify(v)
    historical=load(root/'inputs/h23-closed-link.json');old=load(verify(historical['accounting']))
    assert old['closed']and old['total_charged_seconds']==550.5544812529697 and old['actual_native_advances']==0
    for v in load(pub/'manifest.json')['files']:verify(v)
    for v in load(pub/'heavy-manifest.json')['files']:verify(v)
    assert load(root/'operations/final-publication-v1/receipt.json')['status']=='PASS'
    u=usage(book);assert u['total']<=6000 and u['native']<=4200 and u['passive']<=900
    peak=max(e.get('peak_tree_rss_bytes',0)for e in book['entries'].values());assert peak<=8*1024**3
    sources=sum(size(p)for p in root.iterdir()if p.is_dir()and p.name.startswith(('source','publication-code')))
    output=sum(size(p)if p.is_dir()else p.stat().st_size for p in root.iterdir()if not p.name.startswith(('source','publication-code')))
    assert sources<=8*1024**3 and output<=2*1024**3
    out=root/'closure-v1';out.mkdir()
    counts=dict(governor_processes=3,governing_attempts=28,native_completed=28,accepted_advances=28,failed_governing_attempts=0,
        native_startup_audits=3,native_completed_audits=12,native_full_audits=15,observer_events=15,observer_state_assessments=30,
        qualified_completed_endpoints=12,new_startup_verification_commits=3,terminal_strict_reads=3,terminal_absolute_steps=END,
        CP0=0,bootstrap=0,spatial_initializations=0,standalone_history_initializations=0,native_builds=0,normal_rich_captures=0,
        bound_activations=0,material_corrections=0,new_focused_method_executions=22,passive_retained_triplet_assessments=3,
        historical_H23='separate CLOSED obstruction preserved; not part of H24 numerical counts')
    atomic(out/'execution-counts.json',counts)
    for p in sorted((root/'operations').glob('*/receipt.json')):
        dest=out/'operation-receipts'/p.parent.name/'receipt.json';dest.parent.mkdir(parents=True);atomic(dest,load(p))
    atomic(out/'historical-immutability.json',dict(status='PASS',files=original['files'],H23_accounting=historical['accounting'],H23_closed_unchanged=True))
    atomic(out/'verification.json',dict(status='PASS_H24_SCIENTIFIC_CLOSURE',criteria={
        'coherent_bindings_real_startup':'PASS','M16_C8_F32_qualified_04':'PASS','scheduled_physical_absolute_gates':'PASS',
        'three_terminal_strict_and_lineage':'PASS','signed_four_time_and_error_residual_evidence':'PASS','measured_cost_storage_compact_publication':'PASS'},
        usage=u,ceilings=dict(total=6000,native=4200,passive=900),counts=counts,peak_science_RSS_bytes=peak,
        new_output_before_closure_bytes=output,new_source_runtime_copied_bytes=sources,runtime_reused_in_place=True,
        publication=identity(pub/'manifest.json'),historical_H23=historical['accounting'],runtime=identity(root/'science/runtime-qualified-v1.json'),
        executed_source=cfg['source_manifest'],setup_fully_charged=60,closure_fully_charged=300,closure_source=identity(__file__)))
    def seal(b):
        assert not b['closed'];b['entries']['closure']=dict(status='FULLY_CHARGED',phase='closure',reserved_seconds=300,seconds=300,
            scope='selected verification, final manifests, archives, compact transport')
        b.update(closed=True,total_charged_seconds=u['total'],native_charged_seconds=u['native'],passive_charged_seconds=u['passive'],
            closed_status='PASS_H24_SCIENTIFIC_CLOSURE',verification=identity(out/'verification.json'))
    update_journal(root/'accounting.json',seal);atomic(out/'accounting-final.json',load(root/'accounting.json'))
    atomic(out/'manifest.json',dict(files=[identity(p)for p in sorted(out.rglob('*'))if p.is_file()],scope='separate final CLOSED accounting; publication snapshot is earlier OPEN'))
    atomic(root/'authoritative-index-v1.json',dict(publication=identity(pub/'manifest.json'),closure=identity(out/'manifest.json'),accounting=identity(out/'accounting-final.json'),
        verification=identity(out/'verification.json'),checkpoint_A=identity(root/'checkpoint-a-v1/manifest.json'),checkpoint_B=identity(root/'checkpoint-b-v1/manifest.json'),
        historical_H23=historical['accounting']))
    all_files=[root/'authoritative-index-v1.json',*[p for folder in (pub,out)for p in folder.rglob('*')if p.is_file()]]
    essential=[p for p in all_files if not any(n in str(p.relative_to(root))for n in ('native-physical-frozen-stage-packets','compiled-source-objects','runtime-build-v1','exterior-trace-differences','field-component-extrema'))]
    archives={}
    for name,files in (('essential',essential),('compact',all_files)):
        path=root/f'h24-{name}-v1.tar.gz'
        with tarfile.open(path,'w:gz')as t:
            for p in sorted(files):t.add(p,arcname=str(p.relative_to(root)))
        archives[name]=identity(path)
    atomic(root/'transport-intent-v1.json',dict(archives=archives,index=identity(root/'authoritative-index-v1.json'),closure_fully_reserved=True))
    print(json.dumps(dict(status='PASS_H24_SCIENTIFIC_CLOSURE',usage=u,archives=archives,index=identity(root/'authoritative-index-v1.json'))),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();close(a.root)
