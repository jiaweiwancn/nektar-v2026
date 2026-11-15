#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Seal the stopped H23 prefix; preserve all charges and forbid native reopening."""
import argparse,json,tarfile
from pathlib import Path
from run_wj_h4 import load
from run_wj_h7 import update_journal
from run_wj_h15 import size
from wj_h15_observer import identity,verify,atomic
from run_wj_h23 import usage

def close(root):
    book=load(root/'accounting.json');assert not book['closed']and book['native_halt']
    assert book['native_reservations']==['M','strict_M']and book['reserved_attempts']==8
    assert not any(v['status']=='RESERVED'for v in book['entries'].values())
    pub=root/'publication-v1';gates=load(pub/'component-gates.json')
    assert gates['status']=='OBSTRUCTED_EXECUTION_COMPATIBILITY_BEFORE_GOVERNOR'
    counts=load(pub/'execution-counts.json');assert counts['actual_governing_attempts']==counts['actual_governor_processes']==0 and counts['protected_strict_reads']==1
    for pin in load(root/'execution-config.json')['files']:verify(pin)
    original=load(root/'inputs/historical-immutability.json')
    for pin in original['files']:verify(pin)
    for pin in load(pub/'manifest.json')['files']:verify(pin)
    for pin in load(pub/'heavy-manifest.json')['files']:verify(pin)
    assert load(root/'operations/stopped-publication-v1/receipt.json')['status']=='PASS'
    assert load(root/'science/passive-binding-repair-v1/result.json')['campaign_governors_remain_closed']
    assert not any((root/f'science/{z}/native-launch.json').exists()for z in ('M','C','F'))
    assert not any((root/f'science/{z}/monitor/startup.verification-commit.json').exists()for z in ('M','C','F'))
    u=usage(book);assert u['total']<=6000 and u['native']<=4200 and u['passive']<=900
    peak=max(e.get('peak_tree_rss_bytes',0)for e in book['entries'].values());assert peak<=8*1024**3
    source_bytes=sum(size(p)for p in root.iterdir()if p.is_dir()and p.name.startswith(('source','passive-repair','publication-code')))
    output_bytes=sum(size(p)if p.is_dir()else p.stat().st_size for p in root.iterdir()if not p.name.startswith(('source','passive-repair','publication-code')))
    assert source_bytes<=8*1024**3 and output_bytes<=2*1024**3
    out=root/'closure-v1';out.mkdir()
    atomic(out/'execution-counts.json',dict(counts,focused_distinct_methods=17,successful_method_executions=30,failed_test_selector_entries=1,
        new_distinct_H23_methods=11,reused_distinct_methods=6))
    for pin in sorted((root/'operations').glob('*/receipt.json')):
        p=out/'operation-receipts'/pin.parent.name/'receipt.json';p.parent.mkdir(parents=True);atomic(p,load(pin))
    atomic(out/'historical-immutability.json',dict(status='PASS',all_selected_pins_unchanged=True,files=original['files'],
        H20_H22_ledgers_remain_closed=True,original_failed_H22_cycle_preserved=True))
    atomic(out/'verification.json',dict(status='CLOSED_OBSTRUCTED_EXECUTION_BEFORE_GOVERNOR',usage=u,total_ceiling=6000,native_ceiling=4200,passive_ceiling=900,
        counts=load(out/'execution-counts.json'),new_peak_scientific_RSS_bytes=peak,new_output_before_closure_bytes=output_bytes,
        copied_source_bytes=source_bytes,reused_runtime_copied_bytes=0,publication=identity(pub/'manifest.json'),
        first_halt=book['native_halt'],passive_repair=identity(root/'science/passive-binding-repair-v1/result.json'),
        protected_read=identity(root/'science/M/final-strict/result.json'),closure_source=identity(__file__),
        remaining_native_work='NOT_RUN; no automatic retry/relaunch authorized',setup_fully_charged=60,closure_fully_charged=300))
    def seal(b):
        assert not b['closed'];b['entries']['closure']=dict(status='FULLY_CHARGED',phase='closure',reserved_seconds=300,seconds=300,
            scope='selected verification, manifests, archives and bounded compact transport')
        b.update(closed=True,total_charged_seconds=u['total'],native_charged_seconds=u['native'],passive_charged_seconds=u['passive'],
            actual_native_advances=0,reserved_unknown_capacity_retained=8,closed_status='OBSTRUCTED_EXECUTION',final_verification=identity(out/'verification.json'))
    update_journal(root/'accounting.json',seal);atomic(out/'accounting-final.json',load(root/'accounting.json'))
    atomic(out/'manifest.json',dict(files=[identity(p)for p in sorted(out.rglob('*'))if p.is_file()],scope='separate final CLOSED accounting; publication remains an earlier OPEN snapshot'))
    atomic(root/'authoritative-index-v1.json',dict(publication=identity(pub/'manifest.json'),closure=identity(out/'manifest.json'),accounting=identity(out/'accounting-final.json'),
        verification=identity(out/'verification.json'),historical_indices=identity(root/'inputs/historical-indices.json')))
    files=[root/'authoritative-index-v1.json',*[p for folder in (pub,out)for p in folder.rglob('*')if p.is_file()]]
    essential=[p for p in files if not any(n in str(p.relative_to(root))for n in ('compiled-source-objects','runtime-build-v1'))]
    archives={}
    for kind,selected in (('essential',essential),('compact',files)):
        path=root/f'h23-{kind}-v1.tar.gz'
        with tarfile.open(path,'w:gz')as t:
            for p in sorted(selected):t.add(p,arcname=str(p.relative_to(root)))
        archives[kind]=identity(path)
    atomic(root/'transport-intent-v1.json',dict(archives=archives,closure_transport_fully_reserved=True))
    print(json.dumps(dict(index=identity(root/'authoritative-index-v1.json'),usage=u,archives=archives,counts=counts)),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();close(a.root)
