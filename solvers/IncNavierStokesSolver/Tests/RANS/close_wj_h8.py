#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded H8 closure verification; no native initialization, solve or replay."""
import argparse
import json
from pathlib import Path
from run_wj_h8 import load,save,sha,size


def check(v):
    p=Path(v['path']);assert p.stat().st_size==v['bytes']and sha(p)==v['sha256'],str(p)


def verify(root,out):
    out.mkdir();pub=root/'publication-v2';seal=load(pub/'seal.json')
    for v in seal['files']:check(v)
    gates=load(pub/'component-gates.json')
    assert gates['verification_status']=='PASS'and gates['actual_runtime']=='PASS'and gates['complete_tuple_coverage']and gates['strict_weak']=='PASS'
    assert gates['status']=='UNRESOLVED'and gates['local_timestep_contraction']['status']=='UNRESOLVED','retain measured early noncontraction'
    for v in load(pub/'manifest.json')['heavy']:check(v)
    imports=load(root/'imports-v2.json');ncp=0
    for item in imports['histories']:
        check(item['report']);check(item['xml'])
        for cp in item['checkpoints']:check(cp['copy']);ncp+=1
    for v in load(root/'observer-lineage.json')['files']:check(v)
    rt=load(root/'runtime-equivalence-v2.json')
    for v in rt['compiled_inputs']:assert sha(root/'source'/v['path'])==v['sha256']
    for v in rt['private_libraries']:check(v)
    for v in rt['executables'].values():check(v['ELF']);check(v['wrapper'])
    maps=load(pub/'actual-native-runtime.json');assert len(maps)==9 and all(v['status']=='PASS'for v in maps.values())
    expected={'coarse-1':100,'coarse-3':100,'fine-2':200,'fine-6':200,'quarter-4':400,'quarter-12':400,'quarter-16':400}
    views=load(pub/'strict-weak-stage-views.json');assert set(views)==set(expected)
    for name,v in views.items():assert v['status']=='PASS'and v['context']['alpha']==expected[name];check(v['native']);check(v['context']['producer']);check(v['context']['checkpoint'])
    counts=gates['new_native'];assert counts['attempted']==counts['accepted']==16 and not counts['failure']and counts['initialization_work']==[0,0]
    for key in ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls'):assert counts[key]==16
    assert counts['bound_activations']==0 and load(root/'case/initial-pair.json')['status']=='PASS'
    book=load(root/'accounting.json');assert book['cp0_reservations']==1 and book['strict_view_reservations']==7 and not book.get('coverage_halt')
    assert sum(e.get('units',0)for e in book['entries'].values())==16
    assert size(root)<8*1024**3
    for v in book['entries'].values():assert v.get('peak_tree_rss_bytes',0)<8*1024**3
    for label in ('focused-tests','ctest-focused','focused-h8-final-v2'):assert book['entries'][label]['status']=='PASS'
    text=(root/'operations/ctest-inventory/command.log').read_text();assert 'Total Tests: 4'in text
    logs=out/'logs';logs.mkdir()
    for name in ('focused-tests','header-regression','ctest-inventory','ctest-focused','focused-h8-final','focused-h8-final-v2'):
        p=root/'operations'/name/'command.log';(logs/(name+'.log')).write_bytes(p.read_bytes())
    save(out/'verification.json',dict(status='PASS',science_status=gates['status'],native_attempts=16,accepted=16,CP0=1,strict_views=7,imported_native_checkpoints=ncp,
        tuple_coverage=gates['actual_accepted_records'],loaded_native_process_receipts=len(maps),compiled_inputs=len(rt['compiled_inputs']),qualified_ELFs=len(rt['qualified_ELFs']),
        root_bytes=size(root),selected_CTests=4,unique_focused_methods=36,scope='fresh Linux hash verification of selected heavy bytes; Mac independently verifies compact transport only',
        original_H7_ledger='closed1800 and107.85679090699705 unknown charge unchanged',source=dict(path=str(Path(__file__)),sha256=sha(__file__))))
    print('PASS H8 closure: native/strict histories, original private runtime, hashes and truthful noncontraction',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True,type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args();verify(a.root,a.output)
