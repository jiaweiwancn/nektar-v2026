#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Seal the stopped H22 prefix; never repair an ACK, commit or native cursor."""
import argparse,json,shutil,tarfile
from pathlib import Path
from run_wj_h4 import load
from run_wj_h7 import update_journal
from run_wj_h15 import size
from run_wj_h22 import usage
from wj_h15_observer import identity,verify,atomic

def close(root):
    journal=root/'accounting.json';book=load(journal)
    assert not book['closed']and book['target_halt']and book['integrity_halt']
    assert book['target_levels']==['C','M']and book['target_reserved_attempts']==3
    assert book['strict_reservations']==['C','M']
    assert not any(v['status']=='RESERVED'for v in book['entries'].values())
    assert not(root/'science/F/native-launch.json').exists()
    assert not(root/'science/M/monitor/step-2.qualification-commit.json').exists()
    assert(root/'science/C/monitor/step-1.qualification-commit.json').exists()
    for label in ('focused-tests-v1','prepared-session-test-v1','cursor-schema-test-v1','audit-role-regressions-v1',
                  'sealed-audit-selector-test-v1','retained-role-and-phase-proof-v1','runtime-executed-closure-v1','publication-v1'):
        assert load(root/'operations'/label/'receipt.json')['status']=='PASS'
    assert load(root/'science/C/final-strict/result.json')['selected_absolute_step']==1
    assert load(root/'science/M/final-strict/result.json')['selected_absolute_step']==0
    proof=load(root/'science/passive-role-fix-v1/result.json')
    verify(proof['qualified_cursor_unchanged'])
    pub=root/'publication-v1'
    for v in load(pub/'manifest.json')['files']:verify(v)
    for v in load(pub/'heavy-manifest.json')['files']:verify(v)
    for v in load(root/'execution-config.json')['files']:verify(v)
    for v in load(root/'source-final-v5/manifest.json')['files']:verify(v)
    for a in load(pub/'online-counts.json'):
        assert all(a[k]==1 for k in ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls'))
        assert a['bound_activations']==a['additional_correction_or_restoration_solves']==0
        assert a['raw_k_sha256']==a['corrected_k_sha256']
    used=usage(book)
    assert used['total']<=4500 and used['native']<=1800 and used['build']<=900 and used['analysis']<=600
    source=sum(size(p)for p in root.iterdir()if p.is_dir()and p.name.startswith(('source','build','runtime','passive-code')))
    output=size(root/'science')+size(pub)
    science_peak=max(v.get('peak_tree_rss_bytes',0)for v in book['entries'].values()if v['phase']!='build')
    build_peak=max(v.get('peak_tree_rss_bytes',0)for v in book['entries'].values()if v['phase']=='build')
    assert source<=8*1024**3 and output<=2*1024**3 and science_peak<=8*1024**3 and build_peak<=16*1024**3
    out=root/'closure-v1';out.mkdir()
    final_source=out/'final-source';final_source.mkdir()
    for p in sorted((root/'source-final-v5').glob('*')):
        if p.is_file():shutil.copyfile(p,final_source/p.name)
    shutil.copyfile(__file__,final_source/'executed-close_wj_h22.py')
    counts=dict(governing_process_blocks=2,reserved_target_capacity=3,actual_attempts=3,native_completed=3,
        new_committed_scientific_endpoints=1,completed_unqualified_M_terminal=1,off_cadence_native_completed_M1=1,
        final_strict_reads=2,strict_roles=dict(C='C1 genuine completed callback',M='original M CP0, genuine initial alpha=0'),
        F_attempts_not_run=4,F_final_strict_read='NOT_RUN',native_builds=1,synthetic_native_components=1,
        new_CP0=0,native_history_productions=0,spatial_initializations=0,new_scientific_steps=0,
        vector_packets=3,phase_owner_roles=24,predecessor_owner_roles=6,
        focused_distinct_Python_methods=24,compiled_kernel_requests=0,saved_matrix_operations=0,
        additional_diagnostic_mass_solves=0,corrected_observer_governing_qualification='NOT_RUN')
    atomic(out/'execution-counts.json',counts)
    for p in sorted((root/'operations').glob('*/receipt.json')):
        # These operations have ended; no active log or mutable journal is sealed.
        dest=out/'operation-receipts'/p.parent.name/'receipt.json';dest.parent.mkdir(parents=True,exist_ok=True)
        atomic(dest,load(p))
    atomic(out/'verification.json',dict(status='CLOSED_OBSTRUCTED_OBSERVER_ROLE',counts=counts,
        retained_failure=book['integrity_halt'],qualified_M_cursor_unchanged=proof['qualified_cursor_unchanged'],
        passive_correction=identity(root/'science/passive-role-fix-v1/result.json'),
        source_build_runtime_bytes=source,new_science_and_publication_bytes=output,
        aggregate_science_peak_RSS_bytes=science_peak,build_peak_RSS_bytes=build_peak,
        usage=used,setup_fully_charged=60,closure_fully_charged=300,
        publication=identity(pub/'manifest.json'),runtime=identity(root/'science/runtime-qualified-v1.json'),
        closure_source=identity(__file__),allocation='unused allowance closed; no F reservation or automatic retry'))
    def seal(b):
        assert not b['closed']and b['integrity_halt']==book['integrity_halt']and b['target_halt']==book['target_halt']
        b['entries']['closure']=dict(status='PASS',phase='closure',reserved_seconds=300.,seconds=300.,fully_charged=True,
            scope='final selective verification, archive/hash and compact transport; active closure log excluded')
        b.update(closed=True,total_charged_seconds=used['total'],closure_fully_charged_seconds=300,
                 closed_status='OBSTRUCTED_OBSERVER_ROLE_WITH_PASSIVE_CORRECTION',final_verification=identity(out/'verification.json'))
    update_journal(journal,seal)
    atomic(out/'accounting-final.json',load(journal))
    atomic(out/'manifest.json',dict(files=[identity(p)for p in sorted(out.rglob('*'))if p.is_file()],
        scope='separate CLOSED accounting; earlier publication snapshot stays OPEN'))
    atomic(root/'authoritative-index-v1.json',dict(publication=identity(pub/'manifest.json'),closure=identity(out/'manifest.json'),
        accounting=identity(out/'accounting-final.json'),verification=identity(out/'verification.json')))
    selected=[root/'authoritative-index-v1.json',*[p for folder in (pub,out)for p in folder.rglob('*')if p.is_file()]]
    archive=root/'h22-essential-v1.tar.gz'
    with tarfile.open(archive,'w:gz')as tar:
        for p in sorted(selected):tar.add(p,arcname=str(p.relative_to(root)))
    atomic(root/'essential-transport-intent-v1.json',dict(archive=identity(archive),files=[identity(p)for p in sorted(selected)],
        heavy_science='Linux original files only, selected hashes in publication heavy manifest'))
    print(json.dumps(dict(index=identity(root/'authoritative-index-v1.json'),archive=identity(archive),
        charged_seconds=used['total'],counts=counts)),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args()
    assert a.root.name=='wj-h22-early-phase-replay-20261115';close(a.root)
