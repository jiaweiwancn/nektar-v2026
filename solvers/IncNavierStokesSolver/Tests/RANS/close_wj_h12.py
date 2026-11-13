#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Finalize H12's separately charged closure; never hash a mutable ledger as final."""
import argparse
import json
from pathlib import Path
import shutil
import tarfile
import numpy as np
from run_wj_h4 import load,save
from run_wj_h7 import update_journal
from run_wj_h10 import size
from study_wj_h12 import file_id,serial

def close(root):
    dest=root/'closure-v1';assert load(dest/'verification.json')['status']=='PASS'
    book=load(root/'accounting.json');assert not book['closed']
    charged=lambda phases:sum(v.get('seconds',v['reserved_seconds'])for v in book['entries'].values()if v['phase']in phases)
    operations=charged(('operation','science'));science=charged(('science',));closure=charged(('closure',))
    assert operations<=480 and science<=360 and closure<=90
    total=book['setup_seconds']+operations+book['closure_reserve'];assert total<=600
    assert all(book[k]==0 for k in ('governing_attempts','CP0s','strict_readbacks','builds','initializer_applications','diagnostic_mass_factorizations','diagnostic_mass_solves'))
    result=load(root/'science/results-v2/result.json');roles=load(root/'science/results-v2/local-phase-roles.json')
    witnesses=[]
    for v in load(root/'science/results-v2/cross-window-accounting.json'):
        if v['owner_index']==0 and v['observable']in ('mhat','a','g','G','B','tau','d','nu_eff'):
            if v['observable']in ('G','B'):
                index=lambda a:float(np.asarray(a)[0,2])
            else:index=float
            witnesses.append(dict(window=v['window'],pair=v['pair'],observable=v['observable'],
                incoming=index(v['incoming']),contributions={k:index(x)for k,x in v['contributions'].items()},endpoint=index(v['endpoint'])))
    summary=dict(status='PASS',scientific_temporal='H11 UNRESOLVED, unchanged',result=result,
        same_owner_pairs=[v for v in load(root/'publication-v1/endpoint-pairs.json')if v['owner_index']==0],
        geometry17_upper_window_terms=witnesses,normalization_gap=load(root/'publication-v1/window-normalization-gap.json'),
        maximum_endpoint_replay_error=max(x for v in load(root/'science/results-v2/endpoint-replay.json')['records']for x in v.get('errors',{}).values()),
        compiled_maximum_errors={k:max(v['errors'][k]for v in load(root/'science/results-v2/compiled-checks.json')['records'])for k in ('dimensional','implicit_dN','normalized','R','N')},
        phase_margin_range=[min(v['value']['mhat']for v in roles),max(v['value']['mhat']for v in roles)],
        maximum_implicit_derivative_error=max(abs(v['value']['implicit_derivative_error'])for v in roles),
        matrix_operations=dict(saved_governing_matrix_reads=0,saved_governing_matvecs=0,factorizations=0,solves=0,
            sparse_basis_block_products_initial=20*result['unique_physical_state_evaluations'],sparse_reference_check_block_products=20,
            independent_reference_scalar_dots=60,saved_increment_sparse_block_products=240,
            saved_velocity_vector_scatters=9,scope='3-row native basis dot products and small constitutive tensors, not governing/mass systems'))
    save(dest/'scientific-summary.json',serial(summary))
    versions=[]
    for pattern in ('wj-h12-source-*.json','observer-v*/source.json'):
        for p in sorted(root.glob(pattern)):versions.append(dict(metadata=file_id(p),record=load(p)))
    save(dest/'source-identities.json',versions)
    def finish(b):
        b.update(closed=True,operations_seconds=operations,scientific_seconds=science,closure_operation_seconds=closure,
            charged_seconds=total,closure_fully_charged_seconds=90.,remaining_seconds=600-total,
            peak_process_tree_rss_bytes=max(v['peak_tree_rss_bytes']for v in b['entries'].values()),
            output_bytes_before_transport=size(root)-size(root/'source'),copied_source_bytes=size(root/'source'),
            kernel_requests=result['kernel_requests'],kernel_sessions=1,diagnostic_completion='PASS',H11_temporal='UNRESOLVED',
            passive_failure_resolution='original KeyError and charge retained; aggregation resumed once from unchanged phase records under verified passive-repair receipt')
    update_journal(root/'accounting.json',finish)
    save(dest/'accounting-final.json',load(root/'accounting.json'))
    files=[file_id(p)for p in sorted(dest.iterdir())if p.is_file()]
    receipts=[file_id(p)for p in sorted((root/'operations').glob('*/*'))if p.is_file()]
    save(dest/'manifest.json',dict(files=files,operation_files=receipts,publication_manifest=file_id(root/'publication-v1/manifest.json'),
        publication_seal=file_id(root/'publication-v1/seal.json'),setup=file_id(root/'health-setup.json'),source=file_id(__file__),
        scope='final accounting saved after all operation receipts completed; full 90 seconds closure charge includes this finalizer and compact transport'))
    save(dest/'seal.json',dict(manifest=file_id(dest/'manifest.json'),accounting=file_id(dest/'accounting-final.json')))
    delivery=root/'delivery-v1';delivery.mkdir()
    for name in ('publication-v1','closure-v1','operations'):shutil.copytree(root/name,delivery/name)
    shutil.copyfile(root/'health-setup.json',delivery/'health-setup.json')
    archive=root/'wj-h12-compact-v1.tar.gz'
    with tarfile.open(archive,'w:gz')as f:
        for p in sorted(delivery.iterdir()):f.add(p,arcname=p.name)
    save(root/'compact-transport-v1.json',dict(archive=file_id(archive),original_final_accounting=file_id(dest/'accounting-final.json'),
        closure_seal=file_id(dest/'seal.json'),scope='selected compact/derived outputs only; no heavy coefficients, checkpoints, matrices or native executables'))
    print('CLOSED',json.dumps({k:load(dest/'accounting-final.json')[k]for k in ('charged_seconds','scientific_seconds','peak_process_tree_rss_bytes','output_bytes_before_transport','kernel_requests')}),flush=True)
    print('PRIMARY',json.dumps(serial(witnesses)),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);close(p.parse_args().root)
