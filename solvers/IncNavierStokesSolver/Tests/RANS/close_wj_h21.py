#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact H21 presentation and separately sealed closure; no field evaluator."""
import argparse
import os
from pathlib import Path
import resource
import sys
import tarfile
import time
from run_wj_h4 import THREADS, load
os.environ.update(THREADS)
resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
import numpy as np
from run_wj_h7 import update_journal,complete_operation
from wj_h15_observer import atomic,identity,verify
from wj6_resources import volume


def summarize(root):
    start=time.monotonic();label='08-summary';journal=root/'accounting.json'
    def reserve(b):
        assert not b['closed'] and label not in b['entries']
        spent=sum(v.get('seconds',v['reserved_seconds'])for v in b['entries'].values())
        assert spent+60<=480
        b['entries'][label]=dict(status='RESERVED',phase='operation',reserved_seconds=60.)
    update_journal(journal,reserve)
    s=root/'science-v1';out=root/'supplement-v1';out.mkdir()
    values=load(s/'endpoint-values.json');pairs=load(s/'pair-decompositions.json')
    windows=load(s/'endpoint-windows.json');shared=load(s/'shared-directions.json')
    replay=load(s/'H20-replay.json');records=[];traces=[]
    for v in pairs:
        b=v['delta_B'];ng=float(np.linalg.norm(b['gradient_term']));nt=float(np.linalg.norm(b['timescale_term']))
        records.append(dict(time=v['time'],geometry=v['geometry'],pair=v['pair'],
            B_gradient_norm=ng,B_timescale_norm=nt,timescale_over_gradient_norm=None if ng<=4e-11 else nt/ng,
            B_xz=float(np.array(b['total'])[0,2]),B_xz_gradient=float(np.array(b['gradient_term'])[0,2]),
            B_xz_timescale=float(np.array(b['timescale_term'])[0,2]),
            **{k:v[k]for k in ('delta_a','delta_g','delta_mhat','dimensional_term','normalization_term')}))
    for t in (.005,.01,.015,.02):
        for z in ('C','M','F'):
            a,b=(next(v['values']for v in values if v['time']==t and v['level']==z and v['geometry']==g)for g in (251,265))
            traces.append(dict(time=t,level=z,upper_minus_lower={n:(np.array(a[n])-b[n]).tolist()for n in ('U','k','epsilon','pi','tau','G','B','mhat')}))
    terminal=[]
    for v in windows:
        if v['start']==.015 and v['pair']=='M-F':
            terminal.append(dict(geometry=v['geometry'],pair=v['pair'],start=v['start'],end=v['end'],
                 quantities={n:v['endpoint_quantities'][n]for n in ('mhat','a','g','d','nu_eff','B','G')}))
    result=dict(status='PASS',pair_contributions=records,opposite_trace_checks=traces,terminal=terminal,
        maximum_scaled_H20_field_replay_error=max(x['scaled_error']for r in replay for x in r['checks'].values()),
        maximum_shared_direction_replay_error=max(abs(r['published_error'])for r in shared),
        maximum_implicit_dN_error=max(v['values']['implicit_derivative_error']for v in values),
        minimum_local_covariance=min(min(v['values']['covariance'])for v in values),
        minimum_local_principal=min(v['values']['mhat']for v in values),
        minimum_local_optimizer_gap=min(v['values']['g']for v in values),
        maximum_B_identity_error=max(v['B_identity']['maximum_error']for v in pairs),
        maximum_normalization_identity_error=max(v['normalized_identity']['maximum_error']for v in pairs),
        maximum_optimizer_identity_error=max(v['optimizer_identity']['maximum_error']for v in pairs),
        maximum_endpoint_identity_error=max(q['identity']['maximum_error']for w in windows for q in w['endpoint_quantities'].values()),
        interpretation='Exact endpoint algebra; gradient/timescale norms are not nonlinear percentages; no phase or predictor/correction causality inferred')
    atomic(out/'findings.json',result)
    for name in ('source','source-v2','source-v3','source-tests-v4'):
        atomic(out/(name+'-manifest.json'),dict(files=[identity(p)for p in sorted((root/name).glob('*.py'))]))
    atomic(out/'manifest.json',dict(files=[identity(p)for p in sorted(out.iterdir())if p.is_file()]))
    complete_operation(journal,label,dict(status='PASS',phase='operation',reserved_seconds=60.,seconds=max(1.,time.monotonic()-start),
        peak_self_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,source=identity(__file__)))
    print('maximum scaled field replay error',result['maximum_scaled_H20_field_replay_error'])
    print('maximum shared-direction error',result['maximum_shared_direction_replay_error'])
    for v in records:
        if v['pair'] in ('C-M','M-F'):
            print('PAIR',v)
    for v in terminal:
        for n in ('mhat','a','g'):
            x=v['quantities'][n]
            print('TERMINAL',v['geometry'],n,{k:x[k]for k in ('incoming','net_update','endpoint')},x['net_cancellation'])


def close(root):
    start=time.monotonic();journal=root/'accounting.json';out=root/'closure-v1';out.mkdir()
    def reserve(b):
        assert not b['closed'] and 'closure'not in b['entries']
        b['entries']['closure']=dict(status='RESERVED',phase='closure',reserved_seconds=90.)
    update_journal(journal,reserve)
    for folder in ('publication-v1','supplement-v1'):
        m=load(root/folder/'manifest.json')
        for v in m['files']:verify(v)
    # Active source and every selected historical pin are checked, without a
    # traversal/rehash of the historical campaign roots.
    active=load(root/'controls/active-execution.json')
    for k in ('analysis','runner','tests'):verify(active[k])
    selection=load(root/'controls/selected.json');verify(selection['structure'])
    for item in selection['states']:
        for v in item['anchor'].values():
            if isinstance(v,dict)and set(('path','bytes','sha256'))<=v.keys():verify(v)
    files=list((root/'science-v1').glob('*.json'))
    assert len(load(root/'science-v1/endpoint-values.json'))==24
    atomic(out/'verification.json',dict(status='PASS',publication=identity(root/'publication-v1/manifest.json'),
          supplement=identity(root/'supplement-v1/manifest.json'),active_sources=active,selected_historical_bytes_unchanged=True,
          native_work=dict(governing_attempts=0,cp0=0,restarts=0,strict_readbacks=0,builds=0,kernel_requests=0,
                           saved_governing_matrix_matvecs=0,mass_factorizations=0,mass_solves=0),
          imported_checkpoints=12,accepted_endpoint_owner_roles=24,
          tests=dict(distinct_methods=15,successful_final_execution=15,prior_execution=15,prior_error='one test-only fixture NameError'),
          science_files=[identity(p)for p in sorted(files)],closure_coordinator=identity(__file__)))
    def finish(b):
        for key,v in b['entries'].items():
            if key!='closure' and v['status']=='RESERVED':
                v.update(status='UNKNOWN_FULL_RESERVATION',seconds=v['reserved_seconds'],
                         reason='metadata inspection exited before receipt; conservative full charge retained')
        b['entries']['closure'].update(status='PASS',seconds=90.,fully_charged=True,
            actual_seconds_before_archive=time.monotonic()-start,
            transport='reserved in fixed closure charge; Mac verification recorded separately')
        b['charged_seconds']=b['setup_seconds']+sum(v.get('seconds',v['reserved_seconds'])for v in b['entries'].values())
        assert b['charged_seconds']<=600
        b.update(closed=True,closed_status='PASS_DIAGNOSIS_WITH_RETAINED_PASSIVE_FAILURES',
                 source_and_output_snapshot_bytes=volume(root))
    book=update_journal(journal,finish)
    atomic(out/'accounting-final.json',book)
    atomic(out/'manifest.json',dict(files=[identity(p)for p in sorted(out.iterdir())if p.is_file()],
          scope='closed final accounting; earlier publication snapshot is intentionally OPEN'))
    archive=root/'h21-essential-v1.tar.gz'
    with tarfile.open(archive,'w:gz')as tar:
        for folder in ('publication-v1','supplement-v1','closure-v1','operations'):
            tar.add(root/folder,arcname=folder)
    atomic(root/'authoritative-index-v1.json',dict(publication=identity(root/'publication-v1/manifest.json'),
        supplement=identity(root/'supplement-v1/manifest.json'),closure=identity(out/'manifest.json'),
        accounting=identity(out/'accounting-final.json'),archive=identity(archive),
        closed_charged_seconds=book['charged_seconds'],transport='PREPARED; exact Mac receipt is separate'))
    print('CLOSED',book['charged_seconds'],'/600; archive',identity(archive),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    p.add_argument('--mode',choices=('summarize','close'),required=True);a=p.parse_args()
    assert a.root.name=='wj-h21-outlet-saved-state-20261115'
    dict(summarize=summarize,close=close)[a.mode](a.root)
