#!/usr/bin/env python3
"""Reproduce the final read-only H3 byte/binding checks and compact metrics.

Run with PYTHONPATH pointing at the frozen H3 publication source Tests/RANS.
It neither launches a solver nor edits the original publication or H2 data.
"""
import csv
import json
from pathlib import Path
import sys
import hashlib
import subprocess
import numpy as np
from run_wj_h3 import load,save,sha

root=Path(sys.argv[1]);out=root/'closure-v1';out.mkdir(exist_ok=False)
pub=root/'publication-v1';manifest=load(pub/'manifest.json');heavy=load(pub/'heavy-index.json')
for key,value in manifest['compact'].items():
    path=pub/key
    assert path.stat().st_size==value['bytes'] and sha(path)==value['sha256'],key
for key,value in heavy.items():
    path=Path(key)
    assert path.stat().st_size==value['bytes'] and sha(path)==value['sha256'],key
runs=load(pub/'runs.json');packets=0
for run in runs:
    r=load(run['producer']['path'])
    assert len(r['accepted_states'])==len(r['wj_audit_packets'])
    for state,p in zip(r['accepted_states'],r['wj_audit_packets']):
        digest=lambda v:hashlib.sha256(np.asarray(v,dtype='<f8').tobytes()).hexdigest()
        assert digest(state)==p['state_sha256']
        if p['callback_history']:assert digest(p['callback_history'])==p['stage_history_sha256']
        packets+=1
    for event in r['physical_range_events']:
        p=event['wj_packet'];assert digest(event['state'])==p['state_sha256']
        assert digest(event['history'])==p['stage_history_sha256']
    assert run['solve_counts']==[run['attempts']]*5 and run['raw_corrected_k_exact']
summary=load(root/'analysis-v2/summary.json');diagnosis={}
for key,record in summary['cases'].items():
    d=load(record['path']);s=d['states'][-1];w=s['packet']['principal_witnesses']['in_plane']
    diagnosis[key]=dict(time=s['time'],role=s['role'],minimum=w['minimum_over_nu_eff'],position=w['position'],
        G=w['raw_G_row_major'],k=w['k'],epsilon=w['epsilon'],nu_eff=w['nu_eff'],nu_d=w['nu_d'],N=w['N'],
        covariance=w['covariance_eigenvalues'],Praw=w['Praw'],Pdev=w['Pdev'],F_N=w['F_N'],
        pressure_test_projection=s['pressure_test_projection'],
        largest_sampled_divergence=s['maximum_sampled_divergence'],
        largest_sampled_normal_H_defect=max(abs(x['normal_traction_defect']) for x in s['interface']))
rows=list(csv.DictReader((pub/'paired-local.csv').open()));contrasts={}
for pair in sorted({r['pair'] for r in rows}):
    group=[r for r in rows if r['pair']==pair and abs(float(r['time'])-.08)<1e-12]
    contrasts[pair]={key:max(abs(float(r[key])) for r in group) for key in
        ('delta_Gxz','delta_principal','delta_divergence','delta_pi','delta_Ux','delta_Uz','delta_k','delta_epsilon')}
obs=load(pub/'native-observations.json');boundary={}
for key,record in obs.items():
    p=load(record['native']['path']);b=p['open_boundary']
    boundary[key]=dict(physical_mass_dual=record['stage_norms'],mass_flux=b['mass_flux'],
        maximum_outlet_traction_defect=max(float(np.linalg.norm(np.array(v['actual_traction'])-v['applied_traction']))
            for v in b['boundary_rows'] if v['role']=='outlet'),
        interface_max_normal=max(abs(v['traction_defect'][2]) for v in b['interface']),
        interface_max_tangential=max(abs(v['traction_defect'][0]) for v in b['interface']),
        epsilon_trace_max=max(abs(v['epsilon_trace_target'][0]-v['epsilon_trace_target'][1]) for v in b['interface']))
book=load(root/'accounting.json')
metrics=dict(schema='WJH3ClosureMetricsV1',compact_files=len(manifest['compact']),remote_heavy_files_verified=len(heavy),
    accepted_initial_and_stage_packets_verified=packets,diagnosis=diagnosis,paired_terminal_maxima=contrasts,boundary=boundary,
    current_charged_units=sum(v['units'] for v in book['entries'].values()),
    completed_or_reserved_seconds=sum(v.get('seconds',v['reserved_seconds']) for v in book['entries'].values())+book['prior_profile_setup_seconds'],
    closure_reserve_seconds=book['closure_seconds'],peak_recorded_RSS=max(v.get('peak_tree_rss_bytes',0) for v in book['entries'].values()),
    root_bytes=int(subprocess.check_output(['du','-sb',root],text=True).split()[0]),
    disclaimer='remote heavy verification; only compact artifacts are transported and independently read on Mac')
save(out/'metrics.json',metrics)
save(out/'verification.json',dict(status='PASS',manifest_sha256=sha(pub/'manifest.json'),
    heavy_index_sha256=sha(pub/'heavy-index.json'),compact_files=len(manifest['compact']),heavy_files=len(heavy),
    report_packet_bindings=packets,metrics_sha256=sha(out/'metrics.json'),governing_solves=0))
print(json.dumps({k:v for k,v in metrics.items() if k not in ('diagnosis','boundary','paired_terminal_maxima')},indent=2))
