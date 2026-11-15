#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Close H19 independently of its publication-time accounting snapshot."""
import argparse,json,tarfile
from pathlib import Path
from run_wj_h4 import load
from run_wj_h7 import update_journal
from run_wj_h15 import size
from wj_h15_observer import identity,verify,atomic

def close(root):
    book=load(root/'accounting.json');assert not book['closed']
    assert not any(e['status']=='RESERVED'for e in book['entries'].values())
    assert book['bootstrap_reserved']<=1 and book.get('target_blocks',0)<=1 and book.get('target_charged_attempts',0)<=2
    assert book.get('strict_reservations',[])==['target']
    pub=root/'publication-v1'
    for v in load(pub/'manifest.json')['files']:verify(v)
    gates=load(pub/'component-gates.json');old=load(root/'inputs/selected.json')['H18_closed_accounting'];verify(old)
    phase={k:sum(e.get('seconds',e['reserved_seconds'])for e in book['entries'].values()if e['phase']==k)for k in ('operation','build','qualification','target','target_readback','closure')}
    total=60+300+sum(e.get('seconds',e['reserved_seconds'])for e in book['entries'].values()if e['phase']!='closure')
    assert total<=5400 and phase['build']<=1200 and sum(phase[k]for k in ('qualification','target','target_readback'))<=1800
    source=sum(size(root/n)for n in ('source','build','runtime'))
    output=sum(size(p)if p.is_dir()else p.stat().st_size for p in root.iterdir()if p.name not in ('source','build','runtime'))
    peaks={kind:max([e.get('peak_tree_rss_bytes',0)for e in book['entries'].values()if(e['phase']=='build')==(kind=='build')]+[0])for kind in ('build','science')}
    assert source<=8*1024**3 and output<=2*1024**3 and peaks['build']<=16*1024**3 and peaks['science']<=8*1024**3
    closure=root/'closure-v1';closure.mkdir()
    runtime=load(root/'science/runtime-qualified-v1.json');assert runtime['status']=='PASS'
    for v in runtime['libraries']:verify(v)
    atomic(closure/'runtime-final-identity.json',runtime)
    counts=dict(coarse_native_CP0=book['bootstrap_reserved'],target_blocks=book.get('target_blocks',0),coarse_governing_attempts=book.get('target_charged_attempts',0),fine_governing_attempts=0,
        final_strict_readbacks=1,target_native_history_initializations=1,small_component_history_initializations=1,spatial_initializations=0,compiled_kernel_requests=0,
        focused_python_distinct_methods=18,registered_tests_not_execution=True,source_native='b2365e52',source_frozen_Python='4108e3da',new_coarse_scientific_attempts=book.get('target_charged_attempts',0),maximum_physical_horizon=.01)
    r=load(root/'science/coarse/case.evolution.json')
    counts.update(native_completed=r['native_completed'],completed_independently_qualified=r['independently_qualified'],governor_full_native_audits=r['full_audit_calls'],governor_PhysicalSamples=r['physical_samples_calls'],bootstrap_full_native_audits=1,
        independent_state_role_evaluations=8,independent_Q1_mass_projections=0,strict_diagnostic_mass_algebra='retained in native readback; zero governing solves; no separate diagnostic mass-system campaign')
    atomic(closure/'execution-counts.json',counts)
    atomic(closure/'verification.json',dict(status='PASS',component_gates=gates,phase_seconds=phase,total_charged_seconds=total,ceiling=5400,setup_fully_charged=60,closure_fully_charged=300,
        source_build_runtime_bytes=source,new_output_bytes=output,aggregate_peak_RSS=peaks,H18_closed_ledger_unchanged=old,publication=identity(pub/'manifest.json'),counts=counts,
        limits='provisional two-level screen only; production timestep/duration/application adequacy not admitted'))
    def seal(b):b.update(closed=True,total_charged_seconds=total,closure_fully_charged_seconds=300,final_verification=identity(closure/'verification.json'))
    update_journal(root/'accounting.json',seal)
    atomic(closure/'accounting-final.json',load(root/'accounting.json'))
    atomic(closure/'manifest.json',dict(files=[identity(p)for p in sorted(closure.iterdir())if p.is_file()],scope='separately closed accounting; publication snapshot remains OPEN'))
    atomic(root/'authoritative-index-v1.json',dict(publication=identity(pub/'manifest.json'),closure=identity(closure/'manifest.json'),accounting=identity(closure/'accounting-final.json'),verification=identity(closure/'verification.json')))
    selected=[root/'authoritative-index-v1.json',*closure.glob('*.json')]
    excluded={'final-strict-result.json.gz','reused-fine-terminal-strict.json.gz','copied-evidence.json'}
    selected.extend(p for p in pub.rglob('*')if p.is_file()and p.name not in excluded)
    archive=root/'essential-v1.tar.gz'
    with tarfile.open(archive,'w:gz')as t:
        for p in sorted(set(selected)):t.add(p,arcname=str(p.relative_to(root)))
    atomic(root/'essential-transport-intent-v1.json',dict(archive=identity(archive),files=[identity(p)for p in sorted(set(selected))],scope='compact essentials; optional raw strict views remain Linux; derived compact views preserve signed weak terms'))
    print(json.dumps(dict(index=identity(root/'authoritative-index-v1.json'),archive=identity(archive),resources=load(closure/'verification.json'))),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();close(a.root)
