#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Seal H18 accounting and essential transport without admitting another solve."""
import argparse,json,tarfile,os,shutil,re
from pathlib import Path
from run_wj_h4 import load
from run_wj_h7 import update_journal
from run_wj_h15 import size
from wj_h15_observer import identity,verify,atomic


def close(root):
    book=load(root/'accounting.json');assert not book['closed']
    assert not any(v['status']=='RESERVED'for v in book['entries'].values())
    pub=root/'publication-v1';manifest=load(pub/'manifest.json')
    for v in manifest['files']:verify(v)
    gates=load(pub/'component-gates.json')
    # Partial qualification remains publishable, never relabelled successful.
    if gates['scope_result']=='PASS':
        assert book['small_charged_attempts']==6 and book['target_reserved_attempts']==book['target_charged_attempts']==6
        assert book['target_cases']==['A','B']and book['strict_reservations']==['target']
        assert gates['successful_small_attempts']==5 and gates['deliberate_small_failed_attempts']==1
    old=load(root/'inputs/selected-bootstrap.json')['H16_ledger_unchanged'];verify(old)
    supplement=root/'publication-supplement-v1';supplement.mkdir()
    for name in ('post-science-code-bindings-v1.json','post-science-code-bindings-v2.json','readback-selection-compatibility-v1.json','pending-restored-native-map-parity-v1.json','target-B-measured-admission-v1.json'):
        p=root/'science'/name
        if p.exists():shutil.copyfile(p,supplement/name)
    for label in ('passive-readback-selection-v1','compact-publication-v1'):
        p=root/'operations'/label/'receipt.json'
        if p.exists():shutil.copyfile(p,supplement/(label+'-receipt.json'))
    runs=[load(root/f'science/target/{n}/case.evolution.json')for n in ('A','B')if(root/f'science/target/{n}/case.evolution.json').exists()]
    initial_log=(root/'operations/focused-python-v1/command.log').read_text();post_log=(root/'operations/passive-readback-selection-v1/command.log').read_text()
    initial_count=int(re.search(r'Ran (\d+) tests',initial_log).group(1));post_count=int(re.search(r'Ran (\d+) tests',post_log).group(1))
    starts=sum((root/f'science/target/{n}/monitor/startup.verification-commit.json').exists()for n in ('A','B'))
    atomic(supplement/'execution-counts.json',dict(focused_python_methods=initial_count+post_count,initial_methods=initial_count,post_science_selector_methods=post_count,native_protocol_unit_invocations=1,native_guard_pending_import_component_invocations=1,
        successful_small_governing_attempts=gates['successful_small_attempts'],expected_failed_small_attempts=gates['deliberate_small_failed_attempts'],target_attempts=book.get('target_charged_attempts',0),target_blocks=len(runs),target_completed=sum(r['native_completed']for r in runs),
        target_completed_scientific_native_commits=sum(r['independently_qualified']for r in runs),target_startup_verification_commits=starts,target_full_native_audits=sum(r['full_audit_calls']for r in runs),target_physical_samples=sum(r['physical_samples_calls']for r in runs),target_independent_assessments=len(gates['full_audited_states']),
        skipped_target_completed_full_audits=sum(r['native_completed']-r['independently_qualified']for r in runs),final_strict_readbacks=len(book.get('strict_reservations',[])),new_scientific_timesteps=0,fresh_CP0=0,initialization_projection_lift_calls=0,new_kernel_requests=0,
        scope='actual H18 execution, not CTest registration; intermediate H16 qualification remains historical'))
    atomic(supplement/'retention-and-limits.json',dict(native_resident_limits=dict(pimple=1,accepted=1,WJ=1,stage=1,step=1,explicit=2),
        independent_scope='fixed sealed common support plus complete current checkpoint witnesses; monotone schedule, no historical ready/seen union',
        recovery_rotation='owned/protected/eligible deletion API tested on synthetic files; H18 deliberately retains every short-run recovery',
        late_substep_evidence='captured at available k-bound events/selected steps; discarded operators are not reconstructed',
        native_ACK_validation='existing final-name stale/duplicate checks preserved; no claim of unrestricted duration/namespace cost',
        observation_limits='no new averages; sparse checkpoint samples cannot recover unsaved extrema, exact integrated budgets or nonlinear stress averages',
        accounting='publication snapshot OPEN; final closure is separate; original mutable accounting path in copy records describes its bytes at copy time only'))
    atomic(supplement/'manifest.json',dict(files=[identity(p)for p in sorted(supplement.iterdir())if p.is_file()]))

    phase={k:sum(v.get('seconds',v['reserved_seconds'])for v in book['entries'].values()if v['phase']==k)for k in ('build','qualification','target','target_readback','operation','closure')}
    assert phase['build']<=1200 and phase['qualification']+phase['target']+phase['target_readback']<=3000 and phase['closure']<=300
    total=60+300+sum(v.get('seconds',v['reserved_seconds'])for v in book['entries'].values()if v['phase']!='closure')
    assert total<=7200
    source_names=('source','build','runtime','prequalification-build-v1','prequalification-build-v2')
    src=sum(size(root/n)for n in source_names if(root/n).exists())
    output=sum(size(p)if p.is_dir()else p.stat().st_size for p in root.iterdir()if p.name not in source_names)
    peaks={kind:max([v.get('peak_tree_rss_bytes',0)for v in book['entries'].values()if(v['phase']=='build')==(kind=='build')]+[0])for kind in ('build','science')}
    assert src<=8*1024**3 and output<=2*1024**3 and peaks['build']<=16*1024**3 and peaks['science']<=8*1024**3
    closure=root/'closure-v1';closure.mkdir()
    expected={v['ELF']['path']:v['ELF']for v in load(root/'science/runtime-build-v1.json')['executables']}
    deps={v['actual']['path']:v['actual']for v in load(root/'science/runtime-dependencies-v1.json')['libraries']}
    samples=[];used={}
    for path in sorted((root/'operations').glob('*/native-exec-samples.json')):
        for record in load(path):
            assert record['status']=='PASS'
            for native in record['native']:verify(expected[native])
            for f in record['files']:
                if '.so'not in Path(f).name and not Path(f).name.startswith('ld-linux'):continue
                assert f in deps,f
                if f not in used:used[f]=identity(verify(deps[f]))
            samples.append(dict(operation=path.parent.name,pid=record['pid'],native=record['native'],libc=record['libc']))
    atomic(closure/'runtime-final-identity.json',dict(status='PASS',actual_native_exec_samples=samples,libraries=list(used.values()),new_ELFs=list(expected.values()),pretarget_qualification=identity(root/'science/runtime-qualified-v1.json')))

    atomic(closure/'verification.json',dict(status='PASS',scientific_scope=gates['scope_result'],phase_seconds=phase,total_charged_seconds=total,ceiling=7200,setup_fully_charged=60,closure_fully_charged=300,
        source_build_runtime_bytes=src,new_output_bytes=output,aggregate_peak_RSS=peaks,H16_ledger_unchanged=old,publication=identity(pub/'manifest.json'),new_scientific_steps=0))
    def seal(b):b.update(closed=True,total_charged_seconds=total,closure_fully_charged_seconds=300,final_verification=identity(closure/'verification.json'))
    update_journal(root/'accounting.json',seal)
    atomic(closure/'accounting-final.json',load(root/'accounting.json'))
    atomic(closure/'manifest.json',dict(files=[identity(p)for p in sorted(closure.iterdir())if p.is_file()],scope='final closed accounting; publication-time snapshot retains its original OPEN identity'))
    atomic(root/'authoritative-index-v1.json',dict(publication=identity(pub/'manifest.json'),supplement=identity(supplement/'manifest.json'),closure=identity(closure/'manifest.json'),accounting=identity(closure/'accounting-final.json'),verification=identity(closure/'verification.json')))
    selected=[root/'authoritative-index-v1.json',*closure.glob('*.json'),*supplement.glob('*.json')]
    essentials={'component-gates.json','target-timeline.csv','disjoint-inclusive-timing.csv','checkpoint-storage.csv','conditional-cost-storage.json','compression-measurements.json','state-history-bindings.json','all-signed-Q1-moments.csv','compiled-source-summary.json','checkpoint-cost.png','checkpoint-storage.png','heavy-manifest.json','manifest.json','reused-context.json','final-strict-result.json.gz','audited-state-summaries.json.gz','signed-native-stage-data.json.gz'}
    for p in pub.rglob('*'):
        if p.is_file()and(p.name in essentials or p.parent.name=='controls'or p.parent.name=='small'or'operations'in p.parts and p.name in ('receipt.json','attempt-counts.json')):selected.append(p)
    archive=root/'essential-v1.tar.gz'
    with tarfile.open(archive,'w:gz')as t:
        for p in sorted(set(selected)):t.add(p,arcname=str(p.relative_to(root)))
    atomic(root/'essential-transport-intent-v1.json',dict(archive=identity(archive),files=[identity(p)for p in sorted(set(selected))],closure_fully_charged=300,scope='bounded essential subset; original full publication and all heavy arrays remain Linux'))
    print(json.dumps(dict(index=identity(root/'authoritative-index-v1.json'),archive=identity(archive),total_charged_seconds=total,resources=load(closure/'verification.json'))),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();close(a.root)
