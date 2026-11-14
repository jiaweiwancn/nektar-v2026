#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact publication of H15, including its retained protocol obstruction."""
import argparse
import csv
import gzip
import json
import os
from pathlib import Path
import shlex
import shutil
import time
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,verify,atomic,H14,load_cache

def table(path,records):
    if not records:return
    with path.open('x',newline='')as f:
        w=csv.DictWriter(f,list(records[0]));w.writeheader();w.writerows(records)

def objects(root,build,source_commit):
    records=[];dependencies={}
    for v in load(build/'compile_commands.json'):
        args=shlex.split(v['command']);p=build/args[args.index('-o')+1]
        if not p.exists():continue
        d=Path(str(p)+'.d');records.append(dict(source=identity(v['file']),object=identity(p),command=v['command'],dependency_file=identity(d)if d.exists()else None))
        if d.exists():
            for f in d.read_text().replace('\\\n',' ').split()[1:]:
                if f.startswith(str(root)+'/source')and Path(f).is_file():dependencies[f]=identity(f)
    return dict(native_source_commit=source_commit,cache=identity(build/'CMakeCache.txt'),compile_commands=identity(build/'compile_commands.json'),objects=records,local_dependencies=list(dependencies.values()),static_libraries=[identity(p)for p in sorted(build.glob('*.a'))],compiler=identity('/usr/bin/c++'))

def native_closure(root):
    historical=H14.parent/'wj-h10-flat-inlet-20261113';prior=load(historical/'runtime-qualified.json')['loaded']
    for v in load(historical/'closure/runtime-observed.json')['libraries']:prior[v['path']]=v
    files=set();samples=[]
    for op in ('small-native-protocol-controls-v1','small-genuine-strict-view-v1','single-four-stage-target-v1','last-qualified-CP0-strict-readback-v1'):
        pin=identity(root/'operations'/op/'native-exec-samples.json');records=load(verify(pin));assert records
        for r in records:
            assert r['status']=='PASS'and r['libc']==[str(root/'runtime/glibc/libc.so.6')]
            files.update(f for f in r['files']if '.so'in Path(f).name or Path(f).name.startswith('ld-linux'))
        samples.append(dict(operation=op,saved_maps=pin))
    libraries=[]
    for f in sorted(files):
        key=f
        if '/runtime/glibc/'in f:key=str(historical/'runtime/glibc'/Path(f).name)
        elif '/runtime/geometry/'in f:key=str(historical/'runtime/geometry'/Path(f).name)
        assert key in prior,f;actual=identity(f);p=prior[key];assert actual['sha256']==p['sha256']and actual['bytes']==p['bytes']
        libraries.append(dict(actual=actual,qualified=p))
    return dict(status='PASS',actual_native_exec_maps=samples,libraries=libraries,governor_and_reader=load(root/'science/runtime-build-v1.json')['executables'],scope='original H15 runtime actually executed; no execution of repaired governor/reader')

def main(root,source):
    out=root/'publication-v1';out.mkdir();(out/'controls').mkdir();(out/'operations').mkdir();(out/'analysis').mkdir()
    def copy(p,name):
        p=Path(p);dest=out/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest);assert identity(p)['sha256']==identity(dest)['sha256']
    # Save snapshots before hashing; publication snapshot remains OPEN and is
    # never relabelled as the later final CLOSED journal.
    copy(root/'accounting.json','accounting-at-publication.json')
    aliases=load(root/'source-program-aliases-v1.json')
    copy(root/'source-program-aliases-v1.json','source-program-aliases.json')
    for v in load(root/'target-admission.json')['inputs']:
        if identity(v['path'])==v:continue
        selected=[a for a in aliases['aliases']if a['original']==v]
        assert len(selected)==1,'changed target input without an exact source alias: '+v['path']
        p=verify(selected[0]['retained']);assert identity(p)['sha256']==v['sha256']
        assert v['path'].endswith('/study_wj_h15.py'),'scientific inputs cannot be substituted through source aliases'
    for p in ('target-admission.json','target-final-strict-admission.json','protocol-repair-source-v4.json','bounded-stream-source-v6.json',
        'science/imports.json','science/runtime-build-v1.json','science/runtime-qualified-v1.json','science/small/qualification.json',
        'science/small/pass/view/result.json','science/final-qualified-cp0-view/result.json','science/protocol-repair-checks-v2/result.json',
        'science/observer-cache/cp0-parity.json','science/observer-cache/pilot1-parity.json','science/target/native-outcome.json','science/target/case.evolution.json',
        'science/observer-cache/manifest.json','science/small-cache-v2/manifest.json'):
        copy(root/p,'controls/'+p.replace('/','--'))
    for p in sorted((root/'science/analysis-v1').iterdir()):
        if p.is_file():copy(p,'analysis/'+p.name)
    for p in sorted((root/'science/target/monitor').glob('step-1*')):
        if p.suffix=='.json'and 'inventory'not in p.name and not p.name.endswith(('.report.json','.callback.json')):copy(p,'controls/target-'+p.name)
    began=time.monotonic();cache,cache_manifest=load_cache(root/'science/observer-cache/manifest.json')
    atomic(out/'cache-load-benchmark.json',dict(seconds=time.monotonic()-began,role='fresh passive hash-verification/deserialization benchmark; not a measurement of target observer startup',base_rows=len(cache.base.rows),base_owned_rows=cache.base.point_count,basis_keys=len(cache.index.basis_values),cache=cache_manifest['cache'],rebuilds=0,governing_work=0))
    del cache
    stage=load(root/'science/target/monitor/step-1.report.json');live=load(root/'science/target/monitor/step-1.ready.result.json');summary=load(root/'science/target/case.evolution.json')
    assert summary['native_completed']==1 and summary['independently_qualified']==0 and all(live['replay'].values())
    assert not(root/'science/target/monitor/qualified-cursor.json').exists()
    copy(root/'science/target/execution-diff.json','controls/target-execution-diff.json')
    copy(root/'science/target/case.xml','controls/target-case.xml')
    copy(root/'science/target/intent.json','controls/target-monitor-intent.json')
    atomic(out/'component-gates.json',dict(implementation_build='PASS',bounded_retention_small_controls='PASS',retained_CP0_CP1_observer_parity='PASS',small_native_attempts=10,
        small_strict_view='PASS',target_native_attempts=1,target_reserved_maximum_attempts=4,target_native_completed=1,target_replay_bitwise='PASS',target_native_hard_admission='PASS',target_independent_physical_quality='PASS',
        target_protocol='OBSTRUCTED',target_independently_qualified=0,target_new_scientific_progress=0,target_stages_2_to_4='NOT_RUN',terminal_strict_original_qualified_CP0='PASS',completed_replay_strict_eligibility='NOT_QUALIFIED',
        repaired_protocol_synthetic_regression='PASS',repaired_native_build='PASS',repaired_governor_runtime_execution='NOT_RUN',post_repair_packet_fields_runtime_qualification='NOT_RUN',warm_multistep_route='NOT_QUALIFIED',
        duration_preservation='NOT_RUN',temporal_adequacy='UNRESOLVED',physical_height_selection='NOT_SELECTED',experimental_hill_validation='NOT_RUN'))
    bindings=dict(starting_commit='fa5a48dfe6d82ab6a58efd70e7332636701e20ae',publication_source_commit=source,native_executed_source_commit='2ffd38f061073d17b4b595a8795f31813606862f',
        source_runtime=identity(root/'science/runtime-build-v1.json'),original_qh0=load(root/'science/observer-cache/manifest.json')['reference_identity'],parent_publication=identity(H14/'publication-v2/manifest.json'),
        completed_but_unqualified=dict(producer=identity(root/'science/target/monitor/step-1.report.json'),checkpoint=identity(Path(summary['last_completed_checkpoint'])/'State.bin'),ready=identity(root/'science/target/monitor/step-1.ready.json'),callback=identity(root/'science/target/monitor/step-1.callback.json'),state_sha256=stage['wj_audit_packets'][0]['state_sha256'],callback_history_sha256=stage['wj_audit_packets'][0]['stage_history_sha256']),
        last_qualified=dict(checkpoint=identity(Path(summary['last_qualified_checkpoint'])/'State.bin'),absolute_step=0,time=0.,producer_role='original H14 CP0; no completed callback'),
        failure=summary['failure'],replay=live['replay'],strict_readbacks=dict(small=1,target_last_qualified_CP0=1,completed_H15_checkpoint=0),fresh_CP0s=0,initializer_applications=0)
    atomic(out/'state-bindings.json',bindings);atomic(out/'runtime-verification.json',native_closure(root))
    atomic(out/'compiled-source-objects.json',dict(executed=objects(root,root/'build','2ffd38f061073d17b4b595a8795f31813606862f'),protocol_repair=objects(root,root/'build-protocol-v5','f69dce2b'),complete_bounded_streams=objects(root,root/'build-bounded-v6','90e28ce3'),execution_after_repair='NOT_RUN for governor/reader; synthetic unit only'))
    attempts=[]
    for case in ('legacy','pass','restart','observer-fail','timeout','stale','native-reject','post-callback'):
        p=root/'science/small'/case/'case.evolution.json';r=load(p)
        attempts.append(dict(case=case,governing_attempts=2 if case in ('legacy','pass')else 1,native_completed=2 if case in ('legacy','pass')else 1 if case in ('restart','observer-fail','timeout','stale')else 0,independent_committed=r.get('independently_qualified'),failure=r.get('failure',''),producer=identity(p)))
    attempts.append(dict(case='target_replay_only',governing_attempts=1,native_completed=1,independent_committed=0,failure=summary['failure'],producer=identity(root/'science/target/case.evolution.json')))
    atomic(out/'attempt-counts.json',dict(total_governing_attempts=11,small_verification_attempts=10,target_replay_attempts=1,target_new_scientific_steps=0,reserved_unlaunched_target_stages=3,records=attempts,unit_catalog_calls='synthetic, excluded from governing counts'))
    atomic(out/'replay-linear-and-stage-defects.json',dict(source=identity(root/'science/target/monitor/step-1.report.json'),frozen_subsolve=stage['pimple_attempts'][0],unsplit_callback_stage=stage['stage_audits'],
        stage_columns=['time','lambda','alpha','iterations','maximum_native_scaled_residual','anchor_energy_work','native_scaled_field_residuals...'],
        constraints='existing PIMPLE true linear/Schur/weak gates unchanged; unsplit single-sweep defect is reported, not required to vanish',
        new_H15_completed_stage_strict_mass_dual='NOT_RUN: replay checkpoint unqualified by protocol',reused_H14_CP1_strict=identity(H14/'publication-v2/pilot1-strict.json')))
    # Full native call streams are separately sealed, never inferred from one
    # last in-memory factor allocation or an incidental successful exit code.
    calls=[json.loads(v)for v in (root/'science/target/monitor/factors/calls.jsonl').read_text().splitlines()]
    assert [v['call']for v in calls]==list(range(1,summary['factor_catalog']['calls']+1))
    for v in calls:assert max(0,v['resident_definitions'])<=16 and v['resident_key_bytes']<=32*1024**2
    atomic(out/'retention-and-factor-calls.json',dict(summary=summary,streamed_factor_calls=calls,complete_pattern_definitions=len(list((root/'science/target/monitor/factors').glob('factor-*.json'))),
        stream_roles='implicit attempt remains PROVISIONAL_NOT_TIME_COMMIT; ready packet is completed; prepared PASS is not native qualification commit',post_repair_added_packet_fields='compiled, not retrospectively inserted into original packets'))
    logs=[];timing=[]
    for folder in sorted((root/'operations').iterdir()):
        dest=out/'operations'/folder.name;dest.mkdir()
        for name in ('request.json','receipt.json','receipt-recovered.json','native-exec-samples.json'):
            if(folder/name).exists():copy(folder/name,'operations/'+folder.name+'/'+name)
        if(folder/'command.log').exists():
            with(folder/'command.log').open('rb')as f,gzip.open(dest/'command.log.gz','wb')as g:shutil.copyfileobj(f,g)
            logs.append(dict(raw=identity(folder/'command.log'),compressed=identity(dest/'command.log.gz')))
        if(folder/'timing.jsonl').exists():
            for line in (folder/'timing.jsonl').read_text().splitlines():
                v=json.loads(line);timing.append(dict(operation=folder.name,**v))
    for name in ('native.log','observer.log'):
        p=root/'science/target'/name
        with p.open('rb')as f,gzip.open(out/('target-'+name+'.gz'),'wb')as g:shutil.copyfileobj(f,g)
        logs.append(dict(raw=identity(p),compressed=identity(out/('target-'+name+'.gz'))))
    atomic(out/'compressed-log-bindings.json',logs);atomic(out/'native-timing.json',timing)
    pid=load(root/'science/target/native-launch.json')['pid'];rows=[v for v in timing if v['pid']==pid]
    table(out/'target-native-disjoint-and-inclusive-timing.csv',[{k:v.get(k)for k in ('pid','step','time','name','calls','wall_seconds','exclusive_wall_seconds','process_cpu_seconds','clock_scope')}for v in rows])
    oldrate=live['native_compute_seconds']+live['measured_independent_seconds']
    atomic(out/'conditional-cost-and-limits.json',dict(native_replay_compute_seconds=live['native_compute_seconds'],independent_replay_assessment_seconds=live['measured_independent_seconds'],immutable_cache_build_seconds=load(root/'science/observer-cache/manifest.json')['build_seconds'],
        future_horizon=505.02466521278404,terminal_window=[447.8483832973273,505.02466521278404],dt=.0025,steps=202010,replay_rate_only_conditional_seconds=202010*oldrate,replay_rate_only_conditional_days=202010*oldrate/86400,
        measured_new_warm_scientific_steps=0,warm_rate='UNAVAILABLE: protocol stopped before new progress',duration_admission='CLOSED',additional_costs='cold startup, strict views, any future witness growth, I/O and error/temporal qualification are not established by this extrapolation',
        comparisons='native.step remains inclusive of publication/wait; ready.native_compute_seconds ends before them; accumulated audit buckets are disjoint wall only, CPU not measured'))
    heavy=[identity(p)for p in sorted((root/'science').rglob('*'))if p.is_file()]
    atomic(out/'remote-heavy-manifest.json',dict(files=heavy,bytes=sum(v['bytes']for v in heavy),scope='new H15 science root only; original H14 parent data referenced, not duplicated'))
    atomic(out/'control-commands.json',dict(source_revision=source,root=str(root),runtime='original qualified private loader/library closure mandatory',ledger='OPEN publication-time snapshot; final closure is separate',forbidden='no native retry; historical operation commands are evidence, not permission to repeat them'))
    artifacts=[identity(p)for p in sorted(out.rglob('*'))if p.is_file()]
    atomic(out/'manifest.json',dict(source_revision=source,artifacts=artifacts,heavy_data_copied=False))
    atomic(out/'seal.json',dict(manifest=identity(out/'manifest.json'),accounting_snapshot=identity(out/'accounting-at-publication.json'),status='OBSTRUCTED_PROTOCOL; partial scientific pilot preserved'))
    print('publication-v1 sealed',len(artifacts),'compact artifacts',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--source',required=True);a=p.parse_args();main(a.root,a.source)
