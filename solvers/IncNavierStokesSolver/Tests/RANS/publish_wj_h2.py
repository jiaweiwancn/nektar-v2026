#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact, read-only H2 publication; heavy native reports stay in the run root."""
import argparse
import hashlib
from pathlib import Path
import shlex
import shutil
import subprocess
import time
import numpy as np
from run_wj_h2 import load,save,sha
from compare_cls_epawf_pimple import same_history
from wj3_measure import checkpoint_payload

PRODUCER='c344dc22'  # resolved against the preserved Git object below


def identity(path):
    p=Path(path).resolve();return dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p))


def git(source,*args):
    return subprocess.check_output(['git','-C',str(source),*args],text=True).strip()


def native_inputs(root,source,revision):
    inputs=set()
    for target in ('IncNavierStokesSolver','RANSKernels'):
        for dep in (root/'build-v1/CMakeFiles'/f'{target}.dir').rglob('*.o.d'):
            for token in shlex.split(dep.read_text().replace('\\\n',' '))[1:]:
                p=Path(token)
                try:p.relative_to(source)
                except ValueError:continue
                if p.is_file():inputs.add(p)
    assert inputs,'no native compile dependencies found'
    result=[]
    for p in sorted(inputs):
        rel=str(p.relative_to(source))
        blob=subprocess.check_output(['git','-C',str(source),'show',revision+':'+rel])
        expected=hashlib.sha256(blob).hexdigest();assert expected==sha(p),'changed producing dependency: '+rel
        result.append(dict(path=rel,sha256=expected))
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();start=time.monotonic();r=a.root.resolve();out=a.output.resolve();out.mkdir(exist_ok=False);s=r/'source-v1'
    revision=git(s,'rev-parse',PRODUCER);head=git(s,'rev-parse','HEAD');tree=git(s,'rev-parse','HEAD^{tree}')
    source_inputs=native_inputs(r,s,revision);save(out/'native-source-inputs.json',source_inputs)
    manifest=dict(schema='WJH2CompactPublicationV1',root=str(r),source=dict(publication=head,tree=tree,
        native_producer=revision,native_producer_tree=git(s,'rev-parse',revision+'^{tree}'),
        scientific_driver=git(s,'rev-parse','c2189f8'),half_dt_driver=git(s,'rev-parse','016c25d'),
        observer=git(s,'rev-parse','e62b5f4'),native_input_count=len(source_inputs),native_inputs_unchanged_at_publication=True),
        status='OBSTRUCTED',developed_precursor='NOT_RUN',denied_continuation='NOT_RUN',
        heavy_verification='bytes read and SHA256 independently computed on Linux; transported compact bytes checked separately on Mac',
        arithmetic=dict(compiler=subprocess.check_output(['g++','--version'],text=True).splitlines()[0],
            flags='Release -O3 -DNDEBUG -ffp-contract=off -fno-fast-math; C++17; native targets retain strict warnings',
            RANS_BUILD_SOLVER=True,rank=1,numerical_threads=1),
        files=[],observations={},executions={},runtimes={},checkpoints=[],transactions={})
    specifications={'wj_earsm_neutral_abl_ep_awf_note.tex':'1959e96ea0b61f9deb9ff9ad1f5b26488a932a8451c2de6d3b6f7b7c9dbd9124',
        'wj_earsm_budget_cutoff_awf_appendix.tex':'0ab52e0d813b8fb344973f75ca6d8da97ff389ef9f93476500ca9d1df1fbc47c'}
    manifest['specifications']={}
    for name,expected in specifications.items():
        v=identity(s/'docs/literature'/name);assert v['sha256']==expected;manifest['specifications'][name]=v
    for p in sorted((r/'operations').glob('*/receipt.json')):
        d=load(p);label=p.parent.name
        manifest['executions'][label]={k:d.get(k) for k in ('argv','cwd','status','units','seconds','exit','reason','peak_tree_rss_bytes','peak_tree_threads','output_bytes','log_sha256','inputs','loader')}
        manifest['files'].extend(identity(v) for v in (p,p.parent/'request.json',p.parent/'command.log') if v.exists())
        timing=p.parent/'timing.jsonl'
        if timing.exists():manifest['files'].append(identity(timing))
    for p in sorted((r/'receipts').iterdir()):
        # The current stdout/time receipt is still open. Its closed bytes are
        # bound by the post-publication verification bridge, never prehashed.
        if p.is_file() and p.name!=out.name+'.log':manifest['files'].append(identity(p))
    for p in (r/'build-v1/compile_commands.json',r/'build-v1/CMakeCache.txt',r/'build-v1/CMakeFiles/IncNavierStokesSolver.dir/flags.make',
        r/'build-v1/CMakeFiles/IncNavierStokesSolver.dir/link.txt',r/'loaded-runtime-v2/summary.json'):
        manifest['files'].append(identity(p))
    for role,p in [('governing',r/'runtime-v4/IncNavierStokesSolver'),('initial_observer',r/'runtime-v4/RANSWallResolutionProbe'),
        ('versioned_observer',r/'observer-v6/RANSWallResolutionProbe'),('kernel',r/'build-v1/WJKernelProbe')]:
        manifest['runtimes'][role]=identity(p)
    loaded=load(r/'loaded-runtime-v2/summary.json');assert loaded['status']=='PASS' and loaded['governing_attempts']==0
    assert loaded['solver_sha256']==manifest['runtimes']['governing']['sha256']
    save(out/'loaded-runtime.json',loaded)
    study=load(r/'study-v1/progress.json');assert study['status']=='OBSTRUCTED'
    runs=dict(study['runs']);runs['refined-half-dt-diagnostic']=load(r/'diagnosis-v1/result.json')
    for key,v in runs.items():
        report=load(v['report']);assert sha(v['report'])==v['report_sha256']
        manifest['files'].extend(identity(p) for p in (Path(v['report']),Path(v['xml']),Path(v['xml']).with_name('case-design.json')))
        for cp in sorted(Path(v['xml']).parent.glob('CP*')):
            if not (cp/'State.bin').is_file():continue
            payload=checkpoint_payload(cp);entry=identity(cp/'State.bin')
            entry.update(run=key,step=payload['steps'],time=payload['time'],dt=payload['dt'],reference_sha256=payload['reference_sha256'])
            manifest['checkpoints'].append(entry)
        for attempt in report['pimple_attempts']:
            if attempt.get('k_snapshot'):manifest['files'].append(identity(attempt['k_snapshot']))
    for p in (r/'study-v1/progress.json',r/'study-v1/intent.json',r/'diagnosis-v1/result.json',r/'diagnosis-v1/intent.json',
        r/'observations-v1/intent.json',r/'observations-v2/intent.json',r/'observations-v3/intent.json',r/'observations-v3/index.json',
        r/'observations-v3/summary.json',r/'observations-v3/incident.json',r/'reader-fixtures-v3/summary.json'):
        manifest['files'].append(identity(p))
    for label,record in load(r/'observations-v3/index.json').items():
        assert sha(record['native'])==record['native_sha256'];d=load(record['native'])
        assert d['report_sha256']==record['report_sha256'] and d['checkpoint_sha256']==record['checkpoint_sha256']
        assert d['callback_history_sha256']==record['callback_history_sha256']
        manifest['observations'][label]=record
        manifest['files'].append(identity(Path(record['native']).with_name('request.json')))
        weights=Path(record['native']).with_name('weights.json')
        if weights.exists():manifest['files'].append(identity(weights))
    control=r/'build-v1/wj-h2-open-controls/run-931bb06ba09f/summary.json';c=load(control)
    reports={k:load(v['report']) for k,v in c['runs'].items()}
    for key in ('restart','reuse-on'):same_history(reports['two'],reports[key])
    same_history(reports['one'],reports['frozen']);same_history(reports['one'],reports['rejected'])
    for key in ('physical_samples','state_audits','wj_audit_packets'):
        assert reports['two'][key]==reports['reuse-on'][key]
        assert reports['frozen'][key]==reports['rejected'][key]
    manifest['transactions']=dict(summary=identity(control),signed_zero_preserving_replay='PASS',accepted_native_steps=6,
        rejected_complete_attempts=1,one_U_pressure_k_epsilon_bound_each=True,unmodified_producer=c['unchanged_producer'])
    for v in c['runs'].values():manifest['files'].append(identity(v['report']))
    for pattern in ('wj2-native-Budget/**/summary.json','wj2-legacy-transaction/**/summary.json','wj1-oracle.json'):
        for p in (r/'build-v1').glob(pattern):manifest['files'].append(identity(p))
    for name in ('summary.json','profiles-at-measured-heights.csv','experimental-discrepancy.csv','comparison-table.csv',
        'physical-profiles.png','admission-and-interface.png'):
        shutil.copyfile(r/'analysis-v2'/name,out/name)
    shutil.copyfile(r/'completions-compiled-v1/summary.json',out/'inlet-completion-family.json')
    shutil.copyfile(r/'observations-v3/incident.json',out/'principal-incident.json')
    shutil.copyfile(s/'docs/plans/wj-h2-linux/continuation-denial.json',out/'continuation-denial.json')
    inventory=load(r/'receipts/ctest-inventory-final-v1.json');selection=load(r/'receipts/regression-selection.json')
    manifest['tests']=dict(registrations=len(inventory['tests']),selected=selection['selected'],selected_count=20,
        selected_status=load(r/'operations/selected-regressions-v2/receipt.json')['status'],
        separate_open_boundary_control_status=load(r/'operations/open-boundary-controls-v2/receipt.json')['status'],
        additional_passive_reader_status=load(r/'reader-fixtures-v3/summary.json')['status'])
    book=load(r/'accounting.json');shutil.copyfile(r/'accounting.json',out/'ledger.json')
    manifest['accounting']=dict(charged_units=sum(v['units'] for v in book['entries'].values()),
        charged_seconds=sum(v.get('seconds',v['reserved_seconds']) for v in book['entries'].values())+book['prior_profile_setup_seconds'],
        closure_reserve=book['closure_seconds'],unit_ceiling=book['unit_ceiling'],seconds_ceiling=book['seconds_ceiling'],
        maximum_recorded_global_workers=max(v.get('allocation',{}).get('workers',0) for v in book['entries'].values()),
        ledger=identity(r/'accounting.json'),build_time_scope='separate /usr/bin/time receipts; not numerical-study seconds')
    assert manifest['accounting']['charged_units']<=2200 and manifest['accounting']['charged_seconds']+600<=21600
    save(out/'manifest.json',manifest)
    receipt=dict(schema='WJH2ReceiptV1',source=head,manifest=identity(out/'manifest.json'),
        payloads={str(p.relative_to(out)):identity(p) for p in sorted(out.iterdir()) if p.is_file()},
        publication_seconds_before_final_receipt=time.monotonic()-start,
        root_bytes_before_publication_close=int(subprocess.check_output(['du','-sb',str(r)],text=True).split()[0]),
        status='OBSTRUCTED',continuation='NOT_RUN: automatic review prelaunch denial',stationarity='NOT_REACHED / no developed window')
    save(out/'receipt.json',receipt);print('PUBLISHED',str(out),'OBSTRUCTED')


if __name__=='__main__':main()
