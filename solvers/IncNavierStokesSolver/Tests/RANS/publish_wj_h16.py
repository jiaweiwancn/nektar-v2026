#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Small authoritative H16 publication. Heavy coefficients remain in place."""
import argparse
import csv
import json
import os
from pathlib import Path
import re
import shutil
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import atomic,identity,verify
from publish_wj_h15 import immutable_operation

def native_closure(root):
    build=load(root/'science/runtime-build-v1.json');natives={v['ELF']['path']:v['ELF']for v in build['executables']};deps={v['actual']['path']:v['actual']for v in load(root/'science/runtime-dependencies-v1.json')['libraries']};samples=[];libraries={}
    for p in sorted((root/'operations').glob('*/native-exec-samples.json')):
        for v in load(p):
            assert v['status']=='PASS'and v['libc']==[str(root/'runtime/glibc/libc.so.6')]
            for n in v['native']:verify(natives[n])
            for f in v['files']:
                if '.so'in Path(f).name or Path(f).name.startswith('ld-linux'):libraries[f]=identity(verify(deps[f]))
            samples.append(dict(operation=p.parent.name,pid=v['pid'],native=v['native'],libc=v['libc']))
    return dict(status='PASS',actual_native_exec_samples=samples,libraries=list(libraries.values()),compiled_build=identity(root/'science/runtime-build-v1.json'),scope='all actual H16 governor/reader executions, including target; pre-exec host mappings excluded')

def publish(root,revision):
    out=root/'publication-v1';out.mkdir();book=load(root/'accounting.json')
    atomic(out/'accounting-publication-snapshot.json',book)
    copied=[]
    def copy(p,name):
        p=Path(p)
        if not p.exists():return None
        dest=out/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest);assert identity(p)['sha256']==identity(dest)['sha256'];copied.append(dict(original=identity(p),copy=identity(dest)));return dest
    for n in ('source-intent.json','execution-config.json','target-admission.json','inputs/essential-index.json','inputs/passive-source-v2.json','inputs/h15-final/publication-v3/manifest.json','inputs/h15-final/closure-v3/accounting-final.json'):
        copy(root/n,'provenance/'+n.replace('/','--'))
    for n in ('runtime-build-v1.json','runtime-dependencies-v1.json','runtime-qualified-v1.json','compiled-source-objects.json','target-cache-compatibility.json','small-cache-compatibility.json','scientific-imports-v2.json','retained-parity/result.json','retained-parity/cp0.json','retained-parity/pilot1.json'):
        copy(root/'science'/n,'controls/'+n.replace('/','--'))
    for label in book['entries']:
        copy(root/'operations'/label/'receipt.json','operations/'+label+'.json')
        if immutable_operation(book,label):copy(root/'operations'/label/'command.log','operations/'+label+'.log')
        copy(root/'operations'/label/'native-exec-samples.json','operations/'+label+'-native-maps.json')
    tests=load(root/'operations/focused-python-controls-v1/receipt.json');log=(root/'operations/focused-python-controls-v1/command.log').read_text();count=int(re.search(r'Ran (\d+) tests',log)[1])
    atomic(out/'native-runtime-observed.json',native_closure(root))
    small=[]
    for n in 'ABCDE':
        p=root/f'science/small/{n}/qualification.json'
        if p.exists():small.append(load(p));copy(p,f'small/{n}.json')
        for f in ('outcome.json','execution-diff.json','case.evolution.json'):
            copy(root/f'science/small/{n}'/f,f'small/{n}-'+f.replace('/','--'))
        p=root/f'science/small/{n}/monitor/first-aborted.json'
        if p.exists():
            v=load(p);atomic(out/f'small/{n}-first-aborted-summary.json',dict(source=identity(p),**{k:v[k]for k in ('status','failure','attempted_step','candidate_binary','callback_history_binary','actual_explicit_calls')},provisional_counts={k:len(v[k])for k in ('provisional_audits','provisional_stage_audits','provisional_explicit_audits')},first_native_range_event_present=v['first_native_range_event']is not None))
    copy(root/'science/small/A/view/request.json','small/strict-request.json')
    if(root/'science/small/A/view/result.json').exists():
        v=load(root/'science/small/A/view/result.json');atomic(out/'small/strict-summary.json',{k:v[k]for k in ('status','native','request','selected_absolute_step','physical_mass_dual_norms','callback_history_sha256')})
    target=root/'science/target';outcome=load(target/'target-outcome.json')if(target/'target-outcome.json').exists()else dict(status='NOT_RUN')
    for f in ('target-outcome.json','outcome.json','case.evolution.json','execution-diff.json','view/request.json'):
        copy(target/f,'target/'+f.replace('/','--'))
    if(target/'view/result.json').exists():
        v=load(target/'view/result.json');atomic(out/'target/strict-summary.json',{k:v[k]for k in ('status','native','request','selected_absolute_step','physical_mass_dual_norms','callback_history_sha256')})
        # Full signed weak/interface arrays stay verifiably on Linux.
        atomic(out/'target/strict-heavy-binding.json',dict(result=identity(target/'view/result.json'),native=identity(target/'view/native.json')))
    for p in sorted((root/'science/analysis-v1').glob('*')):
        if p.suffix in ('.csv','.png')or p.name in ('summary.json','native-and-independent.json','witness-history.json'):copy(p,'analysis/'+p.name)
        elif re.fullmatch(r'step-\d+.json',p.name):
            v=load(p);atomic(out/'analysis'/p.name,dict(step=v['step'],time=v['time'],quality=v['final_union_status'],qualified_restart_eligible=v['qualified_restart_eligible'],E_initial=v['E_initial'],E_eq=v['E_eq'],E_drift=v['E_drift'],all_regional_Q1_moments=v['all_regional_Q1_moments']))
    commits=[];interfaces=[];growth=[]
    for p in sorted((target/'monitor').glob('step-*.ready.result.json')):
        v=load(p);step=int(v['binding']['absolute_step']);s=v['summary'];inv=load(target/f'monitor/step-{step}.ready.inventory.json')
        for item in s['owner_interfaces']:interfaces.append(dict(step=step,time=step*.0025,**item))
        growth.append(dict(step=step,physical_points=len(inv['points']),owned_rows=len(inv['owned']),physical_rows=inv['physical_rows'],quadrature_rows=inv['quadrature_rows'],cache=v['cache'],assessment_seconds=v['measured_independent_seconds']))
        for suffix in ('ready.json','qualification-commit.json','ack.consumed.json','report.retention-fields.json','ready.remaining.json'):
            copy(target/f'monitor/step-{step}.{suffix}',f'target/step-{step}.{suffix}')
        c=target/f'monitor/step-{step}.qualification-commit.json'
        if c.exists():commits.append(identity(c))
    first=target/'monitor/step-1.ready.json'
    if first.exists():
        ready=load(first);assert ready['binding']['callback_history_sha256']=='40d7d7178d0ec1cada609e6267ce46b582ba869daf4725e5ed501759e87bed1a'
        assert ready['binding']['state_sha256']=='b9c1c19b6c44513d1778143b695e9720d0463d0a9786a91ec6ae97cc6c469819'
        atomic(out/'target/H14-replay-callback-binding.json',dict(status='PASS',ready=identity(first),genuine_callback=ready['callback'],state_sha256=ready['binding']['state_sha256'],callback_history_sha256=ready['binding']['callback_history_sha256'],scope='direct retained H14 callback/state comparison; native-history/reference/defect/time replay also recorded in independent step1 result'))
    atomic(out/'target/interface-history.json',interfaces);atomic(out/'target/inventory-growth.json',growth)
    gates=dict(repaired_build='PASS'if book['entries']['focused-build-v1']['status']=='PASS'else'FAIL',synthetic_protocol=book['entries']['repaired-native-protocol-fixture-v1']['status'],focused_python_methods=count,focused_python=tests['status'],retained_CP0_CP1_parity='PASS'if(root/'science/retained-parity/result.json').exists()else'NOT_RUN',small_actual_attempts=sum(v['actual_attempts']for v in small),small_governing_qualification='PASS'if len(small)==5 and all(v['status']=='PASS'for v in small)else'NOT_RUN_OR_FAILED',small_strict='PASS'if(root/'science/small/A/view/result.json').exists()else'NOT_RUN',target=outcome,qualified_commits=len(commits),fresh_CP0=0,source_builds=1,target_strict='PASS'if(target/'view/result.json').exists()else'NOT_RUN',strict_readbacks=len(book.get('strict_reservations',[])),historical_H15='OBSTRUCTED_PROTOCOL unchanged; final-v3 essential originals verified; full optional archive not transported',H11_temporal='UNRESOLVED unchanged',H12_audit='PASS unchanged')
    atomic(out/'component-gates.json',gates);atomic(out/'copied-evidence.json',copied)
    heavy=[]
    # Only H16 science is traversed; no old root or historical corpus.
    for p in sorted((root/'science').rglob('*')):
        if p.is_file():heavy.append(identity(p))
    atomic(out/'remote-heavy-manifest.json',dict(root=str(root),files=heavy))
    files=[identity(p)for p in sorted(out.rglob('*'))if p.is_file()]
    atomic(out/'manifest.json',dict(schema='WJH16PublicationV1',source_python_revision=revision,native_source_revision=load(root/'source-intent.json')['native_revision'],files=files,active_operation_logs_excluded=True,accounting_scope='publication OPEN snapshot; closure separately sealed'))
    atomic(out/'seal.json',dict(manifest=identity(out/'manifest.json'),files=len(files)))
    print(json.dumps(dict(status='PASS',publication=identity(out/'manifest.json'),gates=gates)),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--source',required=True);a=p.parse_args();publish(a.root,a.source)
