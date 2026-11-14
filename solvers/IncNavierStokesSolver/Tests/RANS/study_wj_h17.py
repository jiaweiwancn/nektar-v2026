#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Fixed H17 saved-audit and replay orchestration, using existing strict helpers."""
import argparse
import copy
import hashlib
import os
from pathlib import Path
import subprocess
import sys
import time
import uuid
import xml.etree.ElementTree as ET
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,verify,atomic,vector_hash,RULE,H14
from study_wj_h14 import runtime as prepare_runtime
from study_wj_h15 import strict_request,qualified_step
from study_wj_h16 import same_cp,checked_exit
from wj_h17_observer import compatibility,equal,native_load
from wj3_measure import checkpoint_payload
from run_wj_h17 import frozen,halt,reserve_strict
from cls_epawf_m3_cases import set_info,set_parameter,write
HERE=Path(__file__).resolve().parent
OLD=Path('/home/jay/build/wj-h16-repaired-monitor-20261114')
CAP_ENTRIES=196608
CAP_BYTES=256*1024**2
POLICIES=('Off','PatchExactV1')

EXEC_INFO={'RANSWallRestart','RANSWallCheckpoint','RANSWallPimpleSnapshotDirectory','RANSWallRetention','RANSWallMonitorIntent','RANSWallMonitorFault','RANSWallReport','RANSWallImmutableRowReuse','RANSWallImmutableRowEntries','RANSWallImmutableRowBytes'}
def protected(x):
    x=copy.deepcopy(x)
    for p in list(x.find('./CONDITIONS/PARAMETERS')):
        if p.text.split('=')[0].strip()in ('NumSteps','IO_CheckSteps','IO_InfoSteps'):x.find('./CONDITIONS/PARAMETERS').remove(p)
    for p in list(x.find('./CONDITIONS/SOLVERINFO')):
        if p.get('PROPERTY')in EXEC_INFO:x.find('./CONDITIONS/SOLVERINFO').remove(p)
    for p in x.iter():
        if p.text is not None and not p.text.strip():p.text=None
        p.tail=None
    return ET.tostring(x)
def xml_case(original,folder,cp,steps,policy='PatchExactV1'):
    assert policy in POLICIES and steps in (0,1,2,4)
    folder.mkdir(parents=True);x=ET.parse(original).getroot()
    for k,v in (('NumSteps',steps),('IO_CheckSteps',1),('IO_InfoSteps',1)):set_parameter(x,k,v)
    for k,v in (('RANSWallRestart',str(cp)),('RANSWallCheckpoint',str(folder/'CP')),('RANSWallPimpleSnapshotDirectory',str(folder/'snapshots')),
        ('RANSWallRetention','BoundedV1'),('RANSWallMonitorIntent',str(folder/'intent.json')),('RANSWallMonitorFault','Off'),('RANSWallReport',str(folder/'native.json')),
        ('RANSWallImmutableRowReuse',policy),('RANSWallImmutableRowEntries',str(CAP_ENTRIES)),('RANSWallImmutableRowBytes',str(CAP_BYTES))):set_info(x,k,v)
    path=folder/'case.xml';write(x,path)
    assert protected(ET.parse(original).getroot())==protected(x),'scientific XML difference'
    assert next(p.get('VALUE')for p in x.findall('./CONDITIONS/SOLVERINFO/I')if p.get('PROPERTY')=='RANSWallEvaluationReuse')=='ExactStepV1'
    atomic(folder/'execution-diff.json',dict(original=identity(original),as_run=identity(path),protected_sha256=hashlib.sha256(protected(x)).hexdigest(),additional_steps=steps,policy=policy,entry_cap=CAP_ENTRIES,byte_cap=CAP_BYTES,changed_only=sorted(EXEC_INFO|{'NumSteps','IO_CheckSteps','IO_InfoSteps'})))
    return path

def setup(root):
    prepare_runtime(root);b=load(root/'science/runtime-build.json');b.update(package='H17',native_source_revision=load(root/'source-intent.json')['native_revision'],scope='exact native structural-row lifetime only; historical producers remain H14/H16')
    atomic(root/'science/runtime-build-v1.json',b)
    old=load(OLD/'science/runtime-qualified-v1.json');libs=[]
    for rec in old['libraries']:
        verify(rec);p=Path(rec['path']);new=root/'runtime'/p.relative_to(OLD/'runtime')if p.is_relative_to(OLD/'runtime')else p
        v=identity(new);assert v['sha256']==rec['sha256']and v['bytes']==rec['bytes'];libs.append(dict(actual=v,historical=rec))
    atomic(root/'science/runtime-dependencies-v1.json',dict(status='PASS',libraries=libs,historical_qualification=identity(OLD/'science/runtime-qualified-v1.json'),actual_native_maps_required=True))
    compatibility(root,Path('/home/jay/build/wj-h15-persistent-neutral-20261114/science/observer-cache/manifest.json'),root/'science/target-cache-compatibility.json')
    compatibility(root,Path('/home/jay/build/wj-h15-persistent-neutral-20261114/science/small-cache-v2/manifest.json'),root/'science/small-cache-compatibility.json')
    from publish_wj_h15 import objects
    atomic(root/'science/compiled-source-objects.json',objects(root,root/'build',b['native_source_revision']))
    inputs(root)

def expected(root,scope):
    out=root/'inputs'/scope;out.mkdir();records={}
    count=4 if scope=='target'else 2
    for n in range(1,count+1):
        directory=OLD/('science/target'if scope=='target'else'science/small/A')/'monitor';ready=directory/f'step-{n}.ready.json';r=load(ready)
        v={k:r[k]for k in ('checkpoint','callback','report','structure')}
        v.update(ready=identity(ready),commit=identity(directory/f'step-{n}.qualification-commit.json'),result=identity(directory/f'step-{n}.ready.result.json'),inventory=identity(directory/f'step-{n}.ready.inventory.json'))
        for p in v.values():verify(p)
        atomic(out/f'step-{n}.json',v);records[str(n)]=identity(out/f'step-{n}.json')
    atomic(out/'manifest.json',dict(status='PASS',states=records));return records

def cp_request(v):
    return dict(native_report=str(verify(v['producer'])),native_report_sha256=v['producer']['sha256'],native_checkpoint=str(verify(v['checkpoint']).parent),checkpoint_sha256=v['checkpoint']['sha256'])
def completed_request(v):
    for k in ('report','checkpoint','callback','commit'):verify(v[k])
    return dict(native_report=v['report']['path'],native_report_sha256=v['report']['sha256'],native_checkpoint=str(Path(v['checkpoint']['path']).parent),checkpoint_sha256=v['checkpoint']['sha256'],completed_sidecar=v['callback'],qualified_cursor=v['commit'])

def inputs(root):
    target=expected(root,'target');small=expected(root,'small')
    v=next(v for v in load(H14/'publication-v2/state-bindings.json')['states']if v['case']=='cp0')
    zero=cp_request(v);requests=[zero]+[completed_request(load(verify(target[str(n)])))for n in range(1,5)]+[zero]
    intents=[]
    for policy in POLICIES:
        folder=root/'science/audits'/policy;xml=xml_case(verify(v['xml']),folder,verify(v['checkpoint']).parent,0,policy)
        atomic(folder/'request.json',dict(operation='complete_audit_batch',inner_stress_model='Linear',maximum_seconds=900.,states=requests))
        cmd=[str(root/'runtime/RANSWallResolutionProbe'),str(xml),str(folder/'request.json'),str(folder/'native.json')]
        pins=[identity(xml),identity(folder/'request.json'),v['producer'],v['checkpoint']]
        for p in target.values():pins.append(p);pins.extend(load(verify(p)).values())
        intents.append(dict(policy=policy,command=cmd,inputs=pins))
    atomic(root/'audit-admission.json',dict(schema='WJH17FixedAuditAdmissionV1',policies=intents,full_evaluations=12,sequence=[0,1,2,3,4,0],entry_cap=CAP_ENTRIES,byte_cap=CAP_BYTES))
    small_v=load(verify(small['2']));base=OLD/'science/small/A/case.xml';folder=root/'science/row-fixture';xml=xml_case(base,folder,Path(small_v['checkpoint']['path']).parent,0)
    atomic(folder/'request.json',dict(operation='immutable_row_checks',inner_stress_model='Linear',**completed_request(small_v)))
    # Prepared inputs do not authorize governing work; parity remains mandatory.

def freeze(root):
    files=[identity(HERE/n)for n in ('run_wj_h17.py','study_wj_h17.py','wj_h17_observer.py','test_wj_h17_exact_rows.py')]
    for name in ('runtime-build-v1.json','runtime-dependencies-v1.json','compiled-source-objects.json','target-cache-compatibility.json','small-cache-compatibility.json'):files.append(identity(root/'science'/name))
    for v in load(root/'science/runtime-build-v1.json')['executables']:files.extend((v['ELF'],v['wrapper']))
    for v in load(root/'science/target-cache-compatibility.json')['active_sources']:files.append(v['active'])
    files.append(identity(root/'audit-admission.json'))
    atomic(root/'execution-config.json',dict(schema='WJH17FrozenExecutionV1',status='FROZEN_BEFORE_QUALIFICATION',files=files,full_audit_calls_before_freeze=0,small_governing_before_freeze=0,policy_caps=dict(entries=CAP_ENTRIES,bytes=CAP_BYTES)))

def compare(root):
    out=root/'science/audit-comparison.json';old=load(root/'science/audits/Off/native.json');new=load(root/'science/audits/PatchExactV1/native.json')
    assert old['full_audit_calls']==new['full_audit_calls']==6
    keys=('step','time','alpha','initial_history_role','state_sha256','history_sha256','producer_sha256','checkpoint_file_sha256','integrator_primary_sha256','native_history_sha256','native_history_times','reference_sha256','defect_sha256','qh0_sha256','constraint_rows_sha256','constraint_values_sha256','spatial_residual_sha256','completed_residual_sha256','exterior_load_sha256','boundary_exact_native_json_sha256','roughness_budget_sha256','energy_budget_sha256')
    records=[]
    for a,b in zip(old['results'],new['results']):
        for k in keys:equal(a[k],b[k],k)
        for k in ('audit','audit_numerical_sha256','ordered_rows','ordered_observations','canonical','localization'):equal(a['complete'][k],b['complete'][k],k)
        records.append(dict(step=a['step'],status='PASS',audit_sha256=a['complete']['audit_numerical_sha256'],row_stream=a['complete']['ordered_rows']['root_sha256'],row_visits=a['complete']['ordered_rows']['count'],observation_stream=a['complete']['ordered_observations']['root_sha256'],observations=a['complete']['ordered_observations']['count'],baseline_seconds=a['audit_seconds'],new_seconds=b['audit_seconds'],baseline_commitment_seconds=a['complete']['commitment_wall_seconds'],new_commitment_seconds=b['complete']['commitment_wall_seconds']))
    # The repeated state0 must recompute current results with dynamic scopes reset.
    for v in (old,new):equal(v['results'][0]['complete']['audit'],v['results'][5]['complete']['audit'],'A->B->A audit')
    assert new['immutable_rows']['live_bytes']<=CAP_BYTES
    atomic(out,dict(status='PASS',results=records,baseline=identity(root/'science/audits/Off/native.json'),new=identity(root/'science/audits/PatchExactV1/native.json'),scope='all numerical maps, ordered IEEE row/current-value streams and affected outputs exact; metadata/timings/counters excluded'))

def runtime_qualify(root,governor_only=False):
    """Bind actual native-exec maps to the new ELFs and pinned dependency bytes."""
    build=load(root/'science/runtime-build-v1.json');expected={v['ELF']['path']:v['ELF']for v in build['executables']}
    dependencies={v['actual']['path']:v['actual']for v in load(root/'science/runtime-dependencies-v1.json')['libraries']}
    samples=[];seen=set();libraries={}
    for p in sorted((root/'operations').glob('*/native-exec-samples.json')):
        for v in load(p):
            assert v['status']=='PASS'
            for native in v['native']:verify(expected[native]);seen.add(Path(native).name)
            for f in v['files']:
                if '.so'not in Path(f).name and not Path(f).name.startswith('ld-linux'):continue
                assert f in dependencies,'unqualified actual native dependency: '+f
                libraries[f]=identity(verify(dependencies[f]))
            samples.append(dict(operation=p.parent.name,pid=v['pid'],native=v['native'],libc=v['libc']))
    required={'IncNavierStokesSolver.elf'}if governor_only else{'IncNavierStokesSolver.elf','RANSWallResolutionProbe.elf'}
    assert required<=seen,'missing actual governor/reader native-exec sample'
    atomic(root/'science'/('runtime-governor-qualified-v1.json'if governor_only else'runtime-qualified-v1.json'),dict(status='PASS',native_exec_samples=samples,libraries=list(libraries.values()),compiled_build=identity(root/'science/runtime-build-v1.json'),dependencies=identity(root/'science/runtime-dependencies-v1.json'),scope='actual H17 execution; native identity separate from historical cache'))

def coordinate(root,folder,original_cp0,compat,*,behavior='PASS',restart_callback=None,replay=None,deadline=None):
    frozen(root);xml=folder/'case.xml';tree=ET.parse(xml);steps=int(next(p.text.split('=')[1]for p in tree.findall('./CONDITIONS/PARAMETERS/P')if p.text.split('=')[0].strip()=='NumSteps'))
    restart=next(p.get('VALUE')for p in tree.findall('./CONDITIONS/SOLVERINFO/I')if p.get('PROPERTY')=='RANSWallRestart');restored=checkpoint_payload(restart,include_reference=True);qh0=checkpoint_payload(original_cp0,include_reference=True)
    directory=folder/'monitor';directory.mkdir();ip=folder/'intent.json'
    intent=dict(schema='WJMonitorIntentV1',directory=str(directory),nonce=uuid.uuid4().hex,scope_sha256=hashlib.sha256((identity(xml)['sha256']+identity(root/'execution-config.json')['sha256']).encode()).hexdigest(),
        source_sha256=identity(root/'science/runtime-build-v1.json')['sha256'],observer_sha256=identity(HERE/'wj_h17_observer.py')['sha256'],inventory_rule_sha256=hashlib.sha256(RULE.encode()).hexdigest(),
        original_qh0_checkpoint=identity(Path(original_cp0)/'State.bin'),original_qh0_state_sha256=vector_hash(qh0['state']),maximum_absolute_step=restored['steps']+steps,ack_timeout_seconds=90.,strict_reserve_seconds=120)
    if restart_callback:intent['restart_callback']=identity(restart_callback)
    if replay:intent['replay_expected']=replay
    if deadline:intent['deadline_monotonic']=deadline
    native=None;obs=None;error=None;code=None;began=time.monotonic()
    with(folder/'observer.log').open('x')as obslog,(folder/'native.log').open('x')as log:
        try:
            obs=subprocess.Popen([sys.executable,'-B',str(HERE/'wj_h17_observer.py'),'--root',str(root),'--mode','daemon','--intent',str(ip),'--compat',str(compat),'--behavior',behavior],stdout=obslog,stderr=subprocess.STDOUT)
            intent['observer_pid']=obs.pid;atomic(ip,intent)
            while not(directory/'observer-started.json').exists():
                assert obs.poll()is None,'observer failed before native launch'
                assert time.monotonic()-began<30,'observer startup timeout';time.sleep(.025)
            observer_setup=time.monotonic()-began
            native=subprocess.Popen([str(root/'runtime/IncNavierStokesSolver'),str(xml)],stdout=log,stderr=subprocess.STDOUT)
            atomic(folder/'native-launch.json',dict(pid=native.pid,xml=identity(xml),intent=identity(ip),requested_additional_steps=steps,observer_setup_seconds=observer_setup,restored_absolute_step=restored['steps']))
            while native.poll()is None:
                if obs.poll()is not None:raise RuntimeError('observer exited during native run; no recovery/next advance')
                time.sleep(.025)
            code=native.returncode
        except Exception as e:error=repr(e)
        finally:
            if native and native.poll()is None:
                native.terminate()
                try:native.wait(timeout=5)
                except subprocess.TimeoutExpired:native.kill();native.wait()
            atomic(directory/'observer-stop.json',dict(reason='H17 coordinator closure'))
            if obs:
                try:obs.wait(timeout=5)
                except subprocess.TimeoutExpired:obs.terminate();obs.wait(timeout=5)
    report=folder/'case.evolution.json';r=load(report)if report.exists()else None
    atomic(folder/'outcome.json',dict(native_exit=code,error=error,report=identity(report)if r else None,wall_seconds=time.monotonic()-began,restored_state_sha256=vector_hash(restored['state'])))
    if r:
        startup=load(directory/'startup.json')
        for k in ('state','fixed_reference','fixed_defect'):
            assert vector_hash(startup[k])==vector_hash(restored[k]),'restore changed '+k
        same_cp(startup,restored)
        assert startup['original_qh0_state_sha256']==vector_hash(qh0['state'])
        assert r['initialization_work']==[0,0]and r['initializer_calls']==0
        op=os.environ.get('WJ_H17_OPERATION')
        if op:atomic(Path(op)/'attempt-counts.json',dict(kind='target'if folder.name=='target'else'small',actual_attempts=r['implicit_callbacks'],producer=identity(report)))
    if error:raise RuntimeError(error)
    assert r is not None,'native summary unavailable; attempts remain charged unknown'
    return code,r

def strict(root,scope):
    folder=root/'science'/('small/A'if scope=='small'else'target');reserve_strict(root,scope)
    if qualified_step(folder):
        cmd=strict_request(root,folder);out=folder/'view'
    else:
        assert scope=='target';v=next(v for v in load(H14/'publication-v2/state-bindings.json')['states']if v['case']=='cp0');out=folder/'view';out.mkdir()
        request=dict(operation='open_checks',inner_stress_model='Linear',native_report=str(verify(v['producer'])),native_report_sha256=v['producer']['sha256'],native_checkpoint=str(verify(v['checkpoint']).parent),checkpoint_sha256=v['checkpoint']['sha256'],weak=True,weak_compact=True,assembly=False,reuse=False,retain_scalar_observations=False)
        atomic(out/'request.json',request);cmd=[str(root/'runtime/RANSWallResolutionProbe'),str(verify(v['xml'])),str(out/'request.json'),str(out/'native.json')]
    with(out/'native.log').open('x')as log:code=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT).returncode
    checked_exit(code);p=load(out/'native.json');req=load(out/'request.json')
    assert p['initialization_work']==[0,0]and p['initial_velocity_calls']==p['governing_solves']==0 and p['scalar_domain']['admissible']
    assert p['report_sha256']==req['native_report_sha256']and p['checkpoint_sha256']==req['checkpoint_sha256']
    if scope=='small':assert p['selected_absolute_step']==2
    atomic(out/'result.json',dict(status='PASS',native=identity(out/'native.json'),request=identity(out/'request.json'),selected_absolute_step=p['selected_absolute_step'],physical_mass_dual_norms=p['physical_mass_dual_norms'],weak=p['weak_checks'],interface=p['open_boundary'],callback_history_sha256=p['callback_history_sha256']))

def small_case(root,name):
    assert name in ('A','B');frozen(root);assert load(root/'science/audit-comparison.json')['status']=='PASS'
    parent=H14/'science/small/cp0';cp0=Path(load(parent/'case.evolution.json')['checkpoint']);folder=root/'science/small'/name;restart=cp0;side=None
    if name=='B':
        assert load(root/'science/small/A/qualification.json')['status']=='PASS'
        c=load(root/'science/small/A/monitor/step-1.qualification-commit.json');restart=verify(c['checkpoint']).parent;side=verify(c['callback'])
    xml_case(parent/'case.xml',folder,restart,2 if name=='A'else 1)
    expected_states=load(root/'inputs/small/manifest.json')['states']
    code,r=coordinate(root,folder,cp0,root/'science/small-cache-compatibility.json',restart_callback=side,replay=expected_states)
    checked_exit(code);assert r['implicit_callbacks']==(2 if name=='A'else 1)and r['native_completed']==r['independently_qualified']==r['implicit_callbacks']
    assert qualified_step(folder)==2
    same_cp(checkpoint_payload(r['last_qualified_checkpoint'],include_reference=True),checkpoint_payload(verify(load(verify(expected_states['2']))['checkpoint']).parent,include_reference=True))
    atomic(folder/'qualification.json',dict(status='PASS',actual_attempts=r['implicit_callbacks'],completed=r['native_completed'],qualified=r['independently_qualified'],producer=identity(folder/'case.evolution.json'),new_scientific_progress=0))

def admission(root):
    frozen(root)
    assert load(root/'science/audit-comparison.json')['status']=='PASS'
    for n in ('A','B'):assert load(root/f'science/small/{n}/qualification.json')['status']=='PASS'
    assert load(root/'science/small/A/view/result.json')['status']=='PASS'
    assert load(root/'science/runtime-qualified-v1.json')['status']=='PASS'
    book=load(root/'accounting.json');assert not any(book.get(k)for k in ('runtime_halt','qualification_halt','coverage_halt','target_halt'))
    peak=max(load(root/f'operations/{n}/receipt.json')['peak_tree_rss_bytes']for n in ('complete-Off-v1','complete-PatchExactV1-v1'))
    # Conservative native H16 peak + entire independent cache and new row cap.
    memory=6983356416+CAP_BYTES;assert memory<8*1024**3
    forecast=1.25*(4*(95+15)+150+120+60);assert forecast<1800
    charged=sum(v.get('seconds',v['reserved_seconds'])for v in book['entries'].values());assert 60+charged+forecast+300<7200
    v=next(v for v in load(H14/'publication-v2/state-bindings.json')['states']if v['case']=='cp0');folder=root/'science/target';xml=xml_case(verify(v['xml']),folder,verify(v['checkpoint']).parent,4)
    cmd=[sys.executable,'-B',str(HERE/'study_wj_h17.py'),'--root',str(root),'--mode','target']
    atomic(root/'target-admission.json',dict(status='PASS',command=cmd,inputs=[identity(root/'execution-config.json'),identity(xml),v['producer'],v['checkpoint'],v['xml'],identity(root/'science/runtime-qualified-v1.json'),identity(root/'inputs/target/manifest.json')],height_mm=6,state_size=24688,maximum_attempts=4,new_scientific_progress=0,forecast_with25percent_headroom=forecast,memory_upper_estimate=memory,passive_peak_bytes=peak,strict_reserved=120,closure_reserved=300))

def target(root):
    folder=root/'science/target';began=time.monotonic();v=next(v for v in load(H14/'publication-v2/state-bindings.json')['states']if v['case']=='cp0');error=None;code=None;r={}
    try:code,r=coordinate(root,folder,verify(v['checkpoint']).parent,root/'science/target-cache-compatibility.json',replay=load(root/'inputs/target/manifest.json')['states'],deadline=began+1680)
    except Exception as e:error=repr(e)
    if(folder/'case.evolution.json').exists():r=load(folder/'case.evolution.json')
    n=qualified_step(folder);failure=error or r.get('failure')or('incomplete replay'if n!=4 else None)
    if code!=0 or failure:halt(root,'target_halt',dict(first_reason=failure or f'native exit {code}',producer=identity(folder/'case.evolution.json')if r else None,qualified_cursor_absolute_step=n))
    atomic(folder/'target-outcome.json',dict(status='PASS'if code==0 and n==4 and r.get('native_completed')==r.get('independently_qualified')==4 else'OBSTRUCTED',native_exit=code,native_completed=r.get('native_completed'),qualified=n,actual_attempts=r.get('implicit_callbacks'),new_scientific_progress=0,first_failure=failure))
    if error:raise RuntimeError(error)
    checked_exit(code);assert n==4 and r['time']==.01

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('setup','freeze','compare','small','small-strict','target-strict','runtime','admission','target'),required=True);p.add_argument('--case',choices=('A','B'));a=p.parse_args()
    if a.mode=='small':small_case(a.root,a.case)
    elif a.mode.endswith('-strict'):strict(a.root,a.mode.split('-')[0])
    else:globals()[a.mode](a.root)
