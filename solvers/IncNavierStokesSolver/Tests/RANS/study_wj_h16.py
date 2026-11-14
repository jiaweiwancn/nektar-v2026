#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Fixed H16 qualification cases and one source-bound four-stage target."""
import argparse
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
from study_wj_h15 import xml_case,strict_request,qualified_step
from wj_h16_observer import compatibility
from wj3_measure import checkpoint_payload
from run_wj_h16 import frozen,halt,reserve_strict

HERE=Path(__file__).resolve().parent
OLD=Path('/home/jay/build/wj-h15-persistent-neutral-20261114')
CASES=dict(A=(2,'PASS','Off',False),B=(1,'PASS','Off',False),C=(2,'FAIL','Off',False),D=(1,'PASS','Off',True),E=(1,'PASS','AfterImplicitV1',False))

def same_cp(a,b):
    for k in ('state','fixed_reference','fixed_defect'):assert vector_hash(a[k])==vector_hash(b[k]),k
    assert len(a['native_history'])==len(b['native_history'])
    assert all(vector_hash(x)==vector_hash(y)for x,y in zip(a['native_history'],b['native_history']))
    assert a['native_history_times']==b['native_history_times']and a['time']==b['time']and a['steps']==b['steps']and a['dt']==b['dt']and a['signature']==b['signature']

def checked_exit(actual,expected=0):
    if actual!=expected:raise RuntimeError(f'coordinator child exit {actual}, expected {expected}')

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
    atomic(root/'science'/('runtime-governor-qualified-v1.json'if governor_only else'runtime-qualified-v1.json'),dict(status='PASS',native_exec_samples=samples,libraries=list(libraries.values()),compiled_build=identity(root/'science/runtime-build-v1.json'),dependencies=identity(root/'science/runtime-dependencies-v1.json'),scope='actual H16 execution; native identity separate from historical cache'))

def setup(root):
    prepare_runtime(root);b=load(root/'science/runtime-build.json');b.update(package='H16',native_source_revision=load(root/'source-intent.json')['native_revision'],scope='H16 parent-durable repaired governor/reader; original private runtime dependencies')
    atomic(root/'science/runtime-build-v1.json',b)
    prior=load(verify(load(root/'inputs/essential-index.json')['runtime_qualification']));libs=[]
    for v in prior['libraries']:
        old=v['actual'];path=Path(old['path']);new=root/'runtime'/path.relative_to(OLD/'runtime')if path.is_relative_to(OLD/'runtime')else path
        rec=identity(new);assert rec['bytes']==old['bytes']and rec['sha256']==old['sha256'];libs.append(dict(actual=rec,historical=old))
    atomic(root/'science/runtime-dependencies-v1.json',dict(status='PASS',libraries=libs,historical_qualification=identity(OLD/'science/runtime-qualified-v1.json'),new_native=identity(root/'science/runtime-build-v1.json'),actual_native_maps_required=True))
    compatibility(root,OLD/'science/observer-cache/manifest.json',root/'science/target-cache-compatibility.json')
    compatibility(root,OLD/'science/small-cache-v2/manifest.json',root/'science/small-cache-compatibility.json')

def freeze(root):
    verify(load(root/'inputs/essential-index.json')['legacy_producer'])
    assert load(root/'science/retained-parity/result.json')['status']=='PASS'
    files=[identity(HERE/n)for n in ('run_wj_h16.py','study_wj_h16.py','wj_h16_observer.py')]
    for name in ('runtime-build-v1.json','runtime-dependencies-v1.json','compiled-source-objects.json','target-cache-compatibility.json','small-cache-compatibility.json'):files.append(identity(root/'science'/name))
    for v in load(root/'science/runtime-build-v1.json')['executables']:files.extend((v['ELF'],v['wrapper']))
    for v in load(root/'science/target-cache-compatibility.json')['active_sources']:files.append(v['active'])
    atomic(root/'execution-config.json',dict(schema='WJH16FrozenExecutionV1',status='FROZEN_BEFORE_GOVERNING_QUALIFICATION',files=files,native_work_before_freeze=0))

def coordinate(root,folder,original_cp0,compat,*,behavior='PASS',restart_callback=None,replay=None,deadline=None):
    frozen(root);xml=folder/'case.xml';tree=ET.parse(xml);steps=int(next(p.text.split('=')[1]for p in tree.findall('./CONDITIONS/PARAMETERS/P')if p.text.split('=')[0].strip()=='NumSteps'))
    restart=next(p.get('VALUE')for p in tree.findall('./CONDITIONS/SOLVERINFO/I')if p.get('PROPERTY')=='RANSWallRestart');restored=checkpoint_payload(restart,include_reference=True);qh0=checkpoint_payload(original_cp0,include_reference=True)
    directory=folder/'monitor';directory.mkdir();ip=folder/'intent.json'
    intent=dict(schema='WJMonitorIntentV1',directory=str(directory),nonce=uuid.uuid4().hex,scope_sha256=hashlib.sha256((identity(xml)['sha256']+identity(root/'execution-config.json')['sha256']).encode()).hexdigest(),
        source_sha256=identity(root/'science/runtime-build-v1.json')['sha256'],observer_sha256=identity(HERE/'wj_h16_observer.py')['sha256'],inventory_rule_sha256=hashlib.sha256(RULE.encode()).hexdigest(),
        original_qh0_checkpoint=identity(Path(original_cp0)/'State.bin'),original_qh0_state_sha256=vector_hash(qh0['state']),maximum_absolute_step=restored['steps']+steps,ack_timeout_seconds=90.,strict_reserve_seconds=120)
    if restart_callback:intent['restart_callback']=identity(restart_callback)
    if replay:intent['replay_H14_step1']=replay
    if deadline:intent['deadline_monotonic']=deadline
    native=None;obs=None;error=None;code=None;began=time.monotonic()
    with(folder/'observer.log').open('x')as obslog,(folder/'native.log').open('x')as log:
        try:
            obs=subprocess.Popen([sys.executable,'-B',str(HERE/'wj_h16_observer.py'),'--root',str(root),'--mode','daemon','--intent',str(ip),'--compat',str(compat),'--behavior',behavior],stdout=obslog,stderr=subprocess.STDOUT)
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
            atomic(directory/'observer-stop.json',dict(reason='H16 coordinator closure'))
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
        op=os.environ.get('WJ_H16_OPERATION')
        if op:atomic(Path(op)/'attempt-counts.json',dict(kind='target'if folder.name=='target'else'small',actual_attempts=r['implicit_callbacks'],producer=identity(report)))
    if error:raise RuntimeError(error)
    assert r is not None,'native summary unavailable; attempts remain charged unknown'
    return code,r

def legacy_checkpoint(step):
    report=load(OLD/'science/small/legacy/case.evolution.json')
    if step==report['steps']:return checkpoint_payload(report['checkpoint'],include_reference=True)
    candidates=[]
    for p in (OLD/'science/small/legacy').glob('CP*'):
        if p.is_dir()and(p/'State.bin').exists():
            cp=checkpoint_payload(p,include_reference=True)
            if cp['steps']==step:candidates.append(cp)
    assert len(candidates)==1;return candidates[0]

def small_case(root,name):
    assert name in CASES;frozen(root);base=H14/'science/small/cp0';original=load(base/'case.evolution.json');cp0=Path(original['checkpoint']);steps,behavior,fault,hard=CASES[name]
    folder=root/'science/small'/name;restart=cp0;callback=None
    if name=='B':
        assert load(root/'science/small/A/qualification.json')['status']=='PASS'
        commit=load(root/'science/small/A/monitor/step-1.qualification-commit.json');restart=verify(commit['checkpoint']).parent;callback=verify(commit['callback'])
    xml_case(base/'case.xml',folder,restart,steps,fault=fault,hard=hard)
    code,r=coordinate(root,folder,cp0,root/'science/small-cache-compatibility.json',behavior=behavior,restart_callback=callback)
    checked_exit(code,0 if name in ('A','B')else 1)
    if name in ('A','B'):
        assert r['native_completed']==(2 if name=='A'else 1)and r['independently_qualified']==r['native_completed']
        assert qualified_step(folder)==2
        same_cp(checkpoint_payload(r['last_qualified_checkpoint'],include_reference=True),legacy_checkpoint(2))
        if name=='A':
            for n in (1,2):
                c=load(folder/f'monitor/step-{n}.qualification-commit.json');same_cp(checkpoint_payload(verify(c['checkpoint']).parent,include_reference=True),legacy_checkpoint(n))
                assert c['binding']['callback_history_sha256']==load(OLD/'science/small/legacy/case.evolution.json')['wj_audit_packets'][n]['stage_history_sha256']
    elif name=='C':
        assert r['native_completed']==1 and r['independently_qualified']==0 and r['implicit_callbacks']==1
        assert Path(r['last_qualified_checkpoint'])==cp0 and not(folder/'monitor/step-2.ready.json').exists()
        assert not(folder/'monitor/qualified-cursor.json').exists()and not(folder/'monitor/step-1.qualification-commit.json').exists()
        assert Path(r['last_completed_checkpoint'])!=cp0
    else:
        assert r['native_completed']==r['steps']==0 and r['aborted_advances']==1
        original_cp=checkpoint_payload(cp0,include_reference=True)
        assert r['state_sha256']==vector_hash(original_cp['state'])and r['native_history_sha256']==[vector_hash(v)for v in original_cp['native_history']]and r['native_history_times']==original_cp['native_history_times']
        a=load(folder/'monitor/first-aborted.json');verify(a['candidate_binary']);verify(a['callback_history_binary'])
        assert a['status']=='ABORTED_NOT_COMPLETED'and a['attempted_step']==1 and a['provisional_explicit_audits']and a['actual_explicit_calls']>=1
        if name=='D':assert a['first_native_range_event']['native_rejected_audit']['phase']=='rejected_stage'
        else:
            assert a['provisional_stage_audits']and a['provisional_audits'][-1]['phase']=='implicit_provisional'
            assert a['first_native_range_event'] is None
        assert not(folder/'monitor/step-1.ready.json').exists()
    assert r['resident_accepted_states']<=2 and r['resident_wj_packets']<=2 and r['resident_pimple_records']<=1 and r['resident_explicit_records']<=2 and r['factor_catalog']['resident_definitions']<=16
    atomic(folder/'qualification.json',dict(status='PASS',case=name,actual_attempts=r['implicit_callbacks'],native_completed=r['native_completed'],independently_committed=r['independently_qualified'],failure=r['failure'],producer=identity(folder/'case.evolution.json'),scope='expected fixture outcome, not a new physical target'))

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

def admission(root):
    frozen(root)
    for name in CASES:assert load(root/f'science/small/{name}/qualification.json')['status']=='PASS'
    assert load(root/'science/small/A/view/result.json')['status']=='PASS'and load(root/'science/retained-parity/result.json')['status']=='PASS'
    assert load(root/'science/runtime-qualified-v1.json')['status']=='PASS'
    book=load(root/'accounting.json');assert book['small_charged_attempts']==6
    baseline=6266638336;memory=1.25*baseline;forecast=1.25*(165+4*90+4*15+120+60+30)
    assert memory<8*1024**3 and forecast<1800
    v=next(v for v in load(H14/'publication-v2/state-bindings.json')['states']if v['case']=='cp0');cp=checkpoint_payload(verify(v['checkpoint']).parent,include_reference=True);assert vector_hash(cp['state'])=='631caa32a59ff43253f1fcac4effd955eb68b8529d9c28f29e615cf14ca0f077'
    folder=root/'science/target';xml_case(verify(v['xml']),folder,verify(v['checkpoint']).parent,4)
    cmd=[sys.executable,'-B',str(HERE/'study_wj_h16.py'),'--root',str(root),'--mode','target']
    atomic(root/'target-admission.json',dict(status='PASS',command=cmd,inputs=[identity(root/'execution-config.json'),identity(folder/'case.xml'),v['producer'],v['checkpoint'],v['xml'],identity(root/'science/runtime-qualified-v1.json')],height_mm=6,state_size=24688,forecast_with25percent_headroom=forecast,memory_with25percent_headroom=memory,strict_reserved=120,closure_reserved=300,maximum_attempts=4))

def target(root):
    folder=root/'science/target';began=time.monotonic();bindings=load(H14/'publication-v2/state-bindings.json')['states'];cp0=next(v for v in bindings if v['case']=='cp0');cp1=next(v for v in bindings if v['case']=='pilot1');error=None;code=None;r={}
    try:code,r=coordinate(root,folder,verify(cp0['checkpoint']).parent,root/'science/target-cache-compatibility.json',replay=cp1['producer'],deadline=began+1800)
    except Exception as e:error=repr(e)
    if(folder/'case.evolution.json').exists():r=load(folder/'case.evolution.json')
    qualified=qualified_step(folder)
    failure=error or r.get('failure')or('incomplete target'if qualified!=4 else None)
    if code!=0 or failure:halt(root,'target_halt',dict(first_reason=failure or f'native exit {code}',producer=identity(folder/'case.evolution.json')if r else None,qualified_cursor_absolute_step=qualified))
    atomic(folder/'target-outcome.json',dict(status='PASS'if code==0 and qualified==4 and r.get('native_completed')==r.get('independently_qualified')==4 else'OBSTRUCTED',native_exit=code,native_completed=r.get('native_completed'),independently_committed=r.get('independently_qualified'),published_qualified_step=qualified,new_scientific_progress=max(0,qualified-1),first_failure=failure))
    strict(root,'target')  # process/observer are both closed; protected after stop
    if error:raise RuntimeError(error)
    checked_exit(code);assert qualified==4 and r['native_completed']==r['independently_qualified']==4 and r['time']==.01
    print('H16 target four committed stages PASS; three new scientific steps',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('setup','freeze','small','small-strict','runtime-governor','runtime','admission','target'),required=True);p.add_argument('--case',choices=CASES);a=p.parse_args()
    if a.mode=='small':small_case(a.root,a.case)
    elif a.mode=='small-strict':strict(a.root,'small')
    elif a.mode in ('runtime','runtime-governor'):runtime_qualify(a.root,a.mode=='runtime-governor')
    else:globals()[a.mode](a.root)
