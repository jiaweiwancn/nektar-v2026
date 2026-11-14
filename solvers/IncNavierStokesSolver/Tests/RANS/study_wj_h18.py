#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Fixed H18 checkpoint cadence/restart qualification; no initialization route."""
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
from wj_h15_observer import identity,verify,atomic,vector_hash,bits,RULE,H14
from wj_h17_observer import compatibility,native_load,equal
from wj_h18_observer import parity,due
from study_wj_h17 import expected,cp_request,completed_request,runtime_qualify,OLD
from study_wj_h14 import runtime as prepare_runtime
from study_wj_h15 import strict_request,qualified_step
from study_wj_h16 import same_cp,checked_exit
from wj3_measure import checkpoint_payload
from cls_epawf_m3_cases import set_info,set_parameter,write
from run_wj_h18 import frozen,halt,reserve_strict
HERE=Path(__file__).resolve().parent
EXEC={'RANSWallRestart','RANSWallCheckpoint','RANSWallPimpleSnapshotDirectory','RANSWallRetention','RANSWallMonitorIntent','RANSWallMonitorFault','RANSWallReport','RANSWallImmutableRowReuse','RANSWallMonitoringPolicy','RANSWallPimpleSnapshots','RANSWallPimpleArithmetic'}


def protected(x):
    import copy
    x=copy.deepcopy(x)
    for p in list(x.find('./CONDITIONS/PARAMETERS')):
        if p.text.split('=')[0].strip()in ('NumSteps','IO_CheckSteps','IO_InfoSteps'):x.find('./CONDITIONS/PARAMETERS').remove(p)
    for p in list(x.find('./CONDITIONS/SOLVERINFO')):
        if p.get('PROPERTY')in EXEC:x.find('./CONDITIONS/SOLVERINFO').remove(p)
    for p in x.iter():
        if p.text is not None and not p.text.strip():p.text=None
        p.tail=None
    return ET.tostring(x)


def xml_case(original,folder,cp,steps,policy='CheckpointAuditV1',fault='Off'):
    assert policy in ('Off','CheckpointAuditV1')and steps in (1,2,4)
    folder.mkdir(parents=True);x=ET.parse(original).getroot()
    for k,v in (('NumSteps',steps),('IO_CheckSteps',1),('IO_InfoSteps',1)):set_parameter(x,k,v)
    for k,v in (('RANSWallRestart',str(cp)),('RANSWallCheckpoint',str(folder/'CP')),('RANSWallPimpleSnapshotDirectory',str(folder/'snapshots')),
        ('RANSWallRetention','BoundedV1'),('RANSWallMonitorIntent',str(folder/'intent.json')),('RANSWallMonitorFault',fault),('RANSWallReport',str(folder/'native.json')),
        ('RANSWallImmutableRowReuse','Off'),('RANSWallMonitoringPolicy',policy),('RANSWallPimpleSnapshots','Off'),('RANSWallPimpleArithmetic','Off')):set_info(x,k,v)
    assert protected(x)==protected(ET.parse(original).getroot()),'scientific input changed'
    path=folder/'case.xml';write(x,path)
    atomic(folder/'execution-diff.json',dict(original=identity(original),actual=identity(path),protected_sha256=hashlib.sha256(protected(x)).hexdigest(),execution_only=sorted(EXEC|{'NumSteps','IO_CheckSteps','IO_InfoSteps'}),additional_steps=steps,policy=policy,fault=fault))
    return path


def setup(root):
    prepare_runtime(root);b=load(root/'science/runtime-build.json');b.update(package='H18',native_source_revision=load(root/'source-intent.json')['native_revision'])
    atomic(root/'science/runtime-build-v1.json',b)
    old=load(OLD/'science/runtime-qualified-v1.json');libs=[]
    for rec in old['libraries']:
        verify(rec);p=Path(rec['path']);new=root/'runtime'/p.relative_to(OLD/'runtime')if p.is_relative_to(OLD/'runtime')else p
        v=identity(new);assert v['sha256']==rec['sha256'];libs.append(dict(actual=v,historical=rec))
    atomic(root/'science/runtime-dependencies-v1.json',dict(status='PASS',libraries=libs,historical_qualification=identity(OLD/'science/runtime-qualified-v1.json')))
    for scope,manifest in (('target','observer-cache'),('small','small-cache-v2')):
        compatibility(root,root.parent/f'wj-h15-persistent-neutral-20261114/science/{manifest}/manifest.json',root/f'science/{scope}-cache-compatibility.json')
    from publish_wj_h15 import objects
    atomic(root/'science/compiled-source-objects.json',objects(root,root/'build',b['native_source_revision']))
    for scope in ('target','small','small-restart'):expected(root,scope)
    # Target common support is H16 final union, not its evolving live inventory.
    for scope in ('target','small'):
        v=load(root/f'inputs/{scope}/step-{4 if scope=="target"else 2}.json')
        atomic(root/f'inputs/{scope}/common-support.json',load(verify(v['inventory'])))
    files=[identity(HERE/n)for n in ('run_wj_h18.py','study_wj_h18.py','wj_h18_observer.py','test_wj_h18_checkpoint_audits.py','wj_h18_checks.py')]
    for n in ('runtime-build-v1.json','runtime-dependencies-v1.json','compiled-source-objects.json','target-cache-compatibility.json','small-cache-compatibility.json'):files.append(identity(root/'science'/n))
    for v in b['executables']:files.extend((v['ELF'],v['wrapper']))
    for v in load(root/'science/target-cache-compatibility.json')['active_sources']:files.append(v['active'])
    atomic(root/'execution-config.json',dict(schema='WJH18FrozenExecutionV1',status='FROZEN_BEFORE_GOVERNING_QUALIFICATION',files=files,maximum_target_attempts=6,successful_small_attempts=5,deliberate_failed_small_attempts=1,target_policy='CheckpointAuditV1',immutable_rows='Off',evaluation='ExactStepV1'))


def coordinate(root,folder,cp0,scope,steps,*,policy='CheckpointAuditV1',restart_completion=None,fault=False):
    frozen(root);xml=folder/'case.xml';x=ET.parse(xml)
    restart=next(p.get('VALUE')for p in x.findall('./CONDITIONS/SOLVERINFO/I')if p.get('PROPERTY')=='RANSWallRestart')
    restored=checkpoint_payload(restart,include_reference=True);qh0=checkpoint_payload(cp0,include_reference=True)
    directory=folder/'monitor';directory.mkdir();ip=folder/'intent.json'
    intent=dict(schema='WJMonitorIntentV1',directory=str(directory),nonce=uuid.uuid4().hex,policy=policy,
        scope_sha256=hashlib.sha256((identity(xml)['sha256']+identity(root/'execution-config.json')['sha256']).encode()).hexdigest(),
        source_sha256=identity(root/'science/runtime-build-v1.json')['sha256'],observer_sha256=identity(HERE/'wj_h18_observer.py')['sha256'],inventory_rule_sha256=hashlib.sha256(RULE.encode()).hexdigest(),
        original_qh0_checkpoint=identity(Path(cp0)/'State.bin'),original_qh0_state_sha256=vector_hash(qh0['state']),maximum_absolute_step=restored['steps']+steps,
        restored_absolute_step=restored['steps'],recovery_stride=1,scientific_stride=4 if scope=='target'else 2,dt_bits=bits(restored['dt']),ack_timeout_seconds=90.,strict_reserve_seconds=120,substep_evidence_steps=[],preserve_qualification_artifacts=True,
        common_support=identity(root/f'inputs/{scope}/common-support.json'),expected_manifest=identity(root/f'inputs/{scope}/manifest.json'))
    if restart_completion:
        intent['restart_completion']=identity(restart_completion);intent['restart_callback']=load(restart_completion)['callback']
    obs=None;native=None;code=None;error=None;started=time.monotonic()
    with(folder/'observer.log').open('x')as obslog,(folder/'native.log').open('x')as log:
        try:
            obs=subprocess.Popen([sys.executable,'-B',str(HERE/'wj_h18_observer.py'),'--root',str(root),'--intent',str(ip),'--compat',str(root/f'science/{scope}-cache-compatibility.json')],stdout=obslog,stderr=subprocess.STDOUT)
            intent['observer_pid']=obs.pid;atomic(ip,intent)
            while not(directory/'observer-started.json').exists():
                assert obs.poll()is None,'observer setup failed'
                assert time.monotonic()-started<30;time.sleep(.025)
            native=subprocess.Popen([str(root/'runtime/IncNavierStokesSolver'),str(xml)],stdout=log,stderr=subprocess.STDOUT)
            atomic(folder/'native-launch.json',dict(pid=native.pid,xml=identity(xml),intent=identity(ip),requested_additional_steps=steps,restored_absolute_step=restored['steps']))
            while native.poll()is None:
                if obs.poll()is not None:raise RuntimeError('observer exited; no further advancement eligible')
                time.sleep(.025)
            code=native.returncode
        except Exception as e:error=repr(e)
        finally:
            if native and native.poll()is None:native.terminate();native.wait(timeout=5)
            atomic(directory/'observer-stop.json',dict(reason='H18 coordinator exit'))
            if obs:
                try:obs.wait(timeout=5)
                except subprocess.TimeoutExpired:obs.terminate();obs.wait(timeout=5)
    p=folder/'case.evolution.json';r=load(p)if p.exists()else None
    atomic(folder/'outcome.json',dict(native_exit=code,error=error,producer=identity(p)if r else None,wall_seconds=time.monotonic()-started))
    if r:
        op=os.environ.get('WJ_H18_OPERATION')
        if op:atomic(Path(op)/'attempt-counts.json',dict(kind=scope,actual_attempts=r['implicit_callbacks'],producer=identity(p)))
    if error:raise RuntimeError(error)
    assert r is not None,'missing native summary; reservation remains charged unknown'
    assert r['initialization_work']==[0,0]and r['initializer_calls']==0
    if fault:
        assert code!=0 and r['steps']==0 and r['aborted_advances']==1 and r['native_completed']==0
        assert r['state_sha256']==vector_hash(restored['state'])and r['native_history_sha256']==[vector_hash(h)for h in restored['native_history']]
        aborted=load(directory/'first-aborted.json');assert aborted['status']=='ABORTED_NOT_COMPLETED'and aborted['provisional_stage_audits']
        assert not(directory/'latest-native-completed.json').exists()and not(directory/'latest-durable-checkpoint.json').exists()
    else:
        checked_exit(code)
        assert r['implicit_callbacks']==steps and r['native_completed']==steps
        assert r['steps']==restored['steps']+steps and qualified_step(folder)==r['steps']
        for step in range(restored['steps']+1,r['steps']+1):
            report=native_load(directory/f'step-{step}.report.json');base=load(verify(load(root/f'inputs/{scope}/manifest.json')['states'][str(step)]))
            full=policy=='Off'or due(step,intent)
            parity(report,base,full);same_cp(report,checkpoint_payload(load(directory/f'step-{step}.completion.json')['checkpoint']['path'].removesuffix('/State.bin'),include_reference=True))if policy!='Off'else None
            pimple=report['pimple_attempts'][0]
            assert all(pimple[k]==1 for k in ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls'))
            assert pimple['bound_activations']==pimple['additional_correction_or_restoration_solves']==0
            assert pimple['k_bound']['raw_sha256']==pimple['k_bound']['corrected_sha256']
            if policy!='Off':
                assert report['native_full_audit_status']==('PASS'if full else'NOT_SCHEDULED')
                if not full:assert report['native_hard_audit']=={}and report['wj_audit_packets'][0]['principal_witnesses']=={}and not(directory/f'step-{step}.ready.json').exists()
                assert 'k_snapshot'not in pimple
        expected_full=1+(steps if policy=='Off'else sum(due(n,intent)for n in range(restored['steps']+1,r['steps']+1)))
        assert r['full_audit_calls']==r['physical_samples_calls']==expected_full
    atomic(folder/'qualification.json',dict(status='PASS',expected_fault=fault,actual_attempts=r['implicit_callbacks'],native_completed=r['native_completed'],verified=r['independently_qualified'],full_audit_calls=r['full_audit_calls'],physical_samples_calls=r['physical_samples_calls'],producer=identity(p),new_scientific_steps=0))
    return r


def small(root,name):
    assert name in ('legacy','A','B','fault');parent=H14/'science/small/cp0';cp0=Path(load(parent/'case.evolution.json')['checkpoint']);folder=root/'science/small'/name;restart=cp0;completion=None
    if name=='B':
        assert load(root/'science/small/A/qualification.json')['status']=='PASS';completion=root/'science/small/A/monitor/step-1.completion.json';restart=verify(load(completion)['checkpoint']).parent
    steps=2 if name in ('legacy','A')else 1;policy='Off'if name=='legacy'else'CheckpointAuditV1'
    xml_case(parent/'case.xml',folder,restart,steps,policy,fault='AfterImplicitV1'if name=='fault'else'Off')
    return coordinate(root,folder,cp0,'small',steps,policy=policy,restart_completion=completion,fault=name=='fault')


def target(root,name):
    assert name in ('A','B');v=next(v for v in load(H14/'publication-v2/state-bindings.json')['states']if v['case']=='cp0');cp0=verify(v['checkpoint']).parent;restart=cp0;completion=None
    if name=='B':completion=root/'science/target/A/monitor/step-2.completion.json';restart=verify(load(completion)['checkpoint']).parent
    folder=root/'science/target'/name;steps=4 if name=='A'else 2
    xml_case(verify(v['xml']),folder,restart,steps)
    return coordinate(root,folder,cp0,'target',steps,restart_completion=completion)


def strict(root):
    reserve_strict(root,'target')
    folder=next((root/f'science/target/{name}'for name in ('B','A')if(root/f'science/target/{name}/monitor/qualified-cursor.json').exists()),None)
    if folder is not None:cmd=strict_request(root,folder)
    else:
        v=next(v for v in load(H14/'publication-v2/state-bindings.json')['states']if v['case']=='cp0')
        folder=root/'science/target/final-original-CP0';out=folder/'view';out.mkdir(parents=True)
        atomic(out/'request.json',dict(operation='open_checks',inner_stress_model='Linear',weak=True,weak_compact=True,assembly=False,reuse=False,retain_scalar_observations=False,**cp_request(v)))
        cmd=[str(root/'runtime/RANSWallResolutionProbe'),str(verify(v['xml'])),str(out/'request.json'),str(out/'native.json')]
    with(folder/'view/native.log').open('x')as log:code=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT).returncode
    checked_exit(code);p=load(folder/'view/native.json')
    assert p['initialization_work']==[0,0]and p['initial_velocity_calls']==p['governing_solves']==0
    atomic(folder/'view/result.json',dict(status='PASS',native=identity(folder/'view/native.json'),selected_absolute_step=p['selected_absolute_step'],physical_mass_dual_norms=p['physical_mass_dual_norms'],weak=p['weak_checks'],interface=p['open_boundary']))


def admission(root):
    frozen(root)
    for n in ('legacy','A','B','fault'):assert load(root/f'science/small/{n}/qualification.json')['status']=='PASS'
    assert load(root/'science/component-checks.json')['status']=='PASS'
    runtime_qualify(root)
    book=load(root/'accounting.json');assert not any(book.get(k)for k in ('coverage_halt','runtime_halt','qualification_halt','target_halt'))
    # Six advances, two startups with complete audits/independent barriers, strict
    # finish and transport headroom. H16 measured advance/audit loop is the prior.
    forecast=1.25*(6*100+2*90+120+60)
    spent=60+300+sum(v.get('seconds',v['reserved_seconds'])for v in book['entries'].values())
    assert spent+forecast<=7200 and forecast<=3000-sum(v.get('seconds',v['reserved_seconds'])for v in book['entries'].values()if v['phase']=='qualification')
    cases={}
    for name,count in (('A',4),('B',2)):
        cmd=[sys.executable,'-B',str(HERE/'study_wj_h18.py'),'--root',str(root),'--mode','target','--case',name]
        cases[name]=dict(attempts=count,command=cmd,inputs=[identity(root/'execution-config.json'),identity(root/'science/component-checks.json')])
    atomic(root/'target-admission.json',dict(status='PASS',forecast_with_25percent_headroom=forecast,closure_reserved=300,strict_read_reserved=120,cases=cases,remaining_total=7200-spent,science_horizon_advance=0))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('setup','small','admission','target','strict','runtime'),required=True);p.add_argument('--case');a=p.parse_args()
    if a.mode=='small':small(a.root,a.case)
    elif a.mode=='target':target(a.root,a.case)
    elif a.mode=='runtime':runtime_qualify(a.root)
    else:globals()[a.mode](a.root)
