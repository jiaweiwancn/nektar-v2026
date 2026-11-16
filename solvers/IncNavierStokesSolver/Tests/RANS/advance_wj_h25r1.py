#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Genuine C-level continuation and paired sparse/dense restart pilots."""
import argparse,hashlib,importlib,os,subprocess,sys,time,uuid
from pathlib import Path
import xml.etree.ElementTree as ET
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,verify,atomic,vector_hash,bits,RULE
from wj_h16_observer import compatibility
from wj_h17_observer import native_load,equal
from wj_h18_observer import due,rotate_owned
from wj3_measure import checkpoint_payload
from cls_epawf_m3_cases import set_info,set_parameter,write

H24=Path('/home/jay/build/wj-h24-flat-continuation-20261115')


def run(root,mode):
    base=root/'resumption-r1';selected=load(H24/'inputs/selected.json')
    assert load(base/'science/W1-scatter-v1/qualification.json')['status']=='PASS'
    label={'wide':'C-wide-v1','sparse':'C-sparse-restart-v1','dense':'C-dense-restart-v1'}[mode]
    d=base/'science'/label;d.mkdir();monitor=d/'monitor';monitor.mkdir();ip=d/'intent.json'
    if mode=='wide':
        anchor=selected['anchors']['C'];steps=8;stride=2
        expected={n:H24/f'science/C/monitor/step-{n}.report.json' for n in range(5,9)}
    else:
        original=base/'science/C-wide-v1/monitor';q=load(original/'step-8.qualification-commit.json')
        anchor={k:q[k]for k in ('checkpoint','report','callback')}
        anchor.update(qualified_commit=identity(original/'step-8.qualification-commit.json'),completion=identity(original/'step-8.completion.json'))
        assert load(base/'science/C-wide-v1/temporal-qualification.json')['status']=='PASS_SHORT_UNPERTURBED_WINDOW'
        steps=4;stride=4 if mode=='sparse' else 1
        expected={n:original/f'step-{n}.report.json' for n in range(9,13)}
    for pin in anchor.values():verify(pin)
    cp=checkpoint_payload(verify(anchor['checkpoint']).parent,include_reference=True)
    cp0=selected['cp0']['C'];qh0=checkpoint_payload(verify(cp0['checkpoint']).parent,include_reference=True)
    assert cp['dt']==qh0['dt']==.005
    start=int(cp['steps']);terminal=start+steps
    runtime_path=base/'science/runtime-build-scatter-v1.json';runtime=load(runtime_path)
    native=verify(next(e['wrapper']for e in runtime['executables']if Path(e['wrapper']['path']).name=='IncNavierStokesSolver'))
    observer_path=Path(importlib.import_module('wj_h25_observer').__file__).resolve()
    c=compatibility(root,verify(load(root/'science/cache-compatible-v2.json')['historical_manifest']),d/'cache-binding-base.json')
    c.update(new_native_identity=identity(runtime_path),permitted_difference='exact element scatter and tensor lookup; own .005 signature/history; immutable FE geometry/reference/qh0 coefficients unchanged')
    atomic(d/'cache-compatible.json',c)
    x=ET.parse(verify(cp0['xml'])).getroot()
    for k,v in dict(NumSteps=steps,IO_CheckSteps=1,IO_InfoSteps=1).items():set_parameter(x,k,v)
    for k,v in dict(RANSWallRestart=str(verify(anchor['checkpoint']).parent),RANSWallCheckpoint=str(d/'CP'),
        RANSWallRetention='BoundedV1',RANSWallMonitoringPolicy='CheckpointAuditV1',RANSWallMonitorIntent=str(ip),
        RANSWallMonitorFault='Off',RANSWallReport=str(d/'native.json'),RANSWallPimpleSnapshots='Off',
        RANSWallPimpleFlowVectors='Off',RANSWallPimpleArithmetic='Off',RANSWallPimpleVerify='Off',
        RANSWallImmutableRowReuse='Off',RANSWallNativeLookup='TensorExactV1',RANSWallScalarEvaluation='Off',
        RANSWallAssemblyScatter='ElementExactV1').items():set_info(x,k,v)
    write(x,d/'case.xml')
    intent=dict(schema='WJMonitorIntentV1',policy='CheckpointAuditV1',directory=str(monitor),nonce=uuid.uuid4().hex,
        scope_sha256=identity(d/'case.xml')['sha256'],source_sha256=identity(runtime_path)['sha256'],
        observer_sha256=identity(observer_path)['sha256'],inventory_rule_sha256=hashlib.sha256(RULE.encode()).hexdigest(),
        original_qh0_checkpoint=cp0['checkpoint'],original_qh0_state_sha256=vector_hash(qh0['state']),
        spatial_lineage_anchor=cp0['checkpoint'],completed_anchor=anchor,restored_checkpoint=anchor['checkpoint'],
        restart_completion=anchor['completion'],maximum_absolute_step=terminal,restored_absolute_step=start,
        recovery_stride=1,scientific_stride=stride,dt_bits=bits(.005),ack_timeout_seconds=120.,
        substep_evidence_steps=[],flow_capture='Off',snapshot_capture='Off',preserve_qualification_artifacts=True,
        cache_compatibility=identity(d/'cache-compatible.json'),common_support=selected['common_support'])
    governor=observer=None;code=None;error=None;began=time.monotonic()
    with(d/'native.log').open('x')as nl,(d/'observer.log').open('x')as ol:
        try:
            observer=subprocess.Popen([sys.executable,'-B',str(observer_path),'--root',str(root),'--intent',str(ip)],stdout=ol,stderr=subprocess.STDOUT)
            intent['observer_pid']=observer.pid;atomic(ip,intent)
            while not(monitor/'observer-started.json').exists():
                assert observer.poll()is None and time.monotonic()-began<30
                time.sleep(.025)
            governor=subprocess.Popen([str(native),str(d/'case.xml')],stdout=nl,stderr=subprocess.STDOUT)
            while governor.poll()is None:
                assert observer.poll()is None,'observer exited before native'
                time.sleep(.025)
            code=governor.returncode
        except Exception as e:error=repr(e)
        finally:
            if governor and governor.poll()is None:governor.terminate();governor.wait(timeout=5)
            atomic(monitor/'observer-stop.json',dict(native_exit=code,error=error))
            if observer:
                try:observer.wait(timeout=5)
                except subprocess.TimeoutExpired:observer.terminate();observer.wait(timeout=5)
    report=d/'case.evolution.json';r=native_load(report)if report.exists()else None
    actual=int(r['implicit_callbacks'])if r else steps
    atomic(Path(os.environ['WJ_H25_OPERATION'])/'attempt-counts.json',dict(actual_attempts=actual,known=r is not None))
    atomic(d/'outcome.json',dict(native_exit=code,error=error,seconds=time.monotonic()-began,producer=identity(report)if r else None))
    assert code==0 and error is None
    assert r['initialization_work']==[0,0] and r['initializer_calls']==0
    records=[]
    for n in range(start+1,terminal+1):
        current=native_load(monitor/f'step-{n}.report.json')
        if n in expected:
            old=native_load(expected[n])
            for k in ('state','fixed_reference','fixed_defect','native_history','native_history_times','time','steps','dt','signature','last_stage_history','last_stage_alpha','last_stage_time'):
                equal(current[k],old[k],f'{mode} genuine continuation '+k)
        p=current['pimple_attempts'][0]
        assert all(p[k]==1 for k in ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls'))
        assert p['bound_activations']==p['additional_correction_or_restoration_solves']==0
        records.append(dict(step=n,time=current['time'],producer=identity(monitor/f'step-{n}.report.json'),
            oracle=identity(expected[n])if n in expected else None,state=vector_hash(current['state'])))
    cursor=load(monitor/'latest-fully-verified-checkpoint.json');assert int(cursor['binding']['absolute_step'])==terminal
    atomic(d/'qualification.json',dict(status='PASS_NATIVE_AND_PHYSICAL_PENDING_TEMPORAL'if mode=='wide'else'PASS_UNPERTURBED_RESTART_PILOT',
        mode=mode,records=records,source_anchor=anchor,start=start,terminal=terminal,dt=.005,
        native_resource_summary={k:v for k,v in r.items()if k.startswith('resident_')},
        endpoint=identity(monitor/f'step-{terminal}.qualification-commit.json'),
        limitation='unperturbed evolving flat control; prescribed finite perturbation failed initial quality and is not replaced by this pilot'))
    print(label,'native/physical PASS',flush=True)


def strict(root,label):
    d=root/'resumption-r1/science'/label;v=load(d/'qualification.json');terminal=v['terminal']
    q=load(d/f'monitor/step-{terminal}.qualification-commit.json');out=d/'strict';out.mkdir()
    request=dict(operation='open_checks',inner_stress_model='Linear',import_role='QualifiedV1',
        native_report=q['report']['path'],native_report_sha256=q['report']['sha256'],
        native_checkpoint=str(verify(q['checkpoint']).parent),checkpoint_sha256=q['checkpoint']['sha256'],
        completed_sidecar=q['callback'],qualified_cursor=identity(d/f'monitor/step-{terminal}.qualification-commit.json'),
        weak=True,weak_compact=True,assembly=False,reuse=False,retain_scalar_observations=False)
    atomic(out/'request.json',request)
    runtime=load(root/'resumption-r1/science/runtime-build-scatter-v1.json')
    probe=verify(next(e['wrapper']for e in runtime['executables']if Path(e['wrapper']['path']).name=='RANSWallResolutionProbe'))
    with(out/'native.log').open('x')as f:subprocess.run([str(probe),str(d/'case.xml'),str(out/'request.json'),str(out/'native.json')],stdout=f,stderr=subprocess.STDOUT,check=True)
    r=native_load(out/'native.json')
    assert r['initialization_work']==[0,0] and r['initial_velocity_calls']==r['governing_solves']==0
    atomic(out/'qualification.json',dict(status='PASS',request=identity(out/'request.json'),native=identity(out/'native.json'),
        step=terminal,time=native_load(verify(q['report']))['time']))


def rotate(root,label):
    d=root/'resumption-r1/science'/label;assert load(d/'strict/qualification.json')['status']=='PASS'
    v=load(d/'qualification.json');intent=load(d/'intent.json');records=[]
    for n in range(v['start']+1,v['terminal']+1):
        completion=d/f'monitor/step-{n}.completion.json';c=load(completion)
        if due(n,intent):continue
        records.append(dict(id=str(n),step=n,role='recovery',verification='NOT_REQUIRED_SUPERSEDED',
            payloads=[c[k]for k in ('checkpoint','report','callback')],completion=identity(completion)))
    removed=rotate_owned(d,records,protected=[],keep=2)
    atomic(d/'rotation.json',dict(status='PASS',removed=removed,records=records,keep_latest_superseded=2,
        retained_scientific_terminal=v['terminal'],source_anchor_preserved=v['source_anchor'],
        distinction='only owned, superseded off-cadence recovery payloads; immutable completion metadata retained'))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    p.add_argument('--mode',choices=('wide','sparse','dense','strict','rotate'),required=True);p.add_argument('--label')
    a=p.parse_args()
    if a.mode=='strict':strict(a.root,a.label)
    elif a.mode=='rotate':rotate(a.root,a.label)
    else:run(a.root,a.mode)
