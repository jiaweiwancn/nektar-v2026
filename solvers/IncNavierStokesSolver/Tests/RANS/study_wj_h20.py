#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Three fixed continuation blocks; the existing native bootstrap is unchanged."""
import argparse,hashlib,os,subprocess,sys,time,uuid
from pathlib import Path
import xml.etree.ElementTree as ET
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,verify,atomic,vector_hash,bits,RULE
from wj_h17_observer import compatibility,load_compatible,native_load,equal
from wj_h19_analysis import spatial_xml
from wj_h20_analysis import DT,ADDITIONAL,START,END,STRIDE
from wj_h20_observer import read_anchor,evaluate
from study_wj_h16 import same_cp
from study_wj_h17 import cp_request
from wj3_measure import checkpoint_payload
from cls_epawf_m3_cases import set_info,set_parameter,write
from run_wj_h20 import frozen,reserve_strict,usage
HERE=Path(__file__).resolve().parent
NATIVE='b2365e52b0d334d7b61c6be3a4a3bb9029491084'
ACK_TIMEOUT=120.

def validate_monitor_intent(intent):
    # Match ConfigureMonitoring's execution limit before creating a governor.
    # This preflight does not replace any native validation or admission guard.
    assert intent['schema']=='WJMonitorIntentV1'and len(intent['nonce'])>=16
    assert 0<intent['ack_timeout_seconds']<=120,'native bounded ACK timeout is (0,120] seconds'
    level=intent['level']
    assert intent['maximum_absolute_step']==END[level]and intent['restored_absolute_step']==START[level]
    assert intent['recovery_stride']==1 and intent['scientific_stride']==STRIDE[level]
    assert intent['dt_bits']==bits(DT[level])

def executable(root,name):
    b=load(root/'science/runtime-build-v1.json')
    return verify(next(v['wrapper']for v in b['executables']if Path(v['wrapper']['path']).name==name))

def xml_case(original,folder,cp,level,steps):
    assert level in DT and steps in (0,ADDITIONAL[level]);folder.mkdir(parents=True)
    x=ET.parse(original).getroot()
    for k,v in (('TimeStep',DT[level]),('NumSteps',steps),('IO_CheckSteps',1),('IO_InfoSteps',1)):set_parameter(x,k,v)
    for k,v in (('RANSWallRestart',str(cp)),('RANSWallCheckpoint',str(folder/'CP')),('RANSWallRetention','BoundedV1'if steps else'Off'),
        ('RANSWallMonitorIntent',str(folder/'intent.json')),('RANSWallMonitorFault','Off'),('RANSWallReport',str(folder/'native.json')),
        ('RANSWallImmutableRowReuse','Off'),('RANSWallMonitoringPolicy','CheckpointAuditV1'if steps else'Off'),
        ('RANSWallPimpleSnapshots','Off'),('RANSWallPimpleArithmetic','Off'),('RANSWallPimpleSnapshotDirectory',str(folder/'snapshots'))):set_info(x,k,v)
    assert spatial_xml(x)==spatial_xml(ET.parse(original).getroot()),'only timestep and declared execution controls may differ'
    assert next(float(p.text.split('=')[1])for p in x.findall('./CONDITIONS/PARAMETERS/P')if p.text.split('=')[0].strip()=='Time')==0
    p=folder/'case.xml';write(x,p)
    atomic(folder/'execution-diff.json',dict(original=identity(original),actual=identity(p),physical_except_dt_sha256=hashlib.sha256(spatial_xml(x)).hexdigest(),
        TimeStep=DT[level],Time=0,additional_steps=steps,absolute_ceiling=END[level],native_signature_exclusions_changed=False))
    return p

def setup(root):
    inputs=load(root/'inputs/selected.json');b=load(verify(inputs['runtime_build']))
    assert b['native_source_revision']==NATIVE
    for v in b['executables']:verify(v['ELF']);verify(v['wrapper'])
    deps=load(verify(inputs['runtime_dependencies']))
    for v in deps['libraries']:verify(v['actual'])
    compiled=load(verify(inputs['compiled']))
    for v in compiled['dependencies']:verify(v)
    atomic(root/'science/runtime-build-v1.json',dict(b,package='H20',reuse='qualified H19 private closure in place; no build',historical=inputs['runtime_build']))
    atomic(root/'science/runtime-dependencies-v1.json',dict(deps,historical=inputs['runtime_dependencies']))
    historical=load(verify(inputs['cache_manifest']))
    compatibility(root,verify(historical['historical_manifest']),root/'science/cache-base-compatibility.json')
    c=load(root/'science/cache-base-compatibility.json')
    c.update(schema='WJH20DtIndependentCacheCompatibilityV1',permitted_difference='dt .005/.0025/.00125; independent immutable spatial rows/reference only',
        dt_independence=dict(status='PASS',timesteps=DT,cached_objects='geometry, basis/signs/maps, physical quadrature/tests, analytical profile, identical original spatial qh0/reference/defect',
            exclusions='no state, history, alpha, closure, profile/cutoff, boundary reaction or residual reuse',own_signature_checks=True))
    atomic(root/'science/target-cache-compatibility.json',c)
    for level,states in inputs['states'].items():
        for n,v in states.items():read_anchor(v,level,int(n)*DT[level])
    source=checkpoint_payload(verify(inputs['original']['checkpoint']).parent,include_reference=True)
    coarse=checkpoint_payload(verify(inputs['C_CP0']).parent,include_reference=True)
    for k in ('state','fixed_reference','fixed_defect','eta'):equal(source[k],coarse[k],'original/C spatial '+k)
    assert source['steps']==coarse['steps']==0 and bits(source['time'])==bits(coarse['time'])==bits(0.)
    assert source['dt']==DT['M']and coarse['dt']==DT['C']and source['signature']!=coarse['signature']
    atomic(root/'science/initial-spatial-lineage.json',dict(status='PASS',original=inputs['original']['checkpoint'],C=inputs['C_CP0'],
        hashes={k:vector_hash(source[k])for k in ('state','fixed_reference','fixed_defect')},eta=source['eta'],source_signature=source['signature'],C_signature=coarse['signature'],history='each own timestep; no history conversion'))
    files=[identity(HERE/n)for n in ('run_wj_h20.py','study_wj_h20.py','wj_h20_analysis.py','wj_h20_observer.py','test_wj_h20_three_level.py')]
    files.extend(identity(root/'science'/n)for n in ('runtime-build-v1.json','runtime-dependencies-v1.json','target-cache-compatibility.json','initial-spatial-lineage.json'))
    files.append(identity(root/'inputs/selected.json'))
    if(root/'resumption.json').exists():files.append(identity(root/'resumption.json'))
    for v in b['executables']:files.extend((v['ELF'],v['wrapper']))
    for v in c['active_sources']:files.append(v['active'])
    atomic(root/'execution-config.json',dict(schema='WJH20FrozenExecutionV1',files=files,order=['M','C','F'],maximum_CP0=1,maximum_blocks=3,maximum_attempts=22,maximum_final_strict=3,
        native_source=NATIVE,policies='CheckpointAuditV1/BoundedV1/ExactStepV1/ImmutableRowReuse Off',DT=DT,additional=ADDITIONAL,absolute_ceilings=END,scientific_strides=STRIDE,
        matrix='M physical only; C absolute C-M; F all-pair absolute and descriptive contraction',contraction_halts=False))

def bootstrap(root):
    frozen(root);inputs=load(root/'inputs/selected.json');v=inputs['original'];folder=root/'science/F-cp0'
    xml=xml_case(verify(v['xml']),folder,'None','F',0)
    request=dict(inner_stress_model='Linear',source_dt=DT['M'],source_xml=v['xml'],checkpoint_base=str(folder/'CP'),**cp_request(v))
    request['operation']='step_zero_bootstrap';atomic(folder/'request.json',request)
    with(folder/'native.log').open('x')as log:code=subprocess.run([str(executable(root,'RANSWallResolutionProbe')),str(xml),str(folder/'request.json'),str(folder/'native.json')],stdout=log,stderr=subprocess.STDOUT).returncode
    assert code==0,'one F CP0 reservation consumed; no retry'
    r=native_load(folder/'native.json');source=checkpoint_payload(verify(v['checkpoint']).parent,include_reference=True);fine=checkpoint_payload(r['checkpoint'],include_reference=True)
    for k in ('state','fixed_reference','fixed_defect','eta'):equal(source[k],fine[k],'F step-zero spatial '+k)
    assert source['steps']==fine['steps']==0 and bits(source['time'])==bits(fine['time'])==bits(0.)
    assert source['dt']==DT['M']and fine['dt']==DT['F']and fine['signature']!=source['signature']
    assert r['native_history_initializations']==1 and r['initialization_work']==[0,0]and r['initial_velocity_calls']==r['governing_solves']==0
    cache,_=load_compatible(root/'science/target-cache-compatibility.json');fixed=load(verify(inputs['common_support']))['points']
    result,inventory=evaluate(root,cache,fixed,inputs,r,'F',0,folder,'cp0')
    assert result['status']=='PASS'and result['continuation_verdict']=='PROCEED'
    atomic(folder/'lineage.json',dict(status='PASS',source=v,checkpoint=identity(Path(r['checkpoint'])/'State.bin'),report=identity(folder/'native.json'),xml=identity(xml),
        source_dt=source['dt'],dt=fine['dt'],source_signature=source['signature'],signature=fine['signature'],
        hashes={k:vector_hash(fine[k])for k in ('state','fixed_reference','fixed_defect')},eta=fine['eta'],native_history_initializations=1,spatial_initializations=0,
        native_history_times=fine['native_history_times'],native_history_sha256=[vector_hash(h)for h in fine['native_history']],spatial_lineage_anchor=v['checkpoint']))
    atomic(folder/'qualification.json',dict(status='PASS',lineage=identity(folder/'lineage.json'),physical=result,inventory=identity(inventory)))

def anchors(root,level):
    inputs=load(root/'inputs/selected.json')
    if level=='F':
        v=load(root/'science/F-cp0/lineage.json');return v,dict(checkpoint=v['checkpoint'],report=v['report'],xml=v['xml'])
    v=inputs['states'][level][str(START[level])]
    cp0=dict(checkpoint=inputs['C_CP0'],report=inputs['C_CP0_report'],xml=inputs['C_CP0_xml'])if level=='C'else dict(checkpoint=inputs['original']['checkpoint'],report=inputs['original']['producer'],xml=inputs['original']['xml'])
    return v,cp0

def preflight(root):
    """Retained-state evaluator controls, never a new target solution."""
    frozen(root);inputs=load(root/'inputs/selected.json');folder=root/'science/passive-preflight';folder.mkdir()
    cache,_=load_compatible(root/'science/target-cache-compatibility.json');fixed=load(verify(inputs['common_support']))['points']
    m,_=read_anchor(inputs['states']['M']['4'],'M',.01)
    result,_=evaluate(root,cache,fixed,inputs,m,'M',4,folder,'retained-M4')
    assert result['status']=='PASS'and result['comparison']['status']=='NOT_YET_AVAILABLE'
    c,_=read_anchor(inputs['states']['C']['2'],'C',.01)
    result,_=evaluate(root,cache,fixed,inputs,c,'C',2,folder,'retained-C2-M4')
    assert result['status']=='PASS'and result['comparison']['status']=='PASS'
    # Exercise full triplet serialization with actual immutable state aliases.
    # It is an alias fixture, not F evidence or a temporal measurement.
    import numpy as np
    from wj_h20_analysis import triplet_diagnostics
    arrays=np.load(verify(result['signed_arrays']));inv=load(verify(result['inventory']));rows,_,_=cache.view(inv['points'])
    data={z:{k[2:]:arrays[k]for k in arrays.files if k.startswith(z+'_')}for z in ('C','M')};data['F']=data['M']
    diagnostic,_=triplet_diagnostics(rows,data)
    assert len(diagnostic['fixed_probes'])>=4
    atomic(folder/'result.json',dict(status='PASS',retained_pair_result=identity(folder/'retained-C2-M4.result.json')if(folder/'retained-C2-M4.result.json').exists()else result['signed_arrays'],
        retained_M_summary='PASS',retained_C_M_screen='PASS',triplet_alias_fixture='PASS, F=M only here; not scientific F evidence',new_native_work=0))

def admission(root,level):
    frozen(root);b=load(root/'accounting.json')
    assert not any(b.get(k)for k in ('runtime_halt','qualification_halt','target_halt','coverage_halt','screen_halt'))
    assert load(root/'operations/focused-controls-v1/receipt.json')['status']=='PASS'
    previous={'M':None,'C':'M','F':'C'}[level]
    if previous:assert load(root/f'science/{previous}/qualification.json')['status']=='PASS'
    if level=='F':assert load(root/'science/F-cp0/qualification.json')['status']=='PASS'
    # H19 startup.initialize already includes its native full audit. Budget
    # 110s inclusive startup, then 50s online plus 44s audit and 32s/44s observer
    # at NEW due steps; do not count the startup full audit a second time.
    observer=32 if level!='F'else 44
    remaining=['M','C','F'][['M','C','F'].index(level):]
    estimate=sum(110+ADDITIONAL[z]*50+(ADDITIONAL[z]//STRIDE[z])*(44+(44 if z=='F'else 32))+120+20 for z in remaining)+(140 if level!='F'else 0)
    headroom=1.25*estimate
    used=usage(b)
    assert used['total']+headroom<=7200,'global measured remaining admission'
    assert used['native']+headroom<=3600,'combined native/bootstrap/strict remaining admission'
    if(root/'resumption.json').exists():assert load(root/'science/passive-preflight/result.json')['status']=='PASS'
    restored,cp0=anchors(root,level)
    if level!='F':read_anchor(restored,level,.01)
    required=[identity(root/'execution-config.json'),restored['checkpoint'],restored['report'],cp0['checkpoint'],cp0['xml']]
    if 'callback'in restored:required.extend((restored['callback'],restored['completion'],restored['qualified_commit']))
    atomic(root/f'admission-{level}.json',dict(status='PASS',level=level,remaining_estimate=estimate,forecast_with_25percent_headroom=headroom,closure_reserved=300,final_strict_reservations=120*len(remaining),
        command=[sys.executable,'-B',str(HERE/'study_wj_h20.py'),'--root',str(root),'--mode','target','--level',level],inputs=required))

def target(root,level):
    frozen(root);inputs=load(root/'inputs/selected.json');restored,cp0=anchors(root,level)
    folder=root/'science'/level;xml=xml_case(verify(cp0['xml']),folder,verify(restored['checkpoint']).parent,level,ADDITIONAL[level])
    directory=folder/'monitor';directory.mkdir();ip=folder/'intent.json'
    intent=dict(schema='WJMonitorIntentV1',policy='CheckpointAuditV1',level=level,directory=str(directory),nonce=uuid.uuid4().hex,
        scope_sha256=hashlib.sha256((identity(xml)['sha256']+identity(root/'execution-config.json')['sha256']).encode()).hexdigest(),
        source_sha256=identity(root/'science/runtime-build-v1.json')['sha256'],observer_sha256=identity(HERE/'wj_h20_observer.py')['sha256'],inventory_rule_sha256=hashlib.sha256(RULE.encode()).hexdigest(),
        original_qh0_checkpoint=cp0['checkpoint'],original_qh0_state_sha256=inputs['original']['coupled_state_sha256'],spatial_lineage_anchor=inputs['original']['checkpoint'],restored_checkpoint=restored['checkpoint'],
        maximum_absolute_step=END[level],restored_absolute_step=START[level],recovery_stride=1,scientific_stride=STRIDE[level],dt_bits=bits(DT[level]),ack_timeout_seconds=ACK_TIMEOUT,strict_reserve_seconds=120,
        ordinary_step_forecast_seconds=50,observer_forecast_seconds=44 if level=='F'else 32,deadline_monotonic=time.monotonic()+1750,
        substep_evidence_steps=[],preserve_qualification_artifacts=True,common_support=inputs['common_support'],selected_inputs=identity(root/'inputs/selected.json'))
    if START[level]:intent.update(restart_completion=restored['completion'],restart_callback=restored['callback'])
    validate_monitor_intent(intent)
    obs=native=None;code=None;error=None;start=time.monotonic()
    with(folder/'observer.log').open('x')as ol,(folder/'native.log').open('x')as nl:
        try:
            obs=subprocess.Popen([sys.executable,'-B',str(HERE/'wj_h20_observer.py'),'--root',str(root),'--intent',str(ip),'--compat',str(root/'science/target-cache-compatibility.json')],stdout=ol,stderr=subprocess.STDOUT)
            intent['observer_pid']=obs.pid;atomic(ip,intent)
            while not(directory/'observer-started.json').exists():assert obs.poll()is None and time.monotonic()-start<30;time.sleep(.025)
            native=subprocess.Popen([str(executable(root,'IncNavierStokesSolver')),str(xml)],stdout=nl,stderr=subprocess.STDOUT)
            atomic(folder/'native-launch.json',dict(pid=native.pid,xml=identity(xml),intent=identity(ip),additional_steps=ADDITIONAL[level],restored_absolute_step=START[level]))
            while native.poll()is None:
                assert obs.poll()is None,'observer exited; stop governor'
                time.sleep(.025)
            code=native.returncode
        except Exception as e:error=repr(e)
        finally:
            if native and native.poll()is None:native.terminate();native.wait(timeout=5)
            atomic(directory/'observer-stop.json',dict(reason='H20 coordinator exit'))
            if obs:
                try:obs.wait(timeout=5)
                except subprocess.TimeoutExpired:obs.terminate();obs.wait(timeout=5)
    report=folder/'case.evolution.json';r=native_load(report)if report.exists()else None
    atomic(folder/'outcome.json',dict(native_exit=code,error=error,producer=identity(report)if r else None,seconds=time.monotonic()-start))
    if r:
        op=os.environ.get('WJ_H20_OPERATION')
        if op:atomic(Path(op)/'attempt-counts.json',dict(kind=level,actual_attempts=r['implicit_callbacks'],producer=identity(report)))
        assert r['initialization_work']==[0,0]and r['initializer_calls']==0
        for step in range(START[level]+1,r['steps']+1):
            rr=native_load(directory/f'step-{step}.report.json');p=rr['pimple_attempts'][0]
            assert all(p[k]==1 for k in ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls'))
            assert p['bound_activations']==p['additional_correction_or_restoration_solves']==0 and p['k_bound']['raw_sha256']==p['k_bound']['corrected_sha256']
            due=step%STRIDE[level]==0 or step==END[level]
            assert rr['native_full_audit_status']==('PASS'if due else'NOT_SCHEDULED')
            if not due:assert rr['native_hard_audit']=={}and not(directory/f'step-{step}.ready.json').exists()
            same_cp(rr,checkpoint_payload(verify(load(directory/f'step-{step}.completion.json')['checkpoint']).parent,include_reference=True))
    assert error is None,error
    if code:
        reason=r['failure']if r else next((s.strip()for s in (folder/'native.log').read_text().splitlines()if s.startswith('EP-AWF evolution ')),'native nonzero exit; retained native.log')
        atomic(folder/'first-native-failure.json',dict(native_exit=code,reason=reason,producer=identity(report)if r else None,native_log=identity(folder/'native.log')))
        if not(r and reason.startswith('TEMPORAL_REVIEW_STOP')and load(root/'accounting.json').get('screen_halt')):raise RuntimeError(reason)
    else:
        assert r and r['steps']==END[level]and r['implicit_callbacks']==r['native_completed']==ADDITIONAL[level]
        assert load(directory/'latest-fully-verified-checkpoint.json')['binding']['absolute_step']==END[level]
    atomic(folder/'qualification.json',dict(status='PASS'if not code else'PROVISIONAL_SCREEN_STOP',native_exit=code,native_completed=r['native_completed'],independently_qualified=r['independently_qualified'],
        attempts=r['implicit_callbacks'],final_absolute_step=r['steps'],time=r['time'],full_audit_calls=r['full_audit_calls'],physical_samples_calls=r['physical_samples_calls'],producer=identity(report)))

def strict(root,level):
    reserve_strict(root,level);folder=root/'science'/level;directory=folder/'monitor';p=directory/'latest-fully-verified-checkpoint.json'
    if p.exists()and load(p)['schema']=='WJQualificationCommitV1':
        c=load(p);step=int(c['binding']['absolute_step']);commit=directory/f'step-{step}.qualification-commit.json'
        request=dict(operation='open_checks',inner_stress_model='Linear',native_report=c['report']['path'],native_report_sha256=c['report']['sha256'],native_checkpoint=str(verify(c['checkpoint']).parent),
            checkpoint_sha256=c['checkpoint']['sha256'],completed_sidecar=c['callback'],qualified_cursor=identity(commit));xml=folder/'case.xml'
    else:
        v,cp0=anchors(root,level)
        if START[level]:
            read_anchor(v,level,.01);c=load(verify(v['qualified_commit']))
            request=dict(operation='open_checks',inner_stress_model='Linear',**cp_request(dict(v,producer=v['report'])),completed_sidecar=v['callback'],qualified_cursor=v['qualified_commit']);xml=folder/'case.xml'
        else:
            request=dict(operation='open_checks',inner_stress_model='Linear',**cp_request(dict(producer=cp0['report'],checkpoint=cp0['checkpoint'])));xml=verify(cp0['xml'])
    request.update(weak=True,weak_compact=True,assembly=False,reuse=False,retain_scalar_observations=False)
    out=folder/'final-strict';out.mkdir();atomic(out/'request.json',request)
    with(out/'native.log').open('x')as log:code=subprocess.run([str(executable(root,'RANSWallResolutionProbe')),str(xml),str(out/'request.json'),str(out/'native.json')],stdout=log,stderr=subprocess.STDOUT).returncode
    assert code==0;r=native_load(out/'native.json');assert r['initialization_work']==[0,0]and r['initial_velocity_calls']==r['governing_solves']==0
    atomic(out/'result.json',dict(status='PASS',selected_absolute_step=r['selected_absolute_step'],producer=identity(out/'native.json'),physical_mass_dual_norms=r['physical_mass_dual_norms'],weak=r['weak_checks'],interface=r['open_boundary']))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('setup','preflight','bootstrap','admission','target','strict'),required=True);p.add_argument('--level',choices=tuple(DT));a=p.parse_args()
    if a.mode in ('admission','target','strict'):assert a.level;globals()[a.mode](a.root,a.level)
    else:globals()[a.mode](a.root)
