#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""One pending M2 startup verification and the four missing F replay attempts."""
import argparse,copy,hashlib,os,subprocess,sys,time,uuid
from pathlib import Path
import xml.etree.ElementTree as ET
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,verify,atomic,vector_hash,bits,RULE
from wj_h17_observer import native_load,equal,compatibility
from study_wj_h16 import same_cp
from study_wj_h17 import runtime_qualify
from study_wj_h22 import executable,protected_case
from cls_epawf_m3_cases import set_parameter,set_info,write
from wj_h22_audit_roles import select_auxiliary_audits,native_full_map_parity
from wj_h22_replay import DT,END,recovery,read_completed,parity
from wj_h22_observer import prefix
from wj_h22_completion import validate_role,restored_parity,typed_strict,capture_descriptor
from wj3_measure import checkpoint_payload
from run_wj_h22_completion import PARENT,frozen,usage,ROOT

HERE=Path(__file__).resolve().parent

def xml_case(original,folder,cp,level,steps,capture):
    folder.mkdir(parents=True);x=ET.parse(original).getroot()
    for k,v in (('NumSteps',steps),('IO_CheckSteps',1),('IO_InfoSteps',1)):set_parameter(x,k,v)
    for k,v in (('RANSWallRestart',str(cp)),('RANSWallCheckpoint',str(folder/'CP')),('RANSWallRetention','BoundedV1'),
        ('RANSWallMonitorIntent',str(folder/'intent.json')),('RANSWallMonitorFault','Off'),('RANSWallReport',str(folder/'native.json')),
        ('RANSWallImmutableRowReuse','Off'),('RANSWallMonitoringPolicy','CheckpointAuditV1'),('RANSWallPimpleSnapshots','Off'),
        ('RANSWallPimpleArithmetic','Off'),('RANSWallPimpleVerify','Off'),('RANSWallPimpleSnapshotDirectory',str(folder/'snapshots')),
        ('RANSWallPimpleFlowVectors',capture),('RANSWallPimpleFlowVectorIdentity',str(folder.parent/'flow-vector-execution.json'))):set_info(x,k,v)
    assert protected_case(x)==protected_case(ET.parse(original).getroot()),'physical XML/TimeStep changed'
    params={p.text.split('=')[0].strip():p.text.split('=')[1].strip()for p in x.findall('./CONDITIONS/PARAMETERS/P')}
    assert float(params['Time'])==0 and float(params['TimeStep'])==DT[level]and int(params['NumSteps'])==steps
    dest=folder/'case.xml';write(x,dest)
    atomic(folder/'execution-diff.json',dict(original=identity(original),actual=identity(dest),TimeStep='UNCHANGED',Time='UNCHANGED ZERO',
        additional_NumSteps=steps,protected_sha256=hashlib.sha256(protected_case(x)).hexdigest(),physical_signature_definition_changed=False,
        permitted='NumSteps/IO and existing restart/report/monitor/snapshot/passive vector SolverInfo only'))

def setup(root):
    parent_index=load(root/'inputs/parent-index.json')
    for pin in parent_index.values():verify(pin)
    assert load(verify(parent_index['accounting']))['closed']
    selected=load(PARENT/'inputs/selected.json');selected['auxiliary_audits']=select_auxiliary_audits(selected)
    retained={};pins={}
    def remember(pin):verify(pin);pins[pin['path']]=pin
    for z in ('C','M'):
        retained[z]={}
        for n in range(1,END[z]+1):
            anchor=recovery(PARENT/f'science/{z}/monitor',n)
            q=dict(original_execution_qualification='COMMITTED'if z=='C'else'NOT_COMMITTED'if n==2 else'NOT_SCHEDULED',
                   subsequent_startup_verification='NOT_YET_RUN'if z=='M'and n==2 else'NOT_APPLICABLE')
            if z=='C':q['original_commit']=identity(PARENT/'science/C/monitor/step-1.qualification-commit.json')
            d=capture_descriptor(anchor,identity(PARENT/'science/flow-vector-execution.json'),q);retained[z][str(n)]=d
            for p in (*anchor.values(),d['execution_identity'],d['phase_snapshot'],d['sparse_k_matrix'],d['vector_packet'],
                      native_load(verify(d['vector_packet']))['payload']):remember(p)
    selected['retained_attempts']=retained;selected['parent_index']=identity(PARENT/'authoritative-index-v1.json')
    atomic(root/'inputs/selected.json',selected)
    for relative in ('authoritative-index-v1.json','accounting.json','execution-config.json','execution-config-v2.json',
        'publication-v1/manifest.json','closure-v1/manifest.json','closure-v1/accounting-final.json',
        'science/M/case.evolution.json','science/M/monitor/first-replay-capture-failure.json','science/M/monitor/step-2.ack.json',
        'science/M/monitor/latest-fully-verified-checkpoint.json','science/C/monitor/step-1.qualification-commit.json',
        'science/runtime-build-v1.json','science/runtime-dependencies-v1.json','science/compiled-source-objects.json'):
        remember(identity(PARENT/relative))
    atomic(root/'inputs/original-immutability.json',dict(files=list(pins.values()),
        absent=[str(PARENT/'science/M/monitor/step-2.qualification-commit.json')],old_failure='RETAINED_OBSTRUCTED',new_verification='NOT_YET_RUN'))
    for name in ('runtime-build-v1.json','runtime-dependencies-v1.json','compiled-source-objects.json'):
        value=load(PARENT/'science'/name);atomic(root/'science'/name,value)
    b=load(root/'science/runtime-build-v1.json');assert b['native_source_revision']=='057e6b7465660cac5ec8617abf5c876fa8bff71f'
    for v in b['executables']:verify(v['ELF']);verify(v['wrapper'])
    verify(b['compile_commands'])
    for v in b['compiled_inputs']:
        assert identity(v['path'])['sha256']==v['sha256']
        assert '-ffp-contract=off'in v['command']and '-fno-fast-math'in v['command']and '-O3'in v['command']
    objects=load(root/'science/compiled-source-objects.json')
    for v in objects['objects']:
        verify(v['source']);verify(v['object'])
        if v['dependency_file']:verify(v['dependency_file'])
        if v.get('link_command'):verify(v['link_command'])
    for v in objects['local_dependencies']+objects['static_libraries']:verify(v)
    for v in load(root/'science/runtime-dependencies-v1.json')['libraries']:verify(v['actual'])
    previous=load(PARENT/'science/target-cache-compatibility.json')
    c=compatibility(root,verify(previous['historical_manifest']),root/'science/cache-base-compatibility.json')
    c.update(schema='WJH22CompletionCacheCompatibilityV1',permitted_difference='Python producer/audit/verification/source adapters only; same native ELFs and unchanged numerical evaluator/cache',
        native_rebuilds=0,historical_pins_modified=False)
    atomic(root/'science/target-cache-compatibility.json',c)
    atomic(root/'science/flow-vector-execution.json',dict(schema='WJH22CompletionFlowExecutionV1',native_source=b['native_source_revision'],
        runtime_build=identity(root/'science/runtime-build-v1.json'),dependencies=identity(root/'science/runtime-dependencies-v1.json'),
        governor=next(v['ELF']for v in b['executables']if Path(v['wrapper']['path']).name=='IncNavierStokesSolver'),
        Python_sources=identity(root/'source-manifest.json'),original_capture_execution=identity(PARENT/'science/flow-vector-execution.json')))
    m2=retained['M']['2']['anchor'];assert m2['checkpoint']['sha256']=='4d0a8697523a2937526df541cd46b07ff956cbf82d6cbd4cc3a4246ac446b564'
    assert m2['completion']['sha256']=='84bc540a1f57c135674ff2e38143cafc35c93b0cb36e7be9bc58d2fe63b9eb19'
    assert m2['callback']['sha256']=='6fd23abd8605841109b2d11e1bafa61be2c0b144cab64741514fed574cb8c908'
    xml_case(verify(selected['cp0']['M']['xml']),root/'science/M2',verify(m2['checkpoint']).parent,'M',0,'Off')
    xml_case(verify(selected['cp0']['F']['xml']),root/'science/F',verify(selected['cp0']['F']['checkpoint']).parent,'F',4,'MappedIncrementV1')
    atomic(root/'science/runtime-reuse-preflight.json',dict(status='PASS',native_source=b['native_source_revision'],new_builds=0,
        ELFs=[v['ELF']for v in b['executables']],objects=identity(root/'science/compiled-source-objects.json'),cache=identity(root/'science/target-cache-compatibility.json')))

def preflight(root):
    selected=load(root/'inputs/selected.json');offsets=load(verify(selected['sparse_structure']))['offsets'];proof={}
    for z in ('C','M'):
        proof[z]=prefix(root,None,z,END[z],selected,offsets,selected['retained_attempts'][z])
    chain=[];previous=checkpoint_payload(verify(selected['cp0']['F']['checkpoint']).parent,include_reference=True)
    for n in range(1,5):
        r,cp=read_completed(selected['oracles']['F'][str(n)],'F',n)
        assert r['pimple_attempts'][0]['incoming_sha256']==vector_hash(previous['state'])
        assert r['native_full_audit_status']==('PASS'if n==4 else'NOT_SCHEDULED')
        chain.append(dict(step=n,anchor=selected['oracles']['F'][str(n)],state=vector_hash(r['state']),history=vector_hash(r['last_stage_history']),native_full_audit=r['native_full_audit_status']))
        previous=cp
    old,_=read_completed(selected['oracles']['F']['4'],'F',4);native_full_map_parity(old,old)
    m,cp=read_completed(selected['retained_attempts']['M']['2']['anchor'],'M',2)
    assert vector_hash(m['state'])=='6e0daed74de86008984548e9946f9d97bc4d707198f9ecab2a8b315f79224885'
    adversarial=[]
    for key,value in (('last_stage_alpha',0.),('time',0.),('dt',.00125),('signature','wrong'),('last_stage_history',[0.]),
                      ('native_history',[[]]),('fixed_reference',[0.]),('fixed_defect',[1.]),('state',[0.])):
        bad=copy.deepcopy(m);bad[key]=value
        try:parity(bad,m,cp,cp)
        except (AssertionError,ValueError):adversarial.append(key)
        else:raise AssertionError('retained producer parity accepted '+key)
    atomic(root/'science/retained-preflight.json',dict(status='PASS',retained_prefixes=proof,F_historical_chain=chain,
        actual_retained_negative_checks=adversarial,native_operations=0,startup_verification='NOT_YET_RUN'))

def freeze(root):
    assert load(root/'science/retained-preflight.json')['status']=='PASS'
    assert load(root/'operations/focused-completion-tests-v1/receipt.json')['status']=='PASS'
    files=load(root/'source-manifest.json')['files']
    files+= [identity(root/'inputs'/n)for n in ('selected.json','parent-index.json','original-immutability.json')]
    files+=[identity(root/'science'/n)for n in ('runtime-build-v1.json','runtime-dependencies-v1.json','compiled-source-objects.json',
        'target-cache-compatibility.json','flow-vector-execution.json','retained-preflight.json','runtime-reuse-preflight.json')]
    files+=[identity(root/f'science/{z}/{n}')for z in ('M2','F')for n in ('case.xml','execution-diff.json')]
    atomic(root/'execution-config.json',dict(schema='WJH22CompletionFrozenV1',files=files,native_work_before_freeze=0,
        matrix=dict(M2_verification=0,F_advances=4,strict_reads=['M2','F4'],new_CP0=0,new_build=0)))

def intent_for(root,role):
    selected=load(root/'inputs/selected.json');level='M'if role=='M2'else'F';folder=root/'science'/role;xml=folder/'case.xml'
    v=dict(schema='WJMonitorIntentV1',policy='CheckpointAuditV1',completion_role='RetainedM2StartupV1'if role=='M2'else'MissingFReplayV1',
        level=level,directory=str(folder/'monitor'),nonce=uuid.uuid4().hex,
        scope_sha256=hashlib.sha256((identity(xml)['sha256']+identity(root/'execution-config.json')['sha256']).encode()).hexdigest(),
        source_sha256=identity(root/'science/runtime-build-v1.json')['sha256'],observer_sha256=identity(HERE/'wj_h22_completion_observer.py')['sha256'],
        inventory_rule_sha256=hashlib.sha256(RULE.encode()).hexdigest(),original_qh0_checkpoint=selected['cp0'][level]['checkpoint'],
        original_qh0_state_sha256=selected['original']['coupled_state_sha256'],spatial_lineage_anchor=selected['original']['checkpoint'],
        maximum_absolute_step=2 if role=='M2'else 4,restored_absolute_step=2 if role=='M2'else 0,recovery_stride=1,
        scientific_stride=2 if role=='M2'else 4,dt_bits=bits(DT[level]),ack_timeout_seconds=120.,strict_reserve_seconds=120,
        substep_evidence_steps=[]if role=='M2'else[1,2,3,4],flow_capture='Off'if role=='M2'else'MappedIncrementV1',
        preserve_qualification_artifacts=True,selected_inputs=identity(root/'inputs/selected.json'),common_support=selected['common_support'])
    if role=='M2':v['restart_completion']=selected['retained_attempts']['M']['2']['anchor']['completion']
    validate_role(v,0 if role=='M2'else 4);return v

def coordinate(root,role):
    frozen(root);folder=root/'science'/role;directory=folder/'monitor';directory.mkdir();ip=folder/'intent.json'
    intent=intent_for(root,role);obs=native=None;code=None;error=None;start=time.monotonic()
    with(folder/'observer.log').open('x')as ol,(folder/'native.log').open('x')as nl:
        try:
            obs=subprocess.Popen([sys.executable,'-B',str(HERE/'wj_h22_completion_observer.py'),'--root',str(root),'--intent',str(ip)],stdout=ol,stderr=subprocess.STDOUT)
            intent['observer_pid']=obs.pid;atomic(ip,intent)
            while not(directory/'observer-started.json').exists():assert obs.poll()is None and time.monotonic()-start<30;time.sleep(.025)
            native=subprocess.Popen([str(executable(root,'IncNavierStokesSolver')),str(folder/'case.xml')],stdout=nl,stderr=subprocess.STDOUT)
            atomic(folder/'native-launch.json',dict(pid=native.pid,xml=identity(folder/'case.xml'),intent=identity(ip),maximum_new_attempts=0 if role=='M2'else 4))
            while native.poll()is None:
                assert obs.poll()is None,'observer exited before native final budgets/residuals';time.sleep(.025)
            code=native.returncode
        except Exception as e:error=repr(e)
        finally:
            if native and native.poll()is None:native.terminate();native.wait(timeout=5)
            atomic(directory/'observer-stop.json',dict(reason='coordinator observed native exit/stop',native_exit=code))
            if obs:
                try:obs.wait(timeout=5)
                except subprocess.TimeoutExpired:obs.terminate();obs.wait(timeout=5)
    report=folder/'case.evolution.json';r=native_load(report)if report.exists()else None
    atomic(folder/'outcome.json',dict(native_exit=code,error=error,producer=identity(report)if r else None,seconds=time.monotonic()-start))
    if r and os.environ.get('WJ_H22C_OPERATION'):atomic(Path(os.environ['WJ_H22C_OPERATION'])/'attempt-counts.json',dict(actual_attempts=r['implicit_callbacks']))
    assert error is None,error
    assert code==0,r['failure']if r else'original native failure retained'
    assert r['initialization_work']==[0,0]and r['initializer_calls']==0
    if role=='M2':
        assert r['steps']==2 and r['time']==.005 and r['implicit_callbacks']==r['native_completed']==r['aborted_advances']==0
        assert r['state_sha256']=='6e0daed74de86008984548e9946f9d97bc4d707198f9ecab2a8b315f79224885'
        assert r['last_stage_alpha']==400 and r['full_audit_calls']==r['physical_samples_calls']==1
        assert not list(directory.glob('step-*.completion.json'))and not list((folder/'snapshots').glob('*'))
        commit=load(directory/'startup.verification-commit.json');assert commit['schema']=='WJStartupVerificationCommitV1'
        assert commit['completion']==intent['restart_completion']and int(commit['binding']['absolute_step'])==2
        atomic(folder/'verification.json',dict(status='PASS',new_advances=0,completion=commit['completion'],new_startup_commit=identity(directory/'startup.verification-commit.json'),
            producer=identity(report),restored_parity=identity(directory/'startup-replay-parity.json'),old_execution_qualification='NOT_COMMITTED_UNCHANGED'))
    else:
        assert r['steps']==r['implicit_callbacks']==r['native_completed']==4 and r['time']==.005 and r['independently_qualified']==1
        assert r['full_audit_calls']==r['physical_samples_calls']==2
        assert load(directory/'prefix-replay-capture.json')['status']=='PASS'
        assert load(directory/'latest-fully-verified-checkpoint.json')['schema']=='WJQualificationCommitV1'
        atomic(folder/'qualification.json',dict(status='PASS',attempts=4,native_completed=4,scientific_qualified_endpoints=1,
            producer=identity(report),prefix=identity(directory/'prefix-replay-capture.json'),new_scientific_steps=0))

def strict(root,role):
    frozen(root);selected=load(root/'inputs/selected.json');folder=root/'science'/role;directory=folder/'monitor'
    anchor=selected['retained_attempts']['M']['2']['anchor']if role=='M2'else recovery(directory,4)
    request=typed_strict(directory,anchor,'M'if role=='M2'else'F',2 if role=='M2'else 4)
    out=folder/'final-strict';out.mkdir();atomic(out/'request.json',request);start=time.monotonic()
    with(out/'native.log').open('x')as log:code=subprocess.run([str(executable(root,'RANSWallResolutionProbe')),str(folder/'case.xml'),str(out/'request.json'),str(out/'native.json')],stdout=log,stderr=subprocess.STDOUT).returncode
    assert code==0;r=native_load(out/'native.json')
    assert r['initialization_work']==[0,0]and r['initial_velocity_calls']==r['governing_solves']==0
    assert r['selected_absolute_step']==(2 if role=='M2'else 4)
    alpha=r['weak_checks']['alpha_mass'];assert alpha==(400 if role=='M2'else 800)
    expected,_=read_completed(anchor,'M'if role=='M2'else'F',2 if role=='M2'else 4)
    assert r['callback_history_sha256']==vector_hash(expected['last_stage_history'])
    moments={str(b):r['spatial_residual'][load(verify(selected['sparse_structure']))['offsets'][b]:load(verify(selected['sparse_structure']))['offsets'][b+1]]for b in (4,5)}
    assert max(abs(x)for v in moments.values()for x in v)<=2e-15
    atomic(out/'result.json',dict(status='PASS',selected_absolute_step=r['selected_absolute_step'],genuine_alpha=alpha,native_signed_Q1=moments,
        producer=identity(out/'native.json'),request=identity(out/'request.json'),physical_mass_dual_norms=r['physical_mass_dual_norms'],
        weak=r['weak_checks'],interface=r['open_boundary'],seconds=time.monotonic()-start))

def composite(root):
    selected=load(root/'inputs/selected.json');attempts=copy.deepcopy(selected['retained_attempts'])
    assert load(root/'science/M2/verification.json')['status']=='PASS'and load(root/'science/F/qualification.json')['status']=='PASS'
    attempts['M']['2']['qualification'].update(subsequent_startup_verification='PASS',
        startup_commit=identity(root/'science/M2/monitor/startup.verification-commit.json'),strict_read=identity(root/'science/M2/final-strict/result.json'))
    attempts['F']={}
    for n in range(1,5):
        q=dict(original_execution_qualification='COMMITTED'if n==4 else'NOT_SCHEDULED',subsequent_startup_verification='NOT_APPLICABLE')
        if n==4:q['original_commit']=identity(root/'science/F/monitor/step-4.qualification-commit.json')
        attempts['F'][str(n)]=capture_descriptor(recovery(root/'science/F/monitor',n),identity(root/'science/flow-vector-execution.json'),q)
    dest=root/'science/composite-source-index.json'
    atomic(dest,dict(schema='WJH22CompositeSourcesV1',attempts=attempts,auxiliary_audits=selected['auxiliary_audits'],
        parent_accounting=load(root/'inputs/parent-index.json')['accounting'],M2_verification_is_not_an_attempt=True))
    from wj_h22_phase import run
    run(root,source_index=identity(dest))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    p.add_argument('--mode',choices=('setup','preflight','freeze','verify_M2','strict_M2','F','strict_F','analysis','runtime'),required=True);a=p.parse_args()
    assert a.root.name==ROOT
    if a.mode=='verify_M2':coordinate(a.root,'M2')
    elif a.mode=='F':coordinate(a.root,'F')
    elif a.mode.startswith('strict_'):strict(a.root,'M2'if a.mode=='strict_M2'else'F')
    elif a.mode=='analysis':composite(a.root)
    elif a.mode=='runtime':runtime_qualify(a.root)
    else:globals()[a.mode](a.root)
