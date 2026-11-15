#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Seven exact replay attempts, three own-CP0 restores, no initialization route."""
import argparse,hashlib,os,subprocess,sys,time,uuid
from pathlib import Path
import xml.etree.ElementTree as ET
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,verify,atomic,vector_hash,bits,RULE
from wj_h17_observer import compatibility,native_load,equal
from wj_h22_replay import DT,END,read_completed
from study_wj_h14 import runtime as prepare_runtime
from study_wj_h16 import same_cp
from study_wj_h17 import cp_request,runtime_qualify
from study_wj_h18 import protected
from wj3_measure import checkpoint_payload
from cls_epawf_m3_cases import set_info,set_parameter,write
from run_wj_h22 import frozen,usage
HERE=Path(__file__).resolve().parent
EXTRA={'RANSWallPimpleFlowVectors','RANSWallPimpleFlowVectorIdentity'}

def protected_case(x):
    import copy
    x=copy.deepcopy(x)
    for node in list(x.find('./CONDITIONS/SOLVERINFO')):
        if node.get('PROPERTY')in EXTRA:x.find('./CONDITIONS/SOLVERINFO').remove(node)
    return protected(x)

def validate_intent(v):
    z=v['level'];assert v['schema']=='WJMonitorIntentV1'and 0<v['ack_timeout_seconds']<=120
    assert v['maximum_absolute_step']==END[z]and v['restored_absolute_step']==0
    assert v['recovery_stride']==1 and v['scientific_stride']==END[z]and v['dt_bits']==bits(DT[z])
    assert v['substep_evidence_steps']==list(range(1,END[z]+1))

def executable(root,name):
    return verify(next(v['wrapper']for v in load(root/'science/runtime-build-v1.json')['executables']if Path(v['wrapper']['path']).name==name))

def xml_case(original,folder,cp,z):
    folder.mkdir(parents=True);x=ET.parse(original).getroot()
    for k,v in (('NumSteps',END[z]),('IO_CheckSteps',1),('IO_InfoSteps',1)):set_parameter(x,k,v)
    for k,v in (('RANSWallRestart',str(cp)),('RANSWallCheckpoint',str(folder/'CP')),('RANSWallRetention','BoundedV1'),
        ('RANSWallMonitorIntent',str(folder/'intent.json')),('RANSWallMonitorFault','Off'),('RANSWallReport',str(folder/'native.json')),
        ('RANSWallImmutableRowReuse','Off'),('RANSWallMonitoringPolicy','CheckpointAuditV1'),('RANSWallPimpleSnapshots','Off'),
        ('RANSWallPimpleArithmetic','Off'),('RANSWallPimpleSnapshotDirectory',str(folder/'snapshots')),
        ('RANSWallPimpleFlowVectors','MappedIncrementV1'),('RANSWallPimpleFlowVectorIdentity',str(folder.parent/'flow-vector-execution.json'))):set_info(x,k,v)
    assert protected_case(x)==protected_case(ET.parse(original).getroot()),'scientific XML/TimeStep changed'
    params={p.text.split('=')[0].strip():p.text.split('=')[1].strip()for p in x.findall('./CONDITIONS/PARAMETERS/P')}
    assert float(params['Time'])==0 and float(params['TimeStep'])==DT[z]
    path=folder/'case.xml';write(x,path)
    atomic(folder/'execution-diff.json',dict(original=identity(original),actual=identity(path),TimeStep='UNCHANGED',Time='UNCHANGED ZERO',
         permitted_additional_solverinfo=sorted(EXTRA),additional_NumSteps=END[z],protected_sha256=hashlib.sha256(protected_case(x)).hexdigest(),
         physical_signature_definition_changed=False))
    return path

def setup(root):
    inputs=load(root/'inputs/selected.json');old=load(verify(inputs['runtime_build']))
    assert old['native_source_revision']=='b2365e52b0d334d7b61c6be3a4a3bb9029491084'
    for v in old['executables']:verify(v['ELF']);verify(v['wrapper'])
    for v in load(verify(inputs['compiled']))['dependencies']:verify(v)
    prepare_runtime(root);b=load(root/'science/runtime-build.json')
    b.update(package='H22',native_source_revision=load(root/'source-intent.json')['native_revision'],
             prior_native_identity=inputs['runtime_build'],compile_database_scope='registered compile commands; actual built/reused objects recorded separately')
    atomic(root/'science/runtime-build-v1.json',b)
    deps=load(verify(inputs['runtime_dependencies']));libs=[]
    oldroot=Path(inputs['runtime_build']['path']).parents[1]
    for rec in deps['libraries']:
        p=verify(rec['actual']);new=root/'runtime'/p.relative_to(oldroot/'runtime')if p.is_relative_to(oldroot/'runtime')else p
        v=identity(new);assert v['sha256']==rec['actual']['sha256'];libs.append(dict(actual=v,historical=rec['actual']))
    atomic(root/'science/runtime-dependencies-v1.json',dict(status='PASS',libraries=libs,historical=inputs['runtime_dependencies']))
    from publish_wj_h15 import objects
    atomic(root/'science/compiled-source-objects.json',objects(root,root/'build',b['native_source_revision']))
    previous=load(verify(inputs['cache_manifest']))
    compatibility(root,verify(previous['historical_manifest']),root/'science/cache-base-compatibility.json')
    c=load(root/'science/cache-base-compatibility.json')
    c.update(schema='WJH22ImmutableCacheCompatibilityV1',permitted_difference='native passive vector writer only; dt-compatible existing CP0 per branch',
        dt_independence=dict(status='PASS',timesteps=DT,cached_objects='unchanged geometry/basis/signs/tests/quadrature/analytical profile/spatial qh0/reference/defect',
            exclusions='all actual fields, histories, alpha, closure/root, current-K cutoff, reactions and decisions recomputed',historical_manifest_edited=False))
    atomic(root/'science/target-cache-compatibility.json',c)
    atomic(root/'science/flow-vector-execution.json',dict(schema='WJH22FlowCaptureExecutionV1',
        native_source_revision=b['native_source_revision'],compiled=identity(root/'science/compiled-source-objects.json'),
        runtime_build=identity(root/'science/runtime-build-v1.json'),dependencies=identity(root/'science/runtime-dependencies-v1.json'),
        governor=next(v['ELF']for v in b['executables']if Path(v['wrapper']['path']).name=='IncNavierStokesSolver'),
        vector_writer=identity(root/'source/solvers/IncNavierStokesSolver/EquationSystems/RANSWallFlowVectorCapture.h')))
    base=None;initial=[]
    for z in DT:
        anchor=inputs['cp0'][z];cp=checkpoint_payload(verify(anchor['checkpoint']).parent,include_reference=True)
        assert cp['steps']==0 and bits(cp['time'])==bits(0.)and cp['dt']==DT[z]
        producer=native_load(verify(anchor['report']));same_cp(producer,cp)
        assert vector_hash(cp['state'])=='631caa32a59ff43253f1fcac4effd955eb68b8529d9c28f29e615cf14ca0f077'
        if base is None:base=cp
        for k in ('state','fixed_reference','fixed_defect','eta'):equal(cp[k],base[k],'same original spatial '+k)
        initial.append(dict(level=z,anchor=anchor,dt=cp['dt'],native_history_hashes=[vector_hash(v)for v in cp['native_history']],native_history_times=cp['native_history_times'],
                            state_sha256=vector_hash(cp['state']),signature_sha256=hashlib.sha256(cp['signature'].encode()).hexdigest(),native_history_productions=0))
        for n,a in inputs['oracles'][z].items():read_completed(a,z,int(n))
    atomic(root/'science/initial-and-oracle-lineage.json',dict(status='PASS',initial=initial,completed_oracles=7,
          spatial_initializations=0,native_history_productions=0,CP0_productions=0))

def component(root):
    out=root/'science/vector-component';out.mkdir()
    rt=root/'runtime';command=[str(rt/'glibc/ld-linux-x86-64.so.2'),'--library-path',str(rt/'glibc')+':'+str(rt/'geometry')+':/home/jay/.local/opt/nektar-v2026-mpi/lib64',
          str(root/'build/WJFlowVectorChecks'),str(out)]
    with(out/'native.log').open('x')as log:code=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT).returncode
    assert code==0
    packet=native_load(out/'flow.json');payload=verify(packet['payload']).read_bytes()
    assert len(payload)==56 and packet['arrays'][0]['count']==2 and packet['velocity_rows']==[0,1]and packet['pressure_rows']==[4]
    atomic(out/'result.json',dict(status='PASS',command=command,ELF=identity(root/'build/WJFlowVectorChecks'),packet=identity(out/'flow.json'),
        payload=packet['payload'],governing_attempts=0,CP0s=0,history_productions=0))

def freeze(root):
    assert load(root/'science/vector-component/result.json')['status']=='PASS'
    assert load(root/'operations/focused-tests-v1/receipt.json')['status']=='PASS'
    assert load(root/'operations/prepared-session-test-v1/receipt.json')['status']=='PASS'
    files=[identity(HERE/n)for n in ('run_wj_h22.py','study_wj_h22.py','wj_h22_observer.py','wj_h22_replay.py','wj_h22_phase.py','test_wj_h22_phase.py')]
    files.extend(identity(root/'science'/n)for n in ('runtime-build-v1.json','runtime-dependencies-v1.json','compiled-source-objects.json',
         'target-cache-compatibility.json','initial-and-oracle-lineage.json','flow-vector-execution.json','vector-component/result.json'))
    files.extend([identity(root/'inputs/selected.json'),identity(root/'source-intent.json')])
    inputs=load(root/'inputs/selected.json')
    for z in DT:
        cp0=inputs['cp0'][z]
        xml=xml_case(verify(cp0['xml']),root/'science'/z,verify(cp0['checkpoint']).parent,z)
        files.extend((identity(xml),identity(xml.parent/'execution-diff.json')))
    for v in load(root/'science/runtime-build-v1.json')['executables']:files.extend((v['ELF'],v['wrapper']))
    atomic(root/'execution-config.json',dict(schema='WJH22FrozenExecutionV1',files=files,levels=['C','M','F'],DT=DT,attempts=END,
        maximum_blocks=3,maximum_attempts=7,maximum_CP0=0,maximum_native_history_production=0,maximum_strict=3,
        scientific_strides=END,recovery_stride=1,native_arithmetic='Off',capture='selected MappedIncrementV1',
        strict_reserve_seconds=360,closure_reserve_seconds=300,temporal_ratios_are_not_halt_gates=True))

def admission(root,z):
    frozen(root);b=load(root/'accounting.json');assert not any(b.get(k)for k in ('runtime_halt','target_halt','coverage_halt','integrity_halt'))
    order=list(DT);idx=order.index(z)
    for prev in order[:idx]:assert load(root/f'science/{prev}/qualification.json')['status']=='PASS'
    remaining=order[idx:];estimate=sum(110+END[v]*45+45+15+120+15 for v in remaining)
    if idx:
        # Retained C is a measured admission pilot, not a repeatable extra run.
        c=load(root/'science/C/qualification.json')
        actual=c['governing_block_seconds']+c.get('strict_seconds',0)
        estimate=max(estimate,sum(actual+(END[v]-1)*45 for v in remaining))
    u=usage(b);assert u['total']+1.25*estimate<=4500 and u['native']+1.25*estimate<=1800
    atomic(root/f'admission-{z}.json',dict(status='PASS',level=z,remaining_estimate=estimate,headroom=1.25,protected_strict_seconds=120*len(remaining),closure_seconds=300,
        command=[sys.executable,'-B',str(HERE/'study_wj_h22.py'),'--root',str(root),'--mode','target','--level',z]))

def target(root,z):
    frozen(root);inputs=load(root/'inputs/selected.json');cp0=inputs['cp0'][z];folder=root/'science'/z
    xml=folder/'case.xml'  # prepared and sealed for all three levels before C admission
    directory=folder/'monitor';directory.mkdir();ip=folder/'intent.json'
    intent=dict(schema='WJMonitorIntentV1',policy='CheckpointAuditV1',level=z,directory=str(directory),nonce=uuid.uuid4().hex,
        scope_sha256=hashlib.sha256((identity(xml)['sha256']+identity(root/'execution-config.json')['sha256']).encode()).hexdigest(),
        source_sha256=identity(root/'science/runtime-build-v1.json')['sha256'],observer_sha256=identity(HERE/'wj_h22_observer.py')['sha256'],inventory_rule_sha256=hashlib.sha256(RULE.encode()).hexdigest(),
        original_qh0_checkpoint=cp0['checkpoint'],original_qh0_state_sha256=inputs['original']['coupled_state_sha256'],spatial_lineage_anchor=inputs['original']['checkpoint'],
        maximum_absolute_step=END[z],restored_absolute_step=0,recovery_stride=1,scientific_stride=END[z],dt_bits=bits(DT[z]),ack_timeout_seconds=120.,
        strict_reserve_seconds=120,substep_evidence_steps=list(range(1,END[z]+1)),preserve_qualification_artifacts=True,
        selected_inputs=identity(root/'inputs/selected.json'),common_support=inputs['common_support'])
    validate_intent(intent)
    obs=native=None;code=None;error=None;start=time.monotonic()
    with(folder/'observer.log').open('x')as ol,(folder/'native.log').open('x')as nl:
        try:
            obs=subprocess.Popen([sys.executable,'-B',str(HERE/'wj_h22_observer.py'),'--root',str(root),'--intent',str(ip),'--compat',str(root/'science/target-cache-compatibility.json')],stdout=ol,stderr=subprocess.STDOUT)
            intent['observer_pid']=obs.pid;atomic(ip,intent)
            while not(directory/'observer-started.json').exists():assert obs.poll()is None and time.monotonic()-start<30;time.sleep(.025)
            native=subprocess.Popen([str(executable(root,'IncNavierStokesSolver')),str(xml)],stdout=nl,stderr=subprocess.STDOUT)
            atomic(folder/'native-launch.json',dict(pid=native.pid,xml=identity(xml),intent=identity(ip),maximum_attempts=END[z]))
            while native.poll()is None:
                assert obs.poll()is None,'observer exit propagated to coordinator';time.sleep(.025)
            code=native.returncode
        except Exception as e:error=repr(e)
        finally:
            if native and native.poll()is None:native.terminate();native.wait(timeout=5)
            atomic(directory/'observer-stop.json',dict(reason='H22 coordinator exit'))
            if obs:
                try:obs.wait(timeout=5)
                except subprocess.TimeoutExpired:obs.terminate();obs.wait(timeout=5)
    report=folder/'case.evolution.json';r=native_load(report)if report.exists()else None
    atomic(folder/'outcome.json',dict(native_exit=code,error=error,producer=identity(report)if r else None,seconds=time.monotonic()-start))
    if r:
        op=os.environ.get('WJ_H22_OPERATION')
        if op:atomic(Path(op)/'attempt-counts.json',dict(actual_attempts=r['implicit_callbacks'],producer=identity(report)))
    assert error is None,error
    assert code==0, r['failure']if r else'native nonzero; original native.log retained'
    assert r['initialization_work']==[0,0]and r['initializer_calls']==0
    assert r['steps']==END[z]and r['implicit_callbacks']==r['native_completed']==END[z]
    assert int(load(directory/'latest-fully-verified-checkpoint.json')['binding']['absolute_step'])==END[z]
    assert load(directory/'prefix-replay-capture.json')['status']=='PASS'
    assert r['full_audit_calls']==r['physical_samples_calls']==2
    atomic(folder/'qualification.json',dict(status='PASS',attempts=END[z],native_completed=r['native_completed'],
        scientific_qualified_steps=r['independently_qualified'],final_absolute_step=r['steps'],time=r['time'],governing_block_seconds=time.monotonic()-start,
        producer=identity(report),prefix_parity=identity(directory/'prefix-replay-capture.json'),new_scientific_steps=0))

def strict(root,z):
    frozen(root);folder=root/'science'/z;directory=folder/'monitor';cursor=directory/'latest-fully-verified-checkpoint.json'
    if cursor.exists()and load(cursor)['schema']=='WJQualificationCommitV1':
        c=load(cursor);step=int(c['binding']['absolute_step']);commit=directory/f'step-{step}.qualification-commit.json'
        request=dict(operation='open_checks',inner_stress_model='Linear',native_report=c['report']['path'],native_report_sha256=c['report']['sha256'],
             native_checkpoint=str(verify(c['checkpoint']).parent),checkpoint_sha256=c['checkpoint']['sha256'],completed_sidecar=c['callback'],qualified_cursor=identity(commit))
        xml=folder/'case.xml'
    else:
        cp0=load(root/'inputs/selected.json')['cp0'][z]
        request=dict(operation='open_checks',inner_stress_model='Linear',**cp_request(dict(producer=cp0['report'],checkpoint=cp0['checkpoint'])));xml=verify(cp0['xml'])
    request.update(weak=True,weak_compact=True,assembly=False,reuse=False,retain_scalar_observations=False)
    out=folder/'final-strict';out.mkdir();atomic(out/'request.json',request);start=time.monotonic()
    with(out/'native.log').open('x')as log:code=subprocess.run([str(executable(root,'RANSWallResolutionProbe')),str(xml),str(out/'request.json'),str(out/'native.json')],stdout=log,stderr=subprocess.STDOUT).returncode
    assert code==0;r=native_load(out/'native.json');assert r['initialization_work']==[0,0]and r['initial_velocity_calls']==r['governing_solves']==0
    atomic(out/'result.json',dict(status='PASS',selected_absolute_step=r['selected_absolute_step'],producer=identity(out/'native.json'),
        physical_mass_dual_norms=r['physical_mass_dual_norms'],weak=r['weak_checks'],interface=r['open_boundary'],seconds=time.monotonic()-start))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    p.add_argument('--mode',choices=('setup','component','freeze','admission','target','strict','analysis','runtime'),required=True);p.add_argument('--level',choices=tuple(DT));a=p.parse_args()
    if a.mode in ('admission','target','strict'):assert a.level;globals()[a.mode](a.root,a.level)
    elif a.mode=='analysis':
        from wj_h22_phase import run
        run(a.root)
    elif a.mode=='runtime':runtime_qualify(a.root)
    else:globals()[a.mode](a.root)
