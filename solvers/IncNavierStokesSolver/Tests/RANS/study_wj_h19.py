#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""One step-zero native bootstrap and one bounded physical-time comparison."""
import argparse,hashlib,os,subprocess,sys,time,uuid
from pathlib import Path
import xml.etree.ElementTree as ET
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,verify,atomic,vector_hash,bits,RULE,H14
from wj_h17_observer import compatibility,load_compatible,native_load,equal
from wj_h19_analysis import spatial_xml
from wj_h19_observer import assess_pair
from study_wj_h14 import runtime as prepare_runtime
from study_wj_h16 import same_cp
from study_wj_h17 import cp_request,runtime_qualify
from wj3_measure import checkpoint_payload
from cls_epawf_m3_cases import set_info,set_parameter,write
from run_wj_h19 import frozen,halt,reserve_strict
HERE=Path(__file__).resolve().parent
OLD=Path('/home/jay/build/wj-h18-checkpoint-audits-20261114')

def xml_case(original,folder,cp,steps):
    assert steps in (0,2);folder.mkdir(parents=True)
    x=ET.parse(original).getroot()
    for k,v in (('TimeStep',.005),('NumSteps',steps),('IO_CheckSteps',1),('IO_InfoSteps',1)):set_parameter(x,k,v)
    for k,v in (('RANSWallRestart',str(cp)),('RANSWallCheckpoint',str(folder/'CP')),('RANSWallRetention','BoundedV1'if steps else'Off'),
        ('RANSWallMonitorIntent',str(folder/'intent.json')),('RANSWallMonitorFault','Off'),('RANSWallReport',str(folder/'native.json')),
        ('RANSWallImmutableRowReuse','Off'),('RANSWallMonitoringPolicy','CheckpointAuditV1'if steps else'Off'),
        ('RANSWallPimpleSnapshots','Off'),('RANSWallPimpleArithmetic','Off'),('RANSWallPimpleSnapshotDirectory',str(folder/'snapshots'))):set_info(x,k,v)
    assert spatial_xml(x)==spatial_xml(ET.parse(original).getroot()),'only dt and declared execution controls may differ'
    p=folder/'case.xml';write(x,p)
    atomic(folder/'execution-diff.json',dict(original=identity(original),actual=identity(p),physical_except_dt_sha256=hashlib.sha256(spatial_xml(x)).hexdigest(),
        TimeStep=dict(original=.0025 if steps else 'source own dt',actual=.005),additional_steps=steps,Time='unchanged zero',signature_exclusions_changed=False))
    return p

def setup(root):
    prepare_runtime(root);b=load(root/'science/runtime-build.json');b.update(package='H19',native_source_revision=load(root/'source-intent.json')['native_revision'])
    atomic(root/'science/runtime-build-v1.json',b)
    old=load(OLD/'science/runtime-dependencies-v1.json');libs=[]
    for rec in old['libraries']:
        pin=rec['actual'];verify(pin);p=Path(pin['path']);new=root/'runtime'/p.relative_to(OLD/'runtime')if p.is_relative_to(OLD/'runtime')else p
        v=identity(new);assert v['sha256']==pin['sha256'];libs.append(dict(actual=v,historical=pin))
    atomic(root/'science/runtime-dependencies-v1.json',dict(status='PASS',libraries=libs,historical=identity(OLD/'science/runtime-dependencies-v1.json')))
    compatibility(root,Path('/home/jay/build/wj-h15-persistent-neutral-20261114/science/observer-cache/manifest.json'),root/'science/cache-base-compatibility.json')
    c=load(root/'science/cache-base-compatibility.json')
    c.update(schema='WJH19DtIndependentCacheCompatibilityV1',permitted_difference='dt .0025 -> .005; immutable geometry/basis/test/quadrature and analytical profile/qh0/reference/defect only; no history/state-result reuse',
        dt_independence=dict(status='PASS',fine_dt=.0025,coarse_dt=.005,
            cached_objects='native row IDs/signs/basis/geometric Jacobians, projected tests, physical quadrature, exact analytical profile and original qh0 evaluations',
            exclusions='cached report signature retained as historical metadata; no signature checks suppressed; no derivative/history/alpha residuals or current closure/profile values reused',
            spatial_identity_required_before_target=True))
    atomic(root/'science/target-cache-compatibility.json',c)
    from publish_wj_h15 import objects
    atomic(root/'science/compiled-source-objects.json',objects(root,root/'build',b['native_source_revision']))
    selected=load(root/'inputs/selected.json')
    # Bind selected heavy bytes back to the published Mac/Git evidence, not
    # freshly generated hashes alone. Only required files are read.
    h14table=H14/'publication-v2/state-bindings.json'
    assert identity(h14table)['sha256']=='c0f594a286e98db50fd9cddc89de95262ed62c15378e8cb6dab2405a623ca4aa'
    heavy=OLD/'publication-v1/heavy-manifest.json'
    assert identity(heavy)['sha256']=='7354fda2dbf51c02b6a4232207f3d86b5866820ce308b0b4afab8a9d3a12b9d1'
    known={v['path']:v for v in load(heavy)['files']}
    for anchor in selected['fine'].values():
        for key in ('checkpoint','callback','report'):
            assert anchor[key]==known[anchor[key]['path']];verify(anchor[key])
    witness_proof=OLD/'publication-supplement-v1/pending-restored-native-map-parity-v1.json'
    assert identity(witness_proof)['sha256']=='a21a972dc7eab313d11aa1f5fa18759ad9461088c4d2ccc362e11f9a5df1dd95'
    verify(load(witness_proof)['current'])
    for n in (2,4):
        ready=load(OLD/'science/target/B/monitor'/('startup.ready.json'if n==2 else'step-4.ready.json'))
        selected['fine'][str(n)]['audit_report']=ready['report'];verify(ready['report'])
        if n==2:assert ready['report']==load(witness_proof)['current']
    atomic(root/'inputs/anchors-v2.json',selected)
    files=[identity(HERE/n)for n in ('run_wj_h19.py','study_wj_h19.py','wj_h19_analysis.py','wj_h19_observer.py','test_wj_h19_timestep.py')]
    for n in ('runtime-build-v1.json','runtime-dependencies-v1.json','target-cache-compatibility.json','compiled-source-objects.json'):files.append(identity(root/'science'/n))
    for v in b['executables']:files.extend((v['ELF'],v['wrapper']))
    files.extend(identity(root/'inputs'/n)for n in ('selected.json','anchors-v2.json'))
    atomic(root/'execution-config.json',dict(schema='WJH19FrozenExecutionV1',files=files,maximum_bootstrap=1,maximum_target_attempts=2,maximum_target_blocks=1,maximum_strict=1,
        physical='unchanged original 6mm target',monitoring='CheckpointAuditV1 recovery/scientific stride1',comparison_role='H19PhysicalTimeSensitivityV1'))

def component(root):
    original=H14/'science/small/cp0/case.xml';report=H14/'science/small/cp0/case.evolution.json';r=load(report)
    folder=root/'science/component';xml=xml_case(original,folder,'None',0)
    packet=dict(operation='step_zero_bootstrap_checks',inner_stress_model='Linear',source_dt=r['dt'],source_xml=identity(original),
        native_report=str(report),native_report_sha256=identity(report)['sha256'],native_checkpoint=r['checkpoint'],checkpoint_sha256=identity(Path(r['checkpoint'])/'State.bin')['sha256'])
    atomic(folder/'request.json',packet)
    with(folder/'native.log').open('x')as log:code=subprocess.run([str(root/'runtime/RANSWallResolutionProbe'),str(xml),str(folder/'request.json'),str(folder/'result.json')],stdout=log,stderr=subprocess.STDOUT).returncode
    assert code==0;v=native_load(folder/'result.json');assert v['status']=='PASS'and v['fixture']and v['native_history_initializations']==1 and v['governing_solves']==0
    assert v['state_sha256']==vector_hash(r['state']);assert v['dt']==.005 and v['source_dt']==r['dt']
    atomic(folder/'qualification.json',dict(status='PASS',result=identity(folder/'result.json'),fresh_native_CP0=0,governing_attempts=0,component_history_initializations=1))

def bootstrap(root):
    frozen(root);selected=load(root/'inputs/anchors-v2.json');v=selected['original'];folder=root/'science/coarse-cp0'
    xml=xml_case(verify(v['xml']),folder,'None',0)
    request=dict(operation='step_zero_bootstrap',inner_stress_model='Linear',source_dt=.0025,source_xml=v['xml'],checkpoint_base=str(folder/'CP'),**cp_request(v))
    request['operation']='step_zero_bootstrap';atomic(folder/'request.json',request)
    with(folder/'native.log').open('x')as log:code=subprocess.run([str(root/'runtime/RANSWallResolutionProbe'),str(xml),str(folder/'request.json'),str(folder/'native.json')],stdout=log,stderr=subprocess.STDOUT).returncode
    assert code==0,'coarse native CP0 reservation consumed; no retry'
    r=native_load(folder/'native.json');source=checkpoint_payload(verify(v['checkpoint']).parent,include_reference=True);coarse=checkpoint_payload(r['checkpoint'],include_reference=True)
    for key in ('state','fixed_reference','fixed_defect','eta'):equal(coarse[key],source[key],'step-zero physical '+key)
    assert source['steps']==coarse['steps']==0 and bits(source['time'])==bits(coarse['time'])==bits(0.)
    assert source['dt']==.0025 and coarse['dt']==.005 and coarse['signature']!=source['signature']
    assert r['native_history_initializations']==1 and r['initialization_work']==[0,0]and r['initial_velocity_calls']==r['governing_solves']==0
    cache,_=load_compatible(root/'science/target-cache-compatibility.json');fixed=load(verify(selected['common_support']))['points'];fine=native_load(verify(v['producer']))
    result,physical,inventory=assess_pair(cache,fixed,r,fine,fine,0,folder,'cp0')
    assert physical['summary']['status']==result['physical_fine_status']==result['status']=='PASS'
    atomic(folder/'lineage.json',dict(status='PASS',source=v,coarse_checkpoint=identity(Path(r['checkpoint'])/'State.bin'),producer=identity(folder/'native.json'),coarse_xml=identity(xml),
        fine_signature=source['signature'],coarse_signature=coarse['signature'],source_dt=source['dt'],coarse_dt=coarse['dt'],physical_hashes={k:vector_hash(coarse[k])for k in ('state','fixed_reference','fixed_defect')},eta=coarse['eta'],
        native_history_initializations=1,spatial_initializations=0,source_history='never copied/rescaled; native IMEX1 BuildScheme',spatial_lineage_anchor=v['checkpoint']))
    atomic(folder/'qualification.json',dict(status='PASS',lineage=identity(folder/'lineage.json'),result=result,inventory=identity(inventory)))

def admission(root):
    frozen(root);assert load(root/'science/component/qualification.json')['status']=='PASS';assert load(root/'science/coarse-cp0/qualification.json')['status']=='PASS'
    b=load(root/'accounting.json');assert not any(b.get(k)for k in ('runtime_halt','qualification_halt','target_halt','coverage_halt','screen_halt'))
    spent=60+300+sum(v.get('seconds',v['reserved_seconds'])for v in b['entries'].values())
    forecast=1.25*(2*110+100+120+60)
    assert spent+forecast<=5400
    native=sum(v.get('seconds',v['reserved_seconds'])for v in b['entries'].values()if v['phase']in ('qualification','target','target_readback'))
    assert native+forecast<=1800
    for name in ('focused-python-v2','native-protocol-components-v1'):
        assert load(root/f'operations/{name}/receipt.json')['status']=='PASS'
    atomic(root/'target-admission.json',dict(status='PASS',forecast_with_25percent_headroom=forecast,closure_reserved=300,strict_reserved=120,
        command=[sys.executable,'-B',str(HERE/'study_wj_h19.py'),'--root',str(root),'--mode','target'],inputs=[identity(root/'execution-config.json'),identity(root/'science/coarse-cp0/qualification.json'),identity(root/'science/component/qualification.json')]))

def target(root):
    frozen(root);selected=load(root/'inputs/anchors-v2.json');source=selected['original'];cp0=load(root/'science/coarse-cp0/lineage.json')
    folder=root/'science/coarse';xml=xml_case(verify(source['xml']),folder,verify(cp0['coarse_checkpoint']).parent,2)
    directory=folder/'monitor';directory.mkdir();ip=folder/'intent.json'
    intent=dict(schema='WJMonitorIntentV1',policy='CheckpointAuditV1',directory=str(directory),nonce=uuid.uuid4().hex,
        scope_sha256=hashlib.sha256((identity(xml)['sha256']+identity(root/'execution-config.json')['sha256']).encode()).hexdigest(),
        source_sha256=identity(root/'science/runtime-build-v1.json')['sha256'],observer_sha256=identity(HERE/'wj_h19_observer.py')['sha256'],inventory_rule_sha256=hashlib.sha256(RULE.encode()).hexdigest(),
        original_qh0_checkpoint=cp0['coarse_checkpoint'],original_qh0_state_sha256=cp0['physical_hashes']['state'],spatial_lineage_anchor=source['checkpoint'],
        maximum_absolute_step=2,restored_absolute_step=0,recovery_stride=1,scientific_stride=1,dt_bits=bits(.005),ack_timeout_seconds=120.,strict_reserve_seconds=120,
        substep_evidence_steps=[],preserve_qualification_artifacts=True,common_support=selected['common_support'],fine_anchors=identity(root/'inputs/anchors-v2.json'))
    obs=native=None;code=None;error=None;start=time.monotonic()
    with(folder/'observer.log').open('x')as ol,(folder/'native.log').open('x')as nl:
        try:
            obs=subprocess.Popen([sys.executable,'-B',str(HERE/'wj_h19_observer.py'),'--root',str(root),'--intent',str(ip),'--compat',str(root/'science/target-cache-compatibility.json')],stdout=ol,stderr=subprocess.STDOUT)
            intent['observer_pid']=obs.pid;atomic(ip,intent)
            while not(directory/'observer-started.json').exists():assert obs.poll()is None and time.monotonic()-start<30;time.sleep(.025)
            native=subprocess.Popen([str(root/'runtime/IncNavierStokesSolver'),str(xml)],stdout=nl,stderr=subprocess.STDOUT)
            while native.poll()is None:
                assert obs.poll()is None,'observer exited; stop target'
                time.sleep(.025)
            code=native.returncode
        except Exception as e:error=repr(e)
        finally:
            if native and native.poll()is None:native.terminate();native.wait(timeout=5)
            atomic(directory/'observer-stop.json',dict(reason='H19 coordinator exit'))
            if obs:
                try:obs.wait(timeout=5)
                except subprocess.TimeoutExpired:obs.terminate();obs.wait(timeout=5)
    report=folder/'case.evolution.json';r=load(report)if report.exists()else None
    atomic(folder/'outcome.json',dict(native_exit=code,error=error,producer=identity(report)if r else None,seconds=time.monotonic()-start))
    if r:
        op=os.environ.get('WJ_H19_OPERATION')
        if op:atomic(Path(op)/'attempt-counts.json',dict(kind='target',actual_attempts=r['implicit_callbacks'],producer=identity(report)))
        assert r['initialization_work']==[0,0]and r['initializer_calls']==0
        for step in range(1,r['native_completed']+1):
            p=load(directory/f'step-{step}.report.json')['pimple_attempts'][0]
            assert all(p[k]==1 for k in ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls'))
            assert p['bound_activations']==p['additional_correction_or_restoration_solves']==0 and p['k_bound']['raw_sha256']==p['k_bound']['corrected_sha256']
    assert error is None,error
    if code!=0:
        assert r and r['failure'].startswith('TEMPORAL_REVIEW_STOP')and load(root/'accounting.json').get('screen_halt'),'unexpected native failure'
    else:assert r and r['steps']==2 and r['native_completed']==r['independently_qualified']==2
    atomic(folder/'qualification.json',dict(native_exit=code,native_completed=r['native_completed'],independently_qualified=r['independently_qualified'],attempts=r['implicit_callbacks'],
        temporal_status='STOPPED'if code else'PASS',status='PASS'if code==0 else'PROVISIONAL_SCREEN_STOP',new_coarse_scientific_attempts=r['implicit_callbacks'],fine_attempts=0))

def strict(root):
    reserve_strict(root,'target');folder=root/'science/coarse';directory=folder/'monitor';p=directory/'latest-fully-verified-checkpoint.json'
    if p.exists()and load(p)['schema']=='WJQualificationCommitV1':
        c=load(p);step=int(c['binding']['absolute_step']);commit=directory/f'step-{step}.qualification-commit.json'
        request=dict(operation='open_checks',inner_stress_model='Linear',native_report=c['report']['path'],native_report_sha256=c['report']['sha256'],native_checkpoint=str(verify(c['checkpoint']).parent),
            checkpoint_sha256=c['checkpoint']['sha256'],completed_sidecar=c['callback'],qualified_cursor=identity(commit))
        xml=folder/'case.xml'
    else:
        cp=root/'science/coarse-cp0';r=native_load(cp/'native.json');request=dict(operation='open_checks',inner_stress_model='Linear',native_report=str(cp/'native.json'),native_report_sha256=identity(cp/'native.json')['sha256'],
            native_checkpoint=r['checkpoint'],checkpoint_sha256=identity(Path(r['checkpoint'])/'State.bin')['sha256']);xml=cp/'case.xml'
    request.update(weak=True,weak_compact=True,assembly=False,reuse=False,retain_scalar_observations=False)
    out=root/'science/final-strict';out.mkdir();atomic(out/'request.json',request)
    with(out/'native.log').open('x')as log:code=subprocess.run([str(root/'runtime/RANSWallResolutionProbe'),str(xml),str(out/'request.json'),str(out/'native.json')],stdout=log,stderr=subprocess.STDOUT).returncode
    assert code==0;r=load(out/'native.json');assert r['initialization_work']==[0,0]and r['initial_velocity_calls']==r['governing_solves']==0
    atomic(out/'result.json',dict(status='PASS',selected_absolute_step=r['selected_absolute_step'],producer=identity(out/'native.json'),physical_mass_dual_norms=r['physical_mass_dual_norms'],weak=r['weak_checks'],interface=r['open_boundary']))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('setup','component','bootstrap','admission','target','strict'),required=True);a=p.parse_args();globals()[a.mode](a.root)
