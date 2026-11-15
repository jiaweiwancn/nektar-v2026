#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Only three unchanged nonzero restarts and their protected qualified reads."""
import argparse,hashlib,json,os,subprocess,sys,time,uuid
from pathlib import Path
import xml.etree.ElementTree as ET
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,verify,atomic,vector_hash,bits,RULE
from wj_h17_observer import native_load,equal,compatibility,load_compatible
from wj_h22_replay import recovery,parity,integer
from wj_h23_continuation import DT,START,END,ADDITIONAL,STRIDE,ORDER,OBSERVER_FORECAST,branch_forecast,qualified,read_completed,validate_intent,strict_request
from wj_h23_observer import evaluate
from study_wj_h22 import executable,protected_case
from study_wj_h16 import same_cp
from study_wj_h17 import runtime_qualify
from cls_epawf_m3_cases import set_parameter,set_info,write
from wj3_measure import checkpoint_payload
from run_wj_h23 import ROOT,frozen,usage

HERE=Path(__file__).resolve().parent
H20=Path('/home/jay/build/wj-h20-three-level-resume-20261115')
H22=Path('/home/jay/build/wj-h22-completion-20261115')
CP_HASH=dict(C='647044a05d72b72f98526cecd5a731c6a2455a5debe4c10fb775192af44f7adc',
 M='6a1d061c8dee8c396ee588edd3afa86a353ceb58653c8bed2fedfdf291d3188c',F='bf6f580b5a7f3b6af9c88c64186371dd9d67803f4bc76bcdf3901d69e8a61760')

def xml_case(original,folder,anchor,z):
    folder.mkdir(parents=True);x=ET.parse(original).getroot()
    for k,v in (('NumSteps',ADDITIONAL[z]),('IO_CheckSteps',1),('IO_InfoSteps',1)):set_parameter(x,k,v)
    for k,v in (('RANSWallRestart',str(verify(anchor['checkpoint']).parent)),('RANSWallCheckpoint',str(folder/'CP')),
        ('RANSWallRetention','BoundedV1'),('RANSWallMonitorIntent',str(folder/'intent.json')),('RANSWallMonitorFault','Off'),
        ('RANSWallReport',str(folder/'native.json')),('RANSWallImmutableRowReuse','Off'),('RANSWallEvaluationReuse','ExactStepV1'),
        ('RANSWallMonitoringPolicy','CheckpointAuditV1'),('RANSWallPimpleSnapshots','Off'),('RANSWallPimpleArithmetic','Off'),
        ('RANSWallPimpleVerify','Off'),('RANSWallPimpleFlowVectors','Off'),('RANSWallPimpleSnapshotDirectory',str(folder/'snapshots'))):set_info(x,k,v)
    assert protected_case(x)==protected_case(ET.parse(original).getroot()),'changed protected XML/science/TimeStep'
    params={p.text.split('=')[0].strip():p.text.split('=')[1].strip()for p in x.findall('./CONDITIONS/PARAMETERS/P')}
    assert float(params['Time'])==0 and float(params['TimeStep'])==DT[z]and int(params['NumSteps'])==ADDITIONAL[z]
    dest=folder/'case.xml';write(x,dest)
    atomic(folder/'execution-diff.json',dict(original=identity(original),actual=identity(dest),TimeStep=DT[z],Time=0,
        additional_NumSteps=ADDITIONAL[z],restored_step=START[z],maximum_absolute_step=END[z],
        protected_sha256=hashlib.sha256(protected_case(x)).hexdigest(),physical_signature_definition_changed=False,
        allowed='NumSteps/IO and existing restart/checkpoint/report/monitor/capture SolverInfo only; rich captures Off'))

def setup(root):
    indices=load(root/'inputs/historical-indices.json')
    for index in indices.values():
        for pin in index.values():verify(pin)
        assert load(verify(index['accounting']))['closed']
    pub=load(verify(indices['H20']['publication']));published={p['path']:p for p in pub['files']}
    def published_file(relative):return verify(published[str(H20/'publication-v2'/relative)])
    heavy=load(published_file('heavy-manifest.json'));sealed={p['path']:p for p in heavy['files']}
    prior=load(H22/'inputs/selected.json')
    anchors={};pins={}
    for z in ORDER:
        n=START[z];d=H20/f'science/{z}/monitor';a=recovery(d,n)
        for k,p in a.items():
            if p['path']in sealed:assert p==sealed[p['path']];verify(p)
            elif k=='completion':assert p['sha256']==identity(published_file(f'{z}/step-{n}.completion.json'))['sha256']
            else:raise AssertionError('unsealed anchor: '+k)
        a['qualified_commit']=identity(d/f'step-{n}.qualification-commit.json')
        assert a['qualified_commit']['sha256']==identity(published_file(f'{z}/step-{n}.qualification-commit.json'))['sha256']
        assert a['checkpoint']['sha256']==CP_HASH[z]
        r,cp=qualified(a,z,n);anchors[z]=a
        for p in a.values():pins[p['path']]=p
    old_result=sealed[str(H20/'science/F/monitor/step-16.result.json')];result=load(verify(old_result))
    for p in (result['inventory'],result['signed_arrays']):assert p==sealed[p['path']];verify(p)
    selected=dict(schema='WJH23SelectedV1',anchors=anchors,cp0=prior['cp0'],original=prior['original'],common_support=prior['common_support'],
        sparse_structure=prior['sparse_structure'],historical_indices=indices,preflight_oracle=old_result,
        preflight_inventory=result['inventory'],preflight_fields=result['signed_arrays'])
    source=checkpoint_payload(verify(selected['original']['checkpoint']).parent,include_reference=True)
    for z,a in anchors.items():
        p=checkpoint_payload(verify(selected['cp0'][z]['checkpoint']).parent,include_reference=True)
        assert p['steps']==0 and p['time']==0 and p['dt']==DT[z]
        for k in ('state','fixed_reference','fixed_defect','eta'):equal(source[k],p[k],'own CP0 common spatial '+k)
        r,_=qualified(a,z,START[z])
        for k in ('fixed_reference','fixed_defect'):equal(r[k],source[k],'frozen '+k)
        assert r['original_qh0_state_sha256']==selected['original']['coupled_state_sha256']
        for pin in selected['cp0'][z].values():
            if isinstance(pin,dict)and 'path'in pin:verify(pin);pins[pin['path']]=pin
    atomic(root/'inputs/selected.json',selected)
    for name in ('runtime-build-v1.json','runtime-dependencies-v1.json','compiled-source-objects.json'):
        p=H22/'science'/name;value=load(p);pins[str(p)]=identity(p);atomic(root/'science'/name,value)
    build=load(root/'science/runtime-build-v1.json');assert build['native_source_revision']=='057e6b7465660cac5ec8617abf5c876fa8bff71f'
    for v in build['executables']:verify(v['ELF']);verify(v['wrapper'])
    verify(build['compile_commands'])
    for v in build['compiled_inputs']:
        assert identity(v['path'])['sha256']==v['sha256']
        assert all(f in v['command']for f in ('-ffp-contract=off','-fno-fast-math','-O3'))
    objects=load(root/'science/compiled-source-objects.json')
    for v in objects['objects']:
        verify(v['source']);verify(v['object'])
        if v['dependency_file']:verify(v['dependency_file'])
        if v.get('link_command'):verify(v['link_command'])
    for v in objects['local_dependencies']+objects['static_libraries']:verify(v)
    for v in load(root/'science/runtime-dependencies-v1.json')['libraries']:verify(v['actual'])
    previous=load(H22/'science/target-cache-compatibility.json')
    c=compatibility(root,verify(previous['historical_manifest']),root/'science/cache-base-compatibility.json')
    c.update(schema='WJH23ImmutableCompatibilityV1',permitted_difference='continuation roles/time map only; same independent numerical modules and native capsule; captures Off',
        dt_independence='immutable geometry/sign/basis/quadrature/reference only, own signatures checked; no changing fields/profiles/history/decisions reused',
        native_builds=0,historical_pins_modified=False)
    atomic(root/'science/target-cache-compatibility.json',c)
    for z in ORDER:xml_case(verify(selected['cp0'][z]['xml']),root/'science'/z,anchors[z],z)
    for name,index in indices.items():
        for pin in index.values():pins[pin['path']]=pin
    atomic(root/'inputs/historical-immutability.json',dict(files=list(pins.values()),scope='selected original closed inputs only; no old-root traversal'))
    atomic(root/'science/setup.json',dict(status='PASS',anchors=anchors,common_qh0=vector_hash(source['state']),new_builds=0,
        source_revision=build['native_source_revision'],runtime_ELFs=[v['ELF']for v in build['executables']],
        storage_forecast=dict(checkpoint_and_completion_packages=28*6500000,scientific_inventories_fields_results=15*65000000,
            strict_and_reports_logs=150000000,total_upper_bytes=1307000000,cap_bytes=2*1024**3)))

def preflight(root):
    inputs=load(root/'inputs/selected.json');cache,_=load_compatible(root/'science/target-cache-compatibility.json')
    old=load(verify(inputs['preflight_oracle']));inventory=load(verify(inputs['preflight_inventory']))
    report,_=qualified(inputs['anchors']['F'],'F',16);out=root/'science/preflight';out.mkdir()
    result,inv=evaluate(root,cache,load(verify(inputs['common_support']))['points'],inputs,report,'F',16,out,'retained-02-triplet',inventory['points'])
    for z in DT:
        equal(result['physical'][z],old['physical'][z],'retained .02 physical '+z)
        for kind in ('reference','drift'):equal(result['errors'][z][kind],old['errors'][z][kind],'retained error '+z+kind)
    equal(result['comparison'],old['comparison'],'retained all-pair and contraction')
    assert result['status']==result['comparison']['status']=='PASS'
    assert all(all(v['maximum']['value']==v['rms']['value']==0 for values in result['errors'][z]['change_since_02'].values()for v in values.values())for z in DT)
    result.update(binding_role='passive retained H20 .02 triplet, not a new native audit/barrier',oracle=inputs['preflight_oracle'])
    atomic(out/'result.json',result)

def freeze(root):
    assert load(root/'science/preflight/result.json')['status']=='PASS'
    assert load(root/'operations/focused-controls-v2/receipt.json')['status']=='PASS'
    files=list(load(root/'source-manifest.json')['files'])
    files.extend(identity(root/'inputs'/n)for n in ('selected.json','historical-indices.json','historical-immutability.json'))
    files.extend(identity(root/'science'/n)for n in ('runtime-build-v1.json','runtime-dependencies-v1.json','compiled-source-objects.json','target-cache-compatibility.json','setup.json','preflight/result.json'))
    files.extend(identity(root/f'science/{z}/{n}')for z in DT for n in ('case.xml','execution-diff.json'))
    atomic(root/'execution-config.json',dict(schema='WJH23FrozenV1',files=files,order=list(ORDER),start=START,end=END,additional=ADDITIONAL,
        DT=DT,strides=STRIDE,maximum_attempts=28,maximum_governors=3,maximum_strict=3,builds=0,captures='Off',
        forecasts=dict(startup=150,ordinary=45,native_audit=50,observer=OBSERVER_FORECAST,exit=30,strict=120),
        thresholds='unchanged H19/H20; contraction never stops',native_caps=dict(M=1050,C=900,F=1650,strict=120)))

def admission(root,level):
    frozen(root);b=load(root/'accounting.json')
    assert not any(b.get(k)for k in ('native_halt','integrity_halt','coverage_halt','screen_halt','runtime_halt','resource_halt'))
    earlier=ORDER[:ORDER.index(level)]
    for z in earlier:
        assert load(root/f'science/{z}/qualification.json')['status']=='PASS'
        strict=load(root/f'science/{z}/final-strict/result.json')
        assert strict['status']=='PASS'and integer(strict['selected_absolute_step'])==END[z]
    u=usage(b);remaining=ORDER[ORDER.index(level):];costs=forecast_costs(root)
    estimate=sum(branch_forecast(z,**costs[z])for z in remaining)
    assert u['native']+1.25*estimate<=4200 and u['total']+1.25*estimate<=6000,'remaining measured-work headroom'
    path=root/f'admission-{level}.json'
    if not path.exists():atomic(path,dict(status='PASS',level=level,used=u,remaining_estimate=estimate,headroom=1.25,strict_protected=120*len(remaining),closure_protected=300,
        forecasts=costs,forecast_rule='maximum of frozen conservative prior and measured preceding branches; no cap enlargement'))

def forecast_costs(root):
    costs={z:dict(ordinary=45.,audit=50.,observer=float(OBSERVER_FORECAST[z]))for z in ORDER}
    for z in ORDER:
        path=root/f'science/measured-forecast-{z}.json'
        if path.exists():
            v=load(path)
            for level in ORDER:
                for k in costs[level]:costs[level][k]=max(costs[level][k],v['selected_forecasts'][level][k])
    return costs

def measure_forecast(root,level):
    d=root/f'science/{level}/monitor'
    tele=[json.loads(s)for s in(d/'telemetry.jsonl').read_text().splitlines()]
    scopes=[json.loads(s)for s in(root/f'operations/governor-{level}-v1/timing.jsonl').read_text().splitlines()]
    ordinary=[float(t['native_compute_seconds'])for t in tele if integer(t['absolute_step'])%STRIDE[level]]
    def times(name):return [s['wall_seconds']for s in scopes if s['name']==name and START[level]<s['step']<=END[level]]
    audits=times('native.acceptance_audit');waits=times('monitor.independent_wait');costs=forecast_costs(root)
    if ordinary:
        for z in ORDER:costs[z]['ordinary']=max(costs[z]['ordinary'],max(ordinary))
    if audits:
        for z in ORDER:costs[z]['audit']=max(costs[z]['audit'],max(audits))
    if waits:costs[level]['observer']=max(costs[level]['observer'],max(waits))
    atomic(root/f'science/measured-forecast-{level}.json',dict(level=level,ordinary_samples=ordinary,native_audit_samples=audits,
        inclusive_wait_samples=waits,selected_forecasts=costs,rule='conservative max(prior, observed maximum); governing C all-audited is not an ordinary sample',
        timings_source=identity(root/f'operations/governor-{level}-v1/timing.jsonl')))

def target(root,level):
    frozen(root);inputs=load(root/'inputs/selected.json');folder=root/'science'/level;directory=folder/'monitor';directory.mkdir();ip=folder/'intent.json'
    costs=forecast_costs(root)[level]
    intent=dict(schema='WJMonitorIntentV1',policy='CheckpointAuditV1',h23_role='ContinuationFrom02V1',level=level,directory=str(directory),nonce=uuid.uuid4().hex,
        scope_sha256=hashlib.sha256((identity(folder/'case.xml')['sha256']+identity(root/'execution-config.json')['sha256']).encode()).hexdigest(),
        source_sha256=identity(root/'science/runtime-build-v1.json')['sha256'],observer_sha256=identity(HERE/'wj_h23_observer.py')['sha256'],inventory_rule_sha256=hashlib.sha256(RULE.encode()).hexdigest(),
        original_qh0_checkpoint=inputs['cp0'][level]['checkpoint'],original_qh0_state_sha256=inputs['original']['coupled_state_sha256'],spatial_lineage_anchor=inputs['original']['checkpoint'],
        restored_checkpoint=inputs['anchors'][level]['checkpoint'],restart_completion=inputs['anchors'][level]['completion'],
        maximum_absolute_step=END[level],restored_absolute_step=START[level],recovery_stride=1,scientific_stride=STRIDE[level],dt_bits=bits(DT[level]),ack_timeout_seconds=120.,
        substep_evidence_steps=[],flow_capture='Off',snapshot_capture='Off',preserve_qualification_artifacts=True,selected_inputs=identity(root/'inputs/selected.json'),common_support=inputs['common_support'],
        deadline_monotonic=float(os.environ['WJ_H23_DEADLINE']),operation_start_monotonic=float(os.environ['WJ_H23_OPERATION_START']),
        ordinary_step_forecast_seconds=costs['ordinary'],audit_forecast_seconds=costs['audit'],observer_forecast_seconds=costs['observer'],prospective_forecasts=forecast_costs(root))
    validate_intent(intent,ADDITIONAL[level]);obs=native=None;code=None;error=None;start=time.monotonic()
    with(folder/'observer.log').open('x')as ol,(folder/'native.log').open('x')as nl:
        try:
            obs=subprocess.Popen([sys.executable,'-B',str(HERE/'wj_h23_observer.py'),'--root',str(root),'--intent',str(ip)],stdout=ol,stderr=subprocess.STDOUT)
            intent['observer_pid']=obs.pid;atomic(ip,intent)
            while not(directory/'observer-started.json').exists():assert obs.poll()is None and time.monotonic()-start<30;time.sleep(.025)
            native=subprocess.Popen([str(executable(root,'IncNavierStokesSolver')),str(folder/'case.xml')],stdout=nl,stderr=subprocess.STDOUT)
            atomic(folder/'native-launch.json',dict(pid=native.pid,xml=identity(folder/'case.xml'),intent=identity(ip),additional=ADDITIONAL[level],restored=START[level]))
            while native.poll()is None:assert obs.poll()is None,'observer exited before native final evaluation';time.sleep(.025)
            code=native.returncode
        except Exception as e:error=repr(e)
        finally:
            if native and native.poll()is None:native.terminate();native.wait(timeout=5)
            atomic(directory/'observer-stop.json',dict(native_exit=code,reason='normal coordinator stop after native exit/failure'))
            if obs:
                try:obs.wait(timeout=5)
                except subprocess.TimeoutExpired:obs.terminate();obs.wait(timeout=5)
    path=folder/'case.evolution.json';r=native_load(path)if path.exists()else None
    atomic(folder/'outcome.json',dict(native_exit=code,error=error,producer=identity(path)if r else None,seconds=time.monotonic()-start))
    if r and os.environ.get('WJ_H23_OPERATION'):atomic(Path(os.environ['WJ_H23_OPERATION'])/'attempt-counts.json',dict(actual_attempts=integer(r['implicit_callbacks']),producer=identity(path)))
    assert error is None,error
    assert r is not None,'native report absent, original logs preserved'
    assert r['initialization_work']==[0,0]and r['initializer_calls']==0
    completed=[]
    for n in range(START[level]+1,integer(r['steps'])+1):
        anchor=recovery(directory,n);rr,_=read_completed(anchor,level,n);a=rr['pimple_attempts'][0]
        assert all(a[k]==1 for k in ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls'))
        assert a['bound_activations']==a['additional_correction_or_restoration_solves']==0 and a['k_bound']['raw_sha256']==a['k_bound']['corrected_sha256']
        assert not a.get('k_snapshot')and not a.get('flow_vectors'),'normal rich capture unexpectedly present'
        due=n%STRIDE[level]==0
        assert rr['native_full_audit_status']==('PASS'if due else'NOT_SCHEDULED')
        if not due:assert rr['native_hard_audit']=={}and not(directory/f'step-{n}.ready.json').exists()
        completed.append(dict(step=n,time=rr['time'],state=vector_hash(rr['state']),callback=anchor['callback'],completion=anchor['completion'],native_full=rr['native_full_audit_status']))
    if code:raise RuntimeError(r['failure'])
    assert integer(r['steps'])==END[level]and integer(r['implicit_callbacks'])==integer(r['native_completed'])==ADDITIONAL[level]
    assert integer(r['full_audit_calls'])==integer(r['physical_samples_calls'])==5
    assert integer(load(directory/'latest-fully-verified-checkpoint.json')['binding']['absolute_step'])==END[level]
    atomic(folder/'qualification.json',dict(status='PASS',attempts=ADDITIONAL[level],native_completed=ADDITIONAL[level],new_qualified_endpoints=4,
        startup_commit=identity(directory/'startup.verification-commit.json'),completed=completed,producer=identity(path),native_full_audits=5))

def strict(root,level):
    frozen(root);inputs=load(root/'inputs/selected.json');folder=root/'science'/level
    success=(folder/'qualification.json').exists()and load(folder/'qualification.json')['status']=='PASS'
    request,selection=strict_request(folder/'monitor',inputs['anchors'][level],level,expect_terminal=success)
    out=folder/'final-strict';out.mkdir();atomic(out/'request.json',request)
    with(out/'native.log').open('x')as log:code=subprocess.run([str(executable(root,'RANSWallResolutionProbe')),str(folder/'case.xml'),str(out/'request.json'),str(out/'native.json')],stdout=log,stderr=subprocess.STDOUT).returncode
    assert code==0;r=native_load(out/'native.json')
    assert r['initialization_work']==[0,0]and r['initial_velocity_calls']==r['governing_solves']==0
    assert integer(r['selected_absolute_step'])==selection['step']and r['weak_checks']['alpha_mass']==1/DT[level]
    old,_=read_completed(selection['anchor'],level,selection['step'])
    assert r['callback_history_sha256']==vector_hash(old['last_stage_history'])
    atomic(out/'result.json',dict(status='PASS',selected_absolute_step=selection['step'],time=selection['time'],genuine_alpha=selection['alpha'],selection=selection,
        producer=identity(out/'native.json'),request=identity(out/'request.json'),physical_mass_dual_norms=r['physical_mass_dual_norms'],weak=r['weak_checks'],interface=r['open_boundary']))
    if success:measure_forecast(root,level)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('setup','preflight','freeze','admission',*ORDER,'strict_M','strict_C','strict_F','runtime'),required=True)
    p.add_argument('--level',choices=ORDER);a=p.parse_args();assert a.root.name==ROOT
    if a.mode in ORDER:target(a.root,a.mode)
    elif a.mode.startswith('strict_'):strict(a.root,a.mode.split('_')[1])
    elif a.mode=='admission':admission(a.root,a.level)
    elif a.mode=='runtime':runtime_qualify(a.root)
    else:globals()[a.mode](a.root)
