#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H15 source-bound preparation, small fixtures and the single warm block."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import uuid
import xml.etree.ElementTree as ET
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import PersistentRows,identity,verify,atomic,vector_hash,RULE,H14,source_identities,reference_identity
from study_wj_h14 import runtime as prepare_runtime
from wj3_measure import checkpoint_payload
from cls_epawf_m3_cases import set_info,set_parameter,write

HERE=Path(__file__).resolve().parent

def setup(root):
    prepare_runtime(root)
    build=load(root/'science/runtime-build.json');build['package']='H15';build['native_source_revision']='2ffd38f061073d17b4b595a8795f31813606862f'
    atomic(root/'science/runtime-build-v1.json',build)
    parent=load(H14/'publication-v2/state-bindings.json');selected=[]
    for v in parent['states']:
        for k in ('producer','checkpoint','xml'):verify(v[k]);selected.append(v[k])
    old=load(H14/'closure-v1/accounting-final.json');assert old['closed']and old['total_charged_seconds']==2252.529779124001
    atomic(root/'science/imports.json',dict(parent_manifest=identity(H14/'publication-v2/manifest.json'),selected=selected,aliases=[],prior_ledger_unchanged=identity(H14/'closure-v1/accounting-final.json')))

def small_cache(root):
    out=root/'science/small-cache-v2';out.mkdir()
    control=load(H14/'publication-v2/controls/restart-rollback-checks.json')['counts']['cp0']['producer'];r=load(verify(control))
    old=load(H14/'science/small/cp0/assessment-v2/result.json');inv=load(verify(old['inventory']));points=[v['physical']for v in inv['points']]
    ref=reference_identity(root,r)
    cache=PersistentRows(r,points,r['state'],ref)
    import pickle
    p=out/'immutable-cache.pickle'
    with p.open('xb')as f:pickle.dump(cache,f,protocol=5);f.flush();os.fsync(f.fileno())
    atomic(out/'manifest.json',dict(status='PASS',cache=identity(p),reference_identity=ref,source_files=source_identities(),baseline_points=points,build_seconds=cache.build_seconds,original_report=control))

def xml_case(original,folder,cp,steps,retention=True,fault='Off',hard=False):
    folder.mkdir(parents=True);x=ET.parse(original).getroot()
    for k,v in (('NumSteps',steps),('IO_CheckSteps',1),('IO_InfoSteps',1)):set_parameter(x,k,v)
    for k,v in (('RANSWallRestart',str(cp)),('RANSWallCheckpoint',str(folder/'CP')),('RANSWallPimpleSnapshotDirectory',str(folder/'snapshots')),
        ('RANSWallRetention','BoundedV1'if retention else'Off'),('RANSWallMonitorIntent',str(folder/'intent.json')),('RANSWallMonitorFault',fault),('RANSWallReport',str(folder/'native.json'))):
        set_info(x,k,v)
    if hard:set_parameter(x,'RANSWallPimpleHardAuditRejectStep',1)
    path=folder/'case.xml';write(x,path)
    # Time and dt are physical signature inputs. Only listed execution controls
    # are changed; no Function, expansion, geometry or other parameter can drift.
    def protected(t):
        import copy
        t=copy.deepcopy(t)
        for p in t.findall('./CONDITIONS/PARAMETERS/P'):
            if p.text.split('=')[0].strip()in ('NumSteps','IO_CheckSteps','IO_InfoSteps','RANSWallPimpleHardAuditRejectStep'):p.text='execution'
        parent=t.find('./CONDITIONS/SOLVERINFO')
        for p in list(parent):
            if p.get('PROPERTY')in ('RANSWallRestart','RANSWallCheckpoint','RANSWallPimpleSnapshotDirectory','RANSWallRetention','RANSWallMonitorIntent','RANSWallMonitorFault','RANSWallReport'):parent.remove(p)
        # Deliberate rejection parameter is absent in the original; compare the
        # scientific base separately and document that one existing test hook.
        parent=t.find('./CONDITIONS/PARAMETERS')
        for p in list(parent):
            if p.text=='execution':parent.remove(p)
        for p in t.iter():
            if p.text is not None and not p.text.strip():p.text=None
            p.tail=None
        return ET.tostring(t)
    assert protected(ET.parse(original).getroot())==protected(x)
    atomic(folder/'execution-diff.json',dict(original=identity(original),as_run=identity(path),protected_sha256=hashlib.sha256(protected(x)).hexdigest(),additional_steps=steps,retention=retention,fault=fault,deliberate_hard_rejection=hard))
    return path

def run_observed(root,folder,cache_manifest,original_cp0,*,behavior='PASS',restart_callback=None,replay=None,deadline=None):
    for v in load(root/'science/runtime-build-v1.json')['executables']:
        verify(v['ELF']);verify(v['wrapper'])
    xml=folder/'case.xml';cp0=checkpoint_payload(original_cp0,include_reference=True)
    additional=int(next(p.text.split('=')[1]for p in ET.parse(xml).findall('./CONDITIONS/PARAMETERS/P')if p.text.split('=')[0].strip()=='NumSteps'));maximum=cp0['steps']+additional
    # For restarted small fixtures the original qh0 remains CP0, not CP1.
    if restart_callback:maximum=2
    directory=folder/'monitor';directory.mkdir();intent_path=folder/'intent.json'
    source=identity(root/'science/runtime-build-v1.json')['sha256'];observer=identity(HERE/'wj_h15_observer.py')['sha256']
    intent=dict(schema='WJMonitorIntentV1',directory=str(directory),nonce=uuid.uuid4().hex,scope_sha256=hashlib.sha256((source+identity(xml)['sha256']).encode()).hexdigest(),
        source_sha256=source,observer_sha256=observer,inventory_rule_sha256=hashlib.sha256(RULE.encode()).hexdigest(),
        original_qh0_checkpoint=identity(Path(original_cp0)/'State.bin'),original_qh0_state_sha256=vector_hash(cp0['state']),maximum_absolute_step=maximum,
        ack_timeout_seconds=.25 if behavior=='timeout'else 90.,strict_reserve_seconds=120)
    if restart_callback:intent['restart_callback']=identity(restart_callback)
    if replay:intent['replay_H14_step1']=replay
    if deadline:intent['deadline_monotonic']=deadline
    with(folder/'observer.log').open('x')as obslog,(folder/'native.log').open('x')as log:
        obs=subprocess.Popen([sys.executable,'-B',str(HERE/'wj_h15_observer.py'),'--root',str(root),'--mode','daemon','--intent',str(intent_path),'--cache',str(cache_manifest),'--behavior',behavior],stdout=obslog,stderr=subprocess.STDOUT)
        intent['observer_pid']=obs.pid;atomic(intent_path,intent)
        native=None;started=time.monotonic()
        try:
            while not(directory/'observer-started.json').exists():
                assert obs.poll()is None,'persistent observer failed before native launch'
                assert time.monotonic()-started<30,'observer setup timeout';time.sleep(.025)
            native=subprocess.Popen([str(root/'runtime/IncNavierStokesSolver'),str(xml)],stdout=log,stderr=subprocess.STDOUT)
            atomic(folder/'native-launch.json',dict(pid=native.pid,xml=identity(xml),intent=identity(intent_path),reserved_maximum_attempts=additional))
            while native.poll()is None:
                if obs.poll()is not None:native.terminate();raise RuntimeError('observer exited while native waited; no further advance eligible')
                time.sleep(.025)
            return native.returncode,load(folder/'case.evolution.json')
        finally:
            if native and native.poll()is None:native.terminate();native.wait(timeout=5)
            atomic(directory/'observer-stop.json',dict(reason='coordinator closure'))
            try:obs.wait(timeout=5)
            except subprocess.TimeoutExpired:obs.terminate();obs.wait(timeout=5)

def small(root):
    parent=H14/'science/small/cp0';base=parent/'case.xml';cp=Path(load(parent/'case.evolution.json')['checkpoint']);out=root/'science/small';out.mkdir()
    cache=root/'science/small-cache-v2/manifest.json';records=[]
    xml_case(base,out/'legacy',cp,2,retention=False)
    with(out/'legacy/native.log').open('x')as log:code=subprocess.run([str(root/'runtime/IncNavierStokesSolver'),str(out/'legacy/case.xml')],stdout=log,stderr=subprocess.STDOUT).returncode
    assert code==0;legacy=load(out/'legacy/case.evolution.json')
    for name,steps,behavior,fault,hard in (('pass',2,'PASS','Off',False),('restart',1,'PASS','Off',False),
        ('observer-fail',2,'FAIL','Off',False),('timeout',2,'timeout','Off',False),('stale',2,'stale','Off',False),
        ('native-reject',1,'PASS','Off',True),('post-callback',1,'PASS','AfterImplicitV1',False)):
        restart=cp;side=None
        if name=='restart':
            success=load(out/'pass/case.evolution.json');cursor=load(out/'pass/monitor/step-1.qualification-commit.json');restart=Path(cursor['checkpoint']['path']).parent;side=cursor['callback']['path']
        xml_case(base,out/name,restart,steps,fault=fault,hard=hard)
        code,r=run_observed(root,out/name,cache,cp,behavior=behavior,restart_callback=side)
        expected=name in ('pass','restart');assert(code==0)==expected
        assert r['initialization_work']==[0,0]and r['initializer_calls']==0
        if expected:
            end=checkpoint_payload(r['last_qualified_checkpoint'],include_reference=True)
            assert r['steps']==2 and r['independently_qualified']==(2 if name=='pass'else 1)
            for k in ('state','fixed_reference','fixed_defect'):assert vector_hash(end[k])==vector_hash(legacy[k]),k
            assert all(vector_hash(a)==vector_hash(b)for a,b in zip(end['native_history'],legacy['native_history']))and end['native_history_times']==legacy['native_history_times']
            if name=='pass':
                for step in (1,2):
                    stage=load(out/f'pass/monitor/step-{step}.report.json')
                    assert vector_hash(stage['state'])==vector_hash(legacy['accepted_states'][step])
                    assert stage['wj_audit_packets'][0]['stage_history_sha256']==legacy['wj_audit_packets'][step]['stage_history_sha256']
        elif name in ('observer-fail','timeout','stale'):
            assert r['native_completed']==1 and r['independently_qualified']==0 and r['steps']==1
            assert Path(r['last_qualified_checkpoint'])==cp and Path(r['last_completed_checkpoint'])!=cp
            assert not(out/name/'monitor/step-2.ready.json').exists()and not(out/name/'monitor/step-1.qualification-commit.json').exists()
        else:
            assert r['steps']==0 and r['native_completed']==0 and r['aborted_advances']==1
            aborted=load(out/name/'monitor/first-aborted.json');assert aborted['status']=='ABORTED_NOT_COMPLETED'
            assert not(out/name/'monitor/step-1.ready.json').exists()
            original=checkpoint_payload(cp,include_reference=True)
            assert r['state_sha256']==vector_hash(original['state'])and r['native_history_sha256']==[vector_hash(v)for v in original['native_history']]
            assert r['native_history_times']==original['native_history_times']
        assert r['resident_accepted_states']<=2 and r['resident_wj_packets']<=2 and r['resident_pimple_records']<=1 and r['factor_catalog']['resident_definitions']<=16
        records.append(dict(case=name,native_exit=code,summary=identity(out/name/'case.evolution.json'),status='PASS',native_completed=r['native_completed'],independently_qualified=r['independently_qualified'],failure=r['failure']))
    atomic(out/'qualification.json',dict(status='PASS',small_attempts=10,records=records,legacy_final=identity(out/'legacy/case.evolution.json')))
    print('small native/protocol/restart/rollback qualification PASS',flush=True)

def strict_request(root,folder):
    cursor=load(folder/'monitor/qualified-cursor.json');report=load(verify(cursor['report']));cp=verify(cursor['checkpoint']).parent
    out=folder/'view';out.mkdir()
    request=dict(operation='open_checks',inner_stress_model='Linear',native_report=cursor['report']['path'],native_report_sha256=cursor['report']['sha256'],
        native_checkpoint=str(cp),checkpoint_sha256=cursor['checkpoint']['sha256'],completed_sidecar=cursor['callback'],
        qualified_cursor=identity(folder/f"monitor/step-{report['steps']}.qualification-commit.json"),weak=True,weak_compact=True,assembly=False,reuse=False,retain_scalar_observations=False)
    atomic(out/'request.json',request);return [str(root/'runtime/RANSWallResolutionProbe'),str(folder/'case.xml'),str(out/'request.json'),str(out/'native.json')]

def strict_finish(folder):
    p=load(folder/'view/native.json');req=load(folder/'view/request.json');assert p['initialization_work']==[0,0]and p['initial_velocity_calls']==p['governing_solves']==0
    assert p['report_sha256']==req['native_report_sha256']and p['checkpoint_sha256']==req['checkpoint_sha256']and p['scalar_domain']['admissible']
    atomic(folder/'view/result.json',dict(status='PASS',native=identity(folder/'view/native.json'),request=identity(folder/'view/request.json'),physical_mass_dual_norms=p['physical_mass_dual_norms'],weak=p['weak_checks'],interface=p['open_boundary'],selected_absolute_step=p['selected_absolute_step'],callback_history_sha256=p['callback_history_sha256']))

def runtime_qualify(root):
    """Check actual native-exec maps, excluding unrelated parent mappings."""
    old=H14.parent/'wj-h10-flat-inlet-20261113'
    historical=dict(path=str(old/'closure/runtime-observed.json'),bytes=307438,sha256='f8698dbe2108eb1f597696805b342f051c1ddac5cd2006fd655e5a682b43a502')
    observed=load(verify(historical));assert observed['status']=='PASS'
    prior=load(old/'runtime-qualified.json')['loaded']
    for v in observed['libraries']:prior[v['path']]=v
    build=load(root/'science/runtime-build-v1.json')
    for v in build['executables']:verify(v['ELF']);verify(v['wrapper'])
    for v in build['private_libraries']:verify(v['original']);verify(v['copy'])
    libraries={};samples=[];executables=set()
    for p in sorted((root/'operations').glob('*/native-exec-samples.json')):
        for s in load(p):
            assert s['status']=='PASS'
            for native in s['native']:
                assert native in [v['ELF']['path']for v in build['executables']]
                executables.add(Path(native).name)
            samples.append(dict(operation=p.parent.name,pid=s['pid'],native=s['native'],libc=s['libc']))
            for f in s['files']:
                if '.so'not in Path(f).name and not Path(f).name.startswith('ld-linux'):continue
                key=f
                if '/runtime/glibc/'in f:key=str(old/'runtime/glibc'/Path(f).name)
                elif '/runtime/geometry/'in f:key=str(old/'runtime/geometry'/Path(f).name)
                assert key in prior,'unqualified actual native dependency: '+f
                actual=identity(f);expected=prior[key]
                assert actual['sha256']==expected['sha256']and actual['bytes']==expected['bytes']
                libraries[f]=dict(actual=actual,qualified=expected)
    assert executables=={'IncNavierStokesSolver.elf','RANSWallResolutionProbe.elf'}
    atomic(root/'science/runtime-qualified-v1.json',dict(status='PASS',native_exec_samples=samples,libraries=list(libraries.values()),compiled_build=identity(root/'science/runtime-build-v1.json'),historical_runtime=historical))

def admission(root):
    assert load(root/'science/small/qualification.json')['status']=='PASS'and load(root/'science/observer-cache/manifest.json')['status']=='PASS'
    assert load(root/'science/small/pass/view/result.json')['status']=='PASS'
    assert load(root/'science/runtime-qualified-v1.json')['status']=='PASS'
    # Conservative aggregate memory: use H14 governor peak (strict reader is
    # not concurrent with the observer) plus measured cache qualification RSS.
    book=load(root/'accounting.json');old=load(H14/'closure-v1/accounting-final.json');native=old['entries']['target-cost-gated-pilot-v1']['peak_tree_rss_bytes']
    observer=book['entries']['retained-target-persistent-parity-v1']['peak_tree_rss_bytes'];combined=native+observer+128*1024**2
    forecast=1.25*(165+3*80+4*15+120+60+30)
    assert combined<8*1024**3 and forecast<1800
    cp0=next(v for v in load(H14/'publication-v2/state-bindings.json')['states']if v['case']=='cp0')
    folder=root/'science/target';xml_case(verify(cp0['xml']),folder,verify(cp0['checkpoint']).parent,4)
    command=[sys.executable,'-B',str(HERE/'study_wj_h15.py'),'--root',str(root),'--mode','target']
    inputs=[identity(HERE/'study_wj_h15.py'),identity(HERE/'wj_h15_observer.py'),identity(folder/'case.xml'),identity(root/'science/observer-cache/manifest.json'),identity(root/'science/runtime-build-v1.json'),identity(root/'science/runtime-qualified-v1.json')]
    for executable in load(root/'science/runtime-build-v1.json')['executables']:
        inputs.extend((executable['ELF'],executable['wrapper']))
    inputs.extend((cp0['producer'],cp0['checkpoint'],cp0['xml']))
    atomic(root/'target-admission.json',dict(status='PASS',height_mm=6,state_size=24688,command=command,inputs=inputs,
        cost_with25percent_headroom=forecast,phase_ceiling=1800,combined_RSS_forecast=combined,final_strict_reserved=120,closure_reserved=300))

def qualified_step(folder):
    cursor=folder/'monitor/qualified-cursor.json'
    if not cursor.exists():return 0
    committed=load(cursor);assert committed['schema']=='WJQualificationCommitV1'
    for key in ('checkpoint','callback','report','consumed_ack','prepared_cursor'):verify(committed[key])
    step=int(committed['binding']['absolute_step'])
    assert load(folder/f'monitor/step-{step}.qualification-commit.json')==committed
    return step

def final_cp0_intent(root):
    """One final strict view of the actual last qualified state after stop."""
    folder=root/'science/target';r=load(folder/'case.evolution.json')
    assert r['native_completed']==1 and r['independently_qualified']==0
    assert qualified_step(folder)==0 and r['failure']=='R3 durable progress: stale/duplicate observer ACK'
    cp0=next(v for v in load(H14/'publication-v2/state-bindings.json')['states']if v['case']=='cp0')
    assert Path(r['last_qualified_checkpoint'])==verify(cp0['checkpoint']).parent
    out=root/'science/final-qualified-cp0-view';out.mkdir()
    request=dict(operation='open_checks',inner_stress_model='Linear',native_report=str(verify(cp0['producer'])),native_report_sha256=cp0['producer']['sha256'],
        native_checkpoint=r['last_qualified_checkpoint'],checkpoint_sha256=cp0['checkpoint']['sha256'],weak=True,weak_compact=True,assembly=False,reuse=False,retain_scalar_observations=False)
    atomic(out/'request.json',request)
    command=[str(root/'runtime/RANSWallResolutionProbe'),str(verify(cp0['xml'])),str(out/'request.json'),str(out/'native.json')]
    atomic(root/'target-final-strict-admission.json',dict(status='PASS',command=command,inputs=[identity(out/'request.json'),cp0['producer'],cp0['checkpoint'],cp0['xml'],identity(root/'runtime/RANSWallResolutionProbe.elf'),identity(root/'runtime/RANSWallResolutionProbe')],
        role='last qualified checkpoint remains original H14 CP0; no completed-stage sidecar exists at CP0',completed_but_unqualified=identity(folder/'CP/State.bin'),native_failure=identity(folder/'case.evolution.json'),new_governing_attempts=0))

def target(root):
    began=time.monotonic();folder=root/'science/target';states=load(H14/'publication-v2/state-bindings.json')['states'];cp0=next(v for v in states if v['case']=='cp0');cp1=next(v for v in states if v['case']=='pilot1')
    code=None;failure=None
    try:
        code,r=run_observed(root,folder,root/'science/observer-cache/manifest.json',verify(cp0['checkpoint']).parent,replay=cp1['producer'],deadline=began+1800)
    except Exception as error:
        # A coordinator/protocol failure never authorizes another native block.
        # Still protect the one final read of an already committed qualified
        # cursor; a prepared observer cursor is insufficient.
        failure=repr(error)
    summary=folder/'case.evolution.json';r=load(summary)if summary.exists()else {}
    qualified=qualified_step(folder)
    atomic(folder/'native-outcome.json',dict(native_exit=code,summary=identity(summary)if summary.exists()else None,native_completed=r.get('native_completed'),independently_qualified=qualified,new_progress=max(0,qualified-1),failure=failure or r.get('failure')))
    if qualified:
        command=strict_request(root,folder)
        with(folder/'view/native.log').open('x')as log:strict_code=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT).returncode
        assert strict_code==0,'final strict readback failed';strict_finish(folder)
    print('target block',code,'completed',r.get('native_completed'),'qualified',qualified,'wall',time.monotonic()-began,flush=True)
    if failure:raise RuntimeError(failure)
    if code!=0:raise RuntimeError('native target stopped: '+str(r.get('failure')))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('setup','small-cache','small','runtime-qualify','admission','target','final-cp0-intent','small-strict-intent','small-strict-finish'),required=True);a=p.parse_args()
    if a.mode=='small-strict-intent':print(json.dumps(strict_request(a.root,a.root/'science/small/pass')))
    elif a.mode=='small-strict-finish':strict_finish(a.root/'science/small/pass')
    else:globals()[a.mode.replace('-','_')](a.root)
