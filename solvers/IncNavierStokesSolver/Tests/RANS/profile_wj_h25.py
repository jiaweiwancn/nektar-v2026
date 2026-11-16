#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Four advancing profiling-only steps; exact comparison with retained W0."""
import argparse,hashlib,os,subprocess,sys,time,uuid
from pathlib import Path
import xml.etree.ElementTree as ET
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,verify,atomic,vector_hash,bits,RULE
from wj_h17_observer import native_load,equal
from wj_h16_observer import compatibility
from cls_epawf_m3_cases import set_info,set_parameter,write

def run(root,variant='profiling'):
    assert variant in ('profiling','lookup-off','lookup-on')
    lookup=variant=='lookup-on';steps=2 if variant=='lookup-off'else 4
    tag='profile'if variant=='profiling'else'lookup'
    label='W0-profile-v1'if variant=='profiling'else ('W1-lookup-v1'if lookup else'W0-lookup-off-v1')
    terminal=16+steps
    here=Path(__file__).resolve().parent;inputs=load(root/'inputs/baseline.json')
    old=load(root/'science/cache-compatible-v2.json')
    c=compatibility(root,verify(old['historical_manifest']),root/('science/'+variant+'-cache-binding-base.json'))
    c.update(new_native_identity=identity(root/('science/runtime-build-'+tag+'-v1.json')),
        permitted_difference='native timing and opt-in exact interval/geometry lookup only; independent numerical modules/geometry/profile/initial/reference/defect unchanged')
    atomic(root/('science/'+variant+'-cache-compatible.json'),c)
    folder=root/('science/'+label);folder.mkdir();d=folder/'monitor';d.mkdir();ip=folder/'intent.json'
    x=ET.parse(root/'science/W0/case.xml').getroot();set_parameter(x,'NumSteps',steps);set_info(x,'RANSWallNativeLookup','TensorExactV1'if lookup else'Off')
    for k,v in dict(RANSWallCheckpoint=str(folder/'CP'),RANSWallMonitorIntent=str(ip),RANSWallReport=str(folder/'native.json')).items():set_info(x,k,v)
    write(x,folder/'case.xml')
    intent=dict(schema='WJMonitorIntentV1',policy='CheckpointAuditV1',directory=str(d),nonce=uuid.uuid4().hex,
        scope_sha256=identity(folder/'case.xml')['sha256'],source_sha256=identity(root/('science/runtime-build-'+tag+'-v1.json'))['sha256'],observer_sha256=identity(here/'wj_h25_observer.py')['sha256'],inventory_rule_sha256=hashlib.sha256(RULE.encode()).hexdigest(),
        original_qh0_checkpoint=inputs['own_cp0']['checkpoint'],original_qh0_state_sha256=inputs['original']['coupled_state_sha256'],spatial_lineage_anchor=inputs['original']['checkpoint'],
        restored_checkpoint=inputs['anchor']['checkpoint'],restart_completion=inputs['anchor']['completion'],completed_anchor=inputs['anchor'],
        maximum_absolute_step=terminal,restored_absolute_step=16,recovery_stride=1,scientific_stride=terminal,dt_bits=bits(.0025),ack_timeout_seconds=120.,
        substep_evidence_steps=[],flow_capture='Off',snapshot_capture='Off',preserve_qualification_artifacts=True,
        cache_compatibility=identity(root/('science/'+variant+'-cache-compatible.json')),common_support=inputs['support'])
    native=obs=None;code=None;error=None;began=time.monotonic()
    with(folder/'native.log').open('x')as nl,(folder/'observer.log').open('x')as ol:
        try:
            obs=subprocess.Popen([sys.executable,'-B',str(here/'wj_h25_observer.py'),'--root',str(root),'--intent',str(ip)],stdout=ol,stderr=subprocess.STDOUT)
            intent['observer_pid']=obs.pid;atomic(ip,intent)
            while not(d/'observer-started.json').exists():assert obs.poll()is None and time.monotonic()-began<30;time.sleep(.025)
            native=subprocess.Popen([str(root/('runtime-'+tag+'-v1/IncNavierStokesSolver')),str(folder/'case.xml')],stdout=nl,stderr=subprocess.STDOUT)
            while native.poll()is None:assert obs.poll()is None,'observer exited before final native work';time.sleep(.025)
            code=native.returncode
        except Exception as e:error=repr(e)
        finally:
            if native and native.poll()is None:native.terminate();native.wait(timeout=5)
            atomic(d/'observer-stop.json',dict(code=code,error=error))
            if obs:
                try:obs.wait(timeout=5)
                except subprocess.TimeoutExpired:obs.terminate();obs.wait(timeout=5)
    report=folder/'case.evolution.json';r=native_load(report)if report.exists()else None
    actual=int(r['implicit_callbacks'])if r else steps
    atomic(Path(os.environ['WJ_H25_OPERATION'])/'attempt-counts.json',dict(actual_attempts=actual,known=r is not None))
    atomic(folder/'outcome.json',dict(code=code,error=error,actual_attempts=actual,producer=identity(report)if r else None))
    assert code==0 and error is None,(code,error,r.get('failure')if r else None)
    records=[]
    for step in range(17,terminal+1):
        current=native_load(d/f'step-{step}.report.json');old=native_load(root/f'science/W0/monitor/step-{step}.report.json')
        for k in ('state','fixed_reference','fixed_defect','native_history','native_history_times','time','steps','dt','signature','last_stage_history','last_stage_alpha','last_stage_time'):
            equal(current[k],old[k],'profiling-only '+str(step)+' '+k)
        records.append(dict(step=step,time=current['time'],state=vector_hash(current['state']),producer=identity(d/f'step-{step}.report.json'),baseline=identity(root/f'science/W0/monitor/step-{step}.report.json')))
    assert int(load(d/'latest-fully-verified-checkpoint.json')['binding']['absolute_step'])==terminal
    if lookup:
        current=native_load(d/'step-20.report.json');base=native_load(root/'science/W0-profile-v1/monitor/step-20.report.json')
        equal(current['native_hard_audit'],base['native_hard_audit'],'complete native physical map')
    atomic(folder/'qualification.json',dict(status='PASS',variant=variant,lookup='TensorExactV1'if lookup else'Off',scope='exact numerical state/history/callback/reference/time parity for all four evolving steps; startup/final fresh physical audits',records=records))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--variant',choices=('profiling','lookup-off','lookup-on'),default='profiling');a=p.parse_args();run(a.root,a.variant)
