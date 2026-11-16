#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Real restart of the retained exact-lookup trajectory, with sealed W0 oracle."""
import argparse, hashlib, os, subprocess, sys, time, uuid
from pathlib import Path
import xml.etree.ElementTree as ET
from run_wj_h4 import THREADS, load
os.environ.update(THREADS)
from wj_h15_observer import identity, verify, atomic, bits, RULE
from wj_h16_observer import compatibility
from wj_h17_observer import native_load, equal
from cls_epawf_m3_cases import set_info, set_parameter, write


def run(root):
    here = Path(__file__).resolve().parent
    original = root/'science/W1-lookup-v1/monitor'
    commit = original/'step-20.qualification-commit.json'
    q = load(commit); assert q['schema'] == 'WJQualificationCommitV1'
    anchor = {k:q[k] for k in ('checkpoint', 'report', 'callback')}
    anchor.update(qualified_commit=identity(commit), completion=identity(original/'step-20.completion.json'))
    for pin in anchor.values(): verify(pin)
    baseline = load(root/'inputs/baseline.json')
    d = root/'science/P-lookup-restart-v1'; d.mkdir()
    monitor = d/'monitor'; monitor.mkdir(); ip = d/'intent.json'
    c = compatibility(root, verify(load(root/'science/cache-compatible-v2.json')['historical_manifest']), d/'cache-base.json')
    c.update(new_native_identity=identity(root/'science/runtime-build-scalar-v1.json'),
        permitted_difference='new native capsule, exact lookup On, scalar value path Off; all immutable row/profile/reference identities unchanged')
    atomic(d/'cache-compatible.json', c)
    x = ET.parse(baseline['own_cp0']['xml']['path']).getroot()
    for k,v in dict(NumSteps=4, IO_CheckSteps=1, IO_InfoSteps=1).items(): set_parameter(x,k,v)
    for k,v in dict(RANSWallRestart=str(verify(anchor['checkpoint']).parent), RANSWallCheckpoint=str(d/'CP'),
        RANSWallRetention='BoundedV1', RANSWallMonitoringPolicy='CheckpointAuditV1', RANSWallMonitorIntent=str(ip),
        RANSWallMonitorFault='Off', RANSWallReport=str(d/'native.json'), RANSWallPimpleSnapshots='Off',
        RANSWallPimpleFlowVectors='Off', RANSWallPimpleArithmetic='Off', RANSWallPimpleVerify='Off',
        RANSWallImmutableRowReuse='Off', RANSWallNativeLookup='TensorExactV1', RANSWallScalarEvaluation='Off').items(): set_info(x,k,v)
    write(x,d/'case.xml')
    v = dict(schema='WJMonitorIntentV1', policy='CheckpointAuditV1', directory=str(monitor), nonce=uuid.uuid4().hex,
        scope_sha256=identity(d/'case.xml')['sha256'], source_sha256=identity(root/'science/runtime-build-scalar-v1.json')['sha256'],
        observer_sha256=identity(here/'wj_h25_observer.py')['sha256'], inventory_rule_sha256=hashlib.sha256(RULE.encode()).hexdigest(),
        original_qh0_checkpoint=baseline['own_cp0']['checkpoint'], original_qh0_state_sha256=baseline['original']['coupled_state_sha256'],
        spatial_lineage_anchor=baseline['original']['checkpoint'], completed_anchor=anchor, restored_checkpoint=anchor['checkpoint'],
        restart_completion=anchor['completion'], maximum_absolute_step=24, restored_absolute_step=20,
        recovery_stride=1, scientific_stride=24, dt_bits=bits(.0025), ack_timeout_seconds=120.,
        substep_evidence_steps=[], flow_capture='Off', snapshot_capture='Off', preserve_qualification_artifacts=True,
        cache_compatibility=identity(d/'cache-compatible.json'), common_support=baseline['support'])
    native = observer = None; code = None; error = None; started = time.monotonic()
    with (d/'native.log').open('x') as nl, (d/'observer.log').open('x') as ol:
        try:
            observer = subprocess.Popen([sys.executable,'-B',str(here/'wj_h25_observer.py'),'--root',str(root),'--intent',str(ip)],stdout=ol,stderr=subprocess.STDOUT)
            v['observer_pid'] = observer.pid; atomic(ip,v)
            while not (monitor/'observer-started.json').exists():
                assert observer.poll() is None and time.monotonic()-started < 30
                time.sleep(.025)
            native = subprocess.Popen([str(root/'runtime-scalar-v1/IncNavierStokesSolver'),str(d/'case.xml')],stdout=nl,stderr=subprocess.STDOUT)
            while native.poll() is None:
                assert observer.poll() is None, 'observer exited before native'
                time.sleep(.025)
            code = native.returncode
        except Exception as e: error = repr(e)
        finally:
            if native and native.poll() is None: native.terminate(); native.wait(timeout=5)
            atomic(monitor/'observer-stop.json',dict(native_exit=code,error=error))
            if observer:
                try: observer.wait(timeout=5)
                except subprocess.TimeoutExpired: observer.terminate(); observer.wait(timeout=5)
    report = d/'case.evolution.json'; r = native_load(report) if report.exists() else None
    actual = int(r['implicit_callbacks']) if r else 4
    atomic(Path(os.environ['WJ_H25_OPERATION'])/'attempt-counts.json',dict(actual_attempts=actual,known=r is not None))
    atomic(d/'outcome.json',dict(native_exit=code,error=error,seconds=time.monotonic()-started,producer=identity(report) if r else None))
    assert code == 0 and error is None
    assert r['initialization_work'] == [0,0] and r['initializer_calls'] == 0
    records = []
    for step in range(21,25):
        a = native_load(monitor/f'step-{step}.report.json'); b = native_load(root/f'science/W0/monitor/step-{step}.report.json')
        for key in ('state','fixed_reference','fixed_defect','eta','native_history','native_history_times','time','steps','dt','signature','last_stage_history','last_stage_alpha','last_stage_time'):
            equal(a[key],b[key],'real H25 restart '+key)
        attempt = a['pimple_attempts'][0]
        assert all(attempt[k] == 1 for k in ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls'))
        assert attempt['bound_activations'] == attempt['additional_correction_or_restoration_solves'] == 0
        records.append(dict(step=step,producer=identity(monitor/f'step-{step}.report.json'),baseline=identity(root/f'science/W0/monitor/step-{step}.report.json')))
    equal(a['native_hard_audit'],b['native_hard_audit'],'real restart terminal full native audit')
    cursor = load(monitor/'latest-fully-verified-checkpoint.json'); assert int(cursor['binding']['absolute_step']) == 24
    atomic(d/'qualification.json',dict(status='PASS_RESTART_AND_UNPERTURBED_PILOT',records=records,source_anchor=anchor,
        selected=dict(NativeLookup='TensorExactV1',ScalarEvaluation='Off'),endpoint=identity(monitor/'step-24.qualification-commit.json'),
        limitation='four-step unperturbed exact restart; no transient-dt or production-duration admission'))


def strict(root):
    d = root/'science/P-lookup-restart-v1'; q = load(d/'monitor/step-24.qualification-commit.json')
    assert load(d/'qualification.json')['status'] == 'PASS_RESTART_AND_UNPERTURBED_PILOT'
    out = d/'strict'; out.mkdir()
    request = dict(operation='open_checks', inner_stress_model='Linear', import_role='QualifiedV1',
        native_report=q['report']['path'], native_report_sha256=q['report']['sha256'],
        native_checkpoint=str(verify(q['checkpoint']).parent), checkpoint_sha256=q['checkpoint']['sha256'],
        completed_sidecar=q['callback'], qualified_cursor=identity(d/'monitor/step-24.qualification-commit.json'),
        weak=True, weak_compact=True, assembly=False, reuse=False, retain_scalar_observations=False)
    atomic(out/'request.json',request)
    with (out/'native.log').open('x') as f:
        subprocess.run([str(root/'runtime-scalar-v1/RANSWallResolutionProbe'),str(d/'case.xml'),str(out/'request.json'),str(out/'native.json')],stdout=f,stderr=subprocess.STDOUT,check=True)
    r = native_load(out/'native.json')
    assert r['initialization_work'] == [0,0] and r['initial_velocity_calls'] == r['governing_solves'] == 0
    atomic(out/'qualification.json',dict(status='PASS',request=identity(out/'request.json'),result=identity(out/'native.json'),step=24,time=.06))


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True,type=Path);p.add_argument('--mode',choices=('run','strict'),required=True)
    a=p.parse_args();globals()[a.mode](a.root)
