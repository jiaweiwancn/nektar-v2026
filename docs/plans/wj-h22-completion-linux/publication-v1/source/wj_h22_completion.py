# SPDX-License-Identifier: MIT
"""Role/identity adapters only; numerical evaluation remains the H22 evaluator."""
from pathlib import Path
from run_wj_h4 import load
from wj_h15_observer import identity,verify,bits,vector_hash
from wj_h17_observer import native_load,equal
from wj_h22_replay import DT,END,integer,read_completed,parity
from wj_h22_audit_roles import native_full_map_parity

def validate_role(intent,numsteps):
    assert intent['schema']=='WJMonitorIntentV1'and 0<intent['ack_timeout_seconds']<=120
    assert type(numsteps)is int
    assert all(type(intent[k])is int for k in ('maximum_absolute_step','restored_absolute_step','scientific_stride','recovery_stride'))
    assert all(type(v)is int for v in intent['substep_evidence_steps'])
    assert intent['recovery_stride']==1
    if intent['completion_role']=='RetainedM2StartupV1':
        assert intent['level']=='M'and numsteps==0
        assert intent['restored_absolute_step']==intent['maximum_absolute_step']==2
        assert intent['scientific_stride']==2 and intent['dt_bits']==bits(DT['M'])
        assert intent['substep_evidence_steps']==[]and intent['flow_capture']=='Off'
        c=native_load(verify(intent['restart_completion']))
        assert integer(c['binding']['absolute_step'])==2
        if 'restart_callback'in intent:assert intent['restart_callback']==c['callback'],'cannot override authenticated callback'
    else:
        assert intent['completion_role']=='MissingFReplayV1'and intent['level']=='F'and numsteps==4
        assert intent['restored_absolute_step']==0 and intent['maximum_absolute_step']==4
        assert intent['scientific_stride']==4 and intent['dt_bits']==bits(DT['F'])
        assert intent['substep_evidence_steps']==[1,2,3,4]and intent['flow_capture']=='MappedIncrementV1'
        assert 'restart_completion'not in intent and 'restart_callback'not in intent

def restored_parity(current,descriptor,historical,auxiliary,cp):
    """Current startup is not a completed producer; authenticate the old one first."""
    assert current['schema']=='WJMonitoringStartupV1'and current['steps']==2 and current['dt']==DT['M']
    retained,retained_cp=read_completed(descriptor,'M',2)
    old,old_cp=read_completed(historical,'M',2)
    first=parity(retained,old,retained_cp,old_cp)
    actual=parity(current,retained,cp,retained_cp)
    audited=native_load(verify(auxiliary['audit_report']))
    maps=native_full_map_parity(current,old,audited)
    return dict(status='PASS',retained_producer_parity=first,restored_numerical_parity=actual,full_map_parity=maps,
        genuine_completion=descriptor['completion'],genuine_callback=descriptor['callback'],startup_is_producer=False)

def typed_strict(directory,expected,level,step):
    """No fallback: a request for M2/F4 requires that exact qualified commit."""
    directory=Path(directory);cursor=directory/'latest-fully-verified-checkpoint.json'
    q=load(cursor);kind=q['schema']
    assert kind in ('WJQualificationCommitV1','WJStartupVerificationCommitV1')
    assert integer(q['binding']['absolute_step'])==step and step==END[level]
    if kind=='WJStartupVerificationCommitV1':
        assert level=='M'and q['completion']==expected['completion']
        commit=directory/'startup.verification-commit.json';completion=native_load(verify(q['completion']))
        for key in ('report','checkpoint','callback'):assert completion[key]==expected[key]
        assert q['checkpoint']==expected['checkpoint']
        v=load(verify(q['result']));assert v['status']=='PASS'and v['binding']==q['binding']
        result=load(verify(v['result']));assert result['status']=='PASS'and result['binding']==q['binding']
    else:
        assert level=='F';commit=directory/f'step-{step}.qualification-commit.json'
        for key in ('report','checkpoint','callback'):assert q[key]==expected[key]
    assert identity(commit)['sha256']==identity(cursor)['sha256']
    report,cp=read_completed(expected,level,step)
    assert q['binding']['state_sha256']==vector_hash(report['state'])and q['binding']['time_bits']==bits(report['time'])
    assert q['binding']['callback_history_sha256']==vector_hash(report['last_stage_history'])
    assert q['binding']['alpha_bits']==bits(1/DT[level])
    return dict(operation='open_checks',inner_stress_model='Linear',import_role='QualifiedV1',
        native_report=expected['report']['path'],native_report_sha256=expected['report']['sha256'],
        native_checkpoint=str(verify(expected['checkpoint']).parent),checkpoint_sha256=expected['checkpoint']['sha256'],
        completed_sidecar=expected['callback'],qualified_cursor=identity(commit),weak=True,weak_compact=True,
        assembly=False,reuse=False,retain_scalar_observations=False)

def capture_descriptor(anchor,execution,qualification=None):
    r=native_load(verify(anchor['report']));a=r['pimple_attempts'][0]
    snapshot=identity(a['k_snapshot']);assert snapshot['sha256']==a['k_snapshot_sha256']
    packet=native_load(verify(snapshot));matrix=identity(packet['matrix_file']);assert matrix['sha256']==packet['matrix_sha256']
    verify(a['flow_vectors']);verify(execution)
    return dict(anchor=anchor,execution_identity=execution,phase_snapshot=snapshot,sparse_k_matrix=matrix,
        vector_packet=a['flow_vectors'],qualification=qualification or{},producer_role='genuine original completed attempt')

def require_M2_gate(root):
    selected=load(root/'inputs/selected.json');folder=root/'science/M2'
    verified=load(folder/'verification.json');strict=load(folder/'final-strict/result.json')
    assert verified['status']==strict['status']=='PASS'
    request=typed_strict(folder/'monitor',selected['retained_attempts']['M']['2']['anchor'],'M',2)
    assert verified['new_startup_commit']==request['qualified_cursor']
    assert strict['selected_absolute_step']==2 and strict['genuine_alpha']==400
    assert load(verify(strict['request']))==request
    verify(strict['producer'])
    return dict(status='PASS',new_verification=identity(folder/'verification.json'),strict_read=identity(folder/'final-strict/result.json'))
