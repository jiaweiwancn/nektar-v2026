# SPDX-License-Identifier: MIT
"""Fixed H23 continuation roles; unchanged H19/H20 numerical comparison helpers."""
import hashlib
from pathlib import Path
from run_wj_h4 import load
from wj_h15_observer import identity,verify,vector_hash,bits
from wj_h17_observer import native_load,equal
from wj_h4_points import aligned
from wj_h7_analysis import context
from wj_h22_replay import integer,recovery,parity
from wj3_measure import checkpoint_payload
from study_wj_h16 import same_cp

DT=dict(C=.005,M=.0025,F=.00125)
START=dict(C=4,M=8,F=16)
END=dict(C=8,M=16,F=32)
ADDITIONAL={z:END[z]-START[z]for z in DT}
STRIDE=dict(C=1,M=2,F=4)
TIMES=(.025,.030,.035,.040)
ORDER=('M','C','F')
GOVERNOR_CAP=dict(M=1050,C=900,F=1650)
OBSERVER_FORECAST=dict(M=15,C=25,F=45)

def absolute_step(level,t):
    assert level in DT
    matches=[v for v in (.02,*TIMES)if abs(t-v)<=128*2.220446049250313e-16*max(1.,abs(t),abs(v))]
    assert len(matches)==1,'not an H23 physical comparison time'
    n=round(matches[0]/DT[level]);aligned(n,DT[level],t);return n

def selected(report,level,step):
    step=integer(step)
    assert START[level]<=step<=END[level]and report['dt']==DT[level]
    c=context(report,step)
    assert integer(c['step'])==step==integer(report['steps'])
    aligned(step,DT[level],c['time'])
    assert bits(c['time'])==bits(report['time'])and c['alpha']==1/DT[level]
    assert vector_hash(c['state'])==vector_hash(report['state'])
    return c

def read_completed(anchor,level,step):
    """Authenticate native-completed bytes, independently of qualification role."""
    r=native_load(verify(anchor['report']));c=selected(r,level,step)
    cp=checkpoint_payload(verify(anchor['checkpoint']).parent,include_reference=True);same_cp(r,cp)
    side=native_load(verify(anchor['callback']));completion=native_load(verify(anchor['completion']))
    assert side['schema']=='WJCompletedCallbackV1'and completion['schema']=='WJNativeCompletionV1'
    b=completion['binding'];assert integer(b['absolute_step'])==step
    for k in ('checkpoint','callback','report'):assert completion[k]==anchor[k]
    verify(completion['configuration'])
    assert r['run_configuration_sha256']==completion['configuration']['sha256']
    assert b['physical_signature_sha256']==hashlib.sha256(r['signature'].encode()).hexdigest()
    assert b['checkpoint_sha256']==anchor['checkpoint']['sha256']and b['report_sha256']==anchor['report']['sha256']
    assert b['callback_sidecar_sha256']==anchor['callback']['sha256']
    assert b['state_sha256']==side['state_sha256']==vector_hash(r['state'])
    assert b['callback_history_sha256']==side['callback_history_sha256']==vector_hash(c['history'])
    assert b['time_bits']==bits(c['time'])and b['alpha_bits']==bits(c['alpha'])
    assert bits(side['time'])==bits(c['time'])and side['alpha_mass']==c['alpha']
    equal(side['callback_history'],c['history'],'own completed callback')
    for name in ('fixed_reference','fixed_defect'):assert b[name+'_sha256']==vector_hash(r[name])
    assert b['original_qh0_sha256']==r['original_qh0_state_sha256']
    return r,cp

def qualified(anchor,level,step,commit_pin=None):
    r,cp=read_completed(anchor,level,step);pin=commit_pin or anchor['qualified_commit'];q=load(verify(pin));b=q['binding']
    assert integer(b['absolute_step'])==step and b['time_bits']==bits(r['time'])
    assert b['state_sha256']==vector_hash(r['state'])and b['callback_history_sha256']==vector_hash(r['last_stage_history'])
    assert b['alpha_bits']==bits(1/DT[level])
    if q['schema']=='WJQualificationCommitV1':
        for k in ('report','checkpoint','callback'):assert q[k]==anchor[k]
    else:
        assert q['schema']=='WJStartupVerificationCommitV1'and q['completion']==anchor['completion']
        assert q['checkpoint']==anchor['checkpoint']
        wrapper=load(verify(q['result']));assert wrapper['status']=='PASS'and wrapper['binding']==b
        result=load(verify(wrapper['result']));assert result['status']=='PASS'and result['binding']==b
    return r,cp

def partners(level):
    return dict(M=(),C=('M',),F=('C','M'))[level]

def resolve(root,inputs,level,t):
    n=absolute_step(level,t)
    if n==START[level]:return inputs['anchors'][level]
    directory=Path(root)/'science'/level/'monitor';p=directory/f'step-{n}.qualification-commit.json'
    assert p.exists(),'required computed H23 partner missing: '+level+str(n)
    a=recovery(directory,n);a['qualified_commit']=identity(p);return a

def validate_intent(v,numsteps):
    z=v['level'];assert z in DT and v['schema']=='WJMonitorIntentV1'and v['h23_role']=='ContinuationFrom02V1'
    assert type(numsteps)is int and numsteps==ADDITIONAL[z]
    for k in ('maximum_absolute_step','restored_absolute_step','recovery_stride','scientific_stride'):assert type(v[k])is int
    assert v['maximum_absolute_step']==END[z]and v['restored_absolute_step']==START[z]
    assert v['recovery_stride']==1 and v['scientific_stride']==STRIDE[z]and v['dt_bits']==bits(DT[z])
    assert 0<v['ack_timeout_seconds']<=120 and v['substep_evidence_steps']==[]
    assert v['flow_capture']=='Off'and v['snapshot_capture']=='Off'
    c=native_load(verify(v['restart_completion']));assert integer(c['binding']['absolute_step'])==START[z]
    if 'restart_callback'in v:assert v['restart_callback']==c['callback']

def local_forecast(level,step,ordinary=45.,audit=50.,observer=None):
    observer=OBSERVER_FORECAST[level]if observer is None else observer
    due=sum(n%STRIDE[level]==0 for n in range(step+1,END[level]+1))
    return (END[level]-step)*ordinary+due*(audit+observer)+30

def branch_forecast(level,ordinary=45.,audit=50.,observer=None):
    return 150+local_forecast(level,START[level],ordinary,audit,observer)+120

def deadline_ok(deadline,now,level,step,**costs):
    # A later strict reader is outside this governor's own operation deadline.
    return deadline-now>=1.25*local_forecast(level,step,**costs)

def progress_forecast(book,level,step,elapsed,future_costs=None,**costs):
    """Actual completed spend plus current elapsed plus future work, once each."""
    entries=[e for e in book['entries'].values()if e.get('kind')!=level]
    native=sum(e.get('seconds',e['reserved_seconds'])for e in entries if e['phase']=='native')+elapsed
    total=60+300+sum(e.get('seconds',e['reserved_seconds'])for e in entries)+elapsed
    later=ORDER[ORDER.index(level)+1:]
    remaining=local_forecast(level,step,**costs)+120+sum(branch_forecast(z,**(future_costs[z]if future_costs else{}))for z in later)
    return dict(status='PASS'if native+1.25*remaining<=4200 and total+1.25*remaining<=6000 else'RESOURCE_LIMIT',
        actual_native=native,actual_total_with_closure=total,remaining_estimate=remaining,headroom=1.25,
        protected_strict_seconds=120*(1+len(later)),closure_seconds=300,active_reservation_counted_as_spent=False)

def strict_request(directory,original,level,expect_terminal=True):
    """Typed qualified selection; startup uses the original completed producer."""
    directory=Path(directory);cursor=directory/'latest-fully-verified-checkpoint.json'
    if cursor.exists():
        q=load(cursor);step=integer(q['binding']['absolute_step'])
        if q['schema']=='WJQualificationCommitV1':
            anchor=recovery(directory,step);pin=identity(directory/f'step-{step}.qualification-commit.json')
        else:
            assert q['schema']=='WJStartupVerificationCommitV1';anchor=original;pin=identity(directory/'startup.verification-commit.json')
        assert identity(cursor)['sha256']==pin['sha256']
    else:
        anchor=original;step=START[level];pin=original['qualified_commit']
    if expect_terminal:assert step==END[level],'requested .04 endpoint is not qualified; no fallback'
    r,cp=qualified(anchor,level,step,pin)
    request=dict(operation='open_checks',inner_stress_model='Linear',import_role='QualifiedV1',
        native_report=anchor['report']['path'],native_report_sha256=anchor['report']['sha256'],
        native_checkpoint=str(verify(anchor['checkpoint']).parent),checkpoint_sha256=anchor['checkpoint']['sha256'],
        completed_sidecar=anchor['callback'],qualified_cursor=pin,weak=True,weak_compact=True,
        assembly=False,reuse=False,retain_scalar_observations=False)
    return request,dict(step=step,time=r['time'],alpha=r['last_stage_alpha'],anchor=anchor,commit=pin)
