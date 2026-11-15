#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Own-step recovery reader and actual phase/vector binding; no governing algebra."""
import hashlib,math,re
import numpy as np
from wj_h15_observer import identity,verify,bits,vector_hash
from wj_h17_observer import native_load,equal
from wj_h7_analysis import context
from wj_h4_points import aligned
from wj3_measure import checkpoint_payload
from study_wj_h16 import same_cp
from wj_h9_local import same
from wj_h12_phase import identity as arithmetic_identity,scatter

DT=dict(C=.005,M=.0025,F=.00125)
END=dict(C=1,M=2,F=4)
MATRIX=tuple((z,n)for z in DT for n in range(1,END[z]+1))
PHASES=('P','C','F','K','A')
PAYLOAD=('predictor','bGz','du','dz')

def integer(value):
    if isinstance(value,str):assert re.fullmatch(r'0|[1-9][0-9]*',value)
    else:assert not isinstance(value,bool)and math.isfinite(value)and int(value)==value
    return int(value)

def selected(report,level,step):
    assert (level,step)in MATRIX and report['dt']==DT[level]
    c=context(report,step)
    assert c['step']==step and report['steps']==step
    aligned(step,DT[level],c['time'])
    assert c['alpha']==1/DT[level]and bits(c['time'])==bits(report['time'])
    assert vector_hash(c['state'])==vector_hash(report['state'])
    return c

def read_completed(anchor,level,step):
    r=native_load(verify(anchor['report']));c=selected(r,level,step)
    cp=checkpoint_payload(verify(anchor['checkpoint']).parent,include_reference=True);same_cp(r,cp)
    side=native_load(verify(anchor['callback']));completion=native_load(verify(anchor['completion']))
    assert side['schema']=='WJCompletedCallbackV1'and completion['schema']=='WJNativeCompletionV1'
    b=completion['binding'];assert integer(b['absolute_step'])==step
    for key in ('checkpoint','callback','report'):assert completion[key]==anchor[key]
    verify(completion['configuration'])
    assert r['run_configuration_sha256']==completion['configuration']['sha256']
    assert b['physical_signature_sha256']==hashlib.sha256(r['signature'].encode()).hexdigest()
    assert b['checkpoint_sha256']==anchor['checkpoint']['sha256']and b['report_sha256']==anchor['report']['sha256']
    assert b['callback_sidecar_sha256']==anchor['callback']['sha256']
    assert b['state_sha256']==side['state_sha256']==vector_hash(r['state'])
    assert b['callback_history_sha256']==side['callback_history_sha256']==vector_hash(c['history'])
    assert b['time_bits']==bits(c['time'])and b['alpha_bits']==bits(c['alpha'])
    assert bits(side['time'])==bits(c['time'])and side['alpha_mass']==c['alpha']
    equal(side['callback_history'],c['history'],'own genuine callback')
    assert b['fixed_reference_sha256']==vector_hash(r['fixed_reference'])and b['fixed_defect_sha256']==vector_hash(r['fixed_defect'])
    assert b['original_qh0_sha256']==r['original_qh0_state_sha256']
    return r,cp

def recovery(directory,step):
    completion=identity(directory/f'step-{step}.completion.json');v=native_load(verify(completion))
    assert v['schema']=='WJNativeCompletionV1'and integer(v['binding']['absolute_step'])==step
    return dict(**{k:v[k]for k in ('checkpoint','callback','report')},completion=completion)

def parity(actual,old,actual_cp,old_cp):
    for name in ('state','fixed_reference','fixed_defect','native_history','native_history_times','dt','time','steps','signature',
                 'last_stage_history','last_stage_time','last_stage_alpha','original_qh0_state_sha256'):
        equal(actual[name],old[name],'exact replay '+name)
    equal(actual_cp['eta'],old_cp['eta'],'eta')
    return dict(status='PASS',coupled_state_sha256=vector_hash(actual['state']),
                callback_history_sha256=vector_hash(actual['last_stage_history']),
                native_history_sha256=[vector_hash(v)for v in actual['native_history']],time_bits=bits(actual['time']))

def phases(report,packet,predecessor,offsets):
    states={p:np.asarray(v)for p,v in dict(P=predecessor['state'],C=packet['committed'],
         F=packet['post_flow'],K=packet['raw'],A=report['state']).items()}
    assert same(states['P'],states['C']),'stationary lift changes predecessor'
    assert same(packet['post_flow'],packet['pre_k'])and same(packet['raw'],packet['corrected'])
    equal(packet['fixed_reference'],report['fixed_reference'],'phase fixed reference')
    assert bits(packet['time'])==bits(report['time'])and packet['alpha']==report['last_stage_alpha']
    assert vector_hash(packet['history'])==vector_hash(report['last_stage_history'])
    assert packet['k_bound']['activations']==0 and packet['k_bound']['mass_norm_correction']==0
    masks={}
    for before,after,allowed in (('C','F',(0,1,4,5)),('F','K',(2,)),('K','A',(3,))):
        changed=[]
        for b in range(6):
            sl=slice(offsets[b],offsets[b+1]);diff=not same(states[before][sl],states[after][sl])
            assert not diff or b in allowed,'actual phase changes a forbidden field'
            if diff:changed.append(b)
        masks[before+'->'+after]=changed
    for b in (0,1,4,5):assert same(states['F'][offsets[b]:offsets[b+1]],states['A'][offsets[b]:offsets[b+1]])
    return states,dict(states={p:vector_hash(v)for p,v in states.items()},P_C_bitwise=True,post_flow_pre_k_bitwise=True,
         raw_corrected_bitwise=True,flow_A_U_pi_bitwise=True,field_masks=masks,
         provisional_epsilon_trace='K retains actual pre-epsilon field; no trace repair or accepted-state eligibility')

def flow_vectors(pin,states,report,offsets,expected_identity):
    packet=native_load(verify(pin));assert packet['schema']=='WJMappedFlowVectorsV1'
    assert packet['status']=='FLOW_PROVISIONAL_NOT_TIME_COMMIT'
    assert integer(packet['absolute_step'])==report['steps']
    assert packet['time_bits']==bits(report['time'])and packet['alpha_bits']==bits(report['last_stage_alpha'])
    assert packet['base_sha256']==vector_hash(states['C'])and packet['flow_state_sha256']==vector_hash(states['F'])
    assert packet['callback_history_sha256']==vector_hash(report['last_stage_history'])
    assert packet['execution_identity']==expected_identity;verify(expected_identity)
    assert packet['additional_governing_operations']==0
    assert packet['pressure_velocity_sign']=='actual bGz; signed correction is -bGz'
    payload=verify(packet['payload']).read_bytes();arrays={};at=0
    assert [v['name']for v in packet['arrays']]==list(PAYLOAD)
    for v in packet['arrays']:
        assert integer(v['byte_offset'])==at and integer(v['bytes'])==8*integer(v['count'])
        raw=payload[at:at+int(v['bytes'])];assert len(raw)==v['bytes']and hashlib.sha256(raw).hexdigest()==v['sha256']
        arrays[v['name']]=np.frombuffer(raw,dtype='<f8').copy();at+=len(raw)
    assert at==len(payload)
    assert all(isinstance(v,int)and not isinstance(v,bool)for v in packet['velocity_rows']+packet['pressure_rows'])
    vr=np.asarray(packet['velocity_rows'],int);pr=np.asarray(packet['pressure_rows'],int)
    assert len(set(vr)|set(pr))==len(vr)+len(pr)and len(set(vr))==len(vr)and len(set(pr))==len(pr)
    assert np.all((vr>=0)&(vr<offsets[2]))and np.all((pr>=offsets[4])&(pr<offsets[6]))
    assert all(len(arrays[n])==len(vr)for n in ('predictor','bGz','du'))and len(arrays['dz'])==len(pr)
    checks=dict(native_du=arithmetic_identity(arrays['du']-arrays['predictor']+arrays['bGz'],
                    [arrays['du'],arrays['predictor'],arrays['bGz']],1.),
        velocity_coefficients=arithmetic_identity(states['F'][vr]-states['C'][vr]-arrays['du'],
                    [states['F'][vr],states['C'][vr],arrays['du']],1.),
        pressure_coefficients=arithmetic_identity(states['F'][pr]-states['C'][pr]-arrays['dz'],
                    [states['F'][pr],states['C'][pr],arrays['dz']],.051**2))
    size=len(states['C'])
    increments=dict(predictor=scatter(vr,arrays['predictor'],size),pressure=scatter(vr,-arrays['bGz'],size),
                    total=scatter(vr,arrays['du'],size),pressure_increment=scatter(pr,arrays['dz'],size))
    return increments,dict(status='PASS',packet=pin,payload=packet['payload'],execution_identity=expected_identity,
           vector_hashes={n:vector_hash(v)for n,v in arrays.items()},velocity_rows=packet['velocity_rows'],pressure_rows=packet['pressure_rows'],
           checks=checks,signed_pressure_correction='-actual_saved_bGz',additional_governing_operations=0)
