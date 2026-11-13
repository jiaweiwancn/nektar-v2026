#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H11 absolute producer selection, fixed continuation intent and zero-safe metrics."""
import math
import hashlib
import xml.etree.ElementTree as ET
import numpy as np
from wj_h7_analysis import context
from wj_h4_points import aligned,vector_hash

MATRIX=(('coarse',.01,3,4),('fine',.005,7,8),('quarter',.0025,15,16))
VIEWS={'coarse':(2,3,4),'fine':(2,4,6,8),'quarter':(4,8,12,16)}
FACE18=.31946961248954786
EXEC_INFO={'RANSWallRestart','RANSWallCheckpoint','RANSWallPimpleSnapshotDirectory','RANSWallPimpleArithmetic','RANSWallPimpleArithmeticTime'}

def semantics(xml):
    x=ET.parse(xml).getroot()if not isinstance(xml,ET.Element)else xml
    params={p.text.split('=')[0].strip():p.text.split('=')[1].strip()for p in x.findall('./CONDITIONS/PARAMETERS/P')}
    infos={p.get('PROPERTY'):p.get('VALUE')for p in x.findall('./CONDITIONS/SOLVERINFO/I')}
    return dict(parameters={k:v for k,v in params.items()if k!='NumSteps'},solverinfo={k:v for k,v in infos.items()if k not in EXEC_INFO},
        geometry=hashlib.sha256(ET.tostring(x.find('GEOMETRY'))).hexdigest(),expansions=hashlib.sha256(ET.tostring(x.find('EXPANSIONS'))).hexdigest(),
        functions=[hashlib.sha256(ET.tostring(f)).hexdigest()for f in x.findall('./CONDITIONS/FUNCTION')])

def admit_xml(path,branch,original):
    _,dt,n,target=next(v for v in MATRIX if v[0]==branch);x=ET.parse(path).getroot()
    params={p.text.split('=')[0].strip():p.text.split('=')[1].strip()for p in x.findall('./CONDITIONS/PARAMETERS/P')}
    assert semantics(path)==semantics(original),'physical XML intent changed'
    assert float(params['Time'])==0 and float(params['TimeStep'])==dt
    assert float(params['NumSteps']).is_integer()and int(float(params['NumSteps']))==n,'NumSteps must be additional 3/7/15'
    assert params['IO_CheckSteps']=='1'
    return n

def selected(r,step):
    first=r['steps']+1-len(r['accepted_states']);local=step-first-1
    if not 0<=local<len(r['pimple_attempts']):raise ValueError('restored/rejected/nonproduced step is not a completed attempt')
    c=context(r,step)
    if c['step']!=step:raise ValueError('requested absolute step differs from recorded callback time')
    aligned(step,r['dt'],c['time'])
    if c['alpha']!=1/r['dt']:raise ValueError('wrong completed callback alpha')
    a=r['pimple_attempts'][local]
    return c,a,local

def checkpoint_match(payload,r,step):
    c=context(r,step)
    return(payload['steps']==step and payload['dt']==r['dt']and payload['signature']==r['signature']and
        vector_hash(payload['state'])==vector_hash(c['state'])and abs(payload['time']-c['time'])<=128*np.finfo(float).eps*max(1,abs(c['time'])))

def ratio(d1,d2,floor):
    if not all(math.isfinite(v)for v in (d1,d2,floor))or min(d1,d2,floor)<0:raise ValueError('finite nonnegative comparison required')
    value=dict(D1=d1,D2=d2,floor=floor)
    if d1<=floor:return dict(value,rho=None,status='UNRESOLVED',reason='denominator at reconstruction floor')
    if d2<=floor:return dict(value,rho=None,rho_upper_bound=floor/d1,status='INDISTINGUISHABLE',reason='D2 below reconstruction accuracy; no precise zero/order')
    return dict(value,rho=d2/d1,status='CONTRACTS'if d2<d1 else 'DOES_NOT_CONTRACT')

def native_reserve(book,key,units):
    assert not book.get('closed')and not book.get('coverage_halt')and not book.get('h11_halt'),'monotone halt/closed allocation'
    assert key not in book.setdefault('native_reservations',[]),'native reservation cannot be relaunched'
    assert units in (0,3,7,15)and book.get('reserved_attempts',0)+units<=25
    book['native_reservations'].append(key);book['reserved_attempts']=book.get('reserved_attempts',0)+units

def common_step(dt,t):
    step=round(t/dt);aligned(step,dt,t);return step


def callback_relation(history,predecessor):
    """A GLM input is derived arithmetic, not the serialized restored slot.

    Preserve both exact hashes. Only opposite signed zeros may differ here;
    checkpoint/native-slot equality elsewhere remains strictly bitwise.
    """
    a=np.asarray(history,dtype='<f8');b=np.asarray(predecessor,dtype='<f8')
    assert a.shape==b.shape and np.all(np.isfinite(a))and np.all(np.isfinite(b))
    changed=np.flatnonzero(a.view('<u8')!=b.view('<u8'))
    assert np.array_equal(a,b)and np.all(a[changed]==0)and np.all(b[changed]==0),'nonzero physical callback history mismatch'
    return dict(history_sha256=vector_hash(a),predecessor_sha256=vector_hash(b),bitwise_equal=len(changed)==0,
        signed_zero_only_indices=changed.tolist(),numerically_equal=True,scope='actual derived GLM callback retained verbatim; no assertion of bitwise equality with the predecessor or repair of saved values')
