#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Fixed-case immutable rows and a durable, independently qualifying observer."""
import argparse
import copy
from collections import OrderedDict
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import pickle
import struct
import sys
import time
from run_wj_h4 import THREADS
os.environ.update(THREADS)  # before NumPy/SciPy
import numpy as np
from run_wj_h4 import load
from wj_h14_execution import OwnerIndex
from wj_h7_analysis import Rows,classify
from wj_h4_points import vector_hash,physical,aligned
from study_wj_h14 import stored_reference
from study_wj_h9 import serial
from wj3_measure import checkpoint_payload
sys.modules.setdefault('wj_h15_observer',sys.modules[__name__])

H14=Path('/home/jay/build/wj-h14-mapped-backend-20261113')
RULE='H14 all5x5/all owners+fixed physical observations+all retained/native witness1e-6; physical Gauss8 support; current-state recomputation; H15 immutable cache v1'

def bits(v):return struct.pack('>d',float(v)).hex()
def identity(p):
    p=Path(p);h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1024**2),b''):h.update(b)
    return dict(path=str(p),bytes=p.stat().st_size,sha256=h.hexdigest())
def verify(v):
    assert identity(v['path'])==v,'changed pinned input: '+v['path'];return Path(v['path'])
def atomic(p,v,immutable=True):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=Path(str(p)+'.tmp-'+str(os.getpid()))
    with tmp.open('x')as f:json.dump(serial(v),f,allow_nan=False,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
    if immutable:os.link(tmp,p);tmp.unlink()
    else:os.replace(tmp,p)
    fd=os.open(p.parent,os.O_RDONLY|os.O_DIRECTORY)
    try:os.fsync(fd)
    finally:os.close(fd)
def pkey(p):return tuple(float(v).hex()for v in p)
def okey(p,o):return pkey(p)+(o['geometry'],o['normal_side'],o['x_side'],o['interface_role'])
def source_identities():
    names=('wj_h15_observer','wj_h7_analysis','wj_h14_execution','wj_h4_points','wj_h3_analysis','wj_h5_neutral','wj_h2_profiles','validate_cls_epawf_m3','study_wj_h14')
    return [identity(Path(importlib.import_module(n).__file__))for n in names]
def reference_identity(root,r):
    shape={k:r[k]for k in ('fields','offsets','map_parameters','geometry_signature','parameters','quadrature_orders')}
    return dict(qh0=vector_hash(r['state']),fixed_reference=vector_hash(r['fixed_reference']),fixed_defect=vector_hash(r['fixed_defect']),
        physical_signature_sha256=hashlib.sha256(r['signature'].encode()).hexdigest(),geometry_space_sha256=hashlib.sha256(json.dumps(shape,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
        profile=identity(H14/'science/target-resource-design-v1.json'),native_runtime_manifest=identity(root/'science/runtime-build-v1.json'),rule=RULE)

class BoundedIndex(OwnerIndex):
    def basis(self,degree,x):
        if len(self.basis_values)>=8192:self.basis_values.pop(next(iter(self.basis_values)))
        return super().basis(degree,x)
    def matching(self,owner):
        if len(self.matching_rows)>=2048:self.matching_rows.pop(next(iter(self.matching_rows)))
        return super().matching(owner)

def add_witnesses(points,packet,report,c):
    xmin=min(e['lower'][0]for e in report['fields'][0]['mesh']);xmax=max(e['upper'][0]for e in report['fields'][0]['mesh'])
    for w in packet['principal_witnesses'].values():
        p=w['position'];candidates=[p]
        for axis in (0,2):
            for sign in (-1,1):q=list(p);q[axis]+=sign*1e-6;candidates.append(q)
        for q in candidates:
            if xmin<=q[0]<=xmax and c.d0<=q[2]<=c.top:points.setdefault(pkey(q),list(q))

class PersistentRows:
    """A private hash-qualified cache; no evolved-state values are retained."""
    def __init__(self,report,points,qh0,reference_identity,variable_limit=512):
        self.report={k:report[k]for k in ('fields','offsets','map_parameters','geometry_signature','signature','parameters','quadrature_orders')}
        self.config=stored_reference();self.index=BoundedIndex(self.report);self.qh0=list(qh0)
        self.reference_identity=reference_identity;self.variable_limit=variable_limit;self.variable=OrderedDict();self.evictions=0;self.extensions=0
        began=time.monotonic();self.base=Rows(self.report,points,self.config,owner_lookup=self.index,geometry_cache=self.index)
        self.base.report=self.report;self.base.kmap=[]
        self.initial=self.base.evaluate(qh0,principal=False)
        self.known={pkey(p):[]for p in points}
        for j,(p,o)in enumerate(self.base.rows[:self.base.point_count]):self.known[pkey(p)].append(j)
        self.inverse={}
        self.build_seconds=time.monotonic()-began

    def __getstate__(self):
        # wall() carries local analytic callables. Store immutable numerical
        # data only; recreate callable wrappers from the same pinned H13 record
        # on load, without any fit, inverse, FE projection or basis rebuild.
        state=dict(self.__dict__);state['config']=None
        state['reference_parameters']={k:getattr(self.config,k)for k in ('H','h','z0_e','u_tau','nu','d0','top','d_nu')}
        state['base']=copy.copy(self.base);state['base'].config=None
        state['variable']=OrderedDict()
        for k,(rows,initial)in self.variable.items():
            rows=copy.copy(rows);rows.config=None;state['variable'][k]=(rows,initial)
        return state

    def __setstate__(self,state):
        values=state.pop('reference_parameters');self.__dict__.update(state);self.config=stored_reference()
        assert all(bits(getattr(self.config,k))==bits(v)for k,v in values.items()),'cached reference parameter bytes differ'
        self.base.config=self.config
        for rows,initial in self.variable.values():rows.config=self.config

    def _inverse(self,rows,b):
        # Only the immutable base is permanent. Variable-object inverses belong
        # to their LRU entry and are released with it, not an unbounded id map.
        if not hasattr(rows,'h15_inverse'):rows.h15_inverse={}
        if b not in rows.h15_inverse:
            inv=np.full(len(rows.rows),-1,dtype=int);inv[rows.maps[b][0]]=np.arange(len(rows.maps[b][0]));rows.h15_inverse[b]=inv
        return rows.h15_inverse[b]

    def view(self,points):
        refs=[];initials=[];extension_seconds=0.
        for p in points:
            key=pkey(p)
            if key in self.known:
                for j in self.known[key]:refs.append((self.base,j));initials.append((self.initial,j))
                continue
            if key not in self.variable:
                began=time.monotonic();r=Rows(self.report,[p],self.config,owner_lookup=self.index,geometry_cache=self.index,include_volume=False)
                r.kmap=[];v=r.evaluate(self.qh0,principal=False);self.variable[key]=(r,v);self.extensions+=1
                extension_seconds+=time.monotonic()-began
                while len(self.variable)>self.variable_limit:self.variable.popitem(last=False);self.evictions+=1
            r,v=self.variable.pop(key);self.variable[key]=(r,v)
            for j in range(len(r.rows)):refs.append((r,j));initials.append((v,j))
        count=len(refs)
        for j in range(self.base.point_count,len(self.base.rows)):refs.append((self.base,j));initials.append((self.initial,j))
        view=object.__new__(Rows);view.report=self.report;view.config=self.config;view.point_count=count
        view.rows=[r.rows[j]for r,j in refs];view.weights=np.array([r.weights[j]for r,j in refs]);view.outer=np.array([r.outer[j]for r,j in refs]);view.rphys=np.array([r.rphys[j]for r,j in refs]);view.maps={}
        for b in range(6):
            active=np.flatnonzero(view.outer if b in (3,5)else ~view.outer if b==4 else np.ones(len(refs),bool))
            nc=(self.report['fields'][b]['mesh'][0]['degree']+1)**2
            arr=[np.empty((len(active),nc),dtype=int if k==0 else float)for k in range(4)]
            groups={}
            for dest,j in enumerate(active):
                r,source=refs[j];key=id(r)
                if key not in groups:groups[key]=(r,[],[])
                groups[key][1].append(dest);groups[key][2].append(source)
            for r,dest,source in groups.values():
                indices=self._inverse(r,b)[source];assert np.all(indices>=0)
                for k in range(4):arr[k][dest]=r.maps[b][1][k][indices]
            view.maps[b]=(active,tuple(arr))
        view.kids=np.array([r.kids[j]for r,j in refs]);view.kphi=np.array([r.kphi[j]for r,j in refs]);view.kmap=[]
        view.exact={k:np.array([r.exact[k][j]for r,j in refs])for k in self.base.exact}
        init={k:np.array([v[k][j]for v,j in initials])for k in self.initial}
        return view,init,dict(extension_seconds=extension_seconds,extensions=self.extensions,variable_points=len(self.variable),variable_limit=self.variable_limit,
            evictions=self.evictions,basis_keys=len(self.index.basis_values),matching_rows=len(self.index.matching_rows))

def assess(rows,state,initial,step):
    d=rows.evaluate(state);s=rows.summary(d,initial);active=np.flatnonzero(rows.outer);i=active[np.argmin(d['margin'][active])]
    s.update(min_in_plane=float(d['margin'][i]),min_transverse=float(np.min(d['transverse'][active])),min_covariance=float(np.min(d['covariance'])),
        min_scalar=min(float(np.min(d['k'])),float(np.min(d['epsilon'][active]))),witness=rows.witness(d,i),quality_scope='all explicit owners AND physical Gauss8 support',
        divergence_scaled_max=float(np.max(abs(d['divergence'])/np.maximum(1.,np.linalg.norm(d['G'],axis=(1,2))))))
    s['status']=classify(s['min_in_plane'],s['min_transverse'],s['min_covariance'],s['min_scalar']);groups={}
    for n,(p,o)in enumerate(rows.rows[:rows.point_count]):
        if o['interface_role'].startswith('matching_'):
            key=(float(p[0]).hex(),float(o['lower'][0]).hex(),float(o['upper'][0]).hex());groups.setdefault(key,{})[o['region']]=n
    interfaces=[];c=rows.config
    for pair in groups.values():
        if set(pair)!={'inner','outer'}:continue
        a,b=pair['inner'],pair['outer'];K=float(d['K'][b]);target=(K/(81/(2*math.sqrt(143))))**1.5/(c.H*.3819178635467331)
        interfaces.append(dict(physical=rows.rows[b][0],inner_owner=rows.rows[a][1],outer_owner=rows.rows[b][1],K=K,rc=float(d['rc'][b]),Lambda=float(d['Lambda'][b]),
            shared_U=(d['U'][a]-d['U'][b]).tolist(),shared_k=float(d['k'][a]-d['k'][b]),epsilon=float(d['epsilon'][b]),epsilon_target=target,
            epsilon_trace_defect=float(d['epsilon'][b]-target),total_traction_jump=(d['T'][a,:,2]-d['T'][b,:,2]).tolist(),
            k_flux_jump=float(d['nu_d'][a]*d['grad_k'][a,2]-d['nu_d'][b]*d['grad_k'][b,2]),reduced_pressure_jump=float(d['pi'][a]-d['pi'][b])))
    s['owner_interfaces']=interfaces;moments={}
    for b in (4,5):
        ii,(ids,phi,_,_)=rows.maps[b];rhs=np.zeros(rows.report['fields'][b]['size']);np.add.at(rhs,(ids-rows.report['offsets'][b]).ravel(),(phi*(rows.weights[ii]*d['divergence'][ii])[:,None]).ravel())
        moments[str(b)]=dict(moments=rhs.tolist(),maximum=float(np.max(abs(rhs))),signed_integral=float(np.sum(rhs)))
    if max(v['maximum']for v in moments.values())>2e-15:s['status']='OBSTRUCTED'
    if step==0 and s['divergence_scaled_max']>1e-12:s['status']='OBSTRUCTED'
    return serial(dict(summary=s,all_regional_Q1_moments=moments)),d

def parity(root):
    out=root/'science/observer-cache';out.mkdir();bindings=load(H14/'publication-v2/state-bindings.json')['states'];reports={};points=OrderedDict();expected={};inventories={}
    for v in bindings:
        name=v['case'];reports[name]=load(verify(v['producer']));a=load(H14/f'science/target/{name}'/('assessment-v2'if name=='cp0'else'assessment')/'result.json')
        inv=load(verify(a['inventory']));inventories[name]=inv;expected[name]=a
        for p in inv['points']:points.setdefault(pkey(p['physical']),p['physical'])
    r=reports['cp0'];assert len(r['state'])==24688
    reference=reference_identity(root,r)
    cache=PersistentRows(r,list(points.values()),r['state'],reference);decisions=[]
    for name in ('cp0','pilot1'):
        rr,init,stats=cache.view([v['physical']for v in inventories[name]['points']]);new,_=assess(rr,reports[name]['state'],init,reports[name]['steps']);old=expected[name]['states'][0]
        # Numerical fields and owners are unchanged operations, not a looser
        # favorable-status comparison. Preserve the complete first mismatch.
        ok=new['summary']==old['summary']and new['all_regional_Q1_moments']==old['all_regional_Q1_moments']
        record=dict(case=name,status='PASS'if ok else'UNRESOLVED',new=new,original=old,inventory=expected[name]['inventory'],cache=stats)
        atomic(out/(name+'-parity.json'),record);assert ok,'persistent/original H14 numerical assessment differs: '+name
        decisions.append(dict(case=name,status='PASS',owner_rows=rr.point_count,physical_rows=len(rr.rows),source_report=bindings[0 if name=='cp0'else 1]['producer']))
    path=out/'immutable-cache.pickle'
    with path.open('xb')as f:pickle.dump(cache,f,protocol=5);f.flush();os.fsync(f.fileno())
    atomic(out/'manifest.json',dict(status='PASS',cache=identity(path),reference_identity=reference,source_files=source_identities(),
        baseline_points=list(points.values()),parity=decisions,build_seconds=cache.build_seconds,scope='private verified immutable data only; current states never cached'))
    print('persistent cache parity PASS',cache.build_seconds,flush=True)

def load_cache(manifest):
    m=load(manifest);assert m['status']=='PASS'
    for v in m['source_files']:verify(v)
    verify(m['reference_identity']['native_runtime_manifest']);verify(m['reference_identity']['profile'])
    with verify(m['cache']).open('rb')as f:cache=pickle.load(f)
    assert cache.reference_identity==m['reference_identity'];return cache,m

def daemon(root,intent_path,cache_manifest,behavior='PASS'):
    wait=time.monotonic()
    while not intent_path.exists():
        assert time.monotonic()-wait<10,'missing immutable observer intent';time.sleep(.01)
    intent=load(intent_path);directory=Path(intent['directory']);directory.mkdir(parents=True,exist_ok=True)
    assert intent['observer_pid']==os.getpid(),'observer PID/intent differs'
    cache,manifest=load_cache(cache_manifest);required=OrderedDict((pkey(p),p)for p in manifest['baseline_points']);seen=set()
    atomic(directory/'observer-started.json',dict(pid=os.getpid(),nonce=intent['nonce'],cache=manifest['cache'],memory_scope='immutable rows/reference/qh0 only'))
    while not(directory/'observer-stop.json').exists():
        for ready_path in sorted(directory.glob('step-*.ready.json')):
            if ready_path.name in seen:continue
            seen.add(ready_path.name);ready=load(ready_path);binding=ready['binding'];step=int(binding['absolute_step']);began=time.monotonic()
            assert binding['nonce']==intent['nonce']and binding['source_sha256']==intent['source_sha256']and binding['observer_sha256']==intent['observer_sha256']
            report=load(verify(ready['report']));cp=checkpoint_payload(verify(ready['checkpoint']).parent,include_reference=True);side=load(verify(ready['callback']));structure=load(verify(ready['structure']))
            assert structure['fields']==cache.report['fields']and structure['offsets']==cache.report['offsets']and structure['geometry_signature']==cache.report['geometry_signature']
            for k in ('state','fixed_reference','fixed_defect'):assert vector_hash(cp[k])==vector_hash(report[k])
            assert vector_hash(cp['state'])==binding['state_sha256']and bits(cp['time'])==binding['time_bits']and cp['steps']==step
            aligned(step,report['dt'],cp['time']);assert cp['dt']==report['dt']
            assert vector_hash(side['callback_history'])==binding['callback_history_sha256']and bits(side['alpha_mass'])==binding['alpha_bits']
            assert binding['original_qh0_sha256']==cache.reference_identity['qh0']and binding['fixed_reference_sha256']==cache.reference_identity['fixed_reference']and binding['fixed_defect_sha256']==cache.reference_identity['fixed_defect']
            packet=report['wj_audit_packets'][0];assert packet['phase']=='completed_stage'and packet['state_sha256']==binding['state_sha256']
            assert packet['stage_history_sha256']==binding['callback_history_sha256']and bits(packet['time'])==binding['time_bits']and bits(packet['alpha_mass'])==binding['alpha_bits']
            assert report['native_history_times']==cp['native_history_times']and all(vector_hash(a)==vector_hash(b)for a,b in zip(report['native_history'],cp['native_history']))
            add_witnesses(required,packet,cache.report,cache.config)
            rows,initial,stats=cache.view(list(required.values()));result,data=assess(rows,report['state'],initial,step)
            replay=None
            if intent.get('replay_H14_step1')and step==1:
                original=load(verify(intent['replay_H14_step1']));replay={k:vector_hash(report[k])==vector_hash(original[k])for k in('state','fixed_reference','fixed_defect')}
                replay['native_history']=all(vector_hash(a)==vector_hash(b)for a,b in zip(report['native_history'],original['native_history']))and len(report['native_history'])==len(original['native_history'])
                replay['history_times']=report['native_history_times']==original['native_history_times'];assert all(replay.values()),'target replay parity failed'
            p=report['pimple_attempts'][0]
            assert all(p[k]==1 for k in('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls'))
            assert p['bound_activations']==p['additional_correction_or_restoration_solves']==0
            snapshot=load(p['k_snapshot']);assert identity(p['k_snapshot'])['sha256']==p['k_snapshot_sha256']and vector_hash(snapshot['raw'])==vector_hash(snapshot['corrected'])
            inventory=directory/(ready_path.stem+'.inventory.json');owned=[dict(physical=p,owner=o)for p,o in rows.rows[:rows.point_count]]
            atomic(inventory,dict(rule=RULE,points=list(required.values()),owned=owned,physical_rows=len(rows.rows),quadrature_rows=len(rows.rows)-rows.point_count))
            status=result['summary']['status'];status='OBSTRUCTED'if behavior=='FAIL'else status
            result.update(status=status,binding=binding,inventory_sha256=identity(inventory)['sha256'],replay=replay,cache=stats,
                measured_independent_seconds=time.monotonic()-began,reference_identity=cache.reference_identity,
                native_compute_seconds=ready['native_compute_seconds'],signed_frozen_residuals=p,
                completed_stage=report['stage_audits'],accepted_native=True,qualified=status=='PASS')
            dest=directory/(ready_path.stem+'.result.json');atomic(dest,result)
            if behavior=='timeout':continue
            cursor=directory/(ready_path.stem+'.qualified.json')
            if status=='PASS':atomic(cursor,dict(role='PREPARED_BEFORE_ACK; native qualification commit required for restart eligibility',binding=binding,ready_sha256=identity(ready_path)['sha256'],inventory_sha256=identity(inventory)['sha256'],checkpoint=ready['checkpoint'],callback=ready['callback'],report=ready['report']))
            final=step==intent['maximum_absolute_step'];continuation='TERMINAL'if final else'PROCEED'
            if not final and 'deadline_monotonic'in intent:
                remaining=intent['deadline_monotonic']-time.monotonic();native=max(80.,ready['native_compute_seconds']);assessment=max(10.,result['measured_independent_seconds'])
                forecast=1.25*((intent['maximum_absolute_step']-step)*(native+assessment)+intent.get('strict_reserve_seconds',120)+30)
                if remaining<forecast:continuation='RESOURCE_LIMIT'
                atomic(directory/(ready_path.stem+'.remaining.json'),dict(remaining=remaining,required_with_headroom=forecast,strict_reserved=intent.get('strict_reserve_seconds',120),continuation=continuation))
            ack=dict(schema='WJObserverAckV1',status=status,ready_sha256=identity(ready_path)['sha256'],binding=dict(binding),result=identity(dest),
                inventory=identity(inventory),qualified_cursor=identity(cursor)if status=='PASS'else identity(dest),continuation=continuation)
            if behavior=='stale':ack['binding']['nonce']='stale-'+intent['nonce']
            # Result and cursor have been fsynced before this atomic ACK.
            atomic(directory/ready_path.name.replace('.ready.json','.ack.json'),ack)
        time.sleep(.025)

if __name__=='__main__':
    PersistentRows.__module__='wj_h15_observer';BoundedIndex.__module__='wj_h15_observer'
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('parity','daemon'),required=True)
    p.add_argument('--intent',type=Path);p.add_argument('--cache',type=Path);p.add_argument('--behavior',default='PASS');a=p.parse_args()
    if a.mode=='parity':parity(a.root)
    else:daemon(a.root,a.intent,a.cache,a.behavior)
