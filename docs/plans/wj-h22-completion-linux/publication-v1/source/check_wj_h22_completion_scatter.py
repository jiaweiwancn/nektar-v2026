#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Passive whole-field scatter check, including rows outside the native free map."""
import argparse,os
from pathlib import Path
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
import numpy as np
from wj_h15_observer import identity,verify,atomic
from wj_h17_observer import native_load
from wj_h22_replay import read_completed,flow_vectors
from wj_h12_phase import identity as algebra_identity

def complete_scatter(states,increments,offsets):
    result={}
    for name,sl,inc,scale in (('velocity',slice(0,offsets[2]),'total',1.),
        ('pressure',slice(offsets[4],offsets[6]),'pressure_increment',.051**2)):
        before=states['C'][sl];after=states['F'][sl];delta=increments[inc][sl]
        result[name]=algebra_identity(after-before-delta,[after,before,delta],scale)
    return result

def selftest():
    states=dict(C=np.array([1.,-0.,0.,0.,.009,.009,.002,.002,.001,.001]),F=np.array([1.2,-0.,.2,0.,.009,.009,.002,.002,.0011,.0008]))
    inc=dict(total=np.array([.2,0.,.2,0.,0.,0.,0.,0.,0.,0.]),pressure_increment=np.array([0.,0.,0.,0.,0.,0.,0.,0.,.0001,-.0002]))
    complete_scatter(states,inc,[0,2,4,6,8,9,10])
    bad={k:v.copy()for k,v in states.items()};bad['F'][1]=1e-9
    try:complete_scatter(bad,inc,[0,2,4,6,8,9,10])
    except ValueError:return
    raise AssertionError('an omitted essential-row change passed the full scatter check')

def run(root):
    selftest();selected=load(root/'inputs/selected.json');index=load(root/'science/composite-source-index.json')
    offsets=load(verify(selected['sparse_structure']))['offsets'];out=root/'science/whole-field-scatter-v1';out.mkdir()
    arrays={};rows=[]
    for level,attempts in index['attempts'].items():
        for step,desc in attempts.items():
            r,_=read_completed(desc['anchor'],level,int(step));p=native_load(verify(desc['phase_snapshot']))
            states=dict(C=np.asarray(p['committed']),F=np.asarray(p['post_flow']))
            increments,proof=flow_vectors(desc['vector_packet'],states,r,offsets,desc['execution_identity'])
            checks=complete_scatter(states,increments,offsets);summaries={}
            for name,c in checks.items():
                for key in ('residual','bound'):arrays[level+step+'_'+name+'_'+key]=np.asarray(c[key])
                summaries[name]={k:c[k]for k in ('status','maximum_error','maximum_bound')}
            vectors=native_load(verify(desc['vector_packet']));raw=verify(vectors['payload']).read_bytes();zero={}
            for v in vectors['arrays']:
                values=np.frombuffer(raw[int(v['byte_offset']):int(v['byte_offset']+v['bytes'])],dtype='<f8')
                ids=np.flatnonzero((values==0)&np.signbit(values))
                zero[v['name']]=dict(negative_zero_count=len(ids),first16_indices=ids[:16].tolist(),all_bits_preserved_by=v['sha256'])
            rows.append(dict(level=level,step=int(step),source=desc,checks=summaries,zero_bits=zero,
                all_velocity_rows=offsets[2],all_pressure_rows=offsets[6]-offsets[4]))
    assert len(rows)==7
    np.savez_compressed(out/'componentwise-residuals-and-bounds.npz',**arrays)
    atomic(out/'result.json',dict(status='PASS',attempts=rows,selftest='positive full scatter and rejection of omitted-row change PASS',
        componentwise_arrays=identity(out/'componentwise-residuals-and-bounds.npz'),source=identity(__file__),
        governing_or_matrix_operations=0,scope='whole coefficient blocks, not only free rows; raw vector bits remain in original payloads'))
    print('All seven complete velocity/pressure scatter identities PASS, including nonfree rows',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();run(a.root)
