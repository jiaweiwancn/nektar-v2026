#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Measure immutable Python row preparation; every state dot/physics is fresh."""
import argparse,json,os,resource,sys,time
from pathlib import Path
from functools import lru_cache
from rans_qualification import THREADS
os.environ.update(THREADS)
import numpy as np
import ct_accepted8_native as n
from ct_accepted5_traction import profiles
from ct_accepted10_compare import gradients,outlet
from ct_accepted3_metrics import inventory
from ct_accepted_feasibility import write_json
from ct_accepted_open import require,row_value
from rans_request_manifest import RequestManifest,sha
from ct_accepted2_audit import values as original_values

class PreparedRows:
    """Scope-owned immutable parsed arrays. No state, value or validity cache.

    Callers own one hash-bound immutable observation packet. Strong references
    prevent Python object-ID recycling. Exact native ids/phi/derivatives and
    integer field offset are copied once; every actual input state is indexed
    and contracted again in the identical NumPy arithmetic/order.
    """
    policy='CT15-PreparedNativeRowsV1'
    def __init__(self):self.rows={};self.saved=[];self.calls=0;self.original=row_value
    def __call__(self,row,state,offset):
        key=(id(row),int(offset));cached=self.rows.get(key)
        if cached is None:
            ids=np.asarray(row['ids'],int)+offset;phi=np.asarray(row['phi']);d=np.asarray(row['derivative'])
            for v in (ids,phi,d):v.flags.writeable=False
            cached=(row,ids,phi,d);self.rows[key]=cached
        _,ids,phi,d=cached;self.calls+=1
        return phi@state[ids],d@state[ids]
    def __enter__(self):
        for module in list(sys.modules.values()):
            if module and getattr(module,'__name__','').startswith('ct_accepted') and getattr(module,'row_value',None) is self.original:
                self.saved.append(module);module.row_value=self
        return self
    def __exit__(self,*args):
        for module in self.saved:module.row_value=self.original
    @property
    def bytes(self):return sum(v.nbytes for _,*arrays in self.rows.values() for v in arrays)

class PreparedGrid(PreparedRows):
    """Batch the same complete native rows; qualify changed reduction order."""
    policy='CT15-PreparedNativeRowsV2'
    def __init__(self):super().__init__();self.grids={};self.saved_values=[]
    def values(self,initial,state,count):
        key=(id(initial),count);entry=self.grids.get(key);o=initial['observation'];off=o['offsets']
        if entry is None:
            fields=[]
            for field,name in enumerate(('velocity_k','velocity_k','velocity_k','epsilon','pressure')):
                groups={}
                for i,q in enumerate(initial['test_basis'][:count]):
                    row=q[name]
                    if row is None:continue
                    block=field if field<4 else 4+o['grid_regions'][i]
                    groups.setdefault(len(row['ids']),[]).append((i,row,off[block]))
                batches=[]
                for group in groups.values():
                    indices=np.asarray([v[0] for v in group],int)
                    ids=np.asarray([[j+v[2] for j in v[1]['ids']] for v in group],int)
                    phi=np.asarray([v[1]['phi'] for v in group]);derivative=np.asarray([v[1]['derivative'] for v in group])
                    for v in (indices,ids,phi,derivative):v.flags.writeable=False
                    batches.append((indices,ids,phi,derivative))
                fields.append(batches)
            entry=(initial,fields);self.grids[key]=entry
        s=np.asarray(state);v=np.zeros((count,5));g=np.zeros((count,5,3))
        for field,batches in enumerate(entry[1]):
            for indices,ids,phi,d in batches:
                coefficients=s[ids]
                v[indices,field]=np.einsum('ij,ij->i',phi,coefficients,optimize=False)
                g[indices,field,:]=np.einsum('iaj,ij->ia',d,coefficients,optimize=False)
        return dict(observation=dict(time=0.,grid=o['grid'][:count],grid_sides=o['grid_sides'][:count],grid_regions=o['grid_regions'][:count],grid_values=v.tolist(),grid_gradients=g.tolist()))
    def __enter__(self):
        super().__enter__()
        for module in list(sys.modules.values()):
            if module and getattr(module,'__name__','').startswith('ct_accepted') and getattr(module,'values',None) is original_values:
                self.saved_values.append(module);module.values=self.values
        return self
    def __exit__(self,*args):
        for module in self.saved_values:module.values=original_values
        super().__exit__(*args)

def maximum_difference(a,b):
    if isinstance(a,np.ndarray):a=a.tolist()
    if isinstance(b,np.ndarray):b=b.tolist()
    if isinstance(a,dict):
        require(a.keys()==b.keys(),'prepared observation fields changed');return max((maximum_difference(v,b[k]) for k,v in a.items()),default=0.)
    if isinstance(a,list):
        require(len(a)==len(b),'prepared observation membership changed');return max((maximum_difference(x,y) for x,y in zip(a,b,strict=True)),default=0.)
    if isinstance(a,(bool,str)) or a is None:require(a==b,'prepared scientific decision changed');return 0.
    return abs(float(a)-float(b))

def main():
    p=argparse.ArgumentParser();p.add_argument('--ct14',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));n.reference=lru_cache(None)(n.reference)
    pair=n.read(a.ct14/'pair-v2/results.json');c=pair['case'];ini=n.read(a.ct14/'static-v1/target/frozen.json');weights=n.inventory(c)[2];count=len(weights)
    m=RequestManifest(a.output,dict(schema='CT15IndependentRowPreparationV2',source=sha(Path(__file__)),oracle_pair_sha256=sha(a.ct14/'pair-v2/results.json'),initial_packet_sha256=sha(a.ct14/'static-v1/target/frozen.json'),policies=[PreparedRows.policy,PreparedGrid.policy],scaled_action_tolerance=2e-10,native_advance=False))
    def evaluate(state):return dict(metrics=n.physical_metrics(ini,state,c,weights,count),profiles=profiles(ini,state,c),gradients=gradients(ini,state,c,weights,count),outlet=outlet(ini,state,c),inventory=inventory(ini['test_basis'],ini['observation']['grid'],np.asarray(state),ini['observation']['offsets'],'A95',ini['observation']['grid_sides']))
    out=dict(schema='CT15IndependentRowPreparationV2',status='incomplete',members={})
    for dt in ('0.00025','0.000125'):
        state=n.read(a.ct14/'pair-v2'/dt/'case.evolution.json')['state'];cpu=time.process_time();wall=time.monotonic();before=evaluate(state);old=dict(CPU=time.process_time()-cpu,wall=time.monotonic()-wall)
        cpu=time.process_time();wall=time.monotonic()
        with PreparedRows() as rows:after=evaluate(state)
        new=dict(CPU=time.process_time()-cpu,wall=time.monotonic()-wall,unique_rows=len(rows.rows),array_bytes=rows.bytes,calls=rows.calls)
        require(before==after,'prepared row arrays changed independent scientific observations/decisions')
        cpu=time.process_time();wall=time.monotonic()
        with PreparedGrid() as rows:batched=evaluate(state)
        batch_time=dict(CPU=time.process_time()-cpu,wall=time.monotonic()-wall,prepared_grids=len(rows.grids))
        difference=maximum_difference(before,batched)/min(c['ustar']**2,c['K0'],n.col(c).A/(c['H']-c['displacement']))
        require(difference<=2e-10,'batched observation action differs beyond scaled roundoff bound')
        out['members'][dt]=dict(status='pass',baseline=old,prepared=new,batched=batch_time,maximum_scaled_batched_difference=difference,exact_prepared_results_equal=True,all_decisions_equal=True)
        write_json(a.output/'progress.json',out);print('CT15 independent rows',dt,old,new,batch_time,'scaled difference',difference,flush=True)
    out['status']='pass';write_json(a.output/'results.json',out);m.seal(a.output,['results.json'])
if __name__=='__main__':main()
