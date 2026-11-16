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

def main():
    p=argparse.ArgumentParser();p.add_argument('--ct14',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));n.reference=lru_cache(None)(n.reference)
    pair=n.read(a.ct14/'pair-v2/results.json');c=pair['case'];ini=n.read(a.ct14/'static-v1/target/frozen.json');weights=n.inventory(c)[2];count=len(weights)
    m=RequestManifest(a.output,dict(schema='CT15IndependentRowPreparationV1',source=sha(Path(__file__)),oracle_pair_sha256=sha(a.ct14/'pair-v2/results.json'),initial_packet_sha256=sha(a.ct14/'static-v1/target/frozen.json'),policy=PreparedRows.policy,native_advance=False))
    def evaluate(state):return dict(metrics=n.physical_metrics(ini,state,c,weights,count),profiles=profiles(ini,state,c),gradients=gradients(ini,state,c,weights,count),outlet=outlet(ini,state,c),inventory=inventory(ini['test_basis'],ini['observation']['grid'],np.asarray(state),ini['observation']['offsets'],'A95',ini['observation']['grid_sides']))
    out=dict(schema='CT15IndependentRowPreparationV1',status='incomplete',members={})
    for dt in ('0.00025','0.000125'):
        state=n.read(a.ct14/'pair-v2'/dt/'case.evolution.json')['state'];cpu=time.process_time();wall=time.monotonic();before=evaluate(state);old=dict(CPU=time.process_time()-cpu,wall=time.monotonic()-wall)
        cpu=time.process_time();wall=time.monotonic()
        with PreparedRows() as rows:after=evaluate(state)
        new=dict(CPU=time.process_time()-cpu,wall=time.monotonic()-wall,unique_rows=len(rows.rows),array_bytes=rows.bytes,calls=rows.calls)
        require(before==after,'prepared row arrays changed independent scientific observations/decisions')
        out['members'][dt]=dict(status='pass',baseline=old,prepared=new,exact_results_equal=True)
        write_json(a.output/'progress.json',out);print('CT15 independent rows',dt,old,new,flush=True)
    out['status']='pass';write_json(a.output/'results.json',out);m.seal(a.output,['results.json'])
if __name__=='__main__':main()
