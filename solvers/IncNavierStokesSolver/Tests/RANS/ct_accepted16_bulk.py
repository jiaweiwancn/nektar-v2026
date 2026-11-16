#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded immutable batched rows; all state values and physics remain fresh."""
import argparse,hashlib,json,os,resource,sys,time
from pathlib import Path
from functools import lru_cache
from rans_qualification import THREADS
os.environ.update(THREADS)
import numpy as np
import ct_accepted8_native as n
from ct_accepted15_python import PreparedRows,maximum_difference
from ct_accepted15_audits import bundle
from ct_accepted2_audit import values as original_values
from ct_accepted_open import require
from ct_accepted_feasibility import write_json
from rans_trace_index import TraceIndex
from rans_request_manifest import RequestManifest,sha,canonical


class BoundedGrid(PreparedRows):
    policy='CT16-BoundedBatchedRowsV1'
    maximum_plan_bytes=256*1024**2
    maximum_plan_rows=262144
    def __init__(self,initial,count):
        super().__init__();self.initial=initial;self.count=count;self.saved_values=[];self.plan_bytes=0;self.groups=[];self.value_calls=0
        o=initial['observation'];rows=initial['test_basis'];off=o['offsets']
        require(0<count<=len(rows)<=self.maximum_plan_rows and len(o['grid'])==len(rows)==len(o['grid_regions'])==len(o['grid_sides']),'invalid bounded grid inventory')
        self.digest=self.packet_digest()
        # One immutable velocity/k row plan serves three fields with their
        # separately supplied actual offsets. Pressure groups retain regions.
        for name,fields in [('velocity_k',(0,1,2)),('epsilon',(3,)),('pressure',(4,))]:
            grouped={}
            for i,q in enumerate(rows[:count]):
                row=q[name]
                if row is None:continue
                block=4+o['grid_regions'][i] if name=='pressure' else fields[0]
                grouped.setdefault((len(row['ids']),off[block]),[]).append((i,row))
            for (width,base),group in grouped.items():
                planned=len(group)*(8+width*8*5)
                require(width>0 and self.plan_bytes+planned<=self.maximum_plan_bytes,'bounded batched row allocation exceeded')
                indices=np.asarray([i for i,q in group],dtype=np.int64)
                ids=np.asarray([q['ids'] for i,q in group],dtype=np.int64)
                phi=np.asarray([q['phi'] for i,q in group],dtype=float)
                derivative=np.asarray([q['derivative'] for i,q in group],dtype=float)
                require(phi.shape==ids.shape and derivative.shape==(len(group),3,width) and np.isfinite(phi).all() and np.isfinite(derivative).all() and (ids>=0).all(),'malformed immutable native row')
                for a in (indices,ids,phi,derivative):a.flags.writeable=False
                self.plan_bytes+=sum(a.nbytes for a in (indices,ids,phi,derivative))
                require(self.plan_bytes<=self.maximum_plan_bytes,'bounded plan live bytes')
                self.groups.append((fields,base,indices,ids,phi,derivative))
    def packet_digest(self):
        o=self.initial['observation']
        return hashlib.sha256(canonical(dict(basis=self.initial['test_basis'][:self.count],grid=o['grid'][:self.count],regions=o['grid_regions'][:self.count],sides=o['grid_sides'][:self.count],offsets=o['offsets']))).hexdigest()
    def values(self,initial,state,count):
        if initial is not self.initial or count!=self.count:return original_values(initial,state,count)
        self.value_calls+=1;o=initial['observation'];s=np.asarray(state);off=o['offsets']
        require(s.ndim==1 and len(s)==len(o['state']),'batched actual state shape changed')
        v=np.zeros((count,5));g=np.zeros((count,5,3))
        for fields,base,indices,ids,phi,d in self.groups:
            for field in fields:
                block_offset=base if field==4 else off[field]
                current=s[ids+block_offset]
                v[indices,field]=np.einsum('ij,ij->i',phi,current,optimize=False)
                g[indices,field,:]=np.einsum('iaj,ij->ia',d,current,optimize=False)
        return dict(observation=dict(time=0.,grid=o['grid'][:count],grid_sides=o['grid_sides'][:count],grid_regions=o['grid_regions'][:count],grid_values=v.tolist(),grid_gradients=g.tolist()))
    def __enter__(self):
        super().__enter__()
        for module in list(sys.modules.values()):
            if module and getattr(module,'__name__','').startswith('ct_accepted') and getattr(module,'values',None) is original_values:
                self.saved_values.append(module);module.values=self.values
        return self
    def __exit__(self,*args):
        try:require(self.packet_digest()==self.digest,'immutable batched packet changed inside scope')
        finally:
            for module in self.saved_values:module.values=original_values
            super().__exit__(*args)


def same_decisions(a,b,key=''):
    if isinstance(a,dict):
        require(a.keys()==b.keys(),'batched observation fields changed')
        for k,v in a.items():same_decisions(v,b[k],k)
    elif isinstance(a,(list,tuple)):
        require(len(a)==len(b),'batched observation membership changed')
        if key in ('coordinate','geometries','raw_gradient','pressure_gradient'):return
        for x,y in zip(a,b,strict=True):same_decisions(x,y,key)
    elif isinstance(a,(bool,str)) or a is None or isinstance(a,int):require(a==b,'batched decision/identity changed:'+key)


def controls():
    row=dict(ids=[0],phi=[1.],derivative=[[0.],[0.],[0.]],geometry=7)
    packet=dict(observation=dict(grid=[[0.,0.,.1]],grid_sides=[0],grid_regions=[0],offsets=list(range(7)),state=[1.]*6),test_basis=[dict(velocity_k=row,epsilon=None,pressure=row)])
    plan=BoundedGrid(packet,1);state=np.asarray([1.,2.,3.,4.,5.,6.])
    with plan:
        require(plan.values(packet,state,1)==original_values(packet,state,1),'small exact batch fields')
        changed=state.copy();changed[0]+=1.;changed[2]+=2.
        require(plan.values(packet,changed,1)==original_values(packet,changed,1),'actual changed flow/k state reused')
        foreign=dict(packet)
        require(plan.values(foreign,state,1)==original_values(foreign,state,1),'foreign packet borrowed plan')
    class TooSmall(BoundedGrid):maximum_plan_bytes=1
    try:TooSmall(packet,1)
    except ValueError:pass
    else:raise ValueError('oversized immutable batch plan accepted')
    try:
        with BoundedGrid(packet,1):row['phi'][0]=2.
    except ValueError:pass
    else:raise ValueError('stale immutable batch metadata accepted')
    finally:row['phi'][0]=1.
    require(n.values is original_values,'batch hooks remained installed after failure')
    return dict(status='pass',changed_actual_flow_and_k_fresh=True,foreign_packet_fallback=True,byte_bound_rejected=True,stale_rows_rejected=True,hooks_restored_after_failure=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('ct14','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));n.reference=lru_cache(None)(n.reference)
    scripts=Path(__file__).parent;initial=a.ct14/'static-v1/target/frozen.json';paths=[a.ct14/'pair-v2'/dt/'case.evolution.json' for dt in ('0.00025','0.000125')]
    m=RequestManifest(a.output,dict(schema='CT16BoundedBatchedRowsQualificationV1',policy=BoundedGrid.policy,source_sha256=sha(Path(__file__)),initial_sha256=sha(initial),histories_sha256=[sha(f) for f in paths],bounds=dict(rows=BoundedGrid.maximum_plan_rows,bytes=BoundedGrid.maximum_plan_bytes,scalar_row_bytes=256*1024**2),changed_reduction_scaled_tolerance=2e-10,physical_advance=False))
    ini=n.read(initial);c=n.read(a.ct14/'static-v1/results.json')['case'];weights=n.inventory(c)[2];count=len(weights)
    states=[ini['observation']['state']]+[n.read(f)['state'] for f in paths]
    for q,region in zip(ini['observation']['grid'][:count],ini['observation']['grid_regions'][:count],strict=True):n.reference(n.col(c),q[2],bool(region))
    mapped={}
    for line in Path('/proc/self/maps').read_text().splitlines():
        fields=line.split(maxsplit=5)
        if len(fields)==6 and fields[5].startswith('/') and ('/lib' in fields[5] or '.so' in fields[5]):
            file=Path(fields[5]);require(file.is_file(),'unverifiable loaded Python library');mapped[str(file)]=sha(file)
    out=dict(schema='CT16BoundedBatchedRowsQualificationV1',status='incomplete',policy=BoundedGrid.policy,states=[],physical_advance=False,controls=controls(),
        Python_runtime=dict(executable=str(Path(sys.executable).resolve()),executable_sha256=sha(Path(sys.executable).resolve()),version=sys.version,numpy_version=np.__version__,actual_mapped_library_sha256=mapped))
    def evaluate(state):return bundle(ini,state,c,weights,count)
    for i,state in enumerate(states):
        cpu=time.process_time();wall=time.monotonic()
        with PreparedRows(),TraceIndex(ini,c['H']):
            before=evaluate(state)
            warm_cpu=time.process_time();warm_wall=time.monotonic();again=evaluate(state)
            require(before==again,'original repeated independent audit changed')
            old_warm=dict(CPU=time.process_time()-warm_cpu,wall=time.monotonic()-warm_wall)
        old=dict(CPU=time.process_time()-cpu-old_warm['CPU'],wall=time.monotonic()-wall-old_warm['wall'],warm_bundle=old_warm)
        cpu=time.process_time();wall=time.monotonic()
        with BoundedGrid(ini,count) as plan,TraceIndex(ini,c['H']):
            prepared=dict(CPU=time.process_time()-cpu,wall=time.monotonic()-wall)
            after=evaluate(state)
            warm_cpu=time.process_time();warm_wall=time.monotonic();again=evaluate(state)
            require(after==again,'batched repeated independent audit changed')
            new_warm=dict(CPU=time.process_time()-warm_cpu,wall=time.monotonic()-warm_wall)
        new=dict(CPU=time.process_time()-cpu-new_warm['CPU'],wall=time.monotonic()-wall-new_warm['wall'],immutable_preparation=prepared,warm_bundle=new_warm,plan_bytes=plan.plan_bytes,scalar_row_peak_bytes=plan.peak_bytes,value_calls=plan.value_calls,geometry_basis_sha256=plan.digest)
        same_decisions(before,after)
        for name in ('inventory','profiles','outlet'):require(before[name]==after[name],'unchanged scalar constitutive/traction inventory differs:'+name)
        for name in before['gradients']:
            if 'witness' in before['gradients'][name]:require(before['gradients'][name]['witness']==after['gradients'][name]['witness'],'batched gradient extremum witness changed:'+name)
        defect=maximum_difference(before,after)/min(c['ustar']**2,c['K0'],n.col(c).A/(c['H']-c['displacement']))
        require(defect<=2e-10 and after['metrics']['status']=='pass','batched complete physical observation equivalence/physical gate')
        write_json(a.output/('bundle-'+str(i)+'.json'),after)
        out['states'].append(dict(index=i,status='pass',all_decisions_and_scalar_witnesses_equal=True,maximum_scaled_difference=defect,original=old,batched=new,bundle_sha256=sha(a.output/('bundle-'+str(i)+'.json'))))
        write_json(a.output/'progress.json',out);print('CT16 bounded batch',i,old['CPU'],new['CPU'],defect,flush=True)
    out['status']='pass';write_json(a.output/'results.json',out);m.seal(a.output,['results.json'])
    return 0
if __name__=='__main__':raise SystemExit(main())
