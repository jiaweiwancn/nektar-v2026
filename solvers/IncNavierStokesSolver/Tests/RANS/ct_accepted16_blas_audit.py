#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Exact BLAS field rows, fresh per grid or original scalar physics point."""
import argparse,copy,ctypes,hashlib,json,os,resource,struct,subprocess,sys,time
from pathlib import Path
from functools import lru_cache
from rans_qualification import THREADS
os.environ.update(THREADS)
import numpy as np
import ct_accepted8_native as n
import ct_accepted3_metrics as metrics
from ct_accepted15_python import PreparedRows
from ct_accepted15_audits import bundle
from ct_accepted2_audit import values as original_values
from ct_accepted16_blas import BLAS,BLAS_SHA,load
from ct_accepted_open import require,run
from ct_accepted_feasibility import write_json
from rans_trace_index import TraceIndex
from rans_request_manifest import RequestManifest,sha,canonical,runtime_identity

ORIGINAL_POINT=metrics.point
owned_output=None


def exact_equal(a,b):
    if isinstance(a,np.ndarray):a=a.tolist()
    if isinstance(b,np.ndarray):b=b.tolist()
    if isinstance(a,dict):return isinstance(b,dict) and a.keys()==b.keys() and all(exact_equal(v,b[k]) for k,v in a.items())
    if isinstance(a,(list,tuple)):return isinstance(b,(list,tuple)) and len(a)==len(b) and all(exact_equal(x,y) for x,y in zip(a,b,strict=True))
    if isinstance(a,(bool,str)) or a is None:return a==b
    if isinstance(a,(int,np.integer)):return isinstance(b,(int,np.integer)) and a==b
    return struct.pack('!d',float(a))==struct.pack('!d',float(b))


class ExactAuditRows(PreparedRows):
    policy='CT16-ExactBLASAuditRowsV1'
    maximum_plan_bytes=256*1024**2
    maximum_rows=262144
    def __init__(self,initial,count,library):
        super().__init__();self.initial=initial;self.count=count;self.saved_values=[];self.saved_points=[];self.active=None;self.point_calls=0;self.grid_calls=0;self.plan_bytes=0;self.handle=None
        o=initial['observation'];self.offsets=tuple(o['offsets']);self.observation_rows=tuple(initial['test_basis'])+tuple(o['faces']);self.grid_size=len(initial['test_basis']);self.state_size=len(o['state'])
        regions=list(o['grid_regions'])+[q['region'] for q in o['faces']]
        require(0<count<=self.grid_size and len(self.observation_rows)==len(regions)<=self.maximum_rows and len(self.offsets)==7 and all(isinstance(v,(int,np.integer)) for v in self.offsets),'invalid owned BLAS row packet')
        self.row_indices={}
        for i,(row,region) in enumerate(zip(self.observation_rows,regions,strict=True)):
            require(region in (0,1) and (row.get('epsilon') is not None)==bool(region),'inconsistent physical row region/epsilon')
            prior=self.row_indices.setdefault(id(row),(row,i,region));require(prior[0] is row and prior[2]==region,'ambiguous native row owner')
        self.digest=self.packet_digest();self.arrays=[];self.widths=[]
        for name in ('velocity_k','epsilon','pressure'):
            widths={len(q[name]['ids']) for q in self.observation_rows if q.get(name) is not None}
            require(len(widths)==1,'mixed native row widths require the original route');w=widths.pop();self.widths.append(w)
            planned=len(self.observation_rows)*w*8*5
            require(0<w<=128 and self.plan_bytes+planned<=self.maximum_plan_bytes,'BLAS immutable row allocation bound')
            ids=np.zeros((len(self.observation_rows),w),dtype=np.int64);phi=np.zeros((len(self.observation_rows),w));d=np.zeros((len(self.observation_rows),3,w))
            for i,q in enumerate(self.observation_rows):
                if q.get(name) is None:continue
                require(all(isinstance(v,(int,np.integer)) for v in q[name]['ids']),'noninteger native field ID')
                ids[i]=q[name]['ids'];phi[i]=q[name]['phi'];d[i]=q[name]['derivative']
            for a in (ids,phi,d):a.flags.writeable=False
            self.arrays.extend((ids,phi,d));self.plan_bytes+=sum(a.nbytes for a in (ids,phi,d))
        self.present=np.asarray([q.get('epsilon') is not None for q in self.observation_rows],dtype=np.uint8)
        self.pressure_offsets=np.asarray([self.offsets[4+int(r)] for r in regions],dtype=np.int64)
        self.native_offsets=np.asarray(self.offsets,dtype=np.int64)
        for a in (self.present,self.pressure_offsets,self.native_offsets):a.flags.writeable=False;self.plan_bytes+=a.nbytes
        require(self.plan_bytes<=self.maximum_plan_bytes,'BLAS plan byte bound')
        self.lib,_=load(library);require(self.lib.ct16_plan_abi()==1,'foreign BLAS row-plan ABI')
        i64=np.ctypeslib.ndpointer(dtype=np.int64,flags='C_CONTIGUOUS');f64=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');u8=np.ctypeslib.ndpointer(dtype=np.uint8,flags='C_CONTIGUOUS')
        create=self.lib.ct16_plan_create;create.argtypes=[ctypes.c_int64,ctypes.c_int,ctypes.c_int,ctypes.c_int,i64,f64,f64,i64,f64,f64,u8,i64,f64,f64,i64,i64,ctypes.c_int64];create.restype=ctypes.c_void_p
        self.lib.ct16_plan_free.argtypes=[ctypes.c_void_p];self.lib.ct16_plan_free.restype=None
        self.lib.ct16_plan_point.argtypes=[ctypes.c_void_p,ctypes.c_int64,f64,ctypes.c_int64,f64];self.lib.ct16_plan_point.restype=ctypes.c_int
        self.lib.ct16_plan_grid.argtypes=[ctypes.c_void_p,ctypes.c_int64,f64,ctypes.c_int64,f64,f64];self.lib.ct16_plan_grid.restype=ctypes.c_int
        self.handle=create(len(self.observation_rows),*self.widths,*self.arrays[:6],self.present,*self.arrays[6:],self.pressure_offsets,self.native_offsets,self.state_size)
        require(bool(self.handle),'malformed or foreign field row plan rejected before use')
    def packet_digest(self):
        o=self.initial['observation']
        return hashlib.sha256(canonical(dict(rows=self.initial['test_basis'],faces=o['faces'],grid=o['grid'],sides=o['grid_sides'],regions=o['grid_regions'],offsets=o['offsets'],actual_initial_state=o['state']))).hexdigest()
    def close(self):
        if self.handle:self.lib.ct16_plan_free(self.handle);self.handle=None
    def values(self,initial,state,count):
        if initial is not self.initial or count!=self.count:return original_values(initial,state,count)
        require(bool(self.handle),'BLAS row scope already closed');s=np.ascontiguousarray(state,dtype=np.float64)
        require(s.shape==(self.state_size,),'foreign actual field state shape')
        v=np.empty((count,5));g=np.empty((count,5,3))
        code=self.lib.ct16_plan_grid(self.handle,count,s,len(s),v,g);require(code==0,'BLAS grid input rejected:'+str(code));self.grid_calls+=1
        o=initial['observation']
        # Consumers already call asarray; retaining arrays avoids a list round
        # trip. The exported bundle still contains the original scalar/list data.
        return dict(observation=dict(time=0.,grid=o['grid'][:count],grid_sides=o['grid_sides'][:count],grid_regions=o['grid_regions'][:count],grid_values=v,grid_gradients=g))
    def point(self,row,state,offsets,model):
        entry=self.row_indices.get(id(row))
        if not (self.handle and entry and entry[0] is row and entry[2]==1 and model=='A95' and tuple(offsets)==self.offsets and all(isinstance(v,(int,np.integer)) for v in offsets) and isinstance(state,np.ndarray) and state.dtype==np.float64 and state.flags.c_contiguous and state.shape==(self.state_size,)):
            return ORIGINAL_POINT(row,state,offsets,model)
        output=np.empty((5,4));code=self.lib.ct16_plan_point(self.handle,entry[1],state,len(state),output)
        require(code==0,'BLAS point input rejected:'+str(code));self.point_calls+=1
        # Values are freshly computed for THIS original point call. The original
        # scalar domain, invariants, full signed tensors and witnesses below run
        # unchanged. There is no value/state cache across points or grids.
        fields={}
        for b in range(3):fields[(id(row['velocity_k']),self.offsets[b])]=(row['velocity_k'],b)
        fields[(id(row['epsilon']),self.offsets[3])]=(row['epsilon'],3)
        fields[(id(row['pressure']),self.offsets[5])]=(row['pressure'],4)
        previous=self.active;self.active=(state,fields,output)
        try:return ORIGINAL_POINT(row,state,offsets,model)
        finally:self.active=previous
    def __call__(self,row,state,offset):
        active=self.active
        if active is not None and state is active[0]:
            entry=active[1].get((id(row),int(offset)))
            if entry and entry[0] is row:
                b=entry[1];return active[2][b,0],active[2][b,1:]
        return super().__call__(row,state,offset)
    def __enter__(self):
        require(bool(self.handle) and self.active is None,'invalid BLAS audit scope')
        super().__enter__()
        for module in list(sys.modules.values()):
            if not module or not getattr(module,'__name__','').startswith('ct_accepted'):continue
            if getattr(module,'values',None) is original_values:self.saved_values.append(module);module.values=self.values
            if getattr(module,'point',None) is ORIGINAL_POINT:self.saved_points.append(module);module.point=self.point
        return self
    def __exit__(self,*args):
        try:require(self.active is None and self.packet_digest()==self.digest,'immutable BLAS packet changed inside audit scope')
        finally:
            for module in self.saved_points:module.point=ORIGINAL_POINT
            for module in self.saved_values:module.values=original_values
            super().__exit__(*args);self.close()


def controls(library):
    vk=dict(ids=[0],phi=[1.],derivative=[[0.],[0.],[1.]],geometry=7)
    row=dict(velocity_k=vk,epsilon=copy.deepcopy(vk),pressure=copy.deepcopy(vk))
    packet=dict(observation=dict(grid=[[0.,0.,1.]],grid_sides=[0],grid_regions=[1],offsets=list(range(7)),state=[1.,0.,1.,1.,0.,0.],faces=[]),test_basis=[row])
    state=np.asarray(packet['observation']['state']);off=packet['observation']['offsets']
    before=ORIGINAL_POINT(row,state,off,'A95')
    with ExactAuditRows(packet,1,library) as plan:
        require(exact_equal(plan.point(row,state,off,'A95'),before),'small full signed point differs')
        state[0]=2.;state[1]=.05;state[2]=1.2;state[3]=.8
        require(exact_equal(plan.point(row,state,off,'A95'),ORIGINAL_POINT(row,state,off,'A95')),'changed flow/k/epsilon reused an old point value')
        foreign=dict(row)
        require(exact_equal(plan.point(foreign,state,off,'A95'),ORIGINAL_POINT(foreign,state,off,'A95')),'foreign row borrowed a point')
        require(exact_equal(plan.values(dict(packet),state,1),original_values(dict(packet),state,1)),'foreign packet borrowed a plan')
        require(exact_equal(plan.values(packet,state,1),original_values(packet,state,1)),'owned value/gradient grids differ')
        for index,count in [(-1,6),(0,5)]:
            target=np.full((5,4),123.);require(plan.lib.ct16_plan_point(plan.handle,index,state,count,target)!=0 and (target==123.).all(),'invalid point index/state shape wrote output')
    bad=copy.deepcopy(packet);bad['test_basis'][0]['velocity_k']['ids'][0]=99
    try:ExactAuditRows(bad,1,library)
    except ValueError:pass
    else:raise ValueError('out-of-field native IDs accepted')
    try:
        with ExactAuditRows(packet,1,library):vk['phi'][0]=2.
    except ValueError:pass
    else:raise ValueError('mutated immutable row accepted')
    finally:vk['phi'][0]=1.
    require(metrics.point is ORIGINAL_POINT and n.values is original_values,'audit hooks remained after failure')
    return dict(status='pass',actual_flow_k_epsilon_fresh=True,foreign_row_and_packet_fallback=True,all_IDs_and_shapes_checked=True,immutable_mutation_rejected=True,scope_cleanup_after_failure=True)


def main():
    global owned_output
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('ct14','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));n.reference=lru_cache(None)(n.reference)
    source=Path(__file__).with_name('ct16_independent_rows.cpp');initial=a.ct14/'static-v1/target/frozen.json';paths=[a.ct14/'pair-v2'/dt/'case.evolution.json' for dt in ('0.00025','0.000125')]
    require(sha(BLAS)==BLAS_SHA,'actual NumPy BLAS bytes changed')
    names=['ct_accepted16_blas_audit.py','ct16_independent_rows.cpp','ct_accepted16_blas.py','ct_accepted15_python.py','ct_accepted15_audits.py','ct_accepted3_metrics.py','ct_accepted2_diagnosis.py','ct_accepted2_audit.py','ct_accepted8_native.py','ct_accepted8_reference.py','ct_accepted10_compare.py','ct_accepted5_traction.py','ct_accepted4_run.py','ct_accepted4_reference.py','ct_accepted_open.py','validate_a93_a95_e94e95.py','rans_trace_index.py']
    m=RequestManifest(a.output,dict(schema='CT16ExactBLASAuditQualificationV1',policy=ExactAuditRows.policy,source_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=Path(__file__).resolve().parents[4],text=True).strip(),source_sha256={f:sha(Path(__file__).with_name(f)) for f in names},initial_sha256=sha(initial),case_sha256=sha(a.ct14/'static-v1/results.json'),actual_histories_sha256=[sha(f) for f in paths],full_bundle_binary64_required=True,no_native_advance=True));owned_output=a.output
    library=a.output/'CT16IndependentRows.so';execution=run(['c++','-std=c++17','-O3','-fPIC','-shared','-ffp-contract=off','-fno-fast-math',str(source),str(BLAS),'-o',str(library)],a.output,'build.log');require(execution['exit_code']==0,'private exact audit helper build failed')
    ini=n.read(initial);c=n.read(a.ct14/'static-v1/results.json')['case'];weights=n.inventory(c)[2];count=len(weights);states=[ini['observation']['state']]+[n.read(f)['state'] for f in paths]
    for q,region in zip(ini['observation']['grid'][:count],ini['observation']['grid_regions'][:count],strict=True):n.reference(n.col(c),q[2],bool(region))
    out=dict(schema='CT16ExactBLASAuditQualificationV1',status='incomplete',policy=ExactAuditRows.policy,build_execution=execution,library=runtime_identity(library),controls=controls(library),states=[],physical_advance=False)
    for i,state in enumerate(states):
        cpu=time.process_time();wall=time.monotonic()
        with PreparedRows(),TraceIndex(ini,c['H']):
            before=bundle(ini,state,c,weights,count);wc=time.process_time();ww=time.monotonic();again=bundle(ini,state,c,weights,count)
            require(exact_equal(before,again),'original repeated audit changed');old_warm=dict(CPU=time.process_time()-wc,wall=time.monotonic()-ww)
        old=dict(CPU=time.process_time()-cpu-old_warm['CPU'],wall=time.monotonic()-wall-old_warm['wall'],warm_bundle=old_warm)
        cpu=time.process_time();wall=time.monotonic()
        with ExactAuditRows(ini,count,library) as rows,TraceIndex(ini,c['H']):
            prepared=dict(CPU=time.process_time()-cpu,wall=time.monotonic()-wall);after=bundle(ini,state,c,weights,count);wc=time.process_time();ww=time.monotonic();again=bundle(ini,state,c,weights,count)
            require(exact_equal(after,again),'exact BLAS repeated audit changed');new_warm=dict(CPU=time.process_time()-wc,wall=time.monotonic()-ww)
        new=dict(CPU=time.process_time()-cpu-new_warm['CPU'],wall=time.monotonic()-wall-new_warm['wall'],warm_bundle=new_warm,immutable_preparation=prepared,plan_bytes=rows.plan_bytes,point_calls=rows.point_calls,grid_calls=rows.grid_calls,metadata_sha256=rows.digest)
        write_json(a.output/('original-bundle-'+str(i)+'.json'),before);write_json(a.output/('candidate-bundle-'+str(i)+'.json'),after)
        require(exact_equal(before,after) and after['metrics']['status']=='pass' and rows.point_calls>0 and rows.grid_calls>0,'exact BLAS full physical observations/decisions/witnesses differ')
        out['states'].append(dict(index=i,status='pass',all_observations_decisions_witnesses_exact_equal=True,original=old,exact_BLAS=new,bundle_sha256=sha(a.output/('candidate-bundle-'+str(i)+'.json'))));write_json(a.output/'progress.json',out);print('CT16 full exact BLAS audit',i,old['CPU'],new['CPU'],old_warm['CPU'],new_warm['CPU'],flush=True)
    mapped={}
    for line in Path('/proc/self/maps').read_text().splitlines():
        f=line.split(maxsplit=5)
        if len(f)==6 and f[5].startswith('/') and ('/lib' in f[5] or '.so' in f[5]):file=Path(f[5]);require(file.is_file(),'unverifiable loaded audit library');mapped[str(file)]=sha(file)
    out['Python_runtime']=dict(executable_sha256=sha(Path(sys.executable).resolve()),version=sys.version,numpy_version=np.__version__,actual_mapped_library_sha256=mapped)
    out['status']='pass';write_json(a.output/'results.json',out);m.seal(a.output,['results.json']);return 0
if __name__=='__main__':
    try:raise SystemExit(main())
    except Exception as error:
        if owned_output is not None:
            with (owned_output/'failure.json').open('x') as f:json.dump(dict(status='fail',type=type(error).__name__,reason=str(error),physical_advance=False),f,indent=2)
        raise
