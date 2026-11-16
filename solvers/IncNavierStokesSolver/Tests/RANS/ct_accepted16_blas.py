#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Qualify exact C gather/BLAS evaluation against actual NumPy row arithmetic."""
import argparse,ctypes,json,os,resource,sys,time
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
import numpy as np
import ct_accepted8_native as n
from ct_accepted16_bulk import BoundedGrid
from ct_accepted_feasibility import write_json
from ct_accepted_open import require,run
from rans_request_manifest import RequestManifest,sha,runtime_identity

POLICY='CT16-IndependentBLASRowsV1'
BLAS=Path('/usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3')
BLAS_SHA='babdd0d2612d7bd2f0099ac8c82714702fb5bb47b1e2e6fdba819c97c79fc0fd'


def load(path):
    lib=ctypes.CDLL(str(path));require(lib.ct16_row_abi()==1,'foreign independent row ABI')
    i64=np.ctypeslib.ndpointer(dtype=np.int64,flags='C_CONTIGUOUS')
    f64=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
    fn=lib.ct16_gather_blas_rows
    fn.argtypes=[ctypes.c_int64,ctypes.c_int,i64,f64,f64,f64,ctypes.c_int64,ctypes.c_int64,f64,f64]
    fn.restype=ctypes.c_int
    return lib,fn


def execute(fn,ids,phi,d,state,offset):
    require(ids.ndim==2 and phi.shape==ids.shape and d.shape==(len(ids),3,ids.shape[1]) and state.ndim==1,'independent row shapes differ')
    v=np.empty(len(ids));g=np.empty((len(ids),3))
    code=fn(len(ids),ids.shape[1],ids,phi,d,state,len(state),offset,v,g)
    require(code==0,'independent row input rejected:'+str(code));return v,g


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('ct14','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3))
    source=Path(__file__).with_name('ct16_independent_rows.cpp');initial=a.ct14/'static-v1/target/frozen.json';paths=[a.ct14/'pair-v2'/dt/'case.evolution.json' for dt in ('0.00025','0.000125')]
    require(sha(BLAS)==BLAS_SHA,'actual qualified NumPy BLAS library changed')
    m=RequestManifest(a.output,dict(schema='CT16IndependentBLASRowsQualificationV1',policy=POLICY,source_sha256={source.name:sha(source),Path(__file__).name:sha(Path(__file__)), 'ct_accepted16_bulk.py':sha(Path(__file__).with_name('ct_accepted16_bulk.py'))},actual_initial_sha256=sha(initial),actual_histories_sha256=[sha(f) for f in paths],BLAS_path=str(BLAS),BLAS_sha256=BLAS_SHA,bitwise_required=True,physical_advance=False))
    library=a.output/'CT16IndependentRows.so'
    command=['c++','-std=c++17','-O3','-fPIC','-shared','-ffp-contract=off','-fno-fast-math',str(source),str(BLAS),'-o',str(library)]
    execution=run(command,a.output,'build.log');require(execution['exit_code']==0,'independent row helper build failed')
    lib,fn=load(library);ini=n.read(initial);c=n.read(a.ct14/'static-v1/results.json')['case'];plan=BoundedGrid(ini,len(n.inventory(c)[2]));states=[ini['observation']['state']]+[n.read(f)['state'] for f in paths]
    out=dict(schema='CT16IndependentBLASRowsQualificationV1',status='incomplete',policy=POLICY,build_execution=execution,library=runtime_identity(library),plan_bytes=plan.plan_bytes,states=[],physical_advance=False)
    for i,state in enumerate(states):
        state=np.asarray(state);cpu=time.process_time();wall=time.monotonic();truth=[]
        for fields,base,index,ids,phi,d in plan.groups:
            for field in fields:
                offset=base if field==4 else ini['observation']['offsets'][field]
                v=np.empty(len(ids));g=np.empty((len(ids),3))
                for j in range(len(ids)):
                    actual=state[ids[j]+offset];v[j]=phi[j]@actual;g[j]=d[j]@actual
                truth.append((ids,phi,d,offset,v,g))
        old=dict(CPU=time.process_time()-cpu,wall=time.monotonic()-wall)
        cpu=time.process_time();wall=time.monotonic();comparisons=[]
        for ids,phi,d,offset,v,g in truth:
            nv,ng=execute(fn,ids,phi,d,state,offset)
            same_v=bool(np.array_equal(v.view(np.uint64),nv.view(np.uint64)));same_g=bool(np.array_equal(g.view(np.uint64),ng.view(np.uint64)))
            comparisons.append(dict(rows=len(ids),width=ids.shape[1],offset=offset,values_binary64_equal=same_v,gradients_binary64_equal=same_g,maximum_value_difference=float(np.max(abs(v-nv),initial=0)),maximum_gradient_difference=float(np.max(abs(g-ng),initial=0))))
        new=dict(CPU=time.process_time()-cpu,wall=time.monotonic()-wall)
        result=dict(index=i,status='pass' if all(r['values_binary64_equal'] and r['gradients_binary64_equal'] for r in comparisons) else 'fail',original=old,native_gather_BLAS=new,comparisons=comparisons)
        out['states'].append(result);write_json(a.output/'progress.json',out);print('CT16 exact BLAS rows',i,result['status'],old['CPU'],new['CPU'],flush=True)
        require(result['status']=='pass','C BLAS arithmetic differs from original NumPy rows; do not adopt')
    ids=np.asarray([[0]],dtype=np.int64);phi=np.asarray([[1.]]);d=np.zeros((1,3,1));s=np.asarray([1.]);v=np.asarray([123.]);g=np.asarray([[123.]*3])
    for bad_ids,offset in [(np.asarray([[-1]],dtype=np.int64),0),(np.asarray([[1]],dtype=np.int64),0),(ids,-1),(ids,2)]:
        require(fn(1,1,bad_ids,phi,d,s,1,offset,v,g)!=0 and v[0]==123. and (g==123.).all(),'invalid BLAS IDs/offset wrote output')
    out['controls']=dict(status='pass',bounds_and_all_IDs_checked_before_writing=True)
    mapped={}
    for line in Path('/proc/self/maps').read_text().splitlines():
        fields=line.split(maxsplit=5)
        if len(fields)==6 and fields[5].startswith('/') and ('/lib' in fields[5] or '.so' in fields[5]):
            file=Path(fields[5]);require(file.is_file(),'unverifiable mapped BLAS audit library');mapped[str(file)]=sha(file)
    out['Python_runtime']=dict(executable_sha256=sha(Path(sys.executable).resolve()),version=sys.version,numpy_version=np.__version__,mapped_library_sha256=mapped)
    out['status']='pass';write_json(a.output/'results.json',out);m.seal(a.output,['results.json']);return 0
if __name__=='__main__':raise SystemExit(main())
