#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only actual-coordinate geometry microbenchmark; no native advance."""
import argparse,ctypes,json,os,resource,subprocess,time
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
import numpy as np
from ct_accepted_open import require
from ct_accepted_feasibility import write_json
from rans_request_manifest import RequestManifest,sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('ct14','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));resource.setrlimit(resource.RLIMIT_CPU,(100,105))
    source=Path(__file__).with_name('CT16FlatMetricBench.cpp');root=Path(__file__).parents[4];native=root/'library/SolverUtils/Turbulence/RANSWallNormalMap.cpp';header=native.with_suffix('.h');initial=a.ct14/'static-v1/target/frozen.json'
    m=RequestManifest(a.output,dict(schema='CT16ExactFlatMetricPrototypeV1',source_sha256=sha(source),native_geometry_sha256=sha(native),header_sha256=sha(header),controller_sha256=sha(Path(__file__)),actual_packet_sha256=sha(initial),maximum_coordinates=65536,repetitions=32,timed_pairs=3,physical_advance=False,scope='read-only sampled immutable geometry primitives; no governing speedup/admission claim'))
    libfile=a.output/'flat-metric.so';command=['c++','-std=c++17','-O3','-DNDEBUG','-shared','-fPIC','-I'+str(root/'library'),str(source),str(native),'-o',str(libfile)]
    before=resource.getrusage(resource.RUSAGE_CHILDREN);start=time.monotonic();e=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True);(a.output/'compiler.log').write_text(e.stdout);require(e.returncode==0,'flat metric prototype compile')
    after=resource.getrusage(resource.RUSAGE_CHILDREN);build=dict(command=command,CPU=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,wall=time.monotonic()-start,binary_sha256=sha(libfile))
    packet=json.loads(initial.read_text());o=packet['observation'];groups=[o['grid'],[q['coordinate'] for q in o['volume']],[q['coordinate'] for q in o['faces']]];coordinates=[]
    for group in groups:
        coordinates += [group[i] for i in np.linspace(0,len(group)-1,min(16384,len(group)),dtype=int)]
    coordinates=np.asarray(coordinates,float);require(0<len(coordinates)<=65536,'coordinate bound')
    c=json.loads((a.ct14/'static-v1/results.json').read_text())['case'];bounds=np.asarray([min(c['x']),0.,max(c['x']),c['Z']],float)
    lib=ctypes.CDLL(str(libfile));lib.ct16_flat_metric_error.restype=ctypes.c_char_p;lib.ct16_flat_metric_run.argtypes=[ctypes.c_size_t,ctypes.c_void_p,ctypes.c_double,ctypes.c_void_p,ctypes.c_int,ctypes.c_int,ctypes.c_void_p]
    output=[np.zeros((len(coordinates),3)) for _ in range(2)]
    def run(policy,repeats):
        code=lib.ct16_flat_metric_run(len(coordinates),coordinates.ctypes.data,c['H'],bounds.ctypes.data,repeats,policy,output[policy].ctypes.data);require(code==0,'flat metric primitive:'+str(lib.ct16_flat_metric_error()))
    run(0,1);run(1,1);require(np.array_equal(output[0].view(np.uint64),output[1].view(np.uint64)),'flat metric values changed bits')
    timings={'0':[],'1':[]}
    for repeat in range(3):
        for policy in ((0,1) if repeat%2==0 else (1,0)):
            cpu=time.process_time();wall=time.monotonic();run(policy,32);timings[str(policy)].append(dict(CPU=time.process_time()-cpu,wall=time.monotonic()-wall))
    require(np.array_equal(output[0].view(np.uint64),output[1].view(np.uint64)),'timed metric values changed')
    out=dict(schema='CT16ExactFlatMetricPrototypeV1',status='pass',coordinates=len(coordinates),repetitions=32,timed_pairs=3,build=build,binary64_equal=True,maximum_difference=0.,timing=timings,physical_advance=False,interpretation='Exact sampled immutable geometry values and primitive timing only. Governing opt-in implementation,invalid-input/current-point controls and matched native timing are still required. Actual coordinates are never rounded or substituted.')
    write_json(a.output/'results.json',out);m.seal(a.output,['results.json']);print('CT16 flat metric prototype pass',len(coordinates),flush=True)


if __name__=='__main__':main()
