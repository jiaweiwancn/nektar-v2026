#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Measure a genuine selected fresh start, including spatial/history construction."""
import argparse,os,subprocess,time
from pathlib import Path
import xml.etree.ElementTree as ET
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,verify,atomic,assess,vector_hash
from wj_h16_observer import compatibility
from wj_h17_observer import native_load,equal,load_compatible
from wj_h20_observer import common_points
from wj3_measure import checkpoint_payload
from cls_epawf_m3_cases import set_info,set_parameter,write


def run(root):
    base=root/'resumption-r1';d=base/'science/fresh-C-v1';d.mkdir()
    selected=load(Path('/home/jay/build/wj-h24-flat-continuation-20261115/inputs/selected.json'))
    source=selected['cp0']['C'];x=ET.parse(verify(source['xml'])).getroot()
    for k,v in dict(NumSteps=0,Time=0,TimeStep=.005,IO_CheckSteps=1).items():set_parameter(x,k,v)
    for k,v in dict(RANSWallRestart='None',RANSWallCheckpoint=str(d/'CP'),RANSWallRetention='Off',
        RANSWallMonitoringPolicy='Off',RANSWallReport=str(d/'native.json'),RANSWallNativeLookup='TensorExactV1',
        RANSWallScalarEvaluation='Off',RANSWallAssemblyScatter='ElementExactV1',
        RANSWallPimpleSnapshots='Off',RANSWallPimpleFlowVectors='Off',RANSWallPimpleArithmetic='Off',RANSWallPimpleVerify='Off').items():set_info(x,k,v)
    write(x,d/'case.xml')
    runtime_path=base/'science/runtime-build-scatter-v1.json';runtime=load(runtime_path)
    native=verify(next(e['wrapper']for e in runtime['executables']if Path(e['wrapper']['path']).name=='IncNavierStokesSolver'))
    began=time.monotonic()
    with(d/'native.log').open('x')as f:code=subprocess.run([str(native),str(d/'case.xml')],stdout=f,stderr=subprocess.STDOUT).returncode
    native_seconds=time.monotonic()-began
    atomic(d/'launch-result.json',dict(exit=code,native_complete_seconds=native_seconds))
    assert code==0
    r=native_load(d/'case.evolution.json')
    assert r['steps']==0 and r['time']==0 and r['failure']==''
    cp=checkpoint_payload(d/'CP',include_reference=True)
    old=checkpoint_payload(verify(source['checkpoint']).parent,include_reference=True)
    for k in ('state','fixed_reference','fixed_defect','eta','native_history','native_history_times','dt','signature'):
        equal(cp[k],old[k],'genuine fresh C initial '+k)
    c=compatibility(root,verify(load(root/'science/cache-compatible-v2.json')['historical_manifest']),d/'cache-base.json')
    c.update(new_native_identity=identity(runtime_path),permitted_difference='fresh native same physical reference/space; exact initial state/history parity with genuine original C CP0')
    atomic(d/'cache-compatible.json',c);cache,_=load_compatible(d/'cache-compatible.json')
    fixed=load(verify(selected['common_support']))['points'];points=common_points(cache,fixed,{'fresh':r},{'fresh':r})
    rows,initial,stats=cache.view(points);began=time.monotonic()
    physical,data=assess(rows,r['state'],initial,0);assessment_seconds=time.monotonic()-began
    atomic(d/'independent.json',dict(status=physical['summary']['status'],physical=physical,cache=stats,
        assessment_seconds=assessment_seconds,source=identity(d/'case.evolution.json')))
    assert physical['summary']['status']=='PASS'
    atomic(d/'qualification.json',dict(status='PASS_GENUINE_FRESH_INITIALIZATION',checkpoint=identity(d/'CP/State.bin'),
        producer=identity(d/'case.evolution.json'),original=source['checkpoint'],state=vector_hash(cp['state']),
        native_complete_seconds=native_seconds,independent_assessment_seconds=assessment_seconds,
        work=r.get('initialization_work',r.get('initialization_lift_calls_iterations')),
        initializer_calls=r.get('initializer_calls',r.get('initial_velocity_calls')),native_history_initializations=1,
        governing_advances=0,scope='fresh spatial projection/lift/reference/defect/FlatInletTrace and native own-dt history; no restart or seed-import shortcut'))
    print('genuine fresh C PASS',native_seconds,assessment_seconds,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);run(p.parse_args().root)
