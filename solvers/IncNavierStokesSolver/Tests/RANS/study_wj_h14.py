#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Narrow fixed-six-millimetre H14 qualification and immutable input adapters.

Historical states below are operator fixtures, never target restart states.
"""
import argparse
import json
from pathlib import Path
import shlex
import shutil
import xml.etree.ElementTree as ET

from run_wj_h4 import load,save,sha
from wj_h14_execution import clone,file_id,validate_pair
from wj_h7_analysis import context
from wj_h4_points import vector_hash
from cls_epawf_m3_cases import set_info,set_parameter,write

H10=Path('/home/jay/build/wj-h10-flat-inlet-20261113')
H13=Path('/home/jay/build/wj-h13-fullheight-paired-domain-20261113')

def verify(v):
    p=Path(v['path'])
    assert p.stat().st_size==v['bytes'] and sha(p)==v['sha256'],str(p)
    return p

def runtime(root):
    rt=root/'runtime';rt.mkdir();(rt/'glibc').mkdir();(rt/'geometry').mkdir()
    old=load(H10/'setup.json');aliases=[]
    for rec in old['private_libraries']+[old['geometry']]:
        src=verify(rec);dest=rt/('geometry'if rec==old['geometry']else'glibc')/src.name
        shutil.copyfile(src,dest);dest.chmod(0o755)
        assert sha(dest)==rec['sha256'];aliases.append(dict(original=rec,copy=file_id(dest)))
    executables=[]
    for name in ('IncNavierStokesSolver','RANSWallResolutionProbe'):
        source=root/'build'/name;dest=rt/(name+'.elf');shutil.copyfile(source,dest);dest.chmod(0o755)
        wrapper=rt/name
        wrapper.write_text('#!/bin/sh\nexec '+shlex.quote(str(rt/'glibc/ld-linux-x86-64.so.2'))+' --library-path '+
            shlex.quote(str(rt/'glibc')+':'+str(rt/'geometry')+':/home/jay/.local/opt/nektar-v2026-mpi/lib64')+' '+shlex.quote(str(dest))+' "$@"\n')
        wrapper.chmod(0o755);executables.append(dict(build=file_id(source),ELF=file_id(dest),wrapper=file_id(wrapper)))
    commands=load(root/'build/compile_commands.json')
    assert all('-ffp-contract=off'in x['command']and'-fno-fast-math'in x['command']for x in commands)
    save(root/'science/runtime-build.json',dict(executables=executables,private_libraries=aliases,
        compile_commands=file_id(root/'build/compile_commands.json'),compiled_inputs=[dict(path=x['file'],sha256=sha(x['file']),command=x['command'])for x in commands],
        source_capsules=[file_id(p)for p in sorted(root.glob('*source*.json'))],
        compiler='/usr/bin/c++ GNU11.4.0; Release, serial native rank, compile parallel8',
        scope='new H14 governor/reader, exact H10 original glibc/geometry bytes; actual native-exec maps required separately'))

def inputs(root):
    out=root/'science/small';out.mkdir();bindings=load(H10/'publication-v1/state-bindings.json')
    cp0rec=bindings['coarse-cp0']['producer'];runrec=bindings['coarse-run']['producer']
    cp0=load(verify(cp0rec));r=load(verify(runrec));ctx=context(r,1)
    assert ctx['step']==1 and ctx['alpha']==100 and ctx['time']==.01
    assert vector_hash(cp0['state'])==vector_hash(r['accepted_states'][0])
    xml=H10/'science/cases/coarse/cp0/case.xml'
    assert sha(xml)=='f33d721a2f8bd3ceb6cba40f769b15d10bf6a648a1f9e922e30a02c7dadc3eb3'
    case=clone(xml,out/'operator',dt=.01)
    packet=dict(operation='mapped_backend_checks',inner_stress_model='Linear',state=cp0['state'],history=ctx['history'],
        fixed_reference=cp0['fixed_reference'],fixed_defect=cp0['fixed_defect'],alpha=100.,time=.01,
        state_sha256=vector_hash(cp0['state']),history_sha256=vector_hash(ctx['history']))
    save(case.parent/'request.json',packet)
    original=H10/'science/fixture-v2/case.xml'
    assert sha(original)=='08abd8398c2e682e783a1c434c1d96050b0693fb4997981f38f5d418776bf84c'
    fixture=clone(original,out/'initializer',dt=.01)
    save(fixture.parent/'request.json',dict(operation='initializer_checks',inner_stress_model='Linear'))
    cp=clone(xml,out/'cp0',dt=.01)
    runs={}
    for name,steps,restart in (('two',2,out/'cp0/CP'),('one',1,out/'cp0/CP'),('restart',1,out/'one/CP.1'),('reject',1,out/'one/CP.1')):
        path=clone(xml,out/name,dt=.01,steps=steps,restart=restart)
        if name=='reject':
            x=ET.parse(path).getroot();set_parameter(x,'RANSWallPimpleHardAuditRejectStep',2);write(x,path)
        runs[name]=file_id(path)
    save(out/'intent.json',dict(scope='seven maximum small verification attempts: two parity sweeps, four integration, one deliberate rejection; no target',
        original_cp0_report=cp0rec,original_step1_report=runrec,original_xml=file_id(xml),operator_xml=file_id(case),
        operator_packet=file_id(case.parent/'request.json'),initializer_original=file_id(original),initializer_xml=file_id(fixture),cp0=file_id(cp),runs=runs,
        height_mm=6,steps=7,state=2328,operator_state_role='original CP0 as test vector; original actual step1 callback history; never a checkpoint relabeling'))

def small_checks(root):
    """Compare real new histories; the rejected attempt stays rejected."""
    import numpy as np
    from wj_h9_local import same,phase_bindings
    from wj3_measure import checkpoint_payload
    out=root/'science/small';reports={k:load(out/k/'case.evolution.json')for k in ('cp0','two','one','restart','reject')}
    cp0=reports['cp0'];two=reports['two'];one=reports['one'];restart=reports['restart'];reject=reports['reject']
    assert cp0['initial_velocity_calls']==1 and cp0['initial_velocity_mode']=='FlatInletTraceV1'
    assert cp0['steps']==0 and two['steps']==restart['steps']==2 and one['steps']==1
    err=float(np.max(abs(np.array(two['state'])-restart['state'])))
    assert err<=1e-9*max(1.,max(abs(x)for x in two['state']))
    counts={};phase=[]
    for name,r in reports.items():
        assert r['initial_velocity_mode']=='FlatInletTraceV1'
        for key in ('fixed_reference','fixed_defect'):assert same(r[key],cp0[key])
        if name!='cp0':
            assert r['initial_velocity_calls']==0 and r['initialization_lift_calls_iterations']==[0,0]
        for attempt in r['pimple_attempts']:
            assert attempt['scheduled_complete']
            p=load(attempt['k_snapshot']);assert sha(attempt['k_snapshot'])==attempt['k_snapshot_sha256']
            assert same(p['raw'],p['corrected'])and p['k_bound']['activations']==0 and p['k_bound']['mass_norm_correction']==0
            assert p['matrix_encoding'].startswith('COO native-endian int32 row,int32 column') and p['schema']==2
        counts[name]=dict(steps=r['steps'],time=r['time'],failure=r['failure'],attempts=len(r['pimple_attempts']),
            initialization=r['initialization_lift_calls_iterations'],initializer_calls=r['initial_velocity_calls'],
            producer=file_id(out/name/'case.evolution.json'),checkpoint=r['checkpoint'])
    for key in ('state','native_history_times'):assert same(two[key],restart[key])
    assert len(two['native_history'])==len(restart['native_history'])
    assert all(same(x,y)for x,y in zip(two['native_history'],restart['native_history']))
    assert reject['failure']and reject['steps']==1 and reject['time']==.01 and same(reject['state'],one['state'])
    assert len(reject['physical_range_events'])>0
    for p in(out/'reject').glob('CP*/State.bin'):
        c=checkpoint_payload(p.parent);assert c['steps']<=1
    for name in ('two','one','restart'):
        r=reports[name]
        first=r['steps']+1-len(r['accepted_states'])
        for step in range(first+1,r['steps']+1):
            a=r['pimple_attempts'][step-first-1];p=load(a['k_snapshot']);_,v=phase_bindings(r,p,step);assert v['P_C_bitwise'];phase.append(dict(case=name,**v))
    save(out/'restart-rollback-checks.json',dict(status='PASS',state_max_error=err,state_bitwise=same(two['state'],restart['state']),
        native_history_bitwise=True,reference_defect_bitwise=True,rejected_attempt_rollback=True,counts=counts,phase_bindings=phase))

def view(root,case,reuse=False):
    folder=root/'science/small'/case;rp=folder/'case.evolution.json';r=load(rp);out=folder/'view';out.mkdir()
    save(out/'request.json',dict(operation='open_checks',inner_stress_model='Linear',native_report=str(rp),native_report_sha256=sha(rp),
        native_checkpoint=r['checkpoint'],checkpoint_sha256=sha(Path(r['checkpoint'])/'State.bin'),
        physical_observations=[[0.,0.,.2375],[4.,0.,.2375],[8.,0.,.2375],[8.,0.,.25607517888672315]],
        weak=False,assembly=False,reuse=reuse))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('runtime','inputs','small-checks','view'),required=True)
    p.add_argument('--case',default='two');p.add_argument('--reuse',action='store_true');a=p.parse_args()
    if a.mode=='view':view(a.root,a.case,a.reuse)
    else:globals()[a.mode.replace('-','_')](a.root)
