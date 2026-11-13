#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H8 fixed imports, all-state diagnosis, single new block and strict views."""
import argparse
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import time
import numpy as np
from numpy.testing import assert_allclose
from run_wj_h8 import load,save,sha,record_coverage_halt,coverage_quality
from study_wj_h5 import identity,compare_loads
from study_wj_h6 import configuration
from study_wj_h7 import clone,compatible,accounting
from wj_h7_analysis import Rows,context
from wj_h8_analysis import add_local,extend_witnesses,PressureProjection,local_rows,strict_context,FACE
from wj_h4_points import vector_hash,aligned,evaluate,owners
from wj_h3_analysis import outer_loads,principal
from wj3_measure import checkpoint_payload
from wj_h2_weak import independent_weak
from validate_wj1 import Probe

SOURCE=Path(__file__).resolve().parents[4]
H7=Path('/media/jay/新加卷1/nektar-v2026-runs/wj-h7-neutral-transport-20261113')
H6=H7.parent/'wj-h6-neutral-projection-20261112'

def checked(path,v):
    path=Path(path);a=identity(path);assert a['sha256']==v['sha256']and a['bytes']==v['bytes'],str(path);return a

def copy_checked(path,v,dst):
    checked(path,v);dst=Path(dst);dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,dst);checked(dst,v);return identity(dst)

def stage(root):
    """Selective imports only. No traversal of any old source/run root."""
    before=time.monotonic();old=load(SOURCE/'docs/plans/wj-h7-linux/runtime-equivalence.json');compiled=[]
    for v in old['compiled_inputs']:assert sha(SOURCE/v['path'])==v['sha256'];compiled.append(v)
    elf=[checked(v['path'],v)for v in old['qualified_ELFs']]
    runtime=root/'runtime';runtime.mkdir();private=load(SOURCE/'docs/plans/wj-h7-linux/private-runtime.json');glibs=[]
    for v in private['overrides'].values():
        dst=runtime/'glibc'/Path(v['path']).name;copy_checked(v['path'],v,dst);dst.chmod(0o755);glibs.append(identity(dst))
    geo=next(v for v in elf if Path(v['path']).name=='libSpatialDomains.so.5.9.0');copy_checked(geo['path'],geo,runtime/'geometry/libSpatialDomains.so.5.9.0')
    for alias in ('libSpatialDomains.so.5.9','libSpatialDomains.so'):(runtime/'geometry'/alias).symlink_to('libSpatialDomains.so.5.9.0')
    exes={};targets=load(SOURCE/'docs/plans/wj-h7-linux/closure/executable-map-qualification.json')['records']
    loader=runtime/'glibc/ld-linux-x86-64.so.2';library=f'{runtime}/glibc:{runtime}/geometry:/home/jay/.local/opt/nektar-v2026-mpi/lib64'
    for v in targets:
        source=Path(v['path']);assert sha(source)==v['sha256'];dst=runtime/(source.name+'.elf');shutil.copyfile(source,dst);assert sha(dst)==v['sha256'];dst.chmod(0o755)
        wrapper=runtime/source.name;wrapper.write_text('#!/bin/sh\nexec '+shlex.quote(str(loader))+' --library-path '+shlex.quote(library)+' '+shlex.quote(str(dst))+' "$@"\n');wrapper.chmod(0o755);exes[source.name]=dict(ELF=identity(dst),wrapper=identity(wrapper))
    save(root/'runtime-equivalence.json',dict(status='PASS',compiled_inputs=compiled,qualified_ELFs=elf,private_libraries=glibs,executables=exes,geometry=identity(runtime/'geometry/libSpatialDomains.so.5.9.0'),
        governing_producer='c344dc2209f030e426e2ee5501d5e5b72a6d2c46',observer_producer='e62b5f4d265f068dd5af3f090bfb0e385c3ae53e',native_rebuild=False,
        arithmetic='GNU11.4 C++17 Release -O3 -DNDEBUG -ffp-contract=off -fno-fast-math',source=identity(__file__)))
    gates=load(SOURCE/'docs/plans/wj-h7-linux/component-gates.json')['cases'];manifest=load(SOURCE/'docs/plans/wj-h7-linux/manifest.json');expected={v['path']:v for v in manifest['heavy']};histories=[];initials={}
    for h in (5.,5.5,6.):
        g=gates[str(h)];inv=g['inventory'];copy_checked(inv['path'],inv,root/f'imports/H{h:g}mm/inventory.json')
        parent=H6/'cases'/f'H{h:g}mm';rpath=parent/'case.evolution.json';v=next(v for v in load(SOURCE/'docs/plans/wj-h6-linux/manifest.json')['heavy']if v['path']==str(rpath));copy_checked(rpath,v,root/f'imports/H{h:g}mm/initial.evolution.json')
        r0=load(rpath);initials[str(h)]=dict(report=identity(rpath),state_sha256=vector_hash(r0['state']),reference_sha256=vector_hash(r0['fixed_reference']),dt=r0['dt'])
        for branch in ('coarse','fine'):
            rp=Path(g['native'][branch]['report']['path']);v=g['native'][branch]['report'];out=root/f'imports/H{h:g}mm'/branch
            report=copy_checked(rp,v,out/'case.evolution.json');xml=copy_checked(rp.with_name('case.xml'),expected[str(rp.with_name('case.xml'))],out/'case.xml');r=load(rp)
            assert not r['failure']and len(r['state'])==2328;assert vector_hash(r['fixed_reference'])==vector_hash(r0['fixed_reference'])
            snapshots=[]
            for a in r['pimple_attempts']:
                p=Path(a['k_snapshot']);assert sha(p)==a['k_snapshot_sha256'];snapshot=load(p);assert vector_hash(snapshot['raw'])==vector_hash(snapshot['corrected']);snapshots.append(dict(packet=identity(p),fields_present=list(snapshot),predictor_pressure_separation='UNRESOLVED:no separate predictor/correction vectors'))
            checkpoints=[]
            for step in range(r['steps']+1):
                ctx=context(r,step)
                if step:
                    cp=rp.parent/f'CP.{step}';v=expected[str(cp/'State.bin')]
                elif branch=='coarse':cp=parent/'CP';v=next(v for v in load(SOURCE/'docs/plans/wj-h6-linux/manifest.json')['heavy']if v['path']==str(cp/'State.bin'))
                else:
                    original_intent=load(rp.with_name('intent.json'));cp=Path(original_intent['checkpoint']['path']).parent;v=original_intent['checkpoint']
                copy_checked(cp/'State.bin',v,out/f'CP.{step}/State.bin');payload=checkpoint_payload(out/f'CP.{step}')
                assert payload['steps']==step and payload['dt']==r['dt']and payload['signature']==r['signature']and payload['time']==ctx['time']
                assert vector_hash(payload['state'])==vector_hash(ctx['state'])and payload['reference_sha256']==vector_hash(r['fixed_reference'])
                checkpoints.append(dict(step=step,copy=identity(out/f'CP.{step}/State.bin'),original=identity(cp/'State.bin'),state_sha256=vector_hash(ctx['state']),callback_history_sha256=vector_hash(ctx['history']),alpha=ctx['alpha'],time=ctx['time'],native_history_sha256=[vector_hash(v)for v in payload['native_history']],native_history_times=payload['native_history_times']))
            histories.append(dict(height=h,branch=branch,report=report,original_report=identity(rp),xml=xml,checkpoints=checkpoints,snapshot_inventory=snapshots))
    save(root/'imports.json',dict(histories=histories,initials=initials,original_H7_allocation='closed1800 and107.85679090699705 unknown charge unchanged',seconds=time.monotonic()-before))
    # Frozen old XML bytes, not a regenerated/rounded neutral profile.
    xml=H6/'cases/H6mm/case.xml';expected_h6=load(SOURCE/'docs/plans/wj-h6-linux/manifest.json')['heavy'];v=next(v for v in expected_h6 if v['path']==str(xml));copy_checked(xml,v,root/'imports/H6mm/reference.xml')
    save(root/'stage-receipt.json',dict(status='PASS',source=identity(__file__),imports=identity(root/'imports.json'),runtime=identity(root/'runtime-equivalence.json'),seconds=time.monotonic()-before))
    print('PASS selective H6/H7 imports, native histories/checkpoints and exact private runtime',flush=True)


def selected(root,h):
    return [v for v in load(root/'imports.json')['histories']if v['height']==h]

def inventory(root,h,include_new=False):
    items=selected(root,h);r=load(items[0]['report']['path']);c=configuration(h);inv=add_local(r,load(root/f'imports/H{h:g}mm/inventory.json'),c)
    if include_new:
        rp=root/'case/run/case.evolution.json'
        if rp.exists():inv=extend_witnesses(load(rp),inv,c)
    return inv


def diagnose(root,h,version,include_new=False):
    out=root/'diagnosis'/f'H{h:g}mm'/version;out.mkdir(parents=True);items=selected(root,h);c=configuration(h);r0=load(root/f'imports/H{h:g}mm/initial.evolution.json');inv=inventory(root,h,include_new)
    save(out/'inventory.json',inv);cache=Rows(r0,[v['physical']for v in inv['points']],c);projection=PressureProjection(cache);initial=cache.evaluate(r0['state']);records=[];local=[];arrays=[]
    keys={tuple(v)for v in inv['local_keys']};local_idx=[i for i,(p,o)in enumerate(cache.rows[:cache.point_count])if tuple(float(v).hex()for v in p)in keys]
    if include_new and(root/'case/run/case.evolution.json').exists():items=items+[dict(height=h,branch='quarter',report=identity(root/'case/run/case.evolution.json'))]
    probe=Probe(str(root/'runtime/WJKernelProbe'))
    try:
        for item in items:
            r=load(item['report']['path']);branch=item['branch'];ctxs=[context(r,s)for s in range(r['steps']+1-len(r['accepted_states']),r['steps']+1)]
            if r.get('physical_range_events'):ctxs.append(context(r,rejected=True))
            for ctx in ctxs:
                state=ctx['state'];before=vector_hash(state);data=cache.evaluate(state);summary=cache.summary(data,initial);pd=projection.measure(data)
                # Independent and compiled replay of the actual expanded minimum.
                wi=summary['witness'];G=np.array(wi['G']);p=principal(wi['k'],wi['epsilon'],G);H=np.outer(p['polarization'],p['xi']);v=probe.closure(wi['k'],wi['epsilon'],G);d=probe.closure(wi['k'],wi['epsilon'],G,(0,0,H));compiled=-np.sum(H*d['R'])/v['scalar'][1]
                assert_allclose(compiled,summary['min_in_plane'],rtol=0,atol=2e-11)
                assert vector_hash(state)==before
                key=f'{branch}-{ctx["role"]}-{ctx["step"]}';file=out/(key+'.npz');np.savez_compressed(file,**data)
                lr=local_rows(cache,data,local_idx);local.extend(dict(branch=branch,step=ctx['step'],time=ctx['time'],role=ctx['role'],**v)for v in lr)
                cp=next((v for v in item.get('checkpoints',[])if v['step']==ctx['step']),None)
                if branch=='quarter'and ctx['role']=='accepted':
                    from measure_wj_h2 import checkpoint_for
                    cpdir,payload=(str(root/'case/cp0/CP'),checkpoint_payload(root/'case/cp0/CP'))if ctx['step']==0 else checkpoint_for(dict(report=item['report']['path'],xml=str(Path(item['report']['path']).with_name('case.xml'))),ctx['step']);cp=dict(copy=identity(Path(cpdir)/'State.bin'),native_history_times=payload['native_history_times'])
                record=dict(branch=branch,step=ctx['step'],time=ctx['time'],alpha=ctx['alpha'],role=ctx['role'],producer=item['report'],checkpoint=cp,state_sha256=before,callback_history_sha256=vector_hash(ctx['history']),
                    compiled_minimum=float(compiled),inventory_sha256=inv['sha256'],arrays=identity(file),pressure_projection=pd,**summary)
                records.append(record)
                if ctx['role']=='accepted'and summary['status']in ('OBSTRUCTED','UNRESOLVED'):record_coverage_halt(root/'accounting.json',dict(trajectory=f'{h:g}/{branch}',height=h,branch=branch,role='accepted',time=ctx['time'],step=ctx['step'],status=summary['status'],witness=wi))
            save(out/(branch+'-counts.json'),accounting(r))
        save(out/'records.json',records);save(out/'local.json',local)
        status=coverage_quality(records);save(out/'summary.json',dict(status=status,records=identity(out/'records.json'),local=identity(out/'local.json'),inventory=identity(out/'inventory.json'),factorizations=projection.factorizations,
            diagnostic_mass_solves=sum(v['pressure_projection']['mass_solves']for v in records),native_governing_solves=0,causal_predictor_pressure_attribution='UNRESOLVED:no genuine separated arithmetic vectors in retained H7 snapshots',
            coverage_halt=load(root/'accounting.json').get('coverage_halt'),coverage_earliest=load(root/'accounting.json').get('coverage_earliest'),source=identity(__file__)))
        save(root/f'diagnosis-view-H{h:g}.json',identity(out/'summary.json'));print('H8 diagnosis',h,version,status,'states',len(records),flush=True)
    finally:probe.close()


def prepare(root):
    folder=root/'case/cp0';folder.mkdir(parents=True);xml=clone(root/'imports/H6mm/reference.xml',folder/'case.xml',.0025,0)
    save(root/'case/intent.json',dict(xml=identity(xml),physical_parent=identity(root/'imports/H6mm/reference.xml'),dt=.0025,steps=16,H_w_mm=6.,expected_size=2328,NumSteps_CP0=0))


def pair_initial(root):
    r=load(root/'case/cp0/case.evolution.json');base=load(root/'imports/H6mm/initial.evolution.json');compatible(base,r,.0025,same_dt=False);assert not r['failure']
    for a,b in zip(base['native_history'],r['native_history']):assert vector_hash(a)==vector_hash(b)
    assert len(r['native_history'])==len(base['native_history'])
    p=checkpoint_payload(r['checkpoint']);assert p['reference_sha256']==vector_hash(base['fixed_reference'])and vector_hash(p['state'])==vector_hash(r['state'])and p['signature']==r['signature']and p['dt']==.0025
    folder=root/'case/run';folder.mkdir();xml=clone(root/'case/cp0/case.xml',folder/'case.xml',.0025,16,r['checkpoint'])
    save(root/'case/initial-pair.json',dict(status='PASS',coupled_state_sha256=vector_hash(r['state']),fixed_reference_sha256=vector_hash(r['fixed_reference']),physical_histories_equal=True,source_cp0=identity(root/'case/cp0/case.evolution.json'),checkpoint=identity(Path(r['checkpoint'])/'State.bin'),run_xml=identity(xml),distinct_history_labels=[base['native_history_times'],r['native_history_times']]))


def view_request(root,branch,step):
    if branch=='quarter':rp=root/'case/run/case.evolution.json';r=load(rp);from measure_wj_h2 import checkpoint_for;cp,payload=checkpoint_for(dict(report=str(rp),xml=str(rp.with_name('case.xml'))),step)
    else:
        item=next(v for v in selected(root,6.)if v['branch']==branch);rp=Path(item['report']['path']);r=load(rp);cp=str(Path(next(v for v in item['checkpoints']if v['step']==step)['copy']['path']).parent)
    ctx=strict_context(r,step);out=root/'views'/f'{branch}-{step}';out.mkdir(parents=True)
    packet=dict(operation='open_checks',inner_stress_model='Linear',native_report=str(rp),native_report_sha256=sha(rp),native_checkpoint=cp,checkpoint_sha256=sha(Path(cp)/'State.bin'),allow_failed_producer=bool(r['failure']),physical_observations=[[8.,0.,FACE]],weak=True,weak_outer_order=32,assembly=False,reuse=False)
    save(out/'request.json',packet);save(out/'context.json',dict(producer=identity(rp),xml=identity(rp.with_name('case.xml')),step=step,time=ctx['time'],alpha=ctx['alpha'],state_sha256=vector_hash(ctx['state']),callback_history_sha256=vector_hash(ctx['history']),checkpoint=identity(Path(cp)/'State.bin')))


def finish_view(root,branch,step):
    out=root/'views'/f'{branch}-{step}';ctx=load(out/'context.json');request=load(out/'request.json');r=load(ctx['producer']['path']);p=load(out/'native.json');state,h,alpha,t,packet=__import__('wj_h3_analysis').state_context(r,step)
    assert sha(ctx['producer']['path'])==ctx['producer']['sha256']and sha(ctx['checkpoint']['path'])==ctx['checkpoint']['sha256']
    assert p['initialization_work']==[0,0]and p['governing_solves']==0 and p['selected_absolute_step']==step
    assert p['state_sha256']==ctx['state_sha256']==vector_hash(state)and p['callback_history_sha256']==ctx['callback_history_sha256']==vector_hash(h)
    weak=independent_weak(p['weak_checks'],dict(r,open_boundary=p['open_boundary']),enforce=False)
    loads={q:outer_loads(r,state,h,alpha,q)for q in (12,24,32)};q={f'{a}-32':compare_loads(loads[a],loads[32])for a in (12,24)}
    save(out/'outer-loads.json',loads)
    save(out/'result.json',dict(status='PASS'if weak['status']=='PASS'and all(v['status']=='PASS'for v in q.values())else 'UNRESOLVED',context=ctx,native=identity(out/'native.json'),independent_weak=weak,quadrature=q,norms=p['physical_mass_dual_norms'],open_boundary={k:v for k,v in p['open_boundary'].items()if k not in ('boundary_rows','stations')},role='accepted predecessor diagnostic'if r['failure']else 'genuine accepted completed stage'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',required=True,choices=('stage','diagnose','prepare','pair','view-request','view-finish'));p.add_argument('--height',type=float,choices=(5.,5.5,6.));p.add_argument('--version',default='v1');p.add_argument('--new',action='store_true');p.add_argument('--branch',choices=('coarse','fine','quarter'));p.add_argument('--step',type=int);a=p.parse_args()
    if a.mode=='stage':stage(a.root)
    elif a.mode=='diagnose':diagnose(a.root,a.height,a.version,a.new)
    elif a.mode=='prepare':prepare(a.root)
    elif a.mode=='pair':pair_initial(a.root)
    elif a.mode=='view-request':view_request(a.root,a.branch,a.step)
    else:finish_view(a.root,a.branch,a.step)
