#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H7 fixed inputs, native histories, every-state scans and strict endpoints."""
import argparse
import json
import math
from pathlib import Path
import subprocess
import time
import xml.etree.ElementTree as ET
import numpy as np
from numpy.testing import assert_allclose
from run_wj_h7 import execute,load,save,sha,MATRIX,record_coverage_halt,coverage_quality
from study_wj_h5 import identity,compare_loads
from study_wj_h6 import configuration
from wj_h6_projection import Functions
from wj_h4_points import owners,evaluate,vector_hash,aligned
from wj_h3_analysis import state_context,principal,outer_loads
from wj_h7_analysis import Rows,context
from wj_h2_weak import independent_weak
from measure_wj_h2 import checkpoint_for
from wj3_measure import checkpoint_payload
from cls_epawf_m3_cases import set_parameter,set_info,write
from validate_wj1 import Probe

PARENT='wj-h6-neutral-projection-20261112'

def paths(root):
    base=root/'runtime'
    if not (root/'runtime-equivalence.json').exists():raise ValueError('verify exact H6 private runtime closure first')
    return tuple(base/n for n in ('IncNavierStokesSolver','RANSWallResolutionProbe','WJKernelProbe'))

def case(root,h,branch):return root/'cases'/f'H{h:g}mm'/branch

def clone(src,dst,dt,steps,restart=None):
    tree=ET.parse(src).getroot();set_parameter(tree,'TimeStep',dt);set_parameter(tree,'NumSteps',steps)
    set_info(tree,'RANSWallCheckpoint',str(dst.parent/'CP'))
    if restart:set_info(tree,'RANSWallRestart',str(restart))
    write(tree,dst)
    # Enforce exact semantic differences and byte-identical Function elements.
    old=ET.parse(src).getroot();changed=ET.parse(dst).getroot()
    for x in (old,changed):
        for p in x.findall('./CONDITIONS/PARAMETERS/P'):
            if p.text.split('=')[0].strip()in ('TimeStep','NumSteps'):x.find('./CONDITIONS/PARAMETERS').remove(p)
        for e in list(x.findall('./CONDITIONS/SOLVERINFO/I')):
            if e.get('PROPERTY')in ('RANSWallCheckpoint','RANSWallRestart'):x.find('./CONDITIONS/SOLVERINFO').remove(e)
    assert ET.tostring(old)==ET.tostring(changed),'unpermitted reference/mesh/setting change'
    return dst


def compatible(a,b,dt,*,same_dt=True):
    if b['dt']!=dt:raise ValueError('checkpoint timestep differs; no cross-dt relabeling')
    for key in ('state','fixed_reference'):
        if np.asarray(a[key],dtype='<f8').tobytes()!=np.asarray(b[key],dtype='<f8').tobytes():raise ValueError(key+' bytes differ')
    if same_dt and a['signature']!=b['signature']:raise ValueError('physical signature differs')
    if b['steps']!=0 or b['time']!=0:raise ValueError('CP0 required')


def verify_runtime(root):
    source=Path(__file__).resolve().parents[4];sealed=source/'docs/plans/wj-h2-linux'
    checks=[]
    for v in load(sealed/'native-source-inputs.json'):
        assert sha(source/v['path'])==v['sha256'],'changed native source input'
        checks.append(v)
    private=load(root/'private-runtime.json');overrides=private['overrides'];elf=[]
    for v in load(sealed/'loaded-runtime.json')['loaded_files']:
        original=v['path'].replace('/media/jay/新加卷/','/media/jay/新加卷1/',1)
        path=Path(overrides.get(original,{}).get('path',original))
        if not path.is_file()or path.open('rb').read(4)!=b'\x7fELF':continue
        value=identity(path);assert value['sha256']==v['sha256'],'changed qualified ELF: '+original
        elf.append(dict(value,original=original))
    assert len(checks)==340 and len(elf)==133
    base=root.parent/'wj-h2-20261111';executables={}
    for role,path,expected in (
        ('governing',base/'runtime-v4/IncNavierStokesSolver','9ae37a1aba62755b10ad99c9484c498cef9d927a75fb950f00344f4ffc5d06fe'),
        ('observer',base/'observer-v6/RANSWallResolutionProbe','6f35b837f0dd13eb155e0fe34f0dc670c609ee6cbe739609a5de288dd49ed06d'),
        ('kernel',base/'build-v1/WJKernelProbe','3d4edf79a2718a281175be1a575776a33dd5f988c48d5a1b594a9b5f2dd7844c')):
        executables[role]=identity(path);assert executables[role]['sha256']==expected,'changed H6 executable '+role
    save(root/'runtime-equivalence.json',dict(status='PASS',compiled_inputs=checks,qualified_ELFs=elf,executables=executables,
        private_runtime=identity(root/'private-runtime.json'),generating_source=identity(__file__),
        governing_producer='c344dc2209f030e426e2ee5501d5e5b72a6d2c46',observer_producer='e62b5f4d265f068dd5af3f090bfb0e385c3ae53e',
        native_rebuild=False,arithmetic_policy='GNU11.4 C++17 Release -O3 -DNDEBUG -ffp-contract=off -fno-fast-math; RANS_BUILD_SOLVER=ON',
        reason='six host glibc files upgraded; private original-version files exactly match every H6 hash'))


def imports(root):
    verify_runtime(root);parent=root.parent/PARENT;manifest=load(Path(__file__).resolve().parents[4]/'docs/plans/wj-h6-linux/manifest.json')
    expected={v['path']:v for v in manifest['heavy']};records=[]
    for h in (5.5,6.,5.):
        folder=parent/'cases'/f'H{h:g}mm';rp=folder/'case.evolution.json';r=load(rp);cp=Path(r['checkpoint']);payload=checkpoint_payload(cp)
        for p in (rp,folder/'case.xml',cp/'State.bin'):assert identity(p)==expected[str(p)]
        assert not r['failure']and r['dt']==.01 and r['steps']==0 and len(r['state'])==2328
        compatible(r,r,.01);assert payload['signature']==r['signature'] and vector_hash(payload['state'])==vector_hash(r['state'])
        assert payload['reference_sha256']==vector_hash(r['fixed_reference'])
        assert r['native_history']==payload['native_history']and r['native_history_times']==payload['native_history_times']
        records.append(dict(height=h,report=identity(rp),xml=identity(folder/'case.xml'),checkpoint=identity(cp/'State.bin'),
            cp_directory=str(cp),state_sha256=vector_hash(r['state']),reference_sha256=vector_hash(r['fixed_reference']),
            history_hashes=[vector_hash(v)for v in r['native_history']],history_times=r['native_history_times'],signature=r['signature']))
    save(root/'imports.json',records);print('PASS three immutable H6 .01 CP0 imports; no governing work',flush=True)


def prepare(root,h):
    record=next(v for v in load(root/'imports.json')if v['height']==h)
    for branch,dt,n in (('coarse',.01,4),('fine-cp0',.005,0)):
        folder=case(root,h,branch);folder.mkdir(parents=True)
        xml=clone(record['xml']['path'],folder/'case.xml',dt,n,record['cp_directory']if branch=='coarse'else None)
        save(folder/'intent.json',dict(height=h,dt=dt,steps=n,xml=identity(xml),parent=record,reference='unchanged H6'))


def initial_pair(root,h):
    rec=next(v for v in load(root/'imports.json')if v['height']==h);old=load(rec['report']['path']);folder=case(root,h,'fine-cp0');r=load(folder/'case.evolution.json')
    assert not r['failure']and r['steps']==0
    cp=checkpoint_payload(r['checkpoint']);compatible(old,r,.005,same_dt=False)
    assert cp['signature']==r['signature'] and cp['reference_sha256']==vector_hash(old['fixed_reference']) and vector_hash(cp['state'])==vector_hash(r['state'])
    assert np.asarray(r['native_history']) .shape==np.asarray(old['native_history']).shape
    for a,b in zip(old['native_history'],r['native_history']):assert vector_hash(a)==vector_hash(b),'physical history differs'
    assert cp['native_history']==r['native_history']and cp['native_history_times']==r['native_history_times']
    out=case(root,h,'fine');out.mkdir();xml=clone(folder/'case.xml',out/'case.xml',.005,8,r['checkpoint'])
    save(out/'intent.json',dict(xml=identity(xml),cp0=identity(folder/'case.evolution.json'),checkpoint=identity(Path(r['checkpoint'])/'State.bin'),
        bitwise_state_and_fixed_reference=True,physical_history_equal=True,coarse_history_times=old['native_history_times'],fine_history_times=r['native_history_times']))
    print('PASS bitwise same initial physical fields/reference, distinct valid dt identities',h,flush=True)


def reports(root,h):
    record=next(v for v in load(root/'imports.json')if v['height']==h);out=[('H6',load(record['report']['path']),Path(record['report']['path']))]
    for b in ('coarse','fine'):
        p=case(root,h,b)/'case.evolution.json'
        if p.exists():out.append((b,load(p),p))
    return out


def inventory(root,h,rs):
    parent=root.parent/PARENT/'comparisons'/f'H{h:g}mm'/'inventory.json';base=load(parent);points={tuple(float(v).hex()for v in p['physical']):p for p in base['points']};c=configuration(h)
    def add(p,tag):
        if not(0<=p[0]<=8 and c.d0<=p[2]<=c.top):return
        key=tuple(float(v).hex()for v in p);points.setdefault(key,dict(physical=list(p),tags=[]))['tags'].append(tag)
    for name,r,rp in rs:
        packets=r['wj_audit_packets']+[e['wj_packet']for e in r.get('physical_range_events',[])if e.get('wj_packet')]
        for packet in packets:
            for w in packet['principal_witnesses'].values():
                p=w['position'];add(p,'new_native_witness')
                for axis in (0,2):
                    for sign in (-1,1):q=p.copy();q[axis]+=sign*1e-6;add(q,'new_native_offset')
    from wj_h4_points import digest
    result=dict(parent=identity(parent),points=list(points.values()),rule='H6 full physical union + all new native limiting points and prescribed neighborhoods/all owners')
    result['sha256']=digest(result);return result


def accounting(r):
    counts={k:sum(a[k]for a in r['pimple_attempts'])for k in ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls','bound_activations')}
    for a in r['pimple_attempts']:
        if a['scheduled_complete']:
            assert all(a[k]==1 for k in ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls'))
        assert a['bound_activations']==0,'material correction; trajectory cannot qualify'
        if a.get('k_snapshot'):
            p=Path(a['k_snapshot']);assert sha(p)==a['k_snapshot_sha256'];s=load(p)
            assert vector_hash(s['raw'])==vector_hash(s['corrected']),'raw/corrected k differs'
    return dict(counts,attempted=len(r['pimple_attempts']),accepted=r['steps'],failure=r['failure'],
        frozen_blocks=[{k:v for k,v in a.items()if 'residual'in k or k in ('time','scheduled_complete','failure','bound_calls','bound_activations')}for a in r['pimple_attempts']])


def scan(root,h,version):
    rs=reports(root,h);c=configuration(h);inv=inventory(root,h,rs);out=root/'scans'/f'H{h:g}mm'/version;out.mkdir(parents=True)
    save(out/'inventory.json',inv);r0=rs[0][1];cache=Rows(r0,[p['physical']for p in inv['points']],c);initial=cache.evaluate(r0['state']);records=[];saved={};cross=[]
    kernel=Probe(str(paths(root)[2]))
    try:
        for name,r,rp in rs:
            if name!='H6':
                assert r['initialization_lift_calls_iterations']==[0,0],'restart performed Project/Lift'
                assert vector_hash(r['fixed_reference'])==vector_hash(r0['fixed_reference'])
            ctxs=[context(r,s)for s in range(r['steps']+1-len(r['accepted_states']),r['steps']+1)]
            if r['physical_range_events']:ctxs.append(context(r,rejected=True))
            for ctx in ctxs:
                state=ctx['state'];data=cache.evaluate(state);summary=cache.summary(data,initial);key=f'{name}-{ctx["role"]}-{ctx["step"]}'
                np.savez_compressed(out/(key+'.npz'),**data);saved[(name,ctx['step'],ctx['role'])]=key
                packet=ctx['packet'];native_replays={}
                for tag,w in packet['principal_witnesses'].items():
                    own=next(o for o in owners(r,w['position'])if o['geometry']==w['geometry']);v=evaluate(r,state,w['position'],own)
                    G=np.array(v['G']);H=np.array(w['H_row_major']).reshape(3,3);val=kernel.closure(v['k'],v['epsilon'],G);d=kernel.closure(v['k'],v['epsilon'],G,(0,0,H));m=-np.sum(H*d['R'])/val['scalar'][1]
                    assert_allclose(m,w['minimum_over_nu_eff'],rtol=0,atol=2e-11);assert_allclose(v['principal'][tag],w['minimum_over_nu_eff'],rtol=0,atol=2e-11)
                    # Check batched arithmetic against the established owned-point path.
                    ix=next(i for i,(p,o)in enumerate(cache.rows)if o['geometry']==w['geometry']and np.max(abs(np.array(p)-w['position']))<1e-14)
                    assert_allclose(data['G'][ix],v['G'],rtol=0,atol=2e-11);assert_allclose(data['margin'][ix],v['principal']['in_plane'],rtol=0,atol=2e-11)
                    native_replays[tag]=dict(native=w,compiled=float(m),independent=v['principal'])
                wi=summary['witness'];G=np.array(wi['G']);mi=principal(wi['k'],wi['epsilon'],G);xi=np.array(mi['xi']);a=np.array(mi['polarization']);H=np.outer(a,xi)
                val=kernel.closure(wi['k'],wi['epsilon'],G);act=kernel.closure(wi['k'],wi['epsilon'],G,(0,0,H));compiled=-np.sum(H*act['R'])/val['scalar'][1]
                assert_allclose(compiled,summary['min_in_plane'],rtol=0,atol=2e-11)
                record=dict(branch=name,**{k:ctx[k]for k in ('step','time','alpha','role')},producer=identity(rp),state_sha256=vector_hash(state),
                    history_sha256=vector_hash(ctx['history']),packet_state_sha256=packet['state_sha256'],native_witness_replays=native_replays,
                    expanded_compiled_replay=float(compiled),inventory_sha256=inv['sha256'],arrays=identity(out/(key+'.npz')),**summary)
                records.append(record)
                if name!='H6'and ctx['role']=='accepted'and summary['status']in ('OBSTRUCTED','UNRESOLVED'):
                    record_coverage_halt(root/'accounting.json',dict(height=h,branch=name,trajectory=f'{h:g}/{name}',role='accepted',step=ctx['step'],time=ctx['time'],status=summary['status'],witness=wi))
            if name!='H6':save(out/(name+'-accounting.json'),accounting(r))
        pairs=[]
        for step in range(5):
            ka=('coarse',step,'accepted');kb=('fine',step*2,'accepted')
            if ka not in saved or kb not in saved:pairs.append(dict(coarse_step=step,status='NOT_RUN'));continue
            a=next(v for v in records if v['branch']=='coarse'and v['step']==step and v['role']=='accepted');b=next(v for v in records if v['branch']=='fine'and v['step']==step*2 and v['role']=='accepted')
            assert aligned(step,.01,a['time'])==aligned(2*step,.005,b['time'])
            with np.load(out/(saved[ka]+'.npz'))as aa,np.load(out/(saved[kb]+'.npz'))as bb:
                metrics=cache.norms(aa,bb);active=cache.outer & (np.arange(len(cache.rows))<cache.point_count);delta=bb['margin'][active]-aa['margin'][active]
                i=np.flatnonzero(active)[np.argmax(abs(delta))];pairs.append(dict(coarse_step=step,time=step*.01,status='PASS',recorded_times=[a['time'],b['time']],
                    difference=metrics,margin_max_same_point_change=float(max(abs(delta))),margin_at_max_change=[float(aa['margin'][i]),float(bb['margin'][i])],
                    position=cache.rows[i],minima=[a['witness'],b['witness']]))
        save(out/'records.json',records);save(out/'pairs.json',pairs)
        save(out/'summary.json',dict(inventory=identity(out/'inventory.json'),records=identity(out/'records.json'),pairs=identity(out/'pairs.json'),
            peak_quality_status=coverage_quality(records)['quality'],expanded_status=coverage_quality(records)['coverage'],coverage_halt=load(root/'accounting.json').get('coverage_halt'),
            source=identity(__file__),scope='retrospective expanded inventory; native decisions preserved'))
        save(root/f'scan-view-H{h:g}.json',identity(out/'summary.json'));print('PASS completed every-state scan',h,version,[(r['branch'],r['step'],r['role'],r['status'],r['min_in_plane'])for r in records],flush=True)
    finally:kernel.close()


def observe(root,h,branch):
    rp=case(root,h,branch)/'case.evolution.json';r=load(rp);out=root/'observations'/f'H{h:g}mm'/branch;out.mkdir(parents=True)
    if r['steps']==0:
        save(out/'result.json',dict(status='NOT_RUN',reason='no accepted governing step; qualified H6 CP0 view retained',completed_stage='NOT_APPLICABLE'));return
    cp,payload=checkpoint_for(dict(report=str(rp),xml=str(rp.with_name('case.xml'))),r['steps']);before={str(p):sha(p)for p in (rp,Path(cp)/'State.bin')}
    request=dict(operation='open_checks',inner_stress_model='Linear',native_report=str(rp),native_report_sha256=sha(rp),native_checkpoint=cp,checkpoint_sha256=sha(Path(cp)/'State.bin'),
        allow_failed_producer=bool(r['failure']),physical_observations=[],weak=True,weak_outer_order=32,assembly=False,reuse=False)
    save(out/'request.json',request);start=time.monotonic()
    with(out/'native.log').open('x')as f:result=subprocess.run([str(paths(root)[1]),str(rp.with_name('case.xml')),str(out/'request.json'),str(out/'native.json')],stdout=f,stderr=subprocess.STDOUT,check=True)
    p=load(out/'native.json');ctx=context(r);assert p['initialization_work']==[0,0]and p['governing_solves']==0
    assert p['state_sha256']==vector_hash(ctx['state'])and p['callback_history_sha256']==vector_hash(ctx['history'])and p['selected_absolute_step']==r['steps']
    assert all(sha(k)==v for k,v in before.items())
    weak=independent_weak(p['weak_checks'],dict(r,open_boundary=p['open_boundary']),enforce=False)
    views={q:outer_loads(r,ctx['state'],ctx['history'],ctx['alpha'],q)for q in (12,24,32)};quad={f'{q}-32':compare_loads(views[q],views[32])for q in (12,24)}
    save(out/'outer-loads.json',views)
    save(out/'result.json',dict(status='PASS'if weak['status']=='PASS'and all(v['status']=='PASS'for v in quad.values())else 'UNRESOLVED',
        role='diagnostic accepted predecessor of failed producer'if r['failure']else 'terminal accepted completed stage',producer=identity(rp),checkpoint=identity(Path(cp)/'State.bin'),native=identity(out/'native.json'),
        state_sha256=p['state_sha256'],callback_history_sha256=p['callback_history_sha256'],time=ctx['time'],alpha=ctx['alpha'],strict_seconds=time.monotonic()-start,
        spatial_and_completed_norms=p['physical_mass_dual_norms'],independent_weak=weak,quadrature=quad,open_boundary=p['open_boundary'],
        scalar_condition=p['scalar_domain']['admissible']if 'admissible'in p['scalar_domain']else None))
    status=load(out/'result.json')['status']
    print(status+' strict endpoint/predecessor lineage, weak and fixed quadrature',h,branch,flush=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',required=True,choices=('import','prepare','initial-pair','scan','observe'))
    p.add_argument('--height',type=float,choices=(5.,5.5,6.));p.add_argument('--branch',choices=('coarse','fine'));p.add_argument('--version',default='v1');a=p.parse_args()
    if a.mode=='import':imports(a.root)
    elif a.mode=='prepare':prepare(a.root,a.height)
    elif a.mode=='initial-pair':initial_pair(a.root,a.height)
    elif a.mode=='scan':scan(a.root,a.height,a.version)
    else:observe(a.root,a.height,a.branch)
if __name__=='__main__':main()
