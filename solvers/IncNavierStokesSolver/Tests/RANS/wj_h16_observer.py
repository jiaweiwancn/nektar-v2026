#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H16 identity/packet adapter; H15 numerical assessment is unchanged."""
import argparse
import hashlib
import importlib
import json
import os
from pathlib import Path
import time
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
import wj_h15_observer as previous
from wj_h15_observer import identity,verify,atomic,vector_hash

OLD=Path('/home/jay/build/wj-h15-persistent-neutral-20261114')
ORIGINAL_LOADER=previous.load_cache
ORIGINAL_VERIFY=previous.verify

def active_source_bindings(historical):
    records=[]
    for old in historical['source_files']:
        module=importlib.import_module(Path(old['path']).stem);active=identity(module.__file__)
        assert active['bytes']==old['bytes']and active['sha256']==old['sha256'],'active numerical observer source differs: '+module.__name__
        records.append(dict(historical=old,active=active))
    return records

def compatibility(root,manifest,destination):
    m=load(manifest);assert m['status']=='PASS'
    for v in m['source_files']:ORIGINAL_VERIFY(v)
    sources=active_source_bindings(m)
    for k in ('profile','native_runtime_manifest'):ORIGINAL_VERIFY(m['reference_identity'][k])
    cache=ORIGINAL_VERIFY(m['cache'])
    if cache.name=='immutable-cache.pickle'and cache.stat().st_size==211582940:
        assert m['cache']['sha256']=='800ac28c73ecb6498e6504b4e14e774db2cff893f2c0aaf418cb7036d50d00b0'
    runtime=load(root/'science/runtime-build-v1.json')
    for v in runtime['executables']:verify(v['ELF']);verify(v['wrapper'])
    dependencies=load(root/'science/runtime-dependencies-v1.json');assert dependencies['status']=='PASS'
    for v in dependencies['libraries']:verify(v['actual'])
    record=dict(status='PASS',schema='WJH16CacheCompatibilityV1',historical_manifest=identity(manifest),historical_cache=m['cache'],historical_reference_identity=m['reference_identity'],active_sources=sources,
        active_adapter=identity(__file__),new_native_identity=identity(root/'science/runtime-build-v1.json'),runtime_dependencies=identity(root/'science/runtime-dependencies-v1.json'),
        geometry_space_signs_quadrature='unchanged cached report; checked against native structure on every completed packet',complete_physics_profile=m['reference_identity'],
        permitted_difference='native monitor/serialization/directory durability only; no mathematical row/profile/model change',historical_pins_modified=False,actual_imports_verified=True)
    atomic(destination,record)
    return record

def load_compatible(record_path):
    c=load(record_path);assert c['status']=='PASS'and c['actual_imports_verified']
    verify(c['active_adapter']);verify(c['new_native_identity']);verify(c['runtime_dependencies'])
    m=load(verify(c['historical_manifest']));sources=active_source_bindings(m)
    assert sources==c['active_sources']
    for v in load(c['runtime_dependencies']['path'])['libraries']:verify(v['actual'])
    return ORIGINAL_LOADER(verify(c['historical_manifest']))

def packet_fields(r):
    assert r['schema']=='WJCompletedStepV1'
    assert len(r['stage_audits'])==len(r['step_audits'])==1
    assert r['stage_audits'][0][0]==r['step_audits'][0][0]==r['time']
    assert r['stage_audits'][0][2]==r['last_stage_alpha']==1/r['dt']
    rows=r['explicit_audits'];assert 0<len(rows)<=2
    assert r['first_retained_explicit_call']+len(rows)-1==r['actual_explicit_calls']
    assert all(len(v)==2 and v[0]<=r['time']for v in rows)
    assert r['native_hard_audit']['admissible'][0]==1
    assert len(r['accepted_states'])==len(r['wj_audit_packets'])==1
    return dict(stage_rows=1,step_rows=1,explicit_rows=len(rows),actual_explicit_calls=r['actual_explicit_calls'],first_retained_explicit_call=r['first_retained_explicit_call'],hard_audit_keys=sorted(r['native_hard_audit']))

def daemon(root,intent,compat,behavior):
    # Thin hooks validate active imports/new native identity and new packet
    # serialization before the unchanged scientific evaluator can issue PASS.
    def loader(unused):
        start=time.monotonic();cache,m=load_compatible(compat)
        directory=Path(load(intent)['directory'])
        atomic(directory/'h16-cache-load.json',dict(seconds=time.monotonic()-start,compatibility=identity(compat),actual_imports_verified=True,immutable_row_rebuilds=0,base_rows=len(cache.base.rows)))
        return cache,m
    def pinned(v):
        p=ORIGINAL_VERIFY(v)
        if p.name.endswith('.report.json')and p.parent==Path(load(intent)['directory']):
            r=load(p);fields=packet_fields(r);dest=p.parent/(p.stem+'.retention-fields.json')
            if not dest.exists():atomic(dest,dict(source=v,fields=fields,status='PASS',schema='WJH16PacketFieldsV1'))
        return p
    previous.load_cache=loader;previous.verify=pinned
    previous.daemon(root,intent,compat,behavior)

def parity(root):
    cache,m=load_compatible(root/'science/target-cache-compatibility.json');out=root/'science/retained-parity';out.mkdir()
    for name,relative in (('cp0','assessment-v2'),('pilot1','assessment')):
        old=load(previous.H14/f'science/target/{name}/{relative}/result.json');inv=load(verify(old['inventory']))
        bindings=load(previous.H14/'publication-v2/state-bindings.json')['states'];v=next(v for v in bindings if v['case']==name);r=load(verify(v['producer']))
        points=[v['physical']for v in inv['points']];rows,initial,stats=cache.view(points)
        # Independently enumerate all owners, preserving exact order/bits.
        index=previous.OwnerIndex(cache.report);reference_rows=[(list(p),o)for p in points for o in index(p)]
        assert rows.rows[:rows.point_count]==reference_rows
        result,_=previous.assess(rows,r['state'],initial,r['steps']);expected=old['states'][0]
        assert result['summary']==expected['summary']and result['all_regional_Q1_moments']==expected['all_regional_Q1_moments']
        atomic(out/(name+'.json'),dict(status='PASS',source_report=v['producer'],original_assessment=identity(previous.H14/f'science/target/{name}/{relative}/result.json'),inventory=old['inventory'],owned_rows=rows.point_count,
            owner_order_bits='EXACT',summary=result['summary'],all_regional_Q1_moments=result['all_regional_Q1_moments'],cache=stats))
    atomic(out/'result.json',dict(status='PASS',cache_compatibility=identity(root/'science/target-cache-compatibility.json'),scope='same original H14 states and original inventory; all owners/summary/interface/moment values exact'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('daemon','parity'),required=True);p.add_argument('--intent',type=Path);p.add_argument('--compat',type=Path);p.add_argument('--behavior',default='PASS');a=p.parse_args()
    if a.mode=='parity':parity(a.root)
    else:daemon(a.root,a.intent,a.compat,a.behavior)
