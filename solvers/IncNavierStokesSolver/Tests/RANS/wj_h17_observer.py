#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H17 provenance/parity gates around the unchanged independent evaluator."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import struct
import time
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
import wj_h15_observer as previous
from wj_h15_observer import identity,verify,atomic,vector_hash
from wj_h16_observer import active_source_bindings,packet_fields
from wj3_measure import checkpoint_payload

OLD=Path('/home/jay/build/wj-h16-repaired-monitor-20261114')
ORIGINAL_LOADER=previous.load_cache
ORIGINAL_VERIFY=previous.verify
ORIGINAL_ATOMIC=previous.atomic

def native_load(path):
    # Native max_digits10 emits -0 as an integer token; preserve its IEEE sign.
    return json.loads(Path(path).read_text(),parse_int=lambda s:-0.0 if s=='-0'else int(s))

def exact(value):
    """Canonical typed tree, retaining binary64 bits and ordered multiplicity."""
    if isinstance(value,bool)or value is None:return ('atom',value)
    if isinstance(value,(int,float)):return ('number',struct.pack('<d',float(value)).hex())
    if isinstance(value,str):return ('text',value)
    if isinstance(value,list):return ('list',tuple(exact(v)for v in value))
    if isinstance(value,dict):return ('map',tuple((k,exact(v))for k,v in sorted(value.items())))
    raise TypeError(type(value))

def equal(a,b,label):
    assert exact(a)==exact(b),'bitwise numerical/ordered parity: '+label

def compatibility(root,manifest,destination):
    m=load(manifest);assert m['status']=='PASS'
    sources=active_source_bindings(m)
    for v in m['source_files']:verify(v)
    for k in ('profile','native_runtime_manifest'):verify(m['reference_identity'][k])
    cache=verify(m['cache'])
    if cache.stat().st_size==211582940:assert m['cache']['sha256']=='800ac28c73ecb6498e6504b4e14e774db2cff893f2c0aaf418cb7036d50d00b0'
    runtime=load(root/'science/runtime-build-v1.json')
    for v in runtime['executables']:verify(v['ELF']);verify(v['wrapper'])
    for v in load(root/'science/runtime-dependencies-v1.json')['libraries']:verify(v['actual'])
    c=dict(status='PASS',schema='WJH17CacheCompatibilityV1',historical_manifest=identity(manifest),historical_cache=m['cache'],historical_reference_identity=m['reference_identity'],active_sources=sources,
        active_adapter=identity(__file__),new_native_identity=identity(root/'science/runtime-build-v1.json'),runtime_dependencies=identity(root/'science/runtime-dependencies-v1.json'),
        permitted_difference='native exact structural-row lifetime/access only; independent geometry/profile/reference/observer operations unchanged',historical_pins_modified=False,actual_imports_verified=True)
    atomic(destination,c);return c

def load_compatible(path):
    c=load(path);assert c['status']=='PASS'and c['actual_imports_verified']
    for k in ('active_adapter','new_native_identity','runtime_dependencies'):verify(c[k])
    m=load(verify(c['historical_manifest']));assert active_source_bindings(m)==c['active_sources']
    for v in load(c['runtime_dependencies']['path'])['libraries']:verify(v['actual'])
    return ORIGINAL_LOADER(verify(c['historical_manifest']))

def phase_parity(current,expected):
    """Only completed packets, with genuine callback and serialized history."""
    packet_fields(current);packet_fields(expected)
    for k in ('state','fixed_reference','fixed_defect','native_history','native_history_times','time','steps','dt','signature','last_stage_history','last_stage_alpha','last_stage_time','native_hard_audit'):
        equal(current[k],expected[k],k)
    for k in ('stage_audits','step_audits','explicit_audits'):
        equal(current[k],expected[k],k)
    # Kernel result/witness ordering is part of native audit parity, not timing.
    for k in ('state_sha256','stage_history_sha256','time','alpha_mass','principal_witnesses'):
        equal(current['wj_audit_packets'][0][k],expected['wj_audit_packets'][0][k],k)

def independent_parity(result,inventory,expected,old_inventory):
    for k in ('rule','points','owned','physical_rows','quadrature_rows'):equal(inventory[k],old_inventory[k],'live inventory '+k)
    for k in ('summary','all_regional_Q1_moments'):equal(result[k],expected[k],'independent '+k)
    assert result['status']==expected['status']=='PASS'

def daemon(root,intent_path,compat,behavior='PASS'):
    # Reuse H15 numerical loops/cache exactly; hooks act before result and ACK.
    def loader(unused):
        began=time.monotonic();cache,m=load_compatible(compat)
        directory=Path(load(intent_path)['directory'])
        ORIGINAL_ATOMIC(directory/'h17-cache-load.json',dict(seconds=time.monotonic()-began,compatibility=identity(compat),immutable_row_rebuilds=0,base_rows=len(cache.base.rows),base_owned_rows=cache.base.point_count))
        return cache,m
    def pinned(v):
        p=ORIGINAL_VERIFY(v);intent=load(intent_path);directory=Path(intent['directory'])
        if p.parent==directory and p.name.endswith('.report.json'):
            r=native_load(p);packet_fields(r);n=str(r['steps'])
            if intent.get('replay_expected'):
                old=load(verify(intent['replay_expected'][n]));orig=native_load(verify(old['report']))
                phase_parity(r,orig)
                # Read actual serialized coefficient/history arrays too.
                ready=load(directory/f'step-{n}.ready.json');cp=checkpoint_payload(verify(ready['checkpoint']).parent,include_reference=True)
                baseline=checkpoint_payload(verify(old['checkpoint']).parent,include_reference=True)
                from study_wj_h16 import same_cp
                same_cp(cp,baseline)
                equal(native_load(verify(ready['callback']))['callback_history'],native_load(verify(old['callback']))['callback_history'],'genuine callback sidecar')
                dest=directory/f'step-{n}.native-replay-parity.json'
                if not dest.exists():ORIGINAL_ATOMIC(dest,dict(status='PASS',current=v,baseline=old['report'],checkpoint=old['checkpoint'],callback=old['callback'],scope='exact state/reference/defect/native and callback history/time/native numerical audit before independent result/ACK'))
        return p
    def durable(p,v,immutable=True):
        p=Path(p);intent=load(intent_path)
        if p.name.endswith('.ready.result.json')and intent.get('replay_expected'):
            step=str(int(v['binding']['absolute_step']));old=load(verify(intent['replay_expected'][step]));expected=load(verify(old['result']))
            inv=load(p.parent/f'step-{step}.ready.inventory.json');oldinv=load(verify(old['inventory']))
            independent_parity(v,inv,expected,oldinv)
            ORIGINAL_ATOMIC(p.parent/f'step-{step}.independent-replay-parity.json',dict(status='PASS',baseline_result=old['result'],baseline_inventory=old['inventory'],owned_rows=len(inv['owned']),points=len(inv['points']),scope='same live inventory/order and all summary/moment bits before final ACK'))
        return ORIGINAL_ATOMIC(p,v,immutable)
    previous.load_cache=loader;previous.verify=pinned;previous.atomic=durable
    previous.daemon(root,intent_path,compat,behavior)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--intent',type=Path,required=True);p.add_argument('--compat',type=Path,required=True);p.add_argument('--behavior',default='PASS');a=p.parse_args();daemon(a.root,a.intent,a.compat,a.behavior)
