#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H24 binding repair and unchanged H23/H20 continuation mathematics."""
import argparse,os,sys
from pathlib import Path
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,verify,atomic
from wj_h17_observer import load_compatible
from wj_h23_continuation import active_execution_bindings,ORDER,DT,START,END,ADDITIONAL,STRIDE
import study_wj_h23 as inherited
from run_wj_h24 import ROOT,frozen

def check_bindings(root):
    return active_execution_bindings(identity(root/'science/target-cache-compatibility.json'))

def setup(root):
    # The capsule is already at its final immutable execution path. All imports,
    # including the future observer's, resolve here before compatibility is made.
    inherited.setup(root)
    old=load(root/'inputs/historical-indices.json')['H23']
    for pin in old.values():verify(pin)
    book=load(verify(old['accounting']));assert book['closed']and book['actual_native_advances']==0
    assert book['total_charged_seconds']==550.5544812529697
    atomic(root/'inputs/h23-closed-link.json',dict(index=old,accounting=old['accounting'],historical_charges=book['total_charged_seconds'],
        historical_failure=book['native_halt'],historical_reservations=book['native_reservations'],H24_budget='new independent6000/4200/900, no transfer'))
    from run_wj_h7 import update_journal
    update_journal(root/'accounting.json',lambda b:b.update(historical_H23_accounting=old['accounting'],historical_H23_charged_seconds=book['total_charged_seconds']))

def bindings(root):
    current=check_bindings(root)
    old=Path('/home/jay/build/wj-h23-flat-continuation-20261115/science/target-cache-compatibility.json');old_pin=identity(old)
    try:active_execution_bindings(old_pin)
    except AssertionError as error:stale=str(error)
    else:raise AssertionError('actual stale H23 path fixture did not fail')
    cache,manifest=load_compatible(root/'science/target-cache-compatibility.json')
    assert identity(old)==old_pin
    atomic(root/'science/actual-launch-binding-proof.json',dict(status='PASS',current=current,actual_namespace=sys.path,
        observer=str(Path(__file__).with_name('wj_h23_observer.py')),historical_stale_rejected=stale,stale_record_unchanged=old_pin,
        immutable_cache=manifest['cache'],reference_identity=cache.reference_identity,variable_limit=cache.variable_limit,
        frozen_source_path=str(Path(__file__).resolve().parent),native_execution='NOT_YET_RUN'))

def freeze(root):
    proof=check_bindings(root)
    assert load(root/'science/preflight/result.json')['status']=='PASS'
    assert load(root/'science/actual-launch-binding-proof.json')['status']=='PASS'
    assert load(root/'operations/focused-controls-v1/receipt.json')['status']=='PASS'
    files=list(load(root/'source-manifest.json')['files'])
    files.extend(identity(root/'inputs'/n)for n in ('selected.json','historical-indices.json','historical-immutability.json','h23-closed-link.json'))
    files.extend(identity(root/'science'/n)for n in ('runtime-build-v1.json','runtime-dependencies-v1.json','compiled-source-objects.json',
        'target-cache-compatibility.json','setup.json','preflight/result.json','actual-launch-binding-proof.json'))
    files.extend(identity(root/f'science/{z}/{n}')for z in DT for n in ('case.xml','execution-diff.json'))
    atomic(root/'execution-config.json',dict(schema='WJH24FrozenExecutionV1',files=files,compatibility=identity(root/'science/target-cache-compatibility.json'),
        active_import_bindings=proof,source_manifest=identity(root/'source-manifest.json'),native_capsule='057e6b7465660cac5ec8617abf5c876fa8bff71f',
        ORDER=ORDER,DT=DT,START=START,END=END,additional=ADDITIONAL,scientific_stride=STRIDE,recovery_stride=1,
        expected_accepted_advances=28,expected_governors=3,terminal_strict=3,normal_captures='Off',
        internal_schema='inherited qualified H23/H20 numerical schemas; H24 root, source, configuration, nonce and accounting identify the new execution',
        adaptive_authorization='versioned routine repair and necessary failed-part successor within unchanged science and6000/4200/900; preserve all failures'))

def admission(root,level):
    frozen(root);check_bindings(root)
    inherited.admission(root,level)

def target(root,level):
    frozen(root);check_bindings(root)
    inherited.target(root,level,prelaunch=lambda:check_bindings(root))

def strict(root,level):
    frozen(root);inherited.strict(root,level)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    p.add_argument('--mode',choices=('setup','bindings','preflight','freeze','admission',*ORDER,'strict_M','strict_C','strict_F'),required=True)
    p.add_argument('--level',choices=ORDER);a=p.parse_args();assert a.root.name==ROOT
    if a.mode in ORDER:target(a.root,a.mode)
    elif a.mode.startswith('strict_'):strict(a.root,a.mode.split('_')[1])
    elif a.mode=='admission':admission(a.root,a.level)
    elif a.mode=='preflight':inherited.preflight(a.root)
    else:globals()[a.mode](a.root)
