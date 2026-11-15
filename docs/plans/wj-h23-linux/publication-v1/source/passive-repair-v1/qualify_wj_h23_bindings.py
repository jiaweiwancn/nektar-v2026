#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only regression for the prelaunch path mismatch; no governor admission."""
import argparse,os,sys
from pathlib import Path

def qualify(root):
    # The repaired role helper is adjacent; immutable numerics remain in v3.
    base=root/'source-v3/solvers/IncNavierStokesSolver/Tests/RANS'
    sys.path.insert(1,str(base))
    from run_wj_h4 import THREADS,load
    os.environ.update(THREADS)
    from wj_h15_observer import identity,verify,atomic
    from wj_h17_observer import compatibility,load_compatible
    from wj_h23_continuation import active_execution_bindings
    from run_wj_h23 import reserve_native
    old=root/'science/target-cache-compatibility.json';old_pin=identity(old)
    try:active_execution_bindings(old_pin)
    except AssertionError as e:reason=str(e)
    else:raise AssertionError('the actual stale source path fixture was not rejected')
    out=root/'science/passive-binding-repair-v1';out.mkdir()
    path=out/'compatibility.json';historical=verify(load(old)['historical_manifest'])
    compatibility(root,historical,path)
    proof=active_execution_bindings(identity(path));cache,manifest=load_compatible(path)
    assert identity(old)==old_pin
    book=load(root/'accounting.json');assert book['native_halt']and book['native_reservations']==['M','strict_M']
    try:reserve_native(book,'C')
    except AssertionError:pass
    else:raise AssertionError('passive repair reopened native admission')
    atomic(out/'result.json',dict(status='PASS_PASSIVE_REPAIR_ONLY',original_failure=book['native_halt'],original_compatibility=old_pin,
        old_active_paths_rejected=reason,new_compatibility=identity(path),new_import_proof=proof,
        immutable_cache_load='PASS',cache=manifest['cache'],base_owned_rows=cache.base.point_count,
        old_records_unchanged=True,native_operations=0,native_requalification='NOT_RUN',campaign_governors_remain_closed=True))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();qualify(a.root)
