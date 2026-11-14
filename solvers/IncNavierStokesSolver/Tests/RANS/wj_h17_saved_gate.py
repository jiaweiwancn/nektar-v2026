#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""No FE evaluation: compare saved baseline payloads and project remaining cost."""
import argparse
from pathlib import Path
from run_wj_h4 import load
from wj_h17_observer import native_load,equal
from wj_h15_observer import atomic,verify

def saved_baseline(root):
    v=native_load(root/'science/audits/Off/native.json');checks=[]
    assert v['full_audit_calls']==6
    for row in v['results']:
        a=row['complete']['audit']
        assert a['admissible'][0]==a['wj_hard_admissible'][0]==1
        assert all(x[0]==1 for k,x in a.items()if k.startswith('wj_condition_'))
        if row['step']:
            old=load(root/('inputs/target/step-'+str(row['step'])+'.json'))
            baseline=native_load(verify(old['report']))
            equal(a,baseline['native_hard_audit'],'retained H16 complete audit')
            checks.append(dict(step=row['step'],status='PASS',producer=old['report']))
    b=load(root/'accounting.json')
    q=sum(x.get('seconds',x['reserved_seconds'])for x in b['entries'].values()if x['phase']=='qualification')
    forecast=1.25*(b['entries']['complete-Off-v1']['seconds']+22.701710678+17.450035246+7.795539722+20)
    assert q+forecast<1800
    atomic(root/'science/after-baseline-admission.json',dict(status='PASS',retained_H16_complete_numerical_audits=checks,baseline_seconds=b['entries']['complete-Off-v1']['seconds'],qualification_spent=q,mandatory_remaining_with25percent_headroom=forecast,qualification_ceiling=1800,closure_reserved=300,target_strict_reserved=120))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();saved_baseline(a.root)
