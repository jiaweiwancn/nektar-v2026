#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""All exact-optimization cost samples, including an unhelpful candidate."""
import argparse,json,statistics
from pathlib import Path
from analyze_wj_h25_profile import scopes
from run_wj_h4 import load
from wj_h15_observer import atomic,identity

def run(root):
    out=root/'science/benefit-v1';out.mkdir();records=[]
    for label,folder in (('W0-profile-four-v1','W0-profile-v1'),('lookup-on-four-v1','W1-lookup-v1'),('scalar-on-four-v1','W1-scalar-v1')):
        q=load(root/('science/'+folder+'/qualification.json'));assert q['status']=='PASS'
        timing=scopes(root,label,[17,18,19]);tele=[json.loads(s)for s in(root/('science/'+folder+'/monitor/telemetry.jsonl')).read_text().splitlines()]
        raw=[float(t['native_compute_seconds'])for t in tele if int(t['absolute_step'])<20]
        records.append(dict(operation=label,exact_parity=identity(root/('science/'+folder+'/qualification.json')),
            disjoint=timing,ordinary_native_compute_seconds=raw,mean_native_compute=statistics.mean(raw),minimum=min(raw),maximum=max(raw),
            operation_receipt=load(root/('operations/'+label+'/receipt.json'))))
    best=min(records,key=lambda r:r['mean_native_compute'])
    atomic(out/'result.json',dict(status='MEASURED',records=records,lowest_measured_ordinary=best['operation'],
        limitation='same advancing states and IEEE parity; short samples and cross-run variability, no cherry-picked rerun or statistically established production speedup',
        scalar_component=identity(root/'science/scalar-value-component-v1/result.json')))
    print(json.dumps([dict(operation=r['operation'],ordinary=r['ordinary_native_compute_seconds'],equation_fraction=r['disjoint']['equation_fraction'],
        nested=r['disjoint']['nested_explanatory_only'])for r in records],indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True,type=Path);run(p.parse_args().root)
