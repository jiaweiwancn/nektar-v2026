#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Select only accuracy-qualified implementation gains; retain every regression."""
import argparse,statistics
from pathlib import Path
from run_wj_h4 import load
from wj_h15_observer import identity,atomic
from analyze_wj_h25_profile import scopes


def run(root):
    base=root/'resumption-r1';records={}
    for label,steps in (('scatter-off-two-r1-v1',[17]),('scatter-on-four-r1-v1',[17,18,19]),
                        ('memo-off-two-r1-v1',[17]),('memo-on-four-r1-v1',[17,18,19])):
        s=scopes(base,label,steps);s['stdev']=statistics.stdev(s['samples'])if len(steps)>1 else None
        s['minimum']=min(s['samples']);s['maximum']=max(s['samples']);records[label]=s
    a=records['scatter-on-four-r1-v1'];b=records['memo-off-two-r1-v1'];c=records['memo-on-four-r1-v1']
    for label in ('W1-scatter-v1','W1-memo-v1'):
        assert load(base/f'science/{label}/qualification.json')['status']=='PASS'
    assert load(base/'science/memo-component-v1/native.json')['status']=='PASS'
    retained=c['per_step']<.98*min(a['per_step'],b['per_step']) and c['attributed_fraction']>=.95
    chosen='memo'if retained else'scatter';q=base/f'science/W1-{chosen}-v1/qualification.json'
    runtime=base/f'science/runtime-build-{chosen}-v1.json'
    selected=dict(status='PASS_EXACT_IMPLEMENTATION_ONLY',runtime=identity(runtime),qualification=identity(q),
        scalar_policy='MemoExactV1'if retained else'Off',assembly_policy='ElementExactV1',lookup_policy='TensorExactV1',
        record='memo-on-four-r1-v1'if retained else'scatter-on-four-r1-v1',
        limitation='no transient h/q/dt, full-duration or production-cadence admission')
    atomic(base/'science/implementation-selected-v1.json',selected)
    atomic(base/'science/implementation-benefit-v1.json',dict(status='PASS_IMPLEMENTATION_COMPARISON',
        records=records,selected=selected,retention_rule='On mean < .98*both new Off and scatter; exact parity; attribution>=.95',
        memo_decision='RETAIN'if retained else'CORRECT_BUT_NOT_RETAINED',ordinary_definition='actual native.step including current recovery/publication, excluding scheduled audit steps',
        scope='short matched evolving comparisons; native root timing only; no scientific horizon added by replays, no production confidence interval'))
    print(chosen,{k:v['per_step']for k,v in records.items()},flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);run(p.parse_args().root)
