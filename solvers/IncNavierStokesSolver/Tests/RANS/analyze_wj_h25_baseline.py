#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Disjoint native baseline scopes; mixed/unresolved work is not called equation work."""
import argparse,collections,json,statistics
from pathlib import Path
from run_wj_h4 import load
from wj_h15_observer import identity,atomic

def analyze(root):
    out=root/'science/phase-a-v1';out.mkdir()
    operation=root/'operations/W0-700-eight-v1';scopes=[json.loads(s)for s in(operation/'timing.jsonl').read_text().splitlines()]
    selected=[s for s in scopes if 17<=s['step']<=23]
    sums=collections.defaultdict(lambda:[0.,0.,0])
    for s in selected:
        v=sums[s['name']];v[0]+=s['wall_seconds'];v[1]+=s['exclusive_wall_seconds'];v[2]+=1
    table=[dict(name=k,wall_per_ordinary_step=v[0]/7,exclusive_per_ordinary_step=v[1]/7,calls=v[2])for k,v in sorted(sums.items(),key=lambda kv:-kv[1][1])]
    tele=[json.loads(s)for s in(root/'science/W0/monitor/telemetry.jsonl').read_text().splitlines()]
    ordinary=[float(r['native_compute_seconds'])for r in tele if int(r['absolute_step'])<24]
    stock=load(root/'science/S0-explicit-opt/result.json');warm=[r['seconds']for r in stock['samples'][1:]]
    result=dict(status='MEASURED_BASELINES_PASS',W0=dict(ordinary_seconds=ordinary,mean=statistics.mean(ordinary),minimum=min(ordinary),maximum=max(ordinary),
        two_warm_segments=[ordinary[:3],ordinary[3:]],native_completion_count=len(tele),timing=identity(operation/'timing.jsonl'),disjoint_scopes=table),
        S0=dict(samples=stock['samples'],two_warm_segments=[warm[:3],warm[3:]],warm_mean=statistics.mean(warm),full_process_seconds=stock['seconds'],
            result=identity(root/'science/S0-explicit-opt/result.json'),SVV=False),
        comparison=dict(mesh='same700 physical quads',dt=.0025,interval=[.04,.06],workers=1,velocity='Q3',pressure='WJ regionalQ1 versus stock global continuousQ2',
            meaning='elapsed difference combines physics, pressure, boundary and implementation work; not identified as pure RANS overhead'),
        mixed_purpose_scopes=['momentum.solve_validate','pressure.solve_validate','k.solve_validate','epsilon.solve_validate','native.implicit'],
        attribution_status='raw exclusive scopes retained; refine significant mixed/exclusive remainder before final95% claim',
        source=identity(root/'source-manifest-v2.json'))
    atomic(out/'result.json',result)
    print(json.dumps(dict(W0_mean=result['W0']['mean'],S0_mean=result['S0']['warm_mean'],top_disjoint=table[:16]),indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();analyze(a.root)
