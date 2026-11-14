#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Signed terminal weak data and disjoint timing, supplementing sealed v1."""
import argparse
import gzip
import json
import os
from pathlib import Path
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import atomic,identity,verify
from finish_wj_h15 import csv_rows

def supplement(root):
    verify(load(root/'publication-v1/seal.json')['manifest']);out=root/'publication-supplement-v1';out.mkdir();folder=root/'science/target'
    v=load(folder/'view/result.json');atomic(out/'terminal-signed-weak-and-interface.json',dict(status=v['status'],source=identity(folder/'view/result.json'),weak=v['weak'],interface=v['interface'],mass_dual=v['physical_mass_dual_norms'],scope='native projected1,x/5,interface_hat contractions; unchanged q12 operator; no new quadrature claim'))
    summary=load(folder/'case.evolution.json');calls=[json.loads(v)for v in(folder/'monitor/factors/calls.jsonl').read_text().splitlines()]
    assert[v['call']for v in calls]==list(range(1,summary['factor_catalog']['calls']+1))
    assert all(v['resident_definitions']<=16 and v['resident_key_bytes']<=32*1024**2 for v in calls)
    atomic(out/'retention-and-factor-stream.json',dict(summary=summary,complete_calls=calls,pattern_definition_files=len(list((folder/'monitor/factors').glob('factor-*.json'))),source_calls=identity(folder/'monitor/factors/calls.jsonl'),scope='real target records; full streamed definitions retained; required witness union has a separate four-stage growth policy'))
    pid=load(folder/'native-launch.json')['pid'];timing=[]
    for p in sorted((root/'operations').glob('*/timing.jsonl')):
        for line in p.read_text().splitlines():
            row=json.loads(line)
            if row['pid']==pid:timing.append(dict(operation=p.parent.name,**row))
    assert timing;atomic(out/'target-native-timing.json',timing)
    csv_rows(out/'target-disjoint-and-inclusive-timing.csv',[{k:v.get(k)for k in('pid','step','time','name','calls','wall_seconds','exclusive_wall_seconds','process_cpu_seconds','clock_scope')}for v in timing])
    steps=[v for v in timing if v['name']=='native.step'and v['step']in(2,3,4)]
    assert sorted(v['step']for v in steps)==[2,3,4]
    inclusive=sum(v['wall_seconds']for v in steps)/3
    analysis=load(root/'science/analysis-v1/summary.json')
    atomic(out/'conditional-cost-and-limits.json',dict(**analysis['costs'],observed_inclusive_warm_step_mean_seconds=inclusive,observed_loop_conditional_days=inclusive*202010/86400,clocks=analysis['clocks'],scope='ready.native_compute ends before publication/wait; native.step includes them; nested timers never summed; accumulated audit CPU0 is unmeasured; duration/application/temporal admission remains CLOSED'))
    budgets=[];fetch=[]
    for step in range(5):
        p=root/f'science/analysis-v1/step-{step}.json';v=load(p);b=v['summary']['interface_budget']
        budgets.append(dict(step=step,time=step*.0025,source=identity(p),inner_integrals=b['inner_integrals'],wall_shear_range=b['wall_shear_range'],role=b['note']))
        p=root/f'science/analysis-v1/step-{step}-fetch.json';r=load(p);fields={}
        for field in('U','G','k','pi','K','rc','epsilon'):
            values=[item[field]['range']for item in r if field in item]
            def flat(v):return[float(v)]if not isinstance(v,list)else[x for row in v for x in flat(row)]
            fields[field]=max(abs(x)for value in values for x in flat(value))
        fetch.append(dict(step=step,time=step*.0025,source=identity(p),maximum_component_fetch_range=fields))
    atomic(out/'actual-and-counterfactual-budget-history.json',budgets);atomic(out/'fetch-range-history.json',fetch)
    compressed=[]
    for name in('step-4-fetch.json','step-4-exterior.json'):
        p=root/'science/analysis-v1'/name;dest=out/(name+'.gz')
        with p.open('rb')as f,gzip.open(dest,'xb')as g:
            import shutil
            shutil.copyfileobj(f,g)
        compressed.append(dict(raw=identity(p),compressed=identity(dest)))
    atomic(out/'compressed-observation-bindings.json',compressed)
    atomic(out/'manifest.json',dict(original_publication=identity(root/'publication-v1/manifest.json'),files=[identity(p)for p in sorted(out.iterdir())if p.is_file()],scope='no native/readback work; original publication-v1 immutable'))
    atomic(out/'seal.json',dict(manifest=identity(out/'manifest.json')))
    print('H16 signed budget/timing/retention supplement PASS',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();supplement(a.root)
