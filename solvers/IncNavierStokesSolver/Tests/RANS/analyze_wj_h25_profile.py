#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Exclusive native tree accounting; accumulated domain buckets stay nested."""
import argparse,collections,csv,json,statistics
from pathlib import Path
from run_wj_h4 import load
from wj_h15_observer import identity,atomic

def category(name,ancestors):
    names=[name]+ancestors
    if any(n.startswith(('audit.','physical.samples'))for n in names):return 'full_scientific_audit'
    if any(n in ('monitor.observer_wait','monitor.independent_wait')for n in names):return 'observer_enclosing_wait'
    if any(n.startswith('monitor.')for n in names):return 'recovery_publication_IO'
    if any(n.startswith(('domain.','k.bound','verify.','guard.','residual'))or n in ('flow.correction_validate',)for n in names):return 'guards_bounds_residual_verification'
    if name.startswith('data.'):return 'data_allocation_release'
    if name.startswith(('flow.','momentum.','pressure.','k.','epsilon.','mapped.','mass.','dense.','matching.','constitutive.')):return 'equation_evaluation_assembly_solution_update'
    if name in ('native.advance','native.explicit'):return 'equation_evaluation_assembly_solution_update'
    if name=='pimple.sweep':return 'attempt_evidence_and_orchestration'
    if name=='native.step':return 'loop_telemetry_orchestration'
    return 'UNATTRIBUTED'

def scopes(root,label,steps):
    path=root/('operations/'+label+'/timing.jsonl')
    rows=[json.loads(s)for s in path.read_text().splitlines()]
    tree={s['id']:s for s in rows if 'clock_scope'not in s}
    sums=collections.defaultdict(float);by_name=collections.defaultdict(float);counts=collections.Counter();nested=collections.defaultdict(float)
    roots=[s for s in rows if s['name']=='native.step'and s['step']in steps]
    rootids={s['id']for s in roots};total=sum(s['wall_seconds']for s in roots)
    assert len(roots)==len(steps),(label,len(roots),steps)
    for s in rows:
        if s['step']not in steps:continue
        if 'clock_scope'in s:
            nested[s['name']]+=s['wall_seconds'];continue
        a=[];p=s['parent_id'];ids={s['id']}
        while p:
            assert p in tree and p not in ids,'broken timing tree'
            ids.add(p);a.append(tree[p]['name']);p=tree[p]['parent_id']
        if not ids&rootids:continue
        c=category(s['name'],a);sums[c]+=s['exclusive_wall_seconds'];by_name[s['name']]+=s['exclusive_wall_seconds'];counts[s['name']]+=1
    assert abs(sum(sums.values())-total)<1e-7*max(1,total),'exclusive tree does not close'
    n=len(steps)
    return dict(timing=identity(path),steps=steps,wall_seconds=total,per_step=total/n,
        attributed_fraction=1-sums['UNATTRIBUTED']/total,
        equation_fraction=sums['equation_evaluation_assembly_solution_update']/total,
        disjoint={k:v/n for k,v in sums.items()},
        exclusive_by_name={k:dict(seconds_per_step=v/n,calls=counts[k])for k,v in sorted(by_name.items(),key=lambda kv:-kv[1])},
        nested_explanatory_only={k:v/n for k,v in nested.items()},
        samples=[s['wall_seconds']for s in roots])

def analyze(root):
    out=root/'science/phase-a-v2';out.mkdir()
    a=scopes(root,'W0-profile-four-v1',[17,18,19])
    b=scopes(root,'lookup-off-two-v1',[17])
    # The new release scope identifies the old callback remainder. Keep both
    # measured trees: no retroactive alteration of the first profiling receipt.
    assert b['attributed_fraction']>=.95,'meaningful95% wall attribution not established'
    stock=load(root/'science/S0-explicit-opt/result.json');mass=(root/'operations/common-mass-v1/command.log').read_text().splitlines()[0]
    result=dict(status='PASS_PHASE_A',baseline_profile=a,resolved_profile=b,
        stock=stock,common_operator=json.loads(mass),
        exact_advancing_parity=identity(root/'science/W0-profile-v1/qualification.json'),
        second_exact_advancing_parity=identity(root/'science/W0-lookup-off-v1/qualification.json'),
        ranked_opportunities=[
            dict(rank=1,cost='moving-ray FE containment and row construction',evidence=a['nested_explanatory_only'],
                 next='exact native interval/geometry lookup; preserve all scalar contractions and ray/profile checks'),
            dict(rank=2,cost='scalar inventory materialization and destruction',evidence=b['disjoint'].get('data_allocation_release',0),next='dependency-correct compact inventory or exact positivity certificate with fallback if first candidate remains insufficient'),
            dict(rank=3,cost='flow element assembly and pressure mass applications',evidence=a['exclusive_by_name'],next='element-local scatter/batched consistent mass, only if measured benefit warrants'),
            dict(rank=4,cost='step count and fixed audit burden',next='separate350-cell representation/space screen, quadrature and perturbed common-time dt qualification')],
        taxonomy_note='All exclusive intervals close to native.step; accumulator domain/principal buckets are nested explanations, never additional wall. Guards/audits are not equation work. RANS-vs-stock differences include different equations/pressure/BCs.',
        status_limit='Checkpoint A is a measured cost/attribution checkpoint, not production readiness')
    atomic(out/'result.json',result)
    with(out/'disjoint-cost.csv').open('x',newline='')as f:
        w=csv.writer(f);w.writerow(['category','old_profile_seconds_per_step','resolved_profile_seconds_per_step'])
        for k in sorted(set(a['disjoint'])|set(b['disjoint'])):w.writerow([k,a['disjoint'].get(k,0),b['disjoint'].get(k,0)])
    print(json.dumps(dict(status=result['status'],attributed=b['attributed_fraction'],equation_fraction=b['equation_fraction'],
        per_step=b['per_step'],disjoint=b['disjoint']),indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True,type=Path);analyze(p.parse_args().root)
