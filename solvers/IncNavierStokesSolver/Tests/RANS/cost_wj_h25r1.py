#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Complete conditional flat-control forecast from disjoint measured scopes."""
import argparse,collections,math,statistics
from pathlib import Path
from run_wj_h4 import load
from wj_h15_observer import identity,atomic
from wj_h18_observer import due
from analyze_wj_h25_profile import scopes
from run_wj_h25 import usage

T=505.02466521278404
WINDOW=[447.8483832973273,T]
EQ='equation_evaluation_assembly_solution_update'


def stats(samples):
    return dict(count=len(samples),mean=statistics.mean(samples),minimum=min(samples),maximum=max(samples),
        stdev=statistics.stdev(samples)if len(samples)>1 else None,samples=samples)


def block(base,operation,folder):
    p=base/'operations'/operation;r=load(p/'receipt.json');assert r['status']=='PASS'
    d=base/'science'/folder;q=load(d/'qualification.json');intent=load(d/'intent.json')
    allsteps=list(range(q['start']+1,q['terminal']+1));warm=allsteps[1:]
    ordinary=[n for n in warm if not due(n,intent)];scientific=[n for n in warm if due(n,intent)]
    timing=[__import__('json').loads(s)for s in(p/'timing.jsonl').read_text().splitlines()]
    process=next(s for s in timing if s['name']=='process'and 'clock_scope'not in s)
    allprofile=scopes(base,operation,allsteps)
    fixed=r['seconds']-allprofile['wall_seconds'];assert fixed>0
    started=load(d/'monitor/observer-started.json')
    roots={s['name']:s['wall_seconds']for s in timing if s.get('parent_id')==process['id']and 'clock_scope'not in s}
    return dict(operation=identity(p/'receipt.json'),qualification=identity(d/'qualification.json'),
        ordinary=scopes(base,operation,ordinary)if ordinary else None,
        scientific=scopes(base,operation,scientific)if scientific else None,
        whole_advancing_profile=allprofile,first_step_excluded_from_warm=allsteps[0],
        fixed_seconds=fixed,fixed_disjoint=dict(native_startup=sum(roots.get(n,0)for n in ('startup.object','startup.initialize')),
            terminal_fields_and_report=roots.get('output.fields_and_report',0),
            other_native_lifecycle=process['wall_seconds']-allprofile['wall_seconds']-sum(roots.get(n,0)for n in ('startup.object','startup.initialize','output.fields_and_report')),
            controller_observer_setup_and_validation=r['seconds']-process['wall_seconds']),
        observer_setup_seconds=started['setup_seconds'],native_resource_summary=q['native_resource_summary'])


def run(root,provisional=False):
    base=root/'resumption-r1';chosen=load(base/'science/implementation-selected-v1.json')
    assert load(base/'science/C-wide-v1/temporal-qualification.json')['status']=='PASS_SHORT_UNPERTURBED_WINDOW'
    definitions=[('C-wide-r1-v1','C-wide-v1')]
    if not provisional:definitions += [('C-sparse-restart-r1-v1','C-sparse-restart-v1'),('C-dense-restart-r1-v1','C-dense-restart-v1')]
    blocks={op:block(base,op,folder)for op,folder in definitions}
    ordinary=max((v['ordinary']for v in blocks.values()if v['ordinary']),key=lambda v:v['per_step'])
    scientific=max((v['scientific']for v in blocks.values()if v['scientific']),key=lambda v:v['per_step'])
    fixed=max(v['fixed_seconds']for v in blocks.values())
    for v in blocks.values():assert abs(sum(v['fixed_disjoint'].values())-v['fixed_seconds'])<1e-7
    fresh=load(base/'operations/fresh-C-r1-v1/receipt.json')['seconds']if not provisional else 600.
    strict=load(base/'operations/C-terminal-read-r1-v1/receipt.json')['seconds']if not provisional else 180.
    if not provisional:
        assert load(base/'science/fresh-C-v1/qualification.json')['status']=='PASS_GENUINE_FRESH_INITIALIZATION'
        assert load(base/'science/C-sparse-restart-v1/strict/qualification.json')['status']=='PASS'
        assert load(base/'science/C-sparse-restart-v1/rotation.json')['status']=='PASS'
    dt=.005;end=math.ceil(T/dt);forecasts={}
    for stride in (4,50):
        scenarios={}
        for name,start in (('fresh_two_stage',0),('continuation_C8',8)):
            count=end-start;audits=end//stride-start//stride+int(end%stride!=0);normal=count-audits
            disjoint={k:normal*ordinary['disjoint'].get(k,0)+audits*scientific['disjoint'].get(k,0)
                for k in set(ordinary['disjoint'])|set(scientific['disjoint'])}
            # Two advancing processes: initial restore and one planned restart.
            # Each fixed block excludes *all* native.step intervals. Each block
            # has its own startup, terminal export and controller lifecycle.
            # Fresh construction is the separately executed genuine NumSteps=0
            # projection/history/audit/archive stage, followed by ordinary import.
            disjoint['two_process_startup_terminal_controller']=2*fixed
            disjoint['strict_terminal_read']=strict or 0
            if start==0:disjoint['genuine_fresh_construction_initial_archive']=fresh or 0
            total=sum(disjoint.values());eq=disjoint.get(EQ,0);unknown=disjoint.get('UNATTRIBUTED',0)
            scenarios[name]=dict(status='CONDITIONAL_UNPERTURBED_COST_ONLY',start_step=start,start_time=start*dt,
                accepted_advances=count,scientific_steps=audits,ordinary_steps=normal,dt=dt,
                completed_end=end*dt,required_horizon=T,terminal_window=WINDOW,
                disjoint_seconds=disjoint,seconds=total,hours=total/3600,with_25_percent_seconds=1.25*total,
                meets_24_hour_gate=1.25*total<=86400,equation_fraction=eq/total,
                equation_fraction_note='conservative numerator: fixed-stage equation work is not added',
                attributed_fraction=1-unknown/total,
                fixed_costs_measured=not provisional,
                warning='fresh construction 600s and strict read 180s are prospective operation caps, not measured results'if provisional else None)
        forecasts[str(stride)]=dict(scientific_stride=stride,scientific_interval=stride*dt,recovery_stride=1,
            scope='short unperturbed cadence exercised; finite-perturbation production admission open'if stride==4 else'UNEXECUTED .25 cadence sensitivity only; no accuracy or operational admission',scenarios=scenarios)
    storage={};scientific_sizes=[]
    for _,folder in definitions:
        d=base/'science'/folder;intent=load(d/'intent.json');q=load(d/'qualification.json')
        records={}
        for n in range(q['start']+1,q['terminal']+1):
            if due(n,intent):
                records[n]=sum(p.stat().st_size for p in(d/'monitor').glob(f'step-{n}.*')if p.is_file())
                scientific_sizes.append(records[n])
        storage[folder]=dict(scientific_record_bytes=records,retained_bytes=sum(p.stat().st_size for p in d.rglob('*')if p.is_file()),
            logs_bytes=sum(p.stat().st_size for p in d.rglob('*')if p.is_file()and '.log'in p.name and not p.name.endswith('.json')))
    per_record=max(scientific_sizes)
    for value in forecasts.values():
        count=value['scenarios']['fresh_two_stage']['scientific_steps']
        value['retained_scientific_archive_bytes']=count*per_record
        value['fits_8GiB_science_cap']=count*per_record<=8*1024**3
    account=load(base/'accounting.json');original=load(root/'accounting.json')
    result=dict(schema='WJH25R1CompleteCostV1',status='TARGET_NOT_MET',provisional=provisional,
        selected_implementation=chosen,selected_space=dict(cells=700,velocity_degree=3,pressure_degree=1,governing_quadrature=[12,16],dt=.005),
        accuracy_scope='cheapest actually compared short unperturbed C/M/F flat recipe; prescribed finite perturbation never admitted, production accuracy remains open',
        blocks=blocks,ordinary=ordinary,scientific=scientific,
        pooled_ordinary_variability=stats([x for v in blocks.values()if v['ordinary']for x in v['ordinary']['samples']]),
        chosen_cost_rule='largest measured block mean separately for ordinary, scheduled-scientific and fixed costs; first process advance excluded from warm samples; no fastest-step forecast',
        forecasts=forecasts,storage=storage,maximum_measured_scientific_record_bytes=per_record,
        science_rss_cap=8*1024**3,build_rss_cap=16*1024**3,
        measured_peak_science_rss_bytes=max(v['peak_tree_rss_bytes']for v in account['entries'].values()if v['phase']!='build'and 'peak_tree_rss_bytes'in v),
        campaign_snapshot=dict(original=usage(original),continuation=usage(account),
            note='snapshot may include this operation reservation; final closed-account artifact is authoritative'),
        limitations=['Finite seed principal quality failed before evolution; cadence and h/q/dt transient qualifications remain open.',
            'Measured variability is a short control-window observation, not a production confidence interval or long-horizon growth bound.',
            'Scientific archives extrapolated from complete retained pilot records exceed the existing 8 GiB development science cap; no production storage allocation is established.',
            'Fixed stages conservatively include executed source checks, validation and controller work; native.step, audits and waits are never added twice.'])
    name='provisional-cost-r2-v1.json'if provisional else'complete-cost-v1.json'
    assert not(base/'science'/name).exists();atomic(base/'science'/name,result)
    print(name,ordinary['per_step'],scientific['per_step'],forecasts['4']['scenarios']['fresh_two_stage']['hours'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--provisional',action='store_true');a=p.parse_args();run(a.root,a.provisional)
