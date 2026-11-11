# SPDX-License-Identifier: MIT
"""Streamed, both-order WJ prefix qualification with exact tuple coverage."""
import copy
import gc
import math
from pathlib import Path
import numpy as np
from wj_immutable import read,file_hash,verify_artifacts,digest
from wj_qualification import aggregate,component
from wj3_measure import FIELDS
from wj4_measure import evaluate,comparison,BUDGET
from wj4_qualification import sensitivity_sequence,arithmetic_floors
from wj5_measure import boundary_uncertainty
from wj6_qualification import individual_arithmetic,sum_envelopes
from run_cls_epawf_picard import require

TAGS=('coarse','fine','baseline')
SCOPES=('full','central','constrained_boundary')
NAMES=FIELDS+('velocity',)


def observation_steps(tag,end):
    old=[0,4,16,24,32] if tag=='coarse' else [0,8,32,48,64] if tag=='fine' else [0,4,16,32]
    return old+[v*(2 if tag=='fine' else 1) for v in (64,96,128) if v<=end]


def expected_tuples(end):
    return {(order,group,step,scope,field,norm) for order in (8,10)
        for group in ('horizontal','timestep') for step in (0,4,32,64,96,128) if step<=end
        for scope in SCOPES for field in NAMES for norm in ('L2','sampled_Linf')}


def tuple_gate(rows,end):
    keys=[tuple(r['key']) for r in rows]
    ok=len(keys)==len(set(keys)) and set(keys)==expected_tuples(end)
    return component('PASS' if ok else 'NOT_RUN','explicit unique order/group/time/scope/field/norm coverage',
        expected=len(expected_tuples(end)),actual=len(keys),unique=len(set(keys)))


def joint_ready(rows,end):
    return set(rows)=={'EP','Budget'} and all(r.get('end')==end and r.get('gate',{}).get('status')=='PASS' for r in rows.values())


def final_gate(rows):
    return component('PASS' if joint_ready(rows,128) else 'NOT_RUN','both policies must complete fixed T3 and every prefix gate')


def load_measurement(c):
    require(file_hash(c['path'])==c['sha256'],'stale measurement bytes')
    require(file_hash(c['report'])==c['report_sha256'],'wrong/stale producer')
    producer=read(c['report'])
    require(producer['steps']==c['steps'] and producer['time']==c['time'],'absolute native step/time changed')
    require(file_hash(Path(c['checkpoint'])/'State.bin')==c['checkpoint_sha256'],'stale checkpoint')
    m=read(c['path'])['measurements']
    require(m['state_sha256']==c['state_sha256']==m['state_after_sha256'] and
        m['stage_history_sha256']==c['history_sha256']==m['stage_history_after_sha256'] and
        m['time']==c['time'] and m['governing_solves']==0,'wrong state/history/time or governing measurement solve')
    require(all(r.get('constraint_rows') for r in m['observations']),'missing native observation row maps')
    if c['steps']:
        require(m['completed_stage'] and m['completed_stage']['time']==producer['last_stage_time'] and
            m['completed_stage']['alpha_mass']==producer['last_stage_alpha'] and
            m['completed_stage']['history_sha256']==__import__('run_cls_epawf_picard').state_hash(producer['last_stage_history']),
            'wrong completed-stage producer/history/time/alpha')
    return m


def complete_errors(m,policy,c):
    result=evaluate(m,policy,xorder=c['rule_settings']['xorder'],rule_settings=c['rule_settings'])
    tests={}
    for kind in ('continuum_error','initial_drift'):
        for scope in SCOPES:
            require(set(result[kind][scope]['fields'])==set(NAMES),'incomplete accuracy fields')
            for field,row in result[kind][scope]['fields'].items():
                for norm in ('L2','sampled_Linf'):
                    value=row['normalized_'+norm];limit=BUDGET[norm]
                    tests[kind+'/'+scope+'/'+field+'/'+norm]=component('PASS' if math.isfinite(value) and value<=limit else 'OBSTRUCTED',value=value,limit=limit)
    result['accuracy_gate']=aggregate(dict(original=result['accuracy_gate'],all_scopes=aggregate(tests)))
    return result


def qualify(index,policy,design,end):
    """Only one full payload, or the current comparison pair, is live.

    The index consists of scalar records/immutable identities. Row-map passes
    release every full payload; comparisons release each pair before the next.
    Arithmetic formulas and comparison ordering are unchanged from WJ-6.
    """
    expected={(tag,step,order) for tag in TAGS for step in observation_steps(tag,end) for order in (8,10)}
    require(set(index)==expected,'missing/duplicate required observation order/time')
    require(all(c['steps']==step and c['rule_settings']['order']==order for (tag,step,order),c in index.items()),
            'wrong observation role/order/absolute step')
    parts={};metrics={};envelopes={};sources={}
    for tag in TAGS:
        records=[index[tag,step,8] for step in observation_steps(tag,end)]
        packets=[copy.deepcopy(c['arithmetic']) for c in records if c.get('arithmetic')]
        required_new=[v*(2 if tag=='fine' else 1) for v in (64,96,128) if v<=end]
        require(all(index[tag,v,8].get('arithmetic') for v in required_new),'missing cadence absolute-step arithmetic')
        require(packets,'missing retained arithmetic')
        basis=np.zeros(6);evaluation=np.zeros(6)
        for step in observation_steps(tag,end):
            for order in (8,10):
                c=index[tag,step,order];m=load_measurement(c)
                errors=complete_errors(m,policy,c)
                require(errors==c['errors'],'replayed accuracy/interface/divergence decisions differ')
                metrics[tag,step,order]=errors
                parts['accuracy/'+tag+'/'+str(step)+'/'+str(order)]=errors['accuracy_gate']
                basis=np.maximum(basis,m['evaluation_basis_row_norm_max'])
                evaluation=np.maximum(evaluation,m['arithmetic_evaluation_sensitivity'])
                del m;gc.collect()
            a,b=(metrics[tag,step,o] for o in (8,10))
            changes=[abs(a[k][s]['fields'][f]['L2']-b[k][s]['fields'][f]['L2'])
                for k in ('continuum_error','initial_drift') for s in SCOPES for f in NAMES]
            parts['quadrature/'+tag+'/'+str(step)]=component('PASS' if max(changes)<=1e-10 else 'UNRESOLVED',maximum=max(changes),limit=1e-10)
        minimum={name:min(p['systems'][name]['physical_metric_min_eigenvalue'] for p in packets) for name in ('momentum','k','epsilon')}
        pressure=min(p['systems']['pressure']['physical_mass_operator_bound'] for p in packets)
        factors=[basis[i]/math.sqrt(minimum[name]) for i,name in enumerate(('momentum','momentum','k','epsilon'))]+[basis[4]/math.sqrt(pressure)]*2
        converted=[]
        for p in packets:
            verify_artifacts(p['input_lineage']);verify_artifacts(p['code'])
            p['sampled_conversion_factors']=factors;p['evaluation_sensitivity']=evaluation.tolist()
            converted.append(dict(arithmetic=individual_arithmetic(p)))
        floors=arithmetic_floors(dict(a=dict(trajectory=converted)))
        systems={name:dict(physical_metric_min_eigenvalue=value) for name,value in minimum.items()}
        systems['pressure']=dict(physical_mass_operator_bound=pressure)
        envelope=dict(systems=systems,estimated_comparison_floor={f:v['L2'] for f,v in floors.items()},
            estimated_sampled_comparison_floor={f:v['sampled_Linf'] for f,v in floors.items()})
        boundary={f:dict(L2=0.,sampled_Linf=0.) for f in floors};coverage=[]
        for step in observation_steps(tag,end):
            for order in (8,10):
                c=index[tag,step,order];m=load_measurement(c)
                estimate=boundary_uncertainty(m,envelope,c['rule_settings'])
                for f,row in estimate['fields'].items():
                    for norm in ('L2','sampled_Linf'):boundary[f][norm]=max(boundary[f][norm],row[norm])
                coverage.append(dict(step=step,order=order,state=c['state_sha256'],measurement_sha256=c['sha256']))
                del m;gc.collect()
        envelopes[tag]=dict(floors=floors,boundary=boundary,coefficient_envelope=envelope,coverage=coverage,
            maximum_basis=basis.tolist(),evaluation_sensitivity=evaluation.tolist())
    rows=[];comparisons={}
    for order in (8,10):
        for group,ta,tb in (('horizontal','baseline','coarse'),('timestep','coarse','fine')):
            floors=sum_envelopes(envelopes[ta]['floors'],envelopes[tb]['floors'])
            boundary=sum_envelopes(envelopes[ta]['boundary'],envelopes[tb]['boundary'])
            diffs={}
            for step in (0,4,32,64,96,128):
                if step>end:continue
                ca=index[ta,step,order];cb=index[tb,step*2 if tb=='fine' else step,order]
                a=load_measurement(ca);b=load_measurement(cb)
                require(abs(a['time']-b['time'])<=2e-15,'common physical times differ')
                diffs[str(step)]=comparison(a,b,policy,ca['rule_settings'])
                del a,b;gc.collect()
            gate=sensitivity_sequence(diffs,floors,design,list(diffs),scope_uncertainty={'constrained_boundary':boundary})
            parts[group+'/'+str(order)]=gate;comparisons[group+'/'+str(order)]=diffs
            for step,g in gate['components'].items():
                if not step.isdigit():continue
                for key,v in g['components'].items():
                    scope,field,norm=key.split('/')
                    rows.append(dict(key=[order,group,int(step),scope,field,norm],**v))
    parts['tuple_coverage']=tuple_gate(rows,end)
    return dict(schema='WJ7StreamedPrefixQualificationV1',policy=policy,end=end,gate=aggregate(parts),
        sensitivity=rows,envelopes=envelopes,comparisons=comparisons,
        memory_contract='one full row-map payload or current comparison pair/order; no all-history payload cache',
        record_index={tag+'/'+str(step)+'/'+str(order):dict(path=c['path'],sha256=c['sha256'],report=c['report'],report_sha256=c['report_sha256'])
            for (tag,step,order),c in index.items()})
