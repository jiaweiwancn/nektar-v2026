# SPDX-License-Identifier: MIT
"""Complete WJ-4 metadata gates; every required observation contributes."""
import math
from wj_qualification import component,aggregate,finite

FIELDS=('u','v','k','outer_epsilon','reduced_pressure','mechanical_pressure','velocity')
SCOPES=('full','central','constrained_boundary')

def arithmetic_floors(trajectories):
    floors={f:dict(L2=0.,sampled_Linf=0.) for f in FIELDS}
    for t in trajectories.values():
        for c in t.get('trajectory',[]):
            a=c.get('arithmetic')
            if not a:continue
            if a.get('safety_factor')!=100:raise ValueError('declared arithmetic safety factor changed')
            for f in FIELDS:
                if not finite([a['estimated_comparison_floor'][f],a['estimated_sampled_comparison_floor'][f]]) or min(a['estimated_comparison_floor'][f],a['estimated_sampled_comparison_floor'][f])<0:
                    raise ValueError('invalid arithmetic estimate')
                floors[f]['L2']=max(floors[f]['L2'],a['estimated_comparison_floor'][f])
                floors[f]['sampled_Linf']=max(floors[f]['sampled_Linf'],a['estimated_sampled_comparison_floor'][f])
    return floors

def sensitivity_sequence(comparisons,floors,design,expected,scope_uncertainty=None):
    parts={'coverage':component('PASS' if set(comparisons)==set(expected) else 'NOT_RUN','all declared common-time observations required')}
    try:
        if not finite(floors) or any(v<0 for r in floors.values() for v in r.values()):raise ValueError('invalid arithmetic estimate')
        for key,c in comparisons.items():
            tests={}
            for scope in SCOPES:
                if set(c[scope]['fields'])!=set(FIELDS):raise ValueError('incomplete comparison fields')
                for f,r in c[scope]['fields'].items():
                    amp=design['amplitudes'][f];measure=r['measure']
                    if not finite([amp,measure]) or amp<=0 or measure<=0:raise ValueError('invalid physical normalization')
                    for norm in ('L2','sampled_Linf'):
                        value=r['normalized_'+norm];u=floors[f]['sampled_Linf']/amp if scope=='constrained_boundary' or norm=='sampled_Linf' else floors[f]['L2']/(amp*math.sqrt(measure))
                        if scope_uncertainty is not None and scope in scope_uncertainty:
                            # A separately bound estimator may tighten a scope;
                            # complete physical uncertainty, not a measured-zero
                            # exemption, is required for every field and norm.
                            if set(scope_uncertainty[scope])!=set(FIELDS):raise ValueError('incomplete scope uncertainty fields')
                            physical=scope_uncertainty[scope][f][norm]
                            if not finite(physical) or physical<0:raise ValueError('invalid scope uncertainty')
                            u=physical/(amp*(math.sqrt(measure) if norm=='L2' else 1.))
                        limit=design['sensitivity_fraction']*design['budgets'][norm]
                        if not finite([value,u,limit]) or min(value,u)<0 or limit<=0:raise ValueError('nonfinite/negative comparison evidence')
                        tests[scope+'/'+f+'/'+norm]=component('PASS' if value+u<=limit else 'UNRESOLVED',
                            'difference plus finite uncertainty <=10% budget; below-floor values receive no rate',value=value,uncertainty=u,limit=limit,resolved=value>u)
            parts[str(key)]=aggregate(tests)
    except (KeyError,TypeError,ValueError) as error:parts['metadata']=component('FAIL',str(error))
    return aggregate(parts)

def wj4_gates(s,design,work_totals=None):
    parts={k:s.get('components',{}).get(k,component('NOT_RUN','required component missing')) for k in ('mesh_inventory','resource_design','transactions')}
    for k,v in s.get('components',{}).items():
        if k not in ('horizontal','longer',*parts):parts[k]=v
    if design.get('sensitivity_fraction')!=.1 or design.get('budgets',{}).get('L2')!=1e-4 or design.get('budgets',{}).get('sampled_Linf')!=1e-3:
        parts['frozen_budgets']=component('FAIL','WJ-4 declared field budgets changed')
    try:
        horizontal=s.get('horizontal',{});hd=design['horizontal'];target=s.get('components',{}).get('horizontal',{}).get('components',{}).get('target',component('NOT_RUN'))
        comparisons=s.get('horizontal_comparisons',{})
        parts['horizontal']=aggregate(dict(target=target,sensitivity=sensitivity_sequence(comparisons,arithmetic_floors(horizontal),design,
            [str(i) for i in range(0,int(hd['steps'])+1,32)])))
        trajectories=s.get('longer',{});ld=design['longer'];targets={};arithmetic={}
        for j,tag in enumerate(('coarse','fine')):
            t=trajectories.get(tag,{});rows=t.get('trajectory',[]);n=int(ld['steps'][j]);block=int(ld['cadence'][j]);dt=ld['dt'][j]
            expected=list(range(block,n+1,block));actual=[c['steps'] for c in rows]
            ok=actual==expected and not t.get('hard_obstruction') and not t.get('accuracy_obstruction') and t.get('initial',{}).get('errors',{}).get('accuracy_gate',{}).get('status')=='PASS' and all(abs(c['time']-c['steps']*dt)<=2e-15 and c['errors']['accuracy_gate']['status']=='PASS' for c in rows)
            targets[tag]=component('PASS' if ok else 'OBSTRUCTED' if t.get('hard_obstruction') or t.get('accuracy_obstruction') else 'NOT_RUN','complete frozen target and every measured accuracy gate required')
            for step in ld['arithmetic_steps'][j]:
                found=next((c for c in rows if c['steps']==step and c.get('arithmetic')),None)
                arithmetic[tag+'/'+str(step)]=component('PASS' if found else 'NOT_RUN','declared early/middle/end passive arithmetic packet required')
        paired=s.get('timestep_comparisons',{})
        for c in paired.values():
            if abs(c['time']-c['fine_time'])>2e-15:raise ValueError('common physical times differ')
        count=int(ld['steps'][0])//int(ld['cadence'][0])
        parts['longer']=aggregate(dict(target=aggregate(targets),arithmetic=aggregate(arithmetic),
            sensitivity=sensitivity_sequence({k:v['comparison'] for k,v in paired.items()},arithmetic_floors(trajectories),design,[str(i) for i in range(count+1)])))
        a,b=(trajectories.get(tag,{}) for tag in ('coarse','fine'))
        required=('initial_primary_sha256','initial_history_sha256','native_dt_signature_sha256')
        present=all(k in t for t in (a,b) for k in required)
        equal=present and a['initial_primary_sha256']==b['initial_primary_sha256'] and a['initial_history_sha256']==b['initial_history_sha256'] and a['native_dt_signature_sha256']!=b['native_dt_signature_sha256']
        parts['fresh_startup']=component('PASS' if equal else 'FAIL' if present else 'NOT_RUN','equal physical startup fields/history and distinct native dt signatures required')
        tx=parts['transactions']
        if tx['status']=='PASS' and not all(tx.get(k) is True for k in ('exact_primary_gauge_history_time','exact_accepted_audit_samples')):
            parts['transactions']=component('FAIL','incomplete exact transaction metadata')
    except (KeyError,TypeError,ValueError) as error:parts['complete_metadata']=component('FAIL',str(error))
    if not finite(s):parts['finite_evidence']=component('FAIL','nonfinite evidence')
    totals=work_totals if work_totals is not None else s.get('ledger_totals')
    parts['actual_resource']=component('NOT_RUN','missing actual work ledger') if not totals else component(
        'PASS' if finite(totals) and totals['attempts']<=design['attempt_cap_per_policy'] and totals['seconds']<=design['wall_cap_seconds_per_policy'] else 'RESOURCE_LIMIT',totals=totals)
    return aggregate(parts)
