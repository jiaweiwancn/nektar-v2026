# SPDX-License-Identifier: MIT
"""Sufficient statistics of actual native boundary rows, not a new estimator."""
import math
from pathlib import Path
import numpy as np
from wj_immutable import read,file_identity,write_once,verify_artifacts,digest
from wj4_measure import boundary_rule
from wj5_measure import constraint_estimate
from wj3_measure import FIELDS
from run_cls_epawf_picard import require


def build(m,c):
    ids={int(i) for r in m['dependent_trace_action'] for i,value in r['free_K_action']}
    mean={int(r['offset'])+i:float(v) for r in m['passive_diagnostics']['pressure_regions'] for i,v in enumerate(r['gauge_mean_row'])}
    weights=boundary_rule(**{k:v for k,v in c['rule_settings'].items() if k!='all_faces'})
    rows=[];seen=set()
    for r in m['observations']:
        x,_,z=r['position']
        if not (x in (0.,1.) or z==2):continue
        desc=[]
        for j,d in enumerate(r['constraint_rows'][:4]):
            if j==3 and r['inner']:desc.append(None);continue
            item={k:d[k] for k in ('free_row_L2','gross_row_L2','absolute_accumulation','prescribed_accumulation')}
            item['n']=len(d['row'])
            if j==3:
                item['krow']=math.sqrt(sum(v*v for i,v in d['free_row'] if int(i) in ids))
                item['erow']=math.sqrt(sum(v*v for i,v in d['free_row'] if int(i) not in ids))
            desc.append(item)
        pd=r['constraint_rows'][4];p={i:-v for i,v in mean.items()}
        for i,v in pd['row']:p[int(i)]=p.get(int(i),0.)+v
        key=(x,z);w=weights.get(key,0.) if key not in seen else 0.;seen.add(key)
        rows.append(dict(position=r['position'],inner=r['inner'],geometry=r['geometry'],weight=w,desc=desc,
            owner={k:r[k] for k in ('element','side','requested_x_side','requested_normal_side','face_sides_xz','provenance') if k in r},
            pressure_row_L2=math.sqrt(sum(v*v for v in p.values())),pressure_n=len(p),pressure_accumulation=pd['absolute_accumulation']))
    return dict(schema='WJ7NativeRowSufficientStatisticsV1',source=file_identity(c['path']),producer=file_identity(c['report']),steps=c['steps'],
        state_sha256=m['state_sha256'],history_sha256=m['stage_history_sha256'],time=m['time'],
        completed_stage_context={k:m['completed_stage'][k] for k in ('time','alpha_mass','history_sha256')} if m['completed_stage'] else None,
        basis=m['evaluation_basis_row_norm_max'],evaluation=m['arithmetic_evaluation_sensitivity'],
        full_row_count=len(m['observations']),boundary_rows=rows,errors=c['errors'],
        scope='exact sufficient statistics of actual elimination/trace rows; full-volume pressure gauge retained; interior rows have zero boundary mask')


def estimate(cache,envelope):
    systems=envelope['systems'];floor=envelope['estimated_comparison_floor'];sampled=envelope['estimated_sampled_comparison_floor']
    values=[];maxfree=np.zeros(4)
    for r in cache['boundary_rows']:
        v=np.zeros(6)
        for j,name in enumerate(FIELDS[:4]):
            d=r['desc'][j]
            if d is None:continue
            block='momentum' if j<2 else 'k' if j==2 else 'epsilon'
            coefficient=floor[name]/math.sqrt(systems[block]['physical_metric_min_eigenvalue'])
            free=None
            if j==3:
                kcoefficient=floor['k']/math.sqrt(systems['k']['physical_metric_min_eigenvalue'])
                free=d['krow']*kcoefficient+d['erow']*coefficient
            virtual=dict(d,row=[None]*int(d['n']))
            v[j]=constraint_estimate(virtual,coefficient,evaluation_uncertainty=cache['evaluation'][j],free_state_bound=free)['total']
            maxfree[j]=max(maxfree[j],d['free_row_L2'])
        accum=2*100*np.finfo(float).eps*(32*r['pressure_n']+16)*(1+r['pressure_accumulation'])
        coefficient=floor['reduced_pressure']/math.sqrt(systems['pressure']['physical_mass_operator_bound'])
        v[4]=max(sampled['reduced_pressure'],r['pressure_row_L2']*coefficient+accum);v[5]=v[4]+2*v[2]/3
        values.append(v)
    a=np.array(values);w=np.array([r['weight'] for r in cache['boundary_rows']]);outer=np.array([not r['inner'] for r in cache['boundary_rows']]);fields={}
    for j,f in enumerate(FIELDS+('velocity',)):
        active=outer if f=='outer_epsilon' else np.ones(len(w),bool)
        v=np.hypot(a[:,0],a[:,1]) if f=='velocity' else a[:,j]
        fields[f]=dict(L2=float(np.sqrt(w[active]@(v[active]**2))),sampled_Linf=float(max(v[active])),measure=float(sum(w[active])))
    return dict(fields=fields,maximum_boundary_free_rows=dict(zip(FIELDS[:4],maxfree.tolist())))


def verify(cache,c):
    verify_artifacts({cache['source']['path']:cache['source'],cache['producer']['path']:cache['producer']})
    require(cache['source']['path']==c['path'] and cache['source']['sha256']==c['sha256'] and
        cache['producer']['sha256']==c['report_sha256'] and cache['state_sha256']==c['state_sha256'] and
        cache['history_sha256']==c['history_sha256'] and cache['time']==c['time'] and cache['errors']==c['errors'],
        'stale/wrong native row sufficient statistics')
    report=read(c['report'])
    require(report['steps']==cache['steps']==c['steps'] and report['time']==cache['time'],'cache absolute native step/time mismatch')
    require(file_identity(Path(c['checkpoint'])/'State.bin')['sha256']==c['checkpoint_sha256'],'cache checkpoint changed')
    if c['steps']:
        from run_cls_epawf_picard import state_hash
        stage=cache['completed_stage_context']
        require(stage and stage['time']==report['last_stage_time'] and stage['alpha_mass']==report['last_stage_alpha'] and
            stage['history_sha256']==state_hash(report['last_stage_history']),'cache wrong completed-stage producer')
