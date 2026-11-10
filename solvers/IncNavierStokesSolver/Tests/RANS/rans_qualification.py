#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Strict reporting, coverage and physical-comparison primitives."""
import collections
import importlib.util
import json
import math
import os
from pathlib import Path
import struct

THREADS=dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
                      'BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'),'1')
os.environ.update(THREADS)

def dependencies(*names):
    missing=[name for name in names if importlib.util.find_spec(name) is None]
    if missing:raise RuntimeError('Required Python dependencies unavailable: '+', '.join(missing))

def write_json(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def coverage(rows,expected,key=lambda x:(x['model'],x['kappa_policy'],x['wall'])):
    actual=collections.Counter(key(x) for x in rows);wanted=set(expected)
    return dict(status='pass' if set(actual)==wanted and all(n==1 for n in actual.values()) else 'fail',
        missing=sorted(wanted-set(actual)),unexpected=sorted(set(actual)-wanted),
        duplicates=[k for k,n in actual.items() if n!=1],expected_count=len(wanted),actual_count=len(rows))

def exact_vector_hash(values):
    import hashlib
    if not all(math.isfinite(float(x)) for x in values):raise ValueError('nonfinite coefficient payload')
    return hashlib.sha256(struct.pack('<'+str(len(values))+'d',*values)).hexdigest()

def compare_samples(a,b,weights):
    """A shared physical grid, one global Pi gauge, weighted RMS and true L2."""
    dependencies('numpy');import numpy as np
    if not math.isclose(a['time'],b['time'],rel_tol=1e-12,abs_tol=1e-14):raise ValueError('terminal times differ')
    ga,gb=np.asarray(a['grid'],float),np.asarray(b['grid'],float)
    va,vb=np.asarray(a['grid_values'],float),np.asarray(b['grid_values'],float)
    ra,rb=np.asarray(a['grid_regions'],int),np.asarray(b['grid_regions'],int)
    w=np.asarray(weights,float)
    if ga.shape!=gb.shape or ga.ndim!=2 or ga.shape[1]!=3 or not np.array_equal(ga,gb):raise ValueError('physical coordinates differ')
    if va.shape!=vb.shape or va.shape!=(len(ga),5) or len(w)!=len(ga):raise ValueError('physical sample array lengths differ')
    if ra.shape!=rb.shape or ra.shape!=(len(w),) or not np.array_equal(ra,rb):raise ValueError('sample identities/regions differ')
    if not (np.all(np.isfinite(va)) and np.all(np.isfinite(vb)) and np.all(np.isfinite(w)) and np.all(w>0)):raise ValueError('invalid physical samples/weights')
    if 'grid_weights' in a or 'grid_weights' in b:
        if not (np.array_equal(np.asarray(a.get('grid_weights')),w) and np.array_equal(np.asarray(b.get('grid_weights')),w)):raise ValueError('physical weights differ')
    pa=float(np.dot(w,va[:,4])/w.sum());pb=float(np.dot(w,vb[:,4])/w.sum())
    va=va.copy();vb=vb.copy();va[:,4]-=pa;vb[:,4]-=pb
    result={}
    for i,name in enumerate(('u','v','k','epsilon_outer','Pi')):
        mask=ra==1 if i==3 else np.ones(len(w),dtype=bool)
        if not np.any(mask):raise ValueError('required physical region missing')
        ww=w[mask];delta=va[mask,i]-vb[mask,i];square=float(np.dot(ww,delta*delta))
        result[name]=dict(weighted_RMS=math.sqrt(max(0.,square/float(ww.sum()))),L2=math.sqrt(max(0.,square)))
    return dict(metric='weighted_RMS=sqrt(integral(error^2)/volume); L2=sqrt(integral(error^2))',
        time=a['time'],global_pressure_gauge_shifts=[pa,pb],fields=result)

def audited_extrema(report):
    """Keep the winning receipt's phase, time, location and native side."""
    winners={}
    def offer(name,value,time,phase,location=None,side=None,geometry=None):
        if value is None or not math.isfinite(value):return
        if name not in winners or (value<winners[name]['value'] or (value==winners[name]['value'] and winners[name]['location'] is None and location is not None)):
            winners[name]=dict(value=value,time=time,stage=phase,location=location,side=side,geometry=geometry)
    times=report.get('audit_times',[])
    states=report.get('state_audits',[])
    dt=report.get('dt',0.)
    for n,a in enumerate(states):
        t=times[n] if len(times)==len(states) else (report.get('time',0.)-dt*(len(states)-1-n))
        domain=a.get('accepted_awf_scalar_domain',[])
        for index,name in ((2,'inner_k'),(3,'outer_k'),(4,'outer_epsilon')):
            if len(domain)>index:offer(name,domain[index],t,'accepted',a.get('scalar_minimum_point_'+str(index-2)), 'inner' if index==2 else 'outer')
        for side,name in enumerate(('inner_covariance','outer_covariance')):
            v=a.get('accepted_awf_covariance',[])
            if len(v)>side:
                p=a.get('accepted_awf_'+('inner' if side==0 else 'outer')+'_minimum_point',[])
                provenance=a.get('native_'+('inner' if side==0 else 'outer')+'_minimum_geometry_element_side_block',[])
                offer(name,v[side],t,'accepted',p[:3] or None,provenance[2] if provenance else None,provenance[0] if provenance else None)
        r=a.get('native_A95_radicand_range',[])
        if r:
            provenance=a.get('native_A95_minimum_radicand_geometry_element_side_block',[])
            offer('A95_radicand',r[0],t,'accepted',a.get('native_A95_minimum_radicand_point',[])[:3],provenance[2] if provenance else None,provenance[0] if provenance else None)
        budget=a.get('wall_cutoff_required_budget',[]);activity=a.get('wall_cutoff_activity',[])
        if budget and activity and activity[0]>0:offer('required_cutoff_ratio',budget[1],t,'accepted',a.get('wall_cutoff_limiting_point'),'wall_profile')
    for step in report.get('pimple_attempts',[]):
        for phase,d in step.items():
            if not phase.endswith('_domain') or not isinstance(d,dict):continue
            t=step['time'];points=d.get('scalar_minimum_points_inner_k_outer_k_epsilon',[None]*3)
            for j,(name,key) in enumerate((('inner_k','minimum_inner_k'),('outer_k','minimum_outer_k'),('outer_epsilon','minimum_outer_epsilon'))):
                offer(name,d.get(key),t,phase,points[j],'inner' if j==0 else 'outer')
            for j,name in enumerate(('inner_covariance','outer_covariance')):
                v=d.get('sampled_covariance_minima',[])
                if v:offer(name,v[j],t,phase,d['covariance_minimum_points'][j],d['covariance_minimum_sides'][j],d['covariance_minimum_geometries'][j])
            offer('A95_radicand',d.get('minimum_A95_radicand'),t,phase,d.get('A95_radicand_minimum_point'),d.get('A95_radicand_side'),d.get('A95_radicand_geometry'))
            if d.get('required_cutoff_count',0)>0:offer('required_cutoff_ratio',d.get('minimum_required_budget_ratio'),t,phase,d.get('cutoff_limiting_point'))
    return dict(over_all_audited_states_and_stages=winners,final_state=states[-1] if states else None)

def required_native_checks():
    names={prefix+'-'+wall for wall in ('Accepted','EP') for prefix in
        ('compatible-restart','complete-rollback','cache-on-off','time-refinement','spatial-quadrature-refinement','mild-admissible-perturbation')}
    names.update(label+'-restart-rejection' for label in ('changed-model','changed-calibration','changed-wall'))
    names.update('unsupported-'+label for label in ('nonlinear','inner-stress','regularized','curved','3D','MPI'))
    names.add('native-EP-switch-and-baseline');return names

def native_check_coverage(rows,declared):
    base=required_native_checks()
    extra={'A95-current-gradient-cache-on-off','reviewed-CLS-checkpoint-continuation',
        'legacy-resolution-cls-reviewed-whole4','legacy-resolution-cls-reviewed-restart2',
        'unsupported-standalone-resolution-EP','unsupported-standalone-resolution-Accepted'}
    declared=set(declared or [])
    result=coverage(rows,declared,key=lambda x:x['name'])
    result['missing_required_declarations']=sorted(base-declared);result['unknown_declarations']=sorted(declared-(base|extra))
    if result['missing_required_declarations'] or result['unknown_declarations']:result['status']='fail'
    return result
