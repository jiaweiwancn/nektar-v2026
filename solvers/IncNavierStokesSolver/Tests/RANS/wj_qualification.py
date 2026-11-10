# SPDX-License-Identifier: MIT
"""Qualification is distinct from successful execution; no missing gate passes."""
import math

STATUSES=('PASS','FAIL','OBSTRUCTED','UNRESOLVED','RESOURCE_LIMIT','NOT_RUN')

def component(status,reason='',**evidence):
    if status not in STATUSES:raise ValueError('invalid component status')
    return dict(status=status,reason=reason,**evidence)

def aggregate(parts):
    if not parts:return component('NOT_RUN','no components')
    states=[p['status'] for p in parts.values()]
    for status in ('FAIL','OBSTRUCTED','RESOURCE_LIMIT','UNRESOLVED','NOT_RUN'):
        if status in states:return component(status,'required components are not all PASS',components=parts)
    return component('PASS',components=parts)

def finite(value):
    if isinstance(value,dict):return all(finite(v) for v in value.values())
    if isinstance(value,list):return all(finite(v) for v in value)
    return not isinstance(value,(int,float)) or math.isfinite(value)

def wj3_gates(s):
    """Read-only requalification of old summaries. Original receipts are untouched.

    Residual u/k/full-epsilon must decrease above an absolute 1e-11 floor.
    Resolved common-time differential and pressure differences must decrease
    with observed p in [0.5,1.5]. Constant-field roundoff classification is
    accepted only when explicitly identified and both differences fit its floor.
    """
    parts={}
    rows=s.get('spatial',[])
    if len(rows)!=3:parts['spatial']=component('NOT_RUN','three spatial levels required')
    elif not finite(rows):parts['spatial']=component('FAIL','nonfinite spatial evidence')
    else:
        tests={}
        for j,name in ((0,'u'),(2,'k'),(3,'entire_outer_epsilon')):
            v=[r['raw_reference_norms'][j] for r in rows]
            ok=all((a<=1e-11 and b<=1e-11) or 0<=b<a for a,b in zip(v,v[1:]))
            tests[name]=component('PASS' if ok else 'FAIL','decreasing residual or declared absolute floor',values=v)
        parts['spatial']=aggregate(tests)
    rates=s.get('temporal',{}).get('rates',{})
    required=('velocity','k','outer_epsilon','reduced_pressure','mechanical_pressure')
    tr={}
    for f in required:
        r=rates.get(f)
        if r is None:tr[f]=component('NOT_RUN','missing common-time differences');continue
        if not finite(r):tr[f]=component('FAIL','nonfinite temporal evidence');continue
        a,b,floor=(r.get(k) for k in ('A','B','comparison_floor'))
        if any(v is None for v in (a,b,floor)) or min(a,b,floor)<0:tr[f]=component('FAIL','invalid difference/floor');continue
        if min(a,b)<=floor:
            constant=r.get('classification')=='roundoff_limited_constant' and max(a,b)<=floor
            tr[f]=component('PASS' if constant else 'UNRESOLVED','roundoff-limited constant' if constant else 'differences not resolved',A=a,B=b,floor=floor)
        else:
            p=math.log2(a/b)
            tr[f]=component('PASS' if .5<=p<=1.5 else 'FAIL','declared first-order interval [0.5,1.5]',p=p,A=a,B=b,floor=floor)
    parts['temporal']=aggregate(tr)
    main=s.get('main',[])
    parts['trajectories']=component('OBSTRUCTED' if s.get('main_obstructions') else
        'PASS' if len(main)==3 and all('comparison' in r and not r.get('obstruction') for r in main) else 'NOT_RUN',
        'all three frozen-target trajectories required')
    parts['cross_mesh']=component('OBSTRUCTED' if s.get('cross_obstruction') else 'PASS' if s.get('cross_mesh',{}).get('comparison') else 'NOT_RUN',
        'same-time cross-mesh trajectory required')
    tx=s.get('transactions',{})
    parts['transactions']=component('PASS' if all(tx.get(k) is True for k in
        ('exact_16_plus_16','later_frozen_readback','rollback_exact')) else 'NOT_RUN','all transaction checks required')
    if not finite(s):parts['finite_evidence']=component('FAIL','nonfinite evidence')
    return aggregate(parts)
