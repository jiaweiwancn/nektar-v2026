# SPDX-License-Identifier: MIT
"""Explicit three-level time map; absolute screens and contraction are separate."""
import math
import numpy as np
from wj_h19_analysis import metrics,principal_screen,fixed_action,scales,LIMITS,FLOOR,GATE
from wj_h7_analysis import context
from wj_h4_points import aligned

DT=dict(C=.005,M=.0025,F=.00125)
ADDITIONAL=dict(M=4,C=2,F=16)
START=dict(M=4,C=2,F=0)
END=dict(C=4,M=8,F=16)
STRIDE=dict(C=1,M=2,F=4)
TIMES=(.005,.01,.015,.02)
PROBES=(dict(name='outer_principal',physical=[15.,0.,.25607517888672315],geometry=251,x_side='upper',normal_side='upper'),
        dict(name='inner_stress',physical=[-7.5,0.,.1],geometry=57,x_side='lower',normal_side='upper'))

def absolute_step(level,t):
    assert level in DT
    matches=[v for v in TIMES if abs(t-v)<=128*np.finfo(float).eps*max(1.,abs(t),abs(v))]
    assert len(matches)==1,'time is not a declared physical comparison time'
    n=round(matches[0]/DT[level]);aligned(n,DT[level],t);return n

def selected(report,level,t):
    assert report['dt']==DT[level]
    n=absolute_step(level,t);c=context(report,n)
    assert c['step']==n;aligned(n,DT[level],c['time'])
    assert c['alpha']==1/DT[level]and abs(c['time']-t)<=128*np.finfo(float).eps*max(1.,abs(t))
    return c

def availability(phase):
    assert phase in ('M','C','F')
    return dict(M='NOT_YET_AVAILABLE',C='PAIR_C_M_AVAILABLE',F='TRIPLET_AVAILABLE')[phase]

def ratio(d1,d2,floor=FLOOR):
    assert np.isfinite(d1)and np.isfinite(d2)and d1>=0 and d2>=0
    out=dict(D1=float(d1),D2=float(d2),floor=float(floor),rho=None,upper_bound=None,p_obs=None)
    if d1<=floor:return dict(out,status='UNRESOLVED',reason='D1 at/below reconstruction floor')
    if d2<=floor:return dict(out,status='REDUCED_AT_RESOLUTION',upper_bound=float(floor/d1),reason='D2 below floor; no precise rho/order')
    r=d2/d1
    return dict(out,status='CONTRACTING'if r<1 else'NONCONTRACTING',rho=float(r),p_obs=float(-math.log2(r)),reason='observed diagnostic only')

def continuation(physical,absolute,contraction=None):
    # Contraction is descriptive and deliberately not an admission law.
    if physical!='PASS':return 'PHYSICAL_STOP'
    if absolute not in ('PASS','NOT_YET_AVAILABLE'):return 'REVIEW_STOP'
    return 'PROCEED'

def all_pair_screens(rows,data):
    pairs={};arrays={}
    for a,b in (('C','M'),('M','F'),('C','F')):
        m=metrics(rows,data[a],data[b]);p,v=principal_screen(data[a],data[b],rows.outer)
        failed=any(x[k]['status']=='EXCEEDED'for region in m.values()for x in region.values()for k in ('rms','maximum'))
        status='UNRESOLVED'if p['status']=='UNRESOLVED'else'EXCEEDED'if failed or p['status']=='EXCEEDED'else'PASS'
        pairs[a+'-'+b]=dict(status=status,temporal=m,principal=p)
        arrays.update({a+b+'_'+k:value for k,value in v.items()})
    status='UNRESOLVED'if any(p['status']=='UNRESOLVED'for p in pairs.values())else'EXCEEDED'if any(p['status']=='EXCEEDED'for p in pairs.values())else'PASS'
    return dict(status=status,pairs=pairs),arrays

def norm_pair(rows,e1,e2,scale,outer_only=False):
    result={};n=len(rows.rows)
    for region,mask in (('full',np.ones(n,bool)),('inner',~rows.outer),('outer',rows.outer)):
        if outer_only and region=='inner':continue
        active=mask&(rows.outer if outer_only else True);w=rows.weights[active];measure=float(w.sum())
        s=scale if np.isscalar(scale)else scale[active]
        d1=np.linalg.norm(e1.reshape(n,-1)[active],axis=1)/s;d2=np.linalg.norm(e2.reshape(n,-1)[active],axis=1)/s
        result[region]=dict(RMS=ratio(math.sqrt(float(w@d1**2)/measure),math.sqrt(float(w@d2**2)/measure)),maximum=ratio(float(d1.max()),float(d2.max())),measure=measure)
    return result

def shared_direction_actions(data,ids):
    actions={}
    for source in ('C','M','F'):
        angle=data[source]['angle'][ids]
        for level in ('C','M','F'):
            d=data[level]
            actions[source+level]=fixed_action(d['k'][ids],d['epsilon'][ids],d['G'][ids],angle)
    return actions

def triplet_diagnostics(rows,data):
    scalar=scales(rows);fields={};arrays={}
    for name in LIMITS:
        e1=data['C'][name]-data['M'][name];e2=data['M'][name]-data['F'][name]
        fields[name]=norm_pair(rows,e1,e2,scalar[name],name=='epsilon')
        arrays[name+'_E1']=e1;arrays[name+'_E2']=e2;arrays[name+'_E1_minus_2E2']=e1-2*e2
    ids=np.flatnonzero(rows.outer);principal={};actions=shared_direction_actions(data,ids)
    for source in ('C','M','F'):
        for level in ('C','M','F'):
            a=actions[source+level]
            for key,value in a.items():arrays[f'theta_{source}_{level}_{key}']=value
        e1=actions[source+'C']['normalized']-actions[source+'M']['normalized'];e2=actions[source+'M']['normalized']-actions[source+'F']['normalized']
        x=np.zeros(len(rows.rows));y=x.copy();x[ids]=e1;y[ids]=e2
        principal['theta_'+source]=norm_pair(rows,x,y,1.,True)
        arrays['theta_'+source+'_E1']=e1;arrays['theta_'+source+'_E2']=e2
    e1=data['C']['margin']-data['M']['margin'];e2=data['M']['margin']-data['F']['margin']
    principal['minimized']=norm_pair(rows,e1,e2,1.,True)
    arrays['minimized_E1']=e1;arrays['minimized_E2']=e2;arrays['outer_ids']=ids
    probes=[]
    def finite(v):
        a=np.asarray(v);return a.tolist()if np.all(np.isfinite(a))else None
    for probe in PROBES:
        pbits=np.asarray(probe['physical'],dtype='<f8').tobytes()
        same=[i for i,(p,o)in enumerate(rows.rows[:rows.point_count])if np.asarray(p,dtype='<f8').tobytes()==pbits]
        expected=[i for i in same if all(rows.rows[i][1][k]==probe[k]for k in ('geometry','x_side','normal_side'))]
        assert len(expected)==1,'fixed probe exact owner absent/wrong'
        for i in same:
            values={level:{name:finite(data[level][name][i])for name in ('U','G','k','epsilon','R','pi','production','margin','angle','nu_eff','K','rc','destruction','transverse','covariance','root')}for level in ('C','M','F')}
            if not rows.outer[i]:
                for d in values.values():d['epsilon']=None
            contributions={}
            for name in LIMITS:
                if name=='epsilon'and not rows.outer[i]:continue
                a=arrays[name+'_E1'][i];b=arrays[name+'_E2'][i];s=scalar[name]if np.isscalar(scalar[name])else scalar[name][i]
                total=abs(a)+abs(b);cancel=np.divide(abs(a+b),total,out=np.full_like(np.asarray(total,dtype=float),np.nan),where=total>FLOOR*s)
                contributions[name]=dict(E1=finite(a),E2=finite(b),direct_C_F=finite(a+b),cancellation=finite(1-cancel),component_floor=FLOOR*float(s))
            directions={}
            if rows.outer[i]:
                j=int(np.flatnonzero(ids==i)[0])
                for source in ('C','M','F'):
                    a=actions[source+'C']['normalized'][j];b=actions[source+'M']['normalized'][j];c=actions[source+'F']['normalized'][j]
                    directions['theta_'+source]=dict(values={level:{key:finite(actions[source+level][key][j])for key in ('dimensional','nu_eff','normalized','xi','a')}for level in ('C','M','F')},
                        signed_E1=float(a-b),signed_E2=float(b-c),contraction=ratio(abs(a-b),abs(b-c)))
                a,b,c=(data[level]['margin'][i]for level in ('C','M','F'))
                directions['minimized']=dict(signed_E1=float(a-b),signed_E2=float(b-c),contraction=ratio(abs(a-b),abs(b-c)))
            probes.append(dict(name=probe['name'],primary=i==expected[0],row=i,physical=rows.rows[i][0],owner=rows.rows[i][1],states=values,components=contributions,principal=directions))
    return dict(fields=fields,principal=principal,fixed_probes=probes,scope='same support/weights/scales; rho and p_obs descriptive, never a stop gate'),arrays
