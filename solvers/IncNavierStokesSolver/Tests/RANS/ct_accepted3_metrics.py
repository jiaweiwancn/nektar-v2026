# SPDX-License-Identifier: MIT
"""Read-only model-aware inventories on declared native and physical rows."""
import os
from rans_qualification import THREADS
os.environ.update(THREADS)
import numpy as np
from ct_accepted_open import row_value
from ct_accepted2_diagnosis import at as a95_at, grids as ct2_grids, WITNESS
from validate_a93_a95_e94e95 import oracle


def point(row, state, offsets, model):
    if model=='A95':return a95_at(row,state,offsets)
    G=np.zeros((3,3))
    for b,a in enumerate((0,2)):G[a]=row_value(row['velocity_k'],state,offsets[b])[1]
    k=row_value(row['velocity_k'],state,offsets[2])[0]
    e=row_value(row['epsilon'],state,offsets[3])[0]
    p,dp=row_value(row['pressure'],state,offsets[5])
    S=(G+G.T)/2-np.eye(3)*np.trace(G)/3;O=(G-G.T)/2
    v=dict(raw_gradient=G.tolist(),raw_divergence=float(np.trace(G)),k=float(k),epsilon=float(e),
           tau=float(k/e),strain_eigenvalues=np.linalg.eigvalsh(S).tolist(),strain_norm=float(np.linalg.norm(S)),
           rotation_norm=float(np.linalg.norm(O)),Pi=float(p),pressure_gradient=dp.tolist(),
           isotropic_production_work=float(-2*k*np.trace(G)/3),radicand=None)
    try:
        law=oracle(model,k,e,G)
        if k<=0 or e<=0:raise ValueError('nonpositive actual outer scalar')
        T=-law['R']+np.eye(3)*(2*k/3-p)
        v.update(status='admissible',R=law['R'].tolist(),Q=law['Q'].tolist(),production=float(law['production']),
                 total_stress=T.tolist(),nu_t=float(law['nut']))
    except (ValueError,ZeroDivisionError,FloatingPointError) as err:
        v.update(status='constitutive_domain_failure',reason=str(err),full_stress_and_production='undefined')
    return v


def inventory(rows,coords,state,offsets,model,sides=None):
    extrema={};rad=None;failures=[];d2=0.;count=0;divmax=0.;inner=None
    for i,q in enumerate(rows):
        if q.get('epsilon') is None:
            k=float(row_value(q['velocity_k'],state,offsets[2])[0])
            loc=dict(coordinate=coords[i],geometry=q['velocity_k']['geometry'],side=sides[i] if sides is not None else 0,sample_index=i)
            if inner is None or k<inner['value']:inner=dict(loc,value=k)
            if not np.isfinite(k) or k<0:failures.append(dict(loc,status='actual_inner_k_domain_failure',k=k))
            continue
        v=point(q,state,offsets,model);count+=1
        loc=dict(coordinate=coords[i],geometry=q['velocity_k']['geometry'],side=sides[i] if sides is not None else 0,sample_index=i)
        d2+=v['raw_divergence']**2;divmax=max(divmax,abs(v['raw_divergence']))
        for name,value in [('k',v['k']),('epsilon',v['epsilon'])]+([('covariance',float(np.linalg.eigvalsh(v['R'])[0]))] if 'R' in v else []):
            if name not in extrema or value<extrema[name]['value']:extrema[name]=dict(loc,value=value)
        if v.get('radicand') is not None and (rad is None or v['radicand']<rad['radicand']):rad=dict(v,**loc)
        if v['status']!='admissible':failures.append(dict(v,**loc))
    return dict(count=count,domain_failure_count=len(failures),first_failures=failures[:4],minimum_radicand=rad,minimum_inner_k=inner,
                minima=extrema,divergence_unweighted_RMS=float(np.sqrt(d2/count)) if count else None,
                divergence_maximum=divmax)


def grid(c):
    points,weights,n=ct2_grids(c)
    # Six points integrate squares of represented P4 fields exactly on every
    # union cell. CT2's original P3/four-point packets remain unchanged.
    tail=points[n:]
    from numpy.polynomial.legendre import leggauss
    g,w=leggauss(6);points=[];weights=[]
    xs=[0.,.125,.1875,.25];zs=sorted(set(c['knots']+[5.125,6.5]))
    for a,b in zip(xs[:-1],xs[1:]):
        for d,e in zip(zs[:-1],zs[1:]):
            for i in range(6):
                for j in range(6):
                    points.append([(a+b)/2+(b-a)*g[i]/2,0.,(d+e)/2+(e-d)*g[j]/2])
                    weights.append((b-a)*(e-d)*w[i]*w[j]/4)
    n=len(points);points+=tail
    # CT2's final entry was the original witness; keep it LAST for projection.
    points.pop();sides=[0]*len(points)
    for x in (0.,1e-8,.0625,.125-1e-8,.125+1e-8,.1875-1e-8,.1875+1e-8,.22,.24836915803232323,.25):
        for z in sorted(set(np.linspace(5.8,6.35,65).tolist()+[6.109375,WITNESS[2]])):
            points.append([x,0.,z]);sides.append(0)
        for side in (-1,1):points.append([x,0.,c['H']]);sides.append(side)
    points.append(WITNESS);sides.append(0)
    return points,sides,weights,n


def domain(initial,report,c):
    o=initial['observation'];s=np.array(report['state']);off=o['offsets']
    result=dict(native=inventory(o['volume'],[q['coordinate'] for q in o['volume']],s,off,c['model']),
                fixed=inventory(initial['test_basis'],o['grid'],s,off,c['model'],o['grid_sides']))
    result['status']='pass' if all(x['domain_failure_count']==0 and 'covariance' in x['minima']
        and np.isfinite(x['minima']['covariance']['value']) and x['minima']['covariance']['value']>=-1e-10*c['K0']
        for x in (result['native'],result['fixed'])) else 'fail'
    return result
