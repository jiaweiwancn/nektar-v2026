#!/usr/bin/env python3
"""Independent small convex-QP reference, followed separately by the native full oracle.

Uses SciPy SLSQP in a rank-revealing QR subspace and independently polished
80-digit equality KKT equations. No native dual exchange code is imported.
The returned candidate is only a finite-pool reference until the full original
separation oracle, including continuous wall extrema, verifies it.
"""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import platform
import sys
import time
sys.dont_write_bytecode=True

def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result);return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--epoch',type=Path,required=True)
    p.add_argument('--case',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    assert platform.system()=='Linux' and len(os.sched_getaffinity(0))==1
    assert all(os.environ.get(x)=='1' for x in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'))
    import numpy as np
    import mpmath as mp
    from scipy.linalg import qr,solve_triangular,lu_factor,lu_solve
    from scipy.optimize import minimize
    begin=time.monotonic();e=a.epoch.resolve();case=a.case.resolve()
    evidence=e/'source/docs/plans/cls-epawf-revised-evidence'
    h=module(evidence/'m8-r4/prepare.py','r12_ref')
    events=[json.loads(line) for line in (case/'case.r3-progress.jsonl').read_text().splitlines()]
    record,=[x for x in events if x.get('kind')=='r6_cut_audit' and x['name']=='r12_Lf' and x['status']=='passed_original_cut_scan']
    witness,=[x for x in events if x.get('kind')=='r12_candidate' and x['name']=='L1_witness']
    manifest=h.load(case/'predictor.json')
    n=len(h.vector(manifest['state']))
    J=np.fromfile(h.bound(manifest['frozen_action']['matrix']),dtype='<f8').reshape(n,n)
    R=np.array(h.vector(manifest['origin']['residual']))
    con={int(k):v for k,v in np.fromfile(h.bound(manifest['origin']['constraints']),dtype='<f8').reshape(-1,2)}
    sys.path.insert(0,str(e/'source/solvers/IncNavierStokesSolver/Tests/RANS'))
    from cls_epawf_m8_r3_resolution import ResolutionAudit
    raw=h.load(h.bound(manifest['current_preflight']));metric=ResolutionAudit(raw,e/'source')
    rf=R.copy();free=np.zeros(n,bool);B=np.zeros((n,n))
    for block,field in enumerate(metric.fields[:raw['dimension']+2]):
        local=np.array([i for i in range(field.size) if metric.off[block]+i not in con],int)
        ids=metric.off[block]+local;free[ids]=True
        L=np.linalg.cholesky(metric.mass[block][np.ix_(local,local)])
        scale=metric.state_scale(block)*metric.us/metric.H*np.sqrt(field.volume)
        B[np.ix_(ids,ids)]=scale*L
    rf[free]/=.125
    dseed=np.array(h.vector(manifest['r12']['seed']))
    mixed=rf.astype(np.longdouble)+J.astype(np.longdouble)@dseed.astype(np.longdouble)
    c=np.zeros(n)
    for block,field in enumerate(metric.fields[:raw['dimension']+2]):
        ids=np.array([metric.off[block]+i for i in range(field.size) if metric.off[block]+i not in con],int)
        c[ids]=solve_triangular(B[np.ix_(ids,ids)],np.array(mixed[ids],float),lower=True)
    seed_obj=.5*float(c@c)
    assert abs(seed_obj-272.9590401766774)<=2e-9
    rows=record['rows'];N=np.array([h.vector(x['n']) for x in rows]);beta=np.array([x['transformed_rhs'] for x in rows]);length=np.linalg.norm(N,axis=1)
    original=np.array([h.vector(x['a']) for x in rows]);ell=np.array([x['rhs'] for x in rows])
    N/=length[:,None];beta/=length
    Q,upper,piv=qr(N.T,mode='economic',pivoting=True)
    rank=int(np.sum(abs(np.diag(upper))>1e-10));assert rank==len(rows)
    C=N@Q;u0=Q.T@c
    assert min(C@u0-beta)>=-5e-12/min(length)
    solved=minimize(lambda u:.5*float(u@u),u0,jac=lambda u:u,method='SLSQP',
        constraints=[dict(type='ineq',fun=lambda u:C@u-beta,jac=lambda u:C)],
        options=dict(ftol=1e-12,maxiter=200,disp=False))
    assert solved.success,str(solved.message)
    # The native retained active pool has independent rows. A nonnegative
    # multiplier equality solution verifies its finite-pool optimum independently.
    mp.mp.dps=80
    normals=[[mp.mpf(float(x)) for x in row] for row in N]
    gram=mp.matrix([[sum(x*y for x,y in zip(a,b)) for b in normals] for a in normals])
    rhs=mp.matrix([mp.mpf(float(x)) for x in beta])
    multipliers=mp.lu_solve(gram,rhs)
    assert min(multipliers)>=0
    z=np.array([float(sum(multipliers[i]*normals[i][j] for i in range(len(rows)))) for j in range(n)])
    lam=np.array([float(x) for x in multipliers])
    objective=.5*float(z@z);lower=float(beta@lam-.5*np.linalg.norm(N.T@lam)**2)
    assert abs(objective-.5*float(solved.x@solved.x))<=2e-8
    assert abs(objective-.5*float(np.array(h.vector(record['z']))@np.array(h.vector(record['z']))))<=2e-8
    lu=lu_factor(J);load=B@z-rf
    d=lu_solve(lu,load);Jlong=J.astype(np.longdouble);extended=d.astype(np.longdouble)
    for _ in range(4):extended+=lu_solve(lu,np.array(load.astype(np.longdouble)-Jlong@extended,float)).astype(np.longdouble)
    d=np.array(extended,float)
    gaps=original.astype(np.longdouble)@d.astype(np.longdouble)-ell.astype(np.longdouble)
    transformed=(N.astype(np.longdouble)@z.astype(np.longdouble)-beta.astype(np.longdouble))*length.astype(np.longdouble)
    map_error=float(np.max(abs(gaps-transformed)))
    # Separate verification of physical J mapping, including the nonzero base C.
    mapped=Jlong@d.astype(np.longdouble)+rf.astype(np.longdouble)
    algebraic=float(np.max(abs(mapped[~free])))
    metric_obj=.5*sum(x*x for x in metric.norms(np.array(mapped,float),con)[:raw['dimension']+2])
    assert abs(metric_obj-objective)<=2e-8 and objective<=seed_obj+2e-9
    assert map_error<=5e-12 and max(-gaps)<=5e-12 and algebraic<=5e-12
    assert objective-lower<=2e-8
    state=np.array(h.vector(manifest['state']));effective=state.copy();delta=d.copy()
    for col in con:
        # K wall rows use the actual external lifting. Dependent epsilon rows
        # do not occur in the K polynomial and do not affect this calculation.
        effective[col]=state[col]-R[col];delta[col]=0.
    reserve=delta+.95*effective
    from decimal import Decimal
    extrema=[]
    for precision in (80,110):
        candidates=h.wall_rows(raw,reserve,precision)
        worst=min(candidates,key=lambda x:Decimal(x['minimum']['value']))
        extrema.append(dict(worst['minimum'],geometry=worst['geometry']))
    wall_min=Decimal(extrema[-1]['value'])/Decimal(float(metric.state_scale(raw['dimension'])))
    assert wall_min>=Decimal('-5e-12')
    assert abs(Decimal(extrema[0]['value'])-Decimal(extrema[1]['value']))<Decimal('1e-65')
    a.output.parent.mkdir(parents=True,exist_ok=True)
    destination=a.output.with_suffix('.direction.bin');np.asarray(d,dtype='<f8').tofile(destination)
    h.save(a.output,dict(status='finite_pool_reference_verified_pending_full_oracle',script=h.pin(__file__),
        native_record=dict(progress=h.pin(case/'case.r3-progress.jsonl'),sequence=record['sequence']),
        independent_solver='SciPy SLSQP / pivoted QR; 80-digit mpmath KKT polish; independent SciPy refined J solve',
        pool_size=len(rows),rank=rank,seed_coordinate_objective=seed_obj,seed_Phi=seed_obj*.125**2,
        recentered_objective_at_zero=seed_obj,scipy_success=bool(solved.success),scipy_iterations=int(solved.nit),
        scipy_objective=float(solved.fun),Psi=objective,Phi=objective*.125**2,
        nonnegative_dual_bound=lower,primal_dual_gap=objective-lower,minimum_multiplier=float(min(multipliers)),
        stationarity=float(np.linalg.norm(z-N.T@lam)),complementarity=float(np.max(abs(lam*(N@z-beta)))),
        active_original_violation=float(max(0.,max(-gaps))),transformation_gap_error=map_error,
        full_algebraic_coefficient_error=algebraic,independent_metric_objective=metric_obj,
        independent_continuous_wall_reserve=extrema[-1],normalized_wall_reserve_minimum=str(wall_min),
        direction=h.pin(destination),seconds=time.monotonic()-begin,physical_steps=0,
        qualification='Finite pool is a relaxation. Require the complete native plus extended original oracle and recomputed continuous-wall extrema.'))

if __name__=='__main__':main()
