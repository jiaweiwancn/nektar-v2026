# SPDX-License-Identifier: MIT
"""Equivalent LP equilibration with original-coordinate primal/dual checks."""
import time, warnings
import numpy as np
from scipy import sparse
from scipy.optimize import linprog, OptimizeWarning
from cao_tamura_h2h_scalar import dual_check


def equilibrated_lp(c,A,b,bounds,budget):
    # x = scale*y. Positive diagonal scaling preserves signed coefficients and
    # all finite physical bounds. No box, tolerance or acceptance is enlarged.
    maximum=np.asarray(abs(A).max(axis=0).toarray()).ravel()
    scale=1/np.maximum(maximum,1e-8);scale=np.minimum(scale,1e8)
    As=A@sparse.diags(scale);rowmax=np.asarray(abs(As).max(axis=1).toarray()).ravel();rows=1/np.maximum(1,np.maximum(rowmax,abs(b)));As=sparse.diags(rows)@As;bs=b*rows;cs=c*scale
    limits=[(lo/s,hi/s) for (lo,hi),s in zip(bounds,scale)];attempts=[];started=time.monotonic()
    for method in ('highs-ipm','highs-ds'):
        budget.round();tic=time.monotonic();options=dict(primal_feasibility_tolerance=1e-10,dual_feasibility_tolerance=1e-10,threads=1,time_limit=max(.01,budget.seconds-(time.monotonic()-budget.start)))
        if method=='highs-ipm':options['ipm_optimality_tolerance']=1e-12
        else:options['simplex_scale_strategy']=0
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore',category=OptimizeWarning,message='(Unrecognized options detected.*|Unknown solver options: threads)')
            result=linprog(cs,A_ub=As,b_ub=bs,bounds=limits,method=method,options=options)
        attempts.append(dict(method=method,status=int(result.status),message=result.message,seconds=time.monotonic()-tic,iterations=int(result.nit)))
        if result.status!=4:break
    receipt=dict(status=int(result.status),message=result.message,seconds=time.monotonic()-started,iterations=sum(x['iterations'] for x in attempts),rows=A.shape[0],columns=A.shape[1],attempts=attempts,scaling=dict(version='PositiveDiagonalLPScalingV1',column_min_max=[float(scale.min()),float(scale.max())],row_min_max=[float(rows.min()),float(rows.max())],physical_constraints_unchanged=True))
    if result.status==0:
        result.x=result.x*scale
        result.ineqlin.marginals=result.ineqlin.marginals*rows
        result.lower.marginals=result.lower.marginals/scale;result.upper.marginals=result.upper.marginals/scale
        result.ineqlin.residual=b-A@result.x
        result.lower.residual=result.x-np.array([x[0] for x in bounds]);result.upper.residual=np.array([x[1] for x in bounds])-result.x
        receipt['certificate']=dual_check(c,A,b,bounds,result)
    return result,receipt
