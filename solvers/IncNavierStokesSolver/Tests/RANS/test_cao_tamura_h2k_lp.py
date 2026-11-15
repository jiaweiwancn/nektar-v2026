#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Original-coordinate LP certificates after equivalent numerical scaling."""
import argparse
from types import SimpleNamespace
from unittest.mock import patch
from pathlib import Path
import numpy as np
from scipy import sparse
import cao_tamura_h2k_lp as solver
from cao_tamura_h2h_scalar import Budget
from continue_cls_epawf_pimple import atomic


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();tests={}
    A=sparse.csr_matrix([[1e4,-1.],[-1e4,-1.]]);b=np.array([-.4,-.2]);c=np.array([0.,1.]);bounds=[(-1.,1.),(0.,1.)]
    r,d=solver.equilibrated_lp(c,A,b,bounds,Budget());assert r.status==0 and abs(r.x[0]+1e-5)<1e-12 and abs(r.x[1]-.3)<1e-12
    assert d['certificate']['verified'] and d['certificate']['lower_bound']>.2999999 and np.max(A@r.x-b)<1e-9;tests['signed_variables_and_original_coordinate_primal_dual']=True
    zero=sparse.csr_matrix([[-1.,-1.],[1.,-1.]]);r,d=solver.equilibrated_lp(c,zero,np.array([-.2,.5]),bounds,Budget());assert r.status==0 and abs(r.x[1])<1e-12 and d['certificate']['verified'];tests['known_feasible_problem']=True
    original=solver.linprog;calls=[]
    def fail_once(*args,**kwargs):
        calls.append(kwargs['method'])
        return SimpleNamespace(status=4,message='injected numerical status',nit=0) if len(calls)==1 else original(*args,**kwargs)
    budget=Budget()
    with patch.object(solver,'linprog',side_effect=fail_once):r,d=solver.equilibrated_lp(c,A,b,bounds,budget)
    assert r.status==0 and budget.rounds==2 and calls==['highs-ipm','highs-ds'] and len(d['attempts'])==2;tests['bounded_fallback_counts_both_LP_attempts']=True
    atomic(a.output,dict(version='H2kEquivalentLPTestsV1',passed=True,tests=tests))


if __name__=='__main__':main()
