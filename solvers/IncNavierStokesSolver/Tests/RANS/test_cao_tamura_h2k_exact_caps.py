#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Exact-cap candidate construction and checked obstruction on small fixtures."""
import argparse
from pathlib import Path
import numpy as np
from test_cao_tamura_budget_h2i import fixture
from cao_tamura_h2i_convex import JointProblem
from cao_tamura_h2k_exact_caps import exact_caps_phase
from cao_tamura_h2k_lp import equilibrated_lp
from cao_tamura_h2h_scalar import Budget
from continue_cls_epawf_pimple import atomic


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();tests={}
    P=fixture(repair=True);P.geometry=np.arange(P.n);v=np.array([.05,.3,-.4]);P.div=[np.outer(v,v)]*3
    r=exact_caps_phase(P,2e-4*P.K0,[.501,1.999],Budget(30,80,30),lambda event:None,lp_solver=equilibrated_lp)
    assert r['termination']=='seen_feasible',r;s=P.scan(np.array(r['best_coefficients']),2e-4*P.K0,Budget());assert s['hard_violation']<1e-9 and min(s['inner_min'],s['outer_min'])>=2e-4*P.K0-1e-12;tests['signed_native_repair_true_quadratics_and_all_scalar_rows']=True
    P=fixture(inner_bad=True);P.geometry=np.arange(P.n);r=exact_caps_phase(P,1e-4*P.K0,[.5,2],Budget(30,80,30),lambda event:None,lp_solver=equilibrated_lp)
    assert r['termination']=='verified_relaxation_obstruction' and r['proof']['receipt']['certificate']['verified'];tests['numerical_status_alone_does_not_claim_obstruction']=True
    P=fixture(repair=True);P.geometry=np.arange(P.n);r=exact_caps_phase(P,1e-4*P.K0,[.5,2],Budget(0,20,30),lambda event:None,lp_solver=equilibrated_lp);assert r['termination'].startswith('budget_exhausted');tests['work_limit_remains_inconclusive']=True
    count=401;nv=4;g=np.zeros((count,4));g[:,0]=np.linspace(.5,0,count);g[0,0]=100;k=np.full(count,.02);e=np.full(count,.01);inner=np.zeros(count,bool);inner[-1]=True;nt=np.zeros(count);nt[-1]=1e-4;VB=np.zeros((count,nv,6));VB[0,0,0]=1000
    for j in range(1,nv):VB[396+j,j-1,3]=-1000;VB[396+j,j,3]=2000
    E=np.zeros((count,8));E[:8]=np.eye(8)*.01;M=np.diag([0.]+[1.]*nv);D=np.zeros_like(M);P=JointProblem(g,np.zeros((count,2)),np.zeros((count,2)),k,e,k,e,nt,inner,np.ones(count,bool),VB,E,M,[D,D,D],.02);P.geometry=np.where(inner,1,0)
    r=exact_caps_phase(P,2e-4*P.K0,[.501,1.999],Budget(40,80,30),lambda event:None,lp_solver=equilibrated_lp);assert r['termination']=='seen_feasible',r;assert len(r['rounds'])>=4;tests['zero_relaxation_slack_does_not_stop_new_global_violations']=True
    atomic(a.output,dict(version='H2kExactCapsTestsV1',passed=True,tests=tests))


if __name__=='__main__':main()
