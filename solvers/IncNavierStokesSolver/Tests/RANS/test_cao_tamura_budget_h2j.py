#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Finite divergence identities, support conditioning, complete V2 producer."""
import argparse,copy,tempfile,time
from pathlib import Path
import numpy as np
from scipy import linalg
from cao_tamura_h2d_state import load
from cao_tamura_startup_receipt import complete_receipt,validate_receipt,divergence_gate
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import sha,require

def rejected(f):
    try:f()
    except (ValueError,KeyError):return
    raise AssertionError('incompatible V2 receipt accepted')

def main():
    p=argparse.ArgumentParser();p.add_argument('--receipt',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();start=time.monotonic();tests={}
    receipt=load(a.receipt);validate_receipt(receipt);require(receipt['base_sha256']!='038fc076cebbea41bcb7cc349147eacc2b938a9051827bbb386baf7c3cb3716c','flat fixture falsely uses hill base');tests['real_flat_complete_receipt_positive']=True
    from run_cao_tamura_budget_h2f_startup import StartupRun
    progress=load(a.receipt.parent/'import/progress.json');run=StartupRun.__new__(StartupRun);run.out=progress;report=load(a.receipt.parent/'import/zero/case.evolution.json');run.validate_result('zero',report,0,2e-6);tests['actual_native_empty_success_report_checkpoint_validated']=True
    fake=dict(native_divergence_GL16=dict(order=16,regional_L2=[1.5e-8,0.],L2=1.5e-8),audit=dict(strong_continuity=[1.5e-8,0.]));assert not divergence_gate(fake,1e-8)['passed'];tests['native_GL16_does_not_use_extra_grid_allowance']=True
    for key in ('version','base_sha256','probe_output_file','probe_output_sha256','state_sha256','model_specification'):
        x=copy.deepcopy(receipt);x.pop(key);rejected(lambda:validate_receipt(x))
    for key in ('extra_covariance','strong_divergence','native_GL16'):
        x=copy.deepcopy(receipt);x['gates'].pop(key);rejected(lambda:validate_receipt(x));x=copy.deepcopy(receipt);x['gates'][key]=False;rejected(lambda:validate_receipt(x))
    for key in ('state_sha256','base_sha256','probe_output_sha256'):
        x=copy.deepcopy(receipt);x[key]='0'*64;rejected(lambda:validate_receipt(x))
    tests['missing_gate_report_hash_and_state_base_rejected']=True
    raw=load(Path(receipt['probe_output_file']))
    with tempfile.TemporaryDirectory() as temp:
        path=Path(temp)/'report.json'
        for failed in ('covariance','divergence'):
            x=copy.deepcopy(raw)
            if failed=='covariance':x['audit']['accepted_awf_covariance'][0]=-.001
            else:x['native_divergence_GL16'].update(regional_L2=[1.,1.],L2=np.sqrt(2.))
            atomic(path,x);r=complete_receipt(path,receipt['model_specification'],receipt['divergence_limit'],receipt['validation_probe'],receipt['solver_runtime'],receipt['extra_sampling']);assert not r['passed'];rejected(lambda:validate_receipt(r))
    tests['actual_failed_covariance_and_divergence_cannot_issue']=True
    # A tangent at an active sphere permits an orthogonal direction even
    # though every nonzero finite step violates the cap.
    u=np.array([1.,0.]);d=np.array([0.,1.]);assert u@d==0 and np.linalg.norm(u+.1*d)>1;tests['first_order_divergence_tangent_is_insufficient']=True
    from analyse_cao_tamura_budget_h2i import tensor_check
    G=np.diag([.2,0.,-.1]);k=.02;e=.01;nu=.09*k*k/e;R=2*k/3*np.eye(3)-nu*(G+G.T);ev,V=np.linalg.eigh(R);w=dict(observation=dict(k=k,epsilon=e,inner=False,gradient=G.tolist()),tensor=R.tolist(),critical_eigenvector=V[:,0].tolist(),eigenvalues=ev.tolist());v=tensor_check(w);assert abs(v['production']+np.sum(R*G))<1e-15 and abs(v['isotropic_work'])>1e-4;tests['independent_production_keeps_isotropic_divergence_work']=True
    rng=np.random.default_rng(77);B=rng.normal(size=(17,9));weights=np.exp(rng.normal(size=17));u=rng.normal(size=9);d=rng.normal(size=9);Q=B.T@(weights[:,None]*B)
    for alpha in (-.3,.02,.4):assert abs((u+alpha*d)@Q@(u+alpha*d)-(u@Q@u+2*alpha*u@Q@d+alpha*alpha*d@Q@d))<1e-11
    tests['full_weighted_divergence_constant_cross_and_quadratic_identity']=True
    M=np.diag([2.,3.,4.]);W=linalg.solve_triangular(np.linalg.cholesky(M).T,np.eye(3),lower=False);C=np.array([[1.,0.,0.]])@W;_,sv,V=linalg.svd(C,full_matrices=True);N=W@V[1:].T;assert np.max(abs(N.T@M@N-np.eye(2)))<1e-14 and np.max(abs(C@(np.linalg.inv(W)@N)))<1e-14
    # Blind subtraction destroys a qualified physical direction; retaining its
    # exact old-span companion reconstructs the qualified response.
    d=np.array([0.,1.]);z=np.array([1.,1.])/np.sqrt(2);orth=d-z*(z@d);assert abs(orth[0])>.4 and np.max(abs(orth+z*(z@d)-d))<1e-15;tests['raw_response_or_exact_companion_required']=True
    A=np.array([[3.,.2],[.2,2.]]);Q=np.array([[1.,1.]]);rhs=np.array([1.,-.2]);x0=np.linalg.solve(A,rhs);AiQ=np.linalg.solve(A,Q.T);x=x0-AiQ@np.linalg.solve(Q@AiQ,Q@x0);kkt=np.block([[A,Q.T],[Q,np.zeros((1,1))]]);direct=np.linalg.solve(kkt,np.r_[rhs,0])[:2];assert np.max(abs(x-direct))<1e-14 and abs(Q@x)<1e-14;tests['energy_KKT_elimination_matches_full_system']=True
    atomic(a.output,dict(version='H2jFocusedChecksV1',passed=True,tests=tests,seconds=time.monotonic()-start,scope='mathematical fixtures and actual positive flat receipt; not a hill certificate'))
if __name__=='__main__':main()
