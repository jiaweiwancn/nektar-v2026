#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Exactly one authorized H2e-A velocity / H2e-B permitted-scalar combination."""
import argparse,os
from pathlib import Path
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1')
import numpy as np
from cao_tamura_initializer_context import NativeContext
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import sha,state_hash,require
from search_cao_tamura_budget_h2d import tensors

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('context','stage-a','stage-b','preflight','output'):p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();a.output.mkdir(exist_ok=False);c=NativeContext(a.context);A=load(a.stage_a);B=load(a.stage_b);pre=load(a.preflight);base=np.array(c.array('base')[0]);state=base.copy();vr=np.array(c.meta['rows']['velocity'],int);sr=np.array(c.meta['rows']['k']+c.meta['rows']['epsilon'],int);state[vr]=np.array(A['state'])[vr];state[sr]=np.array(B['state'])[sr];fixed=np.ones(len(base),bool);fixed[np.r_[vr,sr]]=False;require(state[fixed].tobytes()==base[fixed].tobytes(),'hybrid changed fixed bits')
    T,Ts,_=c.transforms();x=np.zeros(T.shape[1]);x[:len(A['conditioned_velocity_coefficients'])]=A['conditioned_velocity_coefficients'];s=np.array(B['conditioned_scalar_coefficients']);raw=c.block('velocity_basis');d=np.einsum('nbc,b->nc',raw,x@T.T);g=c.gradient()
    for col,(i,j) in enumerate(((0,0),(0,2),(2,0),(2,2))):g[:,i,j]+=d[:,col]
    ds=np.column_stack([c.block('k_basis'),c.block('epsilon_basis')])@Ts@s;k=c.block('k').copy();e=c.block('epsilon').copy();nk=c.counts['k'];coeff=Ts@s;k+=c.block('k_basis')@coeff[:nk];e+=c.block('epsilon_basis')@coeff[nk:];inner=c.block('inner').astype(bool);cov=c.block('covariance').astype(bool);R=tensors(g,k,e,inner,c.block('inner_nu_t'),c.block('physical'),pre['u_tau'],pre['kappa_CLS']);ev=np.linalg.eigvalsh(R)[:,0]/pre['heights'][-1]['K0'];Q=np.array(c.array('divergence'));z=np.r_[1,T@x];limit=np.sqrt(Q[0,0])*(1+1e-6)+1e-9
    atomic(a.output/'candidate.json',dict(version='H2fSingleHybridV1',stage='B',state=state.tolist(),state_sha256=state_hash(state),base_sha256=state_hash(base),conditioned_velocity_coefficients=x.tolist(),conditioned_scalar_coefficients=s.tolist(),min_lambda_over_K0=float(ev[cov].min()),inner_min_lambda_over_K0=float(ev[inner&cov].min()),divergence_limit=limit,strong_divergence_L2=float(np.sqrt(max(0,z@Q@z))),source_A_sha256=sha(a.stage_a),source_B_sha256=sha(a.stage_b),context_sha256=sha(a.context)))
if __name__=='__main__':main()
