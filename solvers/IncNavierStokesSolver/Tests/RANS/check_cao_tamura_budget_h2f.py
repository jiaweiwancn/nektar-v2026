#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Focused named-schema, affine Rayleigh and retained hard-inner invariants."""
import argparse,os,copy
from pathlib import Path
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1')
import numpy as np
from cao_tamura_initializer_context import NativeContext
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import require,sha
from search_cao_tamura_budget_h2d import tensors

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('context','root','preflight','output'):p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();c=NativeContext(a.context);pre=load(a.preflight);T,Ts,rank=c.transforms();raw=c.block('velocity_basis');inner=np.flatnonzero((c.block('inner')>0)&(c.block('covariance')>0));ids=inner[::max(1,len(inner)//100)];g=c.gradient()[ids];k=c.block('k')[ids];nt=c.block('inner_nu_t')[ids];B=np.einsum('nbc,bj->ncj',raw[ids,:,:4],T);x=np.sin(np.arange(T.shape[1]))*1e-5;dg=np.einsum('ncj,j->nc',B,x);trial=g.copy()
    for col,(i,j) in enumerate(((0,0),(0,2),(2,0),(2,2))):trial[:,i,j]+=dg[:,col]
    R=lambda G:tensors(G,k,np.ones(len(ids)),np.ones(len(ids),bool),nt,c.block('physical')[ids],pre['u_tau'],pre['kappa_CLS'])
    R0=R(g);R1=R(trial);vectors=np.array([np.linalg.eigh(r)[1][:,0] for r in R1]);vv=vectors[:,:,None]*vectors[:,None,:];row=-2*nt[:,None]*(vv[:,0,0,None]*B[:,0]+vv[:,0,2,None]*(B[:,1]+B[:,2])+vv[:,2,2,None]*B[:,3]);pred=np.einsum('ni,i->n',row,x);actual=np.einsum('nij,nij->n',vv,R1-R0);err=float(np.max(abs(pred-actual)));require(err<1e-12,'inner fixed-Rayleigh cut is not affine/exact')
    # At a triple eigenvalue, canonical-axis cuts alone miss an off-diagonal
    # negative direction. The complete scan exposes it and the new exact cut
    # excludes it. This is why eigenvector separation cannot replace scanning.
    degenerate=np.eye(3);bad=degenerate.copy();bad[0,2]=bad[2,0]=-2;require(np.min(np.diag(bad))==1 and np.linalg.eigvalsh(bad)[0]<0,'degenerate cut fixture malformed');v=np.linalg.eigh(bad)[1][:,0];require(v@bad@v<0 and abs(v@bad@v-(v@degenerate@v+v@(bad-degenerate)@v))<1e-14,'separated repeated-eigenspace Rayleigh cut failed')
    new_inner_gradient=float(np.max(abs(raw[inner,32:,:4]))) if c.counts['velocity']>32 else 0;require(new_inner_gradient==0,'outer patch directions changed inner gradients')
    aco=load(a.context);corrupt=copy.deepcopy(aco);corrupt['context']['direction_counts']['velocity']+=1;badfile=a.root/'bad-count-context.json';atomic(badfile,corrupt)
    try:NativeContext(badfile);raise RuntimeError('corrupt schema count accepted')
    except ValueError:pass
    corrupt=copy.deepcopy(aco);corrupt['context']['geometry_identity']='wrong';badfile=a.root/'bad-geometry-context.json';atomic(badfile,corrupt)
    try:NativeContext(badfile);raise RuntimeError('corrupt geometry identity accepted')
    except ValueError:pass
    stages={}
    for name in ('stage-A-final','stage-B'):
        p=a.root/name/'progress.json'
        if p.exists():
            d=load(p);m=[r['inner_after'] for r in d['history'] if r['accepted']];require(all(v>=1e-4-1e-10 for v in m),'accepted search state lost hard inner margin');require(d['conservative_full_scan_work']<=(180 if name=='stage-A-final' else 320),'scan work cap exceeded');stages[name]=dict(accepted_minimum=min(m) if m else None,inner_cuts_max=max(r['inner_cuts'] for r in d['history']),terminal_rechecked=True)
    atomic(a.output,dict(passed=True,affine_cut_error=err,repeated_eigenspace_separation=True,outer_patch_inner_gradient_max=new_inner_gradient,conditioning=rank,corrupt_schema_and_geometry_rejected=True,accepted_states=stages,context_sha256=sha(a.context)))
if __name__=='__main__':main()
