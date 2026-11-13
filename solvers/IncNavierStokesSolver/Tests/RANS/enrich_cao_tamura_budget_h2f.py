#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded patch-constrained native Riesz responses, with no global projector."""
import argparse,os,time,resource
from pathlib import Path
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1')
import numpy as np
from scipy.linalg import cho_factor,cho_solve,qr
from cao_tamura_initializer_context import NativeContext
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import sha,require

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('native','context','output'):p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();a.output.mkdir(exist_ok=False);began=time.monotonic();doc=load(a.native);f=doc['context'];root=Path(doc['array_directory']);c=NativeContext(a.context)
    def arr(name):d=f['files'][name];require(sha(root/d['file'])==d['sha256'],'native data checksum');return np.memmap(root/d['file'],dtype=d['dtype'],mode='r',shape=tuple(d['shape']))
    M=np.array(arr('mass_free'));Q=np.array(arr('continuity_Q'));grad=np.array(arr('gradients'));globalR=np.array(arr('global_responses'));bounds=arr('support_bounds');free=np.array(f['free_rows'],int);field_n=f['velocity_field_size'];old=np.array(c.array('basis')[:32]);oldfree=old[:,free];T,_,_=c.transforms();oldorth=T.T@oldfree;w=f['witnesses'];outer=np.array([not x['inner'] for x in w]);global_norm=np.einsum('ni,ij,nj->n',globalR,M,globalR)
    patches=[];selected=None
    for label,box in [('primary',(-.1,.1,.0065,.4)),('enlarged',(-.2,.6,.0065,.8))]:
        support=bounds[free%field_n];mask=(support[:,0]>=box[0]-1e-13)&(support[:,1]<=box[1]+1e-13)&(support[:,2]>=box[2]-1e-13)&(support[:,3]<=box[3]+1e-13);ids=np.flatnonzero(mask);Mp=M[np.ix_(ids,ids)];mass=cho_factor(Mp,lower=True,check_finite=False)
        # Include every affected global Q row. Add old-span physical-mass
        # orthogonality here, so the residual directions stay patch-supported.
        H=np.vstack([Q[:,ids],oldorth@M[:,ids]]);norm=np.linalg.norm(H,axis=1);keep=norm>1e-14;Hn=H[keep]/norm[keep,None]
        _,R,piv=qr(Hn.T,pivoting=True,mode='economic');diag=abs(np.diag(R));rank=int(np.sum(diag>diag[0]*1e-10)) if len(diag) else 0;ind=Hn[piv[:rank]];iq=cho_solve(mass,ind.T);schur=cho_factor(ind@iq,lower=True,check_finite=False)
        responses=[];fractions=[]
        for g in grad:
            v=cho_solve(mass,g[ids]);v-=iq@cho_solve(schur,ind@v);r=np.zeros(len(free));r[ids]=v;responses.append(r)
        responses=np.array(responses);norms=np.einsum('ni,ij,nj->n',responses,M,responses);frac=np.divide(norms,global_norm,out=np.zeros_like(norms),where=global_norm>1e-24)
        scale=np.maximum(1,np.sqrt(np.maximum(0,norms)));require(np.max(abs(Q@responses.T)/scale)<1e-8,'patch lost a global weak-continuity row');require(np.max(abs(oldorth@M@responses.T)/scale)<1e-8,'patch lost old-span orthogonality')
        record=dict(patch=label,bounds=list(box),coefficient_count=len(ids),continuity_rows=len(Q),continuity_rank=int(np.linalg.matrix_rank(Q[:,ids],tol=1e-10)),combined_rank=rank,nullspace_dimension=len(ids)-rank,response_fraction=frac.tolist(),Q_max=float(np.max(abs(Q@responses.T))),Q_mass_unit_max=float(np.max(abs(Q@responses.T)/scale)),old_overlap_max=float(np.max(abs(oldorth@M@responses.T))))
        patches.append(record);selected=(label,ids,responses,norms,record)
        if np.sum(frac[outer]<.05)<math_ceil_half(int(outer.sum())) and np.max(norms[outer])>1e-24:break
    label,ids,responses,norms,record=selected;new=[];decisions=[]
    # Target outer covariance first. Inner responses are diagnostics only; a
    # patch excluding Gamma_A has exactly zero inner field and gradient.
    for j in np.flatnonzero(outer):
        v=responses[j].copy();original=math_sqrt(v@M@v)
        for _ in range(2):
            for z in new:v-=z*(z@M@v)
        norm=math_sqrt(v@M@v)
        if original<=1e-12 or norm<=1e-10*max(original,1e-30):decisions.append(dict(witness=int(j),retained=False,residual_norm=norm));continue
        v/=norm;new.append(v);decisions.append(dict(witness=int(j),retained=True,residual_norm=norm,original_norm=original))
        if len(new)==16:break
    require(new,'no useful independent patch freedom; stop enrichment rather than fill quota')
    Z=np.zeros((len(new),old.shape[1]));Z[:,free]=new;all_basis=np.vstack([old,Z]);all_basis.astype('<f8').tofile(a.output/'velocity.f64')
    gram=np.array(new)@M@np.array(new).T;cross=oldorth@M@np.array(new).T;Qerr=Q@np.array(new).T
    require(np.max(abs(gram-np.eye(len(new))))<1e-8 and np.max(abs(cross))<1e-8 and np.max(abs(Qerr))<1e-8,'enrichment constraints failed')
    atomic(a.output/'enrichment.json',dict(version='H2fPatchRieszV1',native_data_sha256=sha(a.native),old_span=32,new_rank=len(new),total_velocity_directions=len(all_basis),patches=patches,selected_patch=label,decisions=decisions,basis_file=str(a.output/'velocity.f64'),basis_sha256=sha(a.output/'velocity.f64'),mass_Gram_error=float(np.max(abs(gram-np.eye(len(new))))),old_span_overlap=float(np.max(abs(cross))),Q_error=float(np.max(abs(Qerr))),global_accessibility=w,enriched_eta_squared=[float(q['eta_squared']+sum((z@M@rr)**2 for z in new)/den) if den>1e-24 else None for q,rr,den in zip(w,globalR,global_norm)],seconds=time.monotonic()-began,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,interpretation='local constrained gradient freedom, not a covariance feasibility certificate'))
def math_sqrt(x):return float(np.sqrt(max(0,x)))
def math_ceil_half(n):return (n+1)//2
if __name__=='__main__':main()
