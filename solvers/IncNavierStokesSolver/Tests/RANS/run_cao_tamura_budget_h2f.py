#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H2f bounded native setup, accessibility and certification orchestration."""
import argparse,os,time,resource
from pathlib import Path
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1')
import numpy as np
from cao_tamura_initializer_context import NativeContext
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic,freeze
from diagnose_cls_accepted_covariance import run_probe
from run_cls_epawf_picard import sha,state_hash,require
from search_cao_tamura_budget_h2d import tensors
from run_cao_tamura_budget_h2e import certify

def select(a):
    a.output.mkdir(exist_ok=False);c=NativeContext(a.context);pre=load(a.preflight);start=load(a.start);T,Ts,_=c.transforms();x=np.array(start['conditioned_velocity_coefficients']);raw=c.block('velocity_basis');d=np.einsum('nbc,b->nc',raw,x@T.T);g=c.gradient()
    for column,(i,j) in enumerate(((0,0),(0,2),(2,0),(2,2))):g[:,i,j]+=d[:,column]
    k=c.block('k');e=c.block('epsilon');inner=c.block('inner').astype(bool);cov=c.block('covariance').astype(bool);R=tensors(g,k,e,inner,c.block('inner_nu_t'),c.block('physical'),pre['u_tau'],pre['kappa_CLS']);ev,V=np.linalg.eigh(R)
    require(ev[inner&cov,0].min()>=1e-4*pre['heights'][-1]['K0'],'Stage A starting inner margin failed')
    outer=np.flatnonzero(~inner&cov);phys=c.block('physical');chart=c.block('chart');chosen=[]
    for idx in outer[np.argsort(ev[outer,0])]:
        if all(np.linalg.norm(phys[idx]-phys[j])>=.025 for j in chosen):chosen.append(int(idx))
        if len(chosen)==6:break
    for mask in ((~inner)&cov&(chart[:,2]>.799999),(~inner)&cov&(abs(chart[:,2]-.0065)<1e-12)):
        idx=np.flatnonzero(mask)
        if len(idx):j=int(idx[np.argmin(ev[idx,0])]);
        else:continue
        if j not in chosen:chosen.append(j)
    idx=np.flatnonzero(inner&cov);chosen.append(int(idx[np.argmin(ev[idx,0])]))
    target=np.array([.23549580558569316,0,.0065]);idx=np.flatnonzero(inner&cov&(abs(chart[:,2]-.0065)<1e-12));chosen.append(int(idx[np.argmin(np.linalg.norm(chart[idx]-target,axis=1))]))
    member=c.array('membership');_,first=np.unique(member[:,1],return_index=True);orig={int(member[j,1]):member[j].tolist() for j in first};w=[]
    K0=pre['heights'][-1]['K0']
    for j in dict.fromkeys(chosen):
        for eig in np.flatnonzero(ev[j]-ev[j,0]<=1e-7*K0):
            require(len(w)<16,'near-degenerate witness RHS exceeds predeclared budget')
            m=orig[j];w.append(dict(index=j,chart=chart[j].tolist(),physical=phys[j].tolist(),geometry=m[3],element=m[4],radial_side=m[5],tangential_side=m[6],dependent_epsilon=bool(m[8]),eigenvector=V[j,:,eig].tolist(),eigenvalue=float(ev[j,eig]),inner=bool(inner[j])))
    basis=np.array(c.array('basis')[:c.counts['velocity']],dtype='<f8');basis.tofile(a.output/'velocity.f64');atomic(a.output/'selection.json',dict(version='H2fWitnessSelectionV1',witnesses=w,start_state=start['state'],start_sha256=start['state_sha256'],base_sha256=c.meta['base_sha256'],basis_file=str(a.output/'velocity.f64'),basis_sha256=sha(a.output/'velocity.f64'),context_sha256=sha(a.context),start_inner_minimum=float(ev[inner&cov,0].min()),start_outer_minimum=float(ev[~inner&cov,0].min())))

def native(a):
    a.output.mkdir(exist_ok=False);ctx=load(a.context);selection=load(a.selection);exe,image=freeze(a.probe,a.output/'runtime');folder=a.output/'arrays';folder.mkdir()
    p=dict(operation='h2f_native',inner_stress_model='Linear',xml_sha256=ctx['xml_sha256'],base_sha256=ctx['context']['base_sha256'],basis_file=selection['basis_file'],basis_sha256=selection['basis_sha256'],velocity_direction_count=32,witnesses=selection['witnesses'],start_state=selection['start_state'],context_directory=str(folder))
    d,receipt=run_probe(exe,Path(ctx['xml']),p,a.output/'native');d=load(Path(receipt['output']));d.update(array_directory=str(folder),runtime=image,receipt=receipt,selection=str(a.selection),peak_RSS_KiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss);atomic(a.output/'native-data.json',d)

def context(a):
    a.output.mkdir(exist_ok=False);old=load(a.context);exe,image=freeze(a.probe,a.output/'runtime');folder=a.output/'arrays';folder.mkdir();basis=np.fromfile(a.basis,dtype='<f8').reshape(-1,len(NativeContext(a.context).array('base')[0]));count=len(basis)
    p=dict(operation='h2f_context',inner_stress_model='Linear',xml_sha256=old['xml_sha256'],base_sha256=old['context']['base_sha256'],basis_file=str(a.basis),basis_sha256=sha(a.basis),velocity_direction_count=count,context_directory=str(folder))
    d,receipt=run_probe(exe,Path(old['xml']),p,a.output/'native');d=load(Path(receipt['output']));d.update(xml=old['xml'],xml_sha256=old['xml_sha256'],array_directory=str(folder),runtime=image,receipt=receipt,preflight=str(a.preflight),peak_RSS_KiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss);atomic(a.output/'context.json',d)

def retained_slope(a):
    a.output.mkdir(exist_ok=False);c=NativeContext(a.context);doc=load(a.selection);start=load(a.start);pre=load(a.preflight);T,_,_=c.transforms();ids=np.flatnonzero((c.block('covariance')>0)&(c.block('inner')==0));q=int(ids[np.argmin(np.linalg.norm(c.block('chart')[ids]-np.array([-.022624687254059063,0,.13591982246103707]),axis=1))]);g=c.gradient()[q:q+1];d=c.block('velocity_basis')[q].T@(T@np.array(start['conditioned_velocity_coefficients']))
    for col,(i,j) in enumerate(((0,0),(0,2),(2,0),(2,2))):g[0,i,j]+=d[col]
    R=tensors(g,c.block('k')[q:q+1],c.block('epsilon')[q:q+1],np.array([False]),np.array([0.]),c.block('physical')[q:q+1],pre['u_tau'],pre['kappa_CLS']);ev,v=np.linalg.eigh(R[0]);m=next(row for row in c.array('membership') if row[1]==q);doc['witnesses']=[dict(index=q,chart=c.block('chart')[q].tolist(),physical=c.block('physical')[q].tolist(),geometry=int(m[3]),element=int(m[4]),radial_side=int(m[5]),tangential_side=int(m[6]),eigenvector=v[:,0].tolist(),eigenvalue=float(ev[0]),inner=False,reason='retained H2e Stage B slope witness reevaluated on Stage A')];atomic(a.output/'selection.json',doc)

def combine(a):
    import copy
    a.output.mkdir(exist_ok=False);first=load(a.native);second=load(a.extra);out=copy.deepcopy(first);out['array_directory']=str(a.output);out['component_receipts']=[first['receipt'],second['receipt']]
    for name,desc in first['context']['files'].items():
        p=Path(first['array_directory'])/desc['file']
        if name in ('gradients','global_responses'):
            other=second['context']['files'][name];values=np.vstack([np.fromfile(p,dtype='<f8').reshape(desc['shape']),np.fromfile(Path(second['array_directory'])/other['file'],dtype='<f8').reshape(other['shape'])]);values.tofile(a.output/desc['file']);out['context']['files'][name].update(shape=list(values.shape),sha256=sha(a.output/desc['file']))
        else:(a.output/desc['file']).symlink_to(p)
    out['context']['witnesses']+=second['context']['witnesses'];out['context']['response_RHS_total']=len(out['context']['witnesses']);atomic(a.output/'native-data.json',out)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('phase',choices=['select','native','context','certify','retained-slope','combine'])
    for name in ('context','preflight','start','selection','probe','basis','candidate','native','extra'):p.add_argument('--'+name,type=Path)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    {'select':select,'native':native,'context':context,'certify':certify,'retained-slope':retained_slope,'combine':combine}[a.phase](a)
if __name__=='__main__':main()
