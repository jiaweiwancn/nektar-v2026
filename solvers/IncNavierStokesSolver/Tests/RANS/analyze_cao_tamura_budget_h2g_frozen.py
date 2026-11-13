#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent six-law full-inventory and physical-GL20 frozen-field reduction."""
import argparse,os,time,resource,subprocess
from pathlib import Path
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
import numpy as np
from cao_tamura_h2d_state import load
from cao_tamura_h2g_constitutive import LAWS,parts
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import require,sha

def decode(a):return a[:,6:15].reshape(-1,3,3),a[:,15],a[:,16],a[:,18].astype(bool),a[:,17],a[:,3:6]
def tangents(g,k,e,p,us,kap,law):
    g=np.array(g);p=np.array(p);S=g+g.T;O=g-g.T;qs=k/e*np.sqrt(np.sum(S*S)/2);qo=k/e*np.sqrt(np.sum(O*O)/2);tie=bool(abs(qs-qo)<1e-7*max(1,qs,qo))
    output=[]
    for branch in (('strain','rotation') if tie and law!='calibrated_standard' else ('physical',)):
        def stress(G,K=k,E=e):return parts(G[None],np.array([K]),np.array([E]),np.array([False]),np.array([0.]),p[None],us,kap,law,branch)['tensor'][0]
        T=np.empty((3,3,3,3));err=0
        for a in range(3):
            for b in range(3):
                h=2e-6*max(1,abs(g[a,b]));d=np.zeros((3,3));d[a,b]=h;full=(stress(g+d)-stress(g-d))/(2*h);half=(stress(g+d/2)-stress(g-d/2))/h;T[:,:,a,b]=half;err=max(err,float(np.max(abs(full-half))))
        samples=[]
        for theta in np.linspace(0,np.pi,65)[:-1]:
            n=np.array([np.cos(theta),0,np.sin(theta)]);A=np.array([[-np.sin(theta),0,np.cos(theta)],[0,1,0]]);Q=np.empty((2,2))
            for j,a in enumerate(A):
                dR=np.einsum('ijab,ab->ij',T,np.outer(a,n));Q[:,j]=-A@dR@n
            samples.append(float(np.linalg.eigvalsh((Q+Q.T)/2)[0]))
        derivatives={}
        for name,value in [('k',k),('epsilon',e)]:
            h=2e-6*value
            if name=='k':d=(stress(g,k+h,e)-stress(g,k-h,e))/(2*h)
            else:d=(stress(g,k,e+h)-stress(g,k,e-h))/(2*h)
            derivatives[name]=d.tolist()
        output.append(dict(branch=branch,step_halving_absolute=err,minimum_symmetric_incompressible_symbol_scan=min(samples),maximum_min_symbol_scan=max(samples),velocity_tangent=T.tolist(),scalar_derivatives=derivatives,scope='64 x-z normal directions; both transverse velocity directions; fixed scalars; symmetric principal symbol only, not a full ellipticity proof'))
    return dict(q_strain=qs,q_rotation=qo,branch_tie=tie,derivatives=output)

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--preflight',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--native-probe',type=Path,required=True);a=p.parse_args();pre=load(a.preflight);us=pre['u_tau'];kap=pre['kappa'];start=time.monotonic();out=dict(version='H2gFrozenComparisonV1',states={},governing_solves=0,certificates=0,common_tuple_policy='no scalar, gradient, state, or reference substitution; six laws independently evaluated; inner law unchanged')
 for state in ('base','A','B'):
  d=load(a.root/(state+'-probe/result.json'));stateout=dict(state_sha256=d['state_sha256'],datasets={})
  for dataset in ('inventory','GL20'):
   meta=d[dataset];folder=a.root/(state+'-arrays');require(sha(folder/meta['file'])==meta['sha256'],'tuple hash');r=np.memmap(folder/meta['file'],dtype='<f8',mode='r',shape=tuple(meta['shape']));members=np.memmap(folder/meta['membership_file'],dtype='<i4',mode='r',shape=tuple(meta['membership_shape']));weights=np.array(r[:,19]);multiplicity=np.bincount(members[:,1],minlength=len(r));rows=[];allwitnesses=sorted(set(int(v['index']) for v in meta['witnesses']));g,k,e,inner,nt,physical=decode(r)
   for j,law in enumerate(LAWS):
    v=parts(g,k,e,inner,nt,physical,us,kap,law);ev,vec=np.linalg.eigh(v['tensor']);neg=ev[:,0]<0;worst=int(np.argmin(ev[:,0]));native=meta['witnesses'][j];error=float(np.max(abs(v['tensor'][native['index']]-np.array(native['tensor']))));mine=abs(float(ev[worst,0])-meta['native_minima'][j]);require(error<2e-11 and mine<2e-11,'native/independent tensor or eigenvalue mismatch')
    common=[]
    for idx in allwitnesses:
     vv=vec[idx,:,0];common.append(dict(index=idx,physical=physical[idx].tolist(),k=float(k[idx]),epsilon=float(e[idx]),gradient=g[idx].tolist(),eddy_viscosity=float(v['nu_t'][idx]),production=float(v['production'][idx]),Cmu=float(v['Cmu'][idx]),invariants=[float(v['q_strain'][idx]),float(v['q_rotation'][idx])],eigenvalues=ev[idx].tolist(),critical_vector=vv.tolist(),critical_contributions={n:float(vv@v[n][idx]@vv) for n in ('isotropic','linear','quadratic','cubic','regularization')}))
    negative_members=int(multiplicity[neg].sum());require(int(neg.sum())==meta['negative_unique'][j] and negative_members==meta['negative_memberships'][j],'negative unique/membership count disagreement')
    result=dict(law=law,minimum=float(ev[worst,0]),normalized_by_CLS_K0=float(ev[worst,0]/pre['models'][0]['K0']),witness=native,independent_critical_vector=vec[worst,:,0].tolist(),negative_unique=int(neg.sum()),negative_memberships=negative_members,negative_physical_volume=float(weights[neg].sum()) if dataset=='GL20' else None,negative_inner_volume=float(weights[neg&inner].sum()) if dataset=='GL20' else None,negative_outer_volume=float(weights[neg&~inner].sum()) if dataset=='GL20' else None,native_tensor_error=error,native_minimum_error=mine,common_witnesses=common)
    if dataset=='inventory' and not inner[worst]:result['frozen_tangent']=tangents(g[worst],k[worst],e[worst],physical[worst],us,kap,law)
    rows.append(result)
   stateout['datasets'][dataset]=dict(unique_count=len(r),membership_count=len(members),physical_volume=float(weights.sum()) if dataset=='GL20' else None,laws=rows)
  out['states'][state]=stateout;atomic(a.output,out)
 K=pre['models'][0]['K0'];z=.1;g=np.zeros((3,3));g[0,2]=us/(kap*z);eps=us**3/(kap*z);out['neutral_reference_tangents']={law:tangents(g,K,eps,[0,0,z],us,kap,law) for law in LAWS};tuples=[];references=[]
 for state in ('base','A','B'):
  for lawid,law in enumerate(LAWS):
   row=out['states'][state]['datasets']['inventory']['laws'][lawid]['witness'];G=np.array(row['gradient']);k=row['k'];e=row['epsilon'];phys=np.array(row['physical']);inn=row['inner']
   if inn:continue
   for branch in (0,1,2):
    for h in (0.,2e-6,-2e-6,1e-6,-1e-6):
     d=np.array([[.4,0,-.3],[0,.2,0],[.1,0,-.1]]);GG=G+h*d;kk=k*(1+.23*h);ee=e*(1-.17*h);tuples.append(dict(law=lawid,branch=branch,k=kk,epsilon=ee,physical_x=float(phys[0]),physical_z=float(phys[2]),flat_reference=False,gradient=GG.tolist()));references.append(parts(GG[None],np.array([kk]),np.array([ee]),np.array([False]),np.array([0.]),phys[None],us,kap,law,{0:'physical',1:'strain',2:'rotation'}[branch]))
 # Exact neutral branch-limit cases are distinct from physical max derivatives.
 for lawid,law in enumerate(LAWS):
  for branch in (1,2):
   for h in (0.,2e-6,-2e-6):
    GG=g+h*np.diag([.4,.1,-.2]);tuples.append(dict(law=lawid,branch=branch,k=K,epsilon=eps,physical_x=0.,physical_z=z,flat_reference=True,gradient=GG.tolist()));references.append(parts(GG[None],np.array([K]),np.array([eps]),np.array([False]),np.array([0.]),np.array([[0,0,z]]),us,kap,law,'strain' if branch==1 else 'rotation'))
 request=a.output.with_suffix('.constitutive-input.json');response=a.output.with_suffix('.constitutive-native.json');atomic(request,dict(u_tau=us,kappa=kap,tuples=tuples));subprocess.run([str(a.native_probe),str(request),str(response)],check=True);native=load(response);error=0.;derivative=0.
 for n,r in zip(native['tuples'],references):
  error=max(error,float(np.max(abs(np.array(n['tensor'])-r['tensor'][0]))),abs(n['production']-r['production'][0]))
 for start in range(0,len(native['tuples'])-36,5):
  nr=native['tuples'][start:start+5];rr=references[start:start+5]
  for plus,minus,step in ((1,2,2e-6),(3,4,1e-6)):
   dn=(np.array(nr[plus]['tensor'])-np.array(nr[minus]['tensor']))/(2*step);di=(rr[plus]['tensor'][0]-rr[minus]['tensor'][0])/(2*step);derivative=max(derivative,float(np.max(abs(dn-di))))
 require(error<2e-10 and derivative<2e-5,'native independent non-equilibrium/branch derivative discrepancy');out['constitutive_native_checks']=dict(tuple_count=len(tuples),tensor_and_production_absolute_error=error,combined_velocity_scalar_directional_error=derivative,physical_and_formal_branch_limits=True,probe_sha256=sha(a.native_probe));out.update(seconds=time.monotonic()-start,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,passed_native_independent=True);atomic(a.output,out)
if __name__=='__main__':main()
