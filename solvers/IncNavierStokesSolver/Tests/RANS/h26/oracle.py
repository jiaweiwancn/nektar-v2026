#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent tiny direct assembly using Jacobi polynomials and analytic geometry."""
import argparse,json,os
from pathlib import Path
from run import THREADS
os.environ.update(THREADS)
import numpy as np
import scipy.linalg as la
from numpy.polynomial.legendre import leggauss
from fixtures import modified


def tensor(degree,xi,eta,dx,dz,slope,sign):
    bx,gx=modified(degree,xi);bz,gz=modified(degree,eta)
    b=(bz[:,:,None]*bx[:,None,:]).reshape(len(xi),-1)*sign
    gz=(gz[:,:,None]*bx[:,None,:]).reshape(len(xi),-1)*sign*2/dz
    gx=(bz[:,:,None]*gx[:,None,:]).reshape(len(xi),-1)*sign*2/dx-slope[:,None]*gz
    return b,gx,gz


def native_sparse(packet):
    m=np.zeros((packet['rows'],packet['cols']))
    for i,j,v in packet['entries']:m[i,j]=v
    return m


def assemble(data,curve,order=24):
    n,p=data['velocity_size'],data['pressure_size'];M=np.zeros((n,n));Q=np.zeros((p,p));D=np.zeros((p,2*n));K=np.zeros((2*n,2*n));Klap=np.zeros_like(K);C=np.zeros((2*n,p));mean=np.zeros(p)
    z,w=leggauss(order);xx,zz=np.meshgrid(z,z);xi=xx.ravel();eta=zz.ravel();ww=np.outer(w,w).ravel()
    for e in data['elements']:
        xl,xr,rl,rr=e['bounds'];dx=xr-xl;dz=rr-rl;x=(xl+xr)/2+xi*dx/2;r=(rl+rr)/2+eta*dz/2
        ground=curve*x*x*(1-x)*(1-x);slope=curve*(2*x-6*x*x+4*x*x*x);physical_z=r+ground
        b,gx,gz=tensor(3,xi,eta,dx,dz,slope,np.array(e['us']));q,_,_=tensor(2,xi,eta,dx,dz,slope,np.array(e['ps']));weight=ww*dx*dz/4;mu=.02*(1+x+physical_z)
        ids=np.array(e['u']);ps=np.array(e['p']);M[np.ix_(ids,ids)]+=b.T@(weight[:,None]*b);Q[np.ix_(ps,ps)]+=q.T@(weight[:,None]*q);mean[ps]+=q.T@weight
        for c,g in enumerate((gx,gz)):D[np.ix_(ps,c*n+ids)]+=q.T@(weight[:,None]*g)
        grad=(gx,gz);lap=gx.T@((weight*mu)[:,None]*gx)+gz.T@((weight*mu)[:,None]*gz)
        for a in range(2):
            for c in range(2):
                block=grad[c].T@((weight*mu)[:,None]*grad[a])
                if a==c:block+=lap;Klap[np.ix_(a*n+ids,c*n+ids)]+=lap
                K[np.ix_(a*n+ids,c*n+ids)]+=block
        upper=min(rr,.10854503542535461);lower=max(rl,0.)
        if upper<=lower:continue
        trace=next(v for v in data['elements']if abs(v['bounds'][0]-xl)<1e-10 and abs(v['bounds'][1]-xr)<1e-10 and v['bounds'][2]<=.15+1e-12 and v['bounds'][3]>.15+1e-12)
        tp=np.array(trace['p']);txl,txr,trl,trr=trace['bounds'];etaH=np.full_like(xi,2*(.15-trl)/(trr-trl)-1)
        # Exact physical tangential gradient of the trace: d_a /(1+h'^2).
        _,px,pz=tensor(2,xi,etaH,dx,trr-trl,slope,np.array(trace['ps']));tangent=(px+slope[:,None]*pz)/(1+slope*slope)[:,None]
        rz=(lower+upper)/2+eta*(upper-lower)/2;vb,_,_=tensor(3,xi,2*(rz-rl)/dz-1,dx,dz,slope,np.array(e['us']));cw=ww*dx*(upper-lower)/4
        C[np.ix_(ids,tp)]+=vb.T@(cw[:,None]*tangent);C[np.ix_(n+ids,tp)]+=vb.T@((cw*slope)[:,None]*tangent)
    return dict(mass=M,pressure_mass=Q,D=D,G=-D.T,C=C,K=K,K_component=Klap,mean=mean)


def analyze(path,curve,out):
    data=json.loads(path.read_text());ref=assemble(data,curve);coarse=assemble(data,curve,20);checks={}
    for name in ('mass','pressure_mass','D','G','C','K'):
        native=native_sparse(data[name]);scale=la.norm(ref[name]);error=la.norm(native-ref[name])/scale
        convergence=la.norm(ref[name]-coarse[name])/scale
        checks[name]=dict(relative_error=float(error),independent_quadrature_error=float(convergence));assert error<=1e-10 and convergence<=1e-11,(name,error,convergence)
    n,p,nd=data['velocity_size'],data['pressure_size'],data['dirichlet_size'];nf=n-nd;keep=np.r_[np.arange(nd,n),np.arange(n+nd,2*n)];g=data['G'];alpha=data['alpha'];awf=data['awf']
    A=ref['K'].copy();A[:n,:n]+=alpha*ref['mass'];A[n:,n:]+=alpha*ref['mass'];G=ref['G']-awf*ref['C'];D=ref['D'];gauge=data['gauge'];size=2*nf+p+int(gauge)
    F=np.zeros((size,size));F[:2*nf,:2*nf]=A[np.ix_(keep,keep)];F[:2*nf,2*nf:2*nf+p]=G[keep];F[2*nf:2*nf+p,:2*nf]=D[:,keep]
    if gauge:F[2*nf:2*nf+p,-1]=ref['mean'];F[-1,2*nf:2*nf+p]=ref['mean']
    raw=native_sparse(data['system']);checks['system']=dict(relative_error=float(la.norm(raw-F)/la.norm(F)));assert checks['system']['relative_error']<=1e-10
    rhs=np.array(data['rhs']);direct=la.solve(F,rhs);actual=np.r_[np.array(data['solution'])[keep],np.array(data['solution'])[2*n:],0.]if gauge else np.r_[np.array(data['solution'])[keep],np.array(data['solution'])[2*n:]]
    relative=la.norm(actual-direct)/la.norm(direct);assert relative<=1e-8,relative
    residual=F@actual-rhs;velocity_mass=ref['mass'][nd:,nd:];momentum=sum(v@la.solve(velocity_mass,v,assume_a='pos')for v in (residual[:nf],residual[nf:2*nf]));continuity=residual[2*nf:2*nf+p]
    if gauge:continuity-=ref['mean']*actual[-1]
    volume=float(np.sum(ref['mass']));# replaced below by physical Jacobian measure, independent of coefficient convention
    volume=sum((e['bounds'][1]-e['bounds'][0])*(e['bounds'][3]-e['bounds'][2])for e in data['elements']);L=np.sqrt(volume);force=1+1/L+(6.423173803526447e-5+.3819178635467331*.051*data['elements'][-1]['bounds'][3]/2)/(L*L)+.01/L
    momentum=float(np.sqrt(momentum/volume)/force);continuity=float(np.sqrt(max(0.,continuity@la.solve(ref['pressure_mass'],continuity,assume_a='pos'))/volume)*L);assert max(momentum,continuity)<=1e-9
    unaug=F[:2*nf+p,:2*nf+p];singular=la.svdvals(unaug);nulls=int(np.sum(singular<singular[0]*1e-12));assert nulls==int(gauge),(nulls,gauge,singular[-4:])
    one=np.zeros(p)
    for e in data['elements']:
        for j,gid in enumerate(e['p']):
            if j%3<2 and j//3<2:one[gid]=1/e['ps'][j]
    nullC=float(la.norm(ref['C']@one));assert nullC<1e-12
    wrong_strain=la.norm(ref['K']-ref['K_component'])/la.norm(ref['K']);wrong_awf=la.norm(2*ref['C'][keep])/la.norm(G[keep]);stale=la.norm(.2*ref['K'][np.ix_(keep,keep)])/la.norm(A[np.ix_(keep,keep)]);assert min(wrong_strain,wrong_awf,stale)>1e-8
    result=dict(status='PASS_TINY_NATIVE_INDEPENDENT_ASSEMBLY',curve=curve,source=str(path),checks=checks,direct_solution_relative_error=float(relative),independent_momentum_residual=momentum,independent_continuity_residual=continuity,unaugmented_null_modes=nulls,gauge=gauge,smallest_singular_values=singular[-5:].tolist(),constant_pressure_C_action=nullC,
        mutations=dict(missing_transpose=float(wrong_strain),wrong_AWF_sign=float(wrong_awf),stale_20percent_viscosity=float(stale)),native=data['solve'],current_coefficient_retained_preconditioner=data['changed_mu_retained_preconditioner'],scope='tiny algebraic/direct operator oracle; does not replace continuous MMS')
    assert data['changed_mu_retained_preconditioner']['rebuilds']==0
    out.write_text(json.dumps(result,indent=2)+'\n');print(result['status'],curve,'operator max',max(v['relative_error']for v in checks.values()),'direct',relative,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--native',type=Path,required=True);p.add_argument('--curve',type=float,default=0);p.add_argument('--out',type=Path,required=True);a=p.parse_args();analyze(a.native,a.curve,a.out)
