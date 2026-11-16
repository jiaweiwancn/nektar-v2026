# SPDX-License-Identifier: MIT
"""Independent GL refinement from exact affine cubic fields, not interpolated stresses.

Native GL8 samples determine the p<=3 fields on each common physical panel.
Reconstruct those fields, recompute the selected constitutive law, and integrate
at GL12. Inner dependent epsilon is not a transported comparison field.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss,legvander
from cao_tamura_h2g_constitutive import parts
from run_cls_epawf_picard import require


def resample(report,order=12):
    a=report['wp2_common_samples'];require(int(a['order'])==8,'qualified source quadrature is GL8')
    n=8;zknots=np.asarray(a['physical_z_partition']);nz=len(zknots)-1;raw=np.asarray(a['values'])
    nx=len(raw)//(n*n*nz);require(nx*n*n*nz==len(raw),'common tensor layout')
    # The six retained streamwise panels are shared by every WP2 member.
    xknots=np.array([-.2,-.1,-.05,0,.05,.1,.6]);require(nx==len(xknots)-1,'changed streamwise partition')
    old,wold=leggauss(n);new,wnew=leggauss(order);V=legvander(old,3)
    projection=(np.arange(4)*2+1)[:,None]/2*V.T*wold[None,:]
    interpolation=legvander(new,3)@projection
    native=raw.reshape(nx,n,nz,n,-1).transpose(0,2,1,3,4)
    selected=[0,1,2,3,4,10,11,12,13]
    base=native[:,:,:,:,selected]
    fit=np.einsum('ai,xyijf,bj->xyabf',V@projection,base,V@projection,optimize=True)
    outer=zknots[:-1]>=report['parameters'][1]
    err=abs(fit-base);err[:,~outer,:,:,3]=0.
    fit_error=float(np.max(err/np.maximum(1.,abs(base))))
    require(fit_error<2e-11,'samples do not represent the declared affine cubic fields')
    sampled=np.einsum('ai,xyijf,bj->xyabf',interpolation,base,interpolation,optimize=True)
    # Restore common output layout (x panel,x point,z panel,z point,field).
    sampled=sampled.transpose(0,2,1,3,4).reshape(-1,len(selected))
    x=xknots[:-1,None]+(new[None,:]+1)*np.diff(xknots)[:,None]/2
    z=zknots[:-1,None]+(new[None,:]+1)*np.diff(zknots)[:,None]/2
    X=np.broadcast_to(x[:,:,None,None],(nx,order,nz,order)).reshape(-1)
    Z=np.broadcast_to(z[None,None,:,:],(nx,order,nz,order)).reshape(-1)
    weights=(np.diff(xknots)[:,None,None,None]/2*wnew[None,:,None,None]*np.diff(zknots)[None,None,:,None]/2*wnew[None,None,None,:]).reshape(-1)
    nu,H,hs=report['parameters'][:3];q,kap=report['parameters'][6:8];us=.07426793792419734
    inner=Z<H;U,W,k,e,pi=[sampled[:,j].copy() for j in range(5)];e[inner]=1.
    G=np.zeros((len(k),3,3));G[:,0,0]=sampled[:,5];G[:,0,2]=sampled[:,6];G[:,2,0]=sampled[:,7];G[:,2,2]=sampled[:,8]
    wall=np.asarray(report['physical_samples']['wall_K']).reshape(nx,n)
    KH=(wall@interpolation.T)[:,:,None,None];KH=np.broadcast_to(KH,(nx,order,nz,order)).reshape(-1)
    root=np.sqrt(KH);ratio=hs*root/(70*nu);exponent=np.maximum(.5-.4*ratio**.7,1-.79*ratio**(-.28))
    intercept=-10.7*np.expm1(exponent*np.log(ratio));zv=nu*intercept/root
    nt=(q**-2)*(kap*q**1.5)*root*np.maximum(0.,Z-zv)
    law='calibrated_standard' if report['model_specification']['name']=='StandardKEpsilon' else 'regularized_CLS'
    physical=np.column_stack((X,np.zeros(len(X)),Z));v=parts(G,k,e,inner,nt,physical,us,kap,law)
    T=-v['tensor']+(2*k/3-pi)[:,None,None]*np.eye(3)
    T[inner]+=nu*(G[inner]+G[inner].transpose(0,2,1))
    values=np.column_stack((U,W,k,np.where(inner,0.,e),pi,pi-2*k/3,v['nu_t'],T[:,0,2],T[:,2,2],v['production'],G[:,0,0],G[:,0,2],G[:,2,0],G[:,2,2]))
    return dict(values=values,weights=weights,regions=(~inner).astype(int),fit_error=fit_error,R=v['tensor'],T=T,G=G,
                physical=physical,endpoint_K=KH,
                scope='independent recomputation on affine polynomial fields; inner nu_A from actual endpoint K; dependent inner epsilon excluded')


def verify_native_tuple_replay(report):
    a=resample(report,8);b=np.asarray(report['wp2_common_samples']['values']);d=abs(a['values']-b)
    # Inner epsilon is dependent and deliberately not interpolated/reported here.
    d[a['regions']==0,3]=0.
    maximum=np.max(d,axis=0);require(max(maximum[6:10])<2e-10,'independent native viscosity/stress/production replay failed')
    return dict(maximum_absolute_by_column=maximum.tolist(),polynomial_fit_error=a['fit_error'],passed=True)


def compare_orders(left,right):
    a,b=resample(left),resample(right);require(np.array_equal(a['weights'],b['weights']) and np.array_equal(a['regions'],b['regions']),'changed physical integration')
    native_a,native_b=left['wp2_common_samples'],right['wp2_common_samples'];v=np.asarray(native_a['values'])-np.asarray(native_b['values']);w=np.asarray(native_a['weights']);ref=np.asarray(native_a['reference_values'])
    out={};us=.07426793792419734
    for region in [0,1]:
        old=np.asarray(native_a['regions'])==region;new=a['regions']==region;measure=sum(w[old]);columns={}
        for j,name in enumerate(native_a['columns'].split(',')):
            if j==3 and region==0:continue
            n8=float(np.sqrt(np.dot(w[old],v[old,j]**2)))
            n12=float(np.sqrt(np.dot(a['weights'][new],(a['values'][new,j]-b['values'][new,j])**2)))
            scale=us if j==1 else us*us if j in [4,5,7,8] else us/.0065 if j in [10,12,13] else None
            denominator=scale*np.sqrt(measure) if scale else np.sqrt(np.dot(w[old],ref[old,j]**2))
            columns[name]=dict(L2_GL8=n8,L2_GL12=n12,normalized_change=abs(n12-n8)/denominator)
        out[str(region)]=columns
    c,d=resample(left,8),resample(right,8)
    tensors={}
    for name in ['R','T']:
        row={}
        for region in [0,1]:
            old=c['regions']==region;new=a['regions']==region;volume=sum(c['weights'][old]);den=us**2*np.sqrt(volume)
            n8=float(np.sqrt(np.sum(c['weights'][old,None,None]*(c[name][old]-d[name][old])**2))/den)
            n12=float(np.sqrt(np.sum(a['weights'][new,None,None]*(a[name][new]-b[name][new])**2))/den)
            row[str(region)]=dict(normalized_GL8=n8,normalized_GL12=n12,normalized_change=abs(n12-n8))
        tensors[name]=row
    out['full_3x3_tensors']=tensors
    return out
