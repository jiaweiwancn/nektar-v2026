# SPDX-License-Identifier: MIT
"""Independent fixed-tuple laws. Initialization model selection is explicit."""
import numpy as np
from run_cls_epawf_picard import require
from search_cao_tamura_budget_h2d import tensors as deployed_cls

LAWS=['regularized_CLS','CLS_beta04_no_regularization','Original_CLS','local_CLS_Cmu_linear','local_CLS_Cmu_linear_regularized','calibrated_standard']

def parts(g,k,e,inner,nt,physical,us,kap,law,branch='physical'):
    require(law in LAWS,'unknown frozen constitutive law');S=g+g.transpose(0,2,1);O=g-g.transpose(0,2,1);I=np.eye(3);safe=np.where(inner,1.,e);ss=np.einsum('nij,nij->n',S,S);oo=np.einsum('nij,nij->n',O,O);q=k/safe*np.sqrt((ss if branch=='strain' else oo if branch=='rotation' else np.maximum(ss,oo))/2)
    with np.errstate(over='ignore',under='ignore'):cm=.3/(1+.35*q**1.5)*(-np.expm1(-.36*np.exp(np.minimum(.75*q,700))))
    if law=='calibrated_standard':cm=np.full_like(k,.09)
    nut=cm*k*k/safe;iso=2*k[:,None,None]*I/3;linear=-nut[:,None,None]*S;Q=np.zeros_like(g);C=np.zeros_like(g);reg=np.zeros_like(g)
    if law in LAWS[:3]:
        Q=(k*nut/safe)[:,None,None]*(-.1*(S@S-ss[:,None,None]*I/3)+.1*(O@S-S@O)+.26*(-O@O-oo[:,None,None]*I/3))
        beta=1 if law=='Original_CLS' else .4
        C=beta*(nut*k*k/safe**2)[:,None,None]*((-10*cm*cm)[:,None,None]*(S@S@O-O@S@S)+(5*cm*cm*(oo-ss))[:,None,None]*S)
    if law in ('regularized_CLS','local_CLS_Cmu_linear_regularized'):
        z=np.where(inner,1.,physical[:,2]);sr=np.zeros_like(g);sr[:,0,2]=sr[:,2,0]=us/(kap*z);reg=-(kap*us*z)[:,None,None]*(S-sr)
    linear[inner]=-nt[inner,None,None]*S[inner];Q[inner]=0;C[inner]=0;reg[inner]=0;nut[inner]=nt[inner]
    R=iso+linear+Q+C+reg
    return dict(tensor=(R+R.transpose(0,2,1))/2,isotropic=iso,linear=linear,quadratic=Q,cubic=C,regularization=reg,nu_t=nut,Cmu=cm,q_strain=k/safe*np.sqrt(ss/2),q_rotation=k/safe*np.sqrt(oo/2),production=-np.einsum('nij,nij->n',R,g))

def tensors_for_spec(g,k,e,inner,nt,physical,us,kap,spec):
    name=spec['name']
    if name=='CLSHighRe':return deployed_cls(g,k,e,inner,nt,physical,us,kap)
    require(name=='StandardKEpsilon' and spec['calibration']=='StandardNeutralCalibratedV1' and not spec['physical_regularization_active'],'CLS context cannot be presented as calibrated standard')
    S=g+g.transpose(0,2,1);nut=.09*k*k/np.where(inner,1.,e);nut=np.where(inner,nt,nut)
    return 2*k[:,None,None]*np.eye(3)/3-nut[:,None,None]*S
