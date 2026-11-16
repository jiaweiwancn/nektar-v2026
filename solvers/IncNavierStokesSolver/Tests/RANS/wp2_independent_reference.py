# SPDX-License-Identifier: MIT
"""Independent flat analytical AWF/log reference, including inner ReferenceK.

The flat sigma_k=1, negative-intercept primitive below is an independent
integration of Q0/Dk. It does not project or change a represented field.
"""
import numpy as np
from scipy.special import spence
from wp2_all_state import US,wall
from cao_tamura_h2g_constitutive import parts
from run_cls_epawf_picard import require

Z0=.00012222264395414748


def evaluate(report,physical):
    p=np.asarray(physical);z=p[:,2];nu,H,hs=report['parameters'][:3];K0,q,kap=report['parameters'][5:8]
    require(report['ustar']==US and report['model_specification']['name'] in ['CLSHighRe','StandardKEpsilon'],'WP2 frozen model/reference scope')
    inner=z<H;r=np.minimum(z,H);wa=wall(np.full(len(z),K0),r,nu,H,hs,q,kap)
    zv=wa['intercept'];require(np.all(zv<0),'this independent primitive is for the frozen negative-intercept WP2 calibration')
    m=kap/np.sqrt(q)*np.sqrt(K0);D0=nu-m*zv;D=D0+m*r;L=np.log1p(m*r/D0)
    A=(K0/q)**1.5/kap;rc=wa['rc'];tau=K0/q
    def Fp(x):
        l=np.log1p(m*x/D0)
        return tau*tau/(m*m)*(.5*l*l-nu/D0*(l+D0/(D0+m*x)-1))
    def Fe(x):
        x=np.asarray(x);l=np.log1p(m*x/D0)
        low=A/(rc*m*m)*(m*x-D0*l)
        high=A/m*((1+np.log(np.maximum(x,rc)/rc))*l+spence(1+m*x/D0))
        lc=np.log1p(m*rc/D0)
        shift=A/(rc*m*m)*(m*rc-D0*lc)-A/m*(lc+spence(1+m*rc/D0))
        return np.where(x<=rc,low,high+shift)
    kr=K0-((Fe(H)-Fe(r))-(Fp(H)-Fp(r)))
    kr[~inner]=K0
    UH=US/kap*np.log(H/Z0)
    U=np.zeros(len(z))
    jH=np.log1p(m*H/D0)/m
    U[inner]=UH*(L[inner]/m)/jH[inner]
    U[~inner]=US/kap*np.log(z[~inner]/Z0)
    du=np.where(inner,UH/(jH*D),US/(kap*np.maximum(z,H)))
    er=wa['E'].copy();er[~inner]=US**3/(kap*z[~inner])
    G=np.zeros((len(z),3,3));G[:,0,2]=du
    law='calibrated_standard' if report['model_specification']['name']=='StandardKEpsilon' else 'regularized_CLS'
    c=parts(G,kr,er,inner,wa['nt'],p,US,kap,law)
    # The outer neutral stress kernel supplies the absolute datum and the
    # regional offset. No pressure mean is removed.
    gh=np.zeros((1,3,3));gh[0,0,2]=US/(kap*H)
    ch=parts(gh,np.array([K0]),np.array([US**3/(kap*H)]),np.array([False]),np.array([0.]),np.array([[0.,0.,H]]),US,kap,law)
    xx,zz=ch['tensor'][0,0,0],ch['tensor'][0,2,2];pout=report['developing_boundary']['definition']['Pout']
    pi=np.where(inner,pout+zz-xx,pout+2*K0/3-xx)
    T=-c['tensor']+(2*kr/3-pi)[:,None,None]*np.eye(3);T[inner]+=nu*(G[inner]+G[inner].transpose(0,2,1))
    values=np.column_stack((U,np.zeros(len(z)),kr,er,pi,pi-2*kr/3,c['nu_t'],T[:,0,2],T[:,2,2],c['production'],G[:,0,0],G[:,0,2],G[:,2,0],G[:,2,2]))
    require(np.all(np.isfinite(values)) and min(kr)>0 and min(er)>0,'independent analytical scalar/reference failure')
    return dict(values=values,G=G,R=c['tensor'],T=T)
