#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent dimensional tensor / breakpoint-quadrature oracle (Linux)."""
import argparse
import json
import math
import os
from pathlib import Path
import platform
import subprocess

THREADS = dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'), '1')
os.environ.update(THREADS)
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from numpy.polynomial.legendre import Legendre, leggauss
from a93_a95_e94e95_note_wall_oracle import Kernel

CHECKS = 0
ERRORS = {}

def require(ok, name):
    global CHECKS
    CHECKS += 1
    if not ok: raise AssertionError(name)

def near(a,b,name,tol=3e-10,scale=1.):
    a,b=np.asarray(a),np.asarray(b)
    error=float(np.max(np.abs(a-b)))/max(scale,float(np.max(np.abs(b))))
    ERRORS[name]=max(ERRORS.get(name,0.),error)
    require(np.isfinite(error) and error<tol,(name,error,tol))

def constants(name,policy):
    def c(q):
        return (2/3)/(1.25+1.9*q) if name=='A93' else 1/((6.5 if name=='A95' else 4)+3*q/math.sqrt(2))
    q=brentq(lambda q:c(q)*q*q-1,2,5,xtol=5e-15)
    delta=.48 if name!='E94E95' else 1.9-.43*q
    original=1.3 if name!='E94E95' else 1.2
    kap=math.sqrt(original*delta/q) if policy=='Natural' else .41
    return q,kap,original if policy=='Natural' else kap*kap*q/delta

def dev(a):return a-np.trace(a)*np.eye(3)/3

def oracle(name,k,e,A):
    # Physical, dimensional tensors. Eigenvalues rather than the production
    # evaluator's invariant algorithm. A93 uses the published A-polynomial.
    S=dev((A+A.T)/2);O=(A-A.T)/2;t=k/e
    s=np.linalg.norm(S);r=np.linalg.norm(O);eta=math.sqrt(2)*s*t
    if name=='A93':
        C=(2/3)/(1.25+eta+.9*math.sqrt(2)*r*t)
        T=-4*dev(A@A+(A@A).T)+13*dev(A@A.T)-2*dev(A.T@A)
        # Use trace-free A for the independent polynomial extension.
        B=S+O
        T=-4*dev(B@B+(B@B).T)+13*dev(B@B.T)-2*dev(B.T@B)
        Q=k*t*t/(1000+eta**3)*T
        Q2=k*t*t/(1000+eta**3)*(3*dev(S@S)-19*dev(O@O)+15*(O@S-S@O))
        near(Q,Q2,'independent_A93_tensor_forms')
    else:
        if s==0 and r!=0:raise ValueError('undefined_pure_rotation')
        As=3*np.linalg.eigvalsh(S)[-1]/s if s else 0.
        C=1/((6.5 if name=='A95' else 4)+As*math.hypot(s,r)*t)
        f=1-9*C*C*s*s*t*t
        if name=='A95' and f<0:raise ValueError('negative_radicand')
        Q=2*math.sqrt(f)/(1+6*s*r*t*t)*k*t*t*(O@S-S@O) if name=='A95' else np.zeros((3,3))
    nut=C*k*t;R=2*k*np.eye(3)/3-2*nut*S+Q
    P=-np.sum(R*A)
    ep=(max(.43,eta/(eta+5))*math.sqrt(2)*s*e if name=='E94E95' else 1.44*e/k*P)
    ed=(1.9 if name=='E94E95' else 1.92)*e*e/k
    return dict(R=R,Q=Q,cmu=C,nut=nut,production=P,k_source=P-e,eps_production=ep,eps_destruction=ed,eps_source=ep-ed)

class Probe:
    def __init__(self,path):self.p=subprocess.Popen([str(path)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True,bufsize=1)
    def query(self,line):
        self.p.stdin.write(line+'\n');self.p.stdin.flush()
        result=self.p.stdout.readline()
        require(bool(result),'probe terminated')
        return json.loads(result)
    def closure(self,name,A,k=1.,e=1.,policy='Natural'):
        return self.query('closure '+name+' '+policy+' EP .000015 '+ ' '.join(map(str,[k,e,*np.asarray(A).ravel()])))
    def wall(self,name,policy,wall,K,H=2.,h=.1,nu=1.5e-5,chi=1.,zs=(0.,2.),U=None,F=(0.,0.),order=64):
        q,kap,_=constants(name,policy);kw=Kernel(K,nu,H,h,.2295 if wall=='EP' else kap/math.sqrt(q),2.55 if wall=='EP' else kap*q**1.5)
        U=(K/q*kw.J0(H),0.) if U is None else U
        values=[K,H,h,H,H,chi,order,*U,*F,len(zs),*zs]
        return self.query('wall '+name+' '+policy+' '+wall+' '+str(nu)+' '+' '.join(map(str,values)))
    def close(self):self.p.stdin.close();require(self.p.wait()==0,'probe exit')

def closure_tests(p):
    rng=np.random.default_rng(93);records=[]
    for name in ('E94E95','A95','A93'):
        for policy in ('Natural','Prescribed'):
            q,kap,sigma=constants(name,policy)
            A=np.zeros((3,3));A[0,2]=q
            c=p.closure(name,A,policy=policy)
            near([c['q'],c['kappa'],c['sigma']],[q,kap,sigma],'equilibrium_descriptor')
            near(c['production'],1.,'simple_shear_root')
            near(kap*kap*q/sigma+c['eps_source'],0.,'outer_epsilon_balance')
            near(c['diffusion'],[c['nut'],c['nut'],c['nut']/sigma],'no_outer_molecular_diffusion')
            records.append(dict(model=name,policy=policy,q=q,kappa=kap,sigma=sigma))
        gradients=[np.zeros((3,3)),np.diag([1,1,-2]),np.diag([1,-.5,-.5])]
        gradients += [dev(rng.normal(size=(3,3)))*.3 for _ in range(12)]
        for A in gradients:
            actual=p.closure(name,A,1.2,.7)
            try:expected=oracle(name,1.2,.7,A)
            except ValueError:
                require(not actual['ok'],'independent domain rejection agrees')
                continue
            require(actual['ok'],('admissible tensor rejected',name,A,actual))
            for key,val in expected.items():near(np.array(actual[key]).reshape(val.shape) if isinstance(val,np.ndarray) else actual[key],val,'tensor_'+key)
            R=np.array(actual['R']).reshape(3,3)
            near(R,R.T,'stress_symmetry');near(np.trace(R),2.4,'full_3D_trace')
            near(actual['cubic'],np.zeros(9),'no_CLS_cubic');near(actual['yap'],0.,'no_Yap')
            rotation=np.linalg.qr(rng.normal(size=(3,3)))[0]
            rotated=p.closure(name,rotation@A@rotation.T,1.2,.7)
            near(np.array(rotated['R']).reshape(3,3),rotation@R@rotation.T,'3D_covariance')
            # Independent strictly trace-free oracle then raw-A production.
            divergent=p.closure(name,A+.37*np.eye(3),1.2,.7)
            near(divergent['R'],R.ravel(),'DeviatoricStrainV1_stress')
            near(divergent['production'],actual['production']-2*1.2*.37,'raw_divergence_pressure_work')
            near(divergent['divergence'],1.11,'raw_divergence_reported')
        O=np.zeros((3,3));O[0,2]=8;O[2,0]=-8
        rotation=p.closure(name,O)
        if name=='A93':require(rotation['ok'] and np.array(rotation['R']).reshape(3,3)[1,1]<0,'A93 covariance negative control')
        else:require(not rotation['ok'],'pure rotation must be rejected')
        require(not p.closure(name,np.zeros((3,3)),k=0)['ok'],'nonpositive k rejection')
        guard=p.query('tangent '+name+' Natural EP .000015 1 1 0 0 1 0 0 0 0 0 0')
        require(not guard['ok'] and 'unsupported' in guard['error'].lower(),'unsupported tangent guard')
    boundary=6.5/(3*(math.sqrt(5.42)-1))
    for sign in (-1,1):
        A=boundary*(1+sign*1e-10)*np.diag([1.,.9,-1.9])
        require(p.closure('A95',A)['ok']==(sign<0),'A95 square-root forward boundary')
    bad=np.diag([2.,1.8,-3.8]);require(not p.closure('A95',bad)['ok'],'A95 negative radicand rejection')
    A=np.diag([1.,1.,-2.]);A[0,1]=5;A[1,0]=-5
    c=p.closure('A93',A);require(c['production']<0 and np.linalg.eigvalsh(np.array(c['R']).reshape(3,3))[0]>0,'signed A93 production control')
    near(c['production'],-.3114625443,'A93 negative production value',tol=2e-9)
    # Differentiation of production samples is compared to exact note acoustic
    # eigenvalues, including all coefficient dependence; no tangent is enabled.
    for name,A,a,n,target in [
        ('A93',np.array([[0.,0,8],[0,0,0],[0,0,0]]),np.array([1,0,-2])/math.sqrt(5),np.array([2,0,1])/math.sqrt(5),-51301/296100),
        ('A95',np.array([[0.,0,1],[0,0,0],[0,0,0]]),np.array([0,1,0]),np.array([1,0,1])/math.sqrt(2),-.044533122876222),
        ('E94E95',np.diag([10.,9.,-19.]),np.array([1,0,-1])/math.sqrt(2),np.array([1,0,1])/math.sqrt(2),-19/68)]:
        H=np.outer(a,n);h=1e-5
        plus=np.array(p.closure(name,A+h*H)['R']).reshape(3,3);minus=np.array(p.closure(name,A-h*H)['R']).reshape(3,3)
        c=p.closure(name,A);principal=-np.sum(H*(plus-minus))/(2*h*c['nut'])
        near(principal,target,'negative_principal_'+name,tol=5e-7)
    threshold=5*.43/(1-.43)
    for eta in (threshold-.01,threshold,threshold+.01,8.):
        A=np.zeros((3,3));A[0,2]=eta;c=p.closure('E94E95',A)
        near(c['eps_production'],max(.43,eta/(eta+5))*eta,'E_distinct_source_branch')
    return records

def wall_tests(p):
    records=[]
    for name in ('E94E95','A95','A93'):
        for policy in ('Natural','Prescribed'):
            q,kap,_=constants(name,policy);K=q*.4**2
            for wall in ('EP','Accepted'):
                H=2.;nu=1.5e-5;h=.1
                w=Kernel(K,nu,H,h,.2295 if wall=='EP' else kap/math.sqrt(q),2.55 if wall=='EP' else kap*q**1.5)
                t=.4**2;amplitude=.4**3/kap
                lam=w.integral(lambda z:t*t*w.nt(z)/w.D(z)**2)/amplitude
                zc=H*math.exp(1-lam)
                breaks=sorted(set([0.,H,zc,h,*[x for x in (w.zv,w.ze,.0025,.01,.03,.1,.3,.7) if 0<x<H]]))
                def integ(f,upper=H):
                    br=[z for z in breaks if z<upper]+[upper]
                    return sum(quad(f,a,b,epsabs=3e-12,epsrel=3e-12)[0] for a,b in zip(br[:-1],br[1:]))
                p0=lambda z:t*t*w.nt(z)/w.D(z)**2
                Ew=lambda z:amplitude/max(z,zc)
                Q=lambda z:integ(lambda s:Ew(s)-p0(s),z)
                kg=lambda z:K-integ(lambda s:Q(s)/w.D(s))+integ(lambda s:Q(s)/w.D(s),z)
                nodes,weights=leggauss(24);zs=[0.,H,.01,.2,zc];wg=[0.]*len(zs)
                for a,b in zip(breaks[:-1],breaks[1:]):
                    zs.extend((a+b)/2+(b-a)*nodes/2);wg.extend((b-a)*weights/2)
                c=p.wall(name,policy,wall,K,zs=zs);require(c['ok'],c)
                near(c['J'],[w.J0(H),w.J1(H),w.J2(H)],'wall_integrals',scale=1e-6)
                near(c['cutoff'],zc,'finite_nu_cutoff',scale=.01)
                samples=np.array(c['samples']);z=samples[:,0];ww=np.array(wg)
                near(samples[:,1],list(map(p0,z)),'wall_p0')
                near(np.dot(ww,samples[:,3]),c['production_budget'],'dissipation_budget',tol=3e-9)
                near(t*c['UH'],c['production_budget']+nu*t*t*w.J2(H),'finite_nu_mechanical_work')
                near(samples[0,9:11],[0,0],'no_slip');near(samples[1,9],t*w.J0(H),'upper_velocity')
                if wall=='EP':
                    near(samples[:,5],np.zeros(len(z)),'EP_plateau_cancellation')
                    near(samples[:,6],np.full(len(z),K),'EP_constant_k')
                else:
                    near(samples[:5,6],[kg(x) for x in z[:5]],'accepted_reference_k',tol=2e-9)
                    near(samples[:,7],list(map(Q,z)),'accepted_reference_flux',tol=3e-9)
                    near(samples[:,4],np.zeros(len(z)),'accepted_no_redistribution')
                    require(np.min(samples[:,6])>=K-1e-10,'rough_reference_positive')
                    for degree in range(7):
                        phi=Legendre.basis(degree);dphi=phi.deriv()
                        source=np.dot(ww,phi(2*z/H-1)*samples[:,5]);diff=-np.dot(ww,2/H*dphi(2*z/H-1)*samples[:,7])
                        near(source+diff,0.,'accepted_weak_mode_'+str(degree),tol=2e-8)
                        if degree==1:require(abs(source)>1e-5,'source moments not separately zero')
                # Actual forcing, including reversed ground shear; both components.
                B=w.J1(H)/w.J0(H);F=np.array([1.5*t/B,.35*t/B]);U=np.array([t*w.J0(H),0.]);tw=U/w.J0(H)-B*F
                forced=p.wall(name,policy,wall,K,zs=zs,U=U,F=F);fs=np.array(forced['samples'])
                work=np.dot(U,tw+F*(H-h))
                loss=np.dot(ww,[w.D(x)*np.dot(g,g) for x,g in zip(z,fs[:,11:13])])
                forcing=np.dot(ww,[np.dot(F,u) if x>h else 0. for x,u in zip(z,fs[:,9:11])])
                near(work,loss+forcing,'forced_reversed_mechanical_work',tol=2e-8)
                require(tw[0]<0,'signed_reversed_wall_shear')
                if wall=='Accepted':
                    zero=p.wall(name,policy,wall,K,zs=(.2,),U=(0.,0.))
                    require(zero['ok'] and zero['samples'][0][5]<0,'no_general_dynamic_positivity_guarantee')
                records.append(dict(model=name,policy=policy,wall=wall,Lambda=lam,cutoff=zc,ground_k=float(samples[0,6]),status='pass'))
            # Budget endpoint: use a bracket, then approach Lambda=1 from above.
            alpha=kap/math.sqrt(q)
            def gate(H):
                w=Kernel(K,nu,H,0.,alpha,kap*q**1.5)
                return alpha*math.sqrt(K)*(w.J0(H)-nu*w.J2(H))-1
            H1=brentq(gate,.0001,.1,xtol=1e-16)
            c=p.wall(name,policy,'Accepted',K,H=H1*(1+1e-12),h=0.,zs=(0.,H1*(1+1e-12)))
            require(c['ok'],'Lambda=1 endpoint accepted from admissible side')
            near(c['cutoff'],H1,'Lambda_one_cutoff_equals_H',tol=1e-9,scale=H1)
            require(not p.wall(name,policy,'Accepted',K,H=H1*.99,h=0.,zs=(0.,))['ok'],'Lambda_less_one_rejection')
            # Smooth reference positivity counterexamples, Lambda=4 exactly.
            y4=brentq(lambda x:math.log1p(x)-x/(1+x)-4,100,200)
            alpha=kap/math.sqrt(q);H=nu*(10.7*alpha+y4)/(alpha*math.sqrt(K))
            c=p.wall(name,policy,'Accepted',K,H=H,h=0.,zs=(0.,H))
            require(c['ok'] and c['samples'][0][6]<0,'smooth budget-admissible negative reference')
            require(not p.wall(name,policy,'Accepted',K,H=1e-6,h=0.,zs=(0.,1e-6))['ok'],'molecular-only accepted cutoff rejection')
            # chi=0 must bypass an inactive, inadmissible corrected branch.
            c=p.wall(name,policy,'EP',K,H=1e-6,h=0.,chi=0.,zs=(0.,1e-6))
            require(c['ok'],'direct inactive EP baseline')
    return records

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--probe',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    a=parser.parse_args();require(platform.system()=='Linux','Linux authoritative')
    p=Probe(a.probe);result=dict(schema=1,scope='independent dimensional tensors and scalar wall quadrature; no native PDE claim')
    try:
        result['closures']=closure_tests(p);result['walls']=wall_tests(p);result['status']='pass'
    except Exception as e:
        result['status']='fail';result['failure']=repr(e);raise
    finally:
        p.close();result.update(checks=CHECKS,maximum_scaled_errors=ERRORS,threads=THREADS)
        a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],checks=CHECKS)))
if __name__=='__main__':main()
