#!/usr/bin/env python3
"""Independent finite-nu SCI/EP-AWF identity checks; no CFD solve.
SCI roughness correlation retained; original and explicitly adjusted wall
coefficients are checked with natural and common kappa. Each reference
calibrates h once at K0; no evolving intercept calibration.
Analytic resistance formulas are compared with adaptive quadrature;
Legendre weak moments use separate breakpoint-aware Gauss quadrature.
"""
import json
import math
from dataclasses import dataclass
import numpy as np
from numpy.polynomial.legendre import Legendre, leggauss
from scipy.integrate import quad
from scipy.optimize import brentq

ALPHA=0.2295
CLW=2.55
VS=10.7
KAPPA=0.41
ROOTS={
 'A93':(57+math.sqrt(6249))/40,
 'A95':(3+math.sqrt(61))/(2*math.sqrt(2)),
 'E94_E95':(3+math.sqrt(41))/(2*math.sqrt(2)),
}
COUNT=0

def close(name, actual, expected, atol=3e-11, rtol=3e-10):
 global COUNT
 aa=np.asarray(actual,dtype=float); ee=np.asarray(expected,dtype=float)
 err=float(np.max(np.abs(aa-ee)))
 limit=atol+rtol*float(np.max(np.abs(ee)))
 if not err<=limit: raise AssertionError((name,actual,expected,err,limit))
 COUNT+=1
 return err

def require(name, condition):
 global COUNT
 if not condition:raise AssertionError(name)
 COUNT+=1

def sci_b(h,K,nu):
 if h==0:return VS
 x=h*math.sqrt(K)/(70*nu)
 exponent=max(0.5-0.4*x**0.7,1-0.79*x**(-0.28))
 return VS*(1-x**exponent)

@dataclass
class Kernel:
 K:float
 nu:float
 H:float
 h:float
 alpha:float=ALPHA
 clw:float=CLW
 def __post_init__(self):
  self.a=self.alpha*math.sqrt(self.K)
  self.b=sci_b(self.h,self.K,self.nu)
  self.zv=self.nu*self.b/math.sqrt(self.K)
  self.ze=2*self.clw*self.nu/math.sqrt(self.K)
 def nt(self,z):return self.a*max(0.,z-self.zv)
 def D(self,z):return self.nu+self.nt(z)
 def active(self,z):
  ell=min(max(self.zv,0.),z)
  return ell,self.D(ell),self.D(z)
 def J0(self,z):
  ell,dl,dz=self.active(z)
  return ell/self.nu+(math.log1p((dz-dl)/dl)/self.a if z>ell else 0.)
 def J2(self,z):
  ell,dl,dz=self.active(z)
  return ell/self.nu**2+((1/dl-1/dz)/self.a if z>ell else 0.)
 def J1(self,z):
  if z<=self.h:return 0.
  # Integral from h to the end of the molecular-only part.
  end=min(z,max(self.zv,0.))
  result=max(end-self.h,0.)**2/(2*self.nu)
  start=max(self.h,self.zv,0.)
  if z>start:
   ds=self.D(start);dz=self.D(z)
   result+=((z-start)-(self.h-self.zv+self.nu/self.a)*math.log1p((dz-ds)/ds))/self.a
  return result
 def Ebase(self,z):return self.K**1.5/(self.clw*max(z,self.ze))
 def breaks(self,extra=()):
  return sorted(set([0.,self.H]+[float(x) for x in (self.zv,self.h,self.ze,*extra) if 0<x<self.H]))
 def integral(self,fun,extra=(),upper=None):
  upper=self.H if upper is None else upper
  points=[p for p in self.breaks(extra) if p<upper]+[upper]
  return sum(quad(fun,a,b,epsabs=2e-12,epsrel=2e-12,limit=100)[0] for a,b in zip(points[:-1],points[1:]) if b>a)
 def Pint(self,z,t0):
  # Direct antiderivative of the turbulent, rather than total, work.
  ell,dl,dz=self.active(z)
  if z<=ell:return 0.
  ratio=(dz-dl)/dl
  return t0*t0/self.a*(math.log1p(ratio)-self.nu/dl*ratio/(1+ratio))

def gauss_integral(kernel,fun,extra=(),degree=72):
 x,w=leggauss(degree)
 total=0.
 br=kernel.breaks(extra)
 for a,b in zip(br[:-1],br[1:]):
  z=(a+b)/2+(b-a)*x/2
  total+=(b-a)/2*np.dot(w,[fun(float(v)) for v in z])
 return float(total)

def calibrate_h(q,u,nu,H,z0,kappa=KAPPA,alpha=ALPHA,clw=CLW):
 K=q*u*u; target=math.log(H/z0)/(kappa*u)
 def defect(h):return Kernel(K,nu,H,h,alpha,clw).J0(H)-target
 left,right=0.,H
 dl,dr=defect(left),defect(right)
 endpoint_tol=32*math.ulp(max(target,Kernel(K,nu,H,0.,alpha,clw).J0(H)))
 require('one-time SCI h calibration bracket',
         dl>=-endpoint_tol and dr<0 and target<H/nu)
 h=0. if abs(dl)<=endpoint_tol else brentq(defect,left,right,xtol=1e-14,rtol=2e-14)
 kernel=Kernel(K,nu,H,h,alpha,clw)
 require('physical calibrated roughness and positive production',
         0<=h<H and max(kernel.zv,0.)<H)
 return kernel

def outer_stress(name,q,u):
 K=q*u*u
 if name=='A93':
  fac=q*q/(1000+q**3)
  diag=np.array([2/3+28*fac/3,2/3-11*fac/3,2/3-17*fac/3])*K
 elif name=='A95':
  cq=math.sqrt(1-4.5/q**2)/(1+3*q*q)
  diag=np.array([2/3+cq*q*q,2/3,2/3-cq*q*q])*K
 else:diag=np.full(3,2*K/3)
 R=np.diag(diag);R[0,2]=R[2,0]=-u*u
 return R

def check_reference(label,kernel,q,u,z0=None,model=None,expect_gate=True,kappa=KAPPA):
 H=kernel.H;t0=u*u;A=u**3/kappa
 # Independently integrate all three resistance functions, including partial heights.
 for z in [.13*H,.47*H,H]:
  close(label+' J0',kernel.J0(z),kernel.integral(lambda x:1/kernel.D(x),upper=z))
  close(label+' J2',kernel.J2(z),kernel.integral(lambda x:1/kernel.D(x)**2,upper=z),atol=1e-9)
  close(label+' J1',kernel.J1(z),kernel.integral(lambda x:max(x-kernel.h,0.)/kernel.D(x),upper=z))
 J0=kernel.J0(H);J2=kernel.J2(H);UH=t0*J0
 z0gen=H*math.exp(-kappa*u*J0)
 if z0 is not None:
  close(label+' prescribed z0',z0gen,z0,atol=2e-13)
  close(label+' logarithmic anchor',UH,u/kappa*math.log(H/z0))
 else:z0=z0gen
 p0=lambda z:kernel.nt(z)*(t0/kernel.D(z))**2
 prod=kernel.integral(p0)
 mol=kernel.integral(lambda z:kernel.nu*(t0/kernel.D(z))**2)
 close(label+' production identity',prod,t0*t0*(J0-kernel.nu*J2),atol=4e-12)
 close(label+' positive production antiderivative',prod,kernel.Pint(H,t0))
 close(label+' molecular work',mol,kernel.nu*t0*t0*J2)
 close(label+' reference mean work',t0*UH,prod+mol)
 # SCI's original finite dissipation cutoff and continuous plateau.
 close(label+' baseline cutoff continuity',2*kernel.nu*kernel.K/kernel.ze**2,kernel.K**1.5/(kernel.clw*kernel.ze))
 expected_base=(kernel.K**1.5/kernel.clw*(1+math.log(H/kernel.ze)) if H>=kernel.ze else H*2*kernel.nu*kernel.K/kernel.ze**2)
 close(label+' baseline integral',kernel.integral(kernel.Ebase),expected_base)
 Lambda=prod/A;muwork=kappa*u*kernel.nu*J2
 close(label+' molecular Lambda correction',Lambda,math.log(H/z0)-muwork)
 gate=Lambda>=1+math.log(2.)
 require(label+' admissibility expectation',gate==expect_gate)
 record=dict(case=label,H=H,h=kernel.h,b=kernel.b,zv=kernel.zv,z0=z0,
             J0=J0,production=prod,molecular_loss=mol,Lambda=Lambda,
             molecular_Lambda_correction=muwork,tail_gate_at_H_over_2=gate)
 if not gate:
  if kernel.zv>=H:
   close(label+' all molecular zero turbulent production',prod,0.,atol=1e-14,rtol=0)
   close(label+' all molecular power',t0*UH,mol)
  return record
 zc=H*math.exp(1-Lambda)
 close(label+' finite-nu cutoff',zc,math.e*z0*math.exp(muwork))
 require(label+' finite-nu cutoff exceeds e*z0',zc>math.e*z0)
 Ec=lambda z:A/max(z,zc)
 def Ecint(z):return A*z/zc if z<=zc else A*(1+math.log(z/zc))
 jc=lambda z:Ecint(z)-kernel.Pint(z,t0)
 close(label+' epsilon budget',kernel.integral(Ec,[zc]),prod)
 close(label+' zero redistribution lower',jc(0),0.,atol=1e-13)
 close(label+' zero redistribution upper',jc(H),0.,atol=1e-11)
 # All nonconstant moments: check integration-by-parts realization, not only
 # cancellation of identical pointwise arithmetic expressions.
 raw=[];corrected=[]
 for degree in range(7):
  phi=Legendre.basis(degree);dphi=phi.deriv()
  raw_m=gauss_integral(kernel,lambda z:phi(2*z/H-1)*(p0(z)-Ec(z)),[zc])
  flux_m=-gauss_integral(kernel,lambda z:(2/H)*dphi(2*z/H-1)*jc(z),[zc])
  corrected.append(raw_m+flux_m);raw.append(raw_m)
  close(label+' corrected weak mode '+str(degree),raw_m+flux_m,0.,atol=2e-10,rtol=0)
 require(label+' mean-only correction leaves nonconstant mode',max(abs(x) for x in raw[1:])>1e-6*max(prod,1e-30))
 close(label+' positive corrected dissipation trace',Ec(H),A/H)
 record.update(zc=zc,max_uncorrected_nonconstant_moment=max(abs(x) for x in raw[1:]),
               max_corrected_weak_moment=max(abs(x) for x in corrected))
 if model:
  # Molecular stress belongs to total stress, not Reynolds covariance.
  G=t0/kernel.D(H)
  S=np.zeros((3,3));S[0,2]=S[2,0]=G/2
  RA=2*kernel.K/3*np.eye(3)-2*kernel.nt(H)*S
  Ro=outer_stress(model,q,u)
  piw=Ro[2,2]-2*kernel.K/3
  Tw=-piw*np.eye(3)+2*kernel.nu*S-RA
  To=-Ro
  ez=np.array([0.,0.,1.])
  close(label+' full interface traction',Tw@ez,To@ez)
  close(label+' inner Reynolds shear excludes molecular',RA[0,2],-t0+kernel.nu*G)
  U=np.array([UH,0.,0.])
  close(label+' internal interface work cancels',U@(Tw@ez)+U@(To@(-ez)),0.)
  require(label+' inner reference covariance positive',np.linalg.eigvalsh(RA)[0]>0)
  record['normal_pressure_offset']=float(piw)
 return record

def check_forced(label,kernel,q,u):
 H=kernel.H;t0=u*u;J0=kernel.J0(H);B=kernel.J1(H)/J0
 require(label+' nonzero force transfer',B>0)
 UH=np.array([t0*J0,0.])
 F=np.array([1.5*t0/B,.35*t0/B])
 tw=UH/J0-B*F
 TH=tw+F*(H-kernel.h)
 require(label+' signed reversed ground shear',tw[0]<0)
 def U(z):return tw*kernel.J0(z)+F*kernel.J1(z)
 def grad(z):return (tw+F*max(z-kernel.h,0.))/kernel.D(z)
 PA=lambda z:kernel.nt(z)*np.dot(grad(z),grad(z))
 ML=lambda z:kernel.nu*np.dot(grad(z),grad(z))
 close(label+' no slip',U(0),np.zeros(2))
 close(label+' matching trace',U(H),UH)
 close(label+' total stress derivative/integral',TH-tw,F*(H-kernel.h))
 close(label+' actual upper shear',kernel.D(H)*grad(H),TH)
 prod=kernel.integral(PA);mol=kernel.integral(ML)
 forcepower=kernel.integral(lambda z:np.dot(F,U(z)) if z>kernel.h else 0.)
 upperpower=float(np.dot(UH,TH))
 close(label+' forced molecular+ turbulent work',upperpower-forcepower,prod+mol)
 p0int=kernel.Pint(H,t0)
 src=kernel.integral(lambda z:PA(z)-kernel.nt(z)*(t0/kernel.D(z))**2)
 close(label+' forced source keeps counterfactual p0',src,prod-p0int)
 close(label+' forced combined mean+k energy',upperpower-forcepower-src,p0int+mol)
 # Wrong use of ground shear at the upper interface is a detectable control.
 omitted=float(np.dot(UH,TH-tw))
 require(label+' ground/upper traction confusion detected',abs(omitted)>1e-3*max(abs(upperpower),1e-30))
 # Exact frozen-input force response.
 delta=np.array([.11,-.07]);t2=UH/J0-B*(F+delta)
 close(label+' force response',t2-tw,-B*delta)
 return dict(case=label,wall_shear=tw.tolist(),interface_shear=TH.tolist(),
             production=prod,molecular_loss=mol,forcing_work=forcepower,
             upper_work=upperpower,net_k_source=src,
             omitted_work_if_ground_shear_used=omitted)


def check_kappa_dependencies(label,kernel,q,u,kappa):
 H=kernel.H
 # Centered finite difference independently checks the derivative with H.
 step=1.e-5*H
 def log_z0(hh):return math.log(hh)-kappa*u*kernel.J0(hh)
 derivative=H*(log_z0(H+step)-log_z0(H-step))/(2*step)
 exact=1-kappa*u*H/kernel.D(H)
 close(label+' H-dependence',derivative,exact,atol=2e-9,rtol=2e-8)
 # Active-branch derivative checks the finite molecular/roughness intercept.
 if H>max(kernel.zv,0.):
  close(label+' finite-intercept diffusivity',kernel.D(H),
        kernel.alpha*math.sqrt(kernel.K)*H+kernel.nu*(1-kernel.alpha*kernel.b))
 # The leading slope is tested independently at a high observation height.
 zlarge=max(H,abs(kernel.zv),kernel.nu/kernel.a)*1e8
 effective=kernel.D(zlarge)/(u*zlarge)
 close(label+' overlap slope',effective,kernel.alpha*math.sqrt(q),atol=3e-8,rtol=0)
 return exact

def main():
 results=[];calibrated={}
 H=2.;z0=.03;u=.4;nu=1.5e-5
 for model,q in ROOTS.items():
  D=.48 if model!='E94_E95' else 1.9-.43*q
  sigma_original=1.3 if model!='E94_E95' else 1.2
  natural=math.sqrt(sigma_original*D/q)
  C=(2/3)/(1.25+1.9*q) if model=='A93' else 1/((6.5 if model=='A95' else 4.)+3*q/math.sqrt(2))
  close(model+' constitutive root',C*q*q,1.)
  for choice,kappa in [('natural',natural),('common041',.41)]:
   sigma=kappa*kappa*q/D
   close(model+' '+choice+' outer balance',kappa*kappa*q/sigma-D,0.)
   if choice=='natural':close(model+' unchanged outer sigma',sigma,sigma_original)
   for variant in ['original','slope','both_tails']:
    alpha=ALPHA if variant=='original' else kappa/math.sqrt(q)
    clw=kappa*q**1.5 if variant=='both_tails' else CLW
    kernel=calibrate_h(q,u,nu,H,z0,kappa,alpha,clw)
    label=model+'_'+choice+'_'+variant
    rec=check_reference(label,kernel,q,u,z0,model,kappa=kappa)
    rec.update(kappa=kappa,alpha=alpha,clw=clw,sigma_epsilon=sigma,
               H_log_derivative=check_kappa_dependencies(label,kernel,q,u,kappa))
    if variant!='original':close(label+' matching leading slope',alpha*math.sqrt(q),kappa)
    if variant=='both_tails':
     close(label+' matching baseline tail',kernel.Ebase(H),u**3/(kappa*H))
     close(label+' matching wall cmu',alpha/clw,1/q**2)
    results.append(rec)
    if choice=='common041' and variant=='original':calibrated[model]=kernel
 q=ROOTS['E94_E95'];K=q*u*u
 smooth=Kernel(K,nu,.05,0.)
 smooth_z0=smooth.H*math.exp(-KAPPA*u*smooth.J0(smooth.H))
 smooth_endpoint=calibrate_h(q,u,nu,smooth.H,smooth_z0)
 close('admissible smooth calibration endpoint',smooth_endpoint.h,0.,atol=0.,rtol=0.)
 require('smooth test contains molecular sublayer',0<smooth.zv<smooth.H)
 results.append(check_reference('smooth_molecular_sublayer',smooth,q,u,model='E94_E95'))
 molecular=Kernel(K,nu,5e-5,0.)
 require('all molecular test lies below viscosity intercept',molecular.zv>molecular.H)
 results.append(check_reference('all_molecular_EP_rejected',molecular,q,u,expect_gate=False))
 results.append(check_forced('rough_forced_reversed',calibrated['E94_E95'],q,u))
 results.append(check_forced('smooth_forced_reversed',smooth,q,u))
 print(json.dumps(dict(status='PASS',assertions=COUNT,
  scope='Analytical identities and independent quadrature; no PDE solver or terrain validation.',
  results=results),indent=2))

if __name__=='__main__':main()
