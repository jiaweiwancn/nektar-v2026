#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Split numerical integrals independent of the neutral profile primitives."""
import math
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from wj_h5_neutral import QE,CD,KAPPA,NE


def check_reference(c,probe=None):
    x=c.h*math.sqrt(c.K)/(70*c.nu);m=max(.5-.4*x**.7,1-.79*x**(-.28))
    zv=c.nu*10.7*(1-x**m)/math.sqrt(c.K);a=KAPPA/math.sqrt(QE)*math.sqrt(c.K)
    nt=lambda r:a*max(r-zv,0);D=lambda r:c.nu+nt(r)
    P=lambda r:c.u_tau**4*nt(r)/D(r)**2
    def integrate(fn,l=0.,u=None,extra=()):
        u=c.H if u is None else u
        points=sorted(set([l,u]+[v for v in (c.h,zv,*extra)if l<v<u]))
        return sum(quad(fn,p,q,epsabs=1e-16,epsrel=2e-12,limit=100)[0]for p,q in zip(points[:-1],points[1:]))
    j0=integrate(lambda r:1/D(r));j2=integrate(lambda r:1/D(r)**2);budget=integrate(P)
    rc=math.exp(brentq(lambda logc:integrate(lambda r:c.A/max(r,math.exp(logc)),extra=(math.exp(logc),))-budget,
        math.log(c.H)-20,math.log(c.H),xtol=1e-14,rtol=1e-14))
    E=lambda r:c.A/max(r,rc)
    Q=lambda r:integrate(lambda s:E(s)-P(s),u=r,extra=(rc,))
    kw=lambda r:c.K-integrate(lambda s:Q(s)/D(s),l=r,extra=(rc,))
    rgrid=sorted(set([0.,.05*c.H,rc,.5*c.H,.9*c.H,c.H]))
    errors=dict(J0=abs(j0-c.w['J0']),J2=abs(j2-(c.H/(c.d_nu*(c.d_nu+a*c.H)))),
        cutoff=abs(rc-c.w['rc']),momentum_anchor=abs(c.u_tau**2*j0-c.U_outer(c.H)),
        integrated_budget=abs(budget-c.A*(1+math.log(c.H/rc))),
        profile=max(abs(kw(r)-c.k_inner(r))for r in rgrid),
        Q=max(abs(Q(r)-c.w['Q'](r))for r in rgrid))
    weak=[]
    for n in range(6):
        moment=integrate(lambda r:-(0 if n==0 else c.w['Q'](r)*n/c.H*(r/c.H)**(n-1))+
            (P(r)-E(r))*(r/c.H)**n,extra=(rc,))
        weak.append(moment)
    assert max(abs(v)for v in weak)<2e-11*c.A
    assert errors['profile']<2e-11*c.K and errors['Q']<2e-11*c.A
    assert errors['momentum_anchor']<2e-12 and errors['cutoff']<2e-12*c.H
    rT=c.top-c.d0;outer_weak=[]
    for n in (1,2,3):
        phi=lambda r:((r-c.H)/(rT-c.H))**n
        dp=lambda r:n/(rT-c.H)*((r-c.H)/(rT-c.H))**(n-1)
        source=lambda r:(1.44-1.82)*(c.A/r)**2/c.K
        inside_flux=lambda r:-c.u_tau**4/(1.3*r)
        # Admissible epsilon tests vanish at H. The lower computed flux is a
        # reaction for other tests, never a second imposed epsilon condition.
        value=quad(lambda r:inside_flux(r)*dp(r)-source(r)*phi(r),c.H,rT,epsabs=1e-16,epsrel=2e-12)[0]
        value+=c.u_tau**4/(1.3*rT)
        outer_weak.append(value)
    assert max(abs(v)for v in outer_weak)<2e-11*c.A/c.H
    compiled={}
    if probe:
        for r in (c.H,1.,rT):
            G=np.zeros((3,3));G[0,2]=c.u_tau/(KAPPA*r);v=probe.closure(c.K,c.A/r,G)
            assert np.max(abs(v['R']-c.Rc))<2e-12*c.K
            assert abs(v['scalar'][3]-c.A/r)<2e-12*c.A/r
            assert max(abs(v['scalar'][1:3]-KAPPA*c.u_tau*r))<2e-12*KAPPA*c.u_tau*r
        cmd=['W',1,1,0,c.K,c.nu,c.H,c.h,c.H,0,0,c.U_inner(c.H),0,0,0]
        w=probe.query(' '.join(map(str,cmd)))
        assert abs(w['reference'][5]-rc)<2e-12*c.H
        assert abs(w['traction'][0]-c.u_tau**2)<2e-12*c.u_tau**2
        assert np.max(abs(w['profile'][:,1]-np.array([kw(z)for z in w['profile'][:,0]])))<3e-11*c.K
        compiled=dict(cutoff=float(w['reference'][5]),wall_profile_rows=len(w['profile']),factory='compiled shared WJ/Budget dispatch, not native coupled evolution')
    return dict(status='PASS',errors=errors,inner_weak_modes=weak,outer_epsilon_weak_modes=outer_weak,
        independent_rc=rc,independent_kw_ground=kw(0),split_points=sorted(set([0.,c.h,rc,c.H])),compiled=compiled,
        scope='continuum reference, independent quadrature/root/profile; native projection is a separate measurement')
