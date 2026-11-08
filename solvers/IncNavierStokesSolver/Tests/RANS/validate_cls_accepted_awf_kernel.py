#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent Linux oracle for the calibrated budget AWF appendix.

The oracle derives the outer constants, metric-weighted production primitive,
cutoff and nonconstant reference independently of the native quadrature. Its
factor-three Gauss panels differ from the kernel's panels. A real-dilogarithm
primitive provides another independent flat rough-reference evaluation.
Only the scalar curved column is checked: no curved PIMPLE support is claimed.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

THREADS = {name: "1" for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS", "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS")}
os.environ.update(THREADS)
sys.dont_write_bytecode = True
import numpy as np
import scipy
from scipy.optimize import brentq
from scipy.special import eval_jacobi, eval_legendre, spence

LD = np.longdouble
SPEC = "docs/literature/accepted_awf_appendix.tex"
SPEC_SHA = "9998548d33bd2ca856bb6cc85551c82b86d3b499fd54ae8dc9ba26453450fdbf"
NOTE = "docs/literature/cls_neutral_abl_equilibrium_awf_note_revised.tex"
NOTE_SHA = "d35e437c2f2f8b25a46f7a3f285c6f1ad8b48d7b554aa1d161315f59eeb6de0f"
CHECKS = 0
ERRORS = {}


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def near(actual, expected, name, scale=1., tolerance=2e-10):
    a, b = np.asarray(actual, dtype=LD), np.asarray(expected, dtype=LD)
    require(a.shape == b.shape, name+": shape mismatch")
    require(np.all(np.isfinite(a)) and np.all(np.isfinite(b)), name+": nonfinite")
    denominator = max(LD(scale), np.max(np.abs(b)), LD("1e-300"))
    error = float(np.max(np.abs(a-b))/denominator)
    ERRORS[name] = max(error, ERRORS.get(name, 0.))
    require(error <= tolerance, f"{name}: scaled error {error:.6g} > {tolerance}")
    return error


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for data in iter(lambda: stream.read(1048576), b""):
            h.update(data)
    return h.hexdigest()


def constants():
    def C(q):
        return .3/(1+.35*q**1.5)*(-math.expm1(-.36*math.exp(.75*q)))
    q = brentq(lambda x: C(x)*x*x-1, 2, 4, xtol=5e-15)
    def slope(kappa):
        r = kappa*q**1.5/2.5
        return q*kappa*kappa/1.3+.83*max(r-1, 0)*r*r-.48
    return LD(q), LD(brentq(slope, .1, 1, xtol=5e-15))


class BudgetInapplicable(ValueError):
    pass


class Reference:
    """Closed J-weighted source primitives, independently integrated k."""
    def __init__(self, H=2., h=.1, nu=1.5e-5, ustar=.4,
                 kappa1=0., kappa2=0., sampled_k=None):
        self.H, self.h, self.nu = map(LD, (H, h, nu))
        self.q, self.kappa = constants()
        self.K = self.q*LD(ustar)**2 if sampled_k is None else LD(sampled_k)
        require(self.K > 0 and 0 <= self.h < self.H and self.nu > 0,
                "oracle geometry/scalar domain")
        self.k1, self.k2 = LD(kappa1), LD(kappa2)
        require(min(1, 1+self.H*self.k1) > 0 and min(1, 1+self.H*self.k2) > 0,
                "oracle principal stretch must be positive")
        self.s1, self.s2 = self.k1+self.k2, self.k1*self.k2
        self.cMu, self.cLength = self.q**-2, self.kappa*self.q**LD("1.5")
        self.a = self.cMu*self.cLength*np.sqrt(self.K)
        self.t0 = self.K/self.q
        self.ustar = np.sqrt(self.t0)
        self.A = self.ustar**3/self.kappa
        if self.h == 0:
            b = LD("10.7")
        else:
            x = self.h*np.sqrt(self.K)/(70*self.nu)
            m = max(LD(".5")-LD(".4")*x**LD(".7"),
                    1-LD(".79")*x**LD("-.28"))
            b = LD("10.7")*(1-x**m)
        self.zv = self.nu*b/np.sqrt(self.K)
        self.d = self.nu-self.a*self.zv
        self.lo = max(LD(0), self.zv)
        self.w0 = self.d+self.a*self.lo
        self.c0 = 1-self.s1*self.d/self.a+self.s2*(self.d/self.a)**2
        self.c1 = self.s1/self.a-2*self.s2*self.d/self.a**2
        self.c2 = self.s2/self.a**2
        self.J0, self.J1, self.J2 = self.integrals(self.H)
        self.UH = self.t0*self.J0
        self.z0 = self.H*np.exp(-self.kappa*self.UH/self.ustar)
        self.Pj = self.production_cumulative(self.H)
        self.minimum = self.A*self.volume(self.H)/self.H
        if self.Pj < self.minimum:
            raise BudgetInapplicable("reference budget below metric minimum")
        self.Lambda = self.Pj/self.A
        high = math.log(float(self.H)); low = high-1
        def residual(logc):
            return float(self.dissipation_cumulative(self.H, np.exp(LD(logc)))-self.Pj)
        while residual(low) < 0:
            low -= 2
        self.cutoff = np.exp(LD(brentq(residual, low, high, xtol=5e-14, rtol=1e-14)))

    def J(self, z):
        z = np.asarray(z, dtype=LD)
        return (1+self.k1*z)*(1+self.k2*z)

    def volume(self, z):
        z = np.asarray(z, dtype=LD)
        return z+self.s1*z*z/2+self.s2*z*z*z/3

    def nt(self, z):
        return self.a*np.maximum(np.asarray(z, dtype=LD)-self.zv, 0)

    def D(self, z):
        return self.nu+self.nt(z)

    def integrals(self, z):
        z = np.asarray(z, dtype=LD)
        molecular = np.minimum(z, self.lo)
        turbulent = np.maximum(z, self.lo)
        delta = self.a*(turbulent-self.lo)
        log = np.log1p(delta/self.w0)
        j0 = molecular/self.nu+log/self.a
        j2 = molecular/self.nu**2+delta/(self.a*self.w0*(self.w0+delta))
        lower = max(self.h, self.lo)
        wlower = self.d+self.a*lower
        length = np.maximum(z-lower, 0)
        j1 = np.maximum(molecular-self.h, 0)**2/(2*self.nu)
        j1 += length/self.a-(self.d+self.a*self.h)/self.a**2*np.log1p(self.a*length/wlower)
        return j0, j1, j2

    def p0(self, z):
        return self.t0**2*self.nt(z)/self.D(z)**2

    def production_cumulative(self, z):
        z = np.maximum(np.asarray(z, dtype=LD), self.lo)
        delta = self.a*(z-self.lo); w = self.w0+delta
        return self.t0**2/self.a*(
            (self.c0-self.nu*self.c1)*np.log1p(delta/self.w0)
            +(self.c1-self.nu*self.c2)*delta
            +self.c2*delta*(w+self.w0)/2
            -self.nu*self.c0*delta/(w*self.w0))

    def dissipation_cumulative(self, z, cutoff=None):
        cutoff = self.cutoff if cutoff is None else cutoff
        z = np.asarray(z, dtype=LD)
        below, above = np.minimum(z, cutoff), np.maximum(z, cutoff)
        return self.A*(self.volume(below)/cutoff+np.log(above/cutoff)
            +self.s1*(above-cutoff)+self.s2*(above-cutoff)*(above+cutoff)/2)

    def explicit_budget(self, cutoff):
        return self.A*(1+np.log(self.H/cutoff)+self.s1*(self.H-cutoff/2)
            +self.s2*(self.H**2/2-cutoff**2/6))

    def E(self, z):
        return self.A/np.maximum(np.asarray(z, dtype=LD), self.cutoff)

    def Q(self, z):
        return self.dissipation_cumulative(z)-self.production_cumulative(z)

    def gradient(self, z):
        return self.Q(z)/(self.J(z)*self.D(z))

    def quadrature(self, lower=0., upper=None, order=64, extra=()):
        lo, hi = LD(lower), self.H if upper is None else LD(upper)
        if hi == lo:
            return np.array([], dtype=LD), np.array([], dtype=LD)
        breaks = {lo, hi}
        for z in (self.zv, self.h, self.cutoff, *extra):
            if lo < z < hi:breaks.add(LD(z))
        z = min(self.nu/self.a, self.cutoff)
        while z < hi:
            if lo < z:breaks.add(z)
            z *= 3
        x, w = (np.asarray(x, dtype=LD) for x in np.polynomial.legendre.leggauss(order))
        points, weights = [], []
        panels = sorted(breaks)
        for a, b in zip(panels[:-1], panels[1:]):
            half = (b-a)/2
            points.extend(a+half*(1+x));weights.extend(half*w)
        return np.asarray(points, dtype=LD), np.asarray(weights, dtype=LD)

    def k_integral(self, z, order=64):
        r, w = self.quadrature(z, order=order)
        return self.K-np.sum(w*self.gradient(r), dtype=LD)

    def k_closed(self, z):
        require(self.zv <= 0 and self.k1 == 0 and self.k2 == 0,
                "closed primitive requires a flat rough branch")
        beta, zc = self.nu/self.d, self.a*self.cutoff/self.d
        def low(x):
            l = np.log1p(x)
            return (x-l)/zc-l*l/2+beta*(l+1/(1+x)-1)
        def high(x):
            l = np.log1p(x)
            # scipy's real Spence function is Li_2(1-x).
            dilog = LD(spence(float(1+x)))
            return (1+np.log(x/zc)+beta)*l+dilog-l*l/2+beta/(1+x)
        def joined(x):
            return low(x) if x <= zc else low(zc)+high(x)-high(zc)
        return self.K+self.A/self.a*(joined(self.a*LD(z)/self.d)-joined(self.a*self.H/self.d))

    def summary(self):
        return {k: float(getattr(self, k)) for k in
            ("H", "h", "nu", "K", "q", "kappa", "cMu", "cLength", "a", "zv", "d",
             "J0", "J1", "J2", "UH", "t0", "A", "Pj", "minimum", "Lambda", "cutoff", "z0")}


def weak_moments(p, order=64):
    r, w = p.quadrature(order=order)
    x = np.asarray(2*r/p.H-1, dtype=float)
    moments = []
    for degree in range(6):
        phi = np.asarray(eval_legendre(degree, x), dtype=LD)
        derivative = np.zeros_like(r) if degree == 0 else np.asarray(
            (degree+1)/2*eval_jacobi(degree-1, 1, 1, x), dtype=LD)*2/p.H
        diffusion = -np.sum(w*p.Q(r)*derivative, dtype=LD)
        source = np.sum(w*p.J(r)*(p.p0(r)-p.E(r))*phi, dtype=LD)
        moments.append(dict(degree=degree, diffusion=float(diffusion), source=float(source),
                            residual=float(diffusion+source)))
    return moments


def native_cubic_partition_moments(p):
    """Continuum weak balance on the pilot's assembled modified-CG rows."""
    panels = sorted({LD(0), p.H, *(LD(x) for x in (.01, .1, .25, .5, 1., 2.) if 0 < x < p.H)})
    n = len(panels);diffusion = np.zeros(n+2*(n-1), dtype=LD);source = diffusion.copy()
    for element, (a, b) in enumerate(zip(panels[:-1], panels[1:])):
        r, w = p.quadrature(a, b)
        x = 2*(r-a)/(b-a)-1
        phi = np.asarray([(1-x)/2, (1+x)/2, (1-x*x)/4, x*(1-x*x)/2])
        derivative = np.asarray([-np.ones_like(x)/2, np.ones_like(x)/2,
                                  -x/2, (1-3*x*x)/2])*2/(b-a)
        ids = [element, element+1, n+2*element, n+2*element+1]
        for row, i in enumerate(ids):
            diffusion[i] -= np.sum(w*p.Q(r)*derivative[row], dtype=LD)
            source[i] += np.sum(w*p.J(r)*(p.p0(r)-p.E(r))*phi[row], dtype=LD)
    return dict(panels=[float(x) for x in panels], diffusion=[float(x) for x in diffusion],
        source=[float(x) for x in source], residual=[float(x) for x in diffusion+source],
        physical_flux_convention="Q_J=J*D*k_prime; upward diffusive flux=-Q_J/J")


def independent_checks(p, name):
    scale = max(float(p.Pj), float(p.minimum))
    near(p.a, p.kappa*p.ustar, name+"_calibrated_slope", float(p.a), 2e-14)
    near(p.K**LD("1.5")/p.cLength, p.A, name+"_calibrated_tail", float(p.A), 2e-14)
    near(p.t0**2/p.a, p.A, name+"_reference_amplitude", float(p.A), 2e-14)
    r, w = p.quadrature()
    near(np.sum(w/p.D(r)), p.J0, name+"_J0_integral", float(p.J0), 2e-12)
    near(np.sum(w*np.maximum(r-p.h, 0)/p.D(r)), p.J1, name+"_J1_integral", float(p.J1), 2e-12)
    near(np.sum(w/p.D(r)**2), p.J2, name+"_J2_integral", float(p.J2), 2e-12)
    near(np.sum(w*p.J(r)*p.p0(r)), p.Pj, name+"_production_primitive", scale, 2e-12)
    near(np.sum(w*p.J(r)*p.E(r)), p.Pj, name+"_dissipation_budget", scale, 2e-12)
    near(p.explicit_budget(p.cutoff), p.Pj, name+"_explicit_metric_budget", scale, 2e-12)
    near(p.Q(0), 0, name+"_ground_flux", scale, 2e-12)
    near(p.Q(p.H), 0, name+"_matching_flux", scale, 2e-12)
    near(p.k_integral(p.H), p.K, name+"_endpoint_K", float(p.K), 2e-12)
    require(p.explicit_budget(p.cutoff/2) > p.Pj,
            name+": budget monotonicity below cutoff")
    require(p.explicit_budget(min(p.H, 2*p.cutoff)) < p.Pj,
            name+": budget monotonicity above cutoff")
    c = p.cutoff;step = c*LD("1e-4")
    numerical = (p.explicit_budget(c+step)-p.explicit_budget(c-step))/(2*step)
    exact = -p.A*p.volume(c)/c**2
    near(numerical, exact, name+"_metric_budget_derivative", abs(float(exact)), 2e-8)
    moments = weak_moments(p)
    refined = weak_moments(p, order=96)
    for a, b in zip(moments, refined):
        near(a['residual'], 0, name+"_all_Legendre_modes", scale, 2e-11)
        near([a['diffusion'],a['source']], [b['diffusion'],b['source']],
             name+"_oracle_quad_refinement", scale, 2e-12)
    rows = native_cubic_partition_moments(p)
    near(rows['residual'], np.zeros(len(rows['residual'])), name+"_all_modified_CG_rows", scale, 2e-11)
    k0 = p.k_integral(0)
    near(k0, p.k_integral(0, order=96), name+"_k_integral_refinement", float(p.K), 2e-12)
    closed_errors = []
    if p.k1 == 0 and p.k2 == 0:
        near(p.cutoff, p.H*np.exp(1-p.Pj/p.A), name+"_flat_exponential", float(p.cutoff), 2e-12)
        near(p.t0*p.UH, p.nu*p.t0**2*p.J2+p.Pj, name+"_mechanical_energy", scale, 2e-12)
        if p.zv <= 0:
            for z in (0., float(p.cutoff/2), float(p.cutoff), .1, float(p.H/2), float(p.H)):
                if z > p.H:continue
                err = near(p.k_integral(z), p.k_closed(z), name+"_closed_k_primitive", float(p.K), 2e-11)
                closed_errors.append(err)
            if p.d >= 2*p.nu:
                grid = np.geomspace(float(p.H)*1e-8, float(p.H)*(1-1e-6), 160)
                require(np.all(p.Q(grid) < 0), name+": rough reference flux sign")
                require(k0 > p.K, name+": elevated rough reference ground k")
    # Counterfactual source comparison on identical fields: actual production
    # differs from p0, while both policies retain the same integral by budget.
    represented_production = p.p0(r)*(1+LD(".07")*r/p.H)**2
    accepted_source = represented_production-p.E(r)
    ep_source = represented_production-p.p0(r)
    counterfactual = []
    for degree in range(6):
        phi = np.asarray(eval_legendre(degree,np.asarray(2*r/p.H-1,dtype=float)),dtype=LD)
        counterfactual.append(dict(degree=degree,
            accepted=float(np.sum(w*p.J(r)*accepted_source*phi,dtype=LD)),
            calibrated_ep=float(np.sum(w*p.J(r)*ep_source*phi,dtype=LD))))
    near(counterfactual[0]['accepted'],counterfactual[0]['calibrated_ep'],
         name+"_same_state_source_integral",scale,2e-12)
    require(max(abs(v['accepted']-v['calibrated_ep']) for v in counterfactual[1:]) > 1e-5*scale,
            name+": higher-mode source-policy distinction disappeared")
    return dict(name=name,reference=p.summary(),ground_k=float(k0),
        k_at_0_1=float(p.k_integral(.1)) if p.H>=.1 else None,
        weak_modes=moments,modified_CG=rows,
        closed_profile_max_scaled_difference=max(closed_errors,default=None),
        source_counterfactual=counterfactual,source_counterfactual_scope="identical retained reference coefficients/represented production; no separate evolution")


def compare_native(native, name):
    g, m = native['geometry'], native['model']; ref = native['reference']
    p = Reference(H=g['height'], h=g['roughness'],nu=g['molecular_nu'],
        kappa1=g['kappa1'], kappa2=g['kappa2'], sampled_k=ref['K'])
    near([m['q'],m['kappa']],[p.q,p.kappa],name+"_native_outer_constants",1,3e-14)
    near([native['wall'][k] for k in ('cMu','cLength','smoothIntercept')],
         [p.cMu,p.cLength,10.7],name+"_native_wall_parameters",1,3e-14)
    require(g['sample_height']==g['height'] and g['tail_start']==g['height'],
            name+": endpoint sampler/tail not enforced")
    require(m['sigma_k']==1, name+": selected primary scalar diffusivity")
    expected={"friction_velocity":p.ustar,"traction":p.t0,"matching_velocity":p.UH,
        "numerator":p.A,"production_budget":p.Pj,"minimum_budget":p.minimum,
        "lambda":p.Lambda,"cutoff":p.cutoff,"budget_residual":0.,"J0":p.J0,"J1":p.J1,
        "J2":p.J2,"a":p.a,"zv":p.zv,"d":p.d,"aerodynamic_z0":p.z0,"ground_k":p.k_integral(0)}
    for key,value in expected.items():
        scale=max(abs(float(value)),1e-14)
        if key=='budget_residual':scale=max(float(p.Pj),float(p.minimum))
        near(ref[key],value,name+"_native_"+key,scale,
             tolerance=2e-12 if key=='budget_residual' else 2e-10)
    scale = max(float(p.Pj),float(p.minimum))
    for sample in native['samples']:
        r=LD(sample['r']);D=p.D(r);Q=p.Q(r)
        target={'k':p.k_integral(r),'k_gradient':Q/(p.J(r)*D),'Q':Q,
            'D':D,'nu_t':p.nt(r),'p0':p.p0(r),'E':p.E(r),'source':p.p0(r)-p.E(r),
            'redistribution':0.,'J':p.J(r),'U':p.t0*p.integrals(r)[0],
            'U_gradient':p.t0/D}
        for key,value in target.items():
            unit_scale = p.K if key=='k' else (p.K/p.H if key=='k_gradient' else scale/p.H)
            if key in ('Q',):unit_scale=scale
            if key in ('D','nu_t'):unit_scale=max(p.nu,p.a*p.H)
            if key in ('J',):unit_scale=1
            if key in ('U',):unit_scale=p.UH
            if key in ('U_gradient',):unit_scale=p.t0/p.D(0)
            near(sample[key],value,name+"_native_sample_"+key,float(unit_scale))
        require(sample['redistribution']==0.,name+": prescribed redistribution must be exactly zero")
    quad=native['quadrature'];r=np.array([x['r'] for x in quad],dtype=LD)
    w=np.array([x['weight'] for x in quad],dtype=LD);J=np.array([x['J'] for x in quad],dtype=LD)
    require(np.all(w>0) and np.all(J>0),name+": positive native quadrature")
    near(J,p.J(r),name+"_native_metric",1,2e-13)
    near(np.sum(w*J),p.volume(p.H),name+"_native_volume",float(p.H),2e-12)
    budgets=native['cutoff_budget_samples']
    require(len(budgets)>=3,name+": incomplete native metric-budget samples")
    for sample in budgets:
        c=LD(sample['cutoff']);r,w=p.quadrature(extra=(c,))
        quadrature_value=np.sum(w*p.J(r)*p.A/np.maximum(r,c),dtype=LD)
        near(sample['integral'],p.explicit_budget(c),name+"_native_explicit_budget",scale)
        near(sample['integral'],quadrature_value,name+"_native_numerical_metric_budget",scale)
    require(all(a['cutoff']<b['cutoff'] and a['integral']>b['integral']
        for a,b in zip(budgets[:-1],budgets[1:])),name+": native budget not strictly decreasing")
    if len(budgets)>=5:
        a,c,b=budgets[:3]
        slope=(b['integral']-a['integral'])/(b['cutoff']-a['cutoff'])
        derivative=-p.A*p.volume(LD(c['cutoff']))/LD(c['cutoff'])**2
        near(slope,derivative,name+"_native_metric_budget_derivative",abs(float(derivative)),2e-8)
    moments=weak_moments(p)
    require([int(x['degree']) for x in native['weak_modes']]==list(range(6)),name+": incomplete weak modes")
    for actual,target in zip(native['weak_modes'],moments):
        near([actual[k] for k in ('diffusion','source','residual')],
             [target[k] for k in ('diffusion','source','residual')],name+"_native_all_weak_modes",scale)
    forced=native['forced'];u=np.array(forced['matching_velocity'],dtype=LD);f=np.array(forced['forcing'],dtype=LD)
    traction=(u-f*p.J1)/p.J0
    near(forced['traction'],traction,name+"_forced_traction",float(p.t0))
    r,w=p.quadrature()
    grad=(traction[:,None]+f[:,None]*np.maximum(r-p.h,0))/p.D(r)
    actual_production=p.nt(r)*np.sum(grad*grad,axis=0)
    near(forced['mean_production'],np.sum(w*actual_production)/p.H,
         name+"_forced_mean_production",scale/float(p.H))
    counterfactual=[]
    for actual in native['weak_modes']:
        phi=np.asarray(eval_legendre(actual['degree'],np.asarray(2*r/p.H-1,dtype=float)),dtype=LD)
        accepted=np.sum(w*p.J(r)*(actual_production-p.E(r))*phi,dtype=LD)
        ep=np.sum(w*p.J(r)*(actual_production-p.p0(r))*phi,dtype=LD)
        near([actual['forced_accepted_source'],actual['forced_calibrated_ep_source']],
             [accepted,ep],name+"_forced_same_state_source_modes",scale)
        counterfactual.append(dict(degree=actual['degree'],accepted=float(accepted),calibrated_ep=float(ep)))
    near(counterfactual[0]['accepted'],counterfactual[0]['calibrated_ep'],
         name+"_forced_source_policy_integral",scale,2e-12)
    for sample in forced['samples']:
        z=LD(sample['r']);j0,j1,_=p.integrals(z)
        velocity=traction*j0+f*j1
        gradient=(traction+f*max(z-p.h,LD(0)))/p.D(z)
        near(sample.get('u',sample.get('velocity')),velocity,name+"_forced_velocity",float(p.UH))
        near(sample['gradient'],gradient,name+"_forced_gradient",float(p.t0/p.D(0)))
    return dict(name=name,reference=p.summary(),native_ground_k=ref['ground_k'],
        weak_modes=native['weak_modes'],cutoff_iterations=ref.get('cutoff_iterations'),native_quadrature_points=len(quad),
        native_forced_source_counterfactual=counterfactual,
        independent_profile_quadrature="GL64 factor-three panels, closed metric source primitives",
        curved_scope="frozen normal column only; not a terrain equilibrium or curved PIMPLE evolution")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--self-check',action='store_true',help='independent analytical checks, no native executable')
    parser.add_argument('--quick',action='store_true',help='native order16 only; independent refinement still checked')
    args=parser.parse_args()
    require(sys.platform.startswith('linux'),'Authoritative numerical validation requires Linux')
    require(args.probe is not None or args.self_check,'supply --probe or explicit --self-check')
    root=Path(__file__).resolve().parents[4]
    require(digest(root/SPEC)==SPEC_SHA and digest(root/NOTE)==NOTE_SHA,'supplied specification changed')
    args.output.mkdir(parents=True,exist_ok=False)
    report=dict(schema=1,status='running',captured_UTC=datetime.now(timezone.utc).isoformat(),
        host=platform.node(),python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        specification_sha256={SPEC:SPEC_SHA,NOTE:NOTE_SHA},thread_environment=THREADS,ranks=1,
        source_revision=subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip(),
        scope='calibrated flat analytical reference and metric-weighted frozen scalar columns; no native time evolution or curved PIMPLE claim',
        commands=[],independent_cases=[],native_cases=[])
    cases=[('flat_H2',dict(H=2.)),('flat_H05',dict(H=.5)),('smooth_H2',dict(H=2.,h=0.)),
        ('curved_positive',dict(H=2.,kappa1=.1)),('curved_negative',dict(H=2.,kappa1=-.1)),
        ('curved_mixed',dict(H=2.,kappa1=.12,kappa2=-.06)),
        ('current_K',dict(H=.5,sampled_k=.7))]
    native_lookup={}
    try:
        for name,parameters in cases:
            p=Reference(**parameters);result=independent_checks(p,name);report['independent_cases'].append(result)
            if name=='smooth_H2':
                require(result['ground_k']<0,'smooth expected reference qualification failure disappeared')
                near(result['ground_k'],-.30407,'rounded_smooth_ground_target',1.,1e-5)
            if name in ('flat_H2','flat_H05'):
                table={'flat_H2':[5.909280151,.01475559464,5.379706230,.005342939599,1.022910439,.5447982810],
                       'flat_H05':[4.531094244,.01463643334,4.128270817,.005300459546,.9942980291,.5307430752]}[name]
                for label,actual,target in zip(('Lambda','cutoff','UH','z0','ground_k','k_0.1'),
                    (p.Lambda,p.cutoff,p.UH,p.z0,result['ground_k'],result['k_at_0_1']),table):
                    near(actual,target,name+'_rounded_appendix_'+label,max(abs(target),1e-8),2e-9)
            if args.probe is None:continue
            for order in ((16,) if args.quick else (16,32)):
                filename=f'{name}-q{order}.json'
                command=[str(args.probe.resolve()),'--output',str((args.output/filename).resolve()),
                    '--height',str(parameters.get('H',2.)),'--roughness',str(parameters.get('h',.1)),
                    '--kappa1',str(parameters.get('kappa1',0.)),'--kappa2',str(parameters.get('kappa2',0.)),
                    '--order',str(order)]
                if 'sampled_k' in parameters:command+=['--sampled-k',str(parameters['sampled_k'])]
                begin=time.monotonic();completed=subprocess.run(command,cwd=args.output,capture_output=True,text=True,timeout=300,env=os.environ)
                (args.output/(filename+'.log')).write_text(completed.stdout+completed.stderr)
                report['commands'].append(dict(command=command,exit_code=completed.returncode,seconds=time.monotonic()-begin))
                require(completed.returncode==0,name+': native kernel failure '+completed.stderr[-500:])
                raw=json.loads((args.output/filename).read_text());require(raw['status']=='ok',name+': native not valid')
                report['native_cases'].append(compare_native(raw,name+f'_q{order}'));native_lookup[name,order]=raw
        if args.probe is not None:
            for name,flags in [('inadmissible_budget',['--height','2','--roughness','0','--nu','.05']),
                ('invalid_stretch',['--height','2','--kappa1','-.6']),
                ('positive_J_invalid_stretches',['--height','2','--kappa1','-.6','--kappa2','-.6'])]:
                output=args.output/(name+'.json');command=[str(args.probe.resolve()),'--output',str(output.resolve()),*flags]
                begin=time.monotonic();completed=subprocess.run(command,cwd=args.output,capture_output=True,text=True,timeout=90,env=os.environ)
                (args.output/(name+'.log')).write_text(completed.stdout+completed.stderr)
                report['commands'].append(dict(command=command,exit_code=completed.returncode,seconds=time.monotonic()-begin))
                require(completed.returncode==2,name+': expected classified applicability rejection')
                raw=json.loads(output.read_text());require(raw['status']=='inapplicable',name+': wrong classification')
                report.setdefault('expected_rejections',{})[name]=raw
            if not args.quick:
                for name,_ in cases:
                    a,b=native_lookup[name,16],native_lookup[name,32]
                    near(a['reference']['cutoff'],b['reference']['cutoff'],name+'_native_quad_cutoff',abs(b['reference']['cutoff']),2e-11)
                    near([[v[k] for k in ('diffusion','source')] for v in a['weak_modes']],
                         [[v[k] for k in ('diffusion','source')] for v in b['weak_modes']],
                         name+'_native_quad_weak',b['reference']['production_budget'],2e-10)
                    near(a['reference']['ground_k'],b['reference']['ground_k'],name+'_native_quad_ground_k',abs(b['reference']['K']),2e-10)
        # The metric minimum is a separate condition from reference positivity.
        try:
            Reference(H=2,h=0,nu=.05)
        except BudgetInapplicable:
            report['independent_budget_gate']='rejected_below_minimum_without_cutoff_clipping'
        else:raise RuntimeError('independent budget-inapplicability control failed')
        report['status']='passed'
        report['smooth_reference']='cutoff exists, but stationary-reference positivity qualification fails; no clipping'
        report['curved_limitations']=[
            'positive local principal stretches do not establish global normal-map injectivity',
            'frozen normal-column balance omits actual terrain production, advection and tangential diffusion',
            'full mapped momentum/continuity, curvature/normal-motion traction, surface normals and conservative fluxes require separate native support',
            'endpoint K must select the actual matching face and local physical H for future varying-height geometry',
            'no general volume-positivity theorem or curved time-dependent validation follows']
    except Exception as error:
        report['status']='failed';report['failure']=str(error)
        raise
    finally:
        report['checks']=CHECKS;report['scaled_errors']=ERRORS
        report['source_status']=subprocess.check_output(['git','-C',str(root),'status','--porcelain'],text=True).splitlines()
        sources=[SPEC,NOTE,'solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_accepted_awf_kernel.py',
            'library/SolverUtils/Turbulence/EquilibriumAWFProfiles.h','library/SolverUtils/Turbulence/EquilibriumAWFProfiles.cpp',
            'library/SolverUtils/Turbulence/AnalyticalWallFunction.cpp',
            'library/UnitTests/SolverUtils/RANSAcceptedAWFProbe.cpp']
        report['source_sha256']={p:digest(root/p) for p in sources if (root/p).exists()}
        if args.probe is not None:report['executable_sha256']=digest(args.probe)
        report['artifacts']={p.name:digest(p) for p in args.output.iterdir() if p.is_file()}
        (args.output/'validation.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status=report['status'],checks=CHECKS,report=str(args.output/'validation.json'))))


if __name__=='__main__':main()
