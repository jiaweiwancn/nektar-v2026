# SPDX-License-Identifier: MIT
"""Independent exact fields and strong sources for the smooth hill fixture."""

import numpy as np


class Manufactured:
    def __init__(self, hill, reference):
        self.hill, self.ref = hill, reference
        self.q = np.pi / hill.L

    def graph(self, a):
        c, q = self.hill.height/2, self.q
        return (c*(1+np.cos(q*a)), -c*q*np.sin(q*a),
                -c*q*q*np.cos(q*a), c*q**3*np.sin(q*a))

    def to_tube(self, e, a, rho):
        if not e["epsilon"]:
            return a, rho
        x, z = self.hill.metrics(e, a, rho)["position"].T
        b = np.array(a,dtype=float,copy=True)
        for _ in range(8):
            f, p, pp, _ = self.graph(b)
            b -= (b - x + (f-z)*p) / (1+p*p+(f-z)*pp)
        f, p, _, _ = self.graph(b)
        if np.max(np.abs(b - x + (f-z)*p)) > 1e-13:
            raise ValueError("independent manufactured foot inversion failed")
        return b, (z-f-(x-b)*p)/np.hypot(1, p)-self.hill.d

    def K(self, a, t):
        return self.ref.Kin*(1+.012*np.sin(.6*t)+.008*a*np.cos(.8*t))**2

    def local(self, a, r, t):
        H, d, R = self.hill.H, self.hill.d, .5
        _, p, pp, ppp = self.graph(a)
        ell = np.sqrt(1+p*p)
        turn = pp/(1+p*p)
        turn_a = ppp/(1+p*p) - 2*p*pp**2/(1+p*p)**2
        B = ell-(d+r)*turn
        Ba = p*pp/ell-(d+r)*turn_a
        tangent = np.column_stack((1/ell, p/ell))
        normal = np.column_stack((-p/ell, 1/ell))
        basis = np.stack((tangent, normal), axis=2)
        S, V = .02*(1+.2*np.sin(.9*t)), .07*(1+.4*np.cos(.85*t))
        X, Xa, Xaa = a*a*(1-a), 2*a-3*a*a, 2-6*a
        Y, Yr, Yrr = r*r*(R-r), 2*R*r-3*r*r, 2*R-6*r
        u, w = S*r+V*X*Yr, -V*Xa*Y/B
        ua, ur = V*Xa*Yr, S+V*X*Yrr
        wa = -V*Xaa*Y/B-w*Ba/B
        wr = -V*Xa*Yr/B+w*turn/B
        G = np.stack(((ua-turn*w)/B, ur, (wa+turn*u)/B, wr), axis=1).reshape(-1,2,2)
        K = self.K(a,t)
        k = K*(1+.005*(1+.2*np.sin(.7*t))*(r*r-H*H)/(H*R))
        A = K**1.5/self.ref.cl
        q = (r-H)/R
        epsilon = A/H*(1+.03*np.cos(.5*t)*q+.015*q*q)
        Pi = .0002*(1+.3*np.sin(.8*t))*((a-.5)**2+r*r/R**2+.15*a*r/R)
        return dict(u=np.column_stack((u,w)), G=G, turn=turn, basis=basis, B=B,
                    K=K, k=k, epsilon=epsilon, Pi=Pi)

    @staticmethod
    def rotate(f, key):
        return np.einsum("qij,qj->qi", f["basis"], f[key])

    def law(self, a, r, t, outer):
        f = self.local(a,r,t)
        ref = self.ref
        s = (r/ref.H-.4)/.35
        chi = np.where(s.real<=0,0.,np.where(s.real>=1,1.,s**3*(10-15*s+6*s*s)))
        k, K, G = f["k"], f["K"], f["G"]
        if outer:
            E = f["epsilon"]
            aux = np.zeros_like(E)
        else:
            xi = ref.h*np.sqrt(K)/(70*ref.nu)
            m0, m1 = .5-.4*xi**.7, 1-.79*xi**(-.28)
            zv = 10.7*ref.nu*(1-xi**np.where(m0.real>m1.real,m0,m1))/np.sqrt(K)
            aux = ref.kappa*np.sqrt(K/ref.qe)*(r-zv)
            cutoff = ref.H*np.exp(np.polynomial.chebyshev.chebval(
                (2*K-ref.lo-ref.hi)/(ref.hi-ref.lo), ref.coefficients))
            E = K**1.5/ref.cl/np.where(r.real>=cutoff.real,r,cutoff)
        Dk, De = ref.nu+ref.cd*k*k/E, ref.nu+ref.cd*k*k/(ref.se*E)
        M = G+G.transpose(0,2,1)
        S = (.5*k/E)[:,None,None]*M
        O = (.5*k/E)[:,None,None]*(G-G.transpose(0,2,1))
        IS = np.trace(S@S,axis1=1,axis2=2)
        IO = np.trace(O@O,axis1=1,axis2=2)
        N = np.full_like(IS,1.8)
        for _ in range(12):
            N -= ((N-1.8)*(N*N-2*IO)-2.7*IS*N)/(N*N-2*IO+2*N*(N-1.8)-2.7*IS)
        # The exact manufactured velocity is solenoidal. Its planar
        # trace-free WJ reduction is independent of the native 3D tangent.
        anisotropy = -1.2*(N[:,None,None]*S+S@O-O@S)/(N*N-2*IO)[:,None,None]
        auxiliary = np.zeros_like(G)
        auxiliary[:,0,1] = auxiliary[:,1,0] = -aux*M[:,0,1]/k
        anisotropy = (1-chi[:,None,None])*auxiliary+chi[:,None,None]*anisotropy
        R = k[:,None,None]*(2*np.eye(2)/3+anisotropy)
        V = ref.nu*M-k[:,None,None]*anisotropy
        f.update(E=E, Dk=Dk, De=De, stress=V, P=-np.einsum("qij,qij->q",R,G))
        return f

    def gradients(self, a, r, t):
        delta = 1e-26j
        f = self.local(a,r,t)
        da = self.local(a+delta,r,t)
        dr = self.local(a,r+delta,t)
        return {key: np.column_stack((da[key].imag/(delta.imag*f["B"]),dr[key].imag/delta.imag))
                for key in ("k","epsilon","Pi")}

    def sources(self, e, a, rho, t, steady):
        b, r = self.to_tube(e,a,rho)
        outer, delta = bool(e["epsilon"]), 1e-26j
        f = self.law(b,r,t,outer)
        da = self.law(b+delta,r,t,outer)
        dr = self.law(b,r+delta,t,outer)
        dt = self.local(b,r,t+delta)
        grad = self.gradients(b,r,t)
        V = f["stress"]
        Va, Vr = da["stress"].imag/delta.imag, dr["stress"].imag/delta.imag
        turn, B = f["turn"], f["B"]
        divergence = np.column_stack(((Va[:,0,0]-turn*(V[:,1,0]+V[:,0,1]))/B+Vr[:,0,1],
                                      (Va[:,1,0]+turn*(V[:,0,0]-V[:,1,1]))/B+Vr[:,1,1]))
        force = np.einsum("qij,qj->qi",f["G"],f["u"])+grad["Pi"]-divergence
        if not steady:
            force += dt["u"].imag/delta.imag
        force = np.einsum("qij,qj->qi",f["basis"],force)
        mask = rho<self.hill.H
        if np.any(mask):
            h = self.local(a[mask],np.full(sum(mask),self.hill.H),t)
            dh = self.local(a[mask]+delta,np.full(sum(mask),self.hill.H),t)
            sampled = (dh["Pi"].imag-2*dh["K"].imag/3+2*h["u"][:,0]*dh["u"][:,0].imag)/(delta.imag*h["B"])
            force[mask] -= (sampled*(rho[mask]<self.ref.h))[:,None]*h["basis"][:,:,0]
        # Scalar Hessians use analytic directional differences of first
        # complex-step derivatives. The coefficient derivatives select the
        # cutoff branch at the original point.
        step = 8e-5
        fa = [self.gradients(b+s*step,r,t) for s in (-2,-1,1,2)]
        fr = [self.gradients(b,r+s*step,t) for s in (-2,-1,1,2)]
        scalar = {}
        for key, diffusion in (("k","Dk"),("epsilon","De")):
            if key=="epsilon" and not outer:
                continue
            g = grad[key]
            dd_a = (fa[0][key][:,0]-8*fa[1][key][:,0]+8*fa[2][key][:,0]-fa[3][key][:,0])/(12*step)
            dd_r = (fr[0][key][:,1]-8*fr[1][key][:,1]+8*fr[2][key][:,1]-fr[3][key][:,1])/(12*step)
            div = f[diffusion]*(dd_a/B+dd_r-turn/B*g[:,1]) + (
                da[diffusion].imag*g[:,0]/B+dr[diffusion].imag*g[:,1])/delta.imag
            reaction = (f["P"]-f["E"] if key=="k" else
                        (1.44*f["P"]-1.82*f["E"])*f["E"]/f["k"])
            scalar[key] = np.sum(f["u"]*g,axis=1)-div-reaction
            if not steady:
                scalar[key] += dt[key].imag/delta.imag
        return force, scalar

    def normal_flux(self,e,a,rho,t,n,key):
        b,r = self.to_tube(e,a,rho)
        law = self.law(b,r,t,bool(e["epsilon"]))
        g = self.gradients(b,r,t)[key]
        physical = np.einsum("qij,qj->qi",law["basis"],g)
        return -law["Dk" if key=="k" else "De"]*np.sum(n*physical,axis=1)
