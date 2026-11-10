#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent WJ-1 oracle against a compiled production factory/adapter.

Reference algebra and the exact rational certificate originate in the reviewed
note's embedded programs. Wall expectations use separate SciPy quadrature split
at physical knots. No production result generates an expected result.
"""
import argparse
import json
import math
from pathlib import Path
import subprocess
import numpy as np
from numpy.testing import assert_allclose
from scipy.integrate import quad
from scipy.optimize import brentq
from wj1_reference_certificate import constitutive, L, Z, pivots

CSTAR = 572/6561
QE = 81/(2*math.sqrt(143))
KAPPA = math.sqrt(1.3*(1.82-1.44)*math.sqrt(CSTAR))
CHECKS = 0
MAX_ERROR = 0.


def near(actual, expected, rtol=3e-12, atol=3e-12):
    global CHECKS, MAX_ERROR
    assert_allclose(actual, expected, rtol=rtol, atol=atol)
    CHECKS += 1
    MAX_ERROR = max(MAX_ERROR, float(np.max(np.abs(np.asarray(actual)-expected))))


def independent(k, e, raw, guess=None):
    theta = np.trace(raw)
    dev = raw-theta*np.eye(3)/3
    B = (k/e)*dev
    a, ceff, N = constitutive(B, guess)
    S, O = (B+B.T)/2, (B-B.T)/2
    s, w, iv = np.trace(S@S), np.trace(O@O), np.trace(S@O@O)
    i5 = np.trace(S@S@O@O)
    Q = (5/6)*(N*N-2*w)*(2*N*N-w)
    beta = np.array([-N*(2*N*N-7*w), -12*iv/N, -2*(N*N-2*w), -6*N, 6])/Q
    R = k*(a+2*np.eye(3)/3)
    nut, nud = ceff*k*k/e, CSTAR*k*k/e
    P, Pd = -np.sum(R*raw), -np.sum(k*a*dev)
    ep, ed = 1.44*e/k*P, 1.82*e*e/k
    return dict(scalar=np.array([ceff, nut, nud, P, Pd, P-e, ep, ed, ep-ed,
        nud, nud/1.3, theta, N, s, w, iv, i5, Q,
        1.8*(N*N-2*w)/N+2*N*(N-1.8)]), R=R, a=a, dev=k*a,
        remainder=k*a+2*nut*S/(k/e), beta=beta)


class Probe:
    def __init__(self, executable):
        self.p = subprocess.Popen([executable, '--stream'], stdin=subprocess.PIPE,
                                  stdout=subprocess.PIPE, text=True)
        self.calls = 0

    def query(self, command):
        self.calls += 1
        self.p.stdin.write(command+'\n')
        self.p.stdin.flush()
        result = json.loads(self.p.stdout.readline())
        assert 'rejected' not in result, command
        return {k: np.asarray(v) for k, v in result.items()}

    def closure(self, k, e, G, direction=None, completion=1):
        fields = ['C' if direction is None else 'D', completion, k, e, 0, *G.ravel()]
        if direction is not None:
            dk, de, dG = direction
            fields.extend([dk, de, *dG.ravel()])
        return self.query(' '.join(map(str, fields)))

    def close(self):
        self.p.stdin.close()
        assert self.p.wait(timeout=10) == 0


def constitutive_checks(probe):
    rng = np.random.default_rng(87131)
    states = [np.zeros((3, 3)), np.array([[0, 2, 0], [-2, 0, 1], [0, -1, 0.]]),
              np.array([[0, 0, QE], [0, 0, 0], [0, 0, 0.]]),
              np.array([[.4, 2, -1], [-3, -.1, .7], [1.1, -.4, -.3]])]
    states += [rng.normal(size=(3, 3)) for _ in range(16)]
    directional_errors = []
    for index, G in enumerate(states):
        G = G-np.trace(G)*np.eye(3)/3
        k, e = ((1., 1.) if index<3 else (1.7, .63))
        value = probe.closure(k, e, G, completion=0)
        expected = independent(k, e, G)
        for key in expected:
            near(value[key], expected[key])
        N, s, w = value['scalar'][[12, 13, 14]]
        residual = (N-1.8)*(N*N-2*w)-2.7*s*N
        near(residual/(1+N**3+2.7*s*N+abs(w*N)), 0, atol=1e-14)
        assert N>=1.8 and w<=0
        near(value['R'], value['R'].T)
        near(np.trace(value['R']), 2*k)
        near(value['scalar'][3], -np.sum(value['R']*G))
        assert value['scalar'][3]>=-1e-13
        if index==0 or index==1:
            near(value['R'], 2*k*np.eye(3)/3)
            near(N, 1.8)
        # Orthogonal changes of coordinates, including reflection.
        O, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        if index%2:
            O[:, 0] *= -1
        rotated = probe.closure(k, e, O@G@O.T, completion=0)
        near(rotated['R'], O@value['R']@O.T)
        near(rotated['scalar'], value['scalar'], rtol=4e-12, atol=1e-10)
        if index==3:
            assert abs(value['scalar'][1]-value['scalar'][2])>.01
        # All independent directions, plus a simultaneous k/epsilon/G action.
        directions = [(1., 0., np.zeros((3, 3))), (0., 1., np.zeros((3, 3)))]
        directions += [(0., 0., H.reshape(3, 3)) for H in np.eye(9)]
        directions += [(.3, -.2, rng.normal(size=(3, 3)))]
        for dk, de, dG in directions:
            action = probe.closure(k, e, G, (dk, de, dG))
            step = 1e-30
            cs = independent(k+step*dk*1j, e+step*de*1j,
                             G.astype(complex)+step*dG*1j, expected['scalar'][12])
            for key in cs:
                near(action[key], np.imag(cs[key])/step, rtol=2e-10, atol=5e-10)
        # A simultaneous direction has a separate multi-step convergence check.
        dk, de, dG = directions[-1]
        exact = probe.closure(k, e, G, (dk, de, dG))
        errors = []
        for step in (1e-3, 3e-4, 1e-4, 3e-5):
            plus = probe.closure(k+step*dk, e+step*de, G+step*dG)
            minus = probe.closure(k-step*dk, e-step*de, G-step*dG)
            # Covariance and physical scalar sources/diffusivities.
            fd = np.r_[((plus['R']-minus['R'])/(2*step)).ravel(),
                       (plus['scalar'][:11]-minus['scalar'][:11])/(2*step)]
            action = np.r_[exact['R'].ravel(), exact['scalar'][:11]]
            errors.append(np.max(np.abs(fd-action))/(1+np.max(np.abs(action))))
        assert min(errors)<1e-6, errors
        assert errors[-1]<errors[0]*.03 or errors[0]<1e-10, errors
        directional_errors.append(errors)
        # Isotropic raw-gradient shifts leave constitutive values unchanged.
        shifted = probe.closure(k, e, G+.37*np.eye(3))
        for key in ('R', 'a', 'beta'):
            near(shifted[key], value[key])
        near(shifted['scalar'][[0, 12]], value['scalar'][[0, 12]])
        near(shifted['scalar'][3]-value['scalar'][3], -2*k*.37)
        isotropic = probe.closure(k, e, G, (0, 0, np.eye(3)))
        near(isotropic['R'], np.zeros((3, 3)))
        near(isotropic['scalar'][3], -2*k)
    negative_raw = probe.closure(1,1,2*np.eye(3))
    near(negative_raw['scalar'][3],-4.)
    # Compare all 8 independent trace-free tangent columns with the exact
    # full reference tangent. The note's extension outside that space differs
    # from DevStrainRawWorkV1; compare the two after trace-free restriction.
    B = np.zeros((3, 3)); B[0, 2] = QE
    tangent = np.empty((9, 9))
    for i, H in enumerate(np.eye(9)):
        action = probe.closure(1, 1, B, (0, 0, H.reshape(3, 3)))
        tangent[:, i] = -action['R'].ravel()/CSTAR
    exactL = np.array(L, dtype=float)
    exactZ = np.array(Z, dtype=float)
    near(tangent@exactZ, exactL@exactZ, atol=2e-12)
    # At the planar reference I4=beta3=0 but their actions must be retained.
    action = probe.closure(1, 1, B, (0, 0, np.diag([1., -2., 1.])))
    assert abs(action['scalar'][15])>1 and abs(action['beta'][1])>.01
    # Required adverse witnesses pass by reproducing failure, without clipping.
    B[0, 2] = 3.5
    adverse = probe.closure(1, 1, B)
    xi = np.array([1., 0, -7])/math.sqrt(50)
    a = np.array([7., 0, 1])/math.sqrt(50)
    H = np.outer(a, xi)
    adverse_action = probe.closure(1, 1, B, (0, 0, H))
    principal = -np.sum(adverse_action['R']*H)/adverse['scalar'][1]
    near(adverse['scalar'][12], 4.130462147295278)
    near(principal, -.00633378103539413, atol=1e-13)
    assert np.linalg.eigvalsh(adverse['R']).min()>0
    near(adverse['scalar'][3], 1.03576095435, atol=1e-11)
    B = np.diag([2., 2., -4.]); B[0, 1] = 1.5*math.sqrt(17); B[1, 0] = -B[0, 1]
    unrealizable = probe.closure(1, 1, B)
    near(unrealizable['scalar'][12], 3.)
    near(unrealizable['R'], np.diag([-2/15, -2/15, 34/15]))
    near(unrealizable['scalar'][3], 48/5)
    return dict(fd_errors=directional_errors, reference_tangent_max_error=float(np.max(
        np.abs(tangent@exactZ-exactL@exactZ))), adverse_principal=principal,
        adverse_production=float(adverse['scalar'][3]), adverse_covariance_eigenvalues=
        np.linalg.eigvalsh(adverse['R']).tolist(), unrealizable=unrealizable['R'].tolist(),
        exact_LDL_pivots=list(map(str, pivots)), certificate_bound='nu_t_eff/250')


def neutral_checks(probe):
    for cd in (CSTAR, .087, .12):
        kap = math.sqrt(1.3*(1.82-1.44)*CSTAR**1.5/cd)
        for z in (1., 1.3, 5.):
            u,H,U = .5,1.,6.825550938766217
            n = probe.query(f'N {cd} {z} {u} {H} {U}')
            k, e = QE*u*u,u**3/(kap*z)
            near(n['neutral'], [QE,kap,U+u/kap*math.log(z/H),k,e,kap*u*z,u/(kap*z),
                                -e/z,-(cd/CSTAR)*u**4/(1.3*z),u*u,cd])
            G = np.zeros((3,3)); G[0,2]=u/(kap*z)
            expected = independent(k,e,G)
            near(n['R'], expected['R'])
            near(n['R'][0,2], -u*u)
            near(n['R'][1,1], 2*k/3)
            near(expected['scalar'][3], e)
            diffusion = (cd/CSTAR)*u**4/(1.3*z*z)
            near((1.44-1.82)*e*e/k+diffusion, 0., atol=1e-13)


def pressure_checks(probe):
    rng=np.random.default_rng(2911)
    for i in range(8):
        k=QE*.25
        raw=rng.normal(size=(3,3))
        M=raw+raw.T # Full strain, retaining trace.
        G=rng.normal(size=(3,3))
        R=independent(1.4,.8,G)['R']
        n=rng.normal(size=3); n/=np.linalg.norm(n)
        p,pi,potential,D=1.7,3.8,1.2,.023
        reply=probe.query(' '.join(map(str,['T',k,D,p,pi,potential,*M.ravel(),*R.ravel(),*n])))
        near(reply['traction'],-p*n-R@n)
        near(reply['pressure_jump'],D*n@M@n-2*k/3+n@R@n)
        near(reply['mechanical_pressure'],pi-2*k/3-potential)
    # Reference special case uses the complete outer R_zz.
    R=QE*.25*np.array([[2/3+20/81,0,-1/QE],[0,2/3,0],[-1/QE,0,2/3-20/81]])
    reply=probe.query(' '.join(map(str,['T',QE*.25,.2,0,0,0,*np.zeros(9),*R.ravel(),0,0,1])))
    near(reply['pressure_jump'],-.20906050250177272)


class WallOracle:
    def __init__(self, K, nu, H, h, zm, k1=0., k2=0.):
        self.K,self.nu,self.H,self.h,self.zm = K,nu,H,h,zm
        self.k1,self.k2 = k1,k2
        self.cl = KAPPA*QE**1.5
        self.a = KAPPA/math.sqrt(QE)*math.sqrt(K)
        self.t0 = K/QE
        self.A = (K/QE)**1.5/KAPPA
        ratio = h*math.sqrt(K)/(70*nu)
        b = 10.7 if h==0 else 10.7*(1-ratio**max(.5-.4*ratio**.7,1-.79*ratio**(-.28)))
        self.zv = nu*b/math.sqrt(K)
        self.ze = 2*self.cl*nu/math.sqrt(K)
        self.knots = [h,self.zv,self.ze,zm]
        self.j0 = self.integrate(lambda z: 1/self.D(z))
        self.j1 = self.integrate(lambda z: max(z-h,0)/self.D(z))
        self.j2 = self.integrate(lambda z: 1/self.D(z)**2)
        self.P = self.integrate(lambda z: self.J(z)*self.p0(z))
        self.minimum = self.integrate(lambda z: self.J(z)*self.A/max(z,zm))
        self.c = None
        if self.P>=self.minimum:
            if k1==0 and k2==0:
                self.c = H*math.exp(1-self.P/self.A)
            else:
                self.c = math.exp(brentq(lambda logc: self.integrate(
                    lambda z: self.J(z)*self.A/max(z,math.exp(logc)),
                    knots=[math.exp(logc)])-self.P, math.log(H)-50, math.log(zm), xtol=1e-14))
            self.knots.append(self.c)

    def D(self,z): return self.nu+self.nut(z)
    def nut(self,z): return self.a*max(z-self.zv,0)
    def J(self,z): return (1+self.k1*z)*(1+self.k2*z)
    def p0(self,z): return self.t0**2*self.nut(z)/self.D(z)**2
    def E(self,z): return self.A/max(z,self.c)
    def EA(self,z): return self.A/max(z,self.ze)

    def integrate(self,fun, lower=0., upper=None, knots=()):
        if upper is None: upper=self.H
        points = sorted(set([lower,upper]+[z for z in self.knots+list(knots) if lower<z<upper]))
        return sum(quad(fun,l,r,epsabs=3e-13,epsrel=3e-13,limit=150)[0]
                   for l,r in zip(points[:-1],points[1:]))

    def flux(self,z):
        return self.integrate(lambda s:self.J(s)*(self.E(s)-self.p0(s)),upper=z)

    def kw(self,z):
        return self.K-self.integrate(lambda s:self.flux(s)/(self.J(s)*self.D(s)),lower=z)


def wall_checks(probe):
    K = QE*.25
    summary = []
    # Canonical EP (chi=1, 0, partial), budget, smooth molecular segment,
    # curved metric, reversal and forcing at frozen K, and failed-gate baseline.
    cases = [(0,1.,0,1.5e-5,1.,.1,.5,0.,0.,None,[0.,0.]),
             (1,1.,0,1.5e-5,1.,.1,1.,0.,0.,None,[0.,0.]),
             (0,0.,0,1.5e-5,1.,.1,.5,0.,0.,None,[0.,0.]),
             (0,.35,0,1.5e-5,1.,.1,.5,0.,0.,[2.,-.7],[.4,-.2]),
             (0,1.,0,1.5e-5,1.,.1,.5,0.,0.,[-1.,.3],[1.,.4]),
             (1,1.,0,1.5e-5,1.,.1,1.,0.,0.,[-1.,.3],[1.,.4]),
             (1,1.,0,1.5e-5,1.,.1,1.,.2,-.1,None,[0.,0.]),
             (0,1.,0,1.5e-5,1.,.1,.5,.2,-.1,None,[0.,0.]),
             (0,1.,0,1.5e-5,1.,0.,.5,0.,0.,None,[0.,0.]),
             (1,1.,0,1.5e-5,1.,0.,1.,0.,0.,None,[0.,0.]),
             (0,0.,0,1.,1.,0.,.5,0.,0.,[.3,0.],[.1,0.]),
             (0,1.,1,1.,1.,0.,.5,0.,0.,[.3,0.],[.1,0.])]
    for number,(form,chi,fallback,nu,H,h,zm,k1,k2,U,F) in enumerate(cases):
        oracle = WallOracle(K,nu,H,h,zm,k1,k2)
        U = np.array([oracle.t0*oracle.j0,0.] if U is None else U)
        F = np.array(F)
        wall = probe.query(' '.join(map(str, ['W',form,chi,fallback,K,nu,H,h,zm,k1,k2,*U,*F])))
        ref = wall['reference']
        weight = chi if oracle.c is not None else 0.
        near(ref[[0,1,2,7,8,9,10,11,12,13]], [oracle.t0,oracle.t0*oracle.j0,oracle.A,
            oracle.zv,oracle.ze,CSTAR,oracle.cl,oracle.j0,oracle.j1,oracle.j2],atol=1e-10)
        if chi>0 and oracle.c is not None:
            near(ref[[3,4,5,6]], [oracle.P,oracle.minimum,oracle.c,weight],atol=1e-11)
        else:
            near(ref[5:7], [0.,0.])
        traction = (U-oracle.j1*F)/oracle.j0
        near(wall['traction'],traction)
        near(wall['interface_traction'],traction+F*(H-h))
        def shear(z): return (traction+F*max(z-h,0))/oracle.D(z)
        def velocity(z): return traction*oracle.integrate(lambda s:1/oracle.D(s),upper=z)+F*oracle.integrate(lambda s:max(s-h,0)/oracle.D(s),upper=z)
        def production(z): return oracle.nut(z)*np.sum(shear(z)**2)
        near(ref[14],oracle.integrate(production)/H,atol=1e-11)
        q = wall['quadrature']
        z,w,J = q[:,0],q[:,1],q[:,2]
        pw = w*J
        expected_source = []
        for row in q:
            z = row[0]
            P,p0,EA = production(z),oracle.p0(z),oracle.EA(z)
            E = oracle.E(z) if weight>0 else 0.
            dissipation = E if form==1 else (1-weight)*EA+weight*E
            redistribution = 0 if form==1 else weight*(E-p0)
            source = P-dissipation if form==1 else P-(1-weight)*EA-weight*p0
            near(row[3:8], [P,p0,dissipation,redistribution,source],atol=1e-10)
            near(row[10], oracle.D(z))
            near(row[11:13], velocity(z),atol=1e-10)
            near(row[13:15], shear(z),atol=1e-10)
            if form==1:
                near(row[8],oracle.flux(z),atol=1e-11)
                near(row[9],row[8]/(J[np.argmin(abs(q[:,0]-z))]*oracle.D(z)),atol=1e-11)
            expected_source.append(source)
        # Independent adaptive integrals of each source mode, vs actual C++
        # SourceLoads; weights in the probe are also verified independently.
        for n in range(6):
            expected_load = oracle.integrate(lambda z:oracle.J(z)*(z/H)**n*(production(z)-
                (oracle.E(z) if form==1 else (1-weight)*oracle.EA(z)+weight*oracle.p0(z))))
            near(wall['loads'][n],expected_load,atol=1e-10)
            near(np.dot(pw,(q[:,0]/H)**n),oracle.integrate(lambda z:oracle.J(z)*(z/H)**n),atol=1e-12)
            if form==1 and np.max(abs(U-[oracle.t0*oracle.j0,0]))<1e-12 and np.max(abs(F))==0:
                diffusion = -np.dot(w,q[:,8]*(0 if n==0 else n/H*(q[:,0]/H)**(n-1)))
                near(diffusion+wall['loads'][n],0.,atol=1e-10)
            if form==0 and chi==1 and fallback==0 and np.max(abs(F))==0:
                near(wall['loads'][n],0.,atol=1e-13)
        # Independent signed forced mechanical power, including molecular loss.
        total_loss = oracle.integrate(lambda z:oracle.D(z)*np.sum(shear(z)**2))
        forcing_work = np.dot(F, [oracle.integrate(lambda z:velocity(z)[i],lower=h) for i in range(2)])
        near(np.dot(traction+F*(H-h),U)-total_loss-forcing_work,0.,atol=1e-10)
        profile=wall['profile']
        for row in profile:
            z,kw,flux,pressure,normal,E = row
            if form==1:
                near(kw,oracle.kw(z),atol=1e-10)
                near(flux,oracle.flux(z),atol=1e-11)
                if h>0:
                    assert kw>0
                elif z==0:
                    # The note's positivity theorem is conditional on a rough
                    # branch. This smooth budget reference is a separate
                    # adverse witness, returned without clipping.
                    assert kw<0
            else:
                near(kw,K)
            if k1==0 and k2==0:
                near(pressure,-K/(81/20)+9.81*(H-z)+2*(K-kw)/3)
                near(normal,-K/(81/20)+9.81*H+2*K/3)
        if form==1:
            near(profile[[0,-1],2], [0,0],atol=1e-12)
        if number==0:
            near(ref[5], .015005743260176303)
            near(ref[11],27.30220375506487,atol=1e-11)
            near(ref[13],4984.456376333159,atol=1e-9)
            near(oracle.P,1.701714806838742)
            near(math.exp(-oracle.a*oracle.j0),.0054420488394214764)
            near(nu*oracle.t0**2*oracle.j2,.004672927852812337)
        if number==1:
            near(profile[:3,1],[1.8887084080086511,.907206835406087,.8495674195161094],atol=1e-11)
            near(profile[-1,3],-.20906050250177272)
        if number in (0,3,4,5):
            near(ref[5],.015005743260176303) # fixed K, forcing never changes cutoff
        if number in (8,9):
            assert oracle.zv>0 and np.any(q[:,3]==0) # real molecular-only segment
        summary.append(dict(case=number,formulation=form,chi=chi,cutoff=float(ref[5]),
                            source_loads=wall['loads'].tolist(),kw0=float(profile[0,1])))
    return summary


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--probe',required=True)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    subprocess.run([args.probe],check=True)
    probe=Probe(args.probe)
    constitutive_report=constitutive_checks(probe)
    neutral_checks(probe)
    pressure_checks(probe)
    wall_report=wall_checks(probe)
    probe.close()
    report=dict(status='PASS',checks=CHECKS,compiled_calls=probe.calls,
                maximum_absolute_comparison_error=MAX_ERROR,constitutive=constitutive_report,
                walls=wall_report,scope='Kernel verification only; no native coupled equilibrium or global realizability/stability proof.')
    Path(args.output).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status='PASS',checks=CHECKS,compiled_calls=probe.calls,
                         tangent_error=constitutive_report['reference_tangent_max_error'],
                         adverse_principal=constitutive_report['adverse_principal'])))


if __name__=='__main__':
    main()
