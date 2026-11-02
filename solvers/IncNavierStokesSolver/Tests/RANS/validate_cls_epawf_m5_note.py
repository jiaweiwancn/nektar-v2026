#!/usr/bin/env python3
"""Reproduce the revised note's exact circular-segment examples independently.

The native executable uses the implemented normal map; this oracle uses closed
physical-(s,r) formulas and independently formed boundary measures. No native
metric, gradient, stress or quadrature export constructs the reference answer.
"""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys

sys.dont_write_bytecode = True
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

CURVATURE = .05
STRESS = np.array([[.7, .1, -.2], [.1, -.3, .05], [-.2, .05, 1.1]])


def field(s, r):
    s, r = np.broadcast_arrays(np.asarray(s, float), np.asarray(r, float))
    c, J = CURVATURE, 1+CURVATURE*r
    zero = np.zeros_like(s)
    tangent = np.stack((np.cos(c*s), zero, -np.sin(c*s)), axis=-1)
    normal = np.stack((np.sin(c*s), zero, np.cos(c*s)), axis=-1)
    a, ap, app = 1+.1*np.cos(.4*s), -.04*np.sin(.4*s), -.016*np.cos(.4*s)
    u, w = a*r, -ap*r*r/(2*J)
    us, ur, ws, wr = ap*r, a, -app*r*r/(2*J), -ap*r*(2+c*r)/(2*J*J)
    U = u[..., None]*tangent+w[..., None]*normal
    ds = ((us+c*w)[..., None]*tangent+(ws-c*u)[..., None]*normal)/J[..., None]
    dr = ur[..., None]*tangent+wr[..., None]*normal
    G = np.einsum('...i,...j->...ij', ds, tangent)+np.einsum('...i,...j->...ij', dr, normal)
    k, D = 2+.1*np.sin(s)+.2*r*r, 1+.1*r
    ks, kss, kr, krr = .1*np.cos(s), -.1*np.sin(s), .4*r, .4
    grad_k = tangent*(ks/J)[..., None]+normal*kr[..., None]
    diffusion = D*krr+(.1+c*D/J)*kr+D*kss/(J*J)
    source = np.sum(U*grad_k, axis=-1)-diffusion
    flux = U*k[..., None]-D[..., None]*grad_k
    point = np.stack(((1/c+r)*np.sin(c*s), zero, (1/c+r)*np.cos(c*s)-1/c), axis=-1)
    return dict(U=U, G=G, k=k, D=D, grad_k=grad_k, scalar_flux=flux, source=source,
                J=J, normal=normal, tangent=tangent, u_t=u, w=w, point=point,
                divergence=np.trace(G, axis1=-2, axis2=-1))


def height(s):
    return 2+.2*np.sin(.3*np.asarray(s)), .06*np.cos(.3*np.asarray(s))


def boundary_flux(s, r, area, weights):
    f = field(s, r)
    mass = np.sum(f['U']*area, axis=-1)
    momentum = f['U']*mass[..., None]-np.einsum('ij,...j->...i', STRESS, area)
    scalar = np.sum(f['scalar_flux']*area, axis=-1)
    return dict(mass=float(weights@mass), momentum=weights@momentum, scalar=float(weights@scalar))


def reference(order=96):
    z, w = np.polynomial.legendre.leggauss(order)
    s, sw = -.25+2.75*z, 2.75*w
    H, Hp = height(s)
    top = field(s, H)
    faces = dict(top=boundary_flux(s, H, top['J'][:, None]*top['normal']-Hp[:, None]*top['tangent'], sw))
    ground = field(s, np.zeros_like(s))
    faces['ground'] = boundary_flux(s, np.zeros_like(s), -ground['normal'], sw)
    for name, side, position in (('left', -1, -3.), ('right', 1, 2.5)):
        h, _ = height(position)
        r, rw = h*(z+1)/2, h*w/2
        f = field(position, r)
        faces[name] = boundary_flux(position, r, side*f['tangent'], rw)
    ss = np.broadcast_to(s[:, None], (order, order))
    rr = H[:, None]*(z[None, :]+1)/2
    f = field(ss, rr)
    weights = sw[:, None]*w[None, :]*H[:, None]/2*f['J']
    acceleration = np.einsum('...ij,...j->...i', f['G'], f['U'])
    volume = dict(mass=float(np.sum(weights*f['divergence'])),
                  momentum=np.sum(weights[..., None]*acceleration, axis=(0, 1)),
                  scalar=float(np.sum(weights*f['source'])))
    total = {key: sum(face[key] for face in faces.values()) for key in ('mass', 'momentum', 'scalar')}
    return dict(order=order, faces=faces, boundary=total, volume=volume)


class Audit:
    def __init__(self):
        self.checks, self.errors = 0, {}

    def require(self, good, message):
        self.checks += 1
        if not good:
            raise RuntimeError(message)

    def near(self, actual, expected, name, tolerance=1e-10):
        a, b = np.asarray(actual, float), np.asarray(expected, float)
        self.require(a.shape == b.shape and np.all(np.isfinite(a)) and np.all(np.isfinite(b)), 'invalid '+name)
        error = float(np.max(np.abs(a-b), initial=0.))
        self.errors[name] = max(error, self.errors.get(name, 0.))
        self.require(error <= tolerance, f'{name}: {error} > {tolerance}')


def profile_check(path, audit):
    data = json.loads(Path(path).read_text())
    cases = [case for case in data['cases'] if case['name'] == 'curved_note' and case['chi'] == 1]
    audit.require(bool(cases), 'missing active curved_note native profile')
    records = []
    for case in cases:
        nu, H, rough, rK, rm = case['geometry']
        K, qstar, kappa = case['k'], case['q'], case['kappa']
        audit.near(case['metric'], [CURVATURE, 0.], 'note_profile_metric', 1e-15)
        root = np.sqrt(K)
        x = rough*root/(70*nu)
        power = max(.5-.4*x**.7, 1-.79*x**(-.28))
        zv = nu*10.7*(1-x**power)/root
        slope = .2295*root
        def p0(r):
            nt = slope*max(r-zv, 0.)
            return (K/qstar)**2*nt/(nu+nt)**2
        J = lambda r: 1+CURVATURE*r
        panels = sorted({0., H, *[r for r in (rough, zv) if 0 < r < H]})
        production = sum(quad(lambda r: J(r)*p0(r), lo, hi, epsabs=1e-13, epsrel=2e-13)[0]
                         for lo, hi in zip(panels[:-1], panels[1:]))
        A = (K/qstar)**1.5/kappa
        cap = lambda rc: A*(1+np.log(H/rc)+CURVATURE*H-CURVATURE*rc/2)
        logrc = brentq(lambda value: cap(np.exp(value))-production, -50., np.log(rm), xtol=2e-14, rtol=1e-14)
        rc = float(np.exp(logrc))
        split = sorted(set(panels+[rc]))
        redistribution = sum(quad(lambda r: J(r)*(A/max(r, rc)-p0(r)), lo, hi,
                                  epsabs=2e-13, epsrel=2e-13)[0] for lo, hi in zip(split[:-1], split[1:]))
        audit.near(rc, .009458614009, 'quoted_curved_cutoff', 5e-13)
        audit.near(case['reference'][8], rc, 'native_curved_cutoff', 5e-13)
        audit.near(case['reference'][4], production, 'curved_production_budget', 2e-12)
        audit.near(redistribution, 0., 'curved_redistribution_integral', 3e-12)
        audit.near(case['reference'][10], 0., 'native_curved_budget_residual', 3e-12)
        records.append(dict(order=case['order'], cutoff=rc, production_budget=production,
                            redistribution_integral=redistribution, rK=rK, rm=rm))
    return dict(path=str(Path(path).resolve()), cases=records)


def validate(native, profiles=None):
    audit = Audit()
    audit.require(native['schema'] == 1 and native['success'], 'native note probe failed')
    audit.near(native['curvature'], CURVATURE, 'curvature', 1e-15)
    audit.near(native['stress'], STRESS, 'constant_cartesian_stress', 1e-15)
    independent = reference()
    audit.require([level['order'] for level in native['levels']] == [24, 48, 64], 'missing finite-quadrature sequence')
    for level in native['levels']:
        for face in ('top', 'ground', 'left', 'right'):
            for key in ('mass', 'momentum', 'scalar'):
                audit.near(level['faces'][face][key], independent['faces'][face][key], 'face_'+face+'_'+key)
        for region in ('boundary', 'volume'):
            for key in ('mass', 'momentum', 'scalar'):
                audit.near(level[region][key], independent[region][key], region+'_'+key)
        audit.near(level['boundary']['mass'], 0., 'closed_mass')
        audit.near(level['boundary']['momentum'], level['volume']['momentum'], 'closed_momentum')
        audit.near(level['boundary']['scalar'], level['volume']['scalar'], 'closed_scalar')
        audit.require(level['maximum_divergence'] < 2e-12, 'native nonzero divergence')
    samples = native['finite_difference_samples']
    audit.require(len(samples) == 9, 'missing nine Cartesian FD samples')
    for sample in samples:
        expected = field(sample['s'], sample['r'])
        for key in ('point', 'normal', 'tangent', 'U', 'G', 'J', 'u_t', 'w', 'k', 'D', 'grad_k', 'scalar_flux', 'source'):
            audit.near(sample['field'][key], expected[key], 'sample_'+key, 3e-12)
        audit.near(sample['G_fd'], expected['G'], 'cartesian_inverse_FD', 2e-10)
        audit.near(sample['field']['divergence'], 0., 'sample_divergence', 2e-12)
    sample = native['matching_sample']; s = sample['s']; H, Hp = height(s); f = field(s, H)
    area = f['J']*f['normal']-Hp*f['tangent']
    mass, naive = float(f['U']@area), float(f['J']*f['w'])
    audit.near(sample['area'], area, 'variable_height_area', 3e-12)
    audit.near([sample['field']['H'], sample['field']['w'], sample['m_H'], sample['naive_Jw']],
               [H, f['w'], mass, naive], 'matching_values', 3e-12)
    for name, actual, quoted, tolerance in (
            ('quoted_H', H, 2.0760376830, 5e-11), ('quoted_wH', f['w'], .03880261243, 6e-12),
            ('quoted_mH', mass, -.08237640447, 6e-12), ('quoted_naive_Jw', naive, .04283039671, 6e-12)):
        audit.near(actual, quoted, name, tolerance)
    audit.require(mass < 0 < naive, 'height-slope omission did not reverse the inferred sign')
    pressure = quad(lambda r: CURVATURE*r*r/(1+CURVATURE*r), 0., 2., epsabs=1e-14, epsrel=1e-13)[0]
    closed = 2-2/CURVATURE+np.log1p(2*CURVATURE)/CURVATURE**2
    audit.near(pressure, closed, 'closed_normal_pressure', 1e-12)
    audit.near(pressure, .12407192173, 'quoted_pressure_increment', 6e-12)
    for level in native['levels']:
        audit.near(level['pressure_increment'], pressure, 'native_normal_pressure', 1e-12)
    optional = profile_check(profiles, audit) if profiles else None
    return dict(success=True, checks=audit.checks, maximum_absolute_errors=audit.errors,
                independent_quadrature_order=independent['order'], pressure_increment=pressure,
                matching=dict(H=float(H), wH=float(f['w']), mH=mass, naive_Jw=naive), profiles=optional,
                scope='Exact revised-note geometry and finite-quadrature flux examples; no CFD time evolution, inner epsilon PDE or stability claim')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--profiles', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    native = args.output/'native-note.json'
    command = [str(args.probe.resolve()), '--output', str(native.resolve())]
    environment = dict(os.environ)
    environment.update({name: '1' for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS')})
    report = dict(schema=1, started_utc=datetime.now(timezone.utc).isoformat(), command=command, success=False)
    try:
        result = subprocess.run(command, env=environment, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
        (args.output/'native-note.log').write_text(result.stdout)
        if result.returncode:
            raise RuntimeError(f'native note probe exited {result.returncode}; see native-note.log')
        report.update(validate(json.loads(native.read_text()), args.profiles))
    except Exception as error:
        report['error'] = str(error)
        raise
    finally:
        report['finished_utc'] = datetime.now(timezone.utc).isoformat()
        (args.output/'validation.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(success=report['success'], checks=report['checks'])))


if __name__ == '__main__':
    main()
