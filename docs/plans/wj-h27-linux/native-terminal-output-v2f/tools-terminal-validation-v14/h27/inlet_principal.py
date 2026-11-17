#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent planar H27 principal diagnostic; execute on Linux with python3 -B.

Loads the exact historical inlet through its unchanged implementation and stored
calibration, without refitting. Differentiates the planar WJ cubic analytically,
then checks witness actions by central differences. No governing PDE is solved.
Finite r/c scans provide witnesses, not a continuum universal certificate.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar
from wj_h5_neutral import NeutralReference, QE, KAPPA

CP = 1.8
Q_PEAK = 2 * CP / .65
H0 = np.array([[1., 0.], [0., -1.]])
H1 = np.array([[0., 1.], [0., 0.]])
H2 = H1.T.copy()
BASIS = (H0, H1, H2)


def root(s, w):
    def f(n):
        # Factored evaluation preserves the lower-endpoint sign as q -> 0.
        return (n*n-2*w)*(n-CP)-2.7*s*n
    if s == 0:
        return CP
    hi = CP + 2*math.sqrt(s-w) + 1
    while f(hi) <= 0:
        hi *= 2
    return brentq(f, CP, hi, xtol=1e-13, rtol=1e-14)


def wj(k, e, g, direction=None):
    """Exact planar reduction, trace-free velocity-gradient directions."""
    tau = k/e
    s = tau*(g+g.T)/2
    o = tau*(g-g.T)/2
    invs, invw = float(np.trace(s@s)), float(np.trace(o@o))
    n = root(invs, invw)
    den = n*n-2*invw
    b1, b4 = -1.2*n/den, -1.2/den
    comm = s@o-o@s
    r = k*(2*np.eye(2)/3+b1*s+b4*comm)
    if direction is None:
        return r
    ds = tau*(direction+direction.T)/2
    do = tau*(direction-direction.T)/2
    dis, diw = 2*float(np.trace(s@ds)), 2*float(np.trace(o@do))
    fn = 3*n*n-2*CP*n-2*invw-2.7*invs
    dn = (2.7*n*dis+2*(n-CP)*diw)/fn
    dd = 2*n*dn-2*diw
    db1, db4 = -1.2*(dn*den-n*dd)/den**2, 1.2*dd/den**2
    dr = k*(db1*s+b1*ds+db4*comm+b4*(ds@o+s@do-do@s-o@ds))
    return r, dr


def blend(y):
    if y <= .4:
        return 0.
    if y >= .75:
        return 1.
    q = (y-.4)/.35
    return q**3*(10-15*q+6*q*q)


def tensor(k, e, g, auxiliary, chi, direction=None):
    if direction is None:
        r = wj(k, e, g)
        r[0, 1] = r[1, 0] = chi*r[0, 1]-(1-chi)*auxiliary*(g[0, 1]+g[1, 0])
        return r
    r, dr = wj(k, e, g, direction)
    dr[0, 1] = dr[1, 0] = chi*dr[0, 1]-(1-chi)*auxiliary*(direction[0, 1]+direction[1, 0])
    return dr


def principal(c, r, k, e, gradient, nangles, refine=False, difference=False):
    g = np.array([[0., gradient], [0., 0.]])
    chi, auxiliary = blend(r/c.H), c.nt(r)
    # Normalize by positive split coefficient; sign does not depend on scale.
    rw = wj(k, e, g)
    scale = c.nu+(1-chi)*auxiliary+chi*(-rw[0, 1]/gradient)
    actions = [c.nu*(h+h.T)-tensor(k, e, g, auxiliary, chi, h) for h in BASIS]
    matrix = np.array([[np.sum(a*h) for a in actions] for h in BASIS])
    def direction(angle):
        cc, ss = math.cos(angle), math.sin(angle)
        return np.array([-cc*ss, -ss*ss, cc*cc])
    def form(angle):
        v = direction(angle)
        return float(v@matrix@v/scale)
    angles = np.linspace(0., math.pi, nangles, endpoint=False)
    values = [form(a) for a in angles]
    j = int(np.argmin(values)); angle = float(angles[j])
    if refine:
        delta = math.pi/nangles
        result = minimize_scalar(form, bounds=(angle-delta, angle+delta), method='bounded',
                                 options={'xatol': 1e-13})
        angle = float(result.x)
    v = direction(angle)
    h = sum(vv*hh for vv, hh in zip(v, BASIS))
    out = dict(r=float(r), r_over_H=float(r/c.H), k=float(k), epsilon=float(e),
               G=float(gradient), q=float(k*gradient/e), chi=chi,
               auxiliary=float(auxiliary), scale=float(scale), angle=angle,
               normal=[math.cos(angle), 0., math.sin(angle)],
               polarization=[-math.sin(angle), 0., math.cos(angle)],
               normalized_principal=form(angle), dimensional_principal=form(angle)*scale,
               minimum_covariance_over_k=float(min(np.linalg.eigvalsh(tensor(k, e, g, auxiliary, chi)).min(), 2*k/3)/k))
    if difference:
        expected = float(np.sum(h*(c.nu*(h+h.T)-tensor(k, e, g, auxiliary, chi, h)))/scale)
        checks = []
        for relative in (1e-3, 3e-4, 1e-4, 3e-5):
            step = relative*e/k
            dr = (tensor(k, e, g+step*h, auxiliary, chi)-tensor(k, e, g-step*h, auxiliary, chi))/(2*step)
            fd = float(np.sum(h*(c.nu*(h+h.T)-dr))/scale)
            checks.append(dict(relative_step=relative, normalized_fd=fd, error=fd-expected))
        out['central_difference_checks'] = checks
    return out


def fixed_inlet_scan(c, radii, cutoffs, nangles):
    # Historical U and k are evaluated once with their unchanged functions.
    inlet = [(r, c.k_inner(r), c.u_tau**2/c.D(r)) for r in radii]
    scans = []
    for cutoff in cutoffs:
        rows = [principal(c, r, k, c.A/max(r, cutoff), g, nangles) for r, k, g in inlet]
        best = min(rows, key=lambda v: v['normalized_principal'])
        witness = principal(c, best['r'], best['k'], best['epsilon'], best['G'], nangles, True, True)
        interior = min((v for v in rows if 0 < v['r'] < c.H), key=lambda v: v['normalized_principal'])
        interior = principal(c, interior['r'], interior['k'], interior['epsilon'], interior['G'], nangles, True, True)
        scans.append(dict(c=float(cutoff), c_over_H=float(cutoff/c.H),
                          min_witness=witness, positive_r_witness=interior,
                          negative_sample_count=sum(v['normalized_principal'] < 0 for v in rows),
                          sampled_radius_count=len(rows)))
    tail = []
    for r, k, g in inlet:
        if r > 0:
            value = principal(c, r, k, c.A/r, g, nangles)
            if value['normalized_principal'] < 0:
                value = principal(c, r, k, c.A/r, g, nangles, True, True)
                value['valid_for_all_cutoffs_leq'] = r
                tail.append(value)
    return scans, tail


def constant_k_state(c, r, cutoff):
    """D_inner -> infinity: k_ref -> K, still enforce reference total shear."""
    e, k, chi = c.A/max(r, cutoff), c.K, blend(r/c.H)
    alpha = c.nu+(1-chi)*c.nt(r)
    if chi == 0:
        gradient = c.u_tau**2/alpha
    else:
        def residual(q):
            gradient = q*e/k
            g = np.array([[0., gradient], [0., 0.]])
            return alpha*gradient-chi*wj(k, e, g)[0, 1]-c.u_tau**2
        if residual(Q_PEAK) < 0:
            raise ValueError('no declared increasing-shear-branch root in constant-k limit')
        q = brentq(residual, 0, Q_PEAK, xtol=1e-13, rtol=1e-13)
        gradient = q*e/k
    g = np.array([[0., gradient], [0., 0.]])
    production = -tensor(k, e, g, c.nt(r), chi)[0, 1]*gradient
    return k, e, gradient, production


def constant_k_limit(c, cutoffs, radii, nangles):
    def budget(cutoff):
        cuts = sorted(set([0., cutoff, .4*c.H, .75*c.H, c.H]))
        production = sum(quad(lambda r: constant_k_state(c, r, cutoff)[3], a, b,
                              epsabs=1e-11*c.A, epsrel=2e-9, limit=80)[0]
                         for a, b in zip(cuts[:-1], cuts[1:]))
        return production/c.A-1-math.log(c.H/cutoff)
    scans = []
    for cutoff in cutoffs:
        try:
            scans.append(dict(c=float(cutoff), budget_over_A=budget(cutoff)))
        except ValueError as error:
            scans.append(dict(c=float(cutoff), failure=str(error)))
    roots = []
    for a, b in zip(scans[:-1], scans[1:]):
        if 'failure' in a or 'failure' in b or a['budget_over_A']*b['budget_over_A'] > 0:
            continue
        cutoff = brentq(budget, a['c'], b['c'], xtol=1e-12, rtol=1e-12)
        rows = []
        for r in radii:
            k, e, g, _ = constant_k_state(c, r, cutoff)
            rows.append(principal(c, r, k, e, g, nangles))
        best = min(rows, key=lambda v: v['normalized_principal'])
        roots.append(dict(c=cutoff, budget_over_A=budget(cutoff),
                          witness=principal(c, best['r'], best['k'], best['epsilon'], best['G'], nangles, True, True)))
    return dict(scope='Zero endpoint reference k flux, k_ref=K limiting control only; does not change retained inlet.',
                scans=scans, roots=roots)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--reference-json', type=Path)
    parser.add_argument('--reference-cutoff', type=float, action='append', default=[])
    parser.add_argument('--cutoff-count', type=int, default=49)
    parser.add_argument('--radius-count', type=int, default=129)
    parser.add_argument('--angles', type=int, default=128)
    parser.add_argument('--skip-constant-k-limit', action='store_true')
    args = parser.parse_args()
    if min(args.cutoff_count, args.radius_count, args.angles) < 8:
        parser.error('bounded scans need at least eight points per coordinate')
    repo = Path(__file__).resolve().parents[5]
    calibration = repo/'docs/plans/wj-h5-linux/calibration.json'
    record = json.loads(calibration.read_text())
    historical = next(v for v in record['candidates'] if v['H_w_mm'] == 6.)
    c = NeutralReference(historical['H_w'], historical['z0_e'])
    for name, value in [('K', c.K), ('A', c.A), ('h_s', c.h), ('nu', c.nu), ('z0_e', c.z0_e)]:
        if float(value).hex() != float(historical[name]).hex():
            raise ValueError('historical calibration identity mismatch '+name)
    chosen = list(args.reference_cutoff)
    if args.reference_json:
        payload = json.loads(args.reference_json.read_text())
        for case in payload['cases']:
            if abs(case.get('parameters', {}).get('K', -1)/c.K-1) < 1e-12:
                chosen.extend(v['audit']['c'] for v in case.get('roots', []) if 'audit' in v)
    if any(not 0 < v <= c.H for v in chosen):
        raise ValueError('reference cutoff outside (0,H]')
    radii = sorted(set([0., c.H, .4*c.H, .75*c.H, c.w['rc'], c.h]
                       +np.linspace(0, c.H, args.radius_count).tolist()
                       +(c.H*np.geomspace(1e-8, .1, 32)).tolist()))
    cutoffs = sorted(set((c.H*np.exp(np.linspace(-12, 0, args.cutoff_count))).tolist()+chosen+[c.w['rc']]))
    scans, tail = fixed_inlet_scan(c, radii, cutoffs, args.angles)
    chosen_rows = [v for v in scans if v['c'] in chosen]
    outer = [principal(c, float(r), c.K, c.A/r, c.u_tau/(KAPPA*r), args.angles)
             for r in np.geomspace(c.H, c.top-c.d0, args.radius_count)]
    outer_best = min(outer, key=lambda v: v['normalized_principal'])
    outer_witness = principal(c, outer_best['r'], outer_best['k'], outer_best['epsilon'],
                              outer_best['G'], args.angles, True, True)
    ground_limit = [principal(c, 0., c.k_inner(0), c.A/(c.H*math.exp(logratio)),
                              c.u_tau**2/c.D(0), args.angles, True, True)
                    for logratio in (-4., -8., -12., -16., -24.)]
    output = dict(schema='WJH27IndependentFixedInletPrincipalV1',
        scope='Exact retained inlet, frozen k/Ew rank-one incompressible planar constitutive tangent; no PDE or reference root replaced.',
        calibration_sha256=hashlib.sha256(calibration.read_bytes()).hexdigest(),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        historical_functions='wj_h5_neutral.NeutralReference unchanged, no refit',
        H=c.H, K=c.K, A=c.A, nu=c.nu, cutoff_range=[min(cutoffs), max(cutoffs)],
        angle_count=args.angles, radius_count=len(radii), scans=scans,
        outer_inlet=dict(scope='Historical high-Re-derived inlet, finite-nu H27 momentum action; separate outer H owner',
                         radius_count=len(outer), witness=outer_witness,
                         inner_H_gradient=c.u_tau**2/c.D(c.H),
                         outer_H_gradient=c.u_tau/(KAPPA*c.H)),
        reference_selected_cutoffs=chosen_rows, cutoff_independent_negative_tail=tail,
        ground_small_cutoff_control=ground_limit,
        caveats=['Finite cutoff/radius sampling is not an all-cutoff proof.',
                 'Each negative tail witness applies to every cutoff c<=its positive radius.',
                 'As c approaches zero the ground principal tends to nu+nu_A*cos(2theta)^2, so ground alone cannot prove universal failure.',
                 'Frozen-scalar principal negativity is a necessary local momentum obstruction; positive values alone do not prove full coupled well-posedness.'])
    if not args.skip_constant_k_limit:
        output['infinite_inner_diffusion_reference_control'] = constant_k_limit(c, cutoffs, radii, args.angles)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        stream.write(json.dumps(output, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(output=str(args.output), scanned_cutoffs=len(scans),
          all_scanned_have_negative_witness=all(v['min_witness']['normalized_principal'] < 0 for v in scans),
          negative_tail_witnesses=len(tail), reference_selected=len(chosen_rows)), sort_keys=True))


if __name__ == '__main__':
    main()
