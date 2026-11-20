#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent Linux-only native finite-perturbation projection audit.

Uses signed native maps as data, exact rational Modified-A Q3 polynomials,
Bernstein global bounds and independently assembled pressure test moments.
Synthetic Kbase=Kin trace compatibility is not physical-reference qualification.
"""
import argparse
from fractions import Fraction as F
import json
import math
from pathlib import Path
import sys
import time
import perturbation_certificate as pc
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from numpy.polynomial.legendre import leggauss

Q3 = ([F(1), F(-1)], [F(0), F(1)], [F(0), F(1), F(-1)],
      [F(0), F(-2), F(6), F(-4)])
EDGES = ((0, 4, 8, 12), (1, 5, 9, 13), (0, 1, 2, 3), (4, 5, 6, 7))
REPRESENTATION_LIMIT = 1e-12
DIVERGENCE_LIMIT = 1e-13
EPSILON_PEAK_LIMIT = 1e-5


def scale(a, factor): return [[factor*v for v in row] for row in a]
def subtract(a, b): return pc.padd(a, scale(b, F(-1)))


def derivative(a, axis, width):
    if axis == 0: return [[i*v/width for v in a[i]] for i in range(1, len(a))] or [[F(0)]*len(a[0])]
    return [[j*row[j]/width for j in range(1, len(row))] or [F(0)] for row in a]


def restrict(a, old, new):
    old, new = list(map(pc.rational, old)), list(map(pc.rational, new))
    ax, bx = (new[0]-old[0])/(old[1]-old[0]), (new[1]-new[0])/(old[1]-old[0])
    ar, br = (new[2]-old[2])/(old[3]-old[2]), (new[3]-new[2])/(old[3]-old[2])
    out = [[F(0)]*len(a[0]) for _ in a]
    for i, row in enumerate(a):
        for j, v in enumerate(row):
            for k in range(i+1):
                for l in range(j+1): out[k][l] += v*math.comb(i, k)*ax**(i-k)*bx**k*math.comb(j, l)*ar**(j-l)*br**l
    return out


def polynomial(e, key, coefficients, offset=0):
    signs = e[dict(u='us', p='ps', k='ks', epsilon='epsilon_signs')[key]]
    width = 3 if key == 'p' else 4
    out = [[F(0)]*width for _ in range(width)]
    for local, (gid, sign) in enumerate(zip(e[key], signs)):
        coefficient = pc.rational(coefficients[offset+gid])*pc.rational(sign)
        a, b = Q3[local % width], Q3[local//width]
        for i, u in enumerate(a):
            for j, v in enumerate(b): out[i][j] += coefficient*u*v
    return out


def value(a, x, r):
    return sum((v*x**i*r**j for i, row in enumerate(a) for j, v in enumerate(row)), F(0))


def coefficient_bound(a):
    return max(abs(v) for row in pc.bernstein_exact(a) for v in row)


def native_shapes(packet, certificate, deadline, *, synthetic_base=True):
    layout = certificate['normalization']['layout']; H = layout['H']; a, b = layout['r_knots'][0], layout['r_knots'][-1]
    n = packet['sizes'][0]; kin, Ain = pc.rational(packet['Kin']), pc.rational(packet['Ain'])
    C = pc.rational(packet['normalization']['streamfunction_coefficient']); sk, se = packet['signs'][1:]
    states = dict(velocity=packet['velocity'], k=packet['k_increment'], epsilon=packet['epsilon_independent_increment'])
    pieces = {key: [] for key in states}; error_energy = {key: F(0) for key in states}
    errors = {key: F(0) for key in states}; leaks = errors.copy(); div_bound = F(0)
    gradient_energy = {key: F(0) for key in states}; gradient_bounds = {key: F(0) for key in states}; base_error = F(0)
    gradient_error_bounds = {key: F(0) for key in states}
    for e in packet['elements']:
        pc.clock(deadline); bounds = e['bounds']; xl, xr, rl, rr = bounds
        up = [polynomial(e, 'u', states['velocity'], c*n) for c in (0, 1)]
        kp = scale(polynomial(e, 'k', states['k']), 1/kin)
        if synthetic_base:
            base_error = max(base_error, coefficient_bound(subtract(scale(polynomial(e, 'k', packet['synthetic_k_base']), 1/kin), [[F(1)]])))
        ep = None
        if e['epsilon']:
            physical = polynomial(e, 'epsilon', states['epsilon'])
            ep = [pc.mul(row, [pc.rational(rl)/Ain, (pc.rational(rr)-pc.rational(rl))/Ain]) for row in physical]
        div = pc.padd(derivative(up[0], 0, pc.rational(xr)-pc.rational(xl)),
                     derivative(up[1], 1, pc.rational(rr)-pc.rational(rl)))
        div_bound = max(div_bound, coefficient_bound(div))
        xc = sorted(set([xl, xr]+[v for v in layout['x_knots'] if xl < v < xr]))
        rc = sorted(set([rl, rr]+[v for v in layout['r_knots']+[H, layout['epsilon_peak_radius']] if rl < v < rr]))
        for l, h in zip(xc[:-1], xc[1:]):
            for lo, hi in zip(rc[:-1], rc[1:]):
                cell = [l, h, lo, hi]; measure = (pc.rational(h)-pc.rational(l))*(pc.rational(hi)-pc.rational(lo))
                bx = pc.spline_piece(layout['x_knots'], l, h); dx = pc.derivative(bx, pc.rational(h)-pc.rational(l))
                br = pc.spline_piece(layout['r_knots'], lo, hi); dr = pc.derivative(br, pc.rational(hi)-pc.rational(lo))
                targets = dict(velocity=[scale(pc.outer(bx, dr), C), scale(pc.outer(dx, br), -C)],
                    k=[pc.outer([pc.rational(.01)*sk*F(3, 2)*v for v in bx], pc.hermite_piece(lo, hi, a, H, b))])
                native = dict(velocity=[restrict(v, bounds, cell) for v in up], k=[restrict(kp, bounds, cell)])
                if ep is not None:
                    targets['epsilon'] = [pc.outer([pc.rational(.01)*se*F(3, 2)*v for v in bx],
                        pc.hermite_piece(lo, hi, H, layout['epsilon_peak_radius'], b))]
                    native['epsilon'] = [restrict(ep, bounds, cell)]
                for key, polynomials in native.items():
                    squared = [[F(0)]]; squared_error = [[F(0)]]; component_error = F(0)
                    gradient_squared = [[F(0)]]; gradient_max_squared = F(0)
                    gradient_error_squared = F(0)
                    for actual, target in zip(polynomials, targets[key]):
                        squared = pc.padd(squared, pc.psquare(actual)); difference = subtract(actual, target)
                        squared_error = pc.padd(squared_error, pc.psquare(difference)); component_error += coefficient_bound(difference)**2
                        for axis, width in ((0, pc.rational(h)-pc.rational(l)), (1, pc.rational(hi)-pc.rational(lo))):
                            gradient = derivative(actual, axis, width)
                            gradient_squared = pc.padd(gradient_squared, pc.psquare(gradient))
                            gradient_max_squared += coefficient_bound(gradient)**2
                            if not synthetic_base:
                                gradient_error_squared += coefficient_bound(derivative(difference, axis, width))**2
                    errors[key] = max(errors[key], component_error)
                    error_energy[key] += pc.integral(dict(power=squared_error, area=measure))
                    gradient_energy[key] += pc.integral(dict(power=gradient_squared, area=measure))
                    gradient_bounds[key] = max(gradient_bounds[key], gradient_max_squared)
                    gradient_error_bounds[key] = max(gradient_error_bounds[key], gradient_error_squared)
                    inside = l >= layout['x_knots'][0] and h <= layout['x_knots'][-1] and lo >= (H if key == 'epsilon' else a) and hi <= b
                    if not inside: leaks[key] = max(leaks[key], coefficient_bound(squared))
                    pieces[key].append(dict(bounds=cell, area=measure, power=squared, inside_support=inside))
    diagnostics = dict(global_representation_error_upper={key: pc.sqrt_bounds(v)[1] for key, v in errors.items()},
        representation_error_energy=error_energy, outside_support_peak_upper={key: pc.sqrt_bounds(v)[1] for key, v in leaks.items()},
        strong_divergence_MAX_upper=pc.directed(div_bound, True), gradient_energy=gradient_energy,
        gradient_MAX_upper={key: pc.sqrt_bounds(v)[1] for key, v in gradient_bounds.items()},
        synthetic_Kbase_constant_error_upper=pc.directed(base_error, True))
    if not synthetic_base:
        del diagnostics['synthetic_Kbase_constant_error_upper']
        diagnostics['gradient_representation_error_upper'] = {
            key: pc.sqrt_bounds(v)[1] for key, v in gradient_error_bounds.items()}
    return pieces, diagnostics


def basis(e, key, x, r):
    width = 3 if key == 'p' else 4; xl, xr, rl, rr = e['bounds']
    tx, tr = (x-xl)/(xr-xl), (r-rl)/(rr-rl)
    def one(t):
        b = np.array([np.polynomial.polynomial.polyval(t, list(map(float, v))) for v in Q3[:width]]).T
        d = np.array([np.polynomial.polynomial.polyval(t, list(map(float, pc.derivative(v)))) for v in Q3[:width]]).T
        return b, d
    bx, dx = one(tx); br, dr = one(tr); sign = np.asarray(e[dict(u='us', p='ps', k='ks', epsilon='epsilon_signs')[key]])
    return ((br[:, :, None]*bx[:, None, :]).reshape(len(x), -1)*sign,
            (br[:, :, None]*dx[:, None, :]).reshape(len(x), -1)*sign/(xr-xl),
            (dr[:, :, None]*bx[:, None, :]).reshape(len(x), -1)*sign/(rr-rl))


def divergence_audit(packet, deadline):
    n, npres = packet['sizes'][:2]; U = np.asarray(packet['velocity']); q, w = leggauss(8)
    rows, cols, vals = [], [], []; action = np.zeros(npres)
    for e in packet['elements']:
        pc.clock(deadline); xl, xr, rl, rr = e['bounds']
        x, r = np.meshgrid((xl+xr)/2+(xr-xl)*q/2, (rl+rr)/2+(rr-rl)*q/2)
        weights = np.outer(w, w).ravel()*(xr-xl)*(rr-rl)/4; x, r = x.ravel(), r.ravel()
        _, dx, dr = basis(e, 'u', x, r); p = basis(e, 'p', x, r)[0]
        ui, pi = np.asarray(e['u']), np.asarray(e['p'])
        action[pi] += p.T@(weights*(dx@U[ui]+dr@U[n+ui]))
        block = p.T@(weights[:, None]*p)
        rows.extend(np.repeat(pi, len(pi))); cols.extend(np.tile(pi, len(pi))); vals.extend(block.ravel())
    mass = sp.csc_matrix((vals, (rows, cols)), shape=(npres, npres)); solve = spla.factorized(mass)
    domain = (packet['domain'][1]-packet['domain'][0])*(packet['domain'][3]-packet['domain'][2])
    norm = lambda v: math.sqrt(max(0., float(v@solve(v)))/domain)
    residual = norm(action)
    result = dict(independent_D_action_dual_RMS=residual, pass_=residual <= DIVERGENCE_LIMIT)
    if 'divergence_action' in packet:
        difference = norm(action-np.asarray(packet['divergence_action']))
        result.update(native_D_action_difference_dual_RMS=difference,
                      native_reported_dual_RMS=packet['native_divergence_dual_RMS'],
                      pass_=result['pass_'] and difference <= DIVERGENCE_LIMIT)
    return result


def trace_k(packet, x):
    result = np.empty(len(x)); H = packet['H']; coeff = np.asarray(packet['synthetic_k_new'])
    for i, xx in enumerate(x):
        owners = [e for e in packet['elements'] if e['epsilon'] and abs(e['bounds'][2]-H) < 1e-12 and
                  e['bounds'][0]-1e-12 <= xx <= e['bounds'][1]+1e-12]
        pc.require(owners, 'missing independent outer-H owner')
        values = [float((basis(e, 'k', np.asarray([xx]), np.asarray([H]))[0]@coeff[e['k']])[0]) for e in owners]
        pc.require(max(values)-min(values) <= 1e-12*packet['Kin'], 'K trace depends on shared owner')
        result[i] = values[0]
    return result


def gamma_projection(packet, order):
    q, w = leggauss(order); H = packet['H']; Ain = packet['Ain']; cl = packet['c_length']
    result = np.zeros(packet['sizes'][3]); seen = {}; target_max = 0.; target_error = 0.
    native = np.asarray(packet['gamma_new_lift'])
    for e in packet['elements']:
        if not e['epsilon']: continue
        xl, xr, rl, rr = e['bounds']
        sides = ([0] if abs(xl-packet['domain'][0]) < 1e-12 else [])+([2] if abs(rl-H) < 1e-12 else [])
        for side in sides:
            if side == 2:
                x, r, weight = (xl+xr)/2+(xr-xl)*q/2, np.full(order, H), w*(xr-xl)/2
                ex, er = np.array([xl, xr]), np.full(2, H)
            else:
                x, r, weight = np.full(order, xl), (rl+rr)/2+(rr-rl)*q/2, w*(rr-rl)/2
                ex, er = np.full(2, xl), np.array([rl, rr])
            target = (trace_k(packet, x)**1.5/cl-Ain)/H if side == 2 else np.zeros(order)
            endpoints = (trace_k(packet, ex)**1.5/cl-Ain)/H if side == 2 else np.zeros(2)
            active = list(EDGES[side]); signs = np.asarray(e['epsilon_signs'])[active]
            b = basis(e, 'epsilon', x, r)[0][:, active]
            coefficient = np.zeros(4); coefficient[:2] = endpoints/signs[:2]
            rhs = target-b[:, :2]@coefficient[:2]
            coefficient[2:] = np.linalg.solve(b[:, 2:].T@(weight[:, None]*b[:, 2:]), b[:, 2:].T@(weight*rhs))
            for local, value_ in zip(active, coefficient):
                gid = e['epsilon'][local]
                pc.require(gid not in seen or abs(seen[gid]-value_) <= 1e-12*Ain/H, 'independent Gamma corner conflict')
                seen[gid] = value_; result[gid] = value_
            if side == 2:
                sx = np.r_[x, np.linspace(xl, xr, 65)]; sr = np.full(len(sx), H)
                physical_target = trace_k(packet, sx)**1.5/(cl*H)
                target_max = max(target_max, float(np.max(abs(physical_target/(Ain/H)-1))))
                target_error = max(target_error, float(np.max(abs(basis(e, 'epsilon', sx, sr)[0]@native[e['epsilon']]-physical_target)))/(Ain/H))
    essential = sorted(seen)
    pc.require(set(packet['epsilon_free_ids']) == set(range(len(result)))-set(essential), 'wrong explicit epsilon essential mask')
    return result, target_max, target_error


def synthetic_gamma_audit(packet, table):
    import coupled_mms_check
    reference = coupled_mms_check.IndependentMMS(table)
    for key, expected in (('Kin', reference.Kin), ('Ain', reference.Ain), ('H', reference.H),
                          ('d0', reference.d0), ('c_length', reference.cl)):
        pc.require(abs(packet[key]-expected) <= 1e-12*abs(expected), 'native reference scale differs: '+key)
    expected24, _, _ = gamma_projection(packet, 24); expected40, peak, point_error = gamma_projection(packet, 40)
    scale = packet['Ain']/packet['H']; actual = np.asarray(packet['gamma_increment_lift'])
    error = float(np.max(abs(actual-expected40), initial=0.))/scale
    quadrature = float(np.max(abs(expected40-expected24), initial=0.))/scale
    point = float(trace_k(packet, np.array([0.]))[0]); K = packet['Kin']; sign = packet['signs'][1]
    K_error = abs((point-K)/K-sign*.01)
    cutoff = reference.cutoff(point); native = packet['synthetic_interface_peak']
    cutoff_error = abs(cutoff-native['cutoff_new'])/packet['H']
    return dict(role='Synthetic Kbase=Kin only; no developed-state qualification',
        dependent_Gamma_increment_coefficient_error=error, Gamma_projection_quadrature_change=quadrature,
        sampled_physical_Gamma_projection_target_MAX=point_error, induced_relative_Gamma_change_MAX=peak,
        interface_K_at_x0=point, one_percent_K_increment_error=K_error,
        cutoff_base=reference.cutoff(K), cutoff_new=cutoff, cutoff_change=cutoff-reference.cutoff(K),
        cutoff_native_error_over_H=cutoff_error,
        pass_=error <= 1e-10 and quadrature <= 1e-12 and peak <= .02 and K_error <= REPRESENTATION_LIMIT and cutoff_error <= 1e-12)


def helper_audit(packet, certificate):
    layout = packet['layout']; errors = dict(unit_velocity=0., unit_gradient=0., k_jet=0., epsilon_jet=0.)
    span = layout['r_knots'][-1]-layout['r_knots'][0]; peak = certificate['normalization']['unitPeakEstimate']
    def jet(knots, coordinate):
        i = next(j for j in range(4) if knots[j] <= coordinate < knots[j+1])
        a, b = knots[i:i+2]; p = pc.spline_piece(knots, a, b); width = pc.rational(b)-pc.rational(a)
        t = (pc.rational(coordinate)-pc.rational(a))/width
        def evaluate(q): return float(sum((v*t**j for j, v in enumerate(q)), F(0)))
        d = pc.derivative(p, width)
        return evaluate(p), evaluate(d), evaluate(pc.derivative(d, width))
    for row in packet['helper_samples']:
        x, r = row['x'], row['z']-packet['d0']; a, ax, axx = jet(layout['x_knots'], x); b, br, brr = jet(layout['r_knots'], r)
        u = np.array([a*br, 0., -ax*b]); g = np.array([[ax*br, 0., a*brr], [0., 0., 0.], [-axx*b, 0., -ax*br]])
        errors['unit_velocity'] = max(errors['unit_velocity'], float(np.max(abs(u-row['unit_velocity'])))/peak)
        errors['unit_gradient'] = max(errors['unit_gradient'], float(np.max(abs(g-row['unit_gradient'])))*span/peak)
        for field in ('k', 'epsilon'):
            lower, center, upper = ((layout['r_knots'][0], layout['H'], layout['r_knots'][-1]) if field == 'k' else
                (layout['H'], layout['epsilon_peak_radius'], layout['r_knots'][-1]))
            if not lower < r < upper: target = np.zeros(3)
            else:
                width = center-lower if r <= center else upper-center
                t = (r-lower)/width if r <= center else (upper-r)/width
                h = t*t*(3-2*t); dh = (1 if r <= center else -1)*6*t*(1-t)/width
                bump = np.array([1.5*a*h, 1.5*ax*h, 1.5*a*dh])
                if field == 'k': target = packet['signs'][1]*.01*packet['Kin']*bump
                else: target = packet['signs'][2]*.01*packet['Ain']*np.array([bump[0]/r, bump[1]/r, bump[2]/r-bump[0]/r**2])
            ref = packet['Kin'] if field == 'k' else packet['Ain']/r
            delta = abs(target-row[field+'_jet'])/np.array([ref, ref/span, ref/span])
            errors[field+'_jet'] = max(errors[field+'_jet'], float(np.max(delta)))
    return dict(reference_scaled_errors=errors, pass_=all(v <= 1e-10 for v in errors.values()))


def represented_norms(packet, certificate, deadline, *, synthetic_base=True):
    domain = (pc.rational(packet['domain'][1])-pc.rational(packet['domain'][0]))*(pc.rational(packet['domain'][3])-pc.rational(packet['domain'][2]))
    outer_domain = (pc.rational(packet['domain'][1])-pc.rational(packet['domain'][0]))*(pc.rational(packet['domain'][3])-pc.rational(packet['H']))
    pieces, diagnostics = native_shapes(packet, certificate, deadline, synthetic_base=synthetic_base)
    norms = {}; all_ok = True
    for field, rows in pieces.items():
        pc.clock(deadline); peak = pc.maximum(rows, deadline, relative_gap=1e-11 if field == 'k' else pc.PEAK_GAP)
        peak_bounds = [pc.sqrt_bounds(pc.rational(peak['lower']))[0], pc.sqrt_bounds(pc.rational(peak['upper']))[1]]
        support_rows = [r for r in rows if r['inside_support']]; support = sum((r['area'] for r in support_rows), F(0))
        full_energy = sum((pc.integral(r) for r in rows), F(0)); support_energy = sum((pc.integral(r) for r in support_rows), F(0))
        area = outer_domain if field == 'epsilon' else domain
        full_rms, support_rms = pc.sqrt_bounds(full_energy/area), pc.sqrt_bounds(support_energy/support)
        half = pc.enclosed_half_area(support_rows, peak['lower'], peak['upper'], deadline)
        pc.require(diagnostics['outside_support_peak_upper'][field] < peak_bounds[0]/2, 'half-peak set extends outside declared support')
        expected = .001 if field == 'velocity' else .01
        peak_limit = 1e-8*expected if field == 'velocity' else (REPRESENTATION_LIMIT if field == 'k' else EPSILON_PEAK_LIMIT)
        peak_ok = max(abs(v-expected) for v in peak_bounds) <= peak_limit
        representation = diagnostics['global_representation_error_upper'][field]
        representation_ok = math.isfinite(representation) if field == 'epsilon' else representation <= REPRESENTATION_LIMIT
        good = (peak_ok and representation_ok and full_rms[0] >= (1e-5 if field == 'velocity' else 1e-4) and
                support_rms[0] >= (.0001 if field == 'velocity' else .001) and half[-1]['lower'] >= .05 and
                diagnostics['outside_support_peak_upper'][field] <= REPRESENTATION_LIMIT)
        all_ok = all_ok and good
        norms[field] = dict(global_peak_enclosure=peak_bounds, peak_subdivision=peak,
            peak_deviation_limit=peak_limit, peak_pass=bool(peak_ok), representation_pass=bool(representation_ok),
            domain_RMS_enclosure=full_rms, support_RMS_enclosure=support_rms,
            energy_enclosure=[pc.directed(full_energy, False), pc.directed(full_energy, True)],
            support_area_enclosure=[pc.directed(support, False), pc.directed(support, True)],
            half_peak_area_refinement=half, representation_error_RMS_enclosure=pc.sqrt_bounds(diagnostics['representation_error_energy'][field]/area),
            normalized_field_gradient_RMS_enclosure=pc.sqrt_bounds(diagnostics['gradient_energy'][field]/area),
            normalized_field_gradient_MAX_upper=diagnostics['gradient_MAX_upper'][field], pass_=bool(good))
    diagnostics.pop('representation_error_energy'); diagnostics.pop('gradient_energy')
    return norms, diagnostics, all_ok


def audit(packet, certificate, table, deadline):
    pc.require(packet['schema'] == 'WJH27PerturbationNativeProbeV1' and packet['mode'] == 'project', 'wrong native projection ABI')
    xs, rs = pc.validate_inventory(packet, required_mode='project')
    pc.require(packet['helper_identity'] == 'WJH27FinitePerturbationB3V1', 'wrong frozen native shape helper')
    pc.require(packet['layout'] == certificate['normalization']['layout'], 'refinement reselected the frozen layout')
    for key in ('unitPeakLower', 'unitPeakEstimate', 'unitPeakUpper'):
        pc.require(packet['normalization'][key] == certificate['normalization'][key], 'native peak normalization differs from certificate')
    pc.require(packet['normalization']['streamfunction_coefficient'] == packet['signs'][0]*.001/certificate['normalization']['unitPeakEstimate'], 'native streamfunction coefficient/sign differs')
    pc.require(packet['normalization_load_twice_rejected'] and packet['Gamma_corner_conflict_rejected'], 'native rejection control failed')
    pc.require(all(v >= 4 for v in packet['support_intervals']) and packet['support_area'] >= 1.5, 'support is too narrow or underresolved')
    for key, size, count, sign in (('u', packet['sizes'][0], 16, 'us'), ('p', packet['sizes'][1], 9, 'ps'),
                                  ('k', packet['sizes'][2], 16, 'ks'), ('epsilon', packet['sizes'][3], 16, 'epsilon_signs')):
        covered = set()
        for e in packet['elements']:
            if key == 'epsilon' and not e[key]:
                pc.require(not e[sign] and e['bounds'][3] <= packet['H']+1e-12, 'wrong epsilon subset'); continue
            pc.require(len(e[key]) == len(e[sign]) == count and len(set(e[key])) == count and
                       all(0 <= v < size for v in e[key]) and all(s in (-1., 1.) for s in e[sign]), 'invalid signed native '+key+' map')
            covered.update(e[key])
        pc.require(covered == set(range(size)), 'incomplete native '+key+' map')
    norms, diagnostics, all_ok = represented_norms(packet, certificate, deadline)
    divergence = divergence_audit(packet, deadline); gamma = synthetic_gamma_audit(packet, table); helper = helper_audit(packet, certificate)
    essential = dict(velocity=float(np.max(abs(np.asarray(packet['velocity_essential_lift'])), initial=0.)),
        k=float(np.max(abs(np.asarray(packet['k_essential_lift'])), initial=0.))/packet['Kin'],
        epsilon=float(np.max(abs(np.asarray(packet['epsilon_independent_essential_lift'])), initial=0.))/(packet['Ain']/packet['H']))
    essential_ok = all(v <= REPRESENTATION_LIMIT for v in essential.values())
    composition = dict(k=float(np.max(abs(np.asarray(packet['synthetic_k_new'])-np.asarray(packet['synthetic_k_base'])-np.asarray(packet['k_increment'])), initial=0.))/packet['Kin'],
        epsilon=float(np.max(abs(np.asarray(packet['epsilon_total_increment'])-np.asarray(packet['gamma_increment_lift'])-np.asarray(packet['epsilon_independent_increment'])), initial=0.))/(packet['Ain']/packet['H']))
    composition_ok = max(composition.values()) <= REPRESENTATION_LIMIT and diagnostics['synthetic_Kbase_constant_error_upper'] <= REPRESENTATION_LIMIT
    pc.require(np.max(abs(np.asarray(packet['velocity'])[2*packet['sizes'][0]:]), initial=0.) == 0., 'perturbation changed pressure')
    passed = all_ok and divergence['pass_'] and gamma['pass_'] and helper['pass_'] and essential_ok and composition_ok and diagnostics['strong_divergence_MAX_upper'] <= DIVERGENCE_LIMIT
    return dict(status='PASS_NATIVE_PERTURBATION_CONSTRUCTION_ONLY' if passed else 'NOT_MET_NATIVE_PERTURBATION_CONSTRUCTION',
        native_construction_pass=bool(passed), norms=norms, polynomial_diagnostics=diagnostics,
        weak_divergence=divergence, synthetic_dependent_Gamma=gamma, helper_formula_comparison=helper,
        essential_increment_MAX=essential, essential_pass=essential_ok,
        synthetic_increment_composition_error=composition, synthetic_composition_pass=composition_ok,
        epsilon_gradient_convention='Gradient of r*delta_e/Ain, not grad(delta_e)/(Ain/r).',
        state_scope='Increment and synthetic Kin interface only. No developed reference, coupled histories, admissibility of a seeded physical state, or pilot qualification.')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('native', 'certificate', 'table', 'out'): p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--seconds', type=float, required=True); args = p.parse_args(); start = time.monotonic()
    out = dict(schema='WJH27IndependentPerturbationCheckV1', status='FAIL', full_coupled_gate_pass=False,
        developed_reference_qualified=False, native_construction_pass=False, script_sha256=pc.sha(__file__),
        polynomial_certificate_code_sha256=pc.sha(pc.__file__), native_sha256=pc.sha(args.native),
        certificate_sha256=pc.sha(args.certificate), reference_sha256=pc.sha(args.table),
        prospective=dict(velocity_relative_peak=1e-8, exact_Q3_representation_abs=REPRESENTATION_LIMIT,
                         divergence_abs=DIVERGENCE_LIMIT, epsilon_normalized_peak_abs=EPSILON_PEAK_LIMIT,
                         half_peak_area_minimum=.05, half_peak_area_refinements=list(pc.AREA_REFINEMENTS)))
    failure = None
    try:
        pc.require(sys.platform.startswith('linux'), 'numerical projection audit requires authoritative Linux')
        pc.require(math.isfinite(args.seconds) and args.seconds > 0, 'explicit finite positive admitted cap required')
        packet, certificate, table = [json.loads(path.read_text()) for path in (args.native, args.certificate, args.table)]
        pc.require(packet['certificate_sha256'] == out['certificate_sha256'] and packet['reference_identity'] == out['reference_sha256'] and
            certificate['reference_identity'] == out['reference_sha256'] and packet['source_identity'] == certificate['source_identity'], 'immutable source/reference/certificate binding differs')
        pc.require(certificate['status'] == 'PASS_MATHEMATICAL_SHAPE_ONLY', 'mathematical shape certificate not qualified')
        out.update(shape=certificate['shape'], layout=packet['layout'], mesh_xml_sha256=packet['mesh_xml_sha256'],
                   base_mesh_xml_sha256=certificate['mesh_xml_sha256'], source_identity=packet['source_identity'])
        out.update(audit(packet, certificate, table, start+args.seconds))
    except Exception as error:
        failure = error; out.update(exception=type(error).__name__, reason=str(error))
    out['elapsed_seconds'] = time.monotonic()-start
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x') as f: f.write(json.dumps(out, indent=2, allow_nan=False)+'\n')
    if failure is not None: raise failure
    print(json.dumps(dict(status=out['status'], shape=out['shape'], elapsed_seconds=out['elapsed_seconds'])))
    if not out['native_construction_pass']: raise SystemExit(2)


if __name__ == '__main__': main()
