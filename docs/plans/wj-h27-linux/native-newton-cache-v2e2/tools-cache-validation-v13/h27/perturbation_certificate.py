#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Linux-only independent polynomial perturbation peak/support certificate.

Initial polynomial coefficients are exact rational numbers representing the
exported binary64 knots. Every subsequent Bernstein operation encloses roundoff
with directed nextafter bounds. No sampled maximum is used as an upper bound.
"""
import argparse
from fractions import Fraction as F
import hashlib
import heapq
import itertools
import json
import math
import os
from pathlib import Path
import sys
import time

sys.dont_write_bytecode = True
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS'):
    os.environ[name] = '1'
import numpy as np

PEAK_GAP = 5e-9
AREA_REFINEMENTS = (.01, .003, .001)
MIN_HALF_AREA = .05
FIXED_X = [-2.5, -1.25, 0., 1.25, 2.5]


def require(ok, message):
    if not ok: raise ValueError(message)


def clock(deadline):
    if time.monotonic() >= deadline: raise TimeoutError('perturbation certificate reached its admitted cap')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rational(value):
    return F.from_float(float(value))


def directed(value, upper):
    """Exact rational to outward binary64, including nearest conversion error."""
    candidate = float(value)
    require(math.isfinite(candidate), 'nonfinite interval conversion')
    return float(np.nextafter(candidate, math.inf if upper else -math.inf))


def sqrt_bounds(value):
    require(value >= 0, 'negative squared norm')
    if value == 0: return [0., 0.]
    bits = 180
    integer = math.isqrt((value.numerator << (2*bits))//value.denominator)
    return [max(0., directed(F(integer, 1 << bits), False)),
            directed(F(integer+1, 1 << bits), True)]


def add(a, b):
    out = [F(0)]*max(len(a), len(b))
    for i, value in enumerate(a): out[i] += value
    for i, value in enumerate(b): out[i] += value
    return out


def mul(a, b):
    out = [F(0)]*(len(a)+len(b)-1)
    for i, u in enumerate(a):
        for j, v in enumerate(b): out[i+j] += u*v
    return out


def derivative(a, width=F(1)):
    return [i*a[i]/width for i in range(1, len(a))] or [F(0)]


def spline_piece(knots, left, right):
    """Cox-de Boor polynomial in t on a single physical interval."""
    knots = list(map(rational, knots)); left, right = rational(left), rational(right)
    midpoint = (left+right)/2
    def recurse(i, degree):
        if degree == 0: return [F(int(knots[i] <= midpoint < knots[i+1]))]
        a = [(left-knots[i])/(knots[i+degree]-knots[i]),
             (right-left)/(knots[i+degree]-knots[i])]
        b = [(knots[i+degree+1]-left)/(knots[i+degree+1]-knots[i+1]),
             -(right-left)/(knots[i+degree+1]-knots[i+1])]
        return add(mul(a, recurse(i, degree-1)), mul(b, recurse(i+1, degree-1)))
    return recurse(0, 3)


def hermite_piece(left, right, lower, peak, upper):
    l, r, a, p, b = map(rational, (left, right, lower, peak, upper))
    midpoint = (l+r)/2
    if not a < midpoint < b: return [F(0)]
    t = [(l-a)/(p-a), (r-l)/(p-a)] if midpoint <= p else [(b-l)/(b-p), -(r-l)/(b-p)]
    return mul(mul(t, t), add([F(3)], [-2*v for v in t]))


def outer(a, b):
    return [[u*v for v in b] for u in a]


def padd(a, b):
    out = [[F(0)]*max(len(a[0]), len(b[0])) for _ in range(max(len(a), len(b)))]
    for matrix in (a, b):
        for i, row in enumerate(matrix):
            for j, value in enumerate(row): out[i][j] += value
    return out


def psquare(a):
    out = [[F(0)]*(2*len(a[0])-1) for _ in range(2*len(a)-1)]
    for i, row in enumerate(a):
        for j, u in enumerate(row):
            for k, line in enumerate(a):
                for l, v in enumerate(line): out[i+k][j+l] += u*v
    return out


def bernstein_exact(power):
    nx, nr = len(power)-1, len(power[0])-1
    intermediate = [[sum((power[a][b]*F(math.comb(i, a), math.comb(nx, a)) for a in range(i+1)), F(0))
                     for b in range(nr+1)] for i in range(nx+1)]
    return [[sum((intermediate[i][b]*F(math.comb(j, b), math.comb(nr, b)) for b in range(j+1)), F(0))
             for j in range(nr+1)] for i in range(nx+1)]


def interval_coefficients(power):
    exact = bernstein_exact(power)
    return (np.asarray([[directed(v, False) for v in row] for row in exact]),
            np.asarray([[directed(v, True) for v in row] for row in exact]))


def average(a, b, upper):
    direction = math.inf if upper else -math.inf
    return np.nextafter(np.nextafter(a+b, direction)*.5, direction)


def split_coefficients(lower, upper, axis):
    def split(matrix, high):
        current = np.moveaxis(matrix, axis, 0).copy(); n = len(current)
        left, right = [current[0].copy()], [current[-1].copy()]
        for _ in range(1, n):
            current = average(current[:-1], current[1:], high)
            left.append(current[0].copy()); right.append(current[-1].copy())
        return np.moveaxis(np.asarray(left), 0, axis), np.moveaxis(np.asarray(right[::-1]), 0, axis)
    ll, lr = split(lower, False); ul, ur = split(upper, True)
    return ((ll, ul), (lr, ur))


def center_lower(coefficients):
    value = coefficients.copy()
    for axis in (0, 1):
        value = np.moveaxis(value, axis, 0)
        while len(value) > 1: value = average(value[:-1], value[1:], False)
        value = np.moveaxis(value, 0, axis)
    return float(value[0, 0])


def split_axis(lower, upper):
    # Only efficiency depends on this choice; both children always enclose.
    mid = .5*(lower+upper)
    variation = [float(np.max(abs(np.diff(mid, axis=i)), initial=0.)) for i in (0, 1)]
    return int(variation[1] > variation[0])


def maximum(pieces, deadline, relative_gap=PEAK_GAP):
    heap = []; serial = itertools.count(); lower_bound = 0.; operations = 0
    def insert(lo, hi):
        nonlocal lower_bound
        lower_bound = max(lower_bound, center_lower(lo), float(lo[0, 0]), float(lo[-1, 0]),
                          float(lo[0, -1]), float(lo[-1, -1]))
        bound = float(np.max(hi))
        if bound > lower_bound: heapq.heappush(heap, (-bound, next(serial), lo, hi))
    for piece in pieces:
        clock(deadline); insert(*interval_coefficients(piece['power']))
    while heap:
        clock(deadline); upper_bound = max(lower_bound, -heap[0][0])
        if lower_bound > 0 and (upper_bound-lower_bound)/lower_bound <= relative_gap: break
        neg, _, lo, hi = heapq.heappop(heap)
        if -neg <= lower_bound: continue
        for child in split_coefficients(lo, hi, split_axis(lo, hi)): insert(*child)
        operations += 1
    upper_bound = max(lower_bound, -heap[0][0] if heap else lower_bound)
    return dict(lower=lower_bound, upper=upper_bound, subdivisions=operations, active_boxes=len(heap),
                relative_gap=(upper_bound-lower_bound)/lower_bound)


def enclosed_half_area(pieces, peak_lower, peak_upper, deadline):
    """Area(|field| >= actual_global_peak/2), using squared polynomials."""
    lower_threshold = rational(peak_lower)/4; upper_threshold = rational(peak_upper)/4
    heap = []; serial = itertools.count(); inside = F(0); unknown = F(0); operations = 0
    support = sum((piece['area'] for piece in pieces), F(0)); levels = []
    def insert(lo, hi, measure):
        nonlocal inside, unknown
        # Compare exact binary64 bounds with exact rational thresholds.
        if rational(float(np.min(lo))) >= upper_threshold: inside += measure
        elif rational(float(np.max(hi))) < lower_threshold: pass
        else:
            unknown += measure
            heapq.heappush(heap, (-measure, next(serial), lo, hi, measure))
    for piece in pieces:
        clock(deadline); insert(*interval_coefficients(piece['power']), piece['area'])
    for fraction in AREA_REFINEMENTS:
        while unknown > rational(fraction)*support:
            clock(deadline); require(heap, 'lost undecided half-peak boxes')
            _, _, lo, hi, measure = heapq.heappop(heap); unknown -= measure
            for child in split_coefficients(lo, hi, split_axis(lo, hi)): insert(*child, measure/2)
            operations += 1
        levels.append(dict(relative_support_uncertainty_target=fraction,
            lower=max(0., directed(inside, False)), upper=directed(inside+unknown, True),
            unresolved_area_upper=directed(unknown, True), subdivisions=operations, active_boxes=len(heap)))
    return levels


def integral(piece):
    return piece['area']*sum((v/F((i+1)*(j+1)) for i, row in enumerate(piece['power']) for j, v in enumerate(row)), F(0))


def bound_controls(deadline):
    # Independent closed-form answers, including a non-grid interior maximum
    # and a half-level area whose boundary is not dyadic.
    interior = [[F(1)-F(37, 100)**2-F(61, 100)**2, F(122, 100), F(-1)],
                [F(74, 100), F(0), F(0)], [F(-1), F(0), F(0)]]
    peak = maximum([dict(power=interior, area=F(1))], deadline)
    require(peak['lower'] <= 1 <= peak['upper'], 'interior exact maximum escaped Bernstein enclosure')
    edge = maximum([dict(power=[[F(0), F(0), F(1)], [F(0), F(0), F(0)],
                               [F(1), F(0), F(0)]], area=F(1))], deadline)
    require(edge['lower'] <= 2 <= edge['upper'], 'edge exact maximum escaped Bernstein enclosure')
    piece = dict(power=[[F(9, 100)], [F(-3, 5)], [F(1)]], area=F(1))
    # Supply outward bounds for the exact rational squared maximum 49/100.
    area = enclosed_half_area([piece], directed(F(49, 100), False), directed(F(49, 100), True), deadline)
    require(rational(area[-1]['lower']) <= F(7, 20) <= rational(area[-1]['upper']), 'exact half-level area escaped enclosure')
    require(integral(piece) == F(37, 300), 'exact polynomial moment differs')
    return dict(interior_maximum=peak, edge_maximum=edge, half_area=area, exact_energy='37/300', status='PASS')


def shape_pieces(layout, field):
    xk, rk, H = layout['x_knots'], layout['r_knots'], layout['H']
    rs = rk if field == 'velocity' else ([rk[0], H, rk[-1]] if field == 'k' else [H, layout['epsilon_peak_radius'], rk[-1]])
    result = []
    for xl, xr in zip(xk[:-1], xk[1:]):
        bx = spline_piece(xk, xl, xr); dx = derivative(bx, rational(xr)-rational(xl))
        for rl, rr in zip(rs[:-1], rs[1:]):
            if field == 'velocity':
                br = spline_piece(rk, rl, rr); dr = derivative(br, rational(rr)-rational(rl))
                power = padd(outer(mul(bx, bx), mul(dr, dr)), outer(mul(dx, dx), mul(br, br)))
            else:
                br = hermite_piece(rl, rr, rs[0], rs[1], rs[2])
                power = psquare(outer([F(3, 2)*v for v in bx], br))
            result.append(dict(bounds=[xl, xr, rl, rr], power=power,
                area=(rational(xr)-rational(xl))*(rational(rr)-rational(rl))))
    return result


def validate_inventory(packet, required_mode='layout'):
    require(packet['schema'] == 'WJH27PerturbationNativeProbeV1' and packet['mode'] == required_mode, 'wrong native geometry export')
    require((len(packet['elements']) == 700 if required_mode == 'layout' else len(packet['elements']) >= 700) and
            packet['tensor_inventory_verified'] is True, 'missing retained native tensor inventory')
    xs, rs = packet['x_edges'], packet['r_edges']; cells = set()
    require(all(math.isfinite(v) for v in xs+rs) and all(a < b for a, b in zip(xs[:-1], xs[1:])) and
            all(a < b for a, b in zip(rs[:-1], rs[1:])), 'invalid native edge ordering')
    require(abs(xs[0]+10) < 1e-12 and abs(xs[-1]-15) < 1e-12 and abs(rs[0]) < 1e-12 and
            abs(rs[-1]-19.9125) < 1e-12 and packet['H'] == .15, 'retained physical domain differs')
    for e in packet['elements']:
        bounds = e['bounds']; indices = []
        for edges, value in ((xs, bounds[0]), (xs, bounds[1]), (rs, bounds[2]), (rs, bounds[3])):
            matches = [i for i, edge in enumerate(edges) if abs(edge-value) < 1e-12]
            require(len(matches) == 1, 'ambiguous native edge'); indices.append(matches[0])
        require(indices[1] == indices[0]+1 and indices[3] == indices[2]+1, 'native cell skips tensor layer')
        require((indices[0], indices[2]) not in cells, 'duplicate native tensor cell')
        cells.add((indices[0], indices[2]))
    require(len(cells) == (len(xs)-1)*(len(rs)-1), 'native tensor inventory is incomplete')
    return xs, rs


def independent_layout(xs, rs, H, shape):
    for x in FIXED_X: require(any(abs(v-x) < 1e-12 for v in xs), 'fixed x knot absent')
    candidates = [r for r in rs if .075 <= r <= .45]
    require(len(candidates) >= 5 and candidates[-1]-candidates[0] >= .3, 'unsupported radial layout')
    a, b = candidates[0], candidates[-1]
    def select(target, exclude=None):
        rows = []
        for knots in itertools.combinations(candidates[1:-1], 3):
            trial = [a, *knots, b]
            if trial == exclude: continue
            rows.append((sum(((v-a)/(b-a)-t)**2 for v, t in zip(knots, target)), trial))
        require(rows, 'no distinct secondary shape'); return min(rows, key=lambda v: v[0])[1]
    primary = select([.25, .5, .75])
    radial = primary if shape == 'PrimaryQuartile' else select([.2, .4, .65], primary)
    fraction = .5 if shape == 'PrimaryQuartile' else .65
    eligible = [r for r in rs if H+1e-12 < r < b-1e-12]
    peak = min(eligible, key=lambda r: abs(r-(H+fraction*(b-H))))
    return dict(shape=shape, x_knots=FIXED_X, r_knots=radial, H=H, epsilon_peak_radius=peak)


def certificate(packet, shape, deadline):
    xs, rs = validate_inventory(packet); layout = packet['layouts'][shape]
    require(layout == independent_layout(xs, rs, packet['H'], shape), 'native layout differs from independent selection')
    def index(edges, value):
        found = [i for i, v in enumerate(edges) if abs(v-value) < 1e-12]
        require(len(found) == 1, 'shape knot is not a unique native edge'); return found[0]
    counts = [index(xs, layout['x_knots'][-1])-index(xs, layout['x_knots'][0]),
              index(rs, layout['r_knots'][-1])-index(rs, layout['r_knots'][0]),
              index(rs, layout['r_knots'][-1])-index(rs, layout['H'])]
    require(min(counts) >= 4 and (layout['x_knots'][-1]-layout['x_knots'][0])*(layout['r_knots'][-1]-layout['r_knots'][0]) >= 1.5,
            'support span or native interval minimum not met')
    domain = (rational(xs[-1])-rational(xs[0]))*(rational(rs[-1])-rational(rs[0]))
    outer_domain = (rational(xs[-1])-rational(xs[0]))*(rational(rs[-1])-rational(layout['H']))
    result = {}; normalization = None
    for field in ('velocity', 'k', 'epsilon'):
        pieces = shape_pieces(layout, field); support = sum((p['area'] for p in pieces), F(0))
        peak = maximum(pieces, deadline); bounds = [sqrt_bounds(rational(peak['lower']))[0], sqrt_bounds(rational(peak['upper']))[1]]
        estimate = (bounds[0]+bounds[1])/2
        if field == 'velocity':
            normalization = dict(layout=layout, unitPeakLower=bounds[0], unitPeakEstimate=estimate, unitPeakUpper=bounds[1])
            scale = rational(.001)/rational(estimate)
            require((bounds[1]-bounds[0])/estimate <= 1e-8, 'outward peak bracket exceeds native contract')
        else: scale = rational(.01)
        energy = sum((integral(p) for p in pieces), F(0))*scale*scale
        rms_domain = sqrt_bounds(energy/(outer_domain if field == 'epsilon' else domain))
        rms_support = sqrt_bounds(energy/support)
        areas = enclosed_half_area(pieces, peak['lower'], peak['upper'], deadline)
        result[field] = dict(unit_squared_peak=peak, unit_peak_enclosure=bounds,
            normalized_peak_enclosure=[directed(rational(bounds[0])*scale, False), directed(rational(bounds[1])*scale, True)],
            energy_enclosure=[directed(energy, False), directed(energy, True)],
            domain_RMS_enclosure=rms_domain, support_RMS_enclosure=rms_support,
            support_area_enclosure=[directed(support, False), directed(support, True)], half_peak_area_refinement=areas,
            domain_RMS_pass=rms_domain[0] >= (1e-5 if field == 'velocity' else 1e-4),
            support_RMS_pass=rms_support[0] >= (1e-4 if field == 'velocity' else .001),
            half_peak_area_pass=areas[-1]['lower'] >= MIN_HALF_AREA)
    return normalization, result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh-native', type=Path, required=True)
    p.add_argument('--shape', choices=('PrimaryQuartile', 'SecondaryLowerSkew'), required=True)
    p.add_argument('--seconds', type=float, required=True)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args(); start = time.monotonic()
    out = dict(schema='WJH27PerturbationPeakCertificateV1', status='FAIL', full_coupled_gate_pass=False,
        developed_reference_qualified=False, native_projected_shape_qualified=False,
        script_sha256=sha(__file__), layout_export_sha256=sha(args.mesh_native), shape=args.shape,
        arithmetic='Exact rational polynomial construction from binary64 knots; outward-rounded interval Bernstein de Casteljau; 180-bit integer sqrt enclosure and outward binary64 serialization.',
        peak_method='Global Bernstein upper bounds and enclosed attained-point lower bounds, never a sampled global upper bound.',
        half_area_relative_support_uncertainty_levels=list(AREA_REFINEMENTS), normalized_epsilon_projection_peak_tolerance=1e-5)
    failure = None
    try:
        require(sys.platform.startswith('linux'), 'numerical certificate requires authoritative Linux')
        require(math.isfinite(args.seconds) and args.seconds > 0, 'explicit finite positive admitted cap required')
        packet = json.loads(args.mesh_native.read_text())
        out['closed_form_algorithm_controls'] = bound_controls(start+args.seconds)
        normalization, norms = certificate(packet, args.shape, start+args.seconds)
        out.update(normalization=normalization, target_shape_norms=norms, mesh_xml_sha256=packet['mesh_xml_sha256'],
                   source_identity=packet['source_identity'], reference_identity=packet['reference_identity'],
                   frozen_layout=normalization['layout'], domain=packet['domain'], x_edges=packet['x_edges'], r_edges=packet['r_edges'])
        passed = all(row[key] for row in norms.values() for key in ('domain_RMS_pass', 'support_RMS_pass', 'half_peak_area_pass'))
        out['status'] = 'PASS_MATHEMATICAL_SHAPE_ONLY' if passed else 'NOT_MET_SHAPE_NORMS'
    except Exception as error:
        failure = error; out.update(exception=type(error).__name__, reason=str(error))
    out['elapsed_seconds'] = time.monotonic()-start
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x') as f: f.write(json.dumps(out, indent=2, allow_nan=False)+'\n')
    if failure is not None: raise failure
    print(json.dumps(dict(status=out['status'], shape=args.shape, elapsed_seconds=out['elapsed_seconds'])))
    if out['status'] != 'PASS_MATHEMATICAL_SHAPE_ONLY': raise SystemExit(2)


if __name__ == '__main__': main()
