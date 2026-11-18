#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Owner-resolved inner k balance and retained-inlet diagnostics for H27 CP2."""
import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import resource
import sys
import time

import physical_oracle as physical  # Sets all numerical-library thread limits.
import numpy as np
from numpy.polynomial.legendre import leggauss

SOURCE = '56c35390fdf092180913f60c3ae348ca10cdcdf4e6b889c3afd589eddb511acf'
require = physical.require


def face_flux(oracle, element, side):
    """Flux along positive x or r, with each face evaluated in its owner."""
    if side < 2:
        x, r, weights = oracle.face(element, side)
    else:
        left, right, lower, upper = element['bounds']
        radius = lower if side == 2 else upper
        points, gauss = oracle.gauss
        edges = [left, *oracle.cutoff_x_splits(element, (left, right, radius, radius)), right]
        x = np.concatenate([(a+b)/2+(b-a)*points/2 for a, b in zip(edges, edges[1:])])
        weights = np.concatenate([gauss*(b-a)/2 for a, b in zip(edges, edges[1:])])
        r = np.full_like(x, radius)
    sample = oracle.sample(element, x, r)
    along = 0 if side < 2 else 1
    velocity = sample['U' if along == 0 else 'W']
    advection = float(weights@(velocity*sample['k']))
    diffusion = float(-weights@(sample['Dk']*sample['gradk'][:, along]))
    return np.array([advection, diffusion])


def field_measure(oracle, element, lower, upper):
    """Measure signed error and production on one original material band."""
    left, right = element['bounds'][:2]
    ref = oracle.ref
    x, r, weights = oracle.rule(element, (left, right, lower, upper),
                                extra_cutoff=lambda x: np.full_like(x, ref.old_c))
    require(len(weights) > 0 and np.all(weights > 0), 'nonpositive inner weights')
    sample = oracle.sample(element, x, r)
    old = ref.inlet(r)
    error = (sample['k']-old['k'])/ref.Kin
    return np.array([np.sum(weights), weights@error, weights@(error*error),
                     weights@sample['P'], weights@sample['E']])


def extrema(oracle, element):
    """Include exact original owners and both material/cutoff endpoints."""
    left, right, lower, upper = element['bounds']
    ref = oracle.ref
    nodes, _ = oracle.gauss
    xs = np.r_[(left+right)/2+(right-left)*nodes/2, left, (left+right)/2, right]
    cutoff = ref.cutoff(oracle.trace(element, 'k', xs))
    sample_x, sample_r = [], []
    for point, cut in zip(xs, cutoff):
        cuts = oracle.cuts(lower, upper, (float(cut), ref.old_c))
        radii = np.asarray(cuts+[(a+b)/2 for a, b in zip(cuts, cuts[1:])])
        sample_x.extend(np.full_like(radii, point))
        sample_r.extend(radii)
    x, r = np.asarray(sample_x), np.asarray(sample_r)
    sample = oracle.sample(element, x, r)
    old = ref.inlet(r)['k']
    error = (sample['k']-old)/ref.Kin
    j = int(np.argmax(abs(error)))
    out = dict(x=float(x[j]), r=float(r[j]), owner=element['geometry'],
               k=float(sample['k'][j]), historical_k=float(old[j]),
               signed_over_Kin=float(error[j]))
    return out, float(np.min(sample['k']))


def measure(data, reference, order, deadline, sample_extrema):
    oracle = physical.PhysicalOracle(data, reference, order, deadline)
    bands = sorted((0., reference.old_c, .4*reference.H, reference.h,
                    .75*reference.H, reference.H))
    columns = defaultdict(list)
    for element in oracle.elements:
        if not element['epsilon']:
            columns[tuple(round(x, 12) for x in element['bounds'][:2])].append(element)
    require(columns, 'no native inner columns')
    all_rows, cell_rows = [], []
    for (left, right), cells in sorted(columns.items()):
        cells.sort(key=lambda element: element['bounds'][2])
        integrated = np.zeros((len(bands)-1, 5))
        x_minus = np.zeros(2)
        x_plus = np.zeros(2)
        jumps = np.zeros(2)
        abs_jumps = np.zeros(2)
        previous = None
        peak = None
        for element in cells:
            oracle.clock()
            lower, upper = element['bounds'][2:]
            row = np.zeros(5)
            for i, (a, b) in enumerate(zip(bands, bands[1:])):
                low, high = max(a, lower), min(b, upper)
                if high <= low:
                    continue
                terms = field_measure(oracle, element, low, high)
                integrated[i] += terms
                row += terms
            left_flux, right_flux, lower_flux, upper_flux = (
                face_flux(oracle, element, face) for face in range(4))
            x_minus += left_flux
            x_plus += right_flux
            if previous is None:
                ground = lower_flux
            else:
                jump = previous-lower_flux
                jumps += jump
                abs_jumps += abs(jump)
            previous = upper_flux
            p_minus_e = row[3]-row[4]
            owner_defect = float(np.sum(right_flux-left_flux+upper_flux-lower_flux)-p_minus_e)
            if sample_extrema:
                witness, minimum = extrema(oracle, element)
                if peak is None or abs(witness['signed_over_Kin']) > abs(peak['signed_over_Kin']):
                    peak = witness
                cell_rows.append(dict(owner=element['geometry'], bounds=element['bounds'],
                                      k_departure_integral_squared=float(row[2]), sampled_peak=witness,
                                      min_k_sampled=minimum, owner_integrated_defect=owner_defect))
        inner_H = previous
        outer_H = face_flux(oracle, oracle.owner[cells[-1]['geometry']], 2)
        area, signed, square, production, sink = np.sum(integrated, axis=0)
        require(abs(area-(right-left)*reference.H) <= 1e-11*area, 'column area does not close')
        width = right-left
        scale = width*reference.A
        outward = x_plus-x_minus+inner_H-ground
        terms = dict(
            x_left=x_minus.tolist(), x_right=x_plus.tolist(),
            x_net=(x_plus-x_minus).tolist(), ground_up=ground.tolist(),
            H_inner_up=inner_H.tolist(), H_outer_up=outer_H.tolist(),
            H_jump=(inner_H-outer_H).tolist(), interior_radial_jump=jumps.tolist(),
            interior_radial_absolute_jump=abs_jumps.tolist())
        # Signed owner residual equals the sum of each cell's integrated
        # flux divergence. Surface-only imbalance leaves out the internal
        # numerical flux jumps and therefore has a different interpretation.
        all_rows.append(dict(
            x=[left, right], cells=len(cells), area=float(area),
            signed_mean_over_Kin=float(signed/area), RMS_over_Kin=math.sqrt(square/area),
            integral_square=float(square), sampled_peak=peak,
            production=float(production), dependent_destruction=float(sink),
            production_minus_destruction=float(production-sink), flux=terms,
            surface_defect_over_A_width=float((np.sum(outward)-production+sink)/scale),
            owner_divergence_defect_over_A_width=float((np.sum(outward+jumps)-production+sink)/scale),
            radial_bands=[dict(r=[a, b], terms=list(map(float, row)))
                          for a, b, row in zip(bands, bands[1:], integrated)]))
    area = sum(row['area'] for row in all_rows)
    square = sum(row['integral_square'] for row in all_rows)
    for row in all_rows:
        row['fraction_global_squared_departure'] = row['integral_square']/square
    return dict(order=order, inner_area=area, inner_RMS_over_Kin=math.sqrt(square/area),
                inner_production=sum(row['production'] for row in all_rows),
                inner_dependent_destruction=sum(row['dependent_destruction'] for row in all_rows),
                columns=all_rows, original_cell_indicators=cell_rows,
                band_term_order=['area', 'integral_signed_k_departure', 'integral_squared_k_departure',
                                 'production', 'dependent_destruction'], face_flux_order=['advection', 'diffusion'])


def parallel_inlet(ref):
    """Use the retained inlet as parallel data under both scalar laws."""
    current_cutoff = float(ref.cutoff(ref.Kin))
    def evaluate(radius):
        inlet = ref.inlet(radius)
        k, gradient = inlet['k'], inlet['G']
        aux = ref.aux(np.full_like(radius, ref.Kin), radius)
        new_E = ref.A/np.maximum(radius, current_cutoff)
        old_E = ref.A/np.maximum(radius, ref.old_c)
        old_D = ref.nu+aux
        new_D = ref.nu+ref.cd*k*k/(ref.sk*new_E)
        old_P = aux*gradient*gradient
        chi = np.clip((radius/ref.H-.4)/.35, 0., 1.)
        chi = chi**3*(10-15*chi+6*chi**2)
        active = chi > 0
        new_P = (1-chi)*old_P
        G = np.zeros((int(active.sum()), 3, 3))
        G[:, 0, 2] = gradient[active]
        if np.any(active):
            stress = physical.full_wj(k[active], new_E[active], G)[0]
            new_P[active] -= chi[active]*stress[:, 0, 2]*gradient[active]
        derivative = ref.old_Q(radius)/old_D
        return old_P, new_P, old_E, new_E, old_D, new_D, derivative

    splits = sorted((0., current_cutoff, ref.old_c, .4*ref.H, .75*ref.H, ref.H))
    results = []
    for order in (32, 48):
        a, w = leggauss(order)
        integrated = np.zeros(4)
        rows = []
        for lower, upper in zip(splits, splits[1:]):
            radius = (lower+upper)/2+(upper-lower)*a/2
            terms = evaluate(radius)
            integrated += (upper-lower)/2*(np.vstack(terms[:4])@w)
            point = evaluate(np.array([upper]))
            actual_Q = float(point[5][0]*point[6][0])
            needed_Q = float(integrated[3]-integrated[1])
            rows.append(dict(radius=upper, old_cumulative_Q_over_A=float((integrated[2]-integrated[0])/ref.A),
                             v2_required_Q_over_A=needed_Q/ref.A, actual_v2_gradient_flux_over_A=actual_Q/ref.A,
                             gradient_minus_required_over_A=(actual_Q-needed_Q)/ref.A))
        results.append(dict(order=order, cuts=rows, integrated_over_A=(integrated/ref.A).tolist()))
    r = np.unique(np.r_[np.linspace(0., ref.H, 257), splits])
    terms = evaluate(r)
    ratio = terms[5]/terms[4]
    return dict(old_cutoff=ref.old_c, v2_cutoff_at_Kin=current_cutoff,
                old_ground_destruction=float(terms[2][0]), v2_ground_destruction=float(terms[3][0]),
                ground_v2_over_historical_diffusion=float(ratio[0]),
                v2_over_historical_diffusion_range=[float(np.min(ratio)), float(np.max(ratio))],
                integrated_term_order=['historical_P', 'v2_P', 'historical_E', 'v2_E'],
                independent_rules=results)


def assess(native, table, audit, deadline):
    data, metadata = physical.load_native(native)
    reference = physical.Reference(json.loads(table.read_text()))
    recorded = json.loads(audit.read_text())
    require(data['source_identity'] == SOURCE and data['step'] == data['time'] == 0
            and metadata.get('steady_mode') is True and metadata['steady_solver']['converged'],
            'localization needs a converged source19 normal-main steady result')
    require(recorded['native']['sha256'] == physical.identity(native)['sha256']
            and recorded['source']['sha256'] == physical.identity(physical.__file__)['sha256']
            and recorded['table']['sha256'] == physical.TABLE_SHA
            and recorded['core_steady_equations_status'] == recorded['boundary_checks']['status'] == 'MET',
            'physical oracle identity, equations or boundaries differ')
    expected = recorded['core']['physical_envelope']['domain/inner/k']
    bulk = recorded['core']['bulk']['inner']
    bound = dict(order=recorded['quadrature'], state=recorded['state_sha256'],
                 inner_k=expected, inner_production=bulk['P'], inner_destruction=bulk['E'],
                 mass=recorded['core']['budgets']['mass'],
                 k=recorded['core']['budgets']['k'],
                 H_owner_norms=recorded['owner_interfaces']['norms'],
                 station_inner_k={k: v for k, v in recorded['profiles']['norms'].items()
                                  if k.endswith('/inner/k')})
    del recorded
    coarse = measure(data, reference, 20, deadline, sample_extrema=False)
    fine = measure(data, reference, 32, deadline, sample_extrema=True)
    errors = dict(RMS_orders=abs(coarse['inner_RMS_over_Kin']-fine['inner_RMS_over_Kin']),
                  RMS_bound=abs(expected['RMS']-fine['inner_RMS_over_Kin']),
                  production_over_A_length=abs(bulk['P']-fine['inner_production'])/(25*reference.A),
                  destruction_over_A_length=abs(bulk['E']-fine['inner_dependent_destruction'])/(25*reference.A))
    require(max(errors.values()) <= 1e-5, 'independent fixed-state global integration differs')
    by_column = []
    for a, b in zip(coarse['columns'], fine['columns']):
        require(a['x'] == b['x'], 'quadrature changed original owners')
        width = a['x'][1]-a['x'][0]
        delta = {k: abs(a[k]-b[k]) for k in (
            'RMS_over_Kin', 'surface_defect_over_A_width', 'owner_divergence_defect_over_A_width')}
        delta.update({k+'_over_A_width': abs(a[k]-b[k])/(width*reference.A)
                      for k in ('production', 'dependent_destruction')})
        by_column.append(dict(x=a['x'], errors=delta))
    peak = max((r['sampled_peak'] for r in fine['columns']), key=lambda row: abs(row['signed_over_Kin']))
    return dict(status='MEASURED_FIXED_SOURCE19_REFERENCE_LOCALIZATION',
                source_identity=SOURCE, fixed_state_reference=bound, independent_rules=[coarse, fine],
                independent_global_differences=errors, independent_column_differences=by_column,
                sampled_peak=peak, parallel_inlet=parallel_inlet(reference),
                scope='Fixed coefficient integration and owner flux localization. Physical mesh convergence is assessed separately.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('native', 'table', 'audit', 'out'):
        parser.add_argument('--'+name, required=True, type=Path)
    parser.add_argument('--seconds', required=True, type=float)
    args = parser.parse_args()
    require(sys.platform.startswith('linux') and math.isfinite(args.seconds) and 0 < args.seconds <= 360,
            'requires a bounded Linux operation')
    require(not args.out.exists(), 'immutable output exists')
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    cap = min([8*1024**3]+[v for v in (soft, hard) if v != resource.RLIM_INFINITY])
    resource.setrlimit(resource.RLIMIT_AS, (cap, hard))
    start = time.monotonic()
    report = dict(schema='WJH27PhysicalReferenceLocalizationV1', native=physical.identity(args.native),
                  audit=physical.identity(args.audit), table=physical.identity(args.table),
                  checker=physical.identity(__file__), oracle=physical.identity(physical.__file__),
                  TASK3_READY=False, reference_qualified=False)
    error = None
    try:
        require(report['table']['sha256'] == physical.TABLE_SHA, 'wrong approved table')
        report.update(assess(args.native, args.table, args.audit, start+args.seconds))
    except Exception as exc:
        error = exc
        report.update(status='FAIL', reason=str(exc), exception=type(exc).__name__)
    report['elapsed_seconds'] = time.monotonic()-start
    with args.out.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({k: report[k] for k in ('status', 'elapsed_seconds', 'reason') if k in report}), flush=True)
    if error:
        raise SystemExit(2)


if __name__ == '__main__':
    main()
