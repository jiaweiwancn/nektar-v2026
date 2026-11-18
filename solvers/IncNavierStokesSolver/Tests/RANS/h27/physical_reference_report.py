#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Summarize source19 refinement and locate owner-resolved polynomial disagreement."""

import sys
import time
from collections import defaultdict

from physical_reference_refinement import OUT, ROOT, SOURCE, TABLE_PATH, audit, native, read, require


def intersections(left, right):
    i = j = 0
    while i < len(left) and j < len(right):
        a, b = left[i], right[j]
        lo, hi = max(a[0], b[0]), min(a[1], b[1])
        if hi-lo > 1e-12:
            yield a, b, lo, hi
        if a[1] <= b[1]+1e-12:
            i += 1
        if b[1] <= a[1]+1e-12:
            j += 1


def columns(oracle):
    groups = defaultdict(list)
    for element in oracle.elements:
        xl, xr, rl, rr = element['bounds']
        groups[round(xl, 12), round(xr, 12)].append((rl, rr, element))
    return [(l, r, sorted(cells)) for (l, r), cells in sorted(groups.items())]


def field_difference(first, second, deadline):
    import numpy as np
    from numpy.polynomial.legendre import leggauss
    from physical_oracle import basis

    g, gw = leggauss(8)
    g, gw = np.r_[-1., g, 1.], np.r_[0., gw, 0.]
    band_edges = [-10., -8.75, -7.50, -5., 0., 4.75, 8.75, 12.5, 15.]
    peak = {}
    k_square = defaultdict(float)
    k_area = defaultdict(float)
    for a, b, l, r in intersections(columns(first), columns(second)):
        require(time.monotonic() < deadline, 'spatial mechanism exceeded shared assessment cap')
        for e, f, z0, z1 in intersections(a[2], b[2]):
            x, z = np.meshgrid((l+r)/2+(r-l)*g/2, (z0+z1)/2+(z1-z0)*g/2)
            x, z = x.ravel(), z.ravel()
            weight = np.outer(gw, gw).ravel()*(r-l)*(z1-z0)/4
            region = 'inner' if z1 <= first.ref.H else 'outer'
            values = []
            for oracle, item in ((first, e[2]), (second, f[2])):
                vb, vx, vz = basis(item, 'u', x, z)
                kb = basis(item, 'k', x, z)[0]
                ids = np.asarray(item['u'])
                u, w = oracle.s['flow'][ids], oracle.s['flow'][oracle.n+ids]
                values.append(np.column_stack((vb@u, vb@w, (kb@oracle.s['k'][item['k']])/first.ref.Kin,
                                                vx@u, vz@u, vx@w, vz@w)))
            difference = values[1]-values[0]
            if region == 'inner':
                # Band edges are retained original x boundaries on every mesh.
                band = next(i for i, edge in enumerate(band_edges[1:]) if (l+r)/2 < edge)
                k_square[band] += float(weight@(difference[:, 2]**2))
                k_area[band] += float(sum(weight))
            for index, name in enumerate(('U', 'W', 'k', 'dUdx', 'dUdr', 'dWdx', 'dWdr')):
                delta = difference[:, index]
                q = int(np.argmax(abs(delta)))
                key = region+'/'+name
                if abs(delta[q]) > abs(peak.get(key, {}).get('second_minus_first', 0.)):
                    peak[key] = dict(x=float(x[q]), r=float(z[q]), second_minus_first=float(delta[q]),
                                     first_owner=e[2]['geometry'], second_owner=f[2]['geometry'])
    bands = [dict(x=[a, b], area=k_area[i], square=k_square[i],
                  RMS_over_Kin=(k_square[i]/k_area[i])**.5) for i, (a, b) in
             enumerate(zip(band_edges[:-1], band_edges[1:]))]
    return dict(spatial_samples='Gauss8 plus both exact owner endpoints in each common rectangle',
                sign='second minus first', peak=peak, inner_k_bands=bands)


def summarize():
    from physical_oracle import PhysicalOracle, Reference, load_native
    from resources import atomic, identity

    deadline = time.monotonic()+73
    index = read(OUT/'run-index.json')
    ref = Reference(read(TABLE_PATH))
    reports = {}
    retained_oracles = {}
    for case in index['audited']:
        independent = read(audit(case))
        data, metadata = load_native(native(case))
        require(data['source_identity'] == SOURCE and independent['native_source_identity'] == SOURCE and
                identity(native(case))['sha256'] == independent['native']['sha256'] and
                independent['state_sha256'] == data['state_sha256'], f'invalid source19 report for {case}')
        core = independent['core']
        profiles = independent['profiles']
        reports[case] = dict(native=identity(native(case)), independent=identity(audit(case)),
            native_acceptance=metadata['steady_solver']['true_norm'], state_sha256=data['state_sha256'],
            quadrature=[data['quadratureX'], data['quadratureR']], numerical_current=core['residuals'],
            numerical_limits=core['residual_limits'], equations=independent['core_steady_equations_status'],
            native_admissibility={k: data[k] for k in ('minimum_k', 'minimum_epsilon', 'minimum_covariance_over_k',
                                  'minimum_principal_over_mu', 'cutoff_radius_range')},
            boundary=independent['boundary_checks']['status'],
            principal_scope=core['principal_admissibility'], owner_flux=independent['owner_interfaces']['norms'],
            physical_development=independent['physical_development_status'], envelope=core['physical_envelope'],
            stations=profiles['stations'], station_norms=profiles['norms'],
            budget_normalized_defects={k: {metric: v[metric]/v['fixed_scale']
                 for metric in ('mixed_physical_balance', 'all_owner_gradient_balance')} for k, v in
                 core['budgets'].items() if k in ('k', 'epsilon', 'mechanical')})
        if case in ('radial', 'streamwise', 'targeted'):
            retained_oracles[case] = PhysicalOracle(data, ref, 20, deadline)
    comparisons = {}
    for path in sorted(OUT.glob('*-to-*-order20.json')):
        comparison = read(path)
        norms = comparison['norms']
        failed = sorted(k for k, v in norms.items() if v['status'] != 'MET')
        examples = {k: norms[k] for k in ('inner/k', 'inner/U', 'inner/W', 'inner/G', 'outer/epsilon',
                                         'whole/R', 'whole/Pi', 'whole/traction')}
        first, second = list(comparison['native_and_audit'])
        row = dict(archive=identity(path), status=comparison['status'], failed=failed, example_norms=examples,
                   checked=len(norms), common_owner_subcells=comparison['common_owner_subcells'],
                   budget_changes=comparison['normalized_budget_changes'])
        if 'targeted' in (first, second) and {first, second} <= retained_oracles.keys():
            row['spatial_mechanism'] = field_difference(retained_oracles[first], retained_oracles[second], deadline)
        comparisons[f'{first}-to-{second}'] = row
    atomic(OUT/'assessment.json', dict(schema='WJH27Source19ReferenceAssessmentV1', source_identity=SOURCE,
        run_index=identity(OUT/'run-index.json'), input_manifest=identity(OUT/'inputs.json'),
        localization=identity(ROOT/'science/cp2-localization-v19/owner-balance.json'),
        candidates=reports, comparisons=comparisons,
        notes='Classify using the frozen field, admissibility, physical and convergence conditions together',
        qualified_reference=None, TASK3_READY=False))
    print('SOURCE19_SPATIAL_REFINEMENT_ASSESSED', flush=True)


if __name__ == '__main__':
    require(sys.platform.startswith('linux') and ROOT.is_dir(), 'requires authoritative Linux capsule')
    sys.path[:0] = [str(ROOT/'tools'), str(ROOT/'source-v19/solvers/IncNavierStokesSolver/Tests/RANS/h27')]
    summarize()
