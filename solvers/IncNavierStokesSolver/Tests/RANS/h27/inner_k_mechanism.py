#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Fixed-profile, four-choice diagnosis of the v2 inner-k tendency."""

import argparse
import math
import json
from pathlib import Path
import sys
import time

from physical_oracle import Reference, TABLE_SHA, full_wj, identity, require
import mpmath as mp
import numpy as np
from numpy.polynomial.legendre import leggauss


ROWS = {'A': ('scalar', 'v2'), 'B': ('auxiliary', 'v2'),
        'C': ('scalar', 'historical'), 'D': ('auxiliary', 'historical')}
TOL = 1e-10
LOCALIZATION_SHA = '3e0b8e3e717c7413e404c3af2e34efc3906d660f5f88790833691fe6127cac19'
OLD_COLUMNS_SHA = '7b48adfd27add3aa685602232e1d7d32ae0e815f4c99b98d561cdf49d8b30a3b'
STENCIL = (-147, 360, -450, 400, -225, 72, -10)


def read(path):
    return json.loads(Path(path).read_text())


class FixedProfile:
    def __init__(self, ref):
        self.ref = ref
        self.cutoffs = dict(v2=float(ref.cutoff(ref.Kin)), historical=ref.old_c)
        self.material = sorted({0., *self.cutoffs.values(), .4*ref.H, ref.h, .75*ref.H, ref.H})
        first = self.material[1]
        self.integration = sorted({*self.material, *(first*f for f in (.125, .25, .5, .75))})

    def evaluate(self, r, side=1):
        ref = self.ref
        inlet = ref.inlet(r)
        k, g, u = (inlet[name] for name in ('k', 'G', 'U'))
        aux = ref.a*(r-ref.zv)
        da = ref.nu+aux
        pa = aux*g*g
        ea = ref.A/np.maximum(r, ref.old_c)
        qa = ref.old_Q(r)
        kr = qa/da
        krr = (ea-pa-ref.a*kr)/da
        chi = np.clip((r/ref.H-.4)/.35, 0., 1.)
        chi = chi**3*(10-15*chi+6*chi*chi)
        collar = np.zeros_like(r)
        active = chi > 0
        if np.any(active):
            G = np.zeros((int(active.sum()), 3, 3))
            G[:, 0, 2] = g[active]
            # Both cutoffs lie strictly below the collar. E there is A/r.
            R = full_wj(k[active], ref.A/r[active], G)[0]
            collar[active] = chi[active]*(-R[:, 0, 2]*g[active]-pa[active])
        p = pa+collar
        rows = {}
        for name, (law, choice) in ROWS.items():
            c = self.cutoffs[choice]
            e = ref.A/np.maximum(r, c)
            tail = (r > c) | ((r == c) & (side > 0))
            er = np.where(tail, -ref.A/np.maximum(r, c)**2, 0.)
            if law == 'scalar':
                d = ref.nu+ref.cd*k*k/(ref.sk*e)
                dr = ref.cd/ref.sk*(2*k*kr/e-k*k*er/(e*e))
            else:
                d, dr = da, np.full_like(r, ref.a)
            div = d*krr+dr*kr
            rho = d/da
            rho_r = (dr-rho*ref.a)/da
            rows[name] = dict(D=d, D_r=dr, P=p, E=e, Q=d*kr, diffusion_divergence=div,
                              source=p-e, rhs=div+p-e, diffusion_ratio=rho,
                              identity=(rho-1)*(ea-pa)+rho_r*qa+(ea-e)+collar)
        return dict(U=u, U_r=g, k=k, k_r=kr, k_rr=krr,
                    P_A=pa, E_A=ea, Q_A=qa, D_A=da, collar=collar), rows


def derivative_check(profile, deadline):
    """Differentiate high-precision reconstructed fluxes on smooth sides."""
    ref = profile.ref
    with mp.workdps(50):
        h, nu, K, a, d, strength = map(mp.mpf, (ref.H, ref.nu, ref.Kin, ref.a, ref.D0, ref.A))
        cd, sigma = map(mp.mpf, (ref.cd, ref.sk))

        def j(x):
            return mp.log1p(a*x/d)+nu/(d+a*x)-nu/d

        ca = h*mp.exp(1-j(h))
        cv = mp.mpf(profile.cutoffs['v2'])

        def qold(x):
            return strength*((x/ca if x <= ca else 1+mp.log(x/ca))-j(x))

        def kr(x):
            return qold(x)/(d+a*x)

        def values(x):
            if time.monotonic() >= deadline:
                raise TimeoutError('fixed profile differentiation exhausted its allowance')
            k = K-mp.quad(kr, [x, ca, h] if x < ca else [x, h])
            flux = dict(k_r=kr(x))
            for name, (law, choice) in ROWS.items():
                e = strength/max(x, cv if choice == 'v2' else ca)
                D = nu+cd*k*k/(sigma*e) if law == 'scalar' else d+a*x
                flux[name] = D*flux['k_r']
            return flux

        checks = []
        for lo, hi in zip(profile.material[:-1], profile.material[1:]):
            sites = ((lo, 1, hi-lo), ((lo+hi)/2, 1, (hi-lo)/2), (hi, -1, hi-lo))
            for x, direction, room in sites:
                historical, rows = profile.evaluate(np.array([x]), direction)
                exact = {name: float(row['diffusion_divergence'][0]) for name, row in rows.items()}
                exact['k_r'] = float(historical['k_rr'][0])
                scales = dict.fromkeys(ROWS, ref.A/ref.H)
                scales['k_r'] = ref.Kin/ref.H**2
                result = []
                for div in (2048, 4096, 8192, 16384, 32768):
                    step = min(h/div, mp.mpf(room)/8)*direction
                    points = [values(mp.mpf(x)+j*step) for j in range(len(STENCIL))]
                    slopes = {name: float(mp.fsum(w*p[name] for w, p in zip(STENCIL, points))/(60*step))
                              for name in exact}
                    error = {name: abs(slopes[name]-v)/scales[name] for name, v in exact.items()}
                    result.append(dict(h=float(step), numerical=slopes, scaled_error=error))
                checks.append(dict(r=x, side=direction, analytic=exact, refinement=result))
        return checks


def integrations(profile, order):
    ref = profile.ref
    nodes, weights = leggauss(order)
    names = ('diffusion_divergence', 'P', 'E', 'source', 'rhs', 'identity')
    cumulative = {name: np.zeros(len(names)) for name in ROWS}
    totals = []
    for lo, hi in zip(profile.integration[:-1], profile.integration[1:]):
        r = lo+(hi-lo)/2*(nodes+1)
        _, rows = profile.evaluate(r)
        _, endpoints = profile.evaluate(np.array([hi]))
        regional = {}
        for name, row in rows.items():
            integral = np.array([weights@row[key] for key in names])*(hi-lo)/(2*ref.A)
            cumulative[name] += integral
            actual = float(endpoints[name]['Q'][0]/ref.A)
            required = cumulative[name][2]-cumulative[name][1]
            regional[name] = dict(integral_over_A=dict(zip(names, integral.tolist())),
                                  cumulative_over_A=dict(zip(names, cumulative[name].tolist())),
                                  actual_Q_over_A=actual, required_Q_over_A=required,
                                  cumulative_flux_defect_over_A=actual-required)
        totals.append(dict(r=[lo, hi], rows=regional))
    return totals


def signed_profiles(profile, integration_panels, order=64):
    g, w = leggauss(order)
    cumulative = {name: 0. for name in ROWS}
    output = []
    for panel in integration_panels:
        lo, hi = panel['r']
        outer = np.linspace(lo, hi, 33)
        historical, terms = profile.evaluate(outer)
        # Preserve the approached side at each material endpoint.
        _, left_limit = profile.evaluate(np.array([hi]), -1)
        r = lo+(outer[:, None]-lo)*(g+1)/2
        _, integral_terms = profile.evaluate(r.ravel())
        for name, row in terms.items():
            for key, val in row.items():
                val[-1] = left_limit[name][key][0]
            integrated = (integral_terms[name]['source'].reshape(len(outer), order)@w)*(outer-lo)/2
            row['cumulative_flux_defect'] = row['Q']+cumulative[name]+integrated
            cumulative[name] += integrated[-1]
        a, b, c, d = (terms[key]['rhs'] for key in 'ABCD')
        effects = dict(diffusion_at_v2_cutoff=b-a, direct_sink_with_auxiliary=d-b,
                       cutoff_at_scalar_diffusion=c-a, interaction=d-c-b+a,
                       collar_production=historical['collar'])
        arrays = lambda records: {k: v.tolist() for k, v in records.items()}
        output.append(dict(r=outer.tolist(), historical=arrays(historical),
                           rows={name: arrays(row) for name, row in terms.items()}, effects=arrays(effects)))
    return output


def assess(args, deadline):
    require(identity(args.table)['sha256'] == TABLE_SHA, 'v2 reference table differs')
    require(identity(args.localization)['sha256'] == LOCALIZATION_SHA, 'CP2 fixed-owner evidence differs')
    require(identity(args.earlier_columns)['sha256'] == OLD_COLUMNS_SHA, 'earlier factorial evidence differs')
    ref = Reference(read(args.table))
    profile = FixedProfile(ref)
    require(0 < min(profile.cutoffs.values()) <= max(profile.cutoffs.values()) < .4*ref.H,
            'fixed cutoffs must lie below the collar')
    loc = read(args.localization)
    fixed = loc['fixed_state_reference']
    columns = loc['independent_rules'][-1]['columns']
    area = math.fsum(c['area'] for c in columns)
    rms = math.sqrt(math.fsum(c['integral_square'] for c in columns)/area)
    P = math.fsum(c['production'] for c in columns)
    E = math.fsum(c['dependent_destruction'] for c in columns)
    cp2_audit = identity(loc['audit']['path'])
    require(cp2_audit == loc['audit'], 'retained CP2 independent audit hash differs')
    audit = read(cp2_audit['path'])
    require(audit['core']['physical_envelope']['domain/inner/k'] == fixed['inner_k'], 'CP2 envelope record differs')
    c1, c2 = integrations(profile, 32), integrations(profile, 64)
    samples = signed_profiles(profile, c2)
    independent = derivative_check(profile, deadline)
    scale = ref.A/ref.H
    cp2 = loc['parallel_inlet']
    hground, ground = profile.evaluate(np.array([0.]))
    oldP, oldE, curP, curE = (c2[-1]['rows'][name]['cumulative_over_A'][quantity]
                              for name, quantity in (('D', 'P'), ('D', 'E'), ('A', 'P'), ('A', 'E')))
    historicalP = float(ref.old_integral_p(ref.H))
    archived = cp2['independent_rules'][-1]['integrated_over_A']
    qH = c2[-1]['rows']['A']
    errors = dict(cp2_RMS=abs(rms-fixed['inner_k']['RMS']),
        cp2_inner_production=abs(P-fixed['inner_production'])/(25*ref.A),
        cp2_inner_destruction=abs(E-fixed['inner_destruction'])/(25*ref.A),
        cp2_historical_cutoff=abs(ref.old_c-cp2['old_cutoff'])/ref.H,
        cp2_v2_cutoff=abs(profile.cutoffs['v2']-cp2['v2_cutoff_at_Kin'])/ref.H,
        cp2_ground_diffusion_ratio=abs(float(ground['A']['diffusion_ratio'][0])-cp2['ground_v2_over_historical_diffusion']),
        cp2_parallel_budgets=max(abs(a-b) for a, b in zip((historicalP, curP, oldE, curE), archived)),
        cp2_current_H_flux=abs(qH['actual_Q_over_A']-cp2['independent_rules'][-1]['cuts'][-1]['actual_v2_gradient_flux_over_A']),
        historical_ground_and_H=max(abs(float(v)) for v in ref.old_Q(np.array([0., ref.H])))/ref.A,
        historical_budget=abs(oldE-historicalP),
        cutoff_invariant_production=abs(oldP-curP),
        ratio_identity=max(np.max(np.abs(np.asarray(row['rhs'])-row['identity']))/scale
                           for panel in samples for row in panel['rows'].values()),
        row_D_lower=max(np.max(np.abs(panel['rows']['D']['rhs']))/scale
                        for panel in samples if panel['r'][-1] <= .4*ref.H),
        row_D_collar=max(np.max(np.abs(np.asarray(panel['rows']['D']['rhs'])-panel['effects']['collar_production']))/scale
                         for panel in samples),
        cumulative_fundamental_theorem=max(abs(row['cumulative_over_A']['rhs']-row['cumulative_flux_defect_over_A'])
                                           for panel in c2 for row in panel['rows'].values()),
        split_quadrature_refinement=max(abs(row['cumulative_over_A'][k]-b['rows'][name]['cumulative_over_A'][k])
                                        for a,b in zip(c1,c2) for name,row in a['rows'].items() for k in row['cumulative_over_A']),
        derivative_refined=max(max(row['refinement'][-1]['scaled_error'].values()) for row in independent))
    ground_out = {name: {key: float(v[0]) for key,v in row.items()} for name,row in ground.items()}
    status = 'PASS' if max(errors.values()) <= TOL and ground_out['A']['rhs'] < 0 < ground_out['A']['source'] else 'FAIL'
    budget = {name: dict(lower=c2[next(i for i,p in enumerate(c2) if p['r'][1] == .4*ref.H)]['rows'][name],
                         whole=c2[-1]['rows'][name]) for name in ROWS}
    summary = dict(schema='WJInnerKFixedProfileMechanismV1', MECHANISM_DIAGNOSIS=status, TASK3_READY=False,
        source19=loc['source_identity'], table=identity(args.table), localization=identity(args.localization),
        cp2_audit=cp2_audit, earlier_factorial=identity(args.earlier_columns), script=identity(__file__),
        rows=ROWS, cutoffs=profile.cutoffs, ground_historical={k: float(v[0]) for k,v in hground.items()},
        ground=ground_out, scales=dict(k=ref.Kin, action=scale, flux=ref.A), check_limit=TOL, checks=errors,
        integrated=budget, historical_integrated_P_over_A=historicalP,
        cp2_2D=dict(inner_k=fixed['inner_k'], reintegrated_RMS=rms, sampled_endpoint_peak=loc['sampled_peak'],
                    first_two_squared_share=sum(x['fraction_global_squared_departure'] for x in columns[:2]),
                    downstream_RMS=[c['RMS_over_Kin'] for c in columns[6:]],
                    P_inner=P, E_inner=E, weak_k_budget_error=fixed['k']['assembly_identity_error_over_scale']))
    detail = dict(material_breaks=profile.material, integrations={32: c1, 64: c2},
                  numerical_derivatives=independent, signed_profiles=samples)
    return summary, detail


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('table', 'localization', 'earlier-columns', 'out'):
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--seconds', type=float, required=True)
    args = p.parse_args()
    require(sys.platform.startswith('linux') and 0 < args.seconds <= 100, 'bounded Linux analysis required')
    args.out.mkdir()
    start = time.monotonic()
    summary, detail = assess(args, start+args.seconds)
    summary['elapsed_seconds'] = time.monotonic()-start
    for name, result in (('mechanism.json', summary), ('signed-profiles.json', detail)):
        with (args.out/name).open('x') as stream:
            json.dump(result, stream, indent=2, allow_nan=False)
            stream.write('\n')
    print(json.dumps(dict(status=summary['MECHANISM_DIAGNOSIS'], checks=summary['checks'],
                          elapsed=summary['elapsed_seconds'])), flush=True)
    if summary['MECHANISM_DIAGNOSIS'] != 'PASS':
        raise SystemExit(2)


if __name__ == '__main__':
    main()
