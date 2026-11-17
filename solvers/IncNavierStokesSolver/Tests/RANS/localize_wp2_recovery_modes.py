#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Zero-solve physical/modal localization of retained flat recovery fields.

Local pressure-polynomial visibility is an L2 orthogonality diagnostic. It is
not a global inf-sup test or a prediction that another pressure space will work.
"""
import argparse, time, resource
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import legvander
from analyze_wp2_approach import native_states, paths
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import require, sha, state_hash
from wp2_all_state import Rows, US


def modal_energy(x, z, weight, values, area):
    values = np.asarray(values)
    vx, vz = legvander(x, 3), legvander(z, 3)
    factor = (2*np.arange(4)+1)[:, None]*(2*np.arange(4)+1)[None, :]
    coefficients = np.einsum('i,ia,ib,if->abf', weight, vx, vz, values)*factor[:, :, None]/area
    reconstructed = np.einsum('ia,abf,ib->if', vx, coefficients, vz)
    require(np.max(abs(values-reconstructed)) < 2e-10*max(1., np.max(abs(values))),
            'native field is not the declared affine cubic polynomial')
    energy = coefficients**2*area/factor[:, :, None]
    integrated = np.einsum('i,if->f', weight, values**2)
    require(np.all(abs(np.sum(energy, axis=(0, 1))-integrated) <= 1e-9*integrated+area*1e-24),
            'element modal energies do not reproduce physical quadrature')
    return energy, integrated


def localize(rows, context, state, cold):
    actual = rows.tuple(state, context)
    divergence = actual['G'][:, 0, 0]+actual['G'][:, 2, 2]
    cold_divergence = cold['G'][:, 0, 0]+cold['G'][:, 2, 2]
    fields = np.column_stack((divergence, actual['values'][:, :3]-cold['values'][:, :3],
                              actual['values'][:, 4]-cold['values'][:, 4]))
    elements = rows.a['elements'][:, 0].astype(int)
    table, aggregate = [], np.zeros((4, 4, 5))
    for element in np.unique(elements):
        mask = elements == element
        definition = context['fields'][0]['mesh'][int(element)]
        lo, hi = np.asarray(definition['lower']), np.asarray(definition['upper'])
        area = (hi[0]-lo[0])*(hi[2]-lo[2])
        require(abs(sum(rows.weight[mask])-area) <= 1e-12*area, 'element physical measure mismatch')
        point = rows.physical[mask]
        x, z = 2*(point[:, 0]-lo[0])/(hi[0]-lo[0])-1, 2*(point[:, 2]-lo[2])/(hi[2]-lo[2])-1
        energy, integrated = modal_energy(x, z, rows.weight[mask], fields[mask], area)
        aggregate += energy
        index = int(np.argmax(abs(divergence[mask])))
        table.append(dict(element=int(element), geometry=definition['geometry'], lower=lo.tolist(), upper=hi.tolist(),
            area=area, divergence_L2_squared=float(integrated[0]),
            maximum_absolute_divergence=float(max(abs(divergence[mask]))),
            limiting_physical=point[index].tolist(), k_departure_L2_squared=float(integrated[3]),
            pressure_departure_L2_squared=float(integrated[4])))
    total = float(np.sum(aggregate[:, :, 0]))
    require(abs(total-np.dot(rows.weight, divergence**2)) < 1e-9*total+1e-24,
            'localization changed the strong-divergence integral')
    fractions = {}
    for degree in (1, 2):
        high = np.logical_or(np.arange(4)[:, None] > degree, np.arange(4)[None, :] > degree)
        fractions['outside_local_Q'+str(degree)] = None if total <= 1e-24 else float(sum(aggregate[:, :, 0][high])/total)
    production_departure = actual['values'][:, 9]-cold['values'][:, 9]
    isotropic_departure = -2/3*(actual['values'][:, 2]*divergence-cold['values'][:, 2]*cold_divergence)
    work = []
    for region in (0, 1):
        mask = rows.inner if region == 0 else ~rows.inner
        weight = rows.weight[mask]
        denominator = float(np.dot(weight, production_departure[mask]**2))
        work.append(dict(region=region,
            total_production_departure_RMS=float(np.sqrt(denominator/sum(weight))),
            isotropic_work_departure_RMS=float(np.sqrt(np.dot(weight, isotropic_departure[mask]**2)/sum(weight))),
            signed_projection_on_total_production_departure=None if denominator <= 1e-24 else
                float(np.dot(weight, production_departure[mask]*isotropic_departure[mask])/denominator)))
    table.sort(key=lambda entry: -entry['divergence_L2_squared'])
    for entry in table:
        entry['fraction_of_physical_divergence_L2_squared'] = None if total <= 1e-24 else entry['divergence_L2_squared']/total
    maximum = int(np.argmax(abs(divergence)))
    return dict(strong_divergence_L2_squared=total,
        strong_divergence_RMS_over_ustar_H=float(np.sqrt(total/sum(rows.weight))/(US/.0065)),
        maximum_absolute_divergence=float(abs(divergence[maximum])),
        maximum_divergence_physical=rows.physical[maximum].tolist(),
        pressure_polynomial_orthogonal_energy_fractions=fractions,
        field_order=['divergence', 'U_minus_native_cold', 'W_minus_native_cold', 'k_minus_native_cold', 'pi_minus_native_cold'],
        global_modal_physical_L2_squared=aggregate.tolist(), largest_divergence_elements=table[:10],
        production_work=work)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    rows = Rows(args.root/'geometry-rows-selected/result.json')
    require(rows.a['v_ids'].shape[1] == 16 and rows.a['p_ids'].shape[1] == 4,
            'modal diagnostic requires the retained Q3 velocity / Q1 pressure fields')
    cases = [('CLS_coarse', 'CLS', args.root/'preservation-recovery-CLS/progress.json'),
             ('standard_coarse', 'standard', args.root/'preservation-recovery-standard/progress.json'),
             ('standard_half', 'standard', args.root/'global-half-step-terminal-input/fine.json')]
    output, bindings = {}, {}
    for label, model, progress_path in cases:
        progress = load(progress_path)
        bindings[str(progress_path)] = sha(progress_path)
        states, _ = native_states(args.root, model, 'recovery', progress_path)
        context_path = paths(args.root, model, 'recovery', progress)[0]
        context = load(context_path);bindings[str(context_path)] = sha(context_path)
        require(context['inner_stress_model'] == 'Linear', 'modal localization changed the inner stress law')
        cold = rows.tuple(context['fixed_reference'], context)
        output[label] = dict(dt=progress['dt'], reference_sha256=state_hash(context['fixed_reference']),
            solver_sha256=progress['solver_sha256'], native_tuple_replay=rows.check_native(context), history=[])
        for time_key in (.3, .6, .9, 1.2):
            actual_time, state = states[time_key]
            observation = next(o for o in progress['observations'] if abs(o['time']-actual_time) < 2e-13)
            report_path = Path(observation['report'])
            bindings[str(report_path)] = sha(report_path)
            result = localize(rows, context, state, cold)
            gates = observation['final_physical_checks']['gates']
            native = next(v['value'] for k, v in gates.items() if k.startswith('strong_divergence'))
            require(abs(native-result['strong_divergence_RMS_over_ustar_H']) < 2e-10,
                    'physical quadrature does not reproduce native strong divergence')
            output[label]['history'].append(dict(actual_time=actual_time, state_sha256=state_hash(state),
                native_weak_continuity_gates={k:v for k,v in gates.items() if 'continuity' in k}, **result))
    for path, expected in bindings.items():
        require(sha(Path(path)) == expected, 'localization input changed')
    atomic(args.output, dict(version='WP2RecoveryModalLocalizationV1', cases=output, bindings=bindings,
        source_sha256=sha(Path(__file__)),
        evaluation_sources={name:sha(Path(__file__).with_name(name)) for name in
            ['wp2_all_state.py','cao_tamura_h2g_constitutive.py','analyze_wp2_approach.py']},
        rows_sha256=sha(rows.path), seconds=time.monotonic()-started,
        peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, governing_solves=0,
        scope='Affine native-element Legendre decomposition with physical GL8 integration, exact for squared cubic fields. '
              'Outside local Q1 is orthogonal to every local Q1 pressure polynomial; continuous global pressure has '
              'additional conformity restrictions. This is not a global pressure-space feasibility or inf-sup proof, '
              'not a curved-metric claim, not a residual percentage, and not a runtime correction. '
              'Production terms retain their signed projection and are not added as norms.'))


if __name__ == '__main__':
    main()
