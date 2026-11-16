#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only native-element replay of a retained WP2 recovery rejection."""
import argparse
import time
import resource
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss, legvander, legder
from cao_tamura_h2d_state import load
from run_cls_epawf_picard import require, sha, state_hash
from continue_cls_epawf_pimple import atomic
from wp2_all_state import Rows, wall


def main():
    parser = argparse.ArgumentParser()
    for name in ['root', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    root = args.root
    failed_path = root / 'preservation-recovery-standard/step390/case.evolution.json'
    parent_path = root / 'preservation-recovery-standard/step368/case.evolution.json'
    failed, parent = load(failed_path), load(parent_path)
    require(len(failed['physical_range_events']) == 1, 'expected first rejected attempt only')
    event = failed['physical_range_events'][0]
    audit = event['audit']
    require(audit['accepted_awf_hard_admissible'] == [0], 'not a hard-audit rejection')
    witness = audit['accepted_awf_inner_minimum_point']
    geometry, element, radial_side, block = [int(i) for i in audit['native_inner_minimum_geometry_element_side_block']]
    state = np.asarray(event['state'])
    rows = Rows(root / 'geometry-rows-selected/result.json')
    values = rows.tuple(state, failed)
    partition = np.asarray(rows.meta['physical_z_partition'])
    panel = int(np.searchsorted(partition, witness[2], side='right') - 1)
    selected = ((np.asarray(rows.a['geometry'][:, 0]) == geometry)
                & (rows.physical[:, 2] > partition[panel])
                & (rows.physical[:, 2] < partition[panel + 1]))
    coords = rows.physical[selected]
    xx, zz = np.unique(coords[:, 0]), np.unique(coords[:, 2])
    require(len(xx) == len(zz) == 8 and sum(selected) == 64, 'native element quadrature identity')
    native_element = failed['fields'][0]['mesh'][element]
    require(native_element['geometry'] == geometry, 'element/geometry provenance differs')
    lo, hi = np.asarray(native_element['lower']), np.asarray(native_element['upper'])
    # The immutable common catalogue can split one native element at several
    # source/mesh breakpoints. Fit on exactly one such GL8 panel, retaining the
    # native element and its one-sided tangential limit.
    fit_lo, fit_hi = lo.copy(), hi.copy()
    fit_lo[2], fit_hi[2] = partition[panel:panel + 2]
    reference_point = 2 * (np.asarray(witness[:3]) - fit_lo) / (fit_hi - fit_lo + (fit_hi == fit_lo)) - 1
    nodes, weights = leggauss(8)
    V = legvander(nodes, 3)
    project = (2 * np.arange(4) + 1)[:, None] / 2 * V.T * weights
    vx = legvander([reference_point[0]], 3)[0]
    vz = legvander([reference_point[2]], 3)[0]
    def polynomial(data):
        data = np.asarray(data)[selected].reshape(8, 8, -1)
        coeff = np.einsum('ai,ijf,bj->abf', project, data, project)
        return np.einsum('a,abf,b->f', vx, coeff, vz), coeff
    primitive, coefficients = polynomial(values['values'])
    gradient, _ = polynomial(values['G'].reshape(-1, 9))
    G = gradient.reshape(3, 3)
    dx = legvander([reference_point[0]], 2)[0] @ legder(np.eye(4), axis=0) * 2 / (fit_hi[0] - fit_lo[0])
    dz = legvander([reference_point[2]], 2)[0] @ legder(np.eye(4), axis=0) * 2 / (fit_hi[2] - fit_lo[2])
    independent_G = np.zeros((3, 3))
    for component, column in [(0, 0), (2, 1)]:
        independent_G[component, 0] = np.einsum('a,ab,b', dx, coefficients[:, :, column], vz)
        independent_G[component, 2] = np.einsum('a,ab,b', vx, coefficients[:, :, column], dz)
    gradient_error = float(np.max(abs(independent_G - G)))
    require(gradient_error < 2e-8, 'independent represented-polynomial gradient replay differs')
    # Endpoint K uses its own native functional. Reconstruct its continuous
    # tangential polynomial, never use k at the inner witness as endpoint K.
    endpoint_grid = values['wall']['nt'] * 0 + np.repeat(values['endpoint_K'], len(values['values']) // len(values['endpoint_K']))
    endpoint, _ = polynomial(endpoint_grid)
    K = float(endpoint[0])
    nu, H, hs = failed['parameters'][:3]
    q, kappa = failed['parameters'][6:8]
    profile = wall(np.asarray([K]), np.asarray([witness[2]]), nu, H, hs, q, kappa)
    nt = float(profile['nt'][0])
    k = float(primitive[2])
    isotropic = np.eye(3) * (2 * k / 3)
    shear = -nt * (G + G.T)
    tensor = isotropic + shear
    eigenvalues, eigenvectors = np.linalg.eigh(tensor)
    require(max(abs(k - witness[3]), abs(K - witness[4]), abs(nt - witness[5])) < 2e-12, 'native scalar/profile witness differs')
    require(abs(eigenvalues[0] - witness[-1]) < 2e-11 and eigenvalues[0] < 0, 'negative native tensor not reproduced')
    checks = {}
    for name in ['state', 'native_history', 'native_history_times', 'fixed_reference', 'fixed_defect']:
        checks[name + '_restored_bitwise'] = np.asarray(failed[name], dtype='<f8').tobytes() == np.asarray(parent[name], dtype='<f8').tobytes()
    checks['accepted_time_preserved'] = failed['time'] == parent['time']
    checks['candidate_retained'] = state_hash(event['state']) == state_hash(failed['rejected_state'])
    checks['scalar_domain_admissible'] = audit['accepted_awf_scalar_domain'][0] == 1
    checks['budget_admissible'] = audit['accepted_awf_scalar_domain'][7] >= 1 and audit['accepted_awf_scalar_domain'][8] <= 1e-10
    attempt = failed['pimple_attempts'][0]
    checks['one_each_governing_solve'] = all(attempt[n] == 1 for n in ['momentum_predictor_solves', 'pressure_system_solves', 'epsilon_transport_solves', 'k_transport_solves', 'bound_calls'])
    checks['no_bound_activation'] = attempt['bound_activations'] == 0
    checks['solved_epsilon_retained'] = attempt['solved_epsilon_retained_bitwise']
    require(all(checks.values()), 'rejection transaction or numerical-path contract differs')
    checkpoint = Path(parent['checkpoint']) / 'State.bin'
    rollback = {key: value for key, value in failed['rollback_receipt'].items() if not isinstance(value, (list, dict))}
    result = dict(version='WP2RejectedRecoveryReplayV1', model='calibrated_standard',
        actual_accepted_time=parent['time'], first_rejected_attempt_time=event['time'], accepted_steps=368, attempted_steps=369,
        failure_stage='final complete accepted-AWF audit after all one-sweep solves and logical bound, before coupled candidate assignment and history commit',
        failure=failed['failure'], parent_state_sha256=state_hash(parent['state']), rejected_state_sha256=state_hash(event['state']),
        attempted_history_sha256=state_hash(event['history']), parent_checkpoint=dict(path=str(checkpoint), sha256=sha(checkpoint)),
        reference_sha256=state_hash(parent['fixed_reference']), defect_sha256=state_hash(parent['fixed_defect']),
        native_history_sha256=[state_hash(s) for s in parent['native_history']], native_history_times=parent['native_history_times'],
        witness=dict(physical=witness[:3], geometry=geometry, element=element, radial_side=radial_side,
            tangential_side=-1 if abs(witness[0] - hi[0]) < 1e-14 else 1 if abs(witness[0] - lo[0]) < 1e-14 else 0,
            native_element_lower=lo.tolist(), native_element_upper=hi.tolist(), block=block, k=k, endpoint_K=K,
            dependent_epsilon=float(profile['E'][0]), eddy_viscosity=nt, gradient=G.tolist(),
            independent_gradient=independent_G.tolist(), gradient_replay_maximum_error=gradient_error,
            isotropic_tensor=isotropic.tolist(), linear_stress_contribution=shear.tolist(), reynolds_tensor=tensor.tolist(),
            eigenvalues=eigenvalues.tolist(), critical_eigenvector=eigenvectors[:, 0].tolist(), native_minimum=witness[-1],
            normalized_minimum_over_K0=float(eigenvalues[0] / failed['parameters'][5])),
        scalar_and_budget_domain=audit['accepted_awf_scalar_domain'], covariance_inner_outer=audit['accepted_awf_covariance'],
        checks=checks, rollback_receipt=rollback,
        bindings={str(p): sha(p) for p in [failed_path, parent_path, rows.path, Path(__file__)]},
        cost=dict(seconds=time.monotonic()-started, peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
        zero_PDE_solves=True, interpretation='genuine inner Linear covariance rejection in this discrete dt=.01 recovery; scalar and cutoff budgets pass. This does not establish global model infeasibility or timestep-independent instability.')
    atomic(args.output, result)


if __name__ == '__main__':
    main()
