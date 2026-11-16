#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Signed frozen-operator attribution of the saved rejected k increment.

The adjoint is diagnostic linear algebra, not another governing solve. Its
individual signed load responses retain cancellation and do not constitute
separate unconstrained physical evolution equations.
"""
import argparse, resource, time
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss, legvander
from scipy.linalg import solve
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from diagnose_cls_accepted_covariance import common_packet, run_probe, STAGES
from run_cls_epawf_picard import require, sha, state_hash
from wp2_all_state import Rows


def main():
    parser = argparse.ArgumentParser()
    for key in ('root', 'probe', 'output'):
        parser.add_argument('--' + key, type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    args.output.mkdir(exist_ok=False)
    report_path = args.root / 'preservation-recovery-standard/step390/case.evolution.json'
    report = load(report_path)
    packet_path = Path(report['attempt_packet_status']['last_written_packet'])
    require(sha(packet_path) == report['attempt_packet_status']['sha256'], 'archived packet changed')
    archived = load(packet_path)
    require(sha(Path(archived['matrix_file'])) == archived['matrix_sha256'], 'archived matrix changed')
    checkpoint = Path(load(args.root / 'preservation-recovery-standard/step368/case.evolution.json')['checkpoint'])
    parent_hash = sha(checkpoint / 'State.bin')
    packet = dict(common_packet(checkpoint), operation='archived_k_balance',
        archived_packet=str(packet_path), archived_packet_sha256=sha(packet_path),
        stage_hashes={name: state_hash(archived[name]) for name in STAGES},
        expected_history_sha256=state_hash(archived['history']),
        state=archived['committed'], committed=archived['committed'], history=archived['history'],
        fixed_reference=archived['fixed_reference'], fixed_external_coefficients=archived['fixed_external_coefficients'],
        alpha=archived['alpha'], time=archived['time'], inner_stress_model=report['inner_stress_model'],
        expected_model='StandardKEpsilon')
    xml = report_path.with_name('case.xml')
    native, receipt = run_probe(args.probe, xml, packet, args.output / 'native')
    require(native['governing_solves'] == native['physical_commits'] == 0 and
            native['initialization_work'] == [0, 0], 'diagnostic advanced or initialized a state')
    require(native['pre_k_sha256'] == state_hash(archived['pre_k']), 'diagnostic changed the pre-k state')
    free = np.asarray(archived['transport_rows'], dtype=int)
    require(np.array_equal(free, np.asarray(native['transport_rows'], dtype=int)), 'transport row ordering changed')
    n = len(free)
    matrix = np.fromfile(archived['matrix_file'], dtype='<f8').reshape(n, n)
    rows = Rows(args.root / 'geometry-rows-selected/result.json')
    audit = report['physical_range_events'][0]['audit']
    point = np.asarray(audit['accepted_awf_inner_minimum_point'][:3])
    geometry = int(audit['native_inner_minimum_geometry_element_side_block'][0])
    bounds = next(e for e in report['fields'][0]['mesh'] if e['geometry'] == geometry)
    partition = np.asarray(rows.meta['physical_z_partition'])
    index = int(np.searchsorted(partition, point[2], side='right') - 1)
    mask = ((rows.a['geometry'][:, 0] == geometry) &
            (rows.physical[:, 2] > partition[index]) & (rows.physical[:, 2] < partition[index + 1]))
    require(sum(mask) == 64, 'native sided GL8 panel differs')
    nodes, weight = leggauss(8)
    V = legvander(nodes, 3)
    projector = (2 * np.arange(4) + 1)[:, None] / 2 * V.T * weight
    x = 2 * (point[0] - bounds['lower'][0]) / (bounds['upper'][0] - bounds['lower'][0]) - 1
    z = 2 * (point[2] - partition[index]) / (partition[index + 1] - partition[index]) - 1
    ix, iz = legvander([x], 3)[0] @ projector, legvander([z], 3)[0] @ projector
    phi = np.einsum('i,ijc,j->c', ix, rows.a['v_phi'][mask].reshape(8, 8, -1), iz)
    ids = rows.a['v_ids'][mask]
    require(np.all(ids == ids[0]), 'panel changed native coefficient identities')
    functional = np.zeros(n)
    lookup = {row: i for i, row in enumerate(free)}
    for row, value in zip(ids[0] + rows.offset[2], phi):
        require(int(row) in lookup, 'limiting interior k functional contains a prescribed row')
        functional[lookup[int(row)]] += value
    adjoint = solve(matrix.T, functional, assume_a='gen', check_finite=True)
    residual = np.asarray(archived['transport_residual'])
    increment = np.asarray(archived['raw_increment'])
    changes = {name: -float(np.dot(adjoint.astype(np.longdouble), np.asarray(load, dtype=np.longdouble)))
               for name, load in native['components'].items()}
    from_loads = float(sum(np.longdouble(v) for v in changes.values()))
    frozen_change = float(functional @ increment)
    def value(state):
        data = rows.tuple(np.asarray(state), report)['values'][mask, 2].reshape(8, 8)
        return float(np.einsum('i,ij,j->', ix, data, iz))
    pre, raw = value(archived['pre_k']), value(archived['raw'])
    backward = float(np.max(abs(matrix.T @ adjoint - functional)) /
        (np.linalg.norm(matrix.T, ord=np.inf) * np.linalg.norm(adjoint, ord=np.inf) + np.linalg.norm(functional, ord=np.inf)))
    defect = float(max(abs(from_loads - frozen_change), abs(raw - pre - frozen_change)))
    require(backward < 2e-13 and defect < 2e-11, 'signed witness decomposition does not reproduce the actual increment')
    require(state_hash(archived['raw']) == state_hash(archived['corrected']), 'bound changed the rejected candidate')
    require(sha(checkpoint / 'State.bin') == parent_hash, 'readback changed the accepted parent')
    atomic(args.output / 'result.json', dict(version='WP2FrozenKWitnessBalanceV1',
        physical=point.tolist(), geometry=geometry, radial_side=0, tangential_side=1,
        pre_k_value=pre, raw_k_value=raw, actual_k_change=raw-pre,
        signed_frozen_load_responses=changes, reconstructed_change=from_loads,
        independent_native_functional_change=frozen_change, reconstruction_absolute_defect=defect,
        adjoint_backward_residual=backward,
        observed_load_rounding_effect=float(abs(adjoint @ (np.asarray(native['native']) - residual))),
        conservative_load_rounding_bound=float(sum(abs(adjoint)) * native['packet_residual_maximum_defect']),
        native_residual_export=receipt, rows_sha256=sha(rows.path), matrix_sha256=archived['matrix_sha256'],
        source_report=dict(path=str(report_path), sha256=sha(report_path)),
        archived_packet=dict(path=str(packet_path), sha256=sha(packet_path)),
        parent_checkpoint=dict(path=str(checkpoint), sha256=parent_hash),
        validation_probe=dict(path=str(args.probe), sha256=sha(args.probe)),
        source_sha256=sha(Path(__file__)), zero_governing_solves=True, zero_lift_projection=True,
        correction_activated=False, seconds=time.monotonic()-started,
        peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        interpretation='Signed responses through the exact saved frozen k matrix and native free-row loads. '
                       'Includes boundary terms; constrained rows are not treated as evolution equations. '
                       'Not additive error percentages, a nonlinear causality decomposition, or a separate physical trajectory.'))


if __name__ == '__main__':
    main()
