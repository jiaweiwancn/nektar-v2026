#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only confirmation-pulse invariance and predeclared GL28 checks."""
import argparse, time
from pathlib import Path
import numpy as np
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import require, sha, state_hash
from wp2_all_state import Rows
from wp2_common_quadrature import resample


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    started = time.monotonic()
    progress = load(a.root / 'recovery-seed-confirmation-standard/progress.json')
    report_path = Path(progress['observations'][0]['report'])
    d = load(report_path)
    rows = Rows(a.root / 'geometry-rows-confirmation/result.json')
    state, cold = (np.asarray(d[key], dtype='<f8') for key in ('state', 'fixed_reference'))
    require(len(state) == 11014 and d['model_specification']['name'] == 'StandardKEpsilon', 'wrong seed scope')
    protected = np.ones(len(state), dtype=bool)
    protected[rows.offset[2]:rows.offset[3]] = False
    inner = np.unique(rows.a['v_ids'][rows.inner].reshape(-1)).astype(int) + rows.offset[2]
    protected[inner] = True
    essential = np.asarray(d['wp2_per_state_monitoring'][0]['physical_checks']['physical_essential_rows'], dtype=int)
    protected[essential] = True
    require(np.array_equal(state.view('<u8')[protected], cold.view('<u8')[protected]), 'pulse changed protected coefficient bits')
    actual, reference = rows.tuple(state, d), rows.tuple(cold, d)
    require(actual['endpoint_K'].tobytes() == reference['endpoint_K'].tobytes(), 'pulse changed endpoint K')
    require(np.array_equal(actual['wall']['nt'], reference['wall']['nt']) and
            np.array_equal(actual['wall']['E'], reference['wall']['E']), 'pulse changed dependent AWF data')
    dense_audit = rows.volume_audit(d, actual, allow_accuracy_failure=True)
    def dense(values):
        return resample(dict(d, wp2_common_samples=dict(d['wp2_common_samples'], values=values['values']),
            physical_samples=dict(wall_K=values['endpoint_K'])), 28)
    left, right = dense(actual), dense(reference)
    supply = float(np.dot(left['weights'], left['values'][:, 2] - right['values'][:, 2]))
    atomic(a.output, dict(version='WP2ConfirmationSeedIndependentV1', passed=True,
        seed_sha256=state_hash(state), fixed_reference_sha256=state_hash(cold),
        protected_rows=int(sum(protected)), protected_rows_bitwise_unchanged=True,
        endpoint_K_bitwise_unchanged=True, dependent_epsilon_and_auxiliary_viscosity_unchanged=True,
        setup_integrated_TKE_increment_per_density_unit_span=supply, GL28=dense_audit,
        report=dict(path=str(report_path), sha256=sha(report_path)), rows_sha256=sha(rows.path),
        source_sha256=sha(Path(__file__)), seconds=time.monotonic()-started, governing_solves=0,
        scope='Retrospective read-only check of the exact saved h0. Predeclared seen GL28 volume grid '
              'complements native sided inventory. Not a fresh audit or an alteration of the running state.'))


if __name__ == '__main__':
    main()
