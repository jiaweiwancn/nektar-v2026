#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Additional old-runtime equivalence and forced late-IVP rollback checks."""
import argparse
from pathlib import Path
import xml.etree.ElementTree as ET

from continue_cls_epawf_pimple import CampaignLock, atomic
from run_cls_awf_late_ivp import native, packet, prepare
from run_cls_epawf_picard import load, require, sha, state_hash, set_parameter, write
from diagnose_cls_accepted_covariance import checkpoint


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('solver', 'reference-solver', 'study', 'output'): p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--timeout', type=float, default=1800); a = p.parse_args(); a.output.mkdir(parents=True, exist_ok=False)
    lock = CampaignLock(a.output); parent = a.study/'parent'; q = checkpoint(parent)
    expected = load(a.study/'import-validation/ordinary/case.evolution.json')
    imported = load(a.study/'import-validation/imported/case.evolution.json')
    old, receipt = native(a.reference_solver, prepare(a.output/'old-runtime', 1, q['dt'], restart=parent), a.timeout, lock)
    checks = dict(old_runtime_state_exact=state_hash(old['state']) == state_hash(expected['state']),
        old_runtime_histories_exact=[state_hash(x) for x in old['native_history']] == [state_hash(x) for x in expected['native_history']],
        old_runtime_times_exact=state_hash(old['native_history_times']) == state_hash(expected['native_history_times']))
    left, right = packet(old), packet(imported)
    checks['external_coefficients_exact'] = left['fixed_external_coefficients'] == right['fixed_external_coefficients']
    xml = prepare(a.output/'forced-rejection', 1, q['dt'], parent=parent)
    root = ET.parse(xml).getroot(); set_parameter(root, 'RANSWallPimpleHardAuditRejectStep', 1); write(root, xml)
    rejected, failure = native(a.solver, xml, a.timeout, lock, expected_success=False)
    initial = load(a.study/'coarse/initial/case.evolution.json')
    checks.update(rejection_state_exact=state_hash(rejected['state']) == state_hash(q['coupled']),
        rejection_time_exact=rejected['time'] == q['time'] and rejected['steps'] == 0,
        rejection_history_exact=[state_hash(x) for x in rejected['native_history']] == [state_hash(x) for x in initial['native_history']],
        rejection_time_labels_exact=state_hash(rejected['native_history_times']) == state_hash(initial['native_history_times']),
        actual_rejected_candidate_retained=bool(rejected['rejected_state']) and state_hash(rejected['rejected_state']) != state_hash(q['coupled']))
    attempt = rejected['pimple_attempts'][0]
    checks['rejected_attempt_one_sweep'] = all(attempt[k] == 1 for k in
        ('outer_sweeps', 'momentum_predictor_solves', 'pressure_system_solves', 'epsilon_transport_solves', 'k_transport_solves', 'bound_calls'))
    require(all(checks.values()), 'late IVP rollback/equivalence contract failed: '+str(checks))
    atomic(a.output/'summary.json', dict(status='passed', checks=checks, old_runtime_receipt=receipt,
        forced_rejection_receipt=failure, expected_rejection=rejected['failure'],
        parent_checkpoint_sha256=sha(parent/'State.bin'), old_solver_sha256=sha(a.reference_solver), solver_sha256=sha(a.solver)))


if __name__ == '__main__': main()
