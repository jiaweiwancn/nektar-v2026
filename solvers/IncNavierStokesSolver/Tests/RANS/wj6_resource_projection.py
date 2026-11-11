#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Freeze measured WJ-6 costs, including parent algebra/verification and I/O.

This passive observer leaves the frozen native/measurement driver untouched.
Only a not-yet-created phase resource projection and progress pointer are set.
"""
import argparse
import json
from pathlib import Path
import time

from wj_immutable import (Ledger, atomic, artifacts, file_identity, read,
                          verify_artifacts, write_once)
from wj6_resources import DEFAULT, volume
from wj_qualification import component
from run_cls_epawf_picard import require


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--study', type=Path, required=True)
    parser.add_argument('--policy', choices=('EP', 'Budget'), required=True)
    parser.add_argument('--phase', choices=('Tx', 'T3'), required=True)
    a = parser.parse_args()
    started = time.time()
    root, folder = a.root.resolve(), a.study.resolve()/a.policy
    s = json.loads((folder/'summary.json').read_text())
    ledger = Ledger(root, a.policy, cap=800)
    archive = root/'resource-projections-v1'/a.policy/a.phase
    archive.mkdir(parents=True, exist_ok=False)
    label = 'resource-projection/'+a.phase
    ledger.reserve(label, dict(code=artifacts([Path(__file__)])), 0, 60)
    try:
        design = read(folder/'design.json')
        target = design['primary' if a.phase == 'Tx' else 'extension']
        headroom = 1.25
        require(set(s['trajectories']) == {'coarse', 'fine', 'baseline'}, 'incomplete measured pilot')
        executions = [r for r in s['execution'].values() if r.get('actual_attempts', 0) and r['exit_code'] == 0]
        rates = {tag: max(r['resource']['seconds']/r['actual_attempts'] for r in executions if r['role'] == tag)*headroom
                 for tag in ('coarse', 'fine', 'baseline')}
        obs, svd = {}, {}
        for tag, t in s['trajectories'].items():
            obs[tag], svd[tag] = [], []
            for c in [t['initial']]+t['trajectory']:
                verify_artifacts({c['qualified_record']['path']: c['qualified_record']})
                pair = [Path(c['path']).parent/'receipt.json',
                        Path(c['quadrature']['independent']['path']).parent/'receipt.json']
                obs[tag].append(sum(read(p)['seconds'] for p in pair))
                if c.get('arithmetic'):
                    f = Path(c['path']).parent
                    svd[tag].append(read(f.with_name(f.name+'-svd')/'receipt.json')['seconds'])
            require(obs[tag] and svd[tag], 'missing actual observation or SVD cost: '+tag)
        commands = read(root/'command-receipts.json')
        completed = [c for c in commands if c['label'] in
                     [name+'-'+a.policy for name in ('pilot', 'primary', 'transactions', 'extend')]]
        require(any(c['label'] == 'pilot-'+a.policy and c['exit_code'] == 0 for c in completed),
                'whole pilot process receipt must complete before projection')
        for c in completed:
            from wj_immutable import file_hash
            require(file_hash(c['log']) == c['log_sha256'], 'changed completed phase log')
        children = list(folder.rglob('receipt.json'))
        child_seconds = sum(read(p)['seconds'] for p in children)
        whole_seconds = sum(c['seconds'] for c in completed)
        observed_pairs = sum(len(t['trajectory'])+1 for t in s['trajectories'].values())
        # Includes parent quadrature/algebra, hashing, output parsing, pairing,
        # and other I/O outside monitored child processes. No closure reserve
        # is substituted for this measured work allowance.
        parent_seconds = max(0., whole_seconds-child_seconds)
        parent_per_pair = parent_seconds/observed_pairs
        remaining = {tag: max(0, int(target[tag])-int(t['trajectory'][-1]['steps'])) for tag, t in s['trajectories'].items()}
        nobs = dict(coarse=3, fine=3, baseline=2) if a.phase == 'Tx' else dict(coarse=3, fine=3, baseline=3)
        nsvd = dict(coarse=2, fine=2, baseline=2)
        tx = 0 if s['components']['transactions']['status'] == 'PASS' else 9
        parent_bound = sum(nobs.values())*parent_per_pair*headroom
        estimate = sum(remaining[tag]*rates[tag]+nobs[tag]*max(obs[tag])*headroom+
                       nsvd[tag]*max(svd[tag])*headroom for tag in remaining)
        estimate += parent_bound+tx*rates['coarse']+(60 if tx else 0)+600
        totals = ledger.totals()
        fits = totals['attempts']+sum(remaining.values())+tx <= 800 and totals['seconds']+estimate <= 10800
        write_once(archive/'ledger-basis.json', read(ledger.file))
        write_once(archive/'completed-commands.json', completed)
        value = dict(schema='WJ6MeasuredResourceProjectionV2', phase=a.phase, status='PASS' if fits else 'RESOURCE_LIMIT',
            parent_resource=file_identity(folder/'resource-design.json'), target=target,
            native_seconds_per_step_bound=max(rates.values()), native_seconds_per_step_bounds=rates,
            measured_observation_seconds=obs, measured_svd_seconds=svd, remaining_by_tag=remaining,
            remaining_native_attempts=sum(remaining.values()), transaction_attempts_held=tx,
            remaining_observation_allowance=nobs, remaining_svd_allowance=nsvd,
            measured_whole_phase_seconds=whole_seconds, measured_child_process_seconds=child_seconds,
            measured_parent_algebra_verification_IO_seconds=parent_seconds,
            parent_seconds_per_observation_pair=parent_per_pair, remaining_parent_seconds_bound=parent_bound,
            headroom=headroom, closure_seconds_held=600, projected_remaining_seconds=estimate,
            work_totals=totals, actual_policy_bytes=volume(folder), resource=DEFAULT,
            source=artifacts([Path(__file__)]), basis=artifacts([archive/'ledger-basis.json', archive/'completed-commands.json', folder/'design.json']))
        write_once(archive/'projection.json', value)
        target_path = folder/(a.phase+'-resource-projection.json')
        write_once(target_path, value)
        s['components']['resource_projection' if a.phase == 'Tx' else 'extension_resource'] = component(value['status'], projection=file_identity(target_path))
        atomic(folder/'summary.json', s)
        print(json.dumps(value))
    finally:
        ledger.finish(label, 0, time.time()-started, str(archive/'projection.json'))


if __name__ == '__main__':
    main()
