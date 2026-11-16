#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Focused rejection checks for immutable archived-k diagnostic replay."""
import argparse
from copy import deepcopy
from pathlib import Path
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from diagnose_cls_accepted_covariance import run_probe
from run_cls_epawf_picard import require, sha


def main():
    parser = argparse.ArgumentParser()
    for key in ('root', 'probe', 'output'):
        parser.add_argument('--' + key, type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False)
    folder = args.root / 'rejected-k-witness-balance/native'
    baseline, result, receipt = (load(folder / name) for name in ('input.json', 'result.json', 'receipt.json'))
    require(result['operation'] == 'archived_k_balance' and result['initialization_work'] == [0, 0] and
            result['governing_solves'] == result['physical_commits'] == 0, 'positive archived replay is unqualified')
    require(receipt['exit_code'] == 0 and receipt['output_sha256'] == sha(folder / 'result.json'),
            'positive replay receipt changed')
    parent = Path(baseline['native_checkpoint']) / 'State.bin'
    parent_hash = sha(parent)
    cases = {}
    for key, reason in (('expected_history_sha256', 'archived callback history SHA256 differs'),
                        ('native_signature', 'native replay physical signature differs from archived report')):
        packet = deepcopy(baseline)
        packet[key] = '0' * 64
        log, check = run_probe(args.probe, Path(receipt['command'][1]), packet, args.output / key, success=False)
        require(reason in log, 'negative replay failed for an unexpected reason: ' + key)
        require(sha(parent) == parent_hash, 'negative replay changed the accepted parent')
        cases[key] = dict(expected_reason=reason, receipt=check, parent_checkpoint_preserved=True)
    atomic(args.output / 'result.json', dict(version='WP2ArchivedKBalanceContractsV1', passed=True,
        positive_result_sha256=sha(folder / 'result.json'), probe_sha256=sha(args.probe), cases=cases,
        parent_checkpoint_sha256=parent_hash, governing_solves=0, source_sha256=sha(Path(__file__))))


if __name__ == '__main__':
    main()
