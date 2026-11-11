#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Seal all field/scope sensitivity margins for the adverse WJ-5 dt control.

This is passive reprocessing of the immutable compact publication. A diagnostic
comparison passing its sensitivity budget does not qualify its adverse states.
"""
import argparse
import csv
from pathlib import Path
import sys
import time

from run_cls_epawf_picard import require
from wj4_qualification import sensitivity_sequence
from wj5_measure import boundary_envelope
from wj_immutable import (Ledger, artifacts, file_hash, file_identity, read,
                          verify_artifacts, write_once)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--compact', type=Path, required=True)
    args = parser.parse_args()
    root, compact = args.root.resolve(), args.compact.resolve()
    output = compact / 'control-sensitivity.json'
    require(not output.exists(), 'control supplement is write-once')
    inputs = artifacts([compact / name for name in
                        ('EP.json', 'Budget.json', 'manifest.json', 'receipt.json')])
    verify_artifacts(read(compact / 'receipt.json')['files'])
    code = artifacts([Path(__file__), *[Path(sys.modules[name].__file__) for name in
        ('wj5_measure', 'wj4_qualification', 'wj_immutable')]])
    label = 'publication/control-sensitivity-v2'
    ledgers = {policy: Ledger(root, policy, cap=800) for policy in ('EP', 'Budget')}
    for ledger in ledgers.values():
        ledger.reserve(label, dict(inputs=inputs, code=code), 0, 240)
    started = time.time()
    packets, rows = {}, []
    try:
        for policy in ledgers:
            data = read(compact / (policy + '.json'))
            control = data['adverse_dt_control']
            startup = control['fresh_startup']
            require(all(startup[key] is True for key in
                        ('exact_primary_gauge', 'exact_physical_histories',
                         'native_dt_labels_differ')), 'incompatible fresh control startup')
            coarse = dict(initial=data['retained']['retained-Nx2-0'],
                          trajectory=[data['retained']['retained-Nx2-32']])
            fine = dict(initial=control['initial'], trajectory=[control['terminal']])
            trajectories = dict(coarse=coarse, fine=fine)
            for trajectory in trajectories.values():
                for record in [trajectory['initial']] + trajectory['trajectory']:
                    for key in ('report', 'checkpoint'):
                        path = Path(record[key])
                        if key == 'checkpoint':
                            path = path / 'State.bin'
                        require(file_hash(path) == record[key + '_sha256'],
                                'changed control producer bytes: ' + str(path))
            a, b = coarse['trajectory'][0], fine['trajectory'][0]
            require(abs(a['time'] - b['time']) <= 2e-15,
                    'control states do not share the frozen physical target')
            floors, boundary, coverage = boundary_envelope(trajectories)
            comparison = {'Tx': control['comparison']}
            gate = sensitivity_sequence(comparison, floors, data['frozen_design'], ['Tx'],
                scope_uncertainty={'constrained_boundary': boundary})
            require('metadata' not in gate['components'], 'incomplete control sensitivity')
            for quantity, result in gate['components']['Tx']['components'].items():
                rows.append(dict(policy=policy, quantity=quantity, status=result['status'],
                    value=result['value'], uncertainty=result['uncertainty'], limit=result['limit'],
                    margin=result['limit'] - result['value'] - result['uncertainty'],
                    fraction=(result['value'] + result['uncertainty']) / result['limit']))
            packets[policy] = dict(gate=gate, boundary_coverage=coverage,
                global_floors=floors, constraint_boundary=boundary,
                engineering_outcome=control['engineering_result'],
                interpretation='adverse diagnostic only; no selected-mesh qualification or order claim')
        table = compact / 'control-sensitivity.csv'
        with table.open('x', newline='') as stream:
            writer = csv.DictWriter(stream, list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        verify_artifacts(inputs)
        write_once(output, dict(schema='WJ5AdverseControlSensitivityV1', policies=packets,
            input_lineage=inputs, code=code, table=file_identity(table),
            governing_attempts=0, measured_seconds=time.time() - started))
    finally:
        for ledger in ledgers.values():
            ledger.finish(label, 0, time.time() - started, str(output))
    print(output)


if __name__ == '__main__':
    main()
