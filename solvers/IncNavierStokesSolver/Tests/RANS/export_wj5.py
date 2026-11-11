#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Export exact WJ-5 seals and compact, hash-bridged compressed observations.

This final transport performs no governing or observation solve. Its I/O is
covered by the explicit passive-inspection reserve in the closed work ledger.
Original receipt bytes and heavy artifacts are never rewritten.
"""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import sys
import time

from wj_immutable import artifacts, file_hash, file_identity, verify_artifacts, write_once


def require(ok, message):
    if not ok:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--compact', type=Path, required=True)
    parser.add_argument('--figures', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seal-label', default='seal-v1')
    args = parser.parse_args()
    started = time.time()
    root, compact, figures = args.root.resolve(), args.compact.resolve(), args.figures.resolve()
    load = lambda path: json.loads(path.read_text())
    publication = load(compact / 'publication.json')
    verify_artifacts(publication['files'])
    verify_artifacts(publication['figures'])
    before = load(compact / 'command-receipts-before-seal.json')
    after = load(root / 'command-receipts.json')
    require(after[:-1] == before and after[-1]['label'] == args.seal_label,
            'unexpected commands after the immutable pre-seal snapshot')
    seal = after[-1]
    require(seal['exit_code'] == 0 and file_hash(seal['log']) == seal['log_sha256'],
            'final seal command did not pass with unchanged log bytes')
    with (compact / 'final-command-receipts.json').open('xb') as stream:
        stream.write((root / 'command-receipts.json').read_bytes())
    write_once(compact / 'finalization.json', dict(schema='WJ5OuterReceiptBridgeV1',
        publication=file_identity(compact / 'publication.json'), seal_command=seal,
        final_command_receipts=file_identity(compact / 'final-command-receipts.json'),
        code=artifacts([Path(__file__), Path(sys.modules['wj_immutable'].__file__)]),
        governing_attempts=0, coverage='post-seal transport I/O covered by120-second conservative passive reserve'))
    runs = list(csv.DictReader((compact / 'run-accounting.csv').open()))
    policies = {}
    for policy in ('EP', 'Budget'):
        data = load(compact / (policy + '.json'))
        attempts = sum(int(float(row['attempts'])) for row in runs if row['policy'] == policy)
        require(attempts == 128, 'unexpected bounded new native-attempt count')
        policies[policy] = dict(outcome=publication['policies'][policy],
            new_governing_attempts=attempts, work=publication['accounting'][policy]['totals'],
            screen_accuracy={mesh:dict(unknowns=entry['unknowns'],
                initial=entry['initial']['errors']['accuracy_gate'],
                terminal=entry['terminal']['errors']['accuracy_gate'])
                for mesh, entry in data['screens'].items()},
            adverse_control_accuracy=data['adverse_dt_control']['engineering_result'],
            frozen_T3=.004096, reached_Tx=.001024, selected_mesh=None,
            selected_T3_transactions='NOT_RUN: neither predeclared common local mesh qualifies',
            stationarity='not targeted; WJ-4 NOT_REACHED remains unchanged')
    write_once(compact / 'bounded-summary.json', dict(schema='WJ5BoundedSummaryV1',
        publication=file_identity(compact / 'publication.json'), policies=policies,
        tests={key:load(compact / 'regression-verification.json')[key]
               for key in ('status', 'registered', 'selected_unique', 'names')}))
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    identities = {}

    def copy(path, relative):
        data = path.read_bytes()
        compressed = len(data) > 500_000 and path.suffix in ('.json', '.csv')
        target = output / (relative + ('.gz' if compressed else ''))
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open('xb') as stream:
            stream.write(gzip.compress(data, compresslevel=9, mtime=0) if compressed else data)
        identities[str(target.relative_to(output))] = dict(original=file_identity(path),
            transport_sha256=file_hash(target), transport_bytes=target.stat().st_size,
            encoding='gzip mtime0' if compressed else 'exact bytes',
            uncompressed_sha256=hashlib.sha256(data).hexdigest())

    for path in sorted(compact.iterdir()):
        if path.is_file() and path.suffix != '.lock':
            copy(path, path.name)
    for path in sorted(figures.iterdir()):
        if path.is_file() and path.suffix != '.lock':
            copy(path, 'figures/' + path.name)
    for path in sorted(root.glob('ctest-*.log')):
        copy(path, 'logs/' + path.name)
    for name in ('wj5-focused-final.log', 'measurement-metadata-fixed.log',
                 'publication-metadata-orders.log'):
        copy(root / name, 'logs/' + name)
    copy(root / 'inventory.log', 'logs/inventory.json')
    write_once(output / 'transport.json', dict(schema='WJ5ExactTransportV1',
        files=identities, argv=sys.argv, code=artifacts([Path(__file__)]),
        exit_code=0, seconds=time.time() - started, governing_attempts=0,
        omitted='empty flock coordination files remain indexed in heavy storage',
        restore='gzip-decode *.gz; decoded SHA-256 must equal original/uncompressed_sha256; all receipts are exact original bytes'))
    print(output / 'transport.json')


if __name__ == '__main__':
    main()
