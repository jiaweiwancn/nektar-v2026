#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Transport exact WJ-6 receipts with compressed/decoded identity bridges."""
import argparse
import gzip
import hashlib
from pathlib import Path
import sys
import time

from wj_immutable import artifacts, file_hash, file_identity, read, verify_artifacts, write_once
from run_cls_epawf_picard import require


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'compact', 'figures', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--seal-label', default='seal-v1')
    args = parser.parse_args()
    started = time.time()
    root, compact = args.root.resolve(), args.compact.resolve()
    publication = read(compact/'publication.json')
    verify_artifacts(publication['files'])
    verify_artifacts(publication['figures'])
    before = read(compact/'command-receipts-before-seal.json')
    after = read(root/'command-receipts.json')
    require(after[:-1] == before and after[-1]['label'] == args.seal_label,
            'unexpected command after pre-seal snapshot')
    seal = after[-1]
    require(seal['exit_code'] == 0 and file_hash(seal['log']) == seal['log_sha256'], 'closing seal receipt changed')
    with (compact/'final-command-receipts.json').open('xb') as stream:
        stream.write((root/'command-receipts.json').read_bytes())
    write_once(compact/'finalization.json', dict(schema='WJ6OuterReceiptBridgeV1',
        publication=file_identity(compact/'publication.json'), seal_command=seal,
        final_command_receipts=file_identity(compact/'final-command-receipts.json'), governing_attempts=0,
        transport_coverage='180-second conservative passive/I/O reserve in closed ledger',
        code=artifacts([Path(__file__)])))
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    identities = {}

    def copy(path, relative):
        data = path.read_bytes()
        compress = len(data) > 500000 and path.suffix in ('.json', '.csv')
        target = out/(relative+('.gz' if compress else ''))
        target.parent.mkdir(parents=True, exist_ok=True)
        encoded = gzip.compress(data, compresslevel=9, mtime=0) if compress else data
        with target.open('xb') as stream:
            stream.write(encoded)
        identities[str(target.relative_to(out))] = dict(original=file_identity(path),
            transport_sha256=file_hash(target), transport_bytes=len(encoded),
            encoding='gzip mtime0' if compress else 'exact bytes',
            uncompressed_sha256=hashlib.sha256(data).hexdigest())

    for path in sorted(compact.iterdir()):
        if path.is_file() and path.suffix != '.lock':copy(path, path.name)
    for path in sorted(args.figures.iterdir()):
        if path.is_file() and path.suffix != '.lock':copy(path, 'figures/'+path.name)
    for policy in ('EP', 'Budget'):
        folder=root/'study'/policy
        for name in ('design.json','resource-design.json','Tx-resource-projection.json',
                     'T3-resource-projection.json','transaction-result.json'):
            path=folder/name
            if path.exists():copy(path,'intents/'+policy+'/'+name)
        for case in sorted(folder.glob('*/case.xml')):
            for name in ('case.xml','requested-intent.json','xml-identity.json','execution-intent.json','receipt.json'):
                path=case.parent/name
                if path.exists():copy(path,'cases/'+policy+'/'+case.parent.name+'/'+name)
    for path in sorted((root/'resource-projections-v1').rglob('*.json')):
        copy(path,'resource-bases/'+str(path.relative_to(root/'resource-projections-v1')))
    for version in sorted(root.glob('row-qualified-v*')):
        for policy in ('EP','Budget'):
            for name in ('resource-plan.json','qualification.json','records.json','failure.json','original-summary.json'):
                path=version/policy/name
                if path.exists():copy(path,'row-qualification/'+version.name+'/'+policy+'/'+name)
    for name in ('producer-source-content.json','prior-seals.json'):
        copy(root/name,'identity/'+name)
    copy(Path(read(compact/'manifest.json')['artifact_index']['path']),'identity/heavy-artifact-index.json')
    for path in root.glob('observer-*-source.json'):copy(path,'identity/'+path.name)
    for path in sorted(root.glob('*-whole-resource.txt')):
        copy(path,'identity/'+path.name)
    for policy in ('EP','Budget'):
        copy(root/'checks-v4/wj6-capacity'/policy/'result.json','controls/'+policy+'/result.json')
    commands = read(compact/'final-command-receipts.json')
    for c in commands:copy(Path(c['log']), 'logs/'+Path(c['log']).name)
    copy(root/'inventory-checks-v4.json', 'logs/inventory.json')
    # Preserve preflight/matching-runtime evidence without copying executables.
    for pattern in ('native-source-equivalence-v*.json', '*-allocation.json', 'allocation-*.json',
                    'verified-runtime-aliases-v4.json'):
        for path in sorted(root.glob(pattern)):
            relative = 'identity/'+path.name
            if relative not in identities:copy(path, relative)
    write_once(out/'transport.json', dict(schema='WJ6ExactTransportV1', files=identities, argv=sys.argv,
        code=artifacts([Path(__file__)]), exit_code=0, seconds=time.time()-started, governing_attempts=0,
        restore='gzip-decode *.gz; decoded SHA-256 equals original/uncompressed_sha256; receipts retain exact original bytes',
        verification_scope='compact bytes can be independently verified after transport; heavy Linux bytes were read and verified remotely, not transported'))
    print(out/'transport.json')


if __name__ == '__main__':
    main()
