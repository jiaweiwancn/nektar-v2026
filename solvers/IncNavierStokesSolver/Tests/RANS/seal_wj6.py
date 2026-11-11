#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Verify WJ-6 completed regressions, heavy identities and closed accounting."""
import argparse
import json
from pathlib import Path
import re
import time

from wj_immutable import Ledger, artifacts, file_hash, file_identity, read, verify_artifacts, write_once
from wj_qualification import aggregate, component
from run_cls_epawf_picard import require
from wj6_resources import DEFAULT,volume


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'compact', 'figures'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    started = time.time()
    root, out = args.root.resolve(), args.compact.resolve()
    commands = read(root/'command-receipts.json')
    by_label = {c['label']: c for c in commands}
    for c in commands:
        require(file_hash(c['log']) == c['log_sha256'], 'command log changed: '+c['label'])
    inventory = read(root/'inventory-checks-v4.json')
    registered = {t['name']: t for t in inventory['tests']}
    selections, unique = [], set()
    for label, count in [('wj6-regression-controls', 24), ('wj6-regression-native', 3),
                         ('wj6-regression-readonly', 2), ('wj6-regression-metadata', 1),
                         ('wj6-regression-composite', 2), ('legacy-local-controls', 3)]:
        inv = read(root/(label+'-inventory.json'))
        c = by_label[label if label == 'legacy-local-controls' else 'ctest-'+label]
        require(len(inv['tests']) == count and c['exit_code'] == 0, 'missing/failed selection: '+label)
        require('100% tests passed' in Path(c['log']).read_text(), 'CTest did not complete: '+label)
        names = [t['name'] for t in inv['tests']]
        unique.update(names)
        selections.append(dict(label=label, names=names, count=count, receipt=c,
                               inventory=file_identity(root/(label+'-inventory.json'))))
    require(len(unique) == 35, 'required WJ-5 regression union is not35')
    for label, names in [('capacity-EP', ['WJDenseCapacityEP']),
                         ('capacity-Budget', ['WJDenseCapacityBudget', 'WJDenseCapacityMetadata'])]:
        c = by_label[label]
        require(c['exit_code'] == 0 and '100% tests passed' in Path(c['log']).read_text(), 'capacity controls failed')
        for name in names:
            require(name in registered, 'missing capacity registration: '+name)
        unique.update(names)
        selections.append(dict(label=label, names=names, count=len(names), receipt=c))
    require(len(unique) == 38, 'wrong actual unique selected count')
    latest = by_label['metadata-final-v7']
    require(latest['exit_code'] == 0 and 'Ran 6 tests' in Path(latest['log']).read_text(), 'final metadata fixtures incomplete')
    supplement=by_label['observation-metadata-v2']
    require(supplement['exit_code']==0 and 'Ran 7 tests' in Path(supplement['log']).read_text(), 'independent-order metadata fixtures incomplete')
    resume=by_label['transaction-resumption-metadata-v3']
    require(resume['exit_code']==0 and 'Ran 8 tests' in Path(resume['log']).read_text(), 'zero-launch resumption metadata fixtures incomplete')
    observation_resume=by_label['observation-resume-metadata-v4']
    require(observation_resume['exit_code']==0 and 'Ran 9 tests' in Path(observation_resume['log']).read_text(), 'observation reuse metadata fixtures incomplete')
    write_once(out/'regression-verification.json', dict(schema='WJ6RegressionSealV1', status='PASS',
        registered=len(registered), selected_unique=len(unique), names=sorted(unique), selections=selections,
        latest_metadata=latest,independent_order_metadata=supplement,resumption_metadata=[resume,observation_resume],
        allocation_denials=[file_identity(f) for f in root.glob('*allocation*.json')
                                                   if read(f).get('allowed') is False]))
    manifest = read(out/'manifest.json')
    verify_artifacts(manifest['prior_seals']['wj4'])
    verify_artifacts(manifest['prior_seals']['wj5'])
    heavy_meta = manifest['artifact_index']
    require(file_hash(heavy_meta['path']) == heavy_meta['sha256'], 'heavy index changed')
    heavy = read(heavy_meta['path'])
    verify_artifacts(heavy)
    # This inventory contains completed immutable cases/observations only.
    # Mutable progress, ledgers and the closing publisher log were excluded.
    write_once(out/'heavy-verification.json', dict(status='PASS', index=heavy_meta,
        verified_artifacts=len(heavy), prior_seals=manifest['prior_seals'],
        meaning='independently read and SHA-256 verified remotely; compact transport verification is separate'))
    ledgers = {p: Ledger(root, p, cap=800) for p in ('EP', 'Budget')}
    accounting = {}
    build = [c for c in commands if c['label'].startswith('build')]
    nonbuild = [c for c in commands if c not in build]
    # Shared controls, observations and setup are conservatively charged to
    # each policy. Per-policy process receipts are charged only to that policy.
    for policy, ledger in ledgers.items():
        data = read(ledger.file)
        require(all(e['status'] != 'RUNNING' for e in data['entries'].values()) and
                all(s['status'] != 'RUNNING' for s in data.get('sessions', {}).values()),
                'unreconciled reservation/session: '+policy)
        selected = [c for c in nonbuild if not any(c['label'].endswith('-'+p) for p in ledgers) or
                    c['label'].endswith('-'+policy)]
        desired = sum(c['seconds'] for c in selected)+180.
        before = ledger.totals()
        extra = max(0., desired-before['seconds'])+max(1., time.time()-started)
        label = 'final/outer-process-and-passive-coverage'
        ledger.reserve(label, dict(commands=file_identity(root/'command-receipts.json'), labels=[c['label'] for c in selected],
            passive_and_transport_seconds=180., definition='outer process/session coverage; build time separate; no existing charge reduced'), 0, max(1., extra))
        ledger.finish(label, 0, extra, str(out/'publication.json'))
        totals = ledger.totals()
        require(totals['attempts'] <= 800 and totals['seconds'] <= 10800, 'RESOURCE_LIMIT: final package cap')
        snapshot = out/('final-ledger-'+policy+'.json')
        write_once(snapshot, read(ledger.file))
        accounting[policy] = dict(totals=totals, ledger=file_identity(snapshot), outer_process_seconds=desired,
                                  added_seconds=extra, build_seconds_excluded=sum(c['seconds'] for c in build))
    write_once(out/'command-receipts-before-seal.json', commands)
    policies = {}
    for policy in ledgers:
        s = read(out/(policy+'.json'))
        parts = s['components']
        scopes=[root/'study'/policy,root/'checks-v4/wj6-capacity'/policy]+[v/policy for v in root.glob('row-qualified-v*')]
        peaks=[]
        for scope in scopes:
            for f in scope.rglob('receipt.json'):
                data=read(f)
                if data.get('resource'):peaks.append(data['resource']['peak_rss_bytes'])
        whole=[]
        for f in root.glob('*-'+policy+'-whole-resource.txt'):
            match=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',f.read_text())
            require(match,'missing completed whole-phase resource receipt: '+str(f))
            peaks.append(1024*int(match.group(1)));whole.append(file_identity(f))
        policy_bytes=sum(volume(scope) for scope in scopes)
        actual=component('PASS' if max(peaks)<=DEFAULT['rss_bytes'] and policy_bytes<=DEFAULT['policy_output_bytes'] else 'RESOURCE_LIMIT',
            totals=accounting[policy]['totals'],peak_measured_rss_bytes=max(peaks),policy_data_bytes=policy_bytes,
            limits=DEFAULT,whole_phase_receipts=whole,scope='measured native/passive child and whole-phase maxima; source/build storage separate')
        primary = aggregate({k: parts.get(k, component('NOT_RUN')) for k in
                             ('pilot', 'resource_projection', 'primary', 'transactions')}|dict(actual_resources=actual))
        policies[policy] = dict(primary=primary,
            extension=parts.get('extension', component('NOT_RUN')),
            extension_resource=parts.get('extension_resource', component('NOT_RUN')),
            actual_resource=actual)
    write_once(out/'publication.json', dict(schema='WJ6FinalPublicationV1', inner_receipt=file_identity(out/'receipt.json'),
        files=artifacts([f for f in out.iterdir() if f.is_file() and f.suffix != '.lock']),
        figures=artifacts([f for f in args.figures.iterdir() if f.is_file()]), policies=policies, accounting=accounting,
        native_build_commands=build, native_build_seconds=sum(c['seconds'] for c in build), governing_attempts=0,
        code=artifacts([Path(__file__)]),
        scope='bounded serial flat2D dense finite reference-Dirichlet strip; no stationarity, exact discrete equilibrium, global stability/realizability or scalable production claim',
        closing_receipt='exact outer seal receipt is bound after this process exits; original bytes are never rewritten'))
    print(out/'publication.json')


if __name__ == '__main__':
    main()
