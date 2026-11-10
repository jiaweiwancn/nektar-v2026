#!/usr/bin/env python3
"""Read-only compact export of this package's completed Linux artifacts.

Use a new output directory. No trajectory, request, seal or historical receipt is
changed. Scientific failure is represented in decision.json, not hidden by the
successful completion of this file-copy operation.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as f:
        return json.load(f)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--export', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)
    inventory = []

    def copy(path, target=None):
        path = Path(path)
        name = target or path.name
        dest = a.output/name
        dest.parent.mkdir(parents=True, exist_ok=True)
        compress = path.suffix == '.json' and path.stat().st_size > 200000
        if compress:
            dest = dest.with_suffix(dest.suffix+'.gz')
            with dest.open('wb') as raw:
                with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as f:
                    f.write(path.read_bytes())
        else:
            shutil.copy2(path, dest)
        inventory.append(dict(source=str(path), source_sha256=sha(path),
                              file=str(dest.relative_to(a.output)), sha256=sha(dest),
                              encoding='gzip' if compress else 'original'))

    summary = load(a.export/'summary.json')
    complete = load(a.export/'complete-results.json.gz')
    assert summary['complete_results']['sha256'] == sha(a.export/'complete-results.json.gz')
    assert summary['status'] == complete['status']
    for path in sorted(a.export.iterdir()):
        if path.is_file():
            copy(path)
    root_files = [
        'workspace.json', 'exclusive-costs-complete.json', 'numerical-overview.json',
        'export-final-v1-command.json', 'resource-finalization.json',
        'targeted-regressions-final-v6.json', 'ctest-budget-accuracy.json',
        'ctest-final-v4.json', 'ctest-final-v5.json', 'ctest-current.json',
        'build-ctest-descriptor.json', 'historical-receipt-byte-verification.json',
        'reused-E-duration-endpoint-coverage.json', 'actual-initial-regression.json',
        'real-spatial-gate-export-control.json', 'partial-export-negative-control.json',
        'neutral-as-perturbed-negative-control.json', 'physical-budget-face-regression.json',
        'physical-budget-accuracy-command.json', 'forward-interface-exact-scope-status.json',
        'forward-interface-tolerance-status.json', 'resource-census-exception.json',
        'A95-Accepted-complete-response-audit-process.json', 'A95-EP-response-selection.json',
        'ancillary-A95-results.json', 'A95-EP-preparation.json', 'A95-EP-preparation-v2.json',
        'A95-response-volume-ownership-process.json', 'A95-EP-response-volume-ownership-process.json',
        'complete-response-budgets-A95.json', 'remote-pre-push-inspection.json',
        'build-probe-exact-scope-control.json', 'build-probe-side-control.json',
        'build-budget-face-control.json', 'common-partition-regression.json',
    ]
    for name in root_files:
        copy(a.root/name, 'controls/'+name)
    for folder in ('aggregate-negative-controls', 'complete-export-negative-controls',
                   'collector-identity-negative-controls'):
        copy(a.root/folder/'results.json', 'controls/'+folder+'.json')
    for path in sorted(a.root.glob('*.py')):
        copy(path, 'linux-controllers/'+path.name)
    for path in sorted(a.root.glob('*.log')):
        if path.stat().st_size < 40000:
            copy(path, 'logs/'+path.name)
    for name in ('case.xml', 'request.json', 'packet.json', 'results.json', 'probe.log'):
        copy(a.root/'minimal-volume-ownership-reproducer'/name, 'minimal-ownership/'+name)
    for folder in ('volume-ownership', 'physical-budget-accuracy', 'common-inner-qualification'):
        for path in sorted((a.root/folder).glob('*/results.json')):
            copy(path, folder+'/'+path.parent.name+'.json')
    for name in ('A93-N-Accepted-loader-v2.json', 'A95-N-Accepted-loader-v2.json',
                 'A95-N-EP-loader-v3.json'):
        copy(a.root/'admissions'/name, 'admissions/'+name)
    for folder in ('half-amplitude-comparison', 'wall-response-comparison'):
        copy(a.root/folder/'results.json', folder+'.json')
    for path in sorted((a.root/'physical-budget-analyses').glob('*/results.json')):
        copy(path, 'physical-budget-analyses/'+path.parent.name+'.json')

    costs = load(a.root/'exclusive-costs-complete.json')
    decision = dict(schema='NonlinearDurationResponseDecisionV1',
        status=complete['status'], scientific_process_exit=load(a.root/'export-final-v1-command.json')['exit_code'],
        starting_revision=complete['starting_revision'],
        numerical_export_source_revision=complete['source_report'],
        export_packager_source_revision=subprocess.check_output(
            ['git', '-C', str(a.root/'source'), 'rev-parse', 'HEAD'], text=True).strip(),
        required_scope_classification=complete['required_scope_classification'],
        cohorts={k:dict(numerical_status=v['numerical_status'], classification=v['classification'])
                 for k,v in complete['cohorts'].items()},
        response_statuses={k:v['status'] for k,v in complete['responses'].items()},
        all_complete_comparison_bodies='complete-results.json.gz',
        failures=complete['failures'],
        costs=dict(segments=len(costs['cases']),
                   accepted_steps=sum(v['accepted_steps_in_segment'] for v in costs['cases']),
                   instrumented_exclusive_process_CPU_seconds=sum(
                       v['instrumented_exclusive_process_cpu_seconds'] for v in costs['cases'])),
        scope='serial flat2D prescribed-reference patch; failed response audits do not qualify a longer response',
        packager_status='completed read-only export; this status does not override scientific gates')
    (a.output/'decision.json').write_text(json.dumps(decision, indent=2)+'\n')
    (a.output/'artifact-inventory.json').write_text(json.dumps(inventory, indent=2)+'\n')
    print(json.dumps(dict(files=len(inventory), status=decision['status'], output=str(a.output))))


if __name__ == '__main__':
    main()
