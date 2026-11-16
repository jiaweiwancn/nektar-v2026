#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Hash-bound CT16 investigation ledger; unknown qualification never passes."""
import argparse, hashlib, json, math, os, subprocess, time
from decimal import Decimal, ROUND_CEILING
from pathlib import Path
from ct_accepted15_collect import compact, logical_size, max_address_space

START = 1794835825
GIB = 1024**3


def read(path):
    return json.loads(path.read_text()) if path.is_file() else None


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024**2), b''):
            h.update(block)
    return h.hexdigest()


def write(path, value):
    with path.open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')


def count(T, dt):
    return int((T / Decimal(dt)).to_integral_value(rounding=ROUND_CEILING))


def fit(values):
    return dict(startup=max(0., values[0]),
                step=max(0., values[1]-values[0], (values[2]-values[0])/2,
                         values[2]-values[1]), training=values)


def forecast(root, training):
    q = read(training / 'results.json')
    if not q or q.get('status') != 'pass' or not all(str(i) in q.get('cases', {}) for i in (0, 1, 2)):
        return dict(status='UNATTEMPTED', training=str(training))
    cases = q['cases']
    if any(cases[str(i)].get('detailed_timing', True) or cases[str(i)].get('broad_diagnostic_timing', True) for i in (0, 1, 2)):
        raise ValueError('diagnostic timings cannot be acceptance training')
    model = {k: fit([cases[str(i)]['execution'][field] for i in (0, 1, 2)])
             for k, field in (('CPU', 'process_cpu_seconds'), ('wall', 'wall_seconds'))}
    output = fit([cases[str(i)]['output_bytes'] for i in (0, 1, 2)])
    model['output_all_retained'] = output
    rows = []
    horizons = [('entry', Decimal('.407'))] + [(str(f)+'transit', Decimal(str(f))*Decimal('40.6216378576')) for f in (.1, 1, 3)]
    for name, T in horizons:
        counts = [count(T, dt) for dt in ('.00025', '.000125')]
        for start in ('fresh_zero', 'own_T0.002_conditional'):
            additional = counts if start == 'fresh_zero' else [counts[0]-8, counts[1]-16]
            segments = [math.ceil(n/32) for n in additional]
            N, S = sum(additional), sum(segments)
            initial = 2 if start == 'fresh_zero' else 0
            native = {k: (S+initial)*v['startup']+N*v['step'] for k, v in model.items() if k in ('CPU', 'wall')}
            allowance = (24000/1.25-(S+initial)*model['CPU']['startup'])/N
            rows.append(dict(horizon=name, requested_T=str(T), from_state=start,
                actual_endpoints=[str(Decimal(n)*Decimal(dt)) for n, dt in zip(counts, ('.00025', '.000125'))],
                total_member_steps=counts, additional_member_steps=additional,
                own_start_is_qualified=False, segment_length=32, native_segment_hard_cap=64,
                segments=segments, separate_zero_initializations=initial,
                conditional_native_only_CPU_with_headroom=1.25*native['CPU'],
                conditional_native_only_wall_with_headroom=1.25*native['wall'],
                optimistic_native_step_CPU_allowance_if_all_other_costs_zero=allowance,
                native_only_remaining_reduction_ratio=1.25*native['CPU']/24000,
                all_retained_native_output_bytes_with_headroom=1.25*((S+initial)*output['startup']+N*output['step']),
                complete_campaign_CPU_with_headroom=None, complete_campaign_wall_with_headroom=None,
                qualified_retention_peak_bytes=None, admitted=False,
                status='UNADMITTED', physical_preservation='UNATTEMPTED'))
    return dict(status='UNVALIDATED_PARTIAL_NATIVE_EXTRAPOLATION', training=str(training),
        training_results_sha256=sha(training/'results.json'), model=model, horizons=rows,
        measured_training_AS_bytes=max_address_space(q),
        headroom=1.25, complete_workflow=False, heldout_validated=False,
        missing=['all independent per-state audits', 'full audit/action milestones',
                 'real retention/recovery and replay cost', 'controllers, collection, transfers and publication',
                 'fresh same-component scientific pair and later-state heldout validation'],
        interpretation='Short 0/1/2 envelope and once-only headroom are conditional extrapolations, not mathematical lower bounds on all future implementations or states. Native-only costs do not admit a complete campaign. Own T=.002 starts require fresh qualification of this component. Physical T=.407 is UNATTEMPTED.')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--active-control', type=Path, required=True)
    a = p.parse_args(); root = a.root.resolve(); out = a.output.resolve()
    if out.parent != root or out.exists():
        raise ValueError('collector output must be a new directly owned directory')
    out.mkdir()
    receipts, unclosed, failures = {}, [], []
    peak = written = 0
    for folder in sorted(root.glob('resources-*')):
        if folder.resolve() == a.active_control.resolve():
            continue
        f = folder/'execution.json'; q = read(f)
        if q:
            receipts[folder.name] = dict(sha256=sha(f), **q)
            if q['resource'].get('aggregate_status') != 'pass' or q['resource'].get('project_status') != 'pass':
                failures.append(dict(family=folder.name, resource=q['resource'], scope='includes deliberate paused-guard negative controls; never silently treat incomplete audits as PASS'))
        elif folder.resolve() != a.active_control.resolve():
            unclosed.append(folder.name)
        watch = read(folder/'package-budget-watch.json') or read(folder/'live-budget.json')
        if watch:
            s = watch.get('last_sample', {})
            peak = max(peak, s.get('peak_sampled_live_bytes', 0))
            written += s.get('sampled_cumulative_process_wchar', 0)
    if unclosed:
        raise ValueError('unclosed numerical/engineering family: '+str(unclosed))
    evidence, inventory = {}, {}
    for folder in sorted(root.iterdir()):
        if not folder.is_dir() or folder.name in ('source', 'build', 'tmp') or folder == out or folder.resolve() == a.active_control.resolve():
            continue
        for name in ('request-manifest.json', 'completed-artifacts.json', 'results.json', 'progress.json'):
            f = folder/name
            if f.is_file():
                evidence[str(f.relative_to(root))] = dict(sha256=sha(f), data=compact(read(f)))
        for path, dirs, files in os.walk(folder):
            for name in sorted(files):
                f = Path(path)/name
                if f.is_symlink():
                    inventory[str(f.relative_to(root))] = dict(symlink=os.readlink(f))
                elif f.is_file():
                    inventory[str(f.relative_to(root))] = dict(bytes=f.stat().st_size, sha256=sha(f))
    forecasts = {d.name: forecast(root, d) for d in sorted(root.glob('training-*')) if (d/'results.json').is_file()}
    history = root/'source/docs/plans/cao-tamura-accepted-15-evidence/closure-v2/closed-publication-accounting.json'
    source = subprocess.check_output(['git', '-c', 'index.threads=1', 'rev-parse', 'HEAD'], cwd=root/'source', text=True).strip()
    CPU = sum(q['process_cpu_seconds'] for q in receipts.values())
    gates = dict(G1='UNRESOLVED_AFFECTED_SELECTED_CONTROLS_NOT_COMPLETE',
        G2='UNATTEMPTED_FRESH_CANDIDATE_PAIR', G3='UNRESOLVED_NO_COMPLETE_HELDOUT_MODEL',
        G4='UNATTEMPTED_SELECTED_NATIVE_RETIREMENT_CYCLES', G5='UNADMITTED',
        G6='PREPARED_PENDING_CLOSED_CURRENT_RECEIPT_AND_FINAL_PUSH')
    result = dict(schema='CT16CompactInvestigationV1', status='UNRESOLVED', gates=gates,
        external_root=str(root), source_revision=source,
        prior_CT15=dict(status='CLOSED_UNRESOLVED', accounting_sha256=sha(history), accounting=read(history)),
        engineering_ledger=dict(fixed_start_unix=START, elapsed_wall_seconds=time.time()-START,
            closed_family_CPU=CPU, publication_CPU_reserve=240,
            upper_CPU_before_current_collector=CPU+240, logical_retained_bytes=logical_size(root),
            peak_sampled_live_bytes=peak, sampled_cumulative_process_wchar_lower_bound=written,
            writes_scope='Per-family sampled process wchar,including pipes; can miss short-lived writes. Not exact physical disk writes.',
            receipts=receipts, unclosed_families=unclosed, incomplete_or_failed_worker_audits=failures,
            current_collector_receipt_pending=str(a.active_control)),
        future_ledger=forecasts, physical_T0_407='UNATTEMPTED', full_transit='UNQUALIFIED',
        qualification_scope='Small controls and local timing benefits cannot close selected-space G1,G2,complete G3/G4 or G5. Baseline CT15 short equivalence/physics remains separately valid by exact scope; never relabel its histories as CT16.')
    write(out/'results.json', result); write(out/'forecasts.json', forecasts)
    write(out/'scientific-evidence.json', evidence); write(out/'artifact-inventory.json', inventory)
    write(out/'seal.json', dict(files={f.name: sha(f) for f in out.iterdir() if f.is_file()}))
    print(json.dumps(dict(status=result['status'], closed_CPU=CPU, forecasts=list(forecasts), payloads=len(inventory))))


if __name__ == '__main__':
    main()
