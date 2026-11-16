#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact, reproducible WP2 native execution/cost index (read only)."""
import argparse
import collections
import os
import statistics
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
from cao_tamura_h2d_state import load
from run_cls_epawf_picard import sha, require, THREADS
from continue_cls_epawf_pimple import atomic
from wp2_execution_index import executions as native_executions


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = args.root
    executions = native_executions(root)
    groups = collections.defaultdict(list)
    commands = []
    for name, entry in sorted(executions.items()):
        report = Path(name)
        require(not report.exists() or sha(report) == entry['report_sha256'], 'native report changed')
        groups[report.relative_to(root).parts[0]].append(entry)
        xml = report.with_name('case.xml')
        require(sha(xml) == entry['xml_sha256'], 'native case changed')
        progress = load(report.parent.parent / 'progress.json')
        solver = progress['solver']
        commands.append(dict(cwd=str(report.parent), environment=dict(THREADS,
            LD_LIBRARY_PATH=str(Path(solver).parent / 'lib'), NEKTAR_RANS_TIMING_FILE=str(report.parent / 'timing.jsonl')),
            argv=['/usr/bin/time', '-f', '%e %U %S %M', '-o', str(report.parent / 'usage.txt'), solver, str(xml)],
            report=str(report), report_sha256=entry['report_sha256'], xml_sha256=entry['xml_sha256'],
            process_returncode=entry['code'], failure=entry['failure'], requested_attempt_work=entry['requested_attempt_work'],
            actual_attempts=entry['observed_physical_attempts'], accepted_steps=entry['accepted_step_executions']))
    costs = {}
    for group, entries in groups.items():
        exclusive = collections.defaultdict(float)
        steps = []
        first_steps = []
        for entry in entries:
            local_steps = []
            timing = Path(entry['report']).parent / 'timing.jsonl'
            for line in timing.read_text().splitlines() if timing.exists() else []:
                import json
                row = json.loads(line)
                exclusive[row['name']] += row['exclusive_process_cpu_seconds']
                if row['name'] == 'native.step':
                    local_steps.append(row['wall_seconds'])
            if local_steps:
                first_steps.append(local_steps[0])
                steps.extend(local_steps[1:])
        usage = [[float(s) for s in e['usage'].splitlines()[-1].split()] for e in entries]
        costs[group] = dict(native_invocations=len(entries), observed_attempts=sum(e['observed_physical_attempts'] for e in entries),
            accepted_step_executions=sum(e['accepted_step_executions'] for e in entries),
            driver_wall_seconds=sum(e['seconds'] for e in entries), process_CPU_seconds=sum(u[1]+u[2] for u in usage),
            peak_RSS_KiB=max(u[-1] for u in usage),
            first_native_step_median_seconds=statistics.median(first_steps) if first_steps else None,
            continued_native_steps=len(steps), continued_step_median_seconds=statistics.median(steps) if steps else None,
            major_exclusive_process_CPU=sorted(exclusive.items(), key=lambda item: -item[1])[:15])
    files = ['origin.json', 'aggregate.json', 'native-source-chain.json', 'runtime-v12/identity.json',
        'probe-runtime-v12/identity.json', 'probe-runtime-analysis-rows/identity.json', 'final-status.json',
        'final-readonly-completion.json', 'recovery-standard-native-rejection.json',
        'unperturbed-preservation-v2.json', 'dense-quarter-budget-summary.json', 'native-boundary-comparison-bounds.json',
        'unperturbed-all-pairs-windows/summary.json', 'terminal-checks-after-transactions.json',
        'recovery-CLS-followup-T4p2/progress.json', 'recovery-standard-followup-T4p2/progress.json',
        'recovery-CLS-followup-T4p2/assessment/result.json', 'recovery-CLS-final-all-pairs/CLS-recovery.json',
        'recovery-timestep-comparison-T03.json', 'recovery-figure-final/provenance.json']
    evidence = {}
    for name in files:
        path = root / name
        if path.exists():
            evidence[name] = dict(path=str(path), sha256=sha(path), bytes=path.stat().st_size)
    atomic(args.output, dict(version='WP2ArtifactCostIndexV1', capture_Asia_Shanghai=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),
        root=str(root), costs=costs, native_commands=commands, evidence=evidence, source_sha256=sha(Path(__file__)),
        cost_semantics='driver wall includes startup, native stepping and report/checkpoint/output; native.step is inclusive; major costs are exclusive process CPU. output.report includes read-only budget construction. Reused prefixes appear once by native report path; independent replay executions remain separate.',
        resources=dict(native_ranks=1, numerical_library_threads=1, build_workers_max=4, aggregate_workers_max=72)))


if __name__ == '__main__':
    main()
