#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Collect checked H2k closure evidence without advancing or changing a state."""
import argparse
from datetime import datetime
import json
from pathlib import Path
import re
import time
from zoneinfo import ZoneInfo

from cao_tamura_h2d_state import load
from compare_cls_epawf_pimple import same_history
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import require, sha


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('root', 'comparison', 'output'):
        p.add_argument('--'+key, type=Path, required=True)
    a = p.parse_args()
    a.output.mkdir(exist_ok=False)
    root = a.root
    analysis = load(root/'post-startup-final/analysis.json')
    progress = load(root/'pilot-requalified/progress.json')
    seen = load(root/'seen-reserve/receipt.json')
    fresh = load(root/'fresh-epoch2/receipt.json')
    cert = load(root/'fresh-epoch2/certificate.json')
    imported = load(root/'fresh-epoch2/hill-import.json')
    require(analysis['passed'] and progress['status'] == 'completed_bounded_pair_revalidated', 'incomplete startup/readback')
    require(seen['passed'] and fresh['passed'] and imported['passed'], 'incomplete seed qualification')
    require(progress['seed_sha256'] == cert['state_sha256'] == seen['state_sha256'] == fresh['state_sha256'], 'different physical seeds')
    records = progress['results']
    coarse, fine, split = [load(Path(records[k]['report'])) for k in ('coarse8', 'fine16', 'coarse4-restart')]
    require(same_history(coarse, split), 'direct/split bits differ')
    require(coarse['steps'] == 8 and fine['steps'] == 16, 'incomplete trajectories')
    require(progress['rollback_validation']['passed'], 'rollback not qualified')
    require(sha(Path(cert['H2k']['aggregate_ledger_file'])) == cert['H2k']['aggregate_ledger_sha256'], 'certified ledger changed')
    native_counts = {name: sum(analysis['runs'][k]['counts'][name] for k in ('coarse8', 'fine16', 'coarse4', 'coarse4-restart'))
                     for name in analysis['runs']['coarse8']['counts']}
    require(all(native_counts[k] == 32 for k in ('outer_sweeps', 'momentum_predictor_solves', 'pressure_system_solves',
                'epsilon_transport_solves', 'k_transport_solves', 'bound_calls')), 'native count sum differs')
    rejected_counts = {name: 0 for name in native_counts}
    for label in ('original-rejection-revalidated', 'post-assignment-reject'):
        rejected = load(Path(records[label]['report']))
        require(rejected['failure'] == 'deliberate post-assignment explicit rejection' and rejected['steps'] == 0,
                'unexpected rejected-attempt status')
        for row in rejected['pimple_attempts']:
            for name in rejected_counts:
                rejected_counts[name] += row[name]
    origin = load(root/'origin.json')
    stamp = datetime.now(ZoneInfo('Asia/Shanghai')).isoformat()
    gates = {name: 'PASS' for name in (
        'frozen_physics_and_protected_initial_data', 'useful_native_control_and_old_span_preservation',
        'complete_native_and_seen_sided_inventory', 'fresh_GL28_seed_audit', 'strict_V2_real_hill_zero_import',
        'coarse8_fine16_same_endpoint', 'bitwise_coarse4_plus_4_restart',
        'post_assignment_rejection_full_rollback', 'one_sweep_counts_and_solved_epsilon',
        'strict_saved_history_budget_readbacks')}
    status = dict(version='H2kFinalStatusV1', capture_Asia_Shanghai=stamp,
        result='qualified_seed_and_passing_bounded_hill_startup', numerical_goal_closure=True,
        physical_application_validation_complete=False, gates=gates, state_sha256=cert['state_sha256'],
        cold_sha256=cert['immutable_base_sha256'], reference_sha256=cert['fixed_reference_sha256'],
        defect_sha256=cert['fixed_defect_sha256'], certificate_file=str(root/'fresh-epoch2/certificate.json'),
        certificate_sha256=sha(root/'fresh-epoch2/certificate.json'),
        actual_frontiers={k: dict(time=v['time'], steps=v['global_accepted_steps'], dt=v['dt'],
                                  checkpoint=v['checkpoint'], checkpoint_sha256=v['checkpoint_sha256'])
                          for k, v in analysis['runs'].items() if k in ('coarse8', 'fine16')},
        attempts=dict(primary_accepted_step_executions=24, split_replay_accepted_step_executions=8,
                      total_ordinary_accepted_executions=32, deliberate_rejected_attempts=2,
                      total_physical_attempts=34, unexpected_numerical_rejections=0,
                      zero_step_imports_not_physical_attempts=True),
        accepted_execution_equation_counts=native_counts,
        deliberate_rejected_attempt_equation_counts=rejected_counts,
        validation_only_repair=dict(original_ledger=str(root/'pilot/progress.json'),
            original_terminal_stop='accepted checkpoint missing', repaired_ledger=str(root/'pilot-requalified/progress.json'),
            reused_accepted_step_executions=32, additional_ordinary_steps=0,
            additional_deliberate_rejections=1, certified_native_sources_and_runtime_unchanged=True),
        aggregate_work=load(root/'aggregate.json')['totals'],
        elapsed_wall_seconds_at_collection=time.time()-origin['started_epoch'],
        H2k_jobs_left_running=[], long_continuation_launched=False,
        open_gates=['common flat spatial accuracy', 'independent CLS hill initialization/startup',
                    'temporal order and accumulated time error', 'developed flow',
                    'spatial/height/boundary sensitivity', 'held-out experimental validation'])
    atomic(a.output/'final-status.json', status)
    compact = dict(analysis)
    # Preserve all values, counts and source receipts; retain the complete timer
    # tree outside Git and only the two primary trajectories in the compact tree.
    compact['runs'] = {k: v for k, v in analysis['runs'].items() if k in ('coarse-h0', 'coarse8', 'fine16')}
    atomic(a.output/'startup.json', compact)
    atomic(a.output/'rollback.json', dict(version='H2kRollbackClosureV1',
        native_integrity=load(root/'rollback-native-integrity.json'),
        requalified_progress=progress, focused_regression=load(root/'rollback-validator-regression.json')))
    calibration = load(a.comparison)
    atomic(a.output/'continuing-comparison.json', dict(version='H2kContinuingComparisonV1',
        capture_Asia_Shanghai=stamp, calibration_file=str(a.comparison), calibration_sha256=sha(a.comparison),
        common={k: calibration[k] for k in ('u_tau', 'nu', 'H', 'z0_eff', 'kappa', 'U_H')},
        model_calibrations=calibration['models'],
        standard=dict(state_sha256=cert['state_sha256'], status=status['result'],
                      primary_accepted_steps=24, maximum_physical_time=max(coarse['time'], fine['time']),
                      seed_and_evolution_distinct=True),
        CLS=dict(status='archived blocked hill initializer; unchanged deployed model; zero hill steps in H2k',
                 evidence='docs/plans/cls-cao-tamura-budget-h2f.md'),
        matched_optimizer_performance_claim=False, open_gates=status['open_gates']))
    usages = {}
    for f in sorted(root.glob('*usage.txt')):
        text = f.read_text()
        rss = re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)', text)
        compact = re.search(r'^([0-9.]+) ([0-9.]+) ([0-9.]+) (\d+)$', text, re.M)
        usages[f.name] = dict(raw_file=str(f), sha256=sha(f),
                             maximum_RSS_KiB=int(rss.group(1)) if rss else int(compact.group(4)) if compact else None,
                             wall_user_system_RSS=[float(compact.group(i)) for i in range(1, 5)] if compact else None)
    atomic(a.output/'costs.json', dict(version='H2kMeasuredCostsV1',
        native_export_seconds={name: load(root/name/'native/receipt.json')['seconds'] for name in
                               ('epoch1-context', 'epoch2-context', 'seen-reserve', 'fresh-epoch2')},
        numerical_work=status['aggregate_work'], measured_process_usages=usages,
        native_pilot_usages={k: dict(seconds=v['seconds'], usage=v['usage']) for k, v in records.items()},
        timed_scope='GNU time includes the driver and its native children; native exclusive timer sums are separate in startup.json; no speedup claim from unlike fixtures'))
    commands = {}
    for f in sorted(root.glob('*usage.txt')):
        match = re.search(r'Command being timed: "(.*)"', f.read_text())
        if match:
            commands[f.name] = match.group(1)
    native = []
    for f in sorted(root.rglob('receipt.json')):
        if any(part in ('source', 'build') for part in f.relative_to(root).parts):
            continue
        d = load(f)
        if 'command' in d:
            native.append(dict(cwd=str(f.parent), command=d['command'], receipt=str(f), sha256=sha(f)))
    solver = Path(cert['runtime']['library_path']).parent/'IncNavierStokesSolver'
    require(sha(solver) == cert['runtime']['sha256'], 'certified solver path/hash differs')
    atomic(a.output/'commands.json', dict(version='H2kExecutedCommandsV1', timed_driver_commands=commands,
        native_probe_commands=native,
        production_invocations=[dict(xml=str(Path(v['report']).with_suffix('').with_suffix('.xml')),
                                     executable=str(solver), report=v['report']) for v in records.values()],
        environment=dict(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', OMP_THREAD_LIMIT='1'),
        rebuild_workers=4, native_ranks=1, note='Output directories are immutable; use unused roots for reproduction. Exact input packets, XML, logs and runtime closures are indexed.'))
    paths = []
    for f in sorted(root.rglob('*')):
        if not f.is_file() or any(part in ('source', 'build', 'tmp') for part in f.relative_to(root).parts):
            continue
        paths.append(f)
    entries = [dict(path=str(f.relative_to(root)), bytes=f.stat().st_size, sha256=sha(f)) for f in paths]
    atomic(a.output/'artifact-index.json', dict(version='H2kHeavyArtifactIndexV1',
        capture_Asia_Shanghai=stamp, root=str(root), storage_UUID='5AB8A669B8A642FD', artifacts=entries,
        included_bytes=sum(x['bytes'] for x in entries),
        source_manifest=str(root/'tested-source.json'), build_logs=['configure.log', 'build.log']))
    print(json.dumps(dict(passed=True, files=[f.name for f in sorted(a.output.iterdir())], artifacts=len(entries))))


if __name__ == '__main__':
    main()
