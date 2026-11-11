#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Reduce immutable status captures and matched readbacks to compact evidence."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import statistics

from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import load, require, sha


COUNTS = ('outer_sweeps', 'momentum_predictor_solves', 'pressure_system_solves',
          'epsilon_transport_solves', 'k_transport_solves', 'bound_calls',
          'additional_governing_equation_solves', 'additional_correction_or_restoration_solves',
          'mass_inverse_applications', 'bound_activations', 'bound_active_set_iterations',
          'bound_mass_factor_applications', 'bound_kkt_solves')


def divergence(row):
    out = {}
    for region in (0, 1):
        cells = [e for e in row['element_values'] if e['values'][0] == region]
        value = math.sqrt(math.fsum(e['values'][2] for e in cells))
        require(abs(value-row['complete_audit']['strong_continuity'][region]) < 1e-12,
                'native element divergence does not reproduce the region integral')
        out['outer' if region else 'inner'] = dict(physical_L2=value,
            physical_measure=math.fsum(e['values'][1] for e in cells),
            largest_divergence_squared_cells=sorted(cells, key=lambda e: (-e['values'][2], e['geometry']))[:6],
            columns=row['element_values_columns'])
    return out


def work(capture):
    """Count each accepted ledger row once, through the immutable frontier."""
    progress = load(Path(capture['immutable_progress'])); end = int(capture['accepted_steps'])
    start = int(progress.get('prefix_steps', 0)); records = {}
    with (Path(capture['source_progress']).parent/'steps.jsonl').open() as stream:
        for line in stream:
            row = json.loads(line); n = int(row['global_step'])
            if n > end: continue
            require(n not in records, 'duplicate accepted ledger record'); records[n] = row
    require(sorted(records) == list(range(start+1, end+1)), 'accepted ledger is not contiguous')
    if records: require(records[end]['checkpoint_sha256'] == capture['checkpoint_sha256'], 'ledger/immutable checkpoint differs')
    counts = Counter(); signed = []; positive = []; negative = []; seconds = []; maxima = {}; minima = {}
    for row in records.values():
        t = row['attempt']; require(t['scheduled_complete'] and not t['failure'], 'invalid accepted attempt')
        require(all(t[k] == 1 for k in COUNTS[:6]) and t['additional_governing_equation_solves'] == 0, 'production schedule differs')
        for k in COUNTS: counts[k] += t.get(k, 0)
        b = t['k_bound']; signed.append(b['integral_total']); positive.append(b['integrals_inner_outer_near_positive_negative'][3]); negative.append(b['integrals_inner_outer_near_positive_negative'][4])
        seconds.append(t['seconds'])
        for k in ('predictor_linear_residual', 'pressure_linear_residual', 'epsilon_frozen_linear_residual',
                  'k_frozen_linear_residual', 'boundary_defect', 'momentum_splitting_identity_defect',
                  'native_full_G_action_defect', 'frozen_momentum_defect', 'refreshed_trace_defect'):
            maxima[k] = max(maxima.get(k, 0.), abs(t[k]))
        maxima['corrected_continuity_gauge'] = max(maxima.get('corrected_continuity_gauge', 0.), *map(abs, t['corrected_continuity_gauge']))
        for k in ('minimum_inner_k', 'minimum_outer_k', 'minimum_outer_epsilon', 'minimum_sampled_K'):
            value = t['k_domain'][k]; minima[k] = min(minima.get(k, value), value)
        value = t['k_ground']['minimum']; minima['continuous_ground_k'] = min(minima.get('continuous_ground_k', value), value)
    callbacks = 0; reported = 0; receipts = []; initialization = [0., 0.]; startup_checks = 0
    for segment in progress['segments']:
        if segment['end_step'] > end or segment['failure'] or segment['exit_code']: continue
        file = Path(segment['checkpoint']).parent/'case.evolution.json.gz'
        payload = gzip.open(file, 'rb').read()
        require(hashlib.sha256(payload).hexdigest() == segment['report_sha256'], 'completed native report differs')
        d = json.loads(payload, parse_int=float); steps = len(d['pimple_attempts'])
        require(d['implicit_callback_count'] == steps == segment['end_step']-segment['start_step'], 'callback/step count differs')
        callbacks += int(d['implicit_callback_count']); reported += steps
        initialization = [a+b for a,b in zip(initialization, d['initialization_lift_calls_iterations'])]
        require(d['initial_k_bound_check']['raw_sha256'] == d['initial_k_bound_check']['corrected_sha256'],
                'startup correction changed the restored state')
        startup_checks += 1
        receipts.append(dict(path=str(file), compressed_sha256=sha(file), native_payload_sha256=segment['report_sha256'], seconds=segment['seconds']))
    return dict(prefix_steps_excluded=start, accepted_new_steps=end-start, actual_counts=dict(counts),
        accepted_attempt_maxima=maxima, accepted_inventory_minima=minima,
        completed_block_callbacks=callbacks, completed_block_steps=reported, partial_block_steps=end-start-reported,
        completed_block_initialization_lift_calls_iterations=initialization,
        completed_block_startup_bound_checks=startup_checks,
        correction_integral=math.fsum(signed), correction_positive=math.fsum(positive), correction_negative=math.fsum(negative),
        correction_absolute=math.fsum(positive)-math.fsum(negative),
        implicit_seconds=math.fsum(seconds), median_implicit_seconds=statistics.median(seconds) if seconds else None,
        completed_block_wall_seconds=math.fsum(s['seconds'] for s in progress['segments'] if s['end_step']<=end and not s['failure'] and not s['exit_code']),
        reports=receipts)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'output'): p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args(); n = a.root
    status = load(n/'current-status.json'); matched = load(n/'matched-final-readbacks/summary.json')
    out = dict(captured_UTC=datetime.now(timezone.utc).isoformat(), experiment=str(n),
        collector_sha256=sha(Path(__file__)), source_factory_revision='b6c65b49cfc192753f053c48bab5ded650415ede',
        root_status_sha256=sha(n/'current-status.json'), readback_summary_sha256=sha(n/'matched-final-readbacks/summary.json'),
        cases={}, observations={}, comparisons=matched['comparisons'], readback_status=matched['status'],
        readback_failure=matched.get('failure'), production_h0_bitwise_verified=matched.get('production_h0_bitwise_verified'),
        baseline_summaries=matched['baseline_summaries'], governing_solves_added=0, physical_commits_added=0,
        note='Legacy refreshed_original names mean the selected uncancelled operator. Common normalization and full .5-to-4 epsilon are unchanged. Partial-block callback totals are explicitly separate from committed per-step solve counts.')
    baseline = {}
    for receipt in matched['baseline_summaries']:
        file = Path(receipt['path']); require(sha(file) == receipt['sha256'], 'baseline readback summary differs')
        for label, row in load(file)['cases'].items():
            if label.startswith('refined-'): baseline[int(row['selection']['native_step'])] = row
    for label, capture in status['cases'].items(): out['cases'][label] = dict(capture=capture, work=work(capture))
    for step, row in matched['observations'].items():
        q = row['physical']
        kept = ('state_sha256', 'native_history_sha256', 'native_history_times_sha256', 'reference_sha256', 'defect_sha256',
                'selection', 'receipt', 'complete_audit', 'dimensional_reference_errors', 'common_normalized_reference_errors',
                'drift_from_actual_h0', 'step', 'native_ground_traction')
        entry = {k: row[k] for k in kept if k in row}
        entry['physical'] = {k: q[k] for k in ('time', 'H', 'h', 'z_K', 'z_m', 'K_sampling', 'reference',
            'native_degrees_of_freedom', 'dimensional_mass_dual_spatial_first4', 'common_normalized_spatial_first4',
            'scaled_native_spatial_blocks', 'representation_L2_velocity_gradient_production_referenceProduction',
            'representation_L2_k_gradient', 'fixed_bands', 'mean_wall_stress', 'matching_trace_defect',
            'common_reporting_queries_valid', 'epsilon_reporting_lower')}
        entry['signed_mass_dual_components'] = {b: {k: v for k, v in row['balance'][b]['retained_state'].items()
            if k not in ('element_residual_localization',)} for b in map(str, range(4))}
        entry['native_element_divergence'] = divergence(row)
        out['comparisons'][step]['P1_native_element_divergence'] = divergence(baseline[int(step)])
        out['observations'][step] = entry
    for name in ('verification-matched/summary.json', 'verification/summary.json', 'matched-closure-equivalence.json',
                 'initial-operator-equivalence.json', 'preliminary-prefix-equivalence.json', 'pilot-matched-launch.json', 'analysis-matched-launch.json'):
        file = n/name; out.setdefault('receipts', {})[name] = dict(path=str(file), sha256=sha(file), data=load(file))
    atomic(a.output, out); print(a.output)


if __name__ == '__main__': main()
