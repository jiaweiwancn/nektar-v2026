#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Assess completed native blocks with the existing fixed-window policy.

No native solver, projection, reference reconstruction or history advancement is
performed. Native reports already contain the uncancelled alpha=0 residual.
Exact report bytes, checkpoint/state/reference and every GLM slot/time label are
checked before consuming their physical common-panel samples. Signed JSON zeros
must survive parsing. All heavier copies remain in the Linux experiment root.
"""
import argparse
from datetime import datetime, timezone, timedelta
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import shutil

from analyse_cls_epawf_height_comparison import differences
from analyse_cls_accepted_awf import normalized
from continue_cls_epawf_pimple import atomic, load_progress
from diagnose_cls_accepted_covariance import checkpoint
from run_cls_epawf_picard import THREADS, require, sha, state_hash
from cls_awf_native_metrics import interface_metrics


def native_report(segment, folder):
    folder.mkdir(parents=True, exist_ok=False)
    source = Path(segment['checkpoint']).parent / 'case.evolution.json.gz'
    if not source.exists():
        source = source.with_suffix('')
    payload = gzip.open(source, 'rb').read() if source.suffix == '.gz' else source.read_bytes()
    require(hashlib.sha256(payload).hexdigest() == segment['report_sha256'], 'native report bytes differ')
    # The native writer emits -0; ordinary integer parsing would lose its sign.
    report = json.loads(payload, parse_int=float)
    selected = report.get('accepted_awf_final') or report['height_comparison_final']
    original = Path(segment['checkpoint'])
    cp = folder / 'checkpoint'
    cp.mkdir()
    shutil.copy2(original / 'State.bin', cp / 'State.bin')
    require(sha(cp / 'State.bin') == sha(original / 'State.bin'), 'immutable checkpoint copy differs')
    native = checkpoint(cp)
    require(native['signature'] == report['signature'], 'physical signature differs')
    require(native['steps'] == report['steps'] == segment['end_step'] and native['time'] == report['time'], 'native time/step differs')
    require(state_hash(native['coupled']) == state_hash(report['state']) == selected['state_sha256'], 'selected state differs')
    require(state_hash(native['reference']) == state_hash(report['fixed_reference']) and
            state_hash(native['defect']) == state_hash(report['fixed_defect']), 'saved reference/defect differs')
    require([state_hash(v) for v in native['history']] == [state_hash(v) for v in report['native_history']] and
            state_hash(native['history_times']) == state_hash(report['native_history_times']), 'native GLM history differs')
    require(not report['failure'] and selected['state_unchanged'] and selected['governing_solves'] == 0, 'invalid native observation')
    require(not selected.get('invalid_production_diagnostic_points'), 'optional native reporting query was inadmissible')
    shutil.copy2(source, folder / source.name)
    receipt = dict(source_report=str(source), compressed_report_sha256=sha(source),
        uncompressed_report_sha256=segment['report_sha256'], checkpoint=str(cp), checkpoint_sha256=sha(cp/'State.bin'),
        state_sha256=selected['state_sha256'], native_history_sha256=[state_hash(v) for v in native['history']],
        native_history_times_sha256=state_hash(native['history_times']), physical_signature_sha256=hashlib.sha256(native['signature'].encode()).hexdigest(),
        audit_revision=report.get('audit_revision', 'LegacyInterfaceV1'), time=native['time'], step=native['steps'])
    atomic(folder/'receipt.json', receipt)
    return report, selected, receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign', action='append', required=True, help='LABEL=/absolute/trajectory')
    parser.add_argument('--policy', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    os.environ.update(THREADS)
    args.output.mkdir(parents=True, exist_ok=False)
    policy = load_progress(args.policy)
    require(policy['required_consecutive_windows'] == 3 and policy['primary_physical_window'] == .5, 'unexpected declared policy')
    out = dict(captured_Asia_Shanghai=datetime.now(timezone(timedelta(hours=8))).isoformat(),
        script_sha256=sha(Path(__file__)), policy=policy, policy_sha256=sha(args.policy),
        governing_solves=0, physical_commits=0, native_histories_advanced=False, cases={},
        selection='Latest completed endpoint with three preceding endpoints separated by exactly .5; selection uses availability only. Also retain the latest completed report.',
        scope='Common Ustar=.4,K0=.5105080131,E0=.0726358333,p0=.16,V=24 physical L2 panels; epsilon changes only on common z>=2. Full spatial epsilon residual uses the actual entire outer transport space. Negative classifications do not depend on omitting the lower outer epsilon region.',
        legacy_context='The older long-time report used a broader 1e-3 exploratory field screen. This comparison explicitly uses the fixed H=2 common normalization and the supplied policy for all six cases; it does not silently substitute H-dependent native residual scales or infer stationarity from endpoint samples.')
    for item in args.campaign:
        label, text = item.split('=', 1)
        require(Path(label).name == label and label not in ('.', '..'), 'invalid case label')
        campaign = Path(text)
        progress = load_progress(campaign/'progress.json')
        segments = {int(s['end_step']): s for s in progress['segments'] if not s['failure'] and s['exit_code'] == 0}
        stride = round(policy['primary_physical_window']/progress['dt'])
        eligible = [n for n in segments if n%stride==0 and all(n-i*stride in segments for i in range(4))]
        require(eligible, 'three complete physical windows are unavailable')
        end = max(eligible)
        steps = [end-i*stride for i in (3, 2, 1, 0)]
        latest = max(segments)
        chosen = sorted(set(steps+[latest]))
        samples = {}; diagnostics = {}
        for n in chosen:
            d, q, receipt = native_report(segments[n], args.output/label/f'step-{n:06d}')
            samples[n] = q
            audit = d['state_audits'][-1]
            attempt = d['pimple_attempts'][-1]
            bound_keys = ('mode','version','inventory_policy','success','activations','raw_ground','corrected_ground',
                'raw_minimum','corrected_minimum','mass_norm_raw_step','mass_norm_correction','mass_norm_accepted_step',
                'correction_to_step_ratio','integral_total','integrals_inner_outer_near_positive_negative',
                'active_set_iterations','mass_factor_applications','kkt_solves','seconds_total','raw_sha256','corrected_sha256')
            diagnostics[str(n)] = dict(time=q['time'], common_normalized_spatial_blocks=q['common_normalized_spatial_first4'],
                native_scaled_spatial_blocks=q['scaled_native_spatial_blocks'], dimensional_spatial_blocks=q['dimensional_mass_dual_spatial_first4'],
                matching_defect=q.get('matching_trace_defect', attempt['refreshed_trace_defect']),
                mean_wall_stress=q['mean_wall_stress'], relative_reference_traction_error=q['mean_wall_stress']/.16-1,
                bound={k:attempt['k_bound'][k] for k in bound_keys if k in attempt['k_bound']},
                strong_continuity=audit['strong_continuity'], interface=audit['interface'], named_interface_metrics=interface_metrics(audit['interface']),
                covariance=audit.get('accepted_awf_covariance', audit.get('r2_covariance')), receipt=receipt)
        records = {}
        with (campaign/'steps.jsonl').open() as stream:
            for line in stream:
                try: row = json.loads(line)
                except json.JSONDecodeError: break
                n = int(row['global_step'])
                if steps[0] <= n <= latest: records[n] = row
        windows = []
        for first, last in zip(steps, steps[1:]):
            left, right = samples[first], samples[last]
            duration = right['time']-left['time']
            require(abs(duration-.5) < 1e-10, 'physical window duration differs')
            changes = normalized(differences(right, left))
            require(all(n in records for n in range(first+1, last+1)), 'accepted correction ledger has a gap')
            bound = [records[n]['attempt']['k_bound'] for n in range(first+1, last+1)]
            signed = math.fsum(v['integral_total'] for v in bound)
            positive = math.fsum(v['integrals_inner_outer_near_positive_negative'][3] for v in bound)
            negative = math.fsum(v['integrals_inner_outer_near_positive_negative'][4] for v in bound)
            scale = policy['common_normalization']['K']*24/(2/.4)
            row = dict(start=left['time'], end=right['time'], dimensional_field_change=differences(right,left), normalized_field_change=changes,
                relative_traction_change=abs(right['mean_wall_stress']-left['mean_wall_stress'])/.16,
                signed_correction_integral=signed, positive_correction_integral=positive, negative_correction_integral=negative,
                absolute_correction_integral=positive-negative, normalized_signed_supply_rate=abs(signed)/duration/scale,
                normalized_absolute_supply_rate=(positive-negative)/duration/scale, bound_activations=sum(v['activations'] for v in bound))
            current = diagnostics[str(last)]
            row['fields_pass'] = all(changes[k] <= limit for k,limit in policy['normalized_window_L2_change_maximum'].items())
            row['residuals_pass'] = max(current['common_normalized_spatial_blocks']) <= policy['normalized_selected_spatial_residual_maximum_each_block']
            row['matching_pass'] = current['matching_defect']/policy['common_normalization']['epsilon'] <= policy['relative_matching_defect_maximum']
            row['traction_pass'] = row['relative_traction_change'] <= policy['relative_wall_stress_window_change_maximum']
            row['bound_supply_pass'] = row['normalized_absolute_supply_rate'] <= policy['bound_supply_per_time_relative_to_K0_volume_per_wall_time_maximum']
            row['common_screen_pass'] = all(row[k] for k in ('fields_pass','residuals_pass','matching_pass','traction_pass','bound_supply_pass'))
            windows.append(row)
        # Require all residual/field conditions to persist; one passing final
        # window does not establish steady convergence.
        pass_count = sum(w['common_screen_pass'] for w in windows)
        classification = 'not_steady_by_declared_screen'
        if pass_count == 3:classification = 'common_panel_steady_candidate_requires_full_outer_epsilon_confirmation'
        last = diagnostics[str(latest)]
        for n in chosen:
            require(records[n]['checkpoint_sha256'] == diagnostics[str(n)]['receipt']['checkpoint_sha256'], 'native checkpoint/accepted ledger differs')
        out['cases'][label] = dict(campaign=str(campaign),production_status=progress['status'],frontier_steps=progress['accepted_steps'],frontier_time=progress['accepted_time'],
            selected_steps=steps,latest_completed_step=latest,windows=windows,observations=diagnostics,
            common_screen_passing_windows=pass_count,classification=classification,
            common_residual_persistence_pass=all(max(diagnostics[str(n)]['common_normalized_spatial_blocks']) <= .001 for n in steps[1:]),
            latest_residual_pass=max(last['common_normalized_spatial_blocks']) <= .001)
        atomic(args.output/'summary.json',out)
    out['completed_Asia_Shanghai'] = datetime.now(timezone(timedelta(hours=8))).isoformat()
    atomic(args.output/'summary.json',out)
    print(args.output/'summary.json')


if __name__ == '__main__':
    main()
