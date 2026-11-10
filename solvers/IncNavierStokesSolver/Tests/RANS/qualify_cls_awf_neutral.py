#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Strict native reference/balance and matched physical-panel readbacks.

Selections are immutable native checkpoint/XML/report/ledger identities, not
reconstructed fields. Analytical loads are explicitly separate from native
represented residuals. Signed mass-dual projections preserve cancellation.
"""
import argparse
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import shutil
import struct

from analyse_cls_accepted_awf import compact, normalized, reference_errors
from analyse_cls_epawf_height_comparison import differences, compact_diagnostic
from continue_cls_epawf_pimple import atomic
from diagnose_cls_accepted_covariance import checkpoint, common_packet, run_probe, read_report
from run_cls_epawf_picard import THREADS, load, require, sha, state_hash


TERMS = ('inner_momentum_stress', 'outer_momentum_stress', 'conservative_advection', 'k_diffusion',
    'epsilon_diffusion', 'inner_k_production', 'inner_k_destruction', 'outer_k_production',
    'outer_k_destruction', 'epsilon_production', 'epsilon_destruction', 'epsilon_yap',
    'prescribed_scalar_source', 'roughness_force', 'prescribed_body', 'continuity', 'gauge')


def balances(d):
    import numpy as np
    from scipy.linalg import cho_factor, cho_solve
    from scipy.sparse import coo_matrix,diags
    from scipy.sparse.linalg import spsolve
    out = {}
    for metric in d['mass_metrics']:
        b = int(metric['block']); indices = np.array(metric['free_local_rows'], dtype=int)+int(metric['offset'])
        mass = np.array(metric['consistent_free_mass']); factor = cho_factor(mass)
        constant=None
        if metric.get('native_GL4_rows'):
            rr=[];cc=[];vv=[];weights=[]
            for i,q in enumerate(metric['native_GL4_rows']):
                for j,phi in zip(q['ids'],q['phi']):rr.append(i);cc.append(int(j));vv.append(phi)
                weights.append(q['weight'])
            evaluation=coo_matrix((vv,(rr,cc))).tocsr();weights=np.array(weights)
            full_mass=evaluation.T@diags(weights)@evaluation
            constant=spsolve(full_mass.tocsc(),evaluation.T@weights)
            require(np.max(abs(evaluation@constant-1))<1e-12,'native test-space constant reconstruction failed')
            local=np.array(metric['free_local_rows'],dtype=int)
            require(np.max(abs(full_mass[local][:,local].toarray()-mass))<1e-12,'native full/free mass consistency failed')
        rows = {}
        for kind in ('analytical_selected', 'analytical_Original', 'archived_reference_representation', 'retained_state'):
            values = d[kind]['terms']; loads = np.column_stack([np.array(values[n])[indices] for n in TERMS])
            total = np.array(values['total_constrained'])[indices]
            error = float(np.max(abs(loads.sum(axis=1)-total), initial=0.))
            require(error < 2e-10, 'free weak physical component sum does not close')
            direction = cho_solve(factor, total); norm = math.sqrt(max(0., float(total@direction)))
            signed = loads.T@direction/norm if norm > 1e-30 else np.zeros(len(TERMS))
            rows[kind] = dict(dimensional_mass_dual_norm=norm, component_sum_max_error=error,
                signed_component_projection_onto_total={n: float(v) for n, v in zip(TERMS, signed)},
                signed_projection_sum_error=abs(float(sum(signed))-norm))
            require(rows[kind]['signed_projection_sum_error'] < 1e-8, 'mass-dual cancellation accounting failed')
            if constant is not None:
                n=len(constant);offset=int(metric['offset']);free=np.array(metric['free_local_rows'],dtype=int)
                physical=np.array(values['spatial_sum'][offset:offset+n]);free_sum=float(constant[free]@physical[free]);total_sum=float(constant@physical)
                pieces={key:float(constant@np.array(values[key][offset:offset+n])) for key in TERMS}
                require(abs(sum(pieces.values())-total_sum)<1e-10,'constant-test weak balance does not close')
                rows[kind]['constant_test_balance']=dict(total_volume_residual=total_sum,free_row_contraction=free_sum,
                    prescribed_row_reaction=total_sum-free_sum,inferred_outward_boundary_supply=free_sum-total_sum,
                    signed_component_integrals=pieces,closure_error=abs(sum(pieces.values())-total_sum),
                    convention='alpha=0 native weak volume load; outward prescribed-boundary flux reaction is minus the prescribed-row contraction. This is a discrete reaction, not an exactly conservative local CG physical-flux estimate.')
            if kind in ('retained_state', 'archived_reference_representation') and metric.get('native_GL4_rows'):
                terms_rate = cho_solve(factor, loads)
                index = {int(j): i for i, j in enumerate(metric['free_local_rows'])}
                elements = {}
                for q in metric['native_GL4_rows']:
                    value = 0.; components = np.zeros(len(TERMS))
                    for j, phi in zip(q['ids'], q['phi']):
                        if int(j) in index:
                            i = index[int(j)]; value += phi*direction[i]; components += phi*terms_rate[i]
                    e = elements.setdefault(str(int(q['geometry'])), dict(residual_rate_L2_squared=0., signed_component_dual_work=np.zeros(len(TERMS))))
                    e['residual_rate_L2_squared'] += q['weight']*value*value
                    e['signed_component_dual_work'] += q['weight']*components*value
                quadrature = math.sqrt(max(0., math.fsum(e['residual_rate_L2_squared'] for e in elements.values())))
                require(abs(quadrature-norm)<1e-10*max(1.,norm), 'native element localization does not reproduce physical mass-dual norm')
                for e in elements.values():
                    e['signed_component_dual_work'] = {n: float(v) for n,v in zip(TERMS,e['signed_component_dual_work'])}
                rows[kind]['element_residual_localization'] = elements
                rows[kind]['mass_quadrature_norm_error'] = abs(quadrature-norm)
        out[str(b)] = rows
    out['semantics'] = 'Signed <R_component,M_ff^-1 R_total>/||R_total|| dual projections add to total norm. Component norms are not additive fractions of a traction error. Full loads and constrained reactions retained separately.'
    return out


def observe(probe, selection, folder):
    folder.mkdir(parents=True, exist_ok=False); original = Path(selection['checkpoint']); cp = folder/'checkpoint'
    shutil.copytree(original, cp); require(sha(cp/'State.bin') == selection['checkpoint_sha256'] == sha(original/'State.bin'), 'native selection checkpoint differs')
    saved = checkpoint(cp); require(saved['time'] == selection['native_time'] and saved['steps'] == selection['native_step'], 'selected native clock/count differs')
    H = selection.get('H', .5)
    packet = dict(common_packet(cp), operation='selected_original_operator', state=saved['coupled'], history=saved['coupled'],
        alpha=0., time=saved['time'], inner_stress_model='Linear', complete_outer_epsilon=True,
        reporting_breakpoints=[.0025, .005, .02, .04], neutral_reference_balance=selection.get('balance', True),
        witnesses=[[-3.5, 0, H], [-2.625, 0, H], [-1., 0, H], [1.6, 0, H]])
    result, receipt = run_probe(probe, Path(selection['xml']), packet, folder/'readback')
    require(result['governing_solves'] == result['physical_commits'] == 0 and not result['native_history_advanced'] and result['initialization_work'] == [0, 0], 'readback changed state/history')
    d = result['diagnostics']; physical = d['physical']; require(physical['epsilon_reporting_lower'] == H, 'incomplete transported epsilon region')
    require(physical['state_sha256'] == state_hash(saved['coupled']) and physical['state_unchanged'], 'retained state changed')
    summary = dict(selection=selection, receipt=receipt, state_sha256=state_hash(saved['coupled']),
        native_history_sha256=[state_hash(v) for v in saved['history']], native_history_times_sha256=state_hash(saved['history_times']),
        reference_sha256=state_hash(saved['reference']), defect_sha256=state_hash(saved['defect']),
        state_equals_fixed_reference=state_hash(saved['coupled']) == state_hash(saved['reference']),
        state_reference_numeric_differences=sum(x!=y for x,y in zip(saved['coupled'],saved['reference'])),
        state_reference_bit_differences=sum(struct.pack('<d',x)!=struct.pack('<d',y) for x,y in zip(saved['coupled'],saved['reference'])),
        physical=compact_diagnostic(compact(physical)), complete_audit=result['complete_audits']['selected_state'],
        selected_Original_blocks=[{k: v for k, v in b.items() if k not in ('free_rows', 'selected_load', 'same_state_Original_load', 'difference_load')} for b in d['operators']['blocks']])
    if 'neutral_reference_balance' in d:
        balance = d['neutral_reference_balance']
        summary['balance'] = balances(balance)
        summary['analytical_stress_difference'] = balance['analytical_maximum_stress_difference']
        summary['analytical_production_difference'] = balance['analytical_maximum_production_difference']
        require(summary['analytical_stress_difference']<1e-12 and summary['analytical_production_difference']<1e-12,
            'selected law does not preserve the analytical neutral stress/production')
        summary['native_outer_constitutive_difference_locations']=balance.get('native_outer_constitutive_difference_locations')
        summary['reference_operator_blocks'] = [{k: v for k, v in b.items() if k not in ('free_rows', 'selected_load', 'same_state_Original_load', 'difference_load')} for b in balance['reference_operators']['blocks']]
        summary['element_values_columns'] = balance['element_values_columns']; summary['element_values'] = balance['element_values']
        summary['native_ground_traction'] = balance.get('native_ground_traction')
    if selection.get('source_report'):
        native = read_report(Path(selection['source_report']))
        summary['source_report_sha256'] = sha(Path(selection['source_report']))
        matching = [t for t in native['pimple_attempts'] if t['time'] == saved['time']]
        if matching:
            t = matching[-1]
            summary['step'] = {k: t.get(k) for k in ('outer_sweeps', 'momentum_predictor_solves', 'pressure_system_solves', 'epsilon_transport_solves', 'k_transport_solves', 'bound_calls', 'additional_governing_equation_solves', 'frozen_momentum_defect', 'momentum_splitting_identity_defect', 'refreshed_trace_defect', 'bound_activations')}
            summary['step']['bound_integral'] = t['k_bound']['integral_total']
            summary['step']['bound_positive_negative'] = t['k_bound']['integrals_inner_outer_near_positive_negative'][3:5]
    atomic(folder/'summary.json', summary)
    return summary, physical


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('probe', 'selections', 'output'): p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--reuse',type=Path,help='reuse immutable identical-probe/checkpoint readbacks for new algebra reductions')
    a = p.parse_args(); os.environ.update(THREADS); a.output.mkdir(parents=True, exist_ok=False)
    selected = load(a.selections); out = dict(status='running', probe_sha256=sha(a.probe), script_sha256=sha(Path(__file__)),
        selections_sha256=sha(a.selections), started_UTC=datetime.now(timezone.utc).isoformat(), cases={}, comparisons={},
        common_domains='whole-domain U/k/pressure, full .5<=z<=4 epsilon, fixed union-of-mesh physical panels, one pressure gauge',
        governing_solves=0, physical_commits=0, native_history_advanced=False)
    previous=load(a.reuse/'summary.json') if a.reuse else None
    if previous: require(previous['probe_sha256']==sha(a.probe),'reused native probe identity differs')
    data = {}
    for label, selection in selected['states'].items():
        if previous:
            row=previous['cases'][label];require(row['selection']==selection,'reused native selection differs')
            file=Path(row['receipt']['output']);require(sha(file)==row['receipt']['output_sha256'],'reused native output differs')
            cp=file.parents[1]/'checkpoint';saved=checkpoint(cp)
            require(sha(cp/'State.bin')==selection['checkpoint_sha256'] and state_hash(saved['coupled'])==row['state_sha256'],'reused immutable checkpoint differs')
            raw=load(file);physical=raw['diagnostics']['physical']
            if 'neutral_reference_balance' in raw['diagnostics']:row['balance']=balances(raw['diagnostics']['neutral_reference_balance'])
            row['readback_reuse']=dict(summary=str(a.reuse/'summary.json'),summary_sha256=sha(a.reuse/'summary.json'))
        else: row, physical = observe(a.probe, selection, a.output/label)
        data[label] = physical
        row['dimensional_reference_errors'] = reference_errors(physical)
        row['common_normalized_reference_errors'] = normalized(row['dimensional_reference_errors'])
        out['cases'][label] = row; atomic(a.output/'summary.json', out)
    for label, (left, right) in selected.get('comparisons', {}).items():
        delta = differences(data[left], data[right])
        out['comparisons'][label] = dict(left=left, right=right, dimensional_field_difference=delta, normalized_field_difference=normalized(delta),
            wall_traction_difference=data[left]['mean_wall_stress']-data[right]['mean_wall_stress'])
    out.update(status='passed', completed_UTC=datetime.now(timezone.utc).isoformat()); atomic(a.output/'summary.json', out)


if __name__ == '__main__': main()
