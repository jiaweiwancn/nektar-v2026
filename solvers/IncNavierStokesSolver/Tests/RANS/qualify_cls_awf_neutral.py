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
from cls_awf_native_metrics import flat_partition, interface_metrics


TERMS = ('inner_momentum_stress', 'outer_momentum_stress', 'conservative_advection', 'k_diffusion',
    'epsilon_diffusion', 'inner_k_production', 'inner_k_destruction', 'outer_k_production',
    'outer_k_destruction', 'epsilon_production', 'epsilon_destruction', 'epsilon_yap',
    'prescribed_scalar_source', 'roughness_force', 'prescribed_body', 'continuity', 'gauge')


def step_storage_reaction(summary,raw):
    """Actual callback-history storage; never reconstruct history from a state."""
    d=raw['diagnostics'];balance=d.get('neutral_reference_balance')
    source=summary['selection'].get('source_report')
    if not balance or not source:return None
    native=read_report(Path(source));attempts=[v for v in native['pimple_attempts'] if v['time']==summary['selection']['native_time']]
    if not attempts:return None # h0 has no physical time derivative supplied
    import numpy as np
    t=attempts[-1];file=Path(t['k_snapshot']);packet=load(file)
    require(state_hash(packet['corrected'])==summary['state_sha256'] and packet['time']==summary['selection']['native_time'],
        'storage diagnostic packet/state/time differs')
    metric=balance['mass_metrics'][0];n=int(summary['physical']['native_degrees_of_freedom'][0]);lift=np.zeros(n)
    elements={int(e['geometry']):e for e in balance['element_values']};seen=set()
    xmin=min(e['lower'][0] for e in elements.values());L=max(e['upper'][0] for e in elements.values())-xmin
    # Exact modified-A coefficients of Q2(x)*P1(z), nonzero only on ground.
    for row in metric['native_GL4_rows']:
        g=int(row['geometry']);e=elements[g]
        if g in seen or e['lower'][2]!=0:continue
        seen.add(g);x0,x1=e['lower'][0],e['upper'][0];a=(x0-xmin)/L;b=(x1-xmin)/L
        for j,value in ((0,4*a*(1-a)),(1,4*b*(1-b)),(2,4*((x1-x0)/L)**2)):
            lift[int(row['ids'][j])]=value
    spatial=np.array(balance['retained_state']['terms']['spatial_sum'][:n]);reaction=balance['ground_reaction']['retained_state']
    require(abs(float(lift@spatial)-reaction['full_spatial_lifting_contraction'])<1e-12,'native ground lifting reconstruction differs')
    delta=np.array(packet['corrected'][:n])-np.array(packet['history'][:n]);storage=np.zeros(n);alpha=packet['alpha']
    for row in metric['native_GL4_rows']:
        ids=np.array(row['ids'],dtype=int);phi=np.array(row['phi']);storage[ids]+=alpha*row['weight']*phi*float(phi@delta[ids])
    moment=float(lift@storage);measure=reaction['weighted_ground_measure']
    from scipy.linalg import cho_factor,cho_solve
    free=np.array(metric['free_local_rows'],dtype=int);stage=spatial[free]+storage[free]
    factor=cho_factor(np.array(metric['consistent_free_mass']));stage_norm=math.sqrt(max(0.,float(stage@cho_solve(factor,stage))))
    return dict(packet=str(file),packet_sha256=sha(file),callback_history_sha256=state_hash(packet['history']),alpha=alpha,
        storage_lifting_moment=moment,spatial_lifting_moment=float(lift@spatial),
        refreshed_free_momentum_stage_residual_mass_dual=stage_norm,
        inferred_mean_shear_including_last_callback_storage=-(reaction['prescribed_row_reaction']+moment)/measure,
        scope='actual last-callback alpha*M*(U-history), physical GL4 consistent mass and exact ground lifting; refreshed spatial operator still contains one-sweep splitting/lag defects, so this is a discrete reaction diagnostic, not exact physical wall stress')


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
        reporting_breakpoints=sorted(set(flat_partition(Path(selection['xml']))[1]+selection.get('reporting_breakpoints',[]))), neutral_reference_balance=selection.get('balance', True),
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
    summary['named_interface_metrics'] = interface_metrics(summary['complete_audit']['interface'])
    constitutive=d['selected_state']
    summary['same_state_covariance_minima']={}
    for name,law in (('minimum_inner','selected'),('minimum_outer','selected'),
        ('minimum_same_state_Original_outer','same_state_Original')):
        witness=constitutive[name]
        summary['same_state_covariance_minima'][name]=dict(coordinate=witness['coordinate'],geometry=witness['geometry'],
            element=witness['element'],eigenvalue=witness[law]['eigenvalues'][0],
            velocity_gradient=witness['velocity_gradient'],epsilon_semantics=witness['epsilon_semantics'])
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
        summary['ground_reaction'] = balance.get('ground_reaction')
        if 'transport_decomposition' in balance:
            import numpy as np
            from scipy.linalg import cho_factor,cho_solve
            transport={}
            for kind in ('archived_reference_representation','retained_state'):
                transport[kind]=[]
                for item in balance['transport_decomposition'][kind]:
                    b=int(item['block']);metric=balance['mass_metrics'][b];free=np.array(metric['free_local_rows'],dtype=int)
                    factor=cho_factor(np.array(metric['consistent_free_mass']))
                    keys=('native_conservative_load','advective_load','dilatation_load','outward_boundary_load','physical_isotropic_production_load')
                    loads=np.column_stack([np.array(item[k])[free] for k in keys]);rates=cho_solve(factor,loads)
                    norms=np.sqrt(np.maximum(0,np.sum(loads*rates,axis=0)))
                    direction=rates[:,0];norm=norms[0]
                    signed=loads.T@direction/norm if norm>1e-30 else np.zeros(len(keys))
                    transport[kind].append(dict(block=b,dimensional_mass_dual_norms=dict(zip(keys,map(float,norms))),
                        signed_projection_onto_conservative_advection=dict(zip(keys,map(float,signed))),
                        integrated_q_div_u=item['integrated_q_div_u'],integrated_physical_isotropic_production=item['integrated_physical_isotropic_production'],
                        load_identity_max_error=item['load_identity_max_error']))
            summary['transport_decomposition']=transport
    if selection.get('source_report'):
        native = read_report(Path(selection['source_report']))
        summary['source_report_sha256'] = sha(Path(selection['source_report']))
        matching = [t for t in native['pimple_attempts'] if t['time'] == saved['time']]
        if matching:
            t = matching[-1]
            summary['step'] = {k: t.get(k) for k in ('outer_sweeps', 'momentum_predictor_solves', 'pressure_system_solves', 'epsilon_transport_solves', 'k_transport_solves', 'bound_calls', 'additional_governing_equation_solves', 'frozen_momentum_defect', 'momentum_splitting_identity_defect', 'refreshed_trace_defect', 'bound_activations')}
            summary['step']['bound_integral'] = t['k_bound']['integral_total']
            summary['step']['bound_positive_negative'] = t['k_bound']['integrals_inner_outer_near_positive_negative'][3:5]
    summary['last_step_ground_storage']=step_storage_reaction(summary,result)
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
            row['last_step_ground_storage']=step_storage_reaction(row,raw)
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
