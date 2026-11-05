#!/usr/bin/env python3
"""Create the compact R11 decision ledger from immutable comparison and audit payloads."""
import argparse
from collections import Counter
import importlib.util
import json
from pathlib import Path
import platform
import sys

sys.dont_write_bytecode=True


def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    obj=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def keep(value,*names):
    return {k:value[k] for k in names if k in value}


def clean(value):
    return {k:v for k,v in value.items() if k not in
            ('schema','policy','action_policy','sequence','solve_number','alpha','time','kind')}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--epoch',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    assert platform.system()=='Linux'
    e=a.epoch.resolve()
    evidence=e/'source/docs/plans/cls-epawf-revised-evidence'
    h=module(evidence/'m8-r4/prepare.py','r11_review')
    C=h.load(e/'comparison-v2.json')
    assert C['native_history_verified']
    receipt=h.load(h.bound(C['receipt']))
    assert receipt['status']=='completed' and receipt['runtime_unchanged']
    manifest=h.load(h.bound(receipt['manifest']))
    events=C['events']
    r10=h.load(evidence/'m8-r10/decision-ledger.json')
    old_anchors=r10['accepted_model_anchors']
    affine=next(p for p in C['policies'] if p['name']=='L1')
    assert affine['direction']['sha256']==old_anchors[0]['direction']['sha256']
    anchors=[]
    replays=[]
    for r in events:
        if r['kind']!='r11_replay':continue
        if r['anchor']:
            index=int(r['name'].removeprefix('R10_anchor'))
            anchors.append(dict(index=index,Phi_1=old_anchors[index]['Phi'],
                Phi_fraction=r['Phi_fraction'],model=r['model'],branch_changes=r['ordered_branch_changes']))
        else:
            replays.append(clean(r))
    replay_outcome=next(r for r in events if r['kind']=='r11_replay_outcome')
    assert replay_outcome['passed']
    policy_audits={}
    for name in ('L1','Lf','N1','Nf'):
        path=e/'audits'/str(name+'.json')
        audit=h.load(path)
        result=keep(audit,'status','state','residual','history','fraction','seconds','new_native_trials')
        result['file']=h.pin(path)
        if audit['status'].startswith('reused'):
            result['prior_audit']=audit['audit']
        elif audit['status']=='independently_audited':
            result.update(keep(audit,'native_blocks','independent_blocks','difference_blocks','component_value_error',
                'trace','ground','common','accounting'))
            result['trace']={k:v for k,v in result['trace'].items() if k not in ('samples','projected_target_coefficients')}
            result['domain']=keep(audit['domain'],'admissible','sign_unresolved_observations','minimum_inner_k','minimum_outer_k',
                'minimum_outer_epsilon','minimum_b','minimum_lambda')
        policy_audits[name]=result
    model=[]
    for r in events:
        if r['kind'] in ('r11_model_anchor','r11_tangent','r11_action','r11_model_trial','r11_model_work','r11_seed_fallback'):
            record=clean(r)
            record['record_type']=r['kind']
            if 'details' in record:
                record['details']=keep(record['details'],'method','trust_radius','stationarity_l2',
                    'active_equality_error','complementarity','scaled_residual_objective','exchange_passes')
            model.append(record)
    cuts=[]
    for r in events:
        if r['kind']!='r6_cut_audit':continue
        failed=r['status']!='passed_original_cut_scan'
        relevant=[x for x in r['rows'] if x.get('incoming') or x.get('general_affine_row')]
        if not relevant:
            relevant=sorted(r['rows'],key=lambda x:abs(x['multiplier_scaled']),reverse=True)[:2]
        cuts.append(keep(r,'name','status','exchange','worst_native_violation','worst_extended_violation',
                         'maximum_scan_difference')|dict(
            active_family_counts=dict(Counter(x['kind'] for x in r['rows'])),
            candidate_direction=r['direction'] if failed else None,
            witnesses=[keep(x,'kind','point','geometry','block','axis','side','projected','general_affine_row',
                'incoming','rhs','native_gap','coefficient_gap','transformed_gap','multiplier_scaled') for x in relevant]))
    trials=[]
    for r in events:
        if r['kind']!='r8_trial':continue
        record=keep(r,'comparison','name','fraction','cache_reused','classification','actual','linear','matching_trace',
                    'state','residual','covariance','seconds','charged_seconds')
        record['wall']=r.get('wall')
        record['box']=r.get('box')
        if record['box'] and r['state']['sha256']!=C['selected']['state']['sha256']:
            record['box']={k:v for k,v in record['box'].items() if k!='witnesses'}
        record['scalar']=keep(r.get('scalar',{}),'admissible','minimum_dimensionless_margin','minimum_b','minimum_lambda')
        trials.append(record)
    log=(e/'regressions.log').read_text()
    assert '100% tests passed, 0 tests failed out of 3' in log and 'Ran 9 tests' in log and 'Ran 5 tests' in log
    h.save(a.output,dict(schema=1,milestone='M8-R11',completed=True,M8_closed=False,
        reviewed_baseline='a44639b99f9edd7f0bb4f8f2134a8065ee03b132',
        producing_revision=receipt['source_commit'],runtime=h.pin(e/'runtime-fingerprint.json'),
        executable=receipt['executable'],XML=receipt['session'],receipt=C['receipt'],comparison=h.pin(e/'comparison-v2.json'),
        collection=C['collection'],driver=h.pin(__file__),budget=h.pin(e/'local-bound/budget.json'),
        fractional_row_verification=h.pin(e/'fractional-row-checks.json'),
        failed_read_only_audit=dict(log=h.pin(e/'audits/N1-initial-error.log'),
            script=h.pin(e/'analysis/audit-selected-before-domain-cache.py'),
            error='KeyError phase in range/domain audit handoff',resolved=True,new_native_trials=0,elapsed_seconds=None),
        target=manifest['target'],bindings=keep(C['initial'],'state','history','residual','constraints')|
            dict(physical_J=C['matrix'],source_value=manifest['r10_source']['value'],source_action=manifest['r10_source']['action']),
        source_compatibility=h.load(e/'local-bound/budget.json')['physical_value_and_source_compatibility'],
        retrospective=dict(outcome=clean(replay_outcome),trials=replays,anchor_rank_comparison=anchors),
        policies=[clean(r) for r in C['policies']],model=model,original_row_checks=cuts,
        general_slope_scans=[clean(r) for r in events if r['kind']=='r10_general_row_scan'],
        native_trials=trials,native_costs=[clean(r) for r in events if r['kind']=='r11_native_cost'],
        selected_predictions=[clean(r) for r in events if r['kind']=='r11_selected_prediction'],
        promotions=[clean(r) for r in events if r['kind']=='r11_promotion'],decision=clean(C['decision']),
        audits=policy_audits,retained=dict(snapshot=C['selected']['snapshot'],state=C['selected']['state'],
            residual=C['selected']['residual'],audit_reuse=h.pin(e/'retained-audit-reuse.json'),
            native_report=C['native_report'],history_verified=True,physical_steps=0),
        continuation=dict(run=False,reason='No policy passed practical promotion'),
        work=dict(native_process_seconds=receipt['commands'][1]['seconds'],ranks=1,numerical_threads=receipt['numerical_threads'],
            compiler_workers=16,build=h.pin(e/'build-run/receipt.json')),
        regressions=dict(native_CTests_passed=3,predictor_binding_tests_passed=9,selected_candidate_tests_passed=5,
            compile_time_fractional_objective_assertions=True,log=h.pin(e/'regressions.log')),
        contract={'unknowns': 2165, 'reference_fraction': 0.125, 'root_gate': 2e-12, 'action_gate': 1e-07, 'original_row_gate': 5e-12, 'inner_stress': 'SaturatedR2', 'rho': 0.8, 'beta': 0.4, 'GammaScale': 1, 'GammaUniform': 0, 'destruction': 'Capped', 'smoothing': 'Linear', 'theta': 0.2, 'delta': 0.1, 'total_radius': 0.125, 'reserve': 0.05},
        outcomes={'A_early_slope_constraint': 'N1 improved the full-scale model objective within three tangents; its actual native gain failed promotion.', 'B_trial_scale_affine_objective': 'Unresolved. Lf exhausted the 96-cut exchange limit before a usable direction.', 'C_nonlinear_value_at_trial_scale': 'Unresolved. Nf exhausted its first tangent exchange; its chosen native result is the exact L1 seed.', 'D_root_and_history': 'No coupled root or accepted physical time step; original startup history verified.'},
        next_decision='Diagnose the archived Lf active-set subproblem and construct one independently verified Lf direction before more nonlinear-model or native trajectory work.',
        remaining_gates=['Curved-wall MMS', 'Active-cap and variable-height coupling', 'Full-element positivity', 'Fully discrete positivity', 'Refinement', 'Long-time robustness', 'Scalability', 'Physical floor/hill validation']))


if __name__=='__main__':
    main()
