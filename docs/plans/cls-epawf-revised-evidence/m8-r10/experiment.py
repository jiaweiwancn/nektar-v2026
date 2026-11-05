#!/usr/bin/env python3
"""Prepare and verify the frozen-base R10 comparison on authoritative Linux."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
STATE = 'eb063e15f16b2610f1c5a832280d87894f567e8dab321560202e6dd9af95194c'
MATRIX = '2194965dee9c5960ba4e1c1604f76e8e85372cabef2d5ced7c8a394f593c8dc6'
H0 = 'b8ba7917424057a9f6f16a026d82d53ede43985871537facea170ef6e1841b71'


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def prepare(e, out, h, evidence):
    r9 = h.load(evidence/'m8-r9/decision-ledger.json')
    original = r9['binding']['initial']
    matrix = r9['binding']['current_matrix']
    assert original['state']['sha256'] == STATE and matrix['sha256'] == MATRIX
    assert r9['target']['original_h0']['sha256'] == H0
    assert r9['terminal']['snapshot']['state']['sha256'] != STATE
    for pin in [*original.values(), matrix]:
        if isinstance(pin, dict) and 'path' in pin:
            h.bound(pin)
    comp = h.load(h.bound(r9['binding']['comparison']))
    previous = h.load(h.bound(h.load(h.bound(r9['binding']['receipt']))['manifest']))
    old_xml = h.bound(r9['xml'])
    previous_runtime = h.load(h.bound(r9['runtime']))
    old_source = Path(previous_runtime['source_root'])
    current_source = e/'source'
    relative = Path('solvers/IncNavierStokesSolver/EquationSystems')
    unchanged = {}
    for filename in ('RANSWallEvolutionNonpressure.cpp', 'RANSWallEvolutionRobust.cpp',
                     'RANSWallEvolutionResolution.cpp', 'RANSWallEvolutionTransition.cpp'):
        old = h.pin(old_source/relative/filename)
        current = h.pin(current_source/relative/filename)
        assert old['sha256'] == current['sha256']
        unchanged[filename] = current['sha256']

    def function_bytes(source, name):
        text=(source/relative/'RANSWallEvolutionPatch.cpp').read_text()
        start=text.index('RANSWallEvolutionPatch::'+name+'(')
        begin=text.index('{', start)
        depth=1
        end=begin+1
        while depth:
            depth += (text[end]=='{')-(text[end]=='}')
            end += 1
        return text[start:end].encode()

    for name in ('Residual', 'Norms', 'DestructionValueR9'):
        old = function_bytes(old_source, name)
        current = function_bytes(current_source, name)
        assert old == current
        unchanged[name] = hashlib.sha256(current).hexdigest()
    runtime_path = e/'runtime-fingerprint.json'
    runtime = h.load(runtime_path)
    exe = h.pin(e/'build/IncNavierStokesSolver')
    assert runtime['success'] and runtime['runtime_artifacts']['build/IncNavierStokesSolver'] == exe['sha256']
    out.mkdir()
    tree = ET.parse(old_xml)
    info = tree.find('.//SOLVERINFO')
    selectors = {r.get('PROPERTY'): r for r in info}
    assert selectors['RANSWallDestructionSmoothing'].get('VALUE') == 'Linear'
    assert selectors['RANSWallR9Experiment'].get('VALUE') == 'On'
    info.remove(selectors['RANSWallR9Experiment'])
    ET.SubElement(info, 'I', PROPERTY='RANSWallR10Experiment', VALUE='On')
    selectors['RANSWallR3PredictorFile'].set('VALUE', str(out/'predictor.json'))
    selectors['RANSWallR3ProgressFile'].set('VALUE', str(out/'case.r3-progress.jsonl'))
    xml = out/'case.xml'
    tree.write(xml, encoding='unicode')
    assert h.predictor.protected_xml(xml) == previous['protected_physical_xml_sha256']
    source = r9['source_component']['component']
    ascent = r9['refinements'][1]['direction']
    for pin in (source['D0'], source['JD0'], ascent):
        h.bound(pin)
    h.save(out/'predictor.json', dict(previous, state=original['state'],
        history=r9['target']['original_h0'], producing_source_commit=runtime['source_commit'],
        current_executable=exe, current_runtime=h.pin(runtime_path), current_target_xml=h.pin(xml),
        origin=dict(alpha=400, time=.0025, converged=False, inner_stress_model='SaturatedR2',
            iteration=10, snapshot=original['snapshot'], residual=original['residual'], constraints=original['constraints']),
        frozen_action=dict(matrix=matrix, state_sha256=STATE, history_sha256=H0,
            residual_sha256=original['residual']['sha256'],
            assembly_json=json.dumps(comp['matrix_details']['assembly'])),
        r10_source=dict(value=source['D0'], action=source['JD0'], ascent=ascent)))
    h.save(out/'budget.json', dict(tangent_qp_limit=4,conditional_slope_retry_limit=1,
        source_Armijo_fractions_per_tangent=10,new_source_value_limit=50,native_trials_per_ray=9,
        comparison_cooperative_seconds=600,stop_source_work_at_seconds=540,
        comparison_reserve_seconds=35,comparison_selector_full_endpoint_tested=False,
        physical_value_and_source_compatibility=unchanged,
        native_external_supervisor_seconds=650,source='Linear smoothing',base=original['state'],history=r9['target']['original_h0']))


def summary(e, out, h, evidence):
    receipt = h.load(e/'local-run/receipt.json')
    manifest = h.load(h.bound(receipt['manifest']))
    collection = h.load(h.bound(receipt['durable_collection']['file']))
    events = module(evidence/'m8-r7/experiment.py', 'r7_events').events(collection)
    base, *updates = collection['completed_iterates']
    assert base['state']['sha256'] == STATE and base['history']['sha256'] == H0
    for key in ('residual', 'constraints'):
        assert base[key]['sha256'] == manifest['origin'][key]['sha256']
    matrix, = [j for j in collection['jacobians'] if j['status'] == 'complete']
    assert matrix['iteration'] == 0 and matrix['event']['details']['loaded_bound_matrix']
    assert matrix['event']['matrix']['sha256'] == MATRIX
    assert collection['terminal']['terminal'] == 'r10_fixed_comparison_completed'
    rows = [x for x in events if x['kind'].startswith('r10_')]
    def event(kind):
        return next(x for x in events if x['kind'] == kind)
    assert event('r10_ascent_regression')['rejected_before_native_trials']
    assert event('r10_scaled_anchor_regression')['rejected_as_restored_anchor']
    initial_feasible = [r for r in rows if r['kind']=='r10_anchor_feasibility' and r['name']=='initial']
    anchors = [r for r in rows if r['kind']=='r10_anchor']
    trials = [r for r in rows if r['kind']=='r10_model_trial']
    if anchors:
        assert initial_feasible[0]['passed']
        for a,b in zip(anchors, anchors[1:]):
            assert b['Phi'] < a['Phi'] and b['native_slope'] <= -2e-4*max(base['scaled_block_norms'])
        for trial in trials:
            if trial['accepted']:
                assert trial['reason']=='accepted' and trial['model_armijo']
            else:
                assert trial['reason']!='accepted'
        outcome=event('r10_model_outcome')
        assert outcome['solves']<=4 and outcome['source_values']<=50
        assert outcome['Phi_best']<=anchors[0]['Phi']
        assert all(sum(t['index']==i for t in trials)<=10 for i in range(4))
    actual = [r for r in events if r['kind']=='r8_trial']
    assert len(actual)<=18
    for label in ('r10_linear','r10_source'):
        assert sum(r['comparison']==label for r in actual)<=9
    assert all(not r['history_advanced'] and r['history']['sha256']==H0 for r in collection['completed_iterates'])
    assert receipt['commands'][1]['exit_code']==0
    assert len(updates)<=1
    if updates:
        selected, = [r for r in actual if r['state']['sha256']==updates[0]['state']['sha256']]
        assert selected['classification']=='accepted'
        assert selected['residual']['sha256']==updates[0]['residual']['sha256']
    h.save(out, dict(receipt=h.pin(e/'local-run/receipt.json'),collection=receipt['durable_collection']['file'],
        initial=base,selected=updates[0] if updates else None,matrix=matrix['event']['matrix'],
        events=rows,native=[r for r in events if r['kind'].startswith('r8_')],
        cuts=[r for r in events if r['kind']=='r6_cut_audit'],focused_regressions_passed=True))



def review(e, out, h, evidence):
    """Collect compact completed evidence and independently check the native slopes."""
    import os
    from collections import Counter
    from decimal import Decimal
    import numpy as np
    from scipy.linalg import solve_triangular

    assert len(os.sched_getaffinity(0))==1
    summary_path=e/'comparison.json'
    C=h.load(summary_path)
    receipt=h.load(h.bound(C['receipt']))
    manifest=h.load(h.bound(receipt['manifest']))
    audit=h.load(e/'terminal-audit.json')
    assert audit['success'] and receipt['status']=='completed' and receipt['runtime_unchanged']
    assert C['focused_regressions_passed']
    base=C['initial']
    def one(kind):
        return next(x for x in C['events'] if x['kind']==kind)
    def keep(obj, *names):
        return {name:obj[name] for name in names if name in obj}
    def clean(obj):
        return {k:v for k,v in obj.items() if k not in
                ('schema','kind','policy','action_policy','solve_number','sequence','alpha','time')}
    raw=h.load(h.bound(manifest['current_preflight']))
    sys.path.insert(0,str(e/'source/solvers/IncNavierStokesSolver/Tests/RANS'))
    from cls_epawf_m8_r3_resolution import ResolutionAudit
    metric=ResolutionAudit(raw,e/'source')
    R=np.array(h.vector(base['residual']))
    constraints={int(i):v for i,v in np.fromfile(h.bound(base['constraints']),dtype='<f8').reshape(-1,2)}
    block=raw['dimension']
    free=[i for i in range(metric.fields[block].size) if metric.off[block]+i not in constraints]
    ids=metric.off[block]+np.array(free)
    L=np.linalg.cholesky(metric.mass[block][np.ix_(free,free)])
    white=solve_triangular(L,R[ids],lower=True,check_finite=False)
    rate=metric.K0*metric.us/metric.H
    denominator=rate*rate*metric.fields[block].volume
    norm=metric.norms(R,constraints)[block]
    J=np.fromfile(h.bound(C['matrix']),dtype='<f8').reshape(len(R),-1)
    control=one('r10_control')
    independently={}
    directions=[('control',control['direction']),('R9_ascent',manifest['r10_source']['ascent'])]
    directions += [('anchor'+str(a['index']),a['direction']) for a in C['events'] if a['kind']=='r10_anchor']
    for name,pin in directions:
        action=J@np.array(h.vector(pin))
        dwhite=solve_triangular(L,action[ids],lower=True,check_finite=False)
        independently[name]=float(white@dwhite/denominator/norm)
    assert abs(independently['control']-control['native_slope'])<1e-9
    assert abs(independently['R9_ascent']-one('r10_ascent_regression')['slope'])<1e-9
    r9=h.load(evidence/'m8-r9/decision-ledger.json')
    old=h.load(h.bound(r9['binding']['comparison']))
    old_control=next(x for x in old['results'] if x['kind']=='r9_linear_direction')
    reproducibility=dict(control_direction_identical=control['direction']['sha256']==old_control['step']['sha256'])
    live=next(t for t in C['native'] if t['kind']=='r8_trial' and t['comparison']=='r10_linear' and t['fraction']==.125)
    previous=r9['terminal']['snapshot']
    reproducibility.update(control_trial_state_identical=live['state']['sha256']==previous['state']['sha256'],
                           control_trial_residual_identical=live['residual']['sha256']==previous['residual']['sha256'])
    assert all(reproducibility.values())
    selected=C['selected'] or base
    state=h.vector(selected['state'])
    wall=[]
    for precision in (80,110):
        rows=h.wall_rows(raw,state,precision)
        row=min(rows,key=lambda r:Decimal(r['minimum']['value']))
        wall.append(dict(row['minimum'],geometry=row['geometry']))
    assert Decimal(wall[-1]['value'])>0
    assert abs(Decimal(wall[0]['value'])-Decimal(wall[1]['value']))<Decimal('1e-65')
    phases=audit['phases']
    physical=phases['finite_ranges_and_accounting']
    common=physical['common_inventory_physical_qualification']
    trace=phases['scalar_and_projected_trace']
    native_report=phases['native_report']
    assert native_report['available'] and native_report['native_integrator_history_slots_verified']
    assert native_report['residual_matches_durable_bytes']
    assert physical['independent_accounting_identities_met']
    assert physical['native_budget_crosschecks']['passed']
    assert not audit['native_history_advanced']
    tests=(e/'regressions.log').read_text()
    assert '100% tests passed, 0 tests failed out of 3' in tests and 'Ran 9 tests' in tests and '\nOK\n' in tests
    native_trials=[x for x in C['native'] if x['kind']=='r8_trial']
    tangents=[]
    for t in C['events']:
        if t['kind']!='r10_tangent':
            continue
        small=clean(t)
        small['details']=keep(t['details'],'trust_radius','stationarity_l2','active_equality_error',
                              'complementarity','scaled_residual_objective','exchange_passes')
        tangents.append(small)
    cuts=[]
    for c in C['cuts']:
        groups=Counter(r['kind'] for r in c['rows'])
        governing=sorted(c['rows'],key=lambda r:abs(r['multiplier_scaled']),reverse=True)[:3]
        governing += [r for r in c['rows'] if r.get('general_affine_row') and r not in governing]
        cuts.append(keep(c,'name','status','worst_native_violation','worst_extended_violation')|
            dict(active_families=dict(groups),governing=[
                keep(r,'kind','point','general_affine_row','rhs','native_gap','coefficient_gap','transformed_gap','multiplier_scaled')
                for r in governing]))
    small_trials=[]
    for t in native_trials:
        small_trials.append(keep(t,'comparison','name','fraction','classification','actual','linear','matching_trace','covariance','seconds',
                                'state','residual')|dict(ground=t.get('wall'),box=t.get('box'),
                    scalar=keep(t.get('scalar',{}),'admissible','minimum_dimensionless_margin','minimum_lambda','minimum_b')))
    decision=one('r10_decision')
    control_trial=next(t for t in native_trials if t['comparison']=='r10_linear' and t['fraction']==.125)
    source_trial=next(t for t in native_trials if t['comparison']=='r10_source' and t['classification']=='accepted')
    gain_gates=dict(at_least_25_percent_more_gain=decision['source_gain']>=1.25*decision['control_gain'],
        source_absolute_gain=decision['source_gain']>=1e-4*max(base['scaled_block_norms']),
        ground_noninferior=source_trial['wall']['minimum']>=control_trial['wall']['minimum'],
        trace_block_limit=source_trial['actual'][-1]<=1.1*control_trial['actual'][-1],
        gain_per_second=decision['source_gain']/decision['source_seconds']>decision['control_gain']/decision['control_seconds'])
    assert decision['promote']==all(gain_gates.values())
    h.save(out,dict(schema=1,milestone='M8-R10',completed=True,M8_closed=False,
        reviewed_baseline='3b906f7fd2d2e657526ae291165ab2e4d6b0e901',
        producing_revision=receipt['source_commit'],runtime=h.pin(e/'runtime-fingerprint.json'),
        analysis_driver=h.pin(__file__),reuse_compatibility=h.pin(e/'reuse-compatibility.json'),
        binary=receipt['executable'],xml=receipt['session'],comparison=h.pin(summary_path),
        receipt=C['receipt'],collection=C['collection'],budget=h.pin(e/'local-bound/budget.json'),
        bindings=keep(base,'state','history','residual','constraints')|dict(matrix=C['matrix'],
            source_value=manifest['r10_source']['value'],source_action=manifest['r10_source']['action']),
        compatibility=h.load(e/'local-bound/budget.json')['physical_value_and_source_compatibility'],
        R9_reproducibility=reproducibility,
        ascent_regression=clean(one('r10_ascent_regression')),independent_native_slopes=independently,
        scaled_anchor_regression=clean(one('r10_scaled_anchor_regression')),
        accepted_model_anchors=[clean(x) for x in C['events'] if x['kind']=='r10_anchor'],
        anchor_checks=[clean(x) for x in C['events'] if x['kind']=='r10_anchor_feasibility'],
        tangents=tangents,source_action_checks=[clean(x) for x in C['events'] if x['kind']=='r10_moved_action'],
        model_trials=[clean(x) for x in C['events'] if x['kind']=='r10_model_trial'],
        model_outcome=clean(one('r10_model_outcome')),
        original_row_checks=cuts,general_slope_row_scans=[
            clean(x) for x in C['events'] if x['kind']=='r10_general_row_scan'],
        native_trials=small_trials,decision=clean(decision),promotion_gates=gain_gates,
        continuation=dict(run=False,reason='Promotion gate failed' if not decision['promote'] else 'Required next step'),
        work=dict(native_seconds=receipt['commands'][1]['seconds'],trial_count=len(native_trials),
            ranks=1,numerical_threads=receipt['numerical_threads'],compiler_workers=16,
            build=h.pin(e/'build-run/receipt.json')),
        regressions=dict(compile_time_policy_assertions=True,native_ctest_passed=3,predictor_binding_passed=9,
                         log=h.pin(e/'regressions.log'),native_summary_regression_assertions=True),
        source_proposal_accounting=h.pin(e/'source-proposal-accounting.json'),
        final=dict(snapshot=selected['snapshot'],audit=h.pin(e/'terminal-audit.json'),
            native_report=native_report['file'],residual=phases['coupled_residual'],
            trace=trace['independent_trace'],scalar=keep(trace,'admissible','minimum_inner_k','minimum_outer_k','minimum_outer_epsilon',
                                                       'minimum_b','minimum_lambda','sign_unresolved_observations'),
            ground=wall[-1],box_covered=common['qualified_on_common_inventory'],outer_range=common['outer_range'],
            covariance=common['covariance_minima'],accounting=physical['accounting_defects'],
            traction=physical['matching_observations'],physical_steps=0,history_advanced=False)))


def candidate_accounting(e, out, h, evidence):
    """Read-only accounting on the accepted source proposal, distinct from the saved control."""
    import os
    import time
    import numpy as np

    assert len(os.sched_getaffinity(0))==1
    started=time.monotonic()
    C=h.load(e/'comparison.json')
    receipt=h.load(h.bound(C['receipt']))
    manifest=h.load(h.bound(receipt['manifest']))
    candidate, = [x for x in C['native'] if x['kind']=='r8_trial' and
                  x['comparison']=='r10_source' and x['classification']=='accepted']
    sys.path.insert(0,str(e/'source/solvers/IncNavierStokesSolver/Tests/RANS'))
    import cls_epawf_m8_r2_stage_analysis as analysis
    raw=h.load(h.bound(manifest['current_preflight']))
    analyzer=analysis.Analyzer(dict(raw,stage_test_target_time=.0025),h.bound(receipt['session']),e/'source')
    q=np.array(h.vector(candidate['state']))
    h0=np.array(h.vector(C['initial']['history']))
    effective=analyzer.external_lift(q)
    constraints=analyzer.constraints(q)
    loads,quad=analysis.weak_loads(raw,analyzer.fields,effective)
    unmasked=loads+analyzer.time_mass(effective,h0,400.)
    physical=analyzer.physical(q,h0,400.,unmasked=unmasked,constraints=constraints)
    _,E,D,redistribution,*_=physical['source_energy']['robust_source']
    source_error=E-D-redistribution
    energy_error=physical['energy']['direct_remainder']
    assert max(abs(source_error),abs(redistribution),abs(energy_error))<3e-10
    h.save(out,dict(state=candidate['state'],native_residual=candidate['residual'],
        fraction=candidate['fraction'],proposal_only=True,selected_as_saved_state=False,
        source_identity_error=source_error,integrated_redistribution=redistribution,
        energy_weak_minus_direct=energy_error,source=physical['source_energy']['robust_source'],
        energy=physical['energy'],traction=physical['full_total_traction_and_k_flux'],
        native_scalar=candidate['scalar'],native_box=candidate['box'],native_covariance=candidate['covariance'],
        continuous_ground=candidate['wall'],seconds=time.monotonic()-started,
        new_native_trials=0,new_direction_construction_values=0,producing_source=receipt['source_commit']))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=('prepare', 'summary', 'review', 'candidate-accounting'))
    p.add_argument('--epoch', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a=p.parse_args()
    if platform.system()!='Linux':
        raise ValueError('Use the designated Linux server')
    e=a.epoch.resolve()
    evidence=e/'source/docs/plans/cls-epawf-revised-evidence'
    h=module(evidence/'m8-r4/prepare.py','r10_prepare')
    functions={'prepare':prepare, 'summary':summary, 'review':review, 'candidate-accounting':candidate_accounting}
    functions[a.mode](e,a.out.resolve(),h,evidence)


if __name__=='__main__':
    main()
