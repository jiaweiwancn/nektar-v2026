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


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=('prepare', 'summary'))
    p.add_argument('--epoch', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a=p.parse_args()
    if platform.system()!='Linux':
        raise ValueError('Use the designated Linux server')
    e=a.epoch.resolve()
    evidence=e/'source/docs/plans/cls-epawf-revised-evidence'
    h=module(evidence/'m8-r4/prepare.py','r10_prepare')
    (prepare if a.mode=='prepare' else summary)(e,a.out.resolve(),h,evidence)


if __name__=='__main__':
    main()
