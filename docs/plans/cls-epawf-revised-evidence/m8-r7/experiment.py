#!/usr/bin/env python3
"""Bind and summarize the bounded R7 comparison on the designated Linux host."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import xml.etree.ElementTree as ET

THREADS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS',
           'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')
sys.dont_write_bytecode = True
for key in THREADS:
    os.environ[key] = '1'


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def events(collection):
    item = collection['progress']
    payload = Path(item['path']).read_bytes()[:item['complete_prefix_bytes']]
    if hashlib.sha256(payload).hexdigest() != item['complete_prefix_sha256']:
        raise ValueError('Changed durable event prefix')
    return [json.loads(line) for line in payload.splitlines()]


def prepare(epoch, cases, h):
    evidence = epoch/'source/docs/plans/cls-epawf-revised-evidence'
    r6 = h.load(evidence/'m8-r6/decision-ledger.json')
    pilot = r6['pilot']
    previous_receipt = h.load(h.bound(pilot['receipt']))
    previous = h.load(h.bound(previous_receipt['manifest']))
    collection = h.load(h.bound(pilot['durable']))
    base = collection['last_snapshot']
    prior_j = pilot['last_completed_matrix']
    h.require(base['iteration'] == 17 and base['snapshot'] == pilot['terminal_state']['snapshot']
              and prior_j['iteration'] == 16 and base['history'] == r6['target']['original_h0']
              and not base['history_advanced'], 'Wrong exact-target terminal snapshot')
    for pin in [base[k] for k in ('state', 'history', 'constraints', 'residual', 'snapshot')] + [prior_j['matrix']]:
        h.bound(pin)
    progress = events(collection)
    direction, = [x for x in progress if x['kind'] == 'r5_direction' and x['iteration'] == 16]
    accepted, = [x for x in progress if x['kind'] == 'r5_trial' and x['iteration'] == 16
                 and x['classification'] == 'accepted_comparison']
    displacement = direction['dimensionless_step'] * accepted['fraction']
    runtime_path = epoch/'runtime-fingerprint.json'
    runtime = h.load(runtime_path)
    exe = h.pin(epoch/'build/IncNavierStokesSolver')
    h.require(runtime['success'] and runtime['runtime_artifacts']['build/IncNavierStokesSolver'] == exe['sha256'],
              'Wrong producing runtime')
    old_xml = h.bound(previous_receipt['session'])
    tree = ET.parse(old_xml)
    info = tree.find('.//SOLVERINFO')
    selectors = {row.get('PROPERTY'): row for row in info}
    info.remove(selectors['RANSWallR5Updates'])
    selectors['RANSWallStageTest'].set('VALUE', 'DirectionExperimentR5')
    selectors['RANSWallR3PredictorFile'].set('VALUE', str(cases/'predictor.json'))
    selectors['RANSWallR3ProgressFile'].set('VALUE', str(cases/'case.r3-progress.jsonl'))
    ET.SubElement(info, 'I', PROPERTY='RANSWallR7Policy', VALUE='Compare')
    ET.SubElement(info, 'I', PROPERTY='RANSWallR7PreviousDisplacement', VALUE=repr(displacement))
    cases.mkdir(parents=True)
    xml = cases/'case.xml'
    tree.write(xml, encoding='unicode')
    h.require(h.predictor.protected_xml(xml) == previous['protected_physical_xml_sha256'],
              'Changed physical target input')
    manifest = dict(previous, state=base['state'], history=base['history'],
        producing_source_commit=runtime['source_commit'], current_executable=exe,
        current_target_xml=h.pin(xml), current_runtime=h.pin(runtime_path),
        origin=dict(alpha=400, time=.0025, converged=False, inner_stress_model='SaturatedR2',
                    iteration=17, snapshot=base['snapshot'], residual=base['residual'],
                    constraints=base['constraints'], prior_nonterminal_matrix=prior_j['matrix']))
    h.save(cases/'predictor.json', manifest)
    probe = ET.parse(old_xml.with_name('case-m4-native-rank0.xml'))
    probe.write(cases/'probe.xml', encoding='unicode')
    h.save(cases/'packet.json', dict(state=h.vector(base['state']), history=h.vector(base['history']),
        alpha=400, time=.0025, only_external_direction=True))


def followup(epoch, root, local, h):
    matched = h.load(local/'comparison.json')
    guard, = [d for d in matched['directions'] if d['name'] == 'r7_positive_ground_control']
    baseline, = [d for d in matched['trials'] if d['name'] == 'r7_r6_radius2'
                 and d['classification'] == 'accepted_comparison']
    accepted, = [t for t in matched['trials'] if t['name'] == guard['name']
                 and t['classification'] == 'accepted_comparison']
    h.require(max(accepted['scaled_block_norms']) < max(baseline['scaled_block_norms']) and
              accepted['box']['maximum'] < baseline['box']['maximum'],
              'Matched guard control does not support a continuation')
    old = h.load(local/'predictor.json')
    runtime_path = epoch/'runtime-fingerprint.json'
    runtime = h.load(runtime_path)
    exe = h.pin(epoch/'build/IncNavierStokesSolver')
    h.require(runtime['success'] and runtime['runtime_artifacts']['build/IncNavierStokesSolver'] == exe['sha256'],
              'Wrong R7 followup executable')
    matrix = matched['matrix']
    h.bound(matrix)
    for kind in ('interior','pilot'):
        case = root/kind
        case.mkdir(parents=True)
        tree = ET.parse(h.bound(old['current_target_xml']))
        info = tree.find('.//SOLVERINFO')
        selectors = {row.get('PROPERTY'):row for row in info}
        selectors['RANSWallR7Policy'].set('VALUE','InteriorProbe' if kind=='interior' else 'Adaptive')
        info.remove(selectors['RANSWallR7PreviousDisplacement'])
        selectors['RANSWallStageTest'].set('VALUE','DirectionExperimentR5' if kind=='interior' else 'ExactTargetR3')
        selectors['RANSWallR3PredictorFile'].set('VALUE',str(case/'predictor.json'))
        selectors['RANSWallR3ProgressFile'].set('VALUE',str(case/'case.r3-progress.jsonl'))
        ET.SubElement(info,'I',PROPERTY='RANSWallR7Radius',VALUE=repr(guard['details']['trust_radius']))
        ET.SubElement(info,'I',PROPERTY='RANSWallR7InteriorTarget',VALUE='On')
        xml = case/'case.xml'
        tree.write(xml,encoding='unicode')
        h.require(h.predictor.protected_xml(xml) == old['protected_physical_xml_sha256'],
                  'Changed R7 followup physical input')
        manifest = dict(old, role='retained_r6_terminal_iterate_only',
            producing_source_commit=runtime['source_commit'], current_executable=exe,
            current_target_xml=h.pin(xml), current_runtime=h.pin(runtime_path))
        if kind=='interior':
            manifest['frozen_action']=dict(matrix=matrix, state_sha256=old['state']['sha256'],
                history_sha256=old['history']['sha256'], residual_sha256=old['origin']['residual']['sha256'],
                assembly_json=json.dumps(matched['assembly']['assembly']))
            (case/'probe.xml').write_bytes((local/'probe.xml').read_bytes())
            (case/'packet.json').write_bytes((local/'packet.json').read_bytes())
        h.save(case/'predictor.json',manifest)


def budget(epoch, root, local, h):
    first = h.load(local/'comparison.json')
    interior = h.load(root/'interior/comparison.json')
    base = max(first['base']['scaled_block_norms'])
    guard, = [x for x in first['trials'] if x['name']=='r7_positive_ground_control'
              and x['classification']=='accepted_comparison']
    inside = next((x for x in interior['trials'] if x['classification']=='accepted_comparison'), None)
    use_margin = bool(inside and inside['box']['maximum']==0. and
                      base-max(inside['scaled_block_norms']) >= .9*(base-max(guard['scaled_block_norms'])))
    case = root/'interior' if use_margin else local
    expected_direction = 'r7_interior_probe' if use_margin else 'r7_positive_ground_control'
    qa = h.load(case/'selected-qualification.json')
    h.require(qa['action_gate'] and qa['direction']==expected_direction, 'Selected direction action gate failed')
    receipt = h.load(h.bound(interior['receipt']))
    pilot = root/'pilot'
    xml = pilot/'case.xml'
    manifest_file = pilot/'predictor.json'
    manifest = h.load(manifest_file)
    tree = ET.parse(xml)
    info = tree.find('.//SOLVERINFO')
    info.find("./I[@PROPERTY='RANSWallR7InteriorTarget']").set('VALUE','On' if use_margin else 'Off')
    ET.SubElement(info,'I',PROPERTY='RANSWallR5Updates',VALUE='32')
    tree.write(xml,encoding='unicode')
    h.require(h.predictor.protected_xml(xml)==manifest['protected_physical_xml_sha256'],
              'Budget changed physical target')
    manifest['current_target_xml']=h.pin(xml)
    manifest_file.write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
    h.save(root/'gate.json',dict(schema=1,source_commit=receipt['source_commit'],
        executable_sha256=receipt['executable']['sha256'],absolute_action_tolerance=1e-7,
        required_state_action_gates_met=True,active_dependencies_covered=True,useful_joint_native_descent=True,
        local_evidence=interior['receipt'],base_state=manifest['state'],maximum_updates=32,
        selected_direction=expected_direction,interior_target=use_margin,comparison=h.pin(local/'comparison.json'),
        selected_qualification=h.pin(case/'selected-qualification.json')))
    t=interior['profiling']
    measured=t['r5_direction']['seconds']+t['factorization']['seconds']+t['line_search']['seconds']
    local_receipt=h.load(h.bound(first['receipt']))
    col=h.load(h.bound(local_receipt['durable_collection']['file']))
    assembly=col['terminal']['profiling']['complete_action_assembly']['seconds']
    (root/'budget.txt').write_text(
        f"One R7 continuation, 32 accepted updates, 900 native seconds, and at most three directions per update. "
        "Full searches have 8+8 trials and partial restoration reserves 16. "
        f"The measured complete J17 took {assembly:.3f}s, and the interior direction, factorization and trials "
        f"took {measured:.3f}s. The independent wall cap governs variable cost. No restart.\n")

def summary(cases, run, h):
    receipt = h.load(run/'receipt.json')
    collection = h.load(h.bound(receipt['durable_collection']['file']))
    manifest = h.load(h.bound(receipt['manifest']))
    base, = collection['completed_iterates']
    h.require(all(base[key]['sha256'] == (manifest[key] if key in ('state','history') else
              manifest['origin'][key])['sha256'] for key in ('state','history','constraints','residual')),
              'Changed exact target snapshot')
    matrix, = collection['jacobians']
    h.require(matrix['iteration'] == 0 and matrix['status'] == 'complete'
        and matrix['event']['matrix']['sha256'] != manifest['origin']['prior_nonterminal_matrix']['sha256'],
        'Missing J17 or reused J16')
    if 'frozen_action' in manifest:
        h.require(matrix['event']['details']['loaded_bound_matrix'] and
            matrix['event']['matrix']['sha256']==manifest['frozen_action']['matrix']['sha256'],
            'Changed J17 in the matched interior probe')
    else:
        h.require(not matrix['event']['details']['loaded_bound_matrix'], 'Missing fresh J17 assembly')
    progress = events(collection)
    base_box = next(x['box'] for x in progress if x['kind'] == 'r5_box')
    h.require(collection['terminal']['terminal'] in ('r7_local_comparison_completed','r7_interior_probe_completed'),
              'Incomplete comparison')
    out = dict(receipt=h.pin(run/'receipt.json'), base=base, base_box=base_box,
        matrix=matrix['event']['matrix'], assembly=matrix['event']['details'],
        directions=[x for x in progress if x['kind'] == 'r5_direction'],
        failures=[x for x in progress if x['kind'] == 'r5_direction_failure'],
        cuts=[x for x in progress if x['kind'] == 'r6_cut_audit'],
        trials=[x for x in progress if x['kind'] == 'r5_trial'],
        audits=[x for x in progress if x['kind'] == 'r5_trial_audit'],
        radius_selection=next((x for x in progress if x['kind'] == 'r7_initial_radius'),None),
        interior_target=next((x for x in progress if x['kind'] == 'r7_fixed_interior_target'),None),
        profiling=collection['terminal']['profiling'])
    h.save(cases/'comparison.json', out)
    for item in out['directions']:
        name = item['name']
        trials = [x for x in out['trials'] if x['name'] == name]
        print(dict(name=name, radius=item['details']['trust_radius'], norm=item['dimensionless_step'],
            guard_work=item['guard_work_over_objective'], linear_E=max(item['model_block_norms']),
            trials=[(t['fraction'], t['classification'], max(t.get('scaled_block_norms',[0])),
                     t.get('box',{}).get('family_maxima'), t.get('model_agreement')) for t in trials]))
    for item in out['failures']:
        print(dict(name=item['name'], radius=item['radius'], failure=item['message']))
    print(out['profiling'])


def qualify(epoch, cases, h, name):
    source = epoch/'source'
    tests = source/'solvers/IncNavierStokesSolver/Tests/RANS'
    sys.path.insert(0, str(tests))
    import numpy as np
    from validate_cls_epawf_m8_r3_actions import validate_case, ValueOracle, CD, LD
    evidence = source/'docs/plans/cls-epawf-revised-evidence'
    r5 = module(evidence/'m8-r5/compare.py', 'r5_compare')
    comparison = h.load(cases/'comparison.json')
    base, matrix = comparison['base'], comparison['matrix']
    direction, = [d for d in comparison['directions'] if d['name'] == name]
    packet = h.load(cases/'packet.json')
    packet['r5_matrix_binding'] = dict(state_sha256=base['state']['sha256'],
        history_sha256=base['history']['sha256'], residual_sha256=base['residual']['sha256'], matrix=matrix,
        assembly_json=[json.dumps(x) for x in comparison['assembly']['assembly']])
    bound = cases/'selected-packet.json'
    h.save(bound, packet)
    probe = epoch/'build/RANSNonpressureActionProbe'
    runtime = h.load(epoch/'runtime-fingerprint.json')
    h.require(h.pin(probe)['sha256'] == runtime['runtime_artifacts']['build/RANSNonpressureActionProbe'],
              'Wrong action probe')
    output = cases/'selected-action.json'
    command = ['taskset','-c','49',str(probe),str(cases/'probe.xml'),str(bound),str(output),
               str(h.bound(direction['step']))]
    with (cases/'selected-action.log').open('x') as log:
        subprocess.run(command, cwd=cases, stdout=log, stderr=subprocess.STDOUT, timeout=120, check=True)
    raw = h.load(output)
    h.require(h.pin(raw['matrix_file']) == matrix and np.asarray(raw['directions'][0]['stage_action'],
        dtype='<f8').tobytes() == h.bound(direction['action']).read_bytes(), 'Changed J17 direction action')
    action = validate_case(raw, source)
    witnesses = r5.bound_actions(raw, direction, ValueOracle, CD, LD, np)
    gate = action['success'] and all(max(r['normalized_action_error'],r['normalized_value_error']) <= 1e-7
                                    for r in witnesses)
    h.save(cases/'selected-qualification.json', dict(action_gate=gate, command=command,
        comparison=h.pin(cases/'comparison.json'), report=h.pin(output), direction=name,
        action=action, box_witnesses=witnesses))
    h.require(gate, 'Selected J17 direction failed independent action gate')


def verify(root, h):
    receipt = h.load(root/'pilot-run/receipt.json')
    collection = h.load(h.bound(receipt['durable_collection']['file']))
    progress = events(collection)
    trials = [x for x in progress if x['kind']=='r5_trial']
    directions = [x for x in progress if x['kind'] in ('r5_direction','r5_direction_failure')]
    accepted = {x['iteration']+1:x for x in trials if x['classification']=='accepted_comparison'}
    updates = {x['details']['iteration']+1:x['details'] for x in progress if x['kind']=='newton_update'}
    states = collection['completed_iterates']
    h.require(len(states)<=33 and collection['protocol']['maximum_Newton_updates']==32 and
              receipt['maximum_seconds']==900, 'Changed pilot hard budgets')
    limits = {'r7_full':8, 'r7_retry_full':8, 'r7_partial':16}
    for i in set(x['iteration'] for x in directions):
        names = [x['name'] for x in directions if x['iteration']==i]
        here = [x for x in trials if x['iteration']==i]
        h.require(len(names)==len(set(names))<=3 and len(here)<=32 and all(name in limits for name in names) and
                  all(sum(t['name']==name for t in here)<=limits[name] for name in names),
                  'R7 exceeded a direction or trial allocation')
    cuts = [x for x in progress if x['kind']=='r6_cut_audit' and x['status']=='passed_original_cut_scan']
    h.require(all(max(x['worst_native_violation'],x['worst_extended_violation'])<=5e-12 for x in cuts),
              'Original direction gate failed')
    previous_box = next(x['box'] for x in progress if x['kind']=='r5_box')
    for before, after in zip(states, states[1:]):
        i=after['iteration']
        trial=accepted[i]
        f=trial['fraction']
        E=max(before['scaled_block_norms'])
        actual=max(after['scaled_block_norms'])
        box=trial['box']
        h.require(trial['state']['sha256']==after['state']['sha256'] and
                  trial['residual']['sha256']==after['residual']['sha256'] and
                  trial['required_relative_decrease']==1e-4*f and actual<(1-1e-4*f)*E and
                  box['maximum']<=(1-1e-4*f)*previous_box['maximum'] and
                  all(x<=y for x,y in zip(box['family_maxima'],previous_box['family_maxima'])) and
                  after['scalar_domain']['admissible'] and after['scalar_domain']['r4_wall_trace']['minimum']>=0.,
                  'Accepted state fails unchanged R6 continuation gates')
        previous_box=box
        if i in updates:
            row=updates[i]
            radius=row['radius']
            rho=row['model_agreement']
            expected=min(2.,1.5*radius) if f==1. and rho>=.75 else max(.125,.5*radius) if f<=.25 or rho<.25 else radius
            h.require(expected==row['next_radius'], 'Adaptive radius transition differs')
    h.save(root/'pilot-policy-verification.json',dict(success=True, receipt=h.pin(root/'pilot-run/receipt.json'),
        durable=receipt['durable_collection']['file'],accepted_updates=len(states)-1,trials=len(trials),
        maximum_trials=max((sum(t['iteration']==i for t in trials) for i in range(len(states))),default=0),
        maximum_directions=max((sum(t['iteration']==i for t in directions) for i in range(len(states))),default=0),
        maximum_cut_violation=max((max(x['worst_native_violation'],x['worst_extended_violation']) for x in cuts),default=0.),
        source_bindings_preserved=True))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('prepare','summary','qualify','followup','budget','verify'))
    parser.add_argument('--epoch', type=Path, required=True)
    parser.add_argument('--cases', type=Path, required=True)
    parser.add_argument('--run', type=Path)
    parser.add_argument('--direction')
    parser.add_argument('--local',type=Path)
    args = parser.parse_args()
    if platform.system() != 'Linux':
        raise ValueError('Use the designated Linux server')
    epoch, cases = args.epoch.resolve(), args.cases.resolve()
    h = module(epoch/'source/docs/plans/cls-epawf-revised-evidence/m8-r4/prepare.py','r4_prepare')
    if args.mode == 'prepare':
        prepare(epoch, cases, h)
    elif args.mode == 'summary':
        summary(cases, args.run.resolve(), h)
    elif args.mode == 'qualify':
        qualify(epoch, cases, h, args.direction)
    elif args.mode == 'followup':
        followup(epoch, cases, args.local.resolve(), h)
    elif args.mode == 'budget':
        budget(epoch, cases, args.local.resolve(), h)
    else:
        verify(cases,h)


if __name__ == '__main__':
    main()
