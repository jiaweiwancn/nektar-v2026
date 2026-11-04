#!/usr/bin/env python3
"""Prepare the exact R4 terminal replay and qualify the matched R5 comparison."""
import argparse
from decimal import Decimal
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import xml.etree.ElementTree as ET

THREADS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
           'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')
sys.dont_write_bytecode = True
for name in THREADS:
    os.environ[name] = '1'


def helpers(source):
    path = source / 'docs/plans/cls-epawf-revised-evidence/m8-r4/prepare.py'
    spec = importlib.util.spec_from_file_location('r4_prepare', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def prepare(epoch, cases, h):
    source = epoch / 'source'
    old = h.load(source / 'docs/plans/cls-epawf-revised-evidence/m8-r4/decision-ledger.json')
    pilot = old['pilot']
    previous = h.load(h.bound(pilot['manifest']))
    raw = h.load(h.bound(pilot['native_report']))
    collection = h.load(h.bound(pilot['durable']))
    base = collection['completed_iterates'][-1]
    h.require(base['iteration'] == 12 and base['snapshot'] == pilot['terminal_snapshot'] and
              base['state'] == pilot['state'] and base['residual'] == pilot['residual'] and
              base['constraints'] == pilot['constraints'] and base['alpha'] == 400 and base['time'] == .0025 and
              not base['history_advanced'], 'Wrong exact terminal numerical base')
    prior_matrix = collection['jacobians'][-1]
    h.require(prior_matrix['iteration'] == 11, 'Historical matrix association differs')
    h.bound(prior_matrix['event']['matrix'])
    for key in ('snapshot', 'constraints', 'residual'):
        h.bound(base[key])
    state, history = h.vector(base['state']), h.vector(base['history'])
    h.require(base['history']['sha256'] == h.predictor.H0, 'Original h0 changed')
    wall, high = (h.wall_rows(raw, state, precision) for precision in (80, 110))
    h.require(max(abs(Decimal(p['value']) - Decimal(q['value'])) for a, b in zip(wall, high, strict=True)
                  for p, q in zip(a['samples'], b['samples'], strict=True)) < Decimal('1e-65'),
              'Wall minima precision ladder failed')
    minimum = min((row['minimum'] for row in high), key=lambda row: Decimal(row['value']))
    h.require(Decimal(minimum['value']) > Decimal('1e-60'), 'R4 terminal wall sign is not positive')
    runtime_path = epoch / 'runtime-fingerprint.json'
    runtime = h.load(runtime_path)
    exe = h.pin(epoch / 'build/IncNavierStokesSolver')
    h.require(runtime['success'] and runtime['runtime_artifacts']['build/IncNavierStokesSolver'] == exe['sha256'],
              'Wrong producing runtime')
    old_xml = h.bound(previous['current_target_xml'])
    for phase in ('local', 'pilot'):
        case = cases / phase
        case.mkdir(parents=True)
        tree = ET.parse(old_xml)
        info = tree.find('.//SOLVERINFO')
        selectors = {row.get('PROPERTY'): row for row in info}
        selectors['RANSWallStageTest'].set('VALUE', 'DirectionExperimentR5' if phase == 'local' else 'ExactTargetR3')
        selectors['RANSWallGlobalization'].set('VALUE', 'JointRestorationR5')
        selectors['RANSWallR3PredictorFile'].set('VALUE', str(case / 'predictor.json'))
        selectors['RANSWallR3ProgressFile'].set('VALUE', str(case / 'case.r3-progress.jsonl'))
        xml = case / 'case.xml'
        tree.write(xml, encoding='unicode')
        digest = h.predictor.protected_xml(xml)
        h.require(digest == previous['protected_physical_xml_sha256'], 'Frozen stage input changed')
        manifest = dict(previous, role='retained_r4_terminal_iterate_only', state=base['state'],
                        history=base['history'], current_executable=exe, current_target_xml=h.pin(xml),
                        current_runtime=h.pin(runtime_path), producing_source_commit=runtime['source_commit'],
                        origin=dict(alpha=400, time=.0025, converged=False, inner_stress_model='SaturatedR2',
                                    iteration=12, snapshot=base['snapshot'], residual=base['residual'],
                                    constraints=base['constraints'], prior_nonterminal_matrix=prior_matrix['event']['matrix']))
        h.save(case / 'predictor.json', manifest)
    probe = ET.parse(old_xml.with_name('case-m4-native-rank0.xml'))
    info = probe.find('.//SOLVERINFO')
    selector = info.find("./I[@PROPERTY='RANSWallGlobalization']")
    if selector is None:
        selector = ET.SubElement(info, 'I', PROPERTY='RANSWallGlobalization')
    selector.set('VALUE', 'JointRestorationR5')
    probe.write(cases / 'probe.xml', encoding='unicode')
    h.save(cases / 'packet.json', dict(state=state, history=history, alpha=400, time=.0025,
                                       only_external_direction=True))
    h.save(cases / 'replay-wall.json', dict(base=base['snapshot'], minimum=minimum, wall=high,
        precision_digits=[80, 110], historical_previous_matrix=prior_matrix['event']['matrix']))


def bound_actions(raw, proposal, oracle_type, complex_type, real_type, np):
    oracle = oracle_type(raw, cutoff_relative_step='1e-4')
    state = np.asarray(raw['state'], dtype=complex_type)
    direction = np.asarray(raw['directions'][0]['direction'], dtype=complex_type)
    state[list(oracle.external)] = [oracle.external[i] for i in oracle.external]
    direction[list(oracle.external)] = 0
    records = []
    witnesses = proposal['base_witness_actions']['witnesses']
    for exponent in ('1e-20', '1e-30'):
        step = real_type(exponent)
        perturbed = state + complex_type(1j)*step*direction
        oracle.cache = {}
        projected = perturbed.copy()
        for index, value in oracle.trace(perturbed).items():
            projected[oracle.off[3]+index] = value
        for w in witnesses:
            family = w['family']
            i, j = w['component']
            field = (0 if i == 0 else 1) if family == 0 else 3 if family == 2 else 2
            values = oracle.fields[field].values(projected if w['projected_trace'] else perturbed,
                oracle.off[field], [w['point']], geometry=w['element_geometry'])
            measured = values[4][0, j] if family == 0 else values[3][0]
            derivative = float(measured.imag / step / real_type(w['reference_scale']))
            value = float(measured.real / real_type(w['reference_scale']) - (w['reference_component'] if family == 0 else 0))
            width = (.001, .01, .02, .01)[family]
            records.append(dict(complex_step=exponent, witness=w, independent_value=value,
                independent_derivative=derivative, normalized_value_error=abs(value-w['value'])/width,
                normalized_action_error=abs(derivative-w['directional_quantity_change'])/width))
    return records


def qualify(epoch, cases, h, cpu):
    source = epoch / 'source'
    tests = source / 'solvers/IncNavierStokesSolver/Tests/RANS'
    sys.path.insert(0, str(tests))
    import numpy as np
    from validate_cls_epawf_m8_r3_actions import validate_case, ValueOracle, CD, LD

    receipt_path = cases / 'local-run/receipt.json'
    receipt = h.load(receipt_path)
    h.require(receipt['status'] == 'completed' and receipt['runtime_unchanged'], 'Incomplete local replay')
    collection = h.load(h.bound(receipt['durable_collection']['file']))
    manifest = h.load(h.bound(receipt['manifest']))
    base, = collection['completed_iterates']
    h.require(all(base[key]['sha256'] == (manifest[key] if key in ('state', 'history') else
              manifest['origin'][key])['sha256'] for key in ('state', 'history', 'residual', 'constraints')),
              'Replay changed the exact retained state/history/residual/constraints')
    matrix, = collection['jacobians']
    h.require(matrix['iteration'] == 0 and matrix['status'] == 'complete' and
              matrix['event']['matrix']['sha256'] != manifest['origin']['prior_nonterminal_matrix']['sha256'],
              'No current independent terminal-state matrix')
    packet = h.load(cases / 'packet.json')
    packet['r5_matrix_binding'] = dict(state_sha256=base['state']['sha256'], history_sha256=base['history']['sha256'],
        residual_sha256=base['residual']['sha256'], matrix=matrix['event']['matrix'],
        assembly_json=[json.dumps(row) for row in matrix['event']['details']['assembly']])
    packet_path = cases / 'bound-packet.json'
    if packet_path.exists():
        h.require(packet_path.read_text() == json.dumps(packet, indent=2, allow_nan=False)+'\n',
                  'Existing action packet differs')
    else:
        h.save(packet_path, packet)
    progress = collection['progress']
    payload = Path(progress['path']).read_bytes()[:progress['complete_prefix_bytes']]
    h.require(hashlib.sha256(payload).hexdigest() == progress['complete_prefix_sha256'],
              'Local progress prefix changed')
    events = [json.loads(line) for line in payload.splitlines()]
    directions = [row for row in events if row['kind'] == 'r5_direction']
    trials = [row for row in events if row['kind'] == 'r5_trial']
    failures = [row for row in events if row['kind'] == 'r5_direction_failure']
    probe = epoch / 'build/RANSNonpressureActionProbe'
    runtime = h.load(epoch / 'runtime-fingerprint.json')
    h.require(h.pin(probe)['sha256'] == runtime['runtime_artifacts']['build/RANSNonpressureActionProbe'],
              'Current action probe identity differs')
    verification = {}
    combined = None
    for proposal in directions:
        name = proposal['name']
        output = cases / (name+'-action.json')
        command = ['taskset', '-c', str(cpu), str(probe), str(cases/'probe.xml'),
                   str(cases/'bound-packet.json'), str(output), str(h.bound(proposal['step']))]
        with (cases/(name+'-action.log')).open('x') as log:
            subprocess.run(command, cwd=cases, stdout=log, stderr=subprocess.STDOUT, timeout=120, check=True)
        raw = h.load(output)
        row = raw['directions'][0]
        row['name'] = name
        h.require(h.pin(raw['matrix_file']) == matrix['event']['matrix'] and
            np.asarray(row['stage_action'], dtype='<f8').tobytes() == h.bound(proposal['action']).read_bytes(),
            'Recorded direction action differs from the current terminal-state matrix product')
        if combined is None:
            combined = dict(raw, directions=[])
        combined['directions'].append(row)
        verification[name] = dict(report=h.pin(output), command=command,
            outer_witness_actions=bound_actions(raw, proposal, ValueOracle, CD, LD, np))
    qualification = validate_case(combined, source) if combined is not None else None
    action_pass = bool(qualification and qualification['success'] and all(
        max(row['normalized_action_error'], row['normalized_value_error']) <= 1e-7
        for value in verification.values() for row in value['outer_witness_actions']))
    local_terminal = collection['terminal']
    picked = next((row for row in trials if row['name'] != 'r4_control' and
                   row['classification'] == 'accepted_comparison'), None)
    h.require(bool(picked) == local_terminal['useful_joint_descent'], 'Local merit evidence differs')
    audit = dict(schema=1, local_receipt=h.pin(receipt_path), matrix=matrix['event']['matrix'],
                 base_snapshot=base['snapshot'], base_box=next(row['box'] for row in events if row['kind'] == 'r5_box'),
                 direction_qualification=qualification, witness_qualification=verification, action_gate=action_pass,
                 direction_failures=failures, chosen_trial=picked,
                 comparison_events=[row for row in events if row['kind'] in ('r5_direction', 'r5_trial', 'r5_trial_audit')],
                 profiling=local_terminal['profiling'])
    h.save(cases / 'comparison.json', audit)
    print(json.dumps(dict(action_gate=action_pass, useful_joint_native_descent=bool(picked),
                          directions=list(verification), profiling=local_terminal['profiling']), indent=2))


def stage_budget(cases, h):
    comparison_path = cases / 'comparison.json'
    audit = h.load(comparison_path)
    h.require(audit['action_gate'] and audit['chosen_trial'], 'Local evidence did not pass the continuation gate')
    receipt = h.load(h.bound(audit['local_receipt']))
    t = audit['profiling']
    trial = t['line_search']['seconds'] / t['line_search']['calls']
    # Count even the R4 control time against each update's direction budget.
    # Six trials is a cost estimate, not a change to the 32-trial search cap.
    per_update = (t['complete_action_assembly']['seconds'] + t['factorization']['seconds']
                  + t['r5_direction']['seconds'] + 6*trial + 1.)
    maximum = min(24, int((900-90)/per_update))
    h.require(maximum >= 4, 'Measured direction and trial cost leaves an impractical pilot budget')
    pilot = cases / 'pilot'
    xml, manifest_file = pilot / 'case.xml', pilot / 'predictor.json'
    manifest = h.load(manifest_file)
    tree = ET.parse(xml)
    info = tree.find('.//SOLVERINFO')
    h.require(info.find("./I[@PROPERTY='RANSWallR5Updates']") is None, 'Pilot budget is already set')
    ET.SubElement(info, 'I', PROPERTY='RANSWallR5Updates', VALUE=str(maximum))
    tree.write(xml, encoding='unicode')
    h.require(h.predictor.protected_xml(xml) == manifest['protected_physical_xml_sha256'],
              'Pilot budget changed the frozen physical target')
    manifest['current_target_xml'] = h.pin(xml)
    # This manifest was staged during prepare and has never been used.
    manifest_file.write_text(json.dumps(manifest, indent=2, allow_nan=False)+'\n')
    h.save(cases/'r5-gate.json', dict(schema=1, source_commit=receipt['source_commit'],
        executable_sha256=receipt['executable']['sha256'], absolute_action_tolerance=1e-7,
        required_state_action_gates_met=True, active_dependencies_covered=True, useful_joint_native_descent=True,
        local_evidence=audit['local_receipt'], comparison=h.pin(comparison_path), base_state=manifest['state'],
        maximum_updates=maximum, maximum_trials_per_update=32, maximum_native_seconds=900,
        estimated_seconds_per_update=per_update, setup_and_audit_allowance_seconds=90))
    (cases/'budget.txt').write_text(
        f"One exact-target continuation, at most {maximum} accepted updates, 32 trials per update, and 900 native seconds. "
        f"Measured local assembly, factorization, all directions, six trials at the observed mean cost, and one second "
        f"for publication give {per_update:.3f}s per update. Reserve 90s for setup and terminal audit. "
        "The wall cap applies independently of the measured estimate. No automatic restart.\n")

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=('prepare', 'qualify', 'budget'))
    p.add_argument('--epoch', type=Path, required=True)
    p.add_argument('--cases', type=Path, required=True)
    p.add_argument('--cpu', type=int, default=0)
    a = p.parse_args()
    if platform.system() != 'Linux':
        raise ValueError('Use the designated Linux host')
    epoch, cases = a.epoch.resolve(), a.cases.resolve()
    h = helpers(epoch / 'source')
    if a.mode == 'prepare':
        prepare(epoch, cases, h)
    elif a.mode == 'qualify':
        qualify(epoch, cases, h, a.cpu)
    else:
        stage_budget(cases, h)


if __name__ == '__main__':
    main()
