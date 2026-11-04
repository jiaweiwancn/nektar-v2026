#!/usr/bin/env python3
"""Bind the saved R8 endpoint and verify the fresh R9 comparison."""
import argparse
import importlib.util
from pathlib import Path
import platform
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def prepare(evidence, epoch, out, h):
    r8 = h.load(evidence/'m8-r8/decision-ledger.json')
    pilot = r8['pilot']
    receipt = h.load(h.bound(pilot['receipt']))
    collection = h.load(h.bound(receipt['durable_collection']['file']))
    base = collection['last_snapshot']
    last = pilot['last_matrix']
    assert base['iteration'] == 10 and last['iteration'] == 9
    assert base['snapshot'] == pilot['terminal_snapshot']['snapshot']
    assert base['history'] == r8['target']['original_h0'] and not base['history_advanced']
    for key in ('snapshot', 'state', 'history', 'residual', 'constraints'):
        h.bound(base[key])
    h.bound(last['matrix'])
    previous = h.load(h.bound(receipt['manifest']))
    old_xml = h.bound(receipt['session'])
    tree = ET.parse(old_xml)
    info = tree.find('.//SOLVERINFO')
    selectors = {row.get('PROPERTY'): row for row in info}
    assert selectors['RANSWallDestructionSmoothing'].get('VALUE') == 'Linear'
    assert selectors['RANSWallR8Experiment'].get('VALUE') == 'Pilot'
    for key in ('RANSWallR8Experiment', 'RANSWallR5Updates'):
        info.remove(selectors[key])
    selectors['RANSWallStageTest'].set('VALUE', 'DirectionExperimentR5')
    selectors['RANSWallR3PredictorFile'].set('VALUE', str(out/'predictor.json'))
    selectors['RANSWallR3ProgressFile'].set('VALUE', str(out/'case.r3-progress.jsonl'))
    ET.SubElement(info, 'I', PROPERTY='RANSWallR9Experiment', VALUE='On')
    out.mkdir()
    xml = out/'case.xml'
    tree.write(xml, encoding='unicode')
    assert h.predictor.protected_xml(xml) == previous['protected_physical_xml_sha256']
    runtime_path = epoch/'runtime-fingerprint.json'
    runtime = h.load(runtime_path)
    exe = h.pin(epoch/'build/IncNavierStokesSolver')
    assert runtime['success'] and runtime['runtime_artifacts']['build/IncNavierStokesSolver'] == exe['sha256']
    del previous['frozen_action'], previous['r8_ray']
    h.save(out/'predictor.json', dict(previous, role='retained_r8_terminal_iterate_only',
        state=base['state'], history=base['history'], producing_source_commit=runtime['source_commit'],
        current_executable=exe, current_target_xml=h.pin(xml), current_runtime=h.pin(runtime_path),
        origin=dict(alpha=400, time=.0025, converged=False, inner_stress_model='SaturatedR2',
            iteration=10, snapshot=base['snapshot'], residual=base['residual'],
            constraints=base['constraints'], prior_nonterminal_matrix=last['matrix'])))


def summarize(evidence, epoch, out, h):
    receipt_path = epoch/'local-run/receipt.json'
    receipt = h.load(receipt_path)
    manifest = h.load(h.bound(receipt['manifest']))
    collection = h.load(h.bound(receipt['durable_collection']['file']))
    initial, *updates = collection['completed_iterates']
    assert len(updates) <= 1 and not any(x['history_advanced'] for x in (initial, *updates))
    for key in ('state', 'history', 'residual', 'constraints'):
        expected = (manifest if key in ('state', 'history') else manifest['origin'])[key]
        assert initial[key]['sha256'] == expected['sha256']
    current, = [row for row in collection['jacobians'] if row['status'] == 'complete']
    matrix = current['event']['matrix']
    assert current['iteration'] == 0 and not current['event']['details']['loaded_bound_matrix']
    assert matrix['sha256'] != manifest['origin']['prior_nonterminal_matrix']['sha256']
    assert collection['terminal']['terminal'] == 'r9_fixed_comparison_completed'
    progress = module(evidence/'m8-r7/experiment.py', 'r7_progress').events(collection)
    trials = [x for x in progress if x['kind'] == 'r8_trial']
    assert len(trials) <= 20
    for name in ('r9_linear', 'r9_source'):
        assert sum(t['comparison'] == name for t in trials) <= 9
    decision, = [x for x in progress if x['kind'] == 'r9_decision']
    if updates:
        ray = 'r9_'+decision['selected']
        accepted = [t for t in trials if t['classification'] == 'accepted' and t['comparison'] == ray]
        chosen, = [t for t in accepted if t['state']['sha256'] == updates[0]['state']['sha256']]
        assert chosen['residual']['sha256'] == updates[0]['residual']['sha256']
        assert max(chosen['actual']) == max(updates[0]['scaled_block_norms'])
    h.save(out, dict(receipt=h.pin(receipt_path), initial=initial, current_matrix=matrix,
        matrix_details=current['event']['details'], selected=updates[0] if updates else None,
        results=[x for x in progress if x['kind'].startswith(('r9_', 'r8_'))],
        cut_audits=[x for x in progress if x['kind'] == 'r6_cut_audit']))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=('prepare', 'summary'))
    p.add_argument('--epoch', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    if platform.system() != 'Linux':
        raise ValueError('Use the designated Linux host')
    epoch = a.epoch.resolve()
    evidence = epoch/'source/docs/plans/cls-epawf-revised-evidence'
    h = module(evidence/'m8-r4/prepare.py', 'r4_prepare')
    (prepare if a.mode == 'prepare' else summarize)(evidence, epoch, a.out.resolve(), h)


if __name__ == '__main__':
    main()
