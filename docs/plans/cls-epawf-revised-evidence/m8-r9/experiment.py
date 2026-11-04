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


def audit_model(epoch, out, h):
    """Reconstruct source-model predictions at two existing native trials."""
    import os
    import time

    threads = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
               'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')
    assert all(os.environ.get(key) == '1' for key in threads)
    assert len(os.sched_getaffinity(0)) == 1
    started = time.monotonic()
    source = epoch/'source'
    sys.path.insert(0, str(source/'solvers/IncNavierStokesSolver/Tests/RANS'))
    import numpy as np
    from cls_epawf_m8_r2_stage_analysis import Analyzer
    from cls_epawf_m8_r3_resolution import ResolutionAudit, independent_terms

    comparison_path = epoch/'comparison.json'
    comparison = h.load(comparison_path)
    receipt = h.load(h.bound(comparison['receipt']))
    manifest = h.load(h.bound(receipt['manifest']))
    xml = h.bound(receipt['session'])
    preflight = h.bound(manifest['current_preflight'])
    raw = h.load(preflight)
    analyzer = Analyzer(dict(raw, stage_test_target_time=.0025), xml, source)
    audit = ResolutionAudit(raw, source)
    base = comparison['initial']
    results = comparison['results']
    component, = [r for r in results if r['kind'] == 'r9_component']
    refinements = [r for r in results if r['kind'] == 'r9_refinement' and r['valid_source_state']]
    chosen = min(refinements, key=lambda r: r['merit'])
    ray = [t for t in results if t['kind'] == 'r8_trial' and
           t['comparison'] == 'r9_source' and 'residual' in t]
    d = np.array(h.vector(chosen['direction']))
    q = np.array(h.vector(base['state']))
    R = np.array(h.vector(base['residual']))
    rows = {int(i): value for i, value in
            np.fromfile(h.bound(base['constraints']), dtype='<f8').reshape(-1, 2)}
    Jd = np.fromfile(h.bound(comparison['current_matrix']), dtype='<f8').reshape(len(q), -1) @ d
    JDd = np.fromfile(h.bound(component['JD0']), dtype='<f8').reshape(len(q), -1) @ d
    D0 = np.array(h.vector(component['D0']))

    def destruction(state):
        terms, _ = independent_terms(raw, analyzer.fields, analyzer.external_lift(state))
        value = terms['inner_k_destruction']
        value[list(rows)] = 0.
        return value

    Dbase = destruction(q)
    error = audit.norms(Dbase-D0, rows)
    assert max(error) <= 1e-8
    records = []
    for t in (max(ray, key=lambda t: t['fraction']), min(ray, key=lambda t: t['fraction'])):
        f = t['fraction']
        state = np.array(h.vector(t['state']))
        assert np.array_equal(state, q+f*d)
        residual = np.array(h.vector(t['residual']))
        linear = R + f*Jd
        model = linear + destruction(state) - Dbase - f*JDd
        record = dict(fraction=f, state=t['state'], residual=t['residual'],
                      actual=audit.norms(residual, rows).tolist(),
                      linear=audit.norms(linear, rows).tolist(),
                      nonlinear_destruction=audit.norms(model, rows).tolist(),
                      linear_error=audit.norms(residual-linear, rows).tolist(),
                      nonlinear_error=audit.norms(residual-model, rows).tolist())
        assert max(abs(a-b) for a,b in zip(record['actual'], t['actual'])) < 2e-10
        records.append(record)
    h.save(out, dict(comparison=h.pin(comparison_path), source_commit=receipt['source_commit'],
        script=h.pin(__file__), preflight=manifest['current_preflight'], base_component_error=error.tolist(),
        selected_refinement=chosen['index'], independent_term_evaluations=3, new_native_trial_states=0,
        seconds=time.monotonic()-started, trials=records))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=('prepare', 'summary', 'audit-model'))
    p.add_argument('--epoch', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    if platform.system() != 'Linux':
        raise ValueError('Use the designated Linux host')
    epoch = a.epoch.resolve()
    evidence = epoch/'source/docs/plans/cls-epawf-revised-evidence'
    h = module(evidence/'m8-r4/prepare.py', 'r4_prepare')
    if a.mode == 'audit-model':
        audit_model(epoch, a.out.resolve(), h)
    else:
        (prepare if a.mode == 'prepare' else summarize)(evidence, epoch, a.out.resolve(), h)


if __name__ == '__main__':
    main()
