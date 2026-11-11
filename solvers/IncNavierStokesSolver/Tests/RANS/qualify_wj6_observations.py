#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Versioned native order10 row-map supplement; no governing trajectory rerun.

The original order10 quadrature-only packets and order8 qualification remain
immutable. Reuse the frozen driver's measurement API at a new code/output
identity, requesting the existing native passive row-map export at order10.
"""
import argparse
import copy
import json
import math
from pathlib import Path
import sys
import time


def complete_orders(original, supplement):
    from wj_qualification import aggregate, component
    return aggregate(dict(original=original, all_observation_orders=supplement or
                          component('NOT_RUN', 'missing independent-order row qualification')))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'study', 'output', 'driver-scripts', 'solver', 'probe'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--policy', choices=('EP', 'Budget'), required=True)
    a = parser.parse_args()
    # Import the exact already-frozen engine path, rather than a new checkout
    # copy whose path would change the immutable original native design.
    sys.path.insert(0, str(a.driver_scripts.resolve()))
    import numpy as np
    from run_wj6_native import Study
    from wj_immutable import artifacts, file_identity, file_hash, read, write_once, digest, verify_artifacts
    from wj4_measure import comparison
    from wj4_qualification import sensitivity_sequence
    from wj6_qualification import paired_envelope
    from wj6_resources import DEFAULT, volume
    from wj_qualification import aggregate, component
    from run_cls_epawf_picard import require
    started = time.time()
    out = a.output.resolve()/a.policy
    out.mkdir(parents=True, exist_ok=False)
    original_path = a.study.resolve()/a.policy/'summary.json'
    original = json.loads(original_path.read_text())
    require(original['components']['primary']['status'] != 'NOT_RUN', 'original Tx must be complete before the row-map supplement')
    write_once(out/'original-summary.json', original)
    original_gate = original['components']['primary']
    s = Study(argparse.Namespace(output=a.study, policy=a.policy, solver=a.solver.resolve(), probe=a.probe.resolve(), mode='primary'))
    # The existing API is read-only. Only observation output and code lineage
    # are extended; native code, designs, checkpoints and producer paths stay.
    s.root = out
    s.measure_code.update(artifacts([Path(__file__)]))
    s.measure_version = digest(s.measure_code)[:12]
    native_root = a.study.resolve()/a.policy
    gate = component('NOT_RUN')
    try:
        count = sum(len(t['trajectory'])+1 for t in original['trajectories'].values())
        projection = read(native_root/'Tx-resource-projection.json')
        estimate = count*(2*max(max(v) for v in projection['measured_observation_seconds'].values())+
                          projection['parent_seconds_per_observation_pair'])*1.25
        totals = s.ledger.totals()
        fits = totals['seconds']+estimate+600+9*projection['native_seconds_per_step_bounds']['coarse'] <= 10800
        plan = dict(schema='WJ6Order10RowPlanV1', status='PASS' if fits else 'RESOURCE_LIMIT', observations=count,
            estimated_seconds=estimate, closure_seconds=600, transaction_attempts_held=9, totals=totals,
            inputs=artifacts([out/'original-summary.json', native_root/'Tx-resource-projection.json']),
            code=artifacts([Path(__file__)]), governing_attempts=0)
        write_once(out/'resource-plan.json', plan)
        require(fits, 'RESOURCE_LIMIT: complete order10 row-map supplement does not fit')
        augmented, results, identities = {}, {}, {}
        for tag, t in original['trajectories'].items():
            records = [t['initial']]+t['trajectory']
            extra = []
            basis = np.zeros(6); evaluation = np.zeros(6)
            for c in records:
                verify_artifacts({c['qualified_record']['path']: c['qualified_record']})
                require(file_hash(c['report']) == c['report_sha256'] and file_hash(c['path']) == c['sha256'] and
                        file_hash(Path(c['checkpoint'])/'State.bin') == c['checkpoint_sha256'], 'original producer/state/observation changed')
                verify_artifacts({c['native_receipt']['path']: c['native_receipt']})
                verify_artifacts(read(c['native_receipt']['path'])['artifacts'])
                require(volume(native_root)+volume(out)+DEFAULT['invocation_output_bytes'] <= DEFAULT['policy_output_bytes'],
                        'RESOURCE_LIMIT: combined original/supplement output reservation')
                oldq = Path(c['quadrature']['independent']['path'])
                require(file_hash(oldq) == c['quadrature']['independent']['sha256'], 'original independent-order bytes changed')
                d = read(c['report'])
                m, q = s.measure('order10-rows-'+tag+'-'+str(int(c['steps'])), d, Path(c['report']).parent,
                                 nx=t['nx'], assembly=False, order=10, xorder=6, diagnostics=True)
                old = read(oldq)['measurements']
                for key in ('state_sha256', 'stage_history_sha256', 'time', 'represented_initial_sha256'):
                    require(m[key] == old[key], 'row-map supplement changed immutable context: '+key)
                strip = lambda row: {k: v for k, v in row.items() if k != 'constraint_rows'}
                require([strip(r) for r in m['observations']] == [strip(r) for r in old['observations']],
                        'row-map export changed original physical observations')
                require(all(r['constraint_rows'] for r in m['observations']), 'missing native order10 row map')
                q.update(tag=tag, dt=t['dt'], unknowns=c['unknowns'])
                extra.append(q)
                results[tag, int(c['steps'])] = q
                identities[tag+'/'+str(int(c['steps']))] = dict(original=file_identity(oldq), supplement=file_identity(q['path']),
                    physical_samples_exact=True, governing_solves=0)
                for current in (m, read(c['path'])['measurements']):
                    basis = np.maximum(basis, current['evaluation_basis_row_norm_max'])
                    evaluation = np.maximum(evaluation, current['arithmetic_evaluation_sensitivity'])
            packets = [c['arithmetic'] for c in t['trajectory'] if c.get('arithmetic')]
            require(packets, 'missing retained native SVD packets')
            for packet in packets:
                verify_artifacts(packet['input_lineage']);verify_artifacts(packet['code'])
            minimum = {name: min(p['systems'][name]['physical_metric_min_eigenvalue'] for p in packets)
                       for name in ('momentum', 'k', 'epsilon')}
            pressure = min(p['systems']['pressure']['physical_mass_operator_bound'] for p in packets)
            factors = [basis[i]/math.sqrt(minimum[name]) for i, name in enumerate(('momentum', 'momentum', 'k', 'epsilon'))]
            factors += [basis[4]/math.sqrt(pressure)]*2
            row = copy.deepcopy(t)
            originals = [row['initial']]+row['trajectory']
            for c in originals:
                if c.get('arithmetic'):
                    c['arithmetic']['sampled_conversion_factors'] = factors
                    c['arithmetic']['evaluation_sensitivity'] = evaluation.tolist()
                    c['arithmetic']['row_coverage'] = 'both actual order8/order10 native row sets; worst retained coefficient metric, 100x reserve'
            row['trajectory'] = originals[1:]+extra
            augmented[tag] = row
        parts = {}
        for tag in augmented:
            parts[tag+'_order10_accuracy'] = aggregate({str(int(q['steps'])): q['errors']['accuracy_gate']
                for (name, step), q in results.items() if name == tag})
        comparisons = {}
        for group, ta, tb in (('horizontal', 'baseline', 'coarse'), ('timestep', 'coarse', 'fine')):
            floors, boundary, coverage = paired_envelope(augmented[ta], augmented[tb])
            old = original[group+'_comparisons']
            new = {}
            for step in (0, 4, 32):
                ca, cb = results[ta, step], results[tb, 2*step if tb == 'fine' else step]
                new[str(step)] = comparison(read(ca['path'])['measurements'], read(cb['path'])['measurements'], a.policy, ca['rule_settings'])
            parts[group+'_order8_sensitivity'] = sensitivity_sequence(old, floors, s.design, ['0','4','32'], scope_uncertainty={'constrained_boundary': boundary})
            parts[group+'_order10_sensitivity'] = sensitivity_sequence(new, floors, s.design, ['0','4','32'], scope_uncertainty={'constrained_boundary': boundary})
            comparisons[group] = dict(order10=new, order8=old, floors=floors, boundary=boundary, coverage=coverage)
        gate = aggregate(parts)
        write_once(out/'records.json', {tag+'/'+str(step): q for (tag, step), q in results.items()})
        value = dict(schema='WJ6AllObservationOrdersV1', gate=gate, original_summary=file_identity(out/'original-summary.json'),
            inputs=identities, comparisons=comparisons, code=artifacts([Path(__file__)]), governing_attempts=0,
            records=file_identity(out/'records.json'),physical_samples_unchanged=True, seconds=time.time()-started)
        write_once(out/'qualification.json', value)
    except Exception as error:
        gate = component('RESOURCE_LIMIT' if 'RESOURCE_LIMIT' in str(error) else 'OBSTRUCTED', str(error))
        write_once(out/'failure.json', dict(gate=gate, code=artifacts([Path(__file__)])))
        raise
    finally:
        s.s['components']['primary'] = complete_orders(original_gate, gate)
        result=out/('qualification.json' if (out/'qualification.json').exists() else 'failure.json')
        s.s['order10_row_supplement'] = dict(output=str(out), gate=gate, code=artifacts([Path(__file__)]),result=file_identity(result))
        s.ledger.finish_session(s.session)
        s.save()
    print(out/'qualification.json')


if __name__ == '__main__':
    main()
