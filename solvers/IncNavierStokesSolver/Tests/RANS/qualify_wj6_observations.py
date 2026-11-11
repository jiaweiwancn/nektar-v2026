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


def reusable_observation(receipt):
    return receipt.get('exit_code') == 0 and receipt.get('actual_attempts') == 0 and bool(receipt.get('artifacts'))


def failure_status(error, receipt=None):
    return 'RESOURCE_LIMIT' if 'RESOURCE_LIMIT' in str(error) or (receipt or {}).get('failure') else 'OBSTRUCTED'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'study', 'output', 'driver-scripts', 'solver', 'probe'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--policy', choices=('EP', 'Budget'), required=True)
    parser.add_argument('--reuse',type=Path,help='previous versioned observer root, only after a known zero-launch resource denial')
    a = parser.parse_args()
    # Import the exact already-frozen engine path, rather than a new checkout
    # copy whose path would change the immutable original native design.
    sys.path.insert(0, str(a.driver_scripts.resolve()))
    import numpy as np
    from run_wj6_native import Study
    from wj_immutable import artifacts, file_identity, file_hash, read, write_once, digest, verify_artifacts
    from wj4_measure import comparison,evaluate,BUDGET
    from wj4_qualification import sensitivity_sequence
    from wj6_qualification import paired_envelope
    from wj6_resources import DEFAULT, volume, scientific_workers
    from wj_qualification import aggregate, component
    from run_cls_epawf_picard import require
    started = time.time()
    out = a.output.resolve()/a.policy
    out.mkdir(parents=True, exist_ok=False)
    original_path = a.study.resolve()/a.policy/'summary.json'
    prior=a.reuse.resolve()/a.policy if a.reuse else None
    original = json.loads(((prior/'original-summary.json') if prior else original_path).read_text())
    cached={};prior_seconds=0.
    if prior:
        require(prior!=out and read(prior/'failure.json')['gate']['status']=='RESOURCE_LIMIT','reuse needs a separate root and a retained resource failure')
        for path in prior.rglob('receipt.json'):
            r=read(path)
            if r['exit_code']:
                require(r['actual_attempts']==0 and r['resource'].get('not_launched') is True,'only known zero-launch resource denials authorize this observer resume')
                continue
            if reusable_observation(r):cached[path.parent.name]=path.parent
        for command in read(a.root/'command-receipts.json'):
            argv=command['command']
            if '--output' in argv and argv[argv.index('--output')+1]==str(a.reuse.resolve()) and '--policy' in argv and argv[argv.index('--policy')+1]==a.policy:
                prior_seconds+=command['seconds']
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
    current_failed_receipt=None
    try:
        count = sum(len(t['trajectory'])+1 for t in original['trajectories'].values())
        projection = read(native_root/'Tx-resource-projection.json')
        pending=count-len(cached)
        estimate = (pending*2*max(max(v) for v in projection['measured_observation_seconds'].values())+
                   count*projection['parent_seconds_per_observation_pair'])*1.25+60*pending
        totals = s.ledger.totals()
        fits = totals['seconds']+estimate+600+9*projection['native_seconds_per_step_bounds']['coarse'] <= 10800
        plan = dict(schema='WJ6Order10RowPlanV1', status='PASS' if fits else 'RESOURCE_LIMIT', observations=count,
            estimated_seconds=estimate,reused_native_readbacks=len(cached),pending_native_readbacks=pending,
            prelaunch_wait_seconds_per_observation=60,closure_seconds=600, transaction_attempts_held=9, totals=totals,
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
                current_failed_receipt=None
                verify_artifacts({c['qualified_record']['path']: c['qualified_record']})
                require(file_hash(c['report']) == c['report_sha256'] and file_hash(c['path']) == c['sha256'] and
                        file_hash(Path(c['checkpoint'])/'State.bin') == c['checkpoint_sha256'], 'original producer/state/observation changed')
                verify_artifacts({c['native_receipt']['path']: c['native_receipt']})
                verify_artifacts(read(c['native_receipt']['path'])['artifacts'])
                require(volume(native_root)+volume(out)+(volume(prior) if prior else 0)+DEFAULT['invocation_output_bytes'] <= DEFAULT['policy_output_bytes'],
                        'RESOURCE_LIMIT: combined original/supplement output reservation')
                oldq = Path(c['quadrature']['independent']['path'])
                require(file_hash(oldq) == c['quadrature']['independent']['sha256'], 'original independent-order bytes changed')
                name='order10-rows-'+tag+'-'+str(int(c['steps']))
                d = read(c['report'])
                if name in cached:
                    folder=cached[name];receipt=read(folder/'receipt.json');intent=json.loads((folder/'execution-intent.json').read_text())
                    require(receipt['intent_sha256']==digest(intent) and intent['runtime']==s.runtime,'cached observation intent/runtime differs')
                    verify_artifacts(receipt['artifacts']);verify_artifacts(intent['inputs']);verify_artifacts(intent['code'])
                    core=lambda code:{p:h for p,h in code.items() if Path(p).name!='qualify_wj6_observations.py'}
                    require(core(intent['code'])==core(s.measure_code),'cached physical measurement engine differs')
                    settings=intent['settings']
                    require(settings['order']==10 and settings['xorder']==6 and settings['diagnostics'] is True and settings['assembly'] is False and settings['nx']==t['nx'],
                            'cached row/geometry role differs')
                    q=json.loads((folder/'comparison.json').read_text());m=read(folder/'result.json')['measurements']
                    require(q['path']==str(folder/'result.json') and q['report']==c['report'] and q['sha256']==file_hash(folder/'result.json') and q['report_sha256']==c['report_sha256'] and q['state_sha256']==c['state_sha256'],
                            'cached observation state/producer differs')
                    # Recompute metrics from the verified native output; do
                    # not trust an unsealed intermediate comparison cache.
                    q['errors']=evaluate(m,a.policy,xorder=6,rule_settings=q['rule_settings'])
                    leaves={}
                    for kind in ('continuum_error','initial_drift'):
                        for scope in ('full','central','constrained_boundary'):
                            for field,row in q['errors'][kind][scope]['fields'].items():
                                for norm in ('L2','sampled_Linf'):
                                    val=row['normalized_'+norm];limit=BUDGET[norm]
                                    leaves[kind+'/'+scope+'/'+field+'/'+norm]=component('PASS' if math.isfinite(val) and val<=limit else 'OBSTRUCTED',value=val,limit=limit)
                    q['errors']['accuracy_gate']=aggregate(dict(original=q['errors']['accuracy_gate'],all_scopes=aggregate(leaves)))
                    q['reused_native_lineage']=artifacts([folder/'receipt.json',folder/'execution-intent.json',folder/'result.json'])
                else:
                    for attempt in (0,1):
                        current_failed_receipt=None
                        actual_name=name+('-worker-resume1' if attempt else '')
                        wait_start=time.monotonic();allocations=[]
                        while True:
                            allocation=scientific_workers();allocations.append(dict(elapsed=time.monotonic()-wait_start,**allocation))
                            if allocation['workers']<=60:break
                            require(time.monotonic()-wait_start<60,'RESOURCE_LIMIT: bounded passive prelaunch allocation wait')
                            time.sleep(2)
                        write_once(out/('allocation-'+actual_name+'.json'),dict(samples=allocations,reserved_workers=12,limit=72))
                        try:
                            m,q=s.measure(actual_name,d,Path(c['report']).parent,nx=t['nx'],assembly=False,order=10,xorder=6,diagnostics=True)
                            break
                        except Exception:
                            failed=out/('measure-'+s.measure_version)/actual_name/'receipt.json'
                            r=read(failed) if failed.exists() else {}
                            current_failed_receipt=r
                            if attempt or r.get('actual_attempts')!=0 or r.get('resource',{}).get('not_launched') is not True:raise
                            write_once(out/(name+'-zero-launch-resume.json'),dict(receipt=file_identity(failed),governing_attempts=0,
                                reason='one bounded pre-Popen observation retry; no physical step or rejected state is retried'))
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
            records=file_identity(out/'records.json'),physical_samples_unchanged=True,seconds=time.time()-started+prior_seconds,
            current_observer_seconds=time.time()-started,prior_observer_seconds=prior_seconds,reused_readbacks=len(cached),
            previous_observer=file_identity(prior/'failure.json') if prior else None)
        write_once(out/'qualification.json', value)
    except Exception as error:
        gate = component(failure_status(error,current_failed_receipt), str(error))
        write_once(out/'failure.json', dict(gate=gate, code=artifacts([Path(__file__)])))
        raise
    finally:
        s.s['components']['primary'] = complete_orders(original_gate, gate)
        result=out/('qualification.json' if (out/'qualification.json').exists() else 'failure.json')
        if s.s.get('order10_row_supplement'):
            s.s.setdefault('superseded_observers',[]).append(s.s['order10_row_supplement'])
        s.s['order10_row_supplement'] = dict(output=str(out), gate=gate, code=artifacts([Path(__file__)]),result=file_identity(result))
        s.ledger.finish_session(s.session)
        s.save()
    print(out/'qualification.json')


if __name__ == '__main__':
    main()
