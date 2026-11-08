#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Index immutable Linux AWF evidence without advancing or repairing a state."""
import argparse
from datetime import datetime, timezone
from pathlib import Path
import platform

from run_cls_epawf_picard import load, sha, state_hash, write_json, require


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();root=args.root.resolve()
    require(platform.system()=='Linux' and root.is_relative_to(
        Path('/media/jay/新加卷/nektar-v2026-runs')),'authoritative Linux artifact storage required')
    result=dict(schema=1,captured_UTC=datetime.now(timezone.utc).isoformat(),
        root=str(root),production_governing_solves_added=0,
        scope='read-only content index; terminal results retained; no checkpoint construction or migration',
        artifacts={},short_cases={})

    def record(path):
        path=Path(path)
        require(path.is_file(),'missing indexed artifact: '+str(path))
        item=dict(path=str(path),sha256=sha(path),bytes=path.stat().st_size)
        result['artifacts'][str(path.relative_to(root))]=item
        return item

    fixed=('source-verified-e0d599b.json','forced-recompile-authoritative.json',
        'build-native-authoritative.log','ctest-authoritative.log','runtime-v4/source.json',
        'kernel-independent-self/validation.json','kernel-native-v4/validation.json',
        'native-transactions-v4/summary.json','legacy-active-prefix/summary.json',
        'native-focused-v4/failure-receipt.json','native-constitutive-v6/summary.json',
        'native-constitutive-v6/algebra-accounting.json','native-constitutive-saturated-v6/summary.json',
        'native-failure-reporting-v2-receipt.json','diagnostic-reporting-v6/summary.json',
        'final-solver-cxx-source-sha256.json','final-solver-equivalence/receipt.json',
        'initial-history-provenance.json','accepted-initial-history-provenance.py',
        'source-counterfactual-native-v2/analysis.json',
        'source-counterfactual-native-v2/mass-dual-analysis.json')
    for name in fixed:
        path=root/name
        if path.exists():record(path)
    # Versioned runtime closures and executed command/source receipts are small;
    # libraries, fields and dense matrices themselves stay on the run volume.
    for pattern in ('runtime-v4/*/identity.json','runtime-diagnostic-v*/probe/identity.json',
        'runtime-final-validation/solver/identity.json','final-solver-equivalence/*.json','focused-python/*/manifest.json',
        '*handoff/launch-receipt.json','analysis-suite-receipts/*',
        'analysis-*/compact.json','analysis-*/model-activation-and-interface.json',
        'initial-history-replication/*/input.json','initial-history-replication/*/result.json',
        'initial-history-replication/*/checkpoint/State.bin',
        'status-capture-*/status.json','status-capture-*/*/progress.json',
        'status-capture-*/*/checkpoint/State.bin',
        'status-*-final-*/status.json','compact/*.json'):
        for path in sorted(root.glob(pattern)):
            if path.is_file() and path.resolve()!=args.output.resolve():record(path)

    suites=('primary-2O','primary-05G','primary-05O','saturated-2O',
        'regularized-saturated-2G','regularized-saturated-2O',
        'regularized-05G','regularized-05O','regularized-saturated-05G',
        'regularized-saturated-05O','perturbation-2G','refined-2G','temporal-half-2G')
    audit_names=('admissible','accepted_awf_scalar_domain','accepted_awf_covariance',
        'accepted_awf_sources','inner_stress','inner_constitutive_modification',
        'constitutive_modification','interface','strong_continuity')
    for suite in suites:
        for path in sorted((root/suite).glob('Accepted-*/result.json')):
            native_path=path.parent/'trajectory/case.evolution.json'
            initial_path=path.parent/'initial/case.evolution.json'
            data=load(path);initial=load(initial_path)
            row=dict(status=data['status'],definition=data['definition'],
                actual_h0_sha256=data.get('actual_h0_sha256'),accepted_steps=data.get('steps',0),
                accepted_time=data.get('time',0),failure=data.get('failure') or data.get('initial_failure'),
                result=record(path),initial_report=record(initial_path),
                initial_native_history_sha256=[state_hash(v) for v in initial['native_history']],
                initial_candidate_sha256=state_hash(initial['state']),
                checkpoint_note='The last periodic checkpoint may precede the last accepted state; complete rollback state/history remain in the native report.')
            for file in (path.parent/'initial/definition.json',path.parent/'initial/case.xml'):
                record(file)
            if native_path.exists():
                native=load(native_path);row['native_report']=record(native_path)
                row.update(accepted_steps=native['steps'],accepted_time=native['time'],
                    implicit_callbacks=native['implicit_callback_count'],
                    accepted_coupled_state_sha256=state_hash(native['state']),
                    native_history_sha256=[state_hash(v) for v in native['native_history']],
                    attempted_scheduled_counts={name:sum(t.get(name,0) for t in native['pimple_attempts'])
                        for name in ('outer_sweeps','momentum_predictor_solves','pressure_system_solves',
                            'epsilon_transport_solves','k_transport_solves','bound_calls','bound_activations',
                            'additional_governing_equation_solves')})
                row['audits']={label:{k:v for k,v in audit.items() if k in audit_names}
                    for label,audit in (('initial',native['state_audits'][0]),('final',native['state_audits'][-1]))}
                row['rejections']=[dict(time=e['time'],alpha=e['alpha'],reasons=e['reasons'],
                    state_sha256=state_hash(e['state']),history_sha256=state_hash(e['history']),
                    audit={k:v for k,v in e['audit'].items() if k in audit_names})
                    for e in native['physical_range_events']]
                for attempt in (native['pimple_attempts'][-1:] if native['failure'] else []):
                    packet=Path(attempt['k_snapshot'])
                    if packet.exists():
                        row['failed_frozen_packet']=record(packet);payload=load(packet)
                        if payload.get('matrix_file'):row['failed_transport_A']=record(payload['matrix_file'])
            result['short_cases'][suite]=row
    # Immutable terminal receipts bind complete accepted and rejected payloads.
    for path in sorted(root.glob('*/observations/terminal-attempt-*/receipt.json')):
        record(path)
    for path in sorted(root.glob('*/trajectory/progress.json')):
        progress=load(path)
        if progress['status']!='numerically_rejected':continue
        terminal=Path(progress['segment'])/'case.evolution.json'
        if terminal.exists():
            payload=load(terminal)
            require(payload.get('failure'),'terminal progress lacks native failure')
            record(terminal)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    write_json(args.output,result)
    print(args.output)


if __name__=='__main__':main()
