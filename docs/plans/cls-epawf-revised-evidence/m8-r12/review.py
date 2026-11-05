#!/usr/bin/env python3
"""Build the compact R12 review from hash-bound frozen and native results."""
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
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result);return result

def keep(record,*keys):return {k:record[k] for k in keys if k in record}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--epoch',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();assert platform.system()=='Linux'
    e=a.epoch.resolve();ev=e/'source/docs/plans/cls-epawf-revised-evidence';h=module(ev/'m8-r4/prepare.py','r12_review')
    def events(case):return [json.loads(x) for x in (e/case/'case.r3-progress.jsonl').read_text().splitlines()]
    replay=events('replay-bound');complete=events('complete-bound');verify=events('verify-bound')
    final_cpp=e/'tmp/RANSWallEvolutionFeasible.cpp'
    producing_cpp=e/'source/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionFeasible.cpp'
    assert final_cpp.read_text().replace('{"family",J::Quote(kind)}','{"kind",J::Quote(kind)}')==producing_cpp.read_text()
    C=h.load(e/'comparison.json');assert not C['decision']['promote'],'Run required conditional continuation first'
    receipt=h.load(h.bound(C['receipt']));assert C['native_history_verified']
    def candidate(rows,name):return next(x for x in rows if x['kind']=='r12_candidate' and x['name']==name)
    qp=next(x for x in complete if x['kind']=='r12_frozen_outcome')
    assert qp['status']=='qualified'
    ref=h.load(e/'reference-result.json')
    v=[x for x in verify if x['kind']=='r12_candidate'];assert all(x['qualified_affine'] and x['full_anchor_source_domain'] for x in v)
    def trace(rows):
        cuts={int(x['id']):x for x in rows if 'identity' in x and 'row' in x and 'transformed_rhs' in x}
        ps=[x for x in rows if x['kind']=='r12_qp_pass']
        counts=Counter(x['identity'] for x in cuts.values())
        working=Counter(tuple(sorted(cuts[int(i)]['identity'] for i in x['active_ids'])) for x in ps)
        last=next(x for x in reversed(rows) if x['kind']=='r6_cut_audit')
        exchanges=last['exchanges'];released=[x for x in exchanges if x['op']=='release']
        return dict(admissions=len(cuts),scans=len(ps),distinct_rows=len(counts),exact_alias_readmissions=len(cuts)-len(counts),
            repeated_working_sets=sum(x-1 for x in working.values()),releases=len(released),
            zero_or_small_dual_steps=sum(x['dual_blocker']<=1e-12 for x in released),
            final_active=len(last['rows']),final_min_QR=min(last['qr_diagonal']),
            dual_bound_first=ps[0]['nonnegative_dual_bound'],dual_bound_last=ps[-1]['nonnegative_dual_bound'],
            largest_dual_decrease=max(x['nonnegative_dual_bound']-y['nonnegative_dual_bound'] for x,y in zip(ps,ps[1:])),
            final_scan=keep(last,'status','worst_native_violation','worst_extended_violation','maximum_scan_difference'),
            final_active_rows=[keep(x,'id','kind','point','geometry','side','block','axis','projected','general_affine_row',
                'coefficient','rhs','native_gap','coefficient_gap','transformed_gap','multiplier_scaled') for x in last['rows']],
            sample_bounds=[keep(x,'oracle_scan','admitted','nonnegative_dual_bound','worst_native_violation')
                for x in ps if x['oracle_scan'] in (1,24,48,72,90,96,97,120,len(ps))])
    ps=[x for x in complete if x['kind']=='r12_qp_pass']
    final_gate=ps[-1];previous=ps[-2]
    assert final_gate['admitted']==123 and final_gate['worst_native_violation']<=5e-12
    assert previous['admitted']==122 and previous['worst_native_violation']>5e-12
    tests=(e/'regressions.log').read_text();assert '0 tests failed out of 3' in tests and 'Ran 2 tests' in tests and 'Ran 9 tests' in tests and 'Ran 5 tests' in tests
    audit={}
    for name in ('L1','Lf','Nf'):
        path=e/'audits'/str(name+'.json');r=h.load(path)
        audit[name]=keep(r,'status','state','residual','history','fraction','seconds','native_blocks','independent_blocks','difference_blocks',
            'component_value_error','ground','common','accounting','audit')|dict(file=h.pin(path))
        if 'common' in audit[name]:audit[name]['common'].pop('worst_points',None)
    trials=[]
    for r in C['events']:
        if r['kind']!='r8_trial':continue
        t=keep(r,'comparison','name','fraction','classification','cache_reused','actual','linear','matching_trace','state','residual','covariance','wall','seconds','charged_seconds')
        t['box']={k:v for k,v in r.get('box',{}).items() if k!='witnesses'}
        t['scalar']=keep(r.get('scalar',{}),'admissible','minimum_dimensionless_margin','minimum_b','minimum_lambda')
        trials.append(t)
    model=[]
    for r in C['events']:
        if r['kind'] not in ('r11_model_anchor','r11_model_trial','r11_tangent','r11_action','r11_model_work','r11_seed_fallback'):continue
        row={k:v for k,v in r.items() if k not in ('schema','policy','action_policy','sequence','alpha','time','solve_number')}
        if 'details' in row:row['details']=keep(row['details'],'trust_radius','stationarity_l2','active_equality_error','complementarity','exchange_passes','scaled_residual_objective')
        model.append(row)
    selected=next(x for x in C['policies'] if x['name']==C['decision']['retained_policy'])
    native_report=h.load(h.bound(C['native_report']))
    level,=[x for x in native_report['stage_test_levels'] if x.get('attempted')]
    frozen_b=h.load(e/'execution-budget.json');closed=h.load(e/'frozen-investigation-closed.json');assert closed['remaining_seconds']>0
    history=h.load(ev/'m8-r11/decision-ledger.json')
    h.save(a.output,dict(schema=1,milestone='M8-R12',completed=True,M8_closed=False,
        reviewed_baseline='6ab80ffe1895811e958e703c2a4a51e3e0496999',producing_revision=receipt['source_commit'],
        executable=receipt['executable'],runtime=h.pin(e/'runtime-fingerprint.json'),XML=receipt['session'],driver=h.pin(__file__),
        target=h.load(h.bound(receipt['manifest']))['target'],contract=history['contract'],
        initial=keep(C['initial'],'state','history','residual','constraints'),physical_J=C['matrix'],
        source_value=history['bindings']['source_value'],source_action=history['bindings']['source_action'],
        budget=dict(frozen=h.pin(e/'execution-budget.json'),closure=h.pin(e/'frozen-investigation-closed.json'),
            extension=h.pin(e/'extension-declaration.json'),comparison=h.pin(e/'comparison-declaration.json')),
        integration_failures=dict(failed_compile=h.pin(e.parent/'epoch-001/build-run/receipt.json'),
            failed_entries=[h.pin(Path(x['receipt'])) for x in frozen_b['prior_failed_entry_commands']],
            initial_development_clock=h.pin(e.parent/'epoch-002/investigation-clock.json'),
            audit_import_error_seconds=None),
        build=h.pin(e/'build-run/receipt.json'),private_MPI=h.pin(e/'private-mpi-launch.json'),
        diagnostic_only_final_check=dict(final_source=h.pin(final_cpp),producing_source=h.pin(producing_cpp),
            only_difference='R12 original-row event family key; no mathematical or direction change',
            compiler_and_writer_commands=h.pin(e/'metadata-regressions.json'),
            log=h.pin(e/'metadata-regressions.log'),strict_event_payload=h.pin(e/'tests/row-event.jsonl')),
        audit_import_errors=[h.pin(e/'audits'/str(n+'-initial-import-error.log')) for n in ('L1','Lf','Nf')],
        frozen=dict(replay=trace(replay),completion=trace(complete),
            witness=candidate(complete,'L1_witness'),native_candidate=candidate(complete,'Lf_completed'),
            native_KKT=qp,reference=ref,reference_file=h.pin(e/'reference-result.json'),
            reference_status_after_full_oracle='qualified_to_unchanged_numerical_gates',
            complete_oracle_verification=v,verification_receipt=h.pin(e/'verify-run/receipt.json'),
            replay_receipt=h.pin(e/'replay-run/receipt.json'),completion_receipt=h.pin(e/'complete-run/receipt.json'),
            objective_difference=abs(ref['Psi']-qp['details']['scaled_residual_objective']),
            real_trace_limit_regression=dict(limit123=dict(scans=124,admissions=123,passes=True),
                limit122=dict(scans=123,admissions=122,passes=False,worst=previous['worst_native_violation']),extra_cut_admitted=False)),
        comparison=dict(file=h.pin(e/'comparison.json'),receipt=C['receipt'],collection=C['collection'],
            policies=C['policies'],model=model,trials=trials,
            selected_predictions=[x for x in C['events'] if x['kind']=='r11_selected_prediction'],
            promotion=[x for x in C['events'] if x['kind']=='r11_promotion'],decision=C['decision']),
        audits=audit,retained=dict(policy=selected['name'],state=C['selected']['state'],residual=C['selected']['residual'],snapshot=C['selected']['snapshot'],
            native_report=C['native_report'],history_verified=True,physical_steps=0,
            native_accounting=keep(level['audit'],'coupled_norms','strong_continuity','interface','robust_source','energy')),
        outcomes=dict(A='Pre-update limit snapshot omitted final qualification; corrected post-limit candidate remained infeasible; useful new cut discovery continued.',
            B='Independent finite-pool QP plus complete original oracle qualified.',
            C='Native Lf qualified after 123 admissions and 124 scans.',
            D='Lf improves actual gain 4.31x; distinct Nf adds 4.05% and improves trace, at higher cost.',
            E='No promotion, coupled root or physical advancement.'),
        next_decision='Benchmark a genuine feasible-start construction for successive Nf tangents, transporting physical incumbents and rebuilding current-A transformations; retain all ground/trace and original-row gates.',
        remaining_gates=history['remaining_gates'],
        continuation=dict(run=False,reason='No policy passed all declared practical promotion gates'),
        regressions=dict(log=h.pin(e/'regressions.log'),commands=h.pin(e/'regressions.json'),native_CTests=3,predictor_binding=9,selected_candidate=5,coordinates=2,cut_budget_native=True),
        resources=dict(compiler_workers=16,ranks=1,numerical_threads=receipt['numerical_threads'],shared_Picard_preserved=True,observed=h.pin(e/'resource-snapshot.json')),
        producing_receipt_status=receipt['status']))

if __name__=='__main__':main()
