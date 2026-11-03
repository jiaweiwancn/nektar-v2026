"""Assemble compact final M8-R1 receipts from immutable curated evidence.

Pure field selection and byte/canonical hashing; no solver or constitutive
re-evaluation. Component ledgers retain the complete independent observations.
"""
from datetime import datetime, timezone
import hashlib, json
from pathlib import Path

here = Path(__file__).resolve().parent
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
load = lambda name: json.loads((here/name).read_text())
a, b1, b2 = (load(name) for name in ('archive-ledger.json', 'nonlinear-b1-ledger.json', 'nonlinear-b2-ledger.json'))
points, preflight, controls, kernel = (load(name) for name in (
    'archive-points-summary.json', 'preflight-summary.json', 'controls-ledger.json', 'kernel-summary.json'))
assert a['success'] and b1['success'] and b2['success'] and points['success'] and kernel['verification_success']
assert points['completed_epoch_count'] == 12 and points['points_file']['records'] == 18180
assert b2['fresh_jacobian_same_stage_retry']['classification_counts'] == dict(
    accepted_descent=0, unevaluable_residual_or_norm=21, valid_non_descent=0)
assert not b2['fresh_jacobian_same_stage_retry']['converged']
assert not b2['fresh_jacobian_same_stage_retry']['history_advanced']
assert b2['observed_Jacobian_accuracy_failures']
runtime = load('runtime-summary.json')
assert all(all(v.values()) for v in runtime['before_after_comparisons'].values())
epochs = []
for item in a['states']:
    inner = [item['sets'][name] for name in ('source_volume_inner', 'physical_sample_K', 'interface_inner')]
    outer = [item['sets'][name] for name in ('source_volume_outer', 'physical_sample_outer', 'interface_outer')]
    epochs.append(dict(time=item['time'], accepted_step=item['accepted_step'],
        state_sha256=item['state_sha256'], history_sha256=item['history_sha256'],
        minimum_inner_k=min(v['minimum_k'] for v in inner), minimum_outer_k=min(v['minimum_k'] for v in outer),
        minimum_inner_covariance=min(v['minimum_covariance'] for v in inner),
        minimum_outer_covariance=min(v['minimum_covariance'] for v in outer),
        minimum_outer_epsilon=min(v['minimum_epsilon'] for v in outer),
        source_minimum_b=item['ray_summary']['minimum_b'], trace_minimum_b=item['trace_summary']['minimum_b'],
        minimum_lambda=item['ray_summary']['minimum_lambda'],
        limiter_active_volume=item['robust_source_native'][4], normalization_capped_volume=item['robust_source_native'][5],
        robust_source=item['robust_source_native'], historical_range_membership=item['constitutive_range_membership']))
events = [
    dict(category='representation_error', time=0, scope='original accepted initialization',
        observation='At the first inner witness k=K but represented Gxz=5.19159 versus analytic1.89599; original mesh425 unknowns.',
        evidence='archive-ledger.json#/events/0/detail', causal_status='Initial covariance loss explained by excessive represented shear, not local depletion.'),
    dict(category='inner_realizability_loss', time=0, scope='archived constitutive finite sets',
        detail=next(e for e in a['events'] if e['type']=='first_negative_inner_covariance')['detail'],
        causal_status='Already present before any time step; not repaired by capped destruction.'),
    dict(category='outer_constitutive_range_exit', time=0, scope='both historical local boxes and new restricted region',
        observation='Initial maximum outer gradient departure .0183323283; improved representation .00341233869 still exceeds frozen .001.',
        evidence=['archive-ledger.json#/certified_ranges', 'preflight-summary.json', 'kernel-summary.json#/representation_coverage'],
        causal_status='Coverage failure; not a proof of principal failure at these particular states.'),
    dict(category='source_limiter_activation', time=.045, scope='archived accepted trajectory',
        observation='First active source volume .001207819098335217; cap volume0/lambda1 throughout accepted prefix.',
        evidence='archive-ledger.json#/states', causal_status='No observed cap-induced trace attenuation.'),
    dict(category='pressure_interface_imbalance', time=.055, scope='last accepted state',
        pressure=a['states'][-1]['pressure'], continuity=a['states'][-1]['continuity'],
        interface=a['states'][-1]['interface'], constraints=a['states'][-1]['constraints'], energy=a['states'][-1]['energy'],
        causal_status='Strong continuity and full traction/flux defects are nonzero despite small weak/gauge residuals; terminal causality unproved. No mechanical-pressure continuity imposed.'),
    dict(category='scalar_inadmissibility', time=.06, scope='unaccepted Newton state, first recorded native-GLL violation at iteration8',
        observation='Native inner k=-1.1713577454207364e-7 at iteration8, then -1.0074658877473182e-6 at terminal base; source GL k remains positive.',
        evidence='nonlinear-b2-ledger.json#/original_policy/iterations',
        causal_status='Internal state observation, never an accepted physical time stage.'),
    dict(category='Jacobian_globalization_failure', time=.06, scope='original policy and same-stage fresh diagnostic',
        original_counts=b2['original_policy']['classification_counts'],
        fresh_counts=b2['fresh_jacobian_same_stage_retry']['classification_counts'],
        terminal=b2['fresh_jacobian_same_stage_retry']['terminal'],
        evidence=['nonlinear-b1-ledger.json', 'nonlinear-b2-ledger.json'],
        causal_status='All21 terminal fractions fail actual complete-ray nonnegative-k evaluation. Fresh assembly/retry does not rescue; refresh-policy causality not established; no factor/solve failure.'),
    dict(category='unresolved_cause', time=.06, scope='independent Jacobian action/solution existence',
        detail=b2['observed_Jacobian_accuracy_failures'],
        observation='Pressure/gauge directional accuracy discrepancy retained; existence of a physically admissible stage root and a positive time map unresolved.'),
    dict(category='bounded_spatial_dt_controls', scope='independent fresh starts, dt.005 and.0025 on1571 unknowns',
        observation='Both fail their first nonlinear stage; accepted physical time remains0. All continuation diagnostic from initial applicability exit.',
        evidence='controls-ledger.json', causal_status='No convergence-order, transition, long-time or terrain conclusion.'),
]
selected = ['archive-ledger.json', 'archive-curation-index.json', 'archive-points-index.json',
    'archive-points-summary.json', 'nonlinear-b1-ledger.json', 'nonlinear-b2-ledger.json',
    'nonlinear-curation-index.json', 'preflight-summary.json', 'preflight-curation-index.json',
    'controls-ledger.json', 'controls-curation-index.json', 'controls-focused-guard-identities.json',
    'controls-analysis-runtime.json', 'controls-notes.md',
    'runtime-summary.json', 'source-binding-receipt.json', 'focused-checks-receipt.json', 'kernel-summary.json',
    'kernel-curation-index.json', 'kernel-contract-preregistered.md', 'candidate-contract.md']
result = dict(schema=1, milestone='M8-R1', status='completed_bounded_diagnosis', M8_status='OPEN',
    physical_model_qualified=False, production_constitutive_adoption=False,
    reviewed_baseline='ab5c3c6fee19761e8d88692e8e6c331e5889234e', initial_HEAD_equal_reviewed_baseline=True,
    branch='codex/cls-epawf-revised', compiled_source_epochs=[b1['producing_runtime'], b2['producing_runtime']],
    committed_helper_source_binding=load('source-binding-receipt.json'),
    historical_producing_runtime=a['producing_runtime'],
    original_XML_log_results=a['inputs'], original_checkpoint_exported=False,
    original_rejected_Newton_vector_exported=False,
    complete_archived_history_hashes=a['history'],
    hash_encoding='File references are SHA256 of exact bytes. Canonical vectors use each linked ledger declared encoding; no array is inferred from its norm.',
    source_runtime_receipt='runtime-summary.json (100 runtime images and3921 numerical inputs per instrumented epoch, exact before/after matches)',
    archived_first_failure=dict(last_accepted_time=.055, first_failed_stage_time=.06, accepted_epoch_count=12,
        first_inner_witness=next(e for e in a['events'] if e['type']=='first_negative_inner_covariance'),
        first_negative_outer_covariance=None, accepted_outer_covariance_scope='finite sets only',
        accepted_epochs=epochs, complete_point_export=points['points_file']),
    failing_iterate_receipts={name: dict(file=file, vectors=data['failed_stage_vectors_canonical_sha256'],
        raw_inputs=data['original_input_receipts'], dense_matrices=data['dense_matrices_retained_on_linux'])
        for name, file, data in [('B1', 'nonlinear-b1-ledger.json', b1), ('B2', 'nonlinear-b2-ledger.json', b2)]},
    bounded_matrix=dict(used_slots=[1,2,3,4,5,6], maximum=6, additional_evolution_slots=0,
        native_evolution_initializations=4, same_stage_diagnostic_allocations=2,
        zero_step_preflight_initializations=4, zero_step_focused_guard_checks=2,
        native_solver_cases_sequential=True,
        definition='Six conservative allocations count both same-stage attempts, including the failed strict fresh assembly; preflight/guard zero-step calls and direct unit probes are identified separately.'),
    events=events, focused_regressions='Three existing native unit probes and two zero-step guards passed; exact commands in focused-checks-receipt.json.',
    next_model=dict(choice='a', contract='candidate-contract.md', restricted_kernel_feasible=True,
        actual_initial_trajectory_coverage=False, original_broad_box_status='FAILED unchanged',
        certified_outer_principal_lower_over_nut0=kernel['outer_proof']['certified_principal_lower_over_nut0'],
        actual_dependent_trace='lambda*Epre(H), includingchi0',
        lambda_to_zero='Excluded by explicit applicability; no admissible singular high-Re continuation supplied.',
        necessary_native_changes='Nonlinear inner AWF weak stress/source/trace and full traction operator; no production adoption in R1.'),
    unresolved=['Pressure/gauge derivative discrepancy', 'Existence of an admissible failed-stage root',
        'Complete initial/trajectory coverage of the restricted region',
        'Continuum and fully discrete scalar/turbulence admissibility',
        'Exact-zero inner-k covariance coverage in the optional finite guard',
        'lambda approaching zero outer continuation'],
    later_gates_open=['Genuine curved MMS for an adopted law', 'Active capped-source curved/variable-H coupling',
        'Outer and fully discrete turbulence admissibility', 'Completed unforced perturbation/refinement and long-time qualification',
        'Scalable feasibility beyond dense guard', 'Paired developed floor/hill and reserved physical scores'],
    components={name: dict(sha256=sha(here/name), bytes=(here/name).stat().st_size) for name in selected},
    assembler_sha256=sha(__file__), assembled_utc=datetime.now(timezone.utc).isoformat(), numerical_re_evaluation=False)
(here/'first-failure-ledger.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
lines = ['| Event category | Time / scope | Observation | Causal status |', '|---|---|---|---|']
for event in events:
    observation = event.get('observation')
    if observation is None:
        if event['category']=='inner_realizability_loss': observation='Minimum covariance '+str(event['detail']['eigenvalue'])+'; full eigensystem/G/S/k/K/nuA in ledger.'
        elif event['category']=='pressure_interface_imbalance': observation='Weak/gauge balance small; full normal/tangential traction, k flux and strong divergence remain nonzero.'
        elif event['category']=='Jacobian_globalization_failure': observation='Original/fresh terminal searches:21 invalid evaluations each; no valid non-descent or linear failure.'
        else: observation='Pressure/gauge action accuracy fails fixed diagnostic reporting threshold.'
    scope = str(event.get('time', ''))+' '+event['scope']
    lines.append('| '+event['category']+' | '+scope+' | '+observation+' | '+event.get('causal_status', 'Unresolved; no physical qualification inferred.')+' |')
(here/'event-table.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(dict(ledger_sha256=sha(here/'first-failure-ledger.json'), bytes=(here/'first-failure-ledger.json').stat().st_size, events=len(events))))
