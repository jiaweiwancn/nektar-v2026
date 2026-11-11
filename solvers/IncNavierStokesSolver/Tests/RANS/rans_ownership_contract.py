# SPDX-License-Identifier: MIT
"""Corrective-package scope; physical-budget reconstruction never implies accuracy."""
import math
from rans_duration_contract import require, identity, MEMBERS
from classify_rans_physical_budget_accuracy import classify_row

OPERATOR='GeneratedVolumeOwnerV1'
SCHEMA='NativeOwnershipRequiredScopeV1'


def budget_gate(analysis,accuracy,binding,expected):
    """Recompute scientific primitives after callers validate sealed file identities."""
    require(binding['schema']=='NativeOwnershipBudgetBindingV1','budget_binding_version')
    identity(binding['configuration'],expected['configuration'])
    identity(analysis['configuration'],expected['configuration']);identity(accuracy['configuration'],expected['configuration'])
    require(binding['numerical_operator']==OPERATOR,'old_operator_budget_history')
    require(binding['analysis_sha256']==expected['analysis_sha256'] and
            binding['accuracy_sha256']==expected['accuracy_sha256'],'stale_budget_accuracy_receipt')
    require(binding['analysis_sha256']==accuracy['source_sha256'],'foreign_budget_accuracy_source')
    require(binding['cohort']==expected['cohort'] and binding['amplitude']==2e-5,'foreign_budget_response')
    require(all(abs(float(v)-expected['time'])<=1e-13 for v in (analysis['time'],accuracy['time'],binding['time'])),
            'budget_horizon_mismatch')
    require(set(binding['members'])==set(MEMBERS),'budget_companion_missing')
    require(binding['controls']==dict(volume_orders=[48,64],physical_face_orders=[16,32]),'budget_quadrature_controls_missing')
    require(all(binding['members'][m]==expected['members'][m] for m in MEMBERS),'foreign_budget_native_reports')
    require(analysis['execution_status']=='completed' and analysis['reconstruction_status']=='pass','physical_reconstruction_failed')
    required={(kind,equation) for kind in ('baseline','perturbed') for equation in ('inner_k','outer_k','outer_epsilon','whole_patch_k')}
    rows=[classify_row(row) for row in analysis['finite_defects']]
    require(len(rows)==len(required) and {(v['kind'],v['equation']) for v in rows}==required,'physical_budget_accuracy_coverage')
    require(analysis['quadrature'],'budget_quadrature_missing')
    quadrature=all(math.isfinite(q['limit']) and q['limit']>=0 and
        all(math.isfinite(v) and v>=0 for v in q['changes'].values()) and
        max(q['changes'].values())<=q['limit'] for q in analysis['quadrature'])
    passed=quadrature and all(v['spatial_status'] in ('pass','unresolved_floor') and v['temporal_status']=='pass' for v in rows)
    require(accuracy['checks']==rows,'budget_primitive_recomputation_mismatch')
    require(accuracy['status']==('pass' if passed else 'failed_or_unresolved'),'pass_label_disagrees_with_budget_primitives')
    return dict(status='pass' if passed else 'failed_or_unresolved',reconstruction_status='pass',
        accuracy_status='pass' if passed else 'failed_or_unresolved',quadrature_status='pass' if quadrature else 'unresolved',checks=rows)


def classify_scope(receipt):
    require(receipt['schema']==SCHEMA,'corrective_scope_version')
    require(receipt['numerical_operator']==OPERATOR,'old_operator_evidence_not_qualification')
    required=['ownership_geometry','history_cache_CLS','natural_short_matrix','selected_prescribed',
              'A93_neutral','A95_neutral','A95_response_fields','A95_response_native_weak',
              'physical_budget_reconstruction','physical_budget_accuracy']
    gates={k:receipt['gates'].get(k,dict(status='not_attempted')) for k in required}
    require(receipt['out_of_scope']['A93_response'] and receipt['response_max_time']==.016,'corrective_response_scope_changed')
    # This receipt reports independent failures, including when A93 fails but A95 proceeds.
    return dict(schema=SCHEMA,numerical_operator=OPERATOR,
        status='pass' if all(v['status']=='pass' for v in gates.values()) else 'failed_unresolved_or_incomplete',
        gates=gates,response_extension_eligible=all(gates[k]['status']=='pass' for k in (
            'ownership_geometry','history_cache_CLS','natural_short_matrix','selected_prescribed',
            'A95_neutral','A95_response_fields','A95_response_native_weak','physical_budget_reconstruction','physical_budget_accuracy')),
        scope='corrected bounded patch only; independent accuracy gates retained; no A93 response or extension beyond .016 in this package')
