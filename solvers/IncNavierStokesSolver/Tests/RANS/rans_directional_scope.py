# SPDX-License-Identifier: MIT
"""Versioned required scope of the bounded directional diagnostic package."""
from rans_duration_contract import require

SCHEMA='DirectionalResolutionRequiredScopeV1'
SHARED=('directional_capability','isotropic_compatibility','frozen_scalar_reconciliation',
        'geometry','history_cache_CLS','required_short_verification')
A95=('A95_neutral_control','A95_response_fields','A95_response_native_weak',
     'physical_budget_reconstruction','physical_budget_accuracy')
REQUIRED=(*SHARED,*A95,'A93_mesh_candidate')


def classify(receipt):
    require(receipt['schema']==SCHEMA,'directional_required_scope_version')
    require(receipt['numerical_operator']=='GeneratedVolumeOwnerV1','directional_old_operator')
    require(receipt['response_max_time']==.016 and receipt['A93_response_out_of_scope'],
            'directional_response_scope_changed')
    require(receipt['layout']==dict(fields=[6,5],pressure=[4,3]),'directional_layout_scope_changed')
    require(receipt['unknown_limit']==3000,'directional_unknown_guard_changed')
    gates={key:receipt['gates'].get(key,dict(status='not_attempted')) for key in REQUIRED}
    for key,row in gates.items():
        require(row['status'] in ('pass','fail','failed_or_unresolved','unresolved','not_attempted','blocked'),
                'invalid_directional_gate_status:'+key)
        require(not (row['status']=='pass' and row.get('primitive_status','pass')!='pass'),
                'directional_pass_label_disagrees_with_primitives:'+key)
    # A93 is a separate mesh question. Its failure cannot suppress independent
    # A95 work, while every shared prerequisite remains mandatory for A95.
    a95_ready=all(gates[key]['status']=='pass' for key in (*SHARED,*A95))
    missing=[key for key in REQUIRED if gates[key]['status']!='pass']
    return dict(schema=SCHEMA,status='pass' if not missing else 'failed_unresolved_or_incomplete',
        gates=gates,unmet_required_gates=missing,A95_future_prerequisites_satisfied=a95_ready,
        response_extension_eligible=False,
        reason='this diagnostic package is capped at .016; the historical B5/B11 accuracy decision and new-layout neutral duration are not superseded',
        historical_status_unchanged=receipt.get('historical_status_unchanged',True),
        scope='qualified directional capability, isotropic compatibility and two independent numerical hypotheses; no longer response or A93 perturbation')
