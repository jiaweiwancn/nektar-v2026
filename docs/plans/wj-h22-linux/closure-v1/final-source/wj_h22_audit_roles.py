# SPDX-License-Identifier: MIT
"""Full-audit parity uses an audited state, never an empty producer map."""
from wj_h17_observer import equal
from run_wj_h4 import load
from wj_h15_observer import verify

def select_auxiliary_audits(inputs):
    """Resolve H18 B's audit through H20's sealed selection, retaining A2's role."""
    index=load(verify(inputs['prior_indices']['20']['publication']))
    pin=next(v for v in index['files']if v['path'].endswith('/controls/inputs--selected.json'))
    original=load(verify(pin))['states']['M']['2']
    for key in ('audit_report','qualified_commit'):verify(original[key])
    return {'M':{'2':dict(audit_report=original['audit_report'],original_verification=original['qualified_commit'],
        role='separate H18 B startup audit; A2 remains the completed producer')}}

def native_full_map_parity(current,producer,separate_audit=None):
    assert current['native_full_audit_status']=='PASS'
    if producer['native_full_audit_status']=='PASS':
        audited=producer;role='genuine producing report full audit'
    else:
        if separate_audit is None:
            raise ValueError('producer audit was NOT_SCHEDULED; a separately pinned audited-state report is required')
        assert producer['native_full_audit_status']=='NOT_SCHEDULED'and producer['native_hard_audit']=={}
        assert separate_audit['native_full_audit_status']=='PASS'
        for name in ('state','fixed_reference','fixed_defect','dt','time','signature','original_qh0_state_sha256'):
            equal(producer[name],separate_audit[name],'separate audit physical '+name)
        audited=separate_audit;role='separate startup/re-audit of bitwise identical physical state'
    equal(current['native_hard_audit'],audited['native_hard_audit'],'full numerical audit with correct source role')
    return dict(status='PASS',map_source_role=role,producer_callback_role='unchanged genuine producing completion',
                history_or_alpha_replaced=False,original_producer_full_audit_status=producer['native_full_audit_status'])
