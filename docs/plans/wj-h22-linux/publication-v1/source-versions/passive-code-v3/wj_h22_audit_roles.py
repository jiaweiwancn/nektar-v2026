# SPDX-License-Identifier: MIT
"""Full-audit parity uses an audited state, never an empty producer map."""
from wj_h17_observer import equal

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
