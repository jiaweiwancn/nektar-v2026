#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Versioned classification of retained H9 scalar loads; no state re-evaluation."""
import argparse
from pathlib import Path
import numpy as np
from run_wj_h9 import load,save
from study_wj_h9 import identity
from wj_h9_weak import compare_anchor


def assess(root):
    source=root/'science-v2';out=root/'assessment-v1';out.mkdir();old=load(source/'summary.json');results=load(source/'weak-budgets.json');expected=load(root/'wj-h9-expected.json');checks={}
    for branch,value in results.items():
        native=expected['anchors'][branch]['independent_weak'];anchor=dict(components=native['native_components'],regional=native['regional_native'])
        checks[branch]=dict(anchor=compare_anchor(value['values']['A']['32'],anchor),phase_quadrature=value['checks'])
    assert all(v['anchor']['status']=='PASS'and all(c['status']=='PASS'for c in v['phase_quadrature'].values())for v in checks.values())
    summary=dict(old,status='PASS',original_summary=identity(source/'summary.json'),original_classification=old['status'],repaired_role='physical additive components versus diagnostic/aggregate buckets; no tolerance or load change',
        native_anchor_checks=checks,governing_defect='none demonstrated by these checks',causal_split='UNRESOLVED:no separate saved predictor/correction vectors')
    original=load(source/'admission.json');pre=original['projected_package_seconds']-original['projected_batch_seconds']-90
    # Original recovery projection conservatively included its own240s
    # reservation and excluded the15.338s failed prefix. Keep original bytes;
    # derive the corrected forecast without changing any actual charge.
    failed=load(root/'operations/science-batch/receipt.json')['seconds'];own=load(root/'operations/science-batch-recovery/request.json')['reserved_seconds']
    corrected=dict(original,correction='exclude only current recovery reservation; retain original failed science charge',prior_charged_seconds=pre-own+failed,
        projected_batch_seconds=original['projected_batch_seconds']+failed,projected_package_seconds=original['projected_package_seconds']-own+failed)
    save(out/'summary.json',summary);save(out/'native-anchor-checks.json',checks);save(out/'admission-correction.json',corrected)
    print('PASS scalar-record reassessment; no field/kernel/weak re-evaluation',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True,type=Path);a=p.parse_args();assess(a.root)
