#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Focused guards for extending a completed native campaign, without a PDE."""
from copy import deepcopy
from pathlib import Path
import tempfile
from unittest.mock import patch

from extend_cls_awf_completed import validate_extension
from run_cls_accepted_awf import WORK, case_xml
from run_cls_epawf_picard import load, require


def main():
    with tempfile.TemporaryDirectory() as temporary:
        xml = case_xml('Accepted-AWF-05G', Path(temporary)/'case', 32,
            companion='regularized', companion_reason='preserved faithful failure')
        definition = load(xml.parent/'definition.json')
    saved = dict(steps=8000, time=5., order=1, dt=.000625, signature='actual-physical-signature',
        coupled=[1., -0., 3.], reference=[2., -0., 4.], defect=[0., -0., 0.],
        history=[[1., -0.], [2., 3.]], history_times=[5., 4.999375])
    native = dict(steps=8000, time=5., order=1, dt=.000625, failure='', signature=saved['signature'],
        state=saved['coupled'], fixed_reference=saved['reference'], fixed_defect=saved['defect'],
        native_history=saved['history'], native_history_times=saved['history_times'], audit_revision='NativeElementSideV2')
    progress = dict(status='targets_completed', failure=None, accepted_steps=8000, accepted_time=5., dt=.000625,
        wall_policy='CLSCalibratedBudgetAWF', case_family='AcceptedAWFV1', accepted_awf_case='Accepted-AWF-05G',
        case_definition=definition, accepted_awf_companion='regularized', accepted_awf_companion_reason='preserved faithful failure',
        accepted_awf_inner_refinement='none', projection_work_controls=WORK, inventory='AcceptedBudgetV1',
        bound_solver='RetainedPool', evaluation_reuse='ExactStepV1', audit_revision='NativeElementSideV2')
    # The existing full native-report checker is independently exercised by
    # the actual completed T=5 reports; isolate the orchestration guards here.
    with patch('extend_cls_awf_completed.check_report') as audit:
        validate_extension(progress, native, saved, [9600])
        require(audit.call_count == 1, 'full native report validation was omitted')
        mutations = [
            ('terminal numerical rejection', lambda p,n,s: p.update(status='numerically_rejected')),
            ('unfinished parent', lambda p,n,s: p.update(accepted_steps=7999)),
            ('changed timestep', lambda p,n,s: s.update(dt=.0003125)),
            ('changed physical signature', lambda p,n,s: n.update(signature='different')),
            ('changed history slot', lambda p,n,s: n.update(native_history=[[1., 0.], [2., 3.]])),
            ('changed reference', lambda p,n,s: n.update(fixed_reference=[2., 0., 4.])),
            ('changed model', lambda p,n,s: p['case_definition'].update(model_choices={'outer_variant':'Original'})),
            ('changed tail', lambda p,n,s: p['case_definition'].update(z_m=1.)),
            ('changed inventory', lambda p,n,s: p.update(inventory='ConstitutiveV3'))]
        for name, mutate in mutations:
            p,n,s = deepcopy(progress),deepcopy(native),deepcopy(saved); mutate(p,n,s)
            try:validate_extension(p,n,s,[9600])
            except ValueError:pass
            else:raise AssertionError('invalid extension allowed: '+name)
        try:validate_extension(progress,native,saved,[8000])
        except ValueError:pass
        else:raise AssertionError('old target allowed no advancement')
    print('completed extension: native endpoint/history/signature, model/tail/inventory, terminal rejection and target guards passed')


if __name__ == '__main__': main()
