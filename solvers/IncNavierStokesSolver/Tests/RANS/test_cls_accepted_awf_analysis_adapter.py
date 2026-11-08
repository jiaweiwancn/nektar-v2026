#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Accepted-record identity and lossless summary checks; no native solves."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile

from build_cls_accepted_awf_analysis_suite import accepted_records,compact_diagnostic
from analyse_cls_accepted_awf import comparison_states
from extract_cls_accepted_awf_audits import observation
from run_cls_epawf_picard import require


def rejects(fn,word):
    try:fn()
    except ValueError as e:
        require(word in str(e),str(e));return
    raise AssertionError('expected rejection: '+word)


def main():
    rows=[dict(global_step=n,checkpoint_sha256='exact-'+str(n),
        attempt=dict(scheduled_complete=True,failure='')) for n in (1,2,3)]
    payload=lambda r:(''.join(json.dumps(x)+'\n' for x in r)).encode()
    selected=accepted_records(payload(rows),2)
    require(selected==rows[:2],'adapter altered source checkpoint identities or attempts')
    require(accepted_records(payload(rows[:2])+b'{"global_step":3',2)==rows[:2],
        'in-flight future append broke an immutable completed prefix')
    rejects(lambda:accepted_records(payload([rows[0],rows[0],rows[1]]),2),'duplicate')
    rejects(lambda:accepted_records(payload([rows[1]]),2),'contiguous')
    failed=deepcopy(rows[:2]);failed[1]['attempt']['failure']='hard rejection'
    rejects(lambda:accepted_records(payload(failed),2),'unaccepted')
    diagnostic=dict(state_sha256='actual-state',mean_wall_stress=.16,profiles=[dict(x=1.6,K=.5,
        functional_ids=[1,2],functional_phi=[.5,.5],profile=[[0.,1.],[.5,2.]])],
        independent_vertical_polynomial_minima=dict(scope='declared vertical anchors',state_sha256='actual-state',
            selected_inventory=dict(admissible=False,profile_branch_pattern='A'),
            anchors=[dict(anchor=[1.6,0.,0.],vertical_cubic_panels=[dict(minimum=-1e-6,z=.01)],
                profiles=[dict(rows=[dict(k=-1e-6,coordinate=[1.6,0.,.01],visited_by_selected_inventory=False)])])]))
    saved=deepcopy(diagnostic);summary=compact_diagnostic(diagnostic)
    require(saved==diagnostic,'compaction modified the underlying diagnostic')
    minimum=summary['independent_vertical_polynomial_minima_summary']
    require(minimum['minimum_vertical_panel']['minimum']==-1e-6 and
        minimum['minimum_profile_ray']['k']==-1e-6 and minimum['unregistered_profile_ray_rows']==1 and
        not minimum['selected_inventory']['admissible'],'summary hid a negative or unprotected witness')
    require(summary['profile_summaries'][0]['ground_row']==diagnostic['profiles'][0]['profile'][0],
        'summary changed a represented physical sample')
    with tempfile.TemporaryDirectory(prefix='accepted-audit-extraction-') as temporary:
        initial=dict(time=0.,common_samples=dict(k=[1.]))
        final=dict(time=.1,common_samples=dict(k=[.9]),common_reporting_queries_valid=False)
        invalid=Path(temporary)/'invalid-query.json'
        invalid.write_text(json.dumps(dict(time=.05,common_samples=dict(k=[-1.]),valid=False)))
        result=dict(initialization=initial,final=final,steps=160,
            snapshots={'80':dict(exit_code=0,path=str(invalid))})
        require(comparison_states(result)=={0:initial},
            'cross-suite comparison accepted an invalid final or snapshot query')
        result['initialization']=dict(initial,common_reporting_queries_valid=False)
        require(comparison_states(result)=={},'invalid initial query entered cross-suite norms')
        source=Path(temporary)/'fixture.json';source.write_text('{}\n')
        native=dict(state_audits=[dict(constitutive_modification=[1.,1.,0.])],audit_times=[.1],
            wall_policy='CLSCalibratedBudgetAWF',model_variant='Regularized',inner_stress_model='Linear',
            state=[1.],accepted_awf_final=dict(valid=True,dimensional_mass_dual_spatial_first4=[1.,2.,3.,4.]))
        extracted=observation(native,source,False)
        require(extracted['selected_operator']['outer_variant']=='Regularized' and
            extracted['selected_operator_alpha_zero_dimensional_blocks']==[1.,2.,3.,4.] and
            extracted['native_audit']['inner_stress'] is None and
            not extracted['inner_saturation_applicable'],'extractor substituted a faithful operator or invented inactive audit values')
        rejected=deepcopy(native);rejected.update(state_audits=[],initialization_complete=False,time=0.,
            physical_range_events=[dict(audit=dict(accepted_awf_covariance=[-1.,1.]))],
            accepted_awf_final=dict(valid=False,failure='preserved invalid query'))
        extracted=observation(rejected,source,False)
        require(extracted['state_classification']=='initialization_rejected_candidate' and
            extracted['native_audit']['accepted_awf_covariance']==[-1.,1.] and
            'selected_operator_alpha_zero_dimensional_blocks' not in extracted,
            'extractor relabelled an invalid candidate or its residual as accepted')
    print('analysis adapter: exact identities, contiguous accepted prefixes, quiet append, rejection and negative-witness retention passed')


if __name__=='__main__':main()
