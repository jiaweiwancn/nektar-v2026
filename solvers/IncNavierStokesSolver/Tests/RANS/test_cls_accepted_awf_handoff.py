#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Focused handoff provenance and native-history checks; no solver work."""
from copy import deepcopy
from pathlib import Path
import tempfile

from follow_cls_accepted_awf_pair import copy_exact, supervisor_command, validate_completed
from run_cls_accepted_awf import case_xml, model_choices
from run_cls_epawf_picard import load, require, sha, state_hash


def rejects(fn, word):
    try:
        fn()
    except ValueError as e:
        require(word in str(e), str(e));return
    raise AssertionError('expected rejection: '+word)


def main():
    with tempfile.TemporaryDirectory(prefix='accepted-handoff-') as temporary:
        root = Path(temporary)
        xml = case_xml('Accepted-AWF-05G', root/'native-initial', 0)
        definition = load(xml.parent/'definition.json')
        h0 = [1., 2., -0.]
        state = [1.1, 2.1, -0.]
        history = [[1.1, 2.1], [.9, 1.9]]
        attempts = []
        for n in range(1,161):
            attempts.append(dict(time=n*.000625, scheduled_complete=True,
                outer_sweeps=1,momentum_predictor_solves=1,pressure_system_solves=1,
                epsilon_transport_solves=1,k_transport_solves=1,bound_calls=1,
                additional_governing_equation_solves=0,mass_inverse_applications=2,
                pressure_gauge_size=2,additional_correction_or_restoration_solves=0,bound_kkt_solves=0,
                corrected_continuity_gauge=[0.,0.,0.],boundary_defect=0.,
                momentum_splitting_identity_defect=0.,imposed_trace_defect=0.,
                trace_mass_column_max=1.,retained_trace_column_max=1.,
                epsilon_frozen_linear_residual=0.,k_frozen_linear_residual=0.,
                k_bound=dict(activations=0,work_controls=dict(finite_qp_solver='RetainedPool'))))
        native = dict(temporal_coupling='PimpleSingleSweep',order=1,wall_policy='CLSCalibratedBudgetAWF',
            eta=0.,k_bound_mode='MassProjection',inner_stress_model='Linear',model_variant='Original',
            parameters=[1.5e-5,.5,.1,.5,.5],robust_parameters=[1.,0.,0.,0.,0.,0.,0.],
            implicit_callback_count=160,stage_audits=[[t['time']] for t in attempts],
            pimple_attempts=attempts,steps=160,time=.1,failure=None,state=state,
            native_history=history,native_history_times=[.1,.1],accepted_states=[h0])
        checkpoint = dict(order=1,steps=160,time=.1,dt=.000625,sha256='exact-binary-sha',
            coupled_state=state,history=history,history_times=[.1,.1])
        result = dict(status='numerical_trajectory_completed',steps=160,time=.1,failure=None,
            trajectory_receipt=dict(exit_code=0),definition=definition,
            accepted_checkpoint=dict(sha256='exact-binary-sha'),actual_h0_sha256=state_hash(h0),
            correction_integral=0.)
        proof = validate_completed(result,native,checkpoint,'Accepted-AWF-05G')
        require(proof['steps']==160 and proof['physical_steps_verified']==160,'incomplete handoff proof')
        broken=deepcopy(checkpoint);broken['history'][0][0]+=.01
        rejects(lambda:validate_completed(result,native,broken,'Accepted-AWF-05G'),'histories differ')
        broken=deepcopy(checkpoint);broken['coupled_state'][2]=0.
        rejects(lambda:validate_completed(result,native,broken,'Accepted-AWF-05G'),'coupled state differs')
        broken=deepcopy(result);broken['definition']['inventory']='ConstitutiveV3'
        rejects(lambda:validate_completed(broken,native,checkpoint,'Accepted-AWF-05G'),'sampling, inventory')
        broken=deepcopy(result);broken['status']='numerically_rejected'
        rejects(lambda:validate_completed(broken,native,checkpoint,'Accepted-AWF-05G'),'fixed horizon')
        companion=deepcopy(result)
        companion['definition'].update(companion='regularized',companion_reason='retained Original outer covariance exit',
            model_choices=model_choices('regularized'))
        regularized=deepcopy(native);regularized['model_variant']='Regularized'
        regularized['robust_parameters'][0]=.4;regularized['robust_parameters'][5]=1.
        validate_completed(companion,regularized,checkpoint,'Accepted-AWF-05G','regularized')
        rejects(lambda:validate_completed(companion,regularized,checkpoint,'Accepted-AWF-05G'),'declared unperturbed')
        rejects(lambda:validate_completed(result,native,checkpoint,'Accepted-AWF-05G','regularized'),'declared unperturbed')
        broken=deepcopy(companion);broken['definition']['companion_reason']=''
        rejects(lambda:validate_completed(broken,regularized,checkpoint,'Accepted-AWF-05G','regularized'),'declared unperturbed')
        saturated=deepcopy(companion)
        saturated['definition'].update(companion='regularized-saturated',
            companion_reason='retained Original outer and subsequent Linear inner covariance exits',
            model_choices=model_choices('regularized-saturated'))
        saturated_native=deepcopy(regularized);saturated_native['inner_stress_model']='SaturatedR2'
        validate_completed(saturated,saturated_native,checkpoint,'Accepted-AWF-05G','regularized-saturated')
        rejects(lambda:validate_completed(saturated,saturated_native,checkpoint,'Accepted-AWF-05G','regularized'),'declared unperturbed')
        rejects(lambda:validate_completed(companion,regularized,checkpoint,'Accepted-AWF-05G','regularized-saturated'),'declared unperturbed')
        broken=deepcopy(native);broken['pimple_attempts'][-1]['k_transport_solves']=2
        rejects(lambda:validate_completed(result,broken,checkpoint,'Accepted-AWF-05G'),'extra scheduled pass')
        source=root/'source.bin';source.write_bytes(b'checkpoint bytes\0\xff')
        target=root/'copies/State.bin';copy_exact(source,target);copy_exact(source,target)
        require(sha(source)==sha(target),'idempotent byte copy changed the native checkpoint')
        target.write_bytes(b'changed')
        rejects(lambda:copy_exact(source,target),'existing handoff input differs')
        command=supervisor_command(Path('/frozen/solver'),root/'selected','Accepted-AWF-05G',root/'inputs',[800,1600],7200)
        require('--resume' not in command and '--from-h0' not in command and
            command[command.index('--bound-max-total')+1]=='20000' and
            command[command.index('--prefix-report')+1]==str(root/'inputs/prefix.evolution.json'),
            'handoff altered restart semantics or controls')
        command=supervisor_command(Path('/frozen/solver'),root/'selected-regularized','Accepted-AWF-05G',
            root/'inputs',[800],7200,'regularized','retained faithful failure')
        require(command[command.index('--accepted-awf-companion')+1]=='regularized' and
            command[command.index('--accepted-awf-companion-reason')+1]=='retained faithful failure',
            'named companion identity/reason was lost in native continuation')
    print('accepted-AWF handoff: exact state/history, strict family, terminal rejection, immutable copies and controls passed')


if __name__ == '__main__':
    main()
