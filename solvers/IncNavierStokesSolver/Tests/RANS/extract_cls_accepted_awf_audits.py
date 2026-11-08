#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Extract actual stress activation/interface audits; no operator substitution.

Alpha-zero residuals describe each selected case operator. Native same-state
stress-work/production differences do not constitute a full faithful-operator
residual. A rejected candidate remains distinct from the last accepted state.
"""
import argparse
from datetime import datetime,timezone
from pathlib import Path

from analyse_cls_accepted_awf import report
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import require,sha,state_hash

KEYS=('inner_stress','inner_stress_max_ratio_point','inner_constitutive_modification',
      'constitutive_modification','interface','strong_continuity',
      'accepted_awf_scalar_domain','accepted_awf_covariance','accepted_awf_sources',
      'accepted_awf_hard_admissible')
COLUMNS=dict(
    inner_stress=['active_physical_volume','positive_k_volume','minimum_a_over_c',
        'maximum_a_over_c','zero_k_volume','minimum_turbulent_radial_tangent'],
    inner_stress_max_ratio_point=['x','y','z','local_k','representative_K','nu_tA',
        'strain_norm','a_over_c','t','radial_tangent'],
    inner_constitutive_modification=['stress_work_difference_from_same_state_Linear',
        'production_difference_from_same_state_Linear','work_minus_production'],
    constitutive_modification=['stress_work_difference_from_same_state_Original_outer',
        'production_difference_from_same_state_Original_outer','work_minus_production'],
    interface=['physical_measure','integrated_pi_minus_p','normal_traction_jump_RMS',
        'full_traction_jump_RMS','normal_traction_jump_max','k_diffusive_flux_jump_RMS',
        'net_normal_flow','integrated_normal_traction_jump'],
    strong_continuity=['inner_divergence_L2','outer_divergence_L2',
        'inner_divergence_integral','outer_divergence_integral','net_boundary_flux'],
    accepted_awf_covariance=['inner_minimum_eigenvalue','outer_minimum_eigenvalue'],
    accepted_awf_sources=['inner_production_integral','E_C_integral','redistribution_integral',
        'molecular_work_integral','production_minus_E_C'])


def observation(native,source,first):
    audits=native.get('state_audits',[]);i=0 if first else -1
    if audits:
        audit=audits[i];time=native['audit_times'][i]
        kind='accepted_native_state'
    else:
        events=native.get('physical_range_events',[])
        audit=events[-1]['audit'] if events else {}
        time=native['time'];kind='initialization_rejected_candidate'
    if native.get('initialization_complete') is False:kind='initialization_rejected_candidate'
    selected=native.get('accepted_awf_initial' if first else 'accepted_awf_final') or {}
    values={key:audit.get(key) for key in KEYS}
    inner=values['inner_stress'];volume=inner[1]+inner[4] if inner else None
    result=dict(time=time,state_classification=kind,native_audit=values,
        inner_saturation_applicable=native['inner_stress_model']=='SaturatedR2',
        evaluated_inner_volume=volume,active_inner_volume_fraction=inner[0]/volume if volume else None,
        selected_operator=dict(wall_policy=native['wall_policy'],outer_variant=native['model_variant'],
            inner_stress=native['inner_stress_model'],chi=None,destruction_cap=None,redistribution='zero'),
        source_native_report=dict(path=str(source.resolve()),sha256=sha(source)),
        alpha_zero_residual_interpretation='selected case operator at retained fields; not a counterfactual faithful Original/Linear operator')
    state=(native.get('accepted_states') or [native.get('state')])[0] if first else native.get('state')
    if state is not None:result['coupled_state_sha256']=state_hash(state)
    if selected.get('valid') is not False and selected.get('common_reporting_queries_valid') is not False:
        result['selected_operator_alpha_zero_dimensional_blocks']=selected.get('dimensional_mass_dual_spatial_first4')
        result['selected_operator_alpha_zero_native_scaled_blocks']=selected.get('scaled_native_spatial_blocks')
    else:result['optional_alpha_zero_diagnostic_invalid']=selected.get('failure','invalid common physical query')
    return result


def mapping(text):
    label,files=text.split('=',1);initial,final=files.split(',',1)
    return label,Path(initial).resolve(),Path(final).resolve()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case',type=mapping,action='append',required=True,metavar='LABEL=INITIAL_NATIVE_REPORT,FINAL_NATIVE_REPORT')
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    require(not a.output.exists(),'use a new immutable evidence output')
    require(len({label for label,_,_ in a.case})==len(a.case),'duplicate case labels')
    out=dict(schema=1,captured_UTC=datetime.now(timezone.utc).isoformat(),production_governing_solves_added=0,
        script_sha256=sha(Path(__file__)),native_vector_columns=COLUMNS,cases={},
        caveat='Native stress-work counterfactuals use identical fields; selected-operator residuals are not faithful-operator residuals for companions. None means inapplicable or absent, never an inferred zero.')
    for label,initial,final in a.case:
        left,right=report(initial),report(final)
        require(left['wall_policy']==right['wall_policy']=='CLSCalibratedBudgetAWF','wrong wall formulation')
        row=dict(initial=observation(left,initial,True),final=observation(right,final,False),
            final_native_failure=right.get('failure'))
        if right.get('failure') and right.get('physical_range_events'):
            audit=right['physical_range_events'][-1]['audit']
            row['rejected_candidate']=dict(time=right['pimple_attempts'][-1]['time'] if right.get('pimple_attempts') else right['time'],
                native_audit={key:audit.get(key) for key in KEYS},
                coupled_state_sha256=state_hash(right['rejected_state']) if right.get('rejected_state') else None,
                state_classification='rejected_candidate; no alpha-zero residual inferred from accepted fields')
        out['cases'][label]=row
    a.output.parent.mkdir(parents=True,exist_ok=True);atomic(a.output,out);print(a.output)


if __name__=='__main__':main()
