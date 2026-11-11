#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact verified reductions of the refined qualification follow-up artifacts."""
import argparse
from datetime import datetime,timezone,timedelta
import math
from pathlib import Path
from collections import defaultdict
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import load,sha,require


def receipt(path):return dict(path=str(path),sha256=sha(path),bytes=path.stat().st_size)


def compact_case(row):
    require(sha(Path(row['receipt']['output']))==row['receipt']['output_sha256'],'immutable native readback changed')
    raw=load(Path(row['receipt']['output']));forensic=raw['diagnostics']['selected_state']
    p=row['physical'];keys=('time','H','h','z_K','z_m','reference','mean_wall_stress','native_degrees_of_freedom',
        'scaled_native_spatial_blocks','dimensional_mass_dual_spatial_first4','common_normalized_spatial_first4',
        'representation_L2_velocity_gradient_production_referenceProduction','representation_L2_k_gradient',
        'fixed_bands','matching_trace_defect','epsilon_reporting_lower','common_reporting_queries_valid')
    result={k:row[k] for k in ('selection','receipt','state_sha256','reference_sha256','defect_sha256','native_history_sha256',
        'native_history_times_sha256','state_equals_fixed_reference','state_reference_numeric_differences','state_reference_bit_differences')}
    result.update(physical={k:p[k] for k in keys},dimensional_reference_errors=row['dimensional_reference_errors'],
        normalized_reference_errors=row['common_normalized_reference_errors'],complete_audit=row['complete_audit'],
        named_interface=row['named_interface_metrics'],selected_Original_blocks=row['selected_Original_blocks'],
        analytical_stress_difference=row['analytical_stress_difference'],analytical_production_difference=row['analytical_production_difference'],
        ground_reaction=row['ground_reaction'],transport_decomposition=row['transport_decomposition'],
        native_constitutive_difference_locations=row['native_outer_constitutive_difference_locations'],last_physical_step=row.get('step'),
        last_step_ground_storage=row.get('last_step_ground_storage'))
    result['signed_mass_dual_balance']={b:{k:{key:value for key,value in data.items() if key not in ('element_residual_localization',)}
        for k,data in row['balance'][b].items()} for b in map(str,range(4))}
    minima={}
    for key,law in [('minimum_inner','selected'),('minimum_outer','selected'),('minimum_same_state_Original_outer','same_state_Original')]:
        w=forensic[key];minima[key]={k:w[k] for k in ('coordinate','element','geometry','inventory','epsilon_semantics')}
        minima[key]['minimum_eigenvalue']=w[law]['eigenvalues'][0]
    result['same_state_covariance_minima']=minima
    transport=raw['diagnostics']['neutral_reference_balance']['transport_decomposition']['retained_state'][2]
    error=max(abs(i+2*d/3) for i,d in zip(transport['physical_isotropic_production_load'],transport['dilatation_load']))
    require(error<1e-13,'actual isotropic production/dilatation load identity differs')
    item=row['transport_decomposition']['retained_state'][2]
    result['explicit_k_dilatation_part']=dict(load_identity_max_error=error,
        mass_dual_norm=(5/3)*item['dimensional_mass_dual_norms']['dilatation_load'],
        signed_physical_integral=(5/3)*item['integrated_q_div_u'],
        definition='R_advection contains +k div(u); actual isotropic production is -(2/3) k div(u), entering residual with a minus sign. Their explicit residual part is +(5/3) k div(u). Remaining strain/production contributions are retained; this is not the entire effect of divergence.')
    panels=defaultdict(lambda:[0.,0.,0.])
    for e in row['element_values']:
        v=e['values'];z=(e['lower'][2],e['upper'][2]);panels[z][0]+=v[4];panels[z][1]+=v[6];panels[z][2]+=v[2]
    result['native_panel_gradient_divergence']=[dict(z_bounds=z,dUdz_reference_L2=math.sqrt(v[0]),dkdz_reference_L2=math.sqrt(v[1]),strong_divergence_L2=math.sqrt(v[2])) for z,v in sorted(panels.items())]
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','output'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--response-analysis',type=Path);p.add_argument('--status-capture',type=Path)
    p.add_argument('--retained',type=Path,help='verified qualify_cls_awf_neutral output (including a pure algebra re-reduction)')
    a=p.parse_args();r=a.root;retained=a.retained or r/'retained-diagnosis';d=load(retained/'summary.json');require(d['status']=='passed','incomplete retained diagnosis')
    manifest=load(r/'runtime-manifest.json')
    out=dict(captured_Asia_Shanghai=datetime.now(timezone(timedelta(hours=8))).isoformat(),root=str(r),
        script_sha256=sha(Path(__file__)),runtime_manifest=receipt(r/'runtime-manifest.json'),
        solver_sha256=manifest['solver']['sha256'],probe_sha256=manifest['probe']['sha256'],
        source_base=manifest['base_revision'],cases={k:compact_case(v) for k,v in d['cases'].items()},comparisons=d['comparisons'])
    window=load(r/'P1-all-field-confirmation/summary.json')
    out['P1_all_field_confirmation']={k:v for k,v in window.items() if k!='observations'}
    out['P1_all_field_confirmation']['receipt']=receipt(r/'P1-all-field-confirmation/summary.json')
    files={'terminal_observer':'observer-terminal-validation-v2/progress.json','zero_time':'first-cell-zero-time/summary.json',
        'seed_contracts':'late-response/validation/summary.json','seed_side_telemetry':'seed-side-report-validation-v2/summary.json',
        'old_runtime_equivalence':'old-runtime-equivalence/summary.json','final_state_observer_contract':'terminal-final-state-contract/summary.json'}
    out['checks']={}
    for key,file in files.items():
        q=load(r/file);out['checks'][key]=dict(receipt=receipt(r/file),data=q)
        if key=='seed_contracts':out['checks'][key]['data']['seed_identity'].pop('actual_increment',None)
        if key=='final_state_observer_contract':out['checks'][key]['data']={k:v for k,v in q.items() if k!='row'}
    if a.status_capture:out['P2_status']=load(a.status_capture);out['P2_status_receipt']=receipt(a.status_capture)
    for extra in ('P2-T2-diagnosis','P2-T2.5-diagnosis','P2-terminal-diagnosis'):
        file=r/extra/'summary.json'
        if file.exists():
            q=load(file)
            if q['status']=='passed':out.setdefault('later_pressure_comparisons',{})[extra]=dict(cases={k:compact_case(v) for k,v in q['cases'].items()},comparisons=q['comparisons'],receipt=receipt(file))
    if a.response_analysis:
        q=load(a.response_analysis/'summary.json');require(q['status']=='passed','incomplete response analysis')
        out['response']={k:v for k,v in q.items() if k!='observations'};out['response']['receipt']=receipt(a.response_analysis/'summary.json')
        out['response']['observations']={}
        for leg,rows in q['observations'].items():
            out['response']['observations'][leg]={}
            for t,row in rows.items():
                item={k:row[k] for k in ('selection','receipt','state_sha256','reference_sha256','native_history_sha256','named_interface_metrics','step','complete_audit') if k in row}
                item['physical']={k:row['physical'][k] for k in ('time','H','h','z_K','z_m','reference','mean_wall_stress',
                    'native_degrees_of_freedom','dimensional_mass_dual_spatial_first4','common_normalized_spatial_first4',
                    'fixed_bands','matching_trace_defect','epsilon_reporting_lower','common_reporting_queries_valid')}
                out['response']['observations'][leg][t]=item
    out['response_campaign']={}
    for leg in ('control','perturbed'):
        q=load(r/'late-response'/leg/'progress.json');segments=q['segments'];usage=[list(map(float,x['receipt']['usage'].split())) for x in segments]
        out['response_campaign'][leg]={k:q[k] for k in ('status','accepted_steps','elapsed_time','physical_time','checkpoint','checkpoint_sha256','dt','physical_t0') if k in q}
        out['response_campaign'][leg].update(end_to_end_seconds=q.get('seconds'),complete_block_seconds=sum(v['receipt']['seconds'] for v in segments),
            native_process_wall_seconds=sum(v[0] for v in usage),
            process_user_seconds=sum(v[1] for v in usage),peak_rss_KiB=max((v[3] for v in usage),default=0),
            bound_activations=sum(v['bound_activations'] for v in segments),correction_integral=sum(v['correction_integral'] for v in segments),
            governing_counts={k:sum(x['governing_counts'][k] for x in segments) for k in ('momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves','bound_calls')},
            additional_governing_solves=sum(v['additional_governing_solves'] for v in segments))
    index=[]
    for file in (r/'runtime-manifest.json',r/'retained-selections.json',retained/'summary.json',r/'observer-upgrade-receipt.json',r/'observer-upgrade-final-receipt.json',r/'late-response/summary.json',r/'late-response/declared-policy.json'):
        index.append(receipt(file))
    for item in d['cases'].values():
        index.append(dict(role='immutable native readback',**receipt(Path(item['receipt']['output']))));index.append(dict(role='selected native checkpoint',**receipt(Path(item['selection']['checkpoint'])/'State.bin')))
    out['artifact_index']=index;a.output.parent.mkdir(parents=True,exist_ok=True);atomic(a.output,out)


if __name__=='__main__':main()
