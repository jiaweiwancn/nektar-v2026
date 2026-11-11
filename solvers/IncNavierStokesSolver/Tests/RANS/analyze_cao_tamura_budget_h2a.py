#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Strict, zero-advance budget readbacks of completed H2a startup checkpoints."""
import argparse
from datetime import datetime,timezone
import math
import os
from pathlib import Path
import platform

from continue_cls_epawf_pimple import atomic,freeze
from diagnose_cls_accepted_covariance import checkpoint,common_packet,run_probe
from run_cls_epawf_picard import THREADS,load,require,sha,state_hash


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('probe','startup','h1','output'):p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();require(platform.system()=='Linux','authoritative Linux readbacks required')
    a.output.mkdir(parents=True,exist_ok=False);os.environ.update(THREADS)
    probe,identity=freeze(a.probe.resolve(),a.output/'runtime')
    progress=load(a.startup/'progress.json');roots=[a.startup]
    if progress.get('reused_temporal_root'):roots.append(Path(progress['reused_temporal_root']))
    reports={}
    for root in roots:
        for receipt in root.glob('*/receipt.json'):
            report=receipt.parent/'case.evolution.json'
            if not report.exists():continue
            q=load(receipt)
            if q['code']==0 and not load(report)['failure']:
                require(q['report_sha256']==sha(report) and q['xml_sha256']==sha(receipt.parent/'case.xml'),
                    'completed startup report/session differs from receipt')
                require(receipt.parent.name not in reports,'duplicate startup identity')
                reports[receipt.parent.name]=(report,q)
    require(reports,'no completed native startup states')
    out=dict(schema=1,status='started',created_UTC=datetime.now(timezone.utc).isoformat(),
        startup=str(a.startup),startup_progress_sha256=sha(a.startup/'progress.json'),probe=identity,
        zero_Lift=True,governing_solves=0,physical_commits=0,observations={},costs={},
        integration='native physical consistent weights; independent source panels split at native/profile/AWF breakpoints',
        budget_units=dict(k_integral='m^4/s^2 per unit span',k_rates='m^4/s^3 per unit span',epsilon_rates='m^4/s^4 per unit span'),
        pressure='absolute selected-traction datum; no independent regional or common mean removal')
    atomic(a.output/'summary.json',out)
    for name,(report,receipt) in sorted(reports.items(),key=lambda v:(load(v[1][0])['time'],v[0])):
        q=load(report);cp=Path(q['checkpoint']);saved=checkpoint(cp);before=sha(cp/'State.bin')
        # Internal 32-step chunks remain indexed but need no expensive extra
        # readback: retain seeds, exact temporal endpoints and primary milestones.
        selected=name.endswith('-h0') or q['time']<=.00040000000001 or any(
            abs(q['time']-t)<1e-13 for t in (.0008,.004,.04))
        attempts=q['pimple_attempts'];usage=receipt['usage'].split()
        out['costs'][name]=dict(driver_wall=receipt['seconds'],process_wall=float(usage[0]),
            user_CPU=float(usage[1]),system_CPU=float(usage[2]),peak_RSS_KiB=int(usage[3]),
            local_steps=len(attempts),report=str(report),report_sha256=sha(report),
            mass_factor_hits=sum(v['velocity_mass_factor_reused'] for v in attempts),
            persisted_complete_packets=sum(v['attempt_packet_selected'] for v in attempts))
        if not selected:continue
        packet=dict(common_packet(cp),operation='finite_inflow_h1',state=saved['coupled'],history=saved['coupled'],
            time=saved['time'],alpha=0.,inner_stress_model='Linear',observation_x=0.,
            native_report_file=str(report),native_report_sha256=sha(report),
            permitted_queries_file=str(a.h1/'input/permitted-queries.txt'),
            permitted_queries_sha256=sha(a.h1/'input/permitted-queries.txt'))
        d,rr=run_probe(probe,report.parent/'case.xml',packet,a.output/name)
        require(d['governing_solves']==d['physical_commits']==0 and d['initialization_work']==[0,0]
            and not d['native_history_advanced'] and d['selected_state_unchanged'],'readback advanced or modified selected state')
        require(d['state_sha256']==state_hash(saved['coupled']) and sha(cp/'State.bin')==before,
            'readback changed immutable checkpoint/state')
        fits={}
        for quantity in ('U','minusRxz','sigmaU'):
            rows=[v for v in d['permitted_flat_native_queries'] if v['quantity']==quantity]
            weight=math.fsum(v['weight'] for v in rows);valid=math.fsum(v['weight'] for v in rows if v['valid'])
            fits[quantity]=dict(original_weight=weight,valid_weight=valid,invalid_queries=[v for v in rows if not v['valid']],
                weighted_RMSE=math.sqrt(math.fsum(v['weight']*(v['native_prediction']-v['target'])**2 for v in rows)/weight) if valid==weight else None)
        samples=d['physical_samples'];weights=samples['wall_weights'];measure=math.fsum(weights)
        integrals={str(b):math.fsum(w*v for w,v in zip(samples['weights_'+str(b)],samples['values_'+str(b)])) for b in (2,3)}
        last=attempts[-1] if attempts else {}
        audit=d['complete_audit'];out['observations'][name]=dict(time=q['time'],tstar=q['time']/.04,
            dt=q['dt'],global_steps=q['steps'],checkpoint=str(cp),checkpoint_sha256=before,receipt=rr,
            state_sha256=d['state_sha256'],reference_sha256=d['reference_sha256'],
            mean_wall_traction=math.fsum(w*v for w,v in zip(weights,samples['wall_stress']))/measure,
            wall_K_min_max=[min(samples['wall_K']),max(samples['wall_K'])],physical_scalar_integrals=integrals,
            complete_audit=audit,native_turbulence_budget=d['native_turbulence_budget'],
            boundary=d['developing_boundary'],selected_spatial_blocks=d['selected_spatial_block_norms'],
            physical_reference_projection_errors=d['physical_projection_metrics'],
            profiles=d['physical_profiles'],profile_columns=d['physical_profile_columns'],
            permitted_flat_native_comparison=fits,provider=d['provider'],
            correction_calls=sum(v['bound_calls'] for v in attempts),
            correction_activations=sum(v['bound_activations'] for v in attempts),
            correction_signed=sum(v['k_bound']['integral_total'] for v in attempts),
            correction_positive=sum(v['k_bound']['integrals_inner_outer_near_positive_negative'][3] for v in attempts),
            correction_negative=sum(v['k_bound']['integrals_inner_outer_near_positive_negative'][4] for v in attempts),
            step_diagnostics={k:last.get(k) for k in ('frozen_momentum_defect','momentum_splitting_identity_defect',
                'corrected_continuity_gauge','imposed_trace_defect','refreshed_trace_defect','k_frozen_linear_residual',
                'trace_mass_column_max','trace_physical_storage_contribution_max','solved_epsilon_retained_bitwise')},
            native_stage_schedule=[{k:v[k] for k in ('time','momentum_predictor_solves','pressure_system_solves',
                'epsilon_transport_solves','k_transport_solves','bound_calls','additional_governing_equation_solves')} for v in attempts])
        atomic(a.output/'summary.json',out)
    out.update(status='passed',completed_UTC=datetime.now(timezone.utc).isoformat())
    atomic(a.output/'summary.json',out);print(a.output/'summary.json')


if __name__=='__main__':main()
