#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact H14 publication; heavy native coefficients and rows stay in place."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import tarfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_wj_h4 import load,save,sha
from wj_h14_execution import file_id,validate_pair,OwnerIndex,compact_open_packet
from study_wj_h14 import H13,stored_reference
from wj_h4_points import evaluate,vector_hash
from wj3_measure import checkpoint_payload
from wj_h7_analysis import context

def table(path,rows):
    if not rows:return
    with path.open('w',newline='')as f:
        w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)

def main(root,source):
    out=root/'publication-v1';out.mkdir();(out/'controls').mkdir();(out/'operations').mkdir()
    def copy(p,name=None):
        p=Path(p);dest=out/(name or p.name);shutil.copyfile(p,dest);assert sha(dest)==sha(p);return file_id(dest)
    book=load(root/'accounting.json');copy(root/'accounting.json','accounting-at-publication.json')
    for p in(root/'target-admission-v1.json',root/'target-admission-v2.json',root/'target-admission-v3.json',
        root/'science/target-resource-design-v1.json',root/'science/remaining-projection-v2.json',
        root/'science/target-source-allocation-v1.json',root/'science/qualification-gates.json',
        root/'science/runtime-build-v4.json',root/'science/runtime-qualified-v4.json',
        root/'science/immutable-cache-equivalence.json',root/'science/completed-assessor-equivalence.json',
        root/'science/small/dense-parity-summary.json',root/'science/small/restart-rollback-checks.json',
        root/'science/small/unsupported/result.json'):
        if p.exists():copy(p,'controls/'+('unsupported.json'if p.parent.name=='unsupported'else p.name))
    parents=[]
    m=load(H13/'publication-v1/manifest.json')
    for p in('paired-case/flat/case.xml','paired-case/paired-domain.json','paired-case/preservation-protocol.json',
             'paired-case/mapped-backend-prerequisite.json','component-gates.json'):
        target=H13/'publication-v1'/p;rec=next(v for v in m['artifacts']if v['path']==str(target));assert file_id(target)==rec;parents.append(rec)
    heavy=[];bindings=[];history=[];q1=[];interfaces=[];attempts=[];factors=[];views=[]
    cases=[root/'science/target/cp0']+sorted((root/'science/target').glob('pilot[1-4]'))
    final=None
    for case in cases:
        rp=case/'case.evolution.json'
        if not rp.exists():continue
        r=load(rp);report_identity=file_id(rp);heavy.append(report_identity)
        cp=checkpoint_payload(r['checkpoint'],include_reference=True)
        assert vector_hash(cp['state'])==vector_hash(r['state'])and cp['signature']==r['signature']and cp['steps']==r['steps']
        state_file_identity=file_id(Path(r['checkpoint'])/'State.bin');heavy.append(state_file_identity)
        for k in('fixed_reference','fixed_defect'):assert vector_hash(cp[k])==vector_hash(r[k])
        bindings.append(dict(case=case.name,producer=report_identity,xml=file_id(case/'case.xml'),checkpoint=state_file_identity,
            checkpoint_absolute_step=cp['steps'],checkpoint_time=cp['time'],dt=cp['dt'],coupled_state_sha256=vector_hash(cp['state']),
            fixed_reference_sha256=vector_hash(cp['fixed_reference']),fixed_defect_sha256=vector_hash(cp['fixed_defect']),
            physical_signature_sha256=hashlib.sha256(cp['signature'].encode()).hexdigest(),native_history_sha256=[vector_hash(v)for v in cp['native_history']],
            native_history_times=cp['native_history_times'],native_initialization_work=r['initialization_lift_calls_iterations'],initializer_calls=r['initial_velocity_calls'],
            role='fresh startup'if r['steps']==0 else'new completed state; restored entry remains alpha0'))
        copy(case/'execution-diff.json','controls/'+case.name+'-execution-diff.json')
        assessment=case/('assessment-v2'if case.name=='cp0'else'assessment')/'result.json'
        if assessment.exists():
            a=load(assessment);copy(assessment,case.name+'-assessment.json')
            for v in a['states']:
                s=v['summary'];history.append(dict(case=case.name,absolute_step=v['absolute_step'],time=v['time'],status=s['status'],
                    native_state_sha256=v['state_sha256'],min_principal=s['min_in_plane'],min_transverse=s['min_transverse'],min_covariance=s['min_covariance'],min_scalar=s['min_scalar'],
                    divergence_max=s['raw_divergence_max'],divergence_scaled_max=s['divergence_scaled_max'],E_G=s['errors']['E_G'],E_epsilon=s['errors']['E_epsilon'],
                    E_U=s['errors']['E_U'],E_k=s['errors']['E_k'],U_reference_L2=s['reference_error']['full']['U']['L2'],
                    U_CP0_drift_L2=s['initial_drift']['full']['U']['L2'],G_reference_L2=s['reference_error']['full']['G']['L2'],G_CP0_drift_L2=s['initial_drift']['full']['G']['L2'],
                    witness_x=s['witness']['physical'][0],witness_z_mm=40*s['witness']['physical'][2],witness_geometry=s['witness']['owner']['geometry'],witness_normal_side=s['witness']['owner']['normal_side']))
                for region,d in v['all_regional_Q1_moments'].items():
                    for j,value in enumerate(d['moments']):q1.append(dict(case=case.name,absolute_step=v['absolute_step'],time=v['time'],region=region,native_pressure_id=j,signed_moment=value))
                for row in s['owner_interfaces']:
                    interfaces.append(dict(case=case.name,absolute_step=v['absolute_step'],time=v['time'],x=row['physical'][0],z_mm=40*row['physical'][2],
                        inner_geometry=row['inner_owner']['geometry'],outer_geometry=row['outer_owner']['geometry'],x_side=row['outer_owner']['x_side'],
                        K=row['K'],rc_mm=40*row['rc'],Lambda=row['Lambda'],shared_U_max=max(abs(x)for x in row['shared_U']),shared_k=row['shared_k'],epsilon=row['epsilon'],epsilon_target=row['epsilon_target'],
                        epsilon_trace_defect=row['epsilon_trace_defect'],traction_x=row['total_traction_jump'][0],traction_z=row['total_traction_jump'][2],k_flux_inner_minus_outer=row['k_flux_jump'],pi_inner_minus_outer=row['reduced_pressure_jump']))
        for p in r['pimple_attempts']:
            snapshot=None
            if p.get('k_snapshot'):
                snapshot=load(p['k_snapshot']);assert sha(p['k_snapshot'])==p['k_snapshot_sha256'];heavy.extend([file_id(p['k_snapshot']),file_id(snapshot['matrix_file'])])
                assert sha(snapshot['matrix_file'])==snapshot['matrix_sha256']
            attempts.append(dict(case=case.name,time=p['time'],alpha=p['alpha'],scheduled_complete=p['scheduled_complete'],failure=p['failure'],
                momentum=p['momentum_predictor_solves'],pressure=p['pressure_system_solves'],k=p['k_transport_solves'],epsilon=p['epsilon_transport_solves'],bound_calls=p['bound_calls'],bound_activations=p['bound_activations'],
                raw_corrected_k_bits=vector_hash(snapshot['raw'])==vector_hash(snapshot['corrected'])if snapshot else None,velocity_size=p.get('velocity_size'),pressure_size=p.get('pressure_gauge_size'),flow_size=p.get('assembled_flow_size'),
                k_size=p.get('k_solve_size'),epsilon_size=p.get('epsilon_solve_size'),retained_BG_columns=p.get('retained_BG_columns'),true_schur_residual=p.get('true_schur_reapplication_residual'),
                predictor_residual=p.get('predictor_linear_residual'),pressure_residual=p.get('pressure_linear_residual'),k_residual=p.get('k_frozen_linear_residual'),epsilon_residual=p.get('epsilon_frozen_linear_residual'),
                splitting_identity=p.get('momentum_splitting_identity_defect'),boundary_defect=p.get('boundary_defect'),retained_epsilon_trace_column_max=p.get('retained_trace_column_max'),epsilon_trace_mass_column_max=p.get('trace_mass_column_max')))
            for name in('momentum','mass','k','epsilon'):
                if name+'_factor_allocation'in p:factors.append(dict(case=case.name,time=p['time'],factor=name,allocation=p[name+'_factor_allocation']))
        for allocation in r['dense_evolution_capacity'].get('band_factor_allocations',[]):
            factors.append(dict(case=case.name,time=r['time'],factor='native constructor/projection/norm factor',allocation=allocation))
        vp=case/'view/result.json'
        if vp.exists():copy(vp,case.name+'-strict.json');views.append(file_id(vp));heavy.append(file_id(case/'view/native.json'))
        if not r['failure']:final=(case,r)
    if not final:raise ValueError('no retained target state')
    save(out/'state-bindings.json',dict(parent_H13=parents,states=bindings,source_python_commit=source,native_compiled_source='ec5fb45d (effective H14 C++ tree; exact compiled file identities in runtime-build-v4)',
        scope='producer reports remain execution reports, matrix files remain matrix roles; checkpoint-file and contained-state hashes are distinct'))
    table(out/'timeline-and-errors.csv',history);table(out/'complete-Q1-moments.csv',q1);table(out/'interface-history.csv',interfaces);table(out/'attempt-counts.csv',attempts);save(out/'factor-allocations.json',factors)
    if(root/'science/target/pilot-result.json').exists():pilot=load(root/'science/target/pilot-result.json');copy(root/'science/target/pilot-result.json','pilot-result.json')
    else:pilot=dict(status='UNRESOLVED',completed_steps=0,actual_attempts=None)
    for p in sorted((root/'science/target').glob('pilot-projection-*.json')):copy(p)
    statuses=dict(implementation_build='PASS',small_dense_and_bound_qualification='PASS',fresh_target_CP0='PASS',target_startup_quality='PASS',
        target_strict_CP0='PASS',pilot=pilot['status'],pilot_completed_steps=pilot['completed_steps'],pilot_actual_attempts=pilot['actual_attempts'],
        further_pilot_steps=pilot.get('remaining_steps_status','NOT_RUN'),duration_preservation='NOT_RUN',hill_evolution='NOT_RUN',experimental_validation='NOT_RUN',raised_height_cases='EXCLUDED_NOT_PREPARED',
        temporal_order='NOT_ESTABLISHED',scalable_production_capacity='NOT_QUALIFIED')
    save(out/'component-gates.json',statuses)
    # Retain actual logs/receipts in compact compressed form. Native heavy JSON
    # packets/matrices/checkpoints are referenced above, never duplicated here.
    receipts=[];timing={}
    for folder in sorted((root/'operations').iterdir()):
        if not(folder/'receipt.json').exists():continue
        dest=out/'operations'/folder.name;dest.mkdir()
        for name in('request.json','receipt.json'):
            if(folder/name).exists():shutil.copyfile(folder/name,dest/name)
        p=folder/'command.log'
        with p.open('rb')as f,gzip.open(dest/'command.log.gz','wb')as g:shutil.copyfileobj(f,g)
        receipts.append(dict(operation=folder.name,raw_log=file_id(p),compressed=file_id(dest/'command.log.gz')))
        if(folder/'timing.jsonl').exists():
            for line in(folder/'timing.jsonl').read_text().splitlines():
                v=json.loads(line);key=(folder.name,v['name']);d=timing.setdefault(key,dict(operation=key[0],scope=key[1],calls=0,wall_inclusive=0.,wall_exclusive=0.))
                d['calls']+=1;d['wall_inclusive']+=v['wall_seconds'];d['wall_exclusive']+=v['exclusive_wall_seconds']
    save(out/'compressed-log-bindings.json',receipts);table(out/'native-timing-scopes.csv',list(timing.values()))
    # Spatially meaningful profiles on the actual original physical coordinates.
    c=stored_reference();profile=[]
    heights=sorted(set([c.d0,c.d0+c.H,.175,.25,.25607517888672315,.31946961248954786,.375,c.top]+np.linspace(c.d0,c.d0+c.H,21).tolist()+
        (c.d0+np.geomspace(c.H,c.top-c.d0,55)).tolist()))
    for role,r in [('CP0',load(root/'science/target/cp0/case.evolution.json')),('terminal',final[1])]:
        index=OwnerIndex(r)
        for x in(-10.,0.,15.):
            for z in heights:
                if z>c.top:continue
                p=[x,0.,z]
                for o in index(p):
                    v=evaluate(r,r['state'],p,o);rr=z-c.d0;ue=c.U_inner(rr)if o['region']=='inner'else c.U_outer(rr)
                    profile.append(dict(role=role,step=r['steps'],time=r['time'],x=x,z_mm=40*z,geometry=o['geometry'],normal_side=o['normal_side'],x_side=o['x_side'],region=o['region'],
                        U=v['U'][0],W=v['U'][2],Gxz=v['G'][0][2],divergence=v['raw_divergence'],k=v['k'],epsilon=v['epsilon'],pi=v['pi'],mechanical_pressure=v['mechanical_pressure'],
                        Rxx=v['R'][0][0],Rzz=v['R'][2][2],Rxz=v['R'][0][2],min_covariance=min(v['covariance']),principal=v['principal']['in_plane']if v['principal']else None,
                        K=v['sampled_K'],U_exact=ue))
    table(out/'physical-profiles.csv',profile)
    fig,axes=plt.subplots(1,3,figsize=(11,3.7))
    for role,color in [('CP0','#3568a8'),('terminal','#c45835')]:
        rows=[v for v in profile if v['role']==role and v['x']==15.]
        axes[0].plot([v['U']for v in rows],[v['z_mm']for v in rows],label=role,color=color)
        low=[v for v in rows if v['z_mm']<=30]
        axes[1].plot([v['Gxz']for v in low],[v['z_mm']for v in low],label=role,color=color)
    ref=[v for v in profile if v['role']=='CP0'and v['x']==15.]
    axes[0].plot([v['U_exact']for v in ref],[v['z_mm']for v in ref],'k--',label='continuum reference',linewidth=.9)
    for a in axes[:2]:a.set_ylabel('Physical height above original floor (mm)');a.axhline(9.5,color='.5',linestyle=':',linewidth=.8);a.grid(alpha=.2)
    axes[0].set_xlabel(r'$U_x/U_\infty$');axes[0].set_ylim(3.5,800);axes[0].legend(fontsize=8)
    axes[1].set_xlabel(r'$\partial U_x/\partial z$ (code units)');axes[1].set_ylim(3.5,30);axes[1].axhspan(3.5,9.5,color='.8',alpha=.3)
    axes[1].axhline(7,color='#8b4264',linestyle='--',linewidth=.8)
    axes[2].plot([v['time']for v in history],[v['min_principal']for v in history],'o-',color='#3568a8')
    axes[2].axhline(.002068763923454,color='#b23b3b',linestyle='--',label='half-reference quality')
    axes[2].axhline(.004137527846908,color='.4',linestyle=':',label='exact neutral')
    axes[2].set_xlabel(r'$t U_\infty/H_h$');axes[2].set_ylabel(r'Finite-inventory $\min(m/\nu_{eff})$');axes[2].grid(alpha=.2);axes[2].legend(fontsize=7)
    fig.suptitle('H14: unchanged6mm flat target; outlet profiles and bounded startup/pilot')
    fig.tight_layout();fig.savefig(out/'target-profiles-and-margin.png',dpi=160);plt.close(fig)
    table(out/'operation-accounting.csv',[dict(operation=k,phase=v['phase'],status=v['status'],charged_or_reserved_seconds=v.get('seconds',v['reserved_seconds']),peak_tree_RSS_bytes=v.get('peak_tree_rss_bytes'),reason=v.get('reason'))for k,v in book['entries'].items()])
    totals=dict(setup_charged_seconds=60,operations_charged_at_publication=sum(v.get('seconds',v['reserved_seconds'])for v in book['entries'].values()),
        closure_reserved_charged_seconds=300,allowance=7200,unknown_external_precheck_charge_seconds=5,
        peak_successful_native_RSS_bytes=max((v.get('peak_tree_rss_bytes',0)for k,v in book['entries'].items()if v['status']=='PASS'and k in('target-full-domain-CP0-v1','target-cost-gated-pilot-v1','target-CP0-strict-view-v1')),default=0),
        failed_consumer_RSS_bytes=book['entries']['target-strict-gate-and-pilot-projection-v4']['peak_tree_rss_bytes'],
        source_python_commit=source,build_flags='-ffp-contract=off -fno-fast-math; RANS_BUILD_SOLVER=ON; Release GNU11.4.0',
        rank=1,numerical_threads=1,compile_parallel=8,hard_server_worker_limit=72)
    save(out/'resources-at-publication.json',totals)
    duration=505.02466521278404;dt=final[1]['dt'];steps=math.ceil(duration/dt)if False else int(np.ceil(duration/dt))
    measured=[v for v in timing.values()if v['operation']=='target-cost-gated-pilot-v1'and v['scope']=='native.step']
    warm=sum(v['wall_inclusive']for v in measured)/sum(v['calls']for v in measured)if measured else None
    save(out/'preservation-forecast.json',dict(reference_clocks=dict(bulk=26.13350407707765,slowest_outer=57.17628191545675,
        inner_diffusion=195.80590488759879,top_K_over_epsilon=505.02466521278404),proposed_horizon=duration,
        terminal_window=[447.8483832973273,505.02466521278404],pilot_dt=dt,steps_to_cover_horizon=steps,
        measured_inclusive_native_step_seconds=warm,conditional_serial_seconds=warm*steps if warm else None,
        scope='same measured rate extrapolation only; cold startup, required monitoring, temporal/spatial adequacy and long-time admissibility remain separate; duration NOT_RUN'))
    # Verify all selected current heavy bytes now; distinguish remote recorded
    # verification from the subsequently independently transported compact bytes.
    unique={v['path']:v for v in heavy}
    for v in unique.values():assert file_id(v['path'])==v
    save(out/'heavy-native-manifest.json',dict(scope='recorded remote-heavy verification; no claim the native payload bytes were transported to Mac',files=list(unique.values())))
    artifacts=[file_id(p)for p in sorted(out.rglob('*'))if p.is_file()]
    save(out/'manifest.json',dict(schema='WJH14CompactV1',artifacts=artifacts))
    save(out/'seal.json',dict(manifest=file_id(out/'manifest.json'),parent_H13_manifest=file_id(H13/'publication-v1/manifest.json'),
        publication_accounting_snapshot=file_id(out/'accounting-at-publication.json'),final_closure='separate subsequent immutable snapshot, never relabel this publication snapshot'))
    with tarfile.open(root/'publication-v1.tar.gz','w:gz')as t:t.add(out,arcname='publication-v1')
    save(root/'publication-transport.json',file_id(root/'publication-v1.tar.gz'))
    print('publication complete',statuses,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--source-commit',required=True);a=p.parse_args();main(a.root,a.source_commit)
