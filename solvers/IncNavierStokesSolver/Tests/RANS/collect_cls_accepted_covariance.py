#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact immutable receipts for covariance diagnostics and live controls.

This is an evidence collector, not a solver, restart or acceptance mechanism.
Dense rows, matrices, field vectors and full logs remain on the run volume.
"""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import platform

from continue_cls_epawf_pimple import atomic,header,process_identity
from run_cls_epawf_picard import load,require,sha,state_hash


def witness(p):
    keys=('coordinate','geometry','element','inner','inventory','gradient_element_region',
        'constitutive_region','gradient_region_matches_law','k','epsilon','epsilon_retained',
        'epsilon_semantics','endpoint_K','cutoff','budget_ratio','velocity_gradient','divergence',
        'epsilon_H_retained','epsilon_H_refreshed_projected','epsilon_H_current_unprojected')
    out={k:p[k] for k in keys if k in p}
    for name in ('selected','same_state_Original','same_state_Regularized'):
        law=p.get(name)
        if law is not None:out[name]={k:v for k,v in law.items() if k!='terms'}
    return out


def physical_metrics(q):
    keys=('time','H','h','z_K','z_m','Z','K_sampling','reference','common_reporting_knots',
        'common_reporting_queries_valid','dimensional_mass_dual_spatial_first4',
        'common_normalized_spatial_first4','scaled_native_spatial_blocks','common_scales_U_K_E_H_V',
        'representation_L2_velocity_gradient_production_referenceProduction','representation_L2_k_gradient',
        'mean_wall_stress','matching_trace_defect','native_degrees_of_freedom','state_unchanged','governing_solves')
    out={k:q[k] for k in keys if k in q}
    out['fixed_bands']={**q['fixed_bands'],'regions':[b for b in q['fixed_bands']['regions'] if b['x_region']=='all']}
    out['full_profiles_and_spatial_regions']='retained in indexed native readback; no field vectors in this compact record'
    return out


def transport_attribution(root,operators):
    import numpy as np
    from scipy.linalg import cho_factor,cho_solve
    first=load(root/'operators-v4/Regularized-0/readback/result.json')['diagnostics']
    rows=np.array(first['native_free_k_mass']['rows'],dtype=int)
    factor=cho_factor(np.array(first['native_free_k_mass']['matrix']),lower=True)
    result=dict(scope='consistent-mass dual vector correlations, not signed physical-band integrals',
        postprocessing_mass_applications=0,governing_solves=0,cases={})
    for name,c in operators['cases'].items():
        if not name.startswith('Regularized-'):continue
        d=load(Path(c['full_free_loads_retained']))['diagnostics'];indices=rows-int(d['selected_state']['offsets'][2])
        w=d['physical']['weak_k_component_loads'];r=np.array(w['total_constrained'])[indices];a=np.array(w['conservative_advection'])[indices]
        inverse=cho_solve(factor,np.column_stack((r,a)));rr=float(r@inverse[:,0]);aa=float(a@inverse[:,1]);ra=float(r@inverse[:,1])
        residual_minus_advection=float(max(0.,rr+aa-2*ra)**.5)
        result['cases'][name]=dict(time=c['time'],residual_norm=rr**.5,advection_norm=aa**.5,
            residual_minus_advection_norm=residual_minus_advection,
            mass_dual_cosine=ra/(rr*aa)**.5 if rr*aa>0 else None)
        result['postprocessing_mass_applications']+=2
    return result


def plot_horizons(data,path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
    colors={'baseline':'C0','half-dt':'C1','pressure-p2':'C2'}
    for name,c in data['horizon_comparison']['cases'].items():
        rows=sorted(c['observations'].values(),key=lambda v:v['time']);color=colors[name]
        axes[0,0].plot([v['time'] for v in rows],[v['covariance_minimum_outer']['selected']['eigenvalues'][0] for v in rows],
            marker='o',color=color,label=name)
        point=[(v['time'],w['velocity_gradient'][0][0]) for v in rows for w in v['fixed_witnesses'] if int(w['geometry'])==75]
        axes[0,1].plot([x for x,y in point],[y for x,y in point],marker='o',color=color,label=name)
        axes[1,1].plot([v['time'] for v in rows],[100*(v['physical']['mean_wall_stress']/.16-1) for v in rows],marker='o',color=color,label=name)
    rejected=[('baseline',data['archived_stages']['refined']['attempted_time'],
        data['archived_stages']['refined']['stages']['corrected']['outer']['selected']['eigenvalues'][0])]
    for name,key in [('half-dt','half-dt'),('pressure-p2','outer-pressure-p2-v2')]:
        c=data['controls'][key]
        if c.get('rejected_attempt') and c.get('archived_stage_diagnostic'):
            rejected.append((name,c['rejected_attempt']['time'],c['archived_stage_diagnostic']['stages']['corrected']['outer']['selected']['eigenvalues'][0]))
    for name,t,value in rejected:axes[0,0].plot(t,value,'x',ms=9,mew=2,color=colors[name])
    axes[0,0].axhline(0,color='k',lw=.6);axes[0,0].set_ylabel('Outer covariance minimum');axes[0,0].set_title('Accepted checkpoints; crosses: rejected attempts')
    axes[0,1].set_ylabel('G_xx on left element');axes[0,1].set_title('Fixed x=-3, z=3.68804, geometry 75')
    for name,rows in data['horizon_comparison']['same_time_differences'].items():
        selected=sorted((r for r in rows.values() if r['time']>0),key=lambda r:r['time'])
        for field,style in [('u','-'),('pressure','--')]:
            axes[1,0].semilogy([r['time'] for r in selected],[r['common_normalized_difference'][field] for r in selected],
                style,marker='o',color=colors[name],label=name+' '+field)
    axes[1,0].set_ylabel('Common-normalized L2 difference from baseline');axes[1,1].set_ylabel('Mean traction reference error (%)')
    for ax in axes.flat:ax.set_xlabel('Physical time');ax.legend(fontsize=8)
    fig.savefig(path,dpi=160);plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--campaign-status',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--implementation-commit',required=True);a=p.parse_args();r=a.root.resolve()
    require(platform.system()=='Linux' and r.parent==Path('/media/jay/新加卷/nektar-v2026-runs'),
        'authoritative isolated Linux evidence required')
    a.output.mkdir(exist_ok=False);index={}
    def record(path):
        path=Path(path);require(path.is_file(),'missing evidence '+str(path))
        item=dict(path=str(path.resolve()),sha256=sha(path),bytes=path.stat().st_size);index[str(path.resolve())]=item;return item
    def read(path):path=Path(path);record(path);return load(path)
    snapshots=a.output/'observed-live-snapshots';snapshots.mkdir()
    def read_live(path,label):
        saved=snapshots/(label+'.json');saved.write_bytes(Path(path).read_bytes())
        item=record(saved);item['observed_live_path']=str(path)
        return json.loads(saved.read_text(),parse_int=float)
    out=dict(schema=1,captured_UTC=datetime.now(timezone.utc).isoformat(),implementation_commit=a.implementation_commit,
        base_commit='310ebe427f575f74c9ef1ea240ad5afa1acc45f7',root=str(r),
        scope='zero-PDE evidence collection; instantaneous diagnostics, time-integrated estimates and exact step corrections remain distinct',
        diagnostic_governing_solves=0,physical_commits=0,collector_script_sha256=sha(Path(__file__)),
        archived_stages={},tests={},controls={})
    stages=read(r/'stages-v4/summary.json');require(stages['status']=='passed','archived validation did not pass')
    for name,c in stages['cases'].items():
        checked=read(Path(c['receipt']['output']));full=checked['diagnostics']
        row=out['archived_stages'][name]={k:c[k] for k in ('checkpoint_sha256','archived_packet','archived_packet_sha256','stage_hashes','history_sha256','receipt')}
        row.update(accepted_time=checked['native_replay_context']['accepted_time'],attempted_time=checked['physical_time'])
        row['stages']={s:dict(inner_minimum=d['minimum_inner']['selected']['eigenvalues'][0],
            inner_coordinate=d['minimum_inner']['coordinate'],outer=witness(d['minimum_outer']),
            one_sided_witnesses=[witness(x) for x in d['one_sided_witnesses']],
            scalar_budget=d['scalar_budget']) for s,d in full.items()}
        row['independent_tensor_checks']=c['independent_tensor_checks']
        row['contract_rejections']={k:dict(exit_code=v['receipt']['exit_code'],log_sha256=v['receipt']['log_sha256'],
            receipt=v['receipt'],reason=v['reason']) for k,v in c['contract_rejections'].items()}
        packet=read(c['archived_packet']);record(packet['matrix_file'])
        record(r/'stages-v4'/name/'checkpoint/State.bin')
    tangent=read(r/'tangent-general-v3/native.json');require(tangent['status']=='passed','tangent oracle failed')
    keep=lambda name:any(s in name for s in ('refined-committed-audit','refined-post_flow-audit','refined-analytical-neutral',
        'coarse-post_flow-audit','one-post_flow-audit','one-post_flow-physical'))
    out['tangent']=dict(scope=tangent['scope'],queries=len(tangent['queries']),native_records=len(tangent['tangent_results']),
        maximum_independent_work_difference=max(x['work_difference'] for x in tangent['tangent_results']),
        results=[x for x in tangent['tangent_results'] if keep(x['id'])],
        gradient_attribution={name:c.get('gradient_attribution') for name,c in tangent['cases'].items()},
        flow_preserves_k_and_epsilon_bits={name:c['flow_preserves_k_and_epsilon_bits'] for name,c in tangent['cases'].items()},
        executable_sha256=tangent['tangent_executable_sha256'],receipt=tangent['receipt'])
    if tangent.get('figure'):record(tangent['figure']['path'])
    operators=read(r/'operators-v4/summary.json');require(operators['status']=='passed','strict operator observations failed')
    out['operators']=dict(scope='selected uncancelled operator versus same-state Original; retained epsilon trace; alpha=0',
        windows=operators['windows'],cases={})
    for name,c in operators['cases'].items():
        d={k:v for k,v in c.items() if k not in ('selected_covariance','Original_covariance')}
        d['physical']=physical_metrics(c['physical'])
        op=d['operators'];left=op['selected_continuity_gauge_loads'];right=op['counterfactual_continuity_gauge_loads']
        d['operators']={k:v for k,v in op.items() if not k.endswith('continuity_gauge_loads')}
        d['operators'].update(constraint_rows=len(left),maximum_continuity_gauge_load=max(map(abs,left)),
            continuity_gauge_exactly_equal=left==right,common_gauge_load=left[-1])
        d['selected_covariance']=witness(c['selected_covariance']);d['Original_covariance']=witness(c['Original_covariance'])
        out['operators']['cases'][name]=d;record(c['full_free_loads_retained']);record(Path(c['checkpoint_copy'])/'State.bin')
    out['component_mass_dual']=read(r/'operators-v4/component-mass-dual.json')
    out['transport_attribution']=transport_attribution(r,operators)
    out['milestone_audits']=read(r/'regularized-milestone-audits.json')
    comparison=r/'horizon-comparison-final/summary.json'
    if comparison.exists():
        data=read(comparison);out['horizon_comparison']={k:v for k,v in data.items() if k!='cases'}
        out['horizon_comparison']['cases']={}
        for name,c in data['cases'].items():
            rows={}
            for step,o in c['observations'].items():
                row={k:v for k,v in o.items() if k not in ('covariance_minimum_inner','covariance_minimum_outer','fixed_witnesses')}
                row['physical']=physical_metrics(o['physical'])
                row.update(covariance_minimum_inner=witness(o['covariance_minimum_inner']),
                    covariance_minimum_outer=witness(o['covariance_minimum_outer']),
                    fixed_witnesses=[witness(w) for w in o['fixed_witnesses']])
                rows[step]=row;record(o['receipt']['output'])
            out['horizon_comparison']['cases'][name]=dict(dt=c['dt'],unavailable_times=c['unavailable_times'],observations=rows)
    for name in ('native-regression-v4/summary.json','stages-explicit-contract-v2/summary.json','python-source-validation.json'):
        data=read(r/name);out['tests'][name]={k:v for k,v in data.items() if k not in ('solver','cases','files')}
    for name in ('ctest-focused.log','python-focused-checks.log','observer-final-test.log',
        'tangent-parallel-regression.log','build-solver-tangent.log','build-validation-final.log','build-focused-checks.log'):
        record(r/name)
    for name in ('runtime-v4/identity.json','source-final-verified.json','campaign-library-verification.json',
        'campaign-library-verification-delivery.json'):
        if (r/name).exists():out[name]=read(r/name)
    status=read(a.campaign_status);out['campaigns']=status
    for c in status['cases'].values():
        if c['status']!='running':continue
        progress=read(c['immutable_progress'])
        recent=[s['seconds']/(s['end_step']-s['start_step']) for s in progress['segments'][-4:]
            if not s['failure'] and s['end_step']>s['start_step']]
        if recent:
            cost=sum(recent)/len(recent);c['recent_wall_seconds_per_accepted_step']=cost
            c['conditional_hours_to_next_target']=(c['next_target']-c['accepted_steps'])*cost/3600 if c['next_target'] else None
            c['conditional_hours_to_T5']=(round(5/c['dt'])-c['accepted_steps'])*cost/3600
    observer=r/'window-observer-v2/summary.json'
    if observer.exists():
        data=read_live(observer,'window-observer')['cases']['Accepted-AWF-2G'];selected={}
        for key,w in data['fixed_physical_windows'].items():
            first,last=map(int,key.split('-'));steps=last-first
            if steps in (160,320,800,1600) and first%steps==0:selected[key]=w
        out['fixed_windows']=dict(captured_frontier=[data['accepted_steps'],data['accepted_time']],
            wall_time_scale=data['wall_time_scale'],initial_flow_through_time=data['initial_flow_through_time'],
            windows=selected,rule='Frozen production-report common quadrature; strict augmented-panel checkpoint changes are in operators.windows. Sampled source amounts are trapezoidal estimates, exact repairs come from the step ledger.')
    for name in ('half-dt','outer-pressure-p2-v2'):
        base=r/name;controller=read_live(base/'status.json',name+'-controller');path=base/'native/Accepted-AWF-2G/trajectory'
        files=list(path.glob('checkpoint*/State.bin'));require(files,'control has no accepted checkpoint')
        cp=max(files,key=lambda x:x.stat().st_mtime);native=header(cp.parent,True)
        row=out['controls'][name]=dict(controller=controller,accepted_steps=native['steps'],accepted_time=native['time'],
            checkpoint=record(cp),state_sha256=state_hash(native['coupled_state']),
            native_history_sha256=[state_hash(v) for v in native['history']],
            native_history_times_sha256=state_hash(native['history_times']))
        for part in ('initial','trajectory'):
            folder=base/'native/Accepted-AWF-2G'/part;record(folder/'case.xml');record(folder/'definition.json')
        row['definition']=read(base/'native/Accepted-AWF-2G/initial/definition.json')
        resultfile=base/'native/Accepted-AWF-2G/result.json'
        if resultfile.exists():
            result=read(resultfile);row['result']={k:v for k,v in result.items() if k not in ('initialization','final','snapshots','definition')}
            report=read(path/'case.evolution.json');attempts=report['pimple_attempts'];accepted={x[0] for x in report['stage_audits']}
            require(native['steps']==report['steps'] and native['time']==report['time'] and
                state_hash(native['coupled_state'])==state_hash(report['state']), 'control accepted state/time differs from checkpoint')
            require([state_hash(x) for x in native['history']]==[state_hash(x) for x in report['native_history']] and
                state_hash(native['history_times'])==state_hash(report['native_history_times']), 'control history rollback differs from checkpoint')
            row['accepted_state_time_all_native_history_bits_verified']=True
            counts=('outer_sweeps','momentum_predictor_solves','pressure_system_solves',
                'epsilon_transport_solves','k_transport_solves','bound_calls')
            require(all(tuple(t[k] for k in counts)==(1,1,1,1,1,1) and
                t['additional_governing_equation_solves']==0 for t in attempts),'control repeated a governing solve/sweep')
            row['native_covariance_trajectory']=[dict(time=t,covariance=d.get('accepted_awf_covariance'),
                inner_witness=d.get('accepted_awf_inner_minimum_point'),outer_witness=d.get('accepted_awf_outer_minimum_point'),
                strong_continuity=d.get('strong_continuity')) for t,d in zip(report['audit_times'],report['state_audits'])]
            row['accepted_bound_activations']=sum(t.get('bound_activations',0) for t in attempts if t['time'] in accepted)
            row['governing_solve_counts']={k:sum(t.get(k,0) for t in attempts) for k in ('outer_sweeps','momentum_predictor_solves',
                'pressure_system_solves','epsilon_transport_solves','k_transport_solves','bound_calls','additional_governing_equation_solves')}
            row['rejected_attempt']={k:v for k,v in attempts[-1].items() if not isinstance(v,(list,dict))} if report['failure'] else None
            for t in attempts[-1:] if report['failure'] else []:
                packet=read(t['k_snapshot']);record(packet['matrix_file'])
        child=controller.get('child');row['driver_identity_now']=process_identity(child['pid']) if child else None
        staged=r/('stages-half-dt' if name=='half-dt' else 'stages-pressure-p2')/'summary.json'
        if staged.exists():
            data=read(staged);require(data['status']=='passed','new control stage diagnostic did not pass')
            case=next(iter(data['cases'].values()));full=read(case['receipt']['output'])['diagnostics']
            values=read(case['archived_packet']);offset=[int(x) for x in full['committed']['offsets']]
            unchanged={str(b):state_hash(values['committed'][offset[b]:offset[b+1]])==
                state_hash(values['post_flow'][offset[b]:offset[b+1]]) for b in (2,3)}
            row['archived_stage_diagnostic']=dict(receipt=case['receipt'],independent_tensor_checks=case['independent_tensor_checks'],
                flow_preserves_k_epsilon_coefficients=unchanged,
                stages={n:dict(outer=witness(d['minimum_outer']),inner_minimum=d['minimum_inner']['selected']['eigenvalues'][0],
                    one_sided_witnesses=[witness(w) for w in d['one_sided_witnesses']],scalar_budget=d['scalar_budget']) for n,d in full.items()})
        process_file=path/'process.json'
        if process_file.exists():
            packet=read(process_file);current=process_identity(packet['pid'])
            actual=[]
            if current:
                actual=(Path('/proc')/str(int(packet['pid']))/'cmdline').read_bytes().rstrip(b'\0').decode().split('\0')
            row.update(native_identity_now=current,native_command_matches=bool(current and
                actual==packet['command'] and current['cwd']==str(path)))
            image=record(packet['command'][0]);require(image['sha256']==controller['solver_sha256'],'control executable changed')
            row['native_executable']=image
        manifest=base/'native/manifest.json'
        if manifest.exists():row['frozen_runtime_identity']=read_live(manifest,name+'-manifest')['solver']
    if out.get('horizon_comparison'):
        figure=a.output/'failure-horizons.png';plot_horizons(out,figure);out['failure_horizon_figure']=record(figure)
    atomic(a.output/'artifact-index.json',dict(schema=1,captured_UTC=out['captured_UTC'],root=str(r),artifacts=index))
    atomic(a.output/'native.json',out);print(a.output/'native.json')


if __name__=='__main__':main()
