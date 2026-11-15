#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""CT12 compact primitive decisions; collection success is not admission."""
import argparse,copy,hashlib,json,os,platform,subprocess,sys
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
import numpy as np
import scipy
from ct_accepted9_pressure_run import sealed
from ct_accepted_feasibility import write_json,require
from ct_accepted_budget_gate import physical_budgets,finite_numbers
from rans_request_manifest import RequestManifest,sha


def read(p):return json.loads(Path(p).read_text())


def pressure_gate(q,require_quadrature=True):
    reasons=[]
    if q.get('schema')!='CT12TargetPressureV1' or set(q.get('members',{}))!={'0.001','0.0005'}:reasons.append('coverage')
    for dt,m in q.get('members',{}).items():
        times=[round(v['time'],12) for v in m.get('steps',[])]
        expected=[round(j*float(dt),12) for j in range(1,round(.002/float(dt))+1)]
        if times!=expected:reasons.append('missing_time:'+dt)
        for step in m.get('steps',[]):
            if not finite_numbers(step):reasons.append('nonfinite:'+dt)
            if step.get('C_maximum')!=0 or step.get('unused_gauge') is not False:reasons.append('pressure_level:'+dt)
            if any(v>2e-10 for v in step.get('block_reconstruction',{}).values()):reasons.append('block_reconstruction:'+dt)
            if set(step.get('action_reconstruction',{}))!={'predictor_increment','mass_velocity_correction','velocity_increment','frozen_defect','splitting_identity','predictor_residual','pressure_residual','corrected_free_state'} or any(v>2e-10 for v in step.get('action_reconstruction',{}).values()):reasons.append('native_action:'+dt)
            if step.get('predictor_factor_nonzeros',40000001)>40000000:reasons.append('fill:'+dt)
            expected_terms={region+'/volume/'+name for region in ('inner','outer') for name in
                ('augmented_pressure','advection','molecular_normal','molecular_shear','linear_normal','linear_shear','nonlinear_normal','nonlinear_shear','rough_pressure','rough_inertia')}
            expected_terms|={'prescribed_body'}|{region+'/'+face+'/natural_'+kind for region,face in
                (('inner','outlet'),('outer','outlet'),('outer','top')) for kind in ('traction','advection')}
            if set(step.get('component_images',{}))!=expected_terms:reasons.append('component_coverage:'+dt)
            if not {'total','total-tighter','continuity','continuity-tighter','storage-lifting','storage-lifting-tighter'}<=set(v['label'] for v in step.get('solve_checks',[])):reasons.append('solve_coverage:'+dt)
            for name in ('reconstruction','tighter_repeat','signed_sum_error','telescope','continuity_tighter','storage_lifting_tighter'):
                v=step.get(name,{}).get('maximum_scaled_value')
                if not isinstance(v,(int,float)) or not 0<=v<1e-8:reasons.append(name+':'+dt)
            if any(abs(v)>2e-10 for v in step.get('stage_pressure',{}).values()):reasons.append('stage_pressure:'+dt)
            for solve in step.get('solve_checks',[]):
                if solve['maximum_true_residual']>2e-10*solve['rhs_maximum']+solve['absolute_tolerance']:reasons.append('true_pressure_residual:'+dt)
            for repeat in step.get('component_tighter_repeats',{}).values():
                if repeat['maximum_scaled_value']>=1e-8:reasons.append('component_repeat:'+dt)
    if set(q.get('common',{}))!={'0.001','0.002'}:reasons.append('common_time')
    for t,row in q.get('common',{}).items():
        if row.get('maximum_regional_reconstruction_scaled',1)>=1e-8:reasons.append('regional_pressure:'+t)
        if set(row.get('mechanisms',{}))!={'small_step_spatial','finite_alpha','continuity','storage_lifting'} or row.get('mechanism_signed_sum_error',{}).get('maximum_scaled_value',1)>=1e-8:reasons.append('mechanism_sum:'+t)
    if q.get('bands',{}).get('signed_sum_error',{}).get('maximum_scaled_value',1)>=1e-8:reasons.append('band_sum')
    if require_quadrature:
        for row in q.get('held_quadrature',[]):
            if row.get('status')!='pass':reasons.append('quadrature_capacity')
            elif row['pressure_image']['maximum_scaled_value']>=1e-8 or row['reconstruction']>2e-10:reasons.append('quadrature')
        if len(q.get('held_quadrature',[]))!=2:reasons.append('quadrature_coverage')
    return dict(status='pass' if not reasons else 'fail',reasons=sorted(set(reasons)),scope='read-only reconstruction and signed attribution; no physical uncertainty reclassification')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','ct10','ct11','pressure','axes','profile','self-control','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();pressure=sealed(a.pressure);axes=sealed(a.axes);profile=sealed(a.profile);old=sealed(a.ct11/'target-v1')
    gate=pressure_gate(pressure,False);quadrature_gate=pressure_gate(pressure)
    require(gate['status']=='pass','incomplete pressure reconstruction primitives:'+str(gate))
    require(profile['status']=='pass' and profile['retention']['all_restored_bytes_equal'] and
            profile['rollback']['status']=='pass' and physical_budgets(profile['budgets'])['status']=='pass','profile/retention prerequisites failed')
    require(axes['status']=='pass' and all(r['signed_sum_error']['maximum_scaled_value']<1e-8 for rows in axes['members'].values() for r in rows),'axis reconstruction failed')
    requests={label:dict(path=str(folder),request_sha256=sha(folder/'request-manifest.json'),result_sha256=sha(folder/'results.json'),seal_sha256=sha(folder/'completed-artifacts.json')) for label,folder in [('pressure',a.pressure),('axes',a.axes),('profile',a.profile)]}
    resources={};project=aggregate=0;cpu=wall=0.
    for f in sorted(a.root.glob('resources-*')):
        if f==a.self_control:continue
        require((f/'execution.json').exists(),'unclosed resource family:'+str(f));e=read(f/'execution.json');r=e['resource']
        # A resource failure is a required-scope failure, never a reason to
        # discard otherwise interpretable physical/diagnostic evidence.
        resources[f.name]=dict(execution=e,sha256=sha(f/'execution.json'),highwater={x.name:read(x) if x.suffix=='.json' else [json.loads(line) for line in x.read_text().splitlines()] for x in f.iterdir() if x.name in ('resource-final.json','project-highwater.jsonl','worker-guard-highwater.jsonl','worker-guard-events.jsonl')})
        pause_events=resources[f.name]['highwater'].get('worker-guard-events.jsonl',[])
        active_upper=max(0.,e['wall_seconds']-e['allocation_wait_seconds'])
        resources[f.name]['wall_accounting']=dict(elapsed_seconds=e['wall_seconds'],allocation_wait_seconds=e['allocation_wait_seconds'],
            active_wall_upper_bound_seconds=active_upper,guard_pause_resume_events=pause_events,
            scope='no subtraction of partially resumed child intervals; conservative active-wall upper bound; process CPU is independently measured')
        project=max(project,r['maximum_project']);aggregate=max(aggregate,r['maximum_aggregate']);cpu+=e['process_cpu_seconds'];wall+=active_upper
    # Account for simultaneously active nonnumerical status helpers too. All
    # CT12 controls are identified by their actual root/guard/family PIDs.
    roots={e['execution'][k] for e in resources.values() for k in ('control_parent_pid','guard_pid','owned_family_pid')}
    union_peak=0;union_witness=None
    for e in resources.values():
        for row in e['highwater'].get('worker-guard-highwater.jsonl',[]):
            processes=row['processes'];owned=set(roots)
            while True:
                add={p['pid'] for p in processes if p['parent'] in owned}-owned
                if not add:break
                owned.update(add)
            value=sum(p['threads'] for p in processes if p['pid'] in owned)
            if value>union_peak:union_peak=value;union_witness=dict(time_unix=row['time_unix'],project_threads=value,processes=[p for p in processes if p['pid'] in owned])
    worker_pass=union_peak<=16 and aggregate<=72 and all(e['execution']['resource']['project_status']==e['execution']['resource']['aggregate_status']=='pass' for e in resources.values())
    output_bytes=sum(f.stat().st_size for f in a.root.rglob('*') if f.is_file() and f.relative_to(a.root).parts[0] not in ('source','tmp'))
    package_pass=cpu<12000 and wall<14400 and output_bytes<12*1024**3
    m=RequestManifest(a.output,dict(schema='CT12CompactEvidenceV1',source=sha(Path(__file__)),inputs=requests,resources={k:v['sha256'] for k,v in resources.items()},scope='collection execution,diagnostic,physical and resource decisions distinct'))
    negatives={}
    for name,word,mutate in [('missing_step','missing_time',lambda q:q['members']['0.001']['steps'].pop()),
            ('missing_component','component_coverage',lambda q:q['members']['0.001']['steps'][0]['component_images'].pop('outer/volume/nonlinear_normal')),
            ('bad_row_pass_label','block_reconstruction',lambda q:q['members']['0.001']['steps'][0]['block_reconstruction'].update(G=1e-7)),
            ('bad_pressure_pass_label','reconstruction',lambda q:q['members']['0.001']['steps'][0]['reconstruction'].update(maximum_scaled_value=1e-6)),
            ('wrong_gauge','pressure_level',lambda q:q['members']['0.001']['steps'][0].update(unused_gauge=True)),
            ('missing_common_time','common_time',lambda q:q['common'].pop('0.001'))]:
        v=copy.deepcopy(pressure);mutate(v);v['status']='pass';r=pressure_gate(v,False)
        require(r['status']=='fail' and any(word in x for x in r['reasons']),'nondiscriminating primitive negative:'+name);negatives[name]=r
    initial=read(a.ct10/'target-preflight-v1/target/readiness.json');mom=initial['ct10_momentum'];free=np.asarray(mom['free_rows'],int)
    for name,part,factor in [('omit_nonlinear_normal','nonlinear_normal',1),('reverse_natural_traction','natural_traction',2),('omit_augmented_pressure','augmented_pressure',1)]:
        load=sum((np.asarray(v) for k,v in mom['actual'].items() if part in k),np.zeros(mom['momentum_size']))
        signal=float(np.max(abs(load[free]),initial=0.))*factor;require(signal>2e-10,'nondiscriminating signed-load negative:'+name)
        negatives[name]=dict(status='pass',expected='wrong signed free-row residual',maximum_wrong_row_defect=signal,limit=2e-10)
    categories={};baseline_accepted=0.;replay_accepted=0.;geometry=0.
    for label,row in profile['members'].items():
        baseline_accepted+=row['legacy_exclusive_CPU'];replay_accepted+=row['accepted_exclusive_CPU']
        for key,v in row['categories'].items():
            q=categories.setdefault(key,dict(calls=0,CPU=0.,wall=0.))
            for k in q:q[k]+=v[k]
            if key.startswith('geometry.'):geometry+=v['CPU']
    costs=dict(categories=categories,baseline_uninstrumented_accepted_CPU=baseline_accepted,detailed_accepted_CPU=replay_accepted,
        detailed_minus_baseline_CPU=replay_accepted-baseline_accepted,geometry_exclusive_CPU=geometry,
        instrumentation_caution='detailed clocks can add substantial cost; do not equate instrumented category fractions with achieved optimization savings',
        governing_optimization='none; no unmeasured reuse or ExactStepV1 enabled')
    olddecision=read(a.ct11/'qualified-v1/decision.json');plan=olddecision['closed_campaign_resource_admission']
    native_per=max(v['native']['retained_bytes'] for v in profile['retention']['accepted_only_categories'].values())
    audit_pair=2*max(v['audit']['retained_bytes'] for v in profile['retention']['accepted_only_categories'].values())
    output_future=1.25*(1215*native_per+4*audit_pair)+plan['measured_current_storage_bytes']+output_bytes
    forecast=dict(horizon=.407,additional_pair_steps=1215,segments=39,CPU_limit=24000,active_wall_limit=28800,output_limit_bytes=24*1024**3,
        retained_CT11_closed_CPU_forecast=plan['total_CPU_spent_plus_remaining_with_headroom'],
        conservative_investigation_plus_campaign_CPU=cpu+plan['total_CPU_spent_plus_remaining_with_headroom'],
        retained_CT11_wall_forecast=plan['total_wall_upper_estimate'],conservative_investigation_plus_campaign_wall=wall+plan['total_wall_upper_estimate'],
        measured_new_investigation_CPU=cpu,measured_new_investigation_active_wall_upper=wall,measured_new_output_bytes=output_bytes,
        measured_policy_output_forecast=output_future,CPU_admitted=False,wall_admitted=False,output_admitted=output_future<=24*1024**3,
        scope='spent historical target plus remaining forecast,25percent headroom; actual accepted-only retention excludes rollback dedup; no achieved CPU optimization or future trajectory assumed')
    common=pressure['common']['0.002'];ranked=sorted(common['components'].items(),key=lambda x:abs(x[1]['signed_alignment']),reverse=True)
    summary=dict(schema='CT12PressurePrimitiveSummaryV1',common=pressure['common'],
        top_signed_components=[dict(name=k,**v) for k,v in ranked[:20]],bands=pressure['bands'],
        stage_chain=pressure.get('stage_chain'),row_quality={dt:[{k:r[k] for k in ('time','C_maximum','action_reconstruction','initial_weak_continuity_maximum','reconstruction','tighter_repeat','signed_sum_error','continuity','continuity_tighter','storage_and_lifting','storage_lifting_tighter','finite_alpha_remainder','telescope','stage_pressure','trace','predictor_factor_nonzeros')} for r in member['steps']] for dt,member in pressure['members'].items()},
        capped_diagnostic_solve_flags=[dict(member=dt,time=r['time'],**v) for dt,member in pressure['members'].items() for r in member['steps'] for v in r['solve_checks'] if not v['requested_scipy_stopping_test_met']],
        held_quadrature=pressure['held_quadrature'],stage_loads=pressure['stage_loads'])
    decision=dict(schema='CT12RequiredScopeV1',diagnostic_reconstruction=gate,held_state_quadrature=quadrature_gate,historical_CT11_pressure=old['qualification'],
        target_physical_horizon=.002,new_physical_continuation='not_attempted_or_admitted',
        profile_and_restart=dict(status='pass',accepted_replay_steps=2,deliberate_rejected_attempts=1,governing_change=False),
        retention=dict(status=profile['retention']['status'],scope='lossless compressed/deduplicated bounded replay payloads; no historical pruning or unqualified rolling eviction'),
        campaign_resources=forecast,resources=dict(status='pass' if worker_pass and package_pass else 'fail',worker_status='pass' if worker_pass else 'fail',package_budget_status='pass' if package_pass else 'fail',maximum_project=max(project,union_peak),maximum_aggregate=aggregate,combined_project_witness=union_witness,
            closed_package_CPU=cpu,closed_family_wall_sum=wall,new_output_bytes=output_bytes,package_CPU_limit=12000,package_wall_limit=14400,package_output_limit_bytes=12*1024**3),
        static_contrast='not_yet_selected; target-specific attribution must precede any separate static manifest',
        terrain='unattempted; mapped numbering,joined gradients/cofactors,transmission and curved temporal/spatial gates separate',
        physical_applicability='central physical reference retained; no precursor or hill validation',historical='CT4–CT11 seals unchanged')
    failures={}
    for name in ('pressure-v1','pressure-v2','pressure-v3','pressure-v4','profile-v1','profile-v2'):
        f=a.root/(name+'-family.log');r=a.root/('resources-'+name)/'execution.json'
        if f.exists() and r.exists():failures[name]=dict(log=str(f),log_sha256=sha(f),witness=f.read_text()[-4500:],execution=read(r),execution_sha256=sha(r),
            scope='preserved diagnostic/controller failure, not a new physical trajectory or a relabelled historical result')
    for name,value in [('failed-attempts',failures),('pressure',pressure),('pressure-summary',summary),('axes-and-tractions',axes),('profile',profile),('cost',costs),('forecast',forecast),('negative-controls',negatives),('resources',resources),('decision',decision)]:write_json(a.output/(name+'.json'),value)
    critical=['library/SolverUtils/Turbulence/EquilibriumAWFProfiles.cpp','library/SolverUtils/Turbulence/RANSModels.cpp','library/UnitTests/SolverUtils/RANSPreservationProbe.cpp',
        'solvers/IncNavierStokesSolver/EquationSystems/RANSWallGlobalField.cpp','solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.cpp',
        'solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionOpenBoundary.cpp',
        'solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPimple.cpp','solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionSparse.cpp',
        'solvers/IncNavierStokesSolver/EquationSystems/RANSWallSerialSparse.cpp','solvers/IncNavierStokesSolver/EquationSystems/RANSWallMomentumAttribution.cpp']
    identities={}
    for name in critical:
        original=subprocess.check_output(['git','show','b62e5936:'+name],cwd=a.root/'source')
        digest=hashlib.sha256(original).hexdigest();require(digest==sha(a.root/'source'/name),'governing/diagnostic CPP changed');identities[name]=digest
    write_json(a.output/'provenance.json',dict(source_revision=subprocess.check_output(['git','-c','core.preloadIndex=false','rev-parse','HEAD'],cwd=a.root/'source',text=True).strip(),unchanged_CPP=identities,
        python=sys.version,python_sha256=sha(Path(sys.executable).resolve()),numpy=np.__version__,scipy=scipy.__version__,platform=platform.uname()._asdict(),threads=THREADS,
        governing='pinned CT9 build-v8,unchanged',probe='pinned CT11 build-v2,unchanged',inputs=requests))
    files=[f.name for f in a.output.iterdir() if f.is_file() and f.name not in ('request-manifest.json','completed-artifacts.json')];m.seal(a.output,files)
    print(json.dumps(dict(decision=decision,top_signed=summary['top_signed_components'][:6],cost=costs,retention=profile['retention']),indent=2),flush=True);return 0
if __name__=='__main__':raise SystemExit(main())
