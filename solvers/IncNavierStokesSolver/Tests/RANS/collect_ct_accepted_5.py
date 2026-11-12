#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact CT5 decisions, recomputed primitives and sealed source bindings."""
import argparse,copy,hashlib,json,os,subprocess
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from ct_accepted5_run import basic_gate,traction_gate,physical
from ct_accepted5_traction import profiles
from ct_accepted4_run import traction
from ct_accepted4_reference import write_json
from ct_accepted_pilot import artifact_hash,budget_classification
from ct_accepted_budget_gate import physical_budgets
from ct_accepted_open import require
from ct_accepted3_costs import scopes
from rans_request_manifest import RequestManifest,sha,runtime_identity,canonical

def read(p):return json.loads(Path(p).read_text())

def validate_profile(p,c):
    g=p['ground'];require(g and all(q['weight']>0 and np.isfinite(q['weight']) and np.isfinite(q['actual_traction']) for q in g),'invalid_traction_points')
    require(all(abs(q['D']*q['S']-q['actual_traction'])<2e-10*c['ustar']**2 for q in g),'inconsistent_point_traction')
    mean=sum(q['weight']*q['actual_traction'] for q in g)/sum(q['weight'] for q in g);err=abs(mean-c['ustar']**2)/c['ustar']**2
    require(abs(mean-p['mean_wall_shear'])<2e-10*c['ustar']**2 and abs(err-p['normalized_wall_traction_error'])<2e-10,'inconsistent_traction_summary')
    pairs={}
    for q in p['interface']:pairs.setdefault(round(q['coordinate'][0],14),[]).append(q)
    num=den=0.
    for qs in pairs.values():
        require(len(qs)==2 and {q['side'] for q in qs}=={-1,1},'missing_interface_traction_side')
        if any(q['status']!='defined' for q in qs):continue
        d=np.array(qs[0]['total_traction'])-qs[1]['total_traction'];num+=qs[0]['weight']*float(d@d);den+=qs[0]['weight']
    if p['normalized_interface_traction_RMS'] is not None:require(den>0 and abs(np.sqrt(num/den)/c['ustar']**2-p['normalized_interface_traction_RMS'])<2e-10,'inconsistent_interface_traction_summary')

def summarize_profile(p):
    names=('state_sha256','frozen_endpoint_state_sha256','mean_wall_shear','normalized_wall_traction_error','normalized_interface_traction_RMS','wall_spatial_standard_deviation','minimum_wall','maximum_wall','undefined_interface')
    return dict(**{k:p[k] for k in names},ground_array_sha256=hashlib.sha256(canonical(p['ground'])).hexdigest(),ground_count=len(p['ground']),
        interface_profiles=[{k:q[k] for k in ('coordinate','weight','geometry','region','side','endpoint_K','Pi','total_traction','status') if k in q} for q in p['interface']],interface_jump=p['interface_jump'])

def negative_controls(q,initial_hash):
    out={}
    def reject(name,mutate,reason):
        x=copy.deepcopy(q);mutate(x);b=basic_gate(x['members'],x['temporal'],x['target'],x['dts'],initial_hash)
        require(b['status']=='fail' and any(reason in r for r in b['reasons']),'negative_not_discriminating:'+name);out[name]=dict(status='pass',expected_reason=reason,observed=b['reasons'])
    dt=str(q['dts'][0]);reject('missing_outer_epsilon',lambda x:x['members'][dt]['budgets']['native']['equations'].pop('outer_epsilon'),'physical_budget_primitives')
    reject('missing_member',lambda x:x['members'].pop(dt),'missing_timestep_member')
    reject('unequal_endpoint',lambda x:x['members'][dt].update(time=.007),'incomplete_or_rejected_endpoint')
    p=copy.deepcopy(q['members'][dt]['traction_history'][-1]['profile']);p['mean_wall_shear']=0.;p['normalized_wall_traction_error']=0.
    try:validate_profile(p,q['case'])
    except (ValueError,RuntimeError) as e:require('inconsistent_traction_summary' in str(e),'wrong_forged_traction_reason');out['pass_label_bad_traction_summary']=dict(status='pass',reason=str(e))
    else:raise RuntimeError('forged_traction_summary_accepted')
    return out

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','ct4','runtime','campaign','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');run=read(a.campaign/'results.json');diag=read(a.root/'diagnosis-v2/results.json');obs=read(a.root/'observations-v1/results.json')
    m=RequestManifest(a.output,dict(schema='CT5RequiredScopeV1',campaign_sha256=sha(a.campaign/'results.json'),diagnosis_sha256=sha(a.root/'diagnosis-v2/results.json'),observation_sha256=sha(a.root/'observations-v1/results.json'),
        runtime=runtime_identity(a.runtime/'IncNavierStokesSolver'),sources={f:sha(Path(__file__).with_name(f)) for f in ('collect_ct_accepted_5.py','ct_accepted5_run.py','ct_accepted5_traction.py')}))
    files={};libraries={};cost=[];failures=[]
    def compact(v):
        if isinstance(v,list):return [compact(x) for x in v]
        if not isinstance(v,dict):return v
        out={}
        for k,x in v.items():
            if k=='profile' and isinstance(x,dict) and 'ground' in x:out[k]=summarize_profile(x)
            elif k=='actual_mapped_files':
                h=hashlib.sha256(canonical(x)).hexdigest();libraries[h]=x;out['actual_mapped_files_sha256']=h
            elif k=='stage_diagnostics':out[k]=dict(count=len(x),full_array_sha256=hashlib.sha256(canonical(x)).hexdigest(),source='full native member receipt retained on Linux')
            elif k=='weights':out[k]=dict(count=len(x),sha256=hashlib.sha256(canonical(x)).hexdigest(),source='actual physical quadrature weights retained in initial receipt')
            else:out[k]=compact(x)
        return out
    def save(name,value):write_json(a.output/name,compact(value))
    decision=dict(schema='CT-Accepted-5-DecisionV1',observation_audit=obs['status'],load_localization=diag['status'],models={},mathematical='CT4 central rough reference evidence reused at unchanged continuous inputs',
        future_terrain_admission=False,historical_pressure_interpretation='CT4 unresolved flag retained',scope='one selected split and6mm models only; no9case repeat/.2/curved production/hill prediction')
    member_seen=set();negative=None
    for model,row in run['models'].items():
        if 'initial' not in row:decision['models'][model]=row;continue
        ini=row['initial'];prep=ini['initialization'];case=prep['case'];initial_hash=next(iter(prep['members'].values()))['state_sha256']
        require(ini['physical_input_before']==ini['physical_input_after']==physical(case),'changed_continuous_input')
        decision['models'][model]=dict(initialization=prep['status'],native_counts={dt:v['slots'] for dt,v in prep['members'].items()},milestones={},status=row['status'])
        for label,q in row['milestones'].items():
            require(q['case']==case and q['target'] in (.008,.02),'foreign_case_or_horizon')
            for dt,member in q['members'].items():
                f=Path(member['directory']);native=read(f/'case.evolution.json');require(sha(f/'case.evolution.json')==member['native_report_sha256'] and artifact_hash(native['checkpoint'])==member['history_sha256'],'changed_native_state_or_history')
                require(exact_vector_hash(native['state'])==member['state_sha256'] and member['identity']['model']==model and member['identity']['runtime']==runtime_identity(a.runtime/'IncNavierStokesSolver'),'foreign_model_state_or_runtime')
                root_initial=a.root/('initial-quarter' if dt=='0.0005' else 'initial')/model/dt;packet=read(root_initial/'initial.json')
                require(exact_vector_hash(packet['observation']['state'])==initial_hash,'changed_actual_initial')
                for record in member['traction_history']:validate_profile(record['profile'],case)
                rebuilt=profiles(packet,native['state'],case);require(rebuilt==member['traction_history'][-1]['profile'],'traction_not_from_actual_native_coefficients')
                for control,b in member['budgets'].items():
                    file=Path(b['execution']['command'][-1]);require(sha(file)==b.get('packet_sha256',b.get('raw_packet_sha256')),'changed_budget_packet')
                    if b.get('status')=='failed_probe':continue
                    raw=read(file);expected=budget_classification(raw,case)
                    require(raw['observation']['state']==native['state'] and all(expected[k]==b[k] for k in ('equations','reconstruction','physical_volume')),'foreign_or_forged_budget_primitives')
                require(not any(x.get('bound_calls',0) or x.get('additional_governing_equation_solves',0) or x.get('additional_correction_or_restoration_solves',0) for x in native['pimple_attempts']),'forbidden_runtime_correction')
                if str(f) not in member_seen:
                    member_seen.add(str(f));timing=f/'timing.jsonl';cost.append(dict(model=model,dt=dt,time=native['time'],execution=member['execution'],timing=scopes(timing,native['audit_times'][0],native['time']),native_sha256=sha(f/'case.evolution.json')))
                    if native['failure']:
                        failures.append(dict(model=model,dt=dt,last_accepted_time=native['time'],failure=native['failure'],first_event=native.get('physical_range_events',[None])[0] if native.get('physical_range_events') else None,
                            rejected_state_sha256=exact_vector_hash(native['rejected_state']) if native.get('rejected_state') else None,last_attempt=native['pimple_attempts'][-1] if native['pimple_attempts'] else None,
                            native_report_sha256=sha(f/'case.evolution.json'),directory=str(f)))
            b=basic_gate(q['members'],q['temporal'],q['target'],q['dts'],initial_hash);t=traction_gate(q['members'],case,q['target'])
            require(b==q['basic'] and t==q['traction_uncertainty'],'primitive_milestone_recomputation_differs')
            decision['models'][model]['milestones'][label]=dict(target=q['target'],dts=q['dts'],basic=b,traction=t,qualification=q['qualification'])
            if negative is None and b['status']=='pass':negative=negative_controls(q,initial_hash)
        save(model+'-complete-comparisons.json',row)
    save('native-costs.json',cost);save('native-failure-witnesses.json',failures);save('negative-controls.json',negative or dict(status='not_attempted_no_complete_primitive_fixture'))
    for row in diag['accepted'].values():row['traction']=summarize_profile(row['traction'])
    for row in diag['stages']:
        row['traction']={k:summarize_profile(v) for k,v in row['traction'].items()}
    save('load-and-stage-diagnosis.json',diag)
    obs['height_audit_inventory']=dict(count=len(obs['height_audit']),sha256=hashlib.sha256(canonical(obs.pop('height_audit'))).hexdigest(),source='full versioned observation results retained on Linux')
    save('observation-correction.json',obs)
    resources={}
    for folder in sorted(a.root.glob('launch-*')):
        if (folder/'resource-final.json').exists():resources[folder.name]=dict(read(folder/'resource-final.json'),execution=read(folder/'execution.json'))
    decision['resource_compliance']=dict(status='pass' if resources and all(v['aggregate_status']==v['project_status']=='pass' for v in resources.values()) else 'fail',
        aggregate_peak=max(v['maximum_aggregate'] for v in resources.values()),project_peak=max(v['maximum_project'] for v in resources.values()),pending_own_finalizer=True)
    decision['bounded_flat_numerical']='pass' if len(decision['models'])==3 and all(q['status']=='pass' for q in decision['models'].values()) else 'failed_or_unresolved'
    decision['next_capability']='conservative curved Accepted manufactured patch, retaining physical roughness-sublayer and thin-layer qualifications' if decision['bounded_flat_numerical']=='pass' else 'one mechanism-supported bounded numerical proposal; no model/mesh search'
    save('decision.json',decision);save('resources.json',resources);save('actual-loaded-libraries.json',libraries)
    source=a.root/'source';g=['git','-c','index.threads=1','-c','core.preloadIndex=false'];old=read(source/'docs/plans/cao-tamura-accepted-4-evidence/compact-v1/source-scope.json')
    require(all(sha(source/f)==v for f,v in old['original_kernel_sources'].items()),'changed_governing_kernel')
    require(all(sha(source/f)==v for f,v in old['specifications'].items()),'changed_specifications')
    require(all(sha(source/'docs/plans/cao-tamura-validation-data'/f)==v for f,v in old['dataset'].items()),'changed_frozen_observations')
    changed=subprocess.check_output(g+['diff','--name-only','b9882937..HEAD','--','library','solvers/IncNavierStokesSolver/EquationSystems'],cwd=source,text=True);require(not changed,'changed_native_code_requires_new_runtime')
    save('provenance.json',dict(starting_revision='b9882937a5096aa9c73219dd3b6288cd8a346983',source_revision=subprocess.check_output(g+['rev-parse','HEAD'],cwd=source,text=True).strip(),
        source_status=subprocess.check_output(g+['status','--porcelain'],cwd=source,text=True),runtime=runtime_identity(a.runtime/'IncNavierStokesSolver'),runtime_probe=runtime_identity(a.runtime/'RANSPreservationProbe'),
        governing_source='123f31c974e0b75e551e43a5203749a89378fade',CT4_specifications=old['specifications'],CT4_dataset=old['dataset'],CT4_unchanged_kernels=old['original_kernel_sources'],
        request_manifests={str(f.relative_to(a.root)):dict(sha256=sha(f),body=read(f)) for f in a.root.glob('**/request-manifest.json') if source not in f.parents and a.output not in f.parents},
        full_campaign=dict(path=str(a.campaign/'results.json'),sha256=sha(a.campaign/'results.json')),full_diagnosis=dict(path=str(a.root/'diagnosis-v2/results.json'),sha256=sha(a.root/'diagnosis-v2/results.json')),full_observation=dict(path=str(a.root/'observations-v1/results.json'),sha256=sha(a.root/'observations-v1/results.json')),
        threads=THREADS,historical_scope='CT1-CT4 sealed records and decisions unchanged; CT4 observation claims corrected separately, not silently rewritten'))
    m.seal(a.output,[f.name for f in a.output.iterdir() if f.is_file() and f.name not in ('request-manifest.json','completed-artifacts.json')]);print('CT5 required scope',decision['bounded_flat_numerical']);return 0
if __name__=='__main__':raise SystemExit(main())
