#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact CT4 decisions recomputed from complete bound primitives."""
import argparse,hashlib,json,os,subprocess,platform,sys
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
import scipy
from ct_accepted4_reference import write_json
from ct_accepted4_run import qualify,expected,traction
from ct_accepted4_jobs import CASES,label
from ct_accepted_budget_gate import physical_budgets
from ct_accepted_pilot import artifact_hash,budget_classification
from ct_accepted_open import require
from rans_request_manifest import RequestManifest,sha,runtime_identity,canonical
from ct_accepted2_audit import values
from ct_accepted3_costs import scopes


def read(p):return json.loads(p.read_text())


def projected_gradients(packet,state,c,n,weights):
    """Physical derivative errors; never replace a gradient in evolution."""
    truth,_=expected(packet,c);o=values(packet,state,n)['observation']
    a=np.array(o['grid_gradients']);b=np.array(truth['observation']['grid_gradients'][:n]);w=np.array(weights)
    result={}
    for j,name in enumerate(('u','w','k','outer_epsilon','Pi')):
        mask=np.array(packet['observation']['grid_regions'][:n],bool) if name=='outer_epsilon' else np.ones(n,bool)
        scale=1. if j in (0,1,4) else c['K0']/c['H'] if j==2 else c['reference_column']['counterfactual']['epsilon_amplitude']/c['H']**2
        d=a[mask,j,:]-b[mask,j,:]
        result[name]=dict(weighted_RMS=float(np.sqrt(w[mask]@np.sum(d*d,axis=1)/w[mask].sum())),scale=scale)
        result[name]['scaled_weighted_RMS']=result[name]['weighted_RMS']/scale
    return dict(fields=result,scope='actual represented gradients minus continuous physical reference; separate from continuous-provider agreement and field-value errors')


def stage_summary(rows):
    """Complete scalar receipts remain intact; bulky stage inventories are hashed."""
    winners={};result=[]
    for row in rows:
        stages={}
        for stage,lists in row['stages'].items():
            stages[stage]={}
            for inventory,q in lists.items():
                small={k:q[k] for k in ('count','domain_failure_count','divergence_unweighted_RMS','divergence_maximum')}
                stages[stage][inventory]=small
                for name,v in q.get('minima',{}).items():
                    key=inventory+'/'+name
                    if key not in winners or v['value']<winners[key]['value']:winners[key]=dict(v,time=row['time'],stage=stage)
                v=q.get('minimum_radicand')
                if v is not None:
                    key=inventory+'/A95_radicand'
                    if key not in winners or v['radicand']<winners[key]['radicand']:winners[key]=dict(v,time=row['time'],stage=stage)
                v=q.get('minimum_inner_k')
                if v is not None:
                    key=inventory+'/inner_k'
                    if key not in winners or v['value']<winners[key]['value']:winners[key]=dict(v,time=row['time'],stage=stage)
        result.append(dict(time=row['time'],stages=stages,**{k:v for k,v in row.items() if k not in ('time','stages')}))
    return dict(all_stage_scalar_metrics=result,all_stage_winning_witnesses=winners,external_full_stage_array_count=len(rows),
        external_full_stage_array_sha256=hashlib.sha256(canonical(rows)).hexdigest())


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','runtime','initial-root','evolution-root','analysis','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();source=a.root/'source';git=['git','-c','index.threads=1','-c','core.preloadIndex=false','-c','pack.threads=1']
    os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
    inputs=[*a.initial_root.glob('*/results.json'),*a.evolution_root.glob('*/*/results.json'),
        a.root/'controls-v1/results.json',a.root/'legacy-controls-v1/results.json',a.root/'reference-v1/results.json',a.analysis/'results.json']
    m=RequestManifest(a.output,dict(schema='CT4RequiredScopeV1',source_sha256=sha(Path(__file__)),runtime=runtime_identity(a.runtime/'IncNavierStokesSolver'),
        complete_input_receipts={str(f):sha(f) for f in inputs}))
    exports={}
    def compact(value):
        if isinstance(value,list):return [compact(v) for v in value]
        if not isinstance(value,dict):return value
        out={}
        for k,v in value.items():
            if k=='stage_diagnostics':out[k]=stage_summary(v);continue
            if k in ('effective_initial','state','accepted_states','rejected_state','live_fixed_boundary') and isinstance(v,list):
                out[k+'_external_array']=dict(count=len(v),canonical_sha256=hashlib.sha256(canonical(v)).hexdigest(),scope='full vector/map retained in hashed Linux source artifact; native raw byte identities retained separately')
            else:out[k]=compact(v)
        return out
    def save(name,value):write_json(a.output/name,compact(value))
    def export(name,f):save(name,read(f));exports[name]=dict(source=str(f),sha256=sha(f))
    reference=read(a.root/'reference-v1/results.json');analysis=read(a.analysis/'results.json');controls=read(a.root/'controls-v1/results.json');legacy=read(a.root/'legacy-controls-v1/results.json')
    decision=dict(schema='CT-Accepted-4-DecisionV1',diagnostic_execution='pass',mean_fit='computed_fixed_membership_not_experimental_validation',
        mathematical_rough='pass' if len(reference['cases'])==9 and all(q['classification']=='admissible_reference' for q in reference['cases'].values()) else 'fail',
        provider_implementation='pass' if controls['implementation']=='pass' and legacy['status']=='pass' and analysis['analytical_weak_status']=='pass' else 'fail',
        spatial_sensitivity=analysis['spatial_sensitivity'],cases={},precursor_admission=False,benchmark_admission=False,
        future_terrain=dict(admission=False,thin_0p1='failed for all5.5/6/6.5mm; historical criterion retained',thin_0p25='separate geometric flag only',physical_roughness_sublayer='unestablished',smooth_preparation='unresolved',native_curved='not_attempted'))
    failures=[];initial_representations={};representations={};costs=[]
    for model,h,degree in CASES:
        key=label(model,h,degree);ini_file=a.initial_root/key/'results.json';ini=read(ini_file)
        require(ini['case']['model']==model and ini['case']['H_mm']==h and ini['case']['degree']==degree,'foreign_initial_case')
        export(key+'-initial.json',ini_file);row=dict(initialization=ini['status'],slots={dt:v.get('slots') for dt,v in ini['members'].items()},milestones={},required_times=[.008,.02] if h==6. and degree==3 else [.008]);decision['cases'][key]=row
        representations[key]={}
        for dt,v in ini['members'].items():
            folder=a.initial_root/key/dt;packet=read(folder/'initial.json')
            representations[key][dt]=dict(initial=projected_gradients(packet,packet['observation']['state'],ini['case'],ini['global_count'],ini['weights']),milestones={})
        for target in row['required_times']:
            file=a.evolution_root/key/('T'+str(target))/'results.json'
            if not file.exists():
                row['milestones'][str(target)]=dict(status='not_admitted_or_incomplete',reason='initial or previous paired gate failed; no substitution from another model')
                if ini['status']=='pass' and (target==.008 or row['milestones'].get('0.008',{}).get('status')=='pass'):decision['diagnostic_execution']='incomplete_required_numerical_receipt'
                continue
            q=read(file);require(q['case']==ini['case'] and q['target']==target,'foreign_case_or_horizon_receipt')
            recomputed=qualify(q);require(recomputed==q['qualification'],'primitive_qualification_differs_from_label')
            for dt,v in q['members'].items():
                folder=Path(v['directory']);native=read(folder/'case.evolution.json')
                require(sha(folder/'case.evolution.json')==v['native_report_sha256'] and exact_vector_hash(native['state'])==v['state_sha256'],'changed_native_state')
                require(artifact_hash(native['checkpoint'])==v['history_sha256'],'changed_native_history')
                require(sha(folder/'case.xml')==v['identity']['configuration_sha256'] and v['identity']['runtime']==runtime_identity(a.runtime/'IncNavierStokesSolver'),'changed_configuration_or_runtime')
                for control,b in v['budgets'].items():
                    f=Path(b['execution']['command'][-1]);require(sha(f)==b.get('packet_sha256',b.get('raw_packet_sha256')),'changed_budget_packet')
                    if b.get('status')=='failed_probe':continue
                    packet=read(f);require(packet['observation']['state']==native['state'],'budget_actual_state_mismatch')
                    derived=budget_classification(packet,ini['case'])
                    require(all(derived[k]==b[k] for k in ('equations','reconstruction','physical_volume')),'forged_budget_primitives')
                require(physical_budgets(v)==v['budget_gate'],'budget_gate_primitives_changed')
                packet=read(a.initial_root/key/dt/'initial.json')
                representations[key][dt]['milestones'][str(target)]=dict(final_gradients=projected_gradients(packet,native['state'],ini['case'],ini['global_count'],ini['weights']),
                    accepted_state_traction=[dict(time=t,**{k:v for k,v in traction(packet,s,ini['case']).items() if k!='interface_rows'})
                        for t,s in zip(native['audit_times'],native['accepted_states'],strict=True)],
                    scope='read-only reconstruction; first accuracy crossing is not a production constitutive rejection; continuation stays stopped at the declared paired milestone')
                timing=folder/'timing.jsonl'
                costs.append(dict(case=key,dt=dt,start=q['start'],target=target,actual_time=native['time'],slots=len(native['state']),
                    execution=v['execution'],timing=scopes(timing,q['start'],native['time']) if timing.exists() else None,
                    timing_sha256=sha(timing) if timing.exists() else None,native_report_sha256=sha(folder/'case.evolution.json')))
                if v['failure']:
                    events=native.get('physical_range_events',[]);first=events[0] if events else None
                    failures.append(dict(case=key,dt=dt,last_accepted_time=native['time'],failure=v['failure'],native_report_sha256=sha(folder/'case.evolution.json'),
                        rejected_state_sha256=exact_vector_hash(native['rejected_state']) if native.get('rejected_state') else None,
                        first_failure_event=first,history_sha256=v['history_sha256'],last_attempt=native['pimple_attempts'][-1] if native.get('pimple_attempts') else None))
            name=key+'-T'+str(target)+'.json';export(name,file)
            row['milestones'][str(target)]=dict(recomputed,complete_body=name,times={dt:v['time'] for dt,v in q['members'].items()})
            if recomputed['status']!='pass':
                failures.append(dict(case=key,target=target,kind='paired_qualification_gate',classification=recomputed,
                    receipt_sha256=sha(file),actual_rejection=any(bool(v['failure']) for v in q['members'].values()),
                    member_traction={dt:v['traction'] for dt,v in q['members'].items()},
                    scope='physical completion and failed numerical accuracy remain distinct; no T=.02 continuation from this receipt'))
        row['bounded_numerical_admission']='pass' if ini['status']=='pass' and all(q['status']=='pass' for q in row['milestones'].values()) else 'failed_or_incomplete'
        for version in ('initial-v1','initial-v2'):
            f=a.root/version/key/'results.json'
            if f.exists():initial_representations[version+'/'+key]=dict(source=str(f),sha256=sha(f),body=read(f))
            else:
                f=a.root/version/key/'0.002/case.evolution.json'
                if f.exists():initial_representations[version+'/'+key]=dict(source=str(f),sha256=sha(f),native_preparation_failure=read(f))
    export('mathematical-fit-uncertainty.json',a.root/'reference-v1/results.json');export('observation-audit.json',a.root/'observation-checks-v1/results.json')
    export('provider-history-controls.json',a.root/'controls-v1/results.json');export('legacy-compatibility.json',a.root/'legacy-controls-v1/results.json')
    export('native-weak-and-degree-analysis.json',a.analysis/'results.json');export('pre-evolution-weak.json',a.root/'pre-evolution-weak-v1/results.json')
    for v in ('build-v1','build-v2','build-v3'):export(v+'.json',a.root/v/'results.json')
    save('prior-preparation-representations.json',initial_representations);save('native-failure-witnesses.json',failures)
    save('represented-gradients-and-traction.json',representations);save('costs.json',costs)
    decision['native_numerical']='pass' if all(q['bounded_numerical_admission']=='pass' for q in decision['cases'].values()) else 'failed_accuracy_or_incomplete'
    decision['scope']='implementation/mathematical passes do not admit failed native accuracy, developed precursor or hill validation; prior CT1-CT3 failures and classifications unchanged'
    resources={};peaks=[];events=[];pending=[]
    for folder in sorted(a.root.glob('launch-*')):
        final=folder/'resource-final.json'
        if final.exists():resources[folder.name]=dict(read(final),execution=read(folder/'execution.json'))
        else:pending.append(str(folder))
        for name,target in [('worker-guard-highwater.jsonl',peaks),('project-highwater.jsonl',peaks),('worker-guard-events.jsonl',events)]:
            if (folder/name).exists():target.extend(dict(control=folder.name,kind=name,**json.loads(line)) for line in (folder/name).read_text().splitlines())
    decision['resource_compliance']=dict(status='pass' if resources and all(q['aggregate_status']==q['project_status']=='pass' for q in resources.values()) else 'fail',
        aggregate_peak=max(q['maximum_aggregate'] for q in resources.values()),project_peak=max(q['maximum_project'] for q in resources.values()),
        aggregate_limit=72,project_limit=16,pending_own_finalizer=pending,historical_incidents='CT1-CT3 unchanged')
    save('resources.json',dict(completed=resources,peaks=peaks,events=events));save('decision.json',decision)
    requests={str(f.relative_to(a.root)):dict(sha256=sha(f),body=read(f)) for f in a.root.glob('**/request-manifest.json') if a.output not in f.parents and source not in f.parents}
    save('requests.json',requests);save('export-provenance.json',exports)
    previous=read(source/'docs/plans/cao-tamura-accepted-2-evidence/final-provenance.json')
    specs={f:sha(source/f) for f in previous['specifications']};dataset={f:sha(source/'docs/plans/cao-tamura-validation-data'/f) for f in previous['frozen_dataset']['hashes']}
    kernels={f:sha(source/f) for f in previous['unchanged_reused_kernel_sources']}
    require(specs==previous['specifications'] and dataset==previous['frozen_dataset']['hashes'] and kernels==previous['unchanged_reused_kernel_sources'],'changed_specifications_dataset_or_original_kernels')
    cache=a.root/'build-v1/build/CMakeFiles/IncNavierStokesSolver.dir/flags.make'
    save('source-scope.json',dict(starting_revision='50bd04c315283e996b53c360d8fa7a67877f02f3',current_source_revision=subprocess.check_output(git+['rev-parse','HEAD'],cwd=source,text=True).strip(),
        source_status=subprocess.check_output(git+['status','--porcelain'],cwd=source,text=True),specifications=specs,dataset=dataset,original_kernel_sources=kernels,
        runtime={n:runtime_identity(a.runtime/n) for n in ('IncNavierStokesSolver','RANSPreservationProbe','RANSNewModelsProbe')},
        actual_build_flags=cache.read_text(),actual_build_flags_sha256=sha(cache),python=sys.version,numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform(),threads={k:os.environ.get(k) for k in THREADS},
        mount=json.loads(subprocess.check_output(['findmnt','-T',str(a.root),'-J'],text=True)),
        scope='original governing equations unchanged; new opt-in continuous physical input and constrained native initial/live boundary projection; no interior restoration after initialization'))
    m.seal(a.output,[f.name for f in a.output.iterdir() if f.is_file() and f.name not in ('request-manifest.json','completed-artifacts.json')])
    print(json.dumps(dict(implementation=decision['provider_implementation'],mathematical=decision['mathematical_rough'],cases={k:v['bounded_numerical_admission'] for k,v in decision['cases'].items()},resources=decision['resource_compliance'])));return 0
if __name__=='__main__':raise SystemExit(main())
