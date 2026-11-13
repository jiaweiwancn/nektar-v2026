#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact CT6 receipts; preserve incomplete/failed scopes and raw identities.

Use the owned Linux launcher after numerical families have closed. Large
matrices, coefficient arrays and packets stay at their original archive paths.
"""
import argparse,hashlib,json,math,os,subprocess
from pathlib import Path
from rans_qualification import THREADS,write_json,exact_vector_hash,audited_extrema
os.environ.update(THREADS)
from rans_request_manifest import RequestManifest,sha,canonical
from ct_accepted_open import require
from ct_accepted_budget_gate import physical_budgets
from ct_accepted6_audit import budget,captured_k_matrix
from ct_accepted6_gate import classify as evolution_gate,PRIMARY
from ct_accepted6_shared_gates import classify as shared_gates
from ct_accepted6_quadrature import extract as quadrature_primitives,classify as quadrature_gate


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();r=a.root;files={};libraries={};missing=[]
    def read(path):return json.loads(Path(path).read_text())
    def compact(x):
        if isinstance(x,list):return [compact(v) for v in x]
        if not isinstance(x,dict):return x
        out={}
        for key,value in x.items():
            if key=='actual_mapped_files':
                h=hashlib.sha256(canonical(value)).hexdigest();libraries[h]=value;out['actual_mapped_files_sha256']=h
            elif key in ('state','initial','accepted_states','rejected_state','history','native_history','reference','fixed_defect','fixed_reference','effective_initial','live_fixed_boundary','matrix','mass','weights','mixed_space'):
                out[key+'_archive_identity']=dict(sha256=hashlib.sha256(canonical(value)).hexdigest(),entries=len(value) if isinstance(value,(list,dict)) else None,location='original sealed Linux receipt/packet')
            else:out[key]=compact(value)
        return out
    def receipt(relative,name=None):
        path=r/relative
        if not path.exists():missing.append(str(path));return None
        payload=read(path);files[str(path)]=sha(path)
        write_json(a.output/(name or relative.replace('/','__')),dict(path=str(path),sha256=sha(path),body=compact(payload)))
        return payload
    g=['git','-c','core.filemode=false','-c','index.threads=1','-c','core.preloadIndex=false']
    revision=subprocess.check_output(g+['rev-parse','HEAD'],cwd=r/'source',text=True).strip()
    manifest=RequestManifest(a.output,dict(schema='CT6CompactRequiredScopeV1',root=str(r),source_revision=revision,collector_sha256=sha(Path(__file__))))
    mathematics=receipt('mathematics-v2/continuous/results.json','continuous-oracle.json')
    kernels=receipt('mathematics-v2/kernels/results.json','curved-kernels.json')
    operator=receipt('native-operator-v1/results.json','native-operator.json')
    frozen=receipt('native-operator-v1/frozen-matrix-classification-v1.json','frozen-matrices.json')
    receipt('primitives-v1/results.json','native-continuous-crosscheck.json')
    initial=receipt('initial-v1/results.json','initial-representation-v1.json')
    static=receipt('static-independent-v2/results.json','static-transmission.json')
    receipt('targeted-v4/results.json','static-native-primitives.json')
    zero=receipt('targeted-v3/results.json','zero-and-ground-controls.json')
    primary_endpoint=receipt('endpoint-classification-v2/results.json','A95-primary-repaired-endpoint.json')
    receipt('first-failure-v1.json','first-flat-guard-failure.json')
    receipt('frozen-boundary-stage-proof-v1.json','wrong-corner-precedence-control.json')
    receipt('frozen-boundary-stage-proof-v2.json','actual-stage-reconstruction-control.json')
    receipt('controls-v1/results.json','first-discriminating-control-attempt.json')
    controls=receipt('controls-v2/results.json','discriminating-controls.json')
    regression=receipt('regression-v1/results.json','regression.json')
    receipt('primary-comparisons-v1/results.json','A95-staged-comparisons.json')
    comparison=receipt('primary-comparisons-v2/results.json','A95-complete-comparisons.json')
    if comparison is None:
        path=r/'primary-comparisons-v1/results.json'
        if path.exists():comparison=read(path)
    builds={}
    for folder in sorted(r.glob('build-v*')):
        request=folder/'request-manifest.json';completed=folder/'results.json';progress=folder/'progress.json'
        if not request.exists():continue
        identity=read(request);files[str(request)]=sha(request)
        source=completed if completed.exists() else progress
        body=read(source) if source.exists() else dict(status='no_completed_build_receipt')
        if source.exists():files[str(source)]=sha(source)
        builds[folder.name]=dict(request=compact(identity),request_sha256=sha(request),
            completed=completed.exists(),body=compact(body),
            source_receipt=str(source),source_receipt_sha256=sha(source) if source.exists() else None)
    write_json(a.output/'build-runtime-identities.json',builds)
    members={};audits={};audit_bodies={};costs={}
    for path in sorted(r.glob('evolution-v*/**/results.json')):
        q=read(path)
        if q.get('schema')!='CT6FreshManufacturedMemberV1':continue
        label=str(path.parent.relative_to(r));folder=path.parent/'trajectory';report_path=folder/'case.evolution.json'
        row=dict(receipt_sha256=sha(path),directory=str(path.parent),model=q['model'],degree=q['degree'],dt=q['dt'],target=q['target'],
                 recorded_status=q['status'],physical_completion=q['physical_completion'],diagnostic_status=q['diagnostic_status'])
        if report_path.exists():
            n=read(report_path);complete=not n['failure'] and math.isfinite(n['time']) and abs(n['time']-q['target'])<1e-13
            corrections=sum(x.get('bound_calls',0)+x.get('additional_governing_equation_solves',0)+x.get('additional_correction_or_restoration_solves',0) for x in n['pimple_attempts'])
            row.update(physical_completion='pass' if complete else 'fail',time=n['time'],failure=n['failure'],steps=n['steps'],stored_slots=len(n['state']),
                actual_report_sha256=sha(report_path),actual_state_sha256=exact_vector_hash(n['state']),actual_initial_sha256=exact_vector_hash(n['initial']),
                actual_history_sha256=sha(Path(n['checkpoint'])/'State.bin') if n['checkpoint'] else None,
                accepted_state_times=n['audit_times'], actual_configuration_sha256=sha(folder/'case.xml'),
                corrections_and_additional_solves=corrections,all_audited_extrema=audited_extrema(n))
            require(corrections==0,'forbidden governing correction or extra solve')
            if 'endpoint_metrics' in q:
                row['independent_endpoint_metrics']=q['endpoint_metrics']
            elif complete and label=='evolution-v2/A95-p4-dt002' and primary_endpoint:
                require(primary_endpoint['initialization_receipt']['raw_initial_sha256']==row['actual_state_sha256'] and
                    abs(primary_endpoint['observation_boundary_time']-n['time'])<1e-13,'foreign repaired endpoint state/time')
                require(sha(Path(primary_endpoint['packet_path']))==primary_endpoint['packet_sha256'],'changed supplemented endpoint packet')
                row['independent_endpoint_metrics']=primary_endpoint['physical_metrics']
                row['diagnostic_supplement']=dict(path=str(r/'endpoint-classification-v2/results.json'),
                    sha256=sha(r/'endpoint-classification-v2/results.json'),original_failed_receipt_preserved=True)
        members[label]=row;files[str(path)]=sha(path)
        for phase in ('zero','trajectory'):
            timing=path.parent/phase/'timing.jsonl'
            if not timing.exists():continue
            exclusive={};roots=[]
            for line in timing.read_text().splitlines():
                t=json.loads(line);name=t['name'];v=exclusive.setdefault(name,dict(process_cpu_seconds=0.,wall_seconds=0.,calls=0))
                v['process_cpu_seconds']+=t['exclusive_process_cpu_seconds'];v['wall_seconds']+=t['exclusive_wall_seconds'];v['calls']+=1
                if t['parent_id']==0:roots.append(t)
            costs[label+'/'+phase]=dict(path=str(timing),sha256=sha(timing),exclusive=exclusive,root_scopes=roots,
                scope='exclusive instrumentation only; nested inclusive totals are not added; child process CPU/wall/RSS remain in execution receipts')
        write_json(a.output/(label.replace('/','__')+'.json'),dict(path=str(path),sha256=sha(path),body=compact(q),actual=row))
    for path in sorted(r.glob('physical-audits-v*/**/results.json')):
        q=read(path);label=str(path.parent.relative_to(r));gates=physical_budgets(q)
        member=next((v for v in members.values() if v.get('actual_report_sha256')==q['budget_binding']['native_report_sha256']),None)
        require(member is not None,'foreign actual-state physical audit')
        integration={}
        for name,row in q['budgets'].items():
            source=Path(row['execution']['command'][-1]);require(sha(source)==row['packet_sha256'],'changed physical budget packet')
            b=read(source);rebuilt=budget(b,q['budget_binding']['model'])
            for key in ('equations','reconstruction','physical_volume'):require(rebuilt[key]==row[key],'forged or stale budget primitives:'+name+':'+key)
            require(exact_vector_hash(b['observation']['state'])==q['budget_binding']['state_sha256'],'foreign consumed budget state')
            integration[name]=quadrature_primitives(b,rebuilt)
            del b
        quad=quadrature_gate(integration,comparison or {},q['budget_binding']['model'])
        audits[label]=dict(directory=str(path.parent),body_sha256=sha(path),recomputed_full_budget_gate=gates,
            recorded_status=q['status'],native_diagnostic=q.get('diagnostic_status'),captured_k_matrix=q.get('captured_k_matrix'),mass=q.get('mass'),quadrature_uncertainty=quad)
        audit_bodies[label]=dict(q,quadrature_uncertainty=quad)
        write_json(a.output/(label.replace('/','__')+'.json'),dict(path=str(path),sha256=sha(path),body=compact(q),recomputed_gate=gates,derived_quadrature_uncertainty=quad))
        files[str(path)]=sha(path)
    resources={}
    for path in sorted(r.glob('launch-*/resource-final.json')):
        q=read(path);row=dict(final=q);files[str(path)]=sha(path)
        launch=path.parent
        execution=launch/'execution.json'
        if execution.exists():
            e=read(execution);files[str(execution)]=sha(execution)
            row['execution']=compact(e)
        for stem in ('worker-guard-highwater','project-highwater','worker-guard-events'):
            source=launch/(stem+'.jsonl')
            if source.exists():
                files[str(source)]=sha(source)
                row[stem]=[json.loads(line) for line in source.read_text().splitlines()]
        # Report per-family-root suspension intervals separately from CPU and
        # elapsed time. Summing every child's stopped duration would count the
        # same wall interval repeatedly.
        stopped={};intervals=[]
        for event in row.get('worker-guard-events',[]):
            for kind,root,pid in event['events']:
                if root!=pid:continue
                if kind=='pause':stopped.setdefault(pid,event['time_unix'])
                elif kind=='resume' and pid in stopped:
                    begin=stopped.pop(pid);intervals.append(dict(pid=pid,begin=begin,end=event['time_unix'],seconds=event['time_unix']-begin))
        row['family_root_suspensions']=dict(closed_intervals=intervals,unclosed_begin_times=stopped,
            total_closed_seconds=math.fsum(v['seconds'] for v in intervals),
            scope='root suspension wall intervals; elapsed and process CPU remain separate; no child-duration double counting')
        resources[str(launch.name)]=row
    shared_input=dict(continuous_oracle=mathematics,curved_kernel=kernels,geometry_operator=operator,
        model_initialization=initial,frozen_pressure_action=frozen,static_transmission=static,
        zero_curvature=zero,negative_controls=controls,native_restart_rollback=regression,resources=resources)
    shared=shared_gates(shared_input)
    selected=[v for v in members.values() if v['model']=='A95' and
        v['physical_completion']=='pass' and (v['degree'],v['dt']) in PRIMARY]
    selected_audits={}
    for label,body in audit_bodies.items():
        report=body['budget_binding']['native_report_sha256']
        if report not in {v.get('actual_report_sha256') for v in selected}:continue
        old=selected_audits.get(report)
        # Prefer a complete four-control receipt. Never combine controls from
        # different requests or reinterpret native-only diagnostic completion.
        if old is None or (set(body['budgets'])=={'native','q48f16','q64f16','q64f32'} and
                           set(old['budgets'])!={'native','q48f16','q64f16','q64f32'}):
            selected_audits[report]=body
    primitive_audits={}
    for report,body in selected_audits.items():
        primitive_audits[report]={k:body[k] for k in ('budget_binding','mass','captured_k_matrix','quadrature_uncertainty')}
        primitive_audits[report]['budgets']={name:{k:row[k] for k in
            ('binding','equations','reconstruction','stage_coefficient_validation','packet_sha256','actual_stage_snapshot_sha256')}
            for name,row in body['budgets'].items()}
    admission_input=dict(schema='CT6EvolutionAdmissionInputV1',members=selected,
        comparison=comparison or {},audits=primitive_audits,
        prerequisites={k:v['status'] for k,v in shared.items()},archive_identities=files)
    admission=evolution_gate(admission_input)
    write_json(a.output/'shared-prerequisites.json',shared)
    write_json(a.output/'evolution-admission-input.json',admission_input)
    write_json(a.output/'evolution-admission.json',admission)
    decision=dict(schema='CT-Accepted-6-DecisionV1',mathematical_continuum=mathematics['status'] if mathematics else 'missing',
        curved_kernel=kernels['status'] if kernels else 'missing',
        conservative_native_spatial_operator='pass' if operator and all(v['status']=='pass' for v in operator.values()) else 'fail_or_missing',
        static_variable_height_transmission=static['status'] if static else 'missing',
        zero_curvature='pass' if zero and all(v['packet']['status']=='pass' for k,v in zero['checks'].items() if '-zero-' in k) else 'fail_or_missing',
        primary_MMS_consistency=comparison['status'] if comparison else 'not_completed',
        discriminating_controls=controls['status'] if controls else 'not_completed',regression=regression['status'] if regression else 'not_completed',
        native_members=members,actual_physical_audits=audits,
        physical_experimental_applicability='not_assessed_or_admitted',curved_production_admission=False,
        pending_or_missing=missing,historical_scope='CT1--CT5 and eight epsilon-order observations unchanged',
        resource_scope='completed launch receipts only; collector closure exported separately after return',
        recomputed_shared_prerequisites=shared,recomputed_evolution_admission=admission,
        required_scopes=dict(curved_Accepted_kernel=shared['curved_kernel'],
            conservative_spatial_operator=shared['geometry_operator'],
            full_traction_and_transported_k_transmission=shared['static_transmission'],
            constant_height_evolution=dict(A95=dict(
                bounded_completion=admission['primary_coverage']['status'],
                evolved_consistency=admission['primary_comparison']['status'],
                full_physical_budgets=admission['complete_bound_audit_coverage']),
                A93='not_attempted_before_complete_A95_admission',
                E94E95='not_attempted_before_complete_A95_admission'),
            static_variable_height=shared['static_transmission'],zero_curvature=shared['zero_curvature'],
            resources=shared['resources'],physical_experimental_applicability='not_assessed'),
        status='not_fully_qualified' if admission['status']!='pass' else 'conditional_model_transfers_required')
    write_json(a.output/'decision.json',decision);write_json(a.output/'completed-resources.json',resources)
    write_json(a.output/'native-costs.json',costs)
    write_json(a.output/'actual-loaded-libraries.json',libraries);write_json(a.output/'archive-identities.json',files)
    write_json(a.output/'source-scope.json',dict(source_revision=revision,starting_revision='7d32b43dc9d176628cfef1dbae2aee02a86ce330',
        source_status=subprocess.check_output(g+['status','--porcelain'],cwd=r/'source',text=True),
        changes=subprocess.check_output(g+['diff','--stat','7d32b43dc9d176628cfef1dbae2aee02a86ce330','HEAD'],cwd=r/'source',text=True),
        guards='named original-Natural Accepted serial2D constant-height Gaussian manufactured evolution only; static variable-height transmission separate; no production/EP/MPI/3D/action support'))
    manifest.seal(a.output,[str(p.relative_to(a.output)) for p in a.output.glob('*.json') if p.name not in ('request-manifest.json','completed-artifacts.json')])
    print('CT6 compact evidence exported; primitive admission:',admission['status'],flush=True)
    return 3  # A completed export is not a completed required-scope qualification.


if __name__=='__main__':raise SystemExit(main())
