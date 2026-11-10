#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Complete configuration/time-aware export; required failed gates exit nonzero."""
import argparse,gzip,json,os,platform,shutil,subprocess,sys
from pathlib import Path
from rans_qualification import THREADS,write_json
os.environ.update(THREADS)
from qualify_rans_preservation import load
from rans_request_manifest import RequestManifest,sha
from rans_duration_contract import require,required_times,prefix,MEMBERS,aggregate,identity
from collect_rans_balanced_spatial import compact_analysis,audit_summaries
from rans_nonlinear_evidence_contract import classify


def receipt(path):
    if not path.exists():return dict(status='not_completed',source=str(path))
    return dict(load(path),source=str(path),source_sha256=sha(path))


def collect():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--admissions',type=Path,nargs='+',required=True);p.add_argument('--regression-receipt',type=Path);p.add_argument('--forward-interface-receipt',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=a.root.resolve();out=a.output;out.mkdir(parents=True,exist_ok=False);result=dict(schema='NonlinearDurationResponseEvidenceV1',starting_revision='dfd94161fd9d97f7f1a8ce2966e9cc9461b6d686',source_report=subprocess.check_output(['git','-C',str(r/'source'),'rev-parse','HEAD'],text=True).strip(),environment=dict(platform=platform.platform(),python=sys.version,threads=THREADS,ranks=1,worker_limit=72),admissions={},cohorts={},cases=[],responses={},failures=[],verification={})
    for name in ('workspace.json','historical-receipt-byte-verification.json','actual-initial-regression.json','resource-census-exception.json','worker-guard-status.json'):
        result['verification'][name]=receipt(r/name)
    result['verification']['aggregate-negative-controls']=receipt(r/'aggregate-negative-controls/results.json')
    result['verification']['complete-export-negative-controls']=receipt(r/'complete-export-negative-controls/results.json')
    for name in ('ctest-current.json','partial-export-negative-control.json','neutral-as-perturbed-negative-control.json','physical-budget-regression.json'):
        result['verification'][name]=receipt(r/name)
    for name in ('ctest-final-v4.json','ctest-final-v5.json','build-ctest-descriptor.json','build-probe-exact-scope-control.json'):
        result['verification'][name]=receipt(r/name)
    for name in ('forward-interface-tolerance-status.json','physical-budget-face-regression.json','build-probe-side-control.json'):
        result['verification'][name]=receipt(r/name)
    result['execution_schedules']=[]
    for folder in sorted((r/'duration-schedules').glob('*')):
        if not folder.is_dir():continue
        row=dict(source=str(folder),result=receipt(folder/'results.json'),request=receipt(folder/'request-manifest.json'),logs=[])
        for path in sorted(folder.glob('*.log')):
            text=path.read_text(errors='replace')
            row['logs'].append(dict(source=str(path),sha256=sha(path),tail=text[-5000:]))
        result['execution_schedules'].append(row)
    result['common_scalar_partitions']={p.parent.name:receipt(p) for p in (r/'common-scalar-partitions').glob('*/results.json')}
    result['common_inner_qualification']={p.parent.name:receipt(p) for p in (r/'common-inner-qualification').glob('*/results.json')}
    result['residual_failure_diagnoses']={p.parent.name:receipt(p) for p in (r/'residual-diagnosis').glob('*/results.json')}
    result['volume_ownership_diagnoses']={p.parent.name:receipt(p) for p in (r/'volume-ownership').glob('*/results.json')}
    result['failed_response_spatial_audits']={}
    for path in (r/'evolved-audits').glob('response-*/results.json'):
        value=receipt(path)
        if value.get('status')=='fail':result['failed_response_spatial_audits'][path.parent.name]=dict(value,checks=[compact_analysis(c) for c in value.get('checks',[])])
    result['physical_budget_diagnostic_controls']={p.parent.name:receipt(p) for p in (r/'physical-budgets').glob('test-*/results.json')}
    for name in ('runtime','runtime-resolution-v2','runtime-resolution-v3','runtime-budget-v4','runtime-budget-v5','runtime-probe-v6','runtime-probe-v7'):
        for leaf in ('source-build.json','actual-loaded-libraries.json','actual-loaded-libraries-current.json','lineage.json'):
            path=r/name/leaf
            if path.exists():shutil.copy2(path,out/(name+'-'+leaf))
    for name in ('worker-bookings.jsonl','worker-guard-events.jsonl','targeted-ctest.log','start-duration.log','short-coverage-repair.log'):
        if (r/name).exists():shutil.copy2(r/name,out/name)
    gates={};scope_input=dict(verification={},neutral={},responses={})
    for name,path in (('aggregate_negative_controls',r/'aggregate-negative-controls/results.json'),
                      ('actual_initial_binding',r/'actual-initial-regression.json'),
                      ('historical_receipt_bytes',r/'historical-receipt-byte-verification.json'),
                      ('targeted_regressions',a.regression_receipt or r/'targeted-regressions-current.json'),
                      ('forward_interface_tolerance',a.forward_interface_receipt or r/'forward-interface-tolerance-status.json'),
                      ('physical_budget_controls',r/'physical-budget-face-regression.json')):
        scope_input['verification'][name]=receipt(path)
    accounting=receipt(r/'build-probe-side-control.json');scope_input['verification']['physical_budget_accounting']=dict(status='pass' if accounting.get('accounting_unit_exit_code')==0 else 'fail',receipt=accounting)
    current=scope_input['verification']['targeted_regressions']
    if current.get('status')=='pass':
        valid=bool(current.get('source_file_sha256')) and all(sha(r/'source/solvers/IncNavierStokesSolver/Tests/RANS'/name)==digest for name,digest in current['source_file_sha256'].items())
        if not valid:current['status']='fail';current['failure']='tested_script_identity_changed'
    for name,folder,binary in (('current_solver_libraries','runtime','IncNavierStokesSolver'),
                               ('current_diagnostic_libraries','runtime-resolution-v3','RANSPreservationProbe')):
        value=receipt(r/folder/'actual-loaded-libraries-current.json')
        valid=value.get('executable')==str(r/folder/binary) and value.get('executable_sha256')==sha(r/folder/binary)
        valid=valid and bool(value.get('libraries')) and all(Path(p).is_file() and sha(Path(p))==h for p,h in value.get('libraries',{}).items())
        scope_input['verification'][name]=dict(status='pass' if valid else 'fail',receipt=value)
    for path in a.admissions:
        ad=receipt(path);cfg=ad['configuration'];tag=prefix(cfg['model'],cfg['policy'],cfg['wall']);result['admissions'][tag]=ad;attempts=[];qualified=False
        selected=None
        for dt in (.002,.001):
            cohort=tag+'-duration-'+str(dt);control=r/'duration-cohorts'/cohort/'request-manifest.json'
            if not control.exists():continue
            request=load(control);RequestManifest(control.parent,request['request'],True);row=dict(controller_request=request,comparisons={},audits={},numerical_status='not_completed')
            recovery_events=[]
            for folder in (r/'duration-schedules').glob(cohort+'-*'):
                if not (folder/'request-manifest.json').exists() or not (folder/'results.json').exists():continue
                recovery=load(folder/'request-manifest.json');q=recovery['request']
                if q.get('schema')!='UnchangedMethodExecutionRecoveryV1' or q.get('cohort')!=cohort:continue
                require(q['cohort_request_sha256']==request['sha256'],'recovery_controller_identity_mismatch')
                require(sha(q['previous_failure'])==q['previous_failure_sha256'],'recovery_failure_changed')
                RequestManifest(folder,q,True);recovery_events.extend(load(folder/'results.json')['milestones'])
            row['execution_recoveries']=recovery_events
            for T in ad['milestones']:
                comparison=receipt(r/'comparisons'/(cohort+'-'+str(T)+'.json'))
                candidates=[receipt(Path(e['comparison'])) for e in recovery_events if abs(e['time']-T)<1e-13]
                passing=[v for v in candidates if v.get('status')=='pass']
                if comparison.get('status')!='pass' and passing:comparison=dict(passing[-1],previous_receipts=[comparison,*[v for v in candidates if v.get('status')!='pass']])
                row['comparisons'][str(T)]=comparison;audit=receipt(r/'duration-audits'/(cohort+'-'+str(T))/'results.json')
                if audit.get('status')=='pass':
                    aq=load(Path(audit['source']).parent/'request-manifest.json')['request']
                    require(RequestManifest(Path(audit['source']).parent,aq,True).completed(Path(audit['source']).parent),'duration_audit_unsealed_or_changed')
                    check=aggregate(dict(cfg,audit_times=[t for t in ad['milestones'] if t<=T+1e-13]),audit);require(check['status']=='pass','aggregate_receipt_invalid:'+cohort)
                if 'spatial_audits' in audit:
                    audit['spatial_audits']={m:dict(value,checks=[compact_analysis(c) for c in value['checks']]) for m,value in audit['spatial_audits'].items()}
                row['audits'][str(T)]=audit
            inner=receipt(r/'common-inner-qualification'/cohort/'results.json');row['fixed_common_inner_k']=inner
            if inner.get('status')=='pass':
                from qualify_rans_duration_inner_k import reduction
                identity(inner['configuration'],cfg);required_times(inner['checks'],ad['milestones'])
                required_times(inner['short_checks'],(.0016,.016))
                require(inner['common_packet_sha256']==ad['common_packet_sha256'],'inner_common_packet_changed')
                iq=load(Path(inner['source']).parent/'request-manifest.json')['request']
                require(RequestManifest(Path(inner['source']).parent,iq,True).completed(Path(inner['source']).parent),'inner_common_receipt_unsealed')
                for check in inner['checks']:
                    require(reduction(*[check['groups'][m]['inner']['mass_dual'] for m in MEMBERS[:2]])['status'] in ('pass','unresolved_floor'),'common_inner_reduction_failed')
                    for m in MEMBERS[:2]:require(inner['source_audit_sha256'][str(check['time'])][m]==sha(r/'evolved-audits'/(cohort+'-'+str(check['time'])+'-'+m)/'results.json'),'inner_common_source_audit_changed')
                for check in inner['short_checks']:require(reduction(*[check['groups'][m]['inner']['mass_dual'] for m in MEMBERS[:2]])['status'] in ('pass','unresolved_floor'),'short_common_inner_reduction_failed')
            final=row['audits']['0.5'];complete=all(x.get('status')=='pass' for x in row['comparisons'].values()) and final.get('status')=='pass' and inner.get('status')=='pass'
            row['classification']=dict(
                process_and_endpoint_coverage='pass' if all(len(v.get('provenance',[]))==4 for v in row['comparisons'].values()) else 'incomplete',
                field_and_timestep='pass' if all(v.get('status')=='pass' for v in row['comparisons'].values()) else 'failed_or_incomplete',
                spatial_reconstruction='pass' if all(final.get('spatial_audits',{}).get(m,{}).get('status')=='pass' for m in MEMBERS[:2]) else 'failed_or_incomplete',
                dominant_spatial_reduction=final.get('aggregate',{}).get('evolved_reduction_status','not_completed'),
                fixed_common_inner_k=inner.get('status','not_completed'),
                scalar_coverage=final.get('aggregate',{}).get('coverage_status','not_completed'),
                frozen_quadrature=final.get('quadrature',{}).get('status','not_completed'),
                native_constitutive='pass' if all(final.get('constitutive',{}).get(m,{}).get('status')=='pass' for m in MEMBERS[:2]) else 'failed_or_incomplete')
            failed_reductions=[v for v in final.get('evolved_reductions',[]) if v['status']!='pass']
            if failed_reductions:
                first=min(failed_reductions,key=lambda v:v['time']);time=first['time']
                witness=dict(cohort=cohort,configuration=cfg,time=time,gate='evolved_spatial_reduction',
                    reduction=first,aggregate=final.get('aggregate'),members={},
                    meaning='post-step qualification failure; retained states are accepted physical states, no rejected timestep is invented')
                for member in MEMBERS:
                    endpoint=r/'campaigns'/(cohort+'-'+member)/('T-'+format(time,'.4f'))
                    native=load(endpoint/'case.evolution.json');manifest=load(endpoint.parent/'request-manifest.json')
                    require(RequestManifest(endpoint.parent,manifest['request'],True).completed(endpoint),'failure_witness_endpoint_unsealed')
                    name=cohort+'-'+member+'-'+str(time)+'-spatial-gate.xml';shutil.copy2(endpoint/'case.xml',out/name)
                    witness['members'][member]=dict(source=str(endpoint),configuration_file=name,
                        configuration_sha256=sha(endpoint/'case.xml'),native_report_sha256=sha(endpoint/'case.evolution.json'),
                        trajectory_request=manifest,source_fields={k:native[k] for k in (
                            'signature','time','dt','order','state','initial','fixed_reference','initialization_receipt',
                            'native_history','native_history_times','last_stage_history','last_stage_alpha',
                            'last_stage_time','failure','rejected_state')})
                name=cohort+'-spatial-gate-failure.json.gz'
                with gzip.open(out/name,'wt') as f:json.dump(witness,f,separators=(',',':'),allow_nan=False)
                result['failures'].append(dict(cohort=cohort,time=time,gate='evolved_spatial_reduction',witness=name,sha256=sha(out/name)))
            if complete:
                required_times(list(row['comparisons'].values()),(.016,.1,.5))
                selected=dict(configuration=cfg,status='pass',comparisons=[dict(time=v['time'],status=v['status']) for v in row['comparisons'].values()],final_audit_status=final['status'],companions={m:v['scalar_coverage'] for m,v in final['companions'].items()})
            row['numerical_status']='pass' if complete else 'failed_unresolved_or_incomplete';qualified|=complete;result['cohorts'][cohort]=row;attempts.append(cohort)
        gates[tag+'-duration']='pass' if qualified else 'failed_unresolved_or_incomplete'
        result['admissions'][tag]['attempted_cohorts']=attempts
        if cfg['wall']=='Accepted':scope_input['neutral'][cfg['model']]=selected or dict(configuration=cfg,status='failed_unresolved_or_incomplete')
    for folder in sorted((r/'campaigns').iterdir()):
        path=folder/'request-manifest.json'
        if not path.exists():continue
        manifest=load(path);request=manifest['request'];seal=RequestManifest(folder,request,True);row=dict(label=folder.name,request=manifest,initial_metrics=receipt(folder/'initial-metrics.json'),milestones=[])
        for path in sorted(folder.glob('T-*/metrics.json')):
            if not seal.completed(path.parent):continue
            metric=load(path);native=load(path.parent/'case.evolution.json');item={k:v for k,v in metric.items() if k not in ('all_audited_accuracy','spatial_alpha_zero_full','extrema','initialization_receipt')};item['initialization']={k:v for k,v in metric.get('initialization_receipt',{}).items() if k not in ('live_fixed_boundary','effective_initial')};item['extrema']={k:v for k,v in metric.get('extrema',{}).items() if k!='final_state'};item['final_state_audit']=metric.get('final_audit');item['all_state_and_stage']=audit_summaries(native);item['native_report_sha256']=sha(path.parent/'case.evolution.json');item['correction_calls']=sum(x.get('bound_calls',0) for x in native.get('pimple_attempts',[]));item['correction_signed_integral']=item['correction_absolute_integral']=0. if item['correction_calls']==0 else None;row['milestones'].append(item)
            if metric['status']!='pass':
                witness=dict(label=folder.name,metrics=metric,configuration_sha256=sha(path.parent/'case.xml'),native_failure=native.get('failure'),last_accepted_state=native.get('state'),rejected_state=native.get('rejected_state'),native_history=native.get('native_history'),range_events=native.get('physical_range_events'),source=str(path.parent))
                dest=out/(folder.name+'-'+str(metric['target'])+'-failure.json.gz')
                with gzip.open(dest,'wt') as f:json.dump(witness,f)
                shutil.copy2(path.parent/'case.xml',out/(folder.name+'-failure.xml'));result['failures'].append(dict(label=folder.name,target=metric['target'],failure=metric.get('failure'),witness=dest.name))
        result['cases'].append(row)
    for path in sorted((r/'response-cohorts').glob('*/results.json')) if (r/'response-cohorts').exists() else []:
        data=receipt(path);controller=load(path.parent/'request-manifest.json')['request'];response_seal=RequestManifest(path.parent,controller,True);identity(data['configuration'],controller['configuration']);bodies={}
        for milestone in data.get('milestones',[]):
            target=Path(milestone['comparison']);body=receipt(target)
            if body.get('status')=='pass':
                require(response_seal.completed(target.parent),'response_comparison_unsealed');identity(body['configuration'],controller['configuration']);require(abs(body['time']-milestone['time'])<1e-13,'response_comparison_time_mismatch')
                require(body['amplitude']==controller['amplitude'] and body['timesteps']==controller['timesteps'],'response_request_physics_mismatch')
            bodies[str(milestone['time'])]=body
        data['complete_comparisons']=bodies;data['audits']={};rows=[]
        for time,body in bodies.items():
            audit=receipt(r/'response-audits'/(path.parent.name+'-'+time)/'results.json')
            if audit.get('status')=='pass':
                aq=load(Path(audit['source']).parent/'request-manifest.json')['request'];identity(aq['configuration'],controller['configuration']);require(aq['admission_sha256']==controller['admission_sha256'],'response_audit_admission_mismatch');require(abs(aq['time']-float(time))<1e-13,'response_audit_time_mismatch');require(RequestManifest(Path(audit['source']).parent,aq,True).completed(Path(audit['source']).parent),'response_audit_unsealed')
            if 'spatial' in audit:audit['spatial']={m:dict(v,checks=[compact_analysis(c) for c in v['checks']]) for m,v in audit['spatial'].items()}
            data['audits'][time]=audit
            supplements={p.parent.name:receipt(p) for p in (r/'response-audits').glob(path.parent.name+'-'+time+'-*/results.json')}
            for value in supplements.values():
                aq=load(Path(value['source']).parent/'request-manifest.json')['request'];identity(aq['configuration'],controller['configuration'])
                require(aq['trajectory_requests']==load(r/'response-audits'/(path.parent.name+'-'+time)/'request-manifest.json')['request']['trajectory_requests'],'response_supplement_trajectory_mismatch')
                require(aq['previous_failure']['sha256']==sha(r/'response-audits'/(path.parent.name+'-'+time)/'results.json'),'response_supplement_failure_changed')
                require(RequestManifest(Path(value['source']).parent,aq,True).completed(Path(value['source']).parent),'response_supplement_unsealed')
                if 'spatial' in value:value['spatial']={m:dict(v,checks=[compact_analysis(c) for c in v['checks']]) for m,v in value['spatial'].items()}
            data.setdefault('complete_audit_supplements',{})[time]=supplements
            if audit.get('status')=='fail':
                name=path.parent.name+'-'+time+'-response-audit-failure.json.gz'
                witness=dict(configuration=controller['configuration'],time=float(time),gate='response_native_spatial_audit',audit=audit,members={},
                    meaning='post-step read-only qualification failure in accepted states; no rejected timestep is fabricated')
                for m in MEMBERS:
                    endpoint=r/'campaigns'/(path.parent.name+'-'+m)/('T-'+format(float(time),'.4f'));native=load(endpoint/'case.evolution.json')
                    q=load(endpoint.parent/'request-manifest.json');require(RequestManifest(endpoint.parent,q['request'],True).completed(endpoint),'response_failure_endpoint_unsealed')
                    config=path.parent.name+'-'+m+'-'+time+'-audit-gate.xml';shutil.copy2(endpoint/'case.xml',out/config)
                    witness['members'][m]=dict(configuration_file=config,configuration_sha256=sha(endpoint/'case.xml'),native_report_sha256=sha(endpoint/'case.evolution.json'),request=q,
                        source_fields={k:native[k] for k in ('signature','time','dt','order','state','initial','fixed_reference','initialization_receipt','native_history','native_history_times','last_stage_history','last_stage_alpha','last_stage_time','failure','rejected_state')})
                with gzip.open(out/name,'wt') as f:json.dump(witness,f,separators=(',',':'),allow_nan=False)
                result['failures'].append(dict(cohort=path.parent.name,time=float(time),gate='response_native_spatial_audit',witness=name,sha256=sha(out/name)))
            budgets={q.parent.name:receipt(q) for q in (r/'physical-budgets').glob(path.parent.name+'-*'+time+'*/results.json')}
            data.setdefault('physical_budgets',{})[time]=budgets
            required_budgets={path.parent.name+'-'+kind+'-'+m+'-'+time+'-Q'+str(q)+suffix for kind in ('perturbed','baseline') for m in MEMBERS for q,suffix in (((1,''),(2,''),(2,'-F32')) if m in MEMBERS[:2] else ((2,'-F32'),))}
            coverage=required_budgets.issubset(budgets) and all(budgets[n]['status']=='pass' for n in required_budgets)
            if coverage:
                for kind in ('perturbed','baseline'):
                    head=path.parent.name if kind=='perturbed' else prefix(**{k:controller['configuration'][k] for k in ('model','policy','wall')})+'-duration-'+str(controller['timesteps'][0])
                    for m in MEMBERS:
                        report=r/'campaigns'/(head+'-'+m)/('T-'+format(float(time),'.4f'))/'case.evolution.json'
                        for q,suffix in (((1,''),(2,''),(2,'-F32')) if m in MEMBERS[:2] else ((2,'-F32'),)):
                            value=budgets[path.parent.name+'-'+kind+'-'+m+'-'+time+'-Q'+str(q)+suffix]
                            identity(value['configuration'],controller['configuration']);require(value['source_report_sha256']==sha(report),'physical_budget_native_report_mismatch');require(abs(value['time']-float(time))<1e-13,'physical_budget_horizon_mismatch')
                for name in required_budgets:
                    folder=Path(budgets[name]['source']).parent;q=load(folder/'request-manifest.json')['request'];require(RequestManifest(folder,q,True).completed(folder),'physical_budget_unsealed');identity(q['configuration'],controller['configuration'])
            budget_analysis=receipt(r/'physical-budget-analyses'/(path.parent.name+'-'+time)/'results.json');data.setdefault('physical_budget_analysis',{})[time]=budget_analysis
            coverage=coverage and budget_analysis.get('execution_status')=='completed' and budget_analysis.get('reconstruction_status')=='pass'
            data.setdefault('physical_budget_coverage',{})[time]=dict(status='reconstructed' if coverage else 'not_completed',required=sorted(required_budgets),observed=sorted(budgets))
            rows.append(dict(time=float(time),status=body.get('status','not_completed'),audit_status=audit.get('status','not_completed'),physical_budget_status='reconstructed' if coverage else 'not_completed',companions={m:dict(last_time=float(time)) for m in MEMBERS} if len(body.get('provenance',[]))==8 else {}))
        result['responses'][path.parent.name]=data
        if controller.get('amplitude')==2e-5:
            value=dict(configuration=data['configuration'],status=data['status'],comparisons=rows,extended_to_point_one=any(v['time']==.1 for v in rows))
            if data['configuration']['wall']=='Accepted':scope_input['responses'][data['configuration']['model']]=value
            elif data['configuration']['wall']=='EP':
                value['wall_comparison_status']=receipt(r/'wall-response-comparison/results.json').get('status','not_completed');scope_input['ep_control']=value
    result['half_amplitude']=receipt(r/'half-amplitude-comparison/results.json');scope_input['half_amplitude']=result['half_amplitude']
    result['wall_response_comparison']=receipt(r/'wall-response-comparison/results.json')
    result['historical_epsilon_orders']='eight observations remain unresolved under their immutable original criterion; new pairs bound uncertainty and do not measure order'
    result['gates']=gates;result['required_scope_input']=scope_input;result['required_scope_classification']=classify(scope_input);result['status']=result['required_scope_classification']['status']
    result['scope']='serial flat2D prescribed-reference patch; unsupported geometry/3D/MPI/actions unchanged; response and independent physical-flux scopes reported separately'
    complete=out/'complete-results.json.gz'
    with gzip.open(complete,'wt',compresslevel=9) as f:json.dump(result,f,separators=(',',':'),allow_nan=False)
    light={k:v for k,v in result.items() if k not in ('cases','cohorts','admissions','responses')};light.update(complete_results=dict(file=complete.name,sha256=sha(complete)),admissions={k:dict(status=v['status'],configuration=v['configuration'],source=v['source'],source_sha256=v['source_sha256']) for k,v in result['admissions'].items()},cohorts={k:dict(numerical_status=v['numerical_status'],classification=v['classification'],comparison_statuses={t:x['status'] for t,x in v['comparisons'].items()},audit_statuses={t:x['status'] for t,x in v['audits'].items()}) for k,v in result['cohorts'].items()},response_statuses={k:v['status'] for k,v in result['responses'].items()})
    write_json(out/'summary.json',light);print(result['status']);return 0 if result['status']=='pass' else 1


def main():
    # Identity rejection must leave a structured receipt too. Never replace a
    # completed destination, including on a malformed reuse request.
    small=argparse.ArgumentParser(add_help=False);small.add_argument('--output',type=Path);small.add_argument('--root',type=Path);small.add_argument('--admissions',type=Path,nargs='+')
    args,_=small.parse_known_args();existed=bool(args.output and args.output.exists())
    try:return collect()
    except Exception as error:
        import time
        body=dict(status='fail',execution_status='failed_before_complete_export',failure=repr(error),
            command=sys.argv,root=str(args.root),existing_destination_preserved=existed,
            admissions={str(p):dict(sha256=sha(p),content_utf8=p.read_text(errors='replace')) for p in (args.admissions or []) if p.is_file()})
        if args.output:
            if existed:path=args.output.with_name(args.output.name+'-rejected-'+str(time.time_ns())+'.json')
            else:
                args.output.mkdir(parents=True,exist_ok=True);path=args.output/'export-failure.json'
            write_json(path,body);print(str(path))
        else:print(json.dumps(body,allow_nan=False))
        return 1


if __name__=='__main__':raise SystemExit(main())
