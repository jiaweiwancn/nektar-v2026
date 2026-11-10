#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Complete configuration/time-aware export; required failed gates exit nonzero."""
import argparse,gzip,json,os,platform,shutil,subprocess,sys
from pathlib import Path
from rans_qualification import THREADS,write_json
os.environ.update(THREADS)
from qualify_rans_preservation import load
from rans_request_manifest import RequestManifest,sha
from rans_duration_contract import require,required_times,prefix,MEMBERS,aggregate
from collect_rans_balanced_spatial import compact_analysis,audit_summaries
from rans_nonlinear_evidence_contract import classify


def receipt(path):
    if not path.exists():return dict(status='not_completed',source=str(path))
    return dict(load(path),source=str(path),source_sha256=sha(path))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--admissions',type=Path,nargs='+',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=a.root.resolve();out=a.output;out.mkdir(parents=True,exist_ok=False);result=dict(schema='NonlinearDurationResponseEvidenceV1',starting_revision='dfd94161fd9d97f7f1a8ce2966e9cc9461b6d686',source_report=subprocess.check_output(['git','-C',str(r/'source'),'rev-parse','HEAD'],text=True).strip(),environment=dict(platform=platform.platform(),python=sys.version,threads=THREADS,ranks=1,worker_limit=72),admissions={},cohorts={},cases=[],responses={},failures=[],verification={})
    for name in ('workspace.json','historical-receipt-byte-verification.json','actual-initial-regression.json','resource-census-exception.json','worker-guard-status.json'):
        result['verification'][name]=receipt(r/name)
    result['verification']['aggregate-negative-controls']=receipt(r/'aggregate-negative-controls/results.json')
    result['verification']['complete-export-negative-controls']=receipt(r/'complete-export-negative-controls/results.json')
    result['execution_schedules']=[]
    for folder in sorted((r/'duration-schedules').glob('*')):
        if not folder.is_dir():continue
        row=dict(source=str(folder),result=receipt(folder/'results.json'),request=receipt(folder/'request-manifest.json'),logs=[])
        for path in sorted(folder.glob('*.log')):
            text=path.read_text(errors='replace')
            row['logs'].append(dict(source=str(path),sha256=sha(path),tail=text[-5000:]))
        result['execution_schedules'].append(row)
    result['common_scalar_partitions']={p.parent.name:receipt(p) for p in (r/'common-scalar-partitions').glob('*/results.json')}
    result['physical_budget_diagnostic_controls']={p.parent.name:receipt(p) for p in (r/'physical-budgets').glob('test-*/results.json')}
    for name in ('runtime','runtime-resolution-v2','runtime-resolution-v3','runtime-budget-v4'):
        for leaf in ('source-build.json','actual-loaded-libraries.json','actual-loaded-libraries-current.json','lineage.json'):
            path=r/name/leaf
            if path.exists():shutil.copy2(path,out/(name+'-'+leaf))
    for name in ('worker-bookings.jsonl','worker-guard-events.jsonl','targeted-ctest.log','start-duration.log','short-coverage-repair.log'):
        if (r/name).exists():shutil.copy2(r/name,out/name)
    gates={};scope_input=dict(verification={},neutral={},responses={})
    for name,path in (('aggregate_negative_controls',r/'aggregate-negative-controls/results.json'),
                      ('actual_initial_binding',r/'actual-initial-regression.json'),
                      ('historical_receipt_bytes',r/'historical-receipt-byte-verification.json'),
                      ('targeted_regressions',r/'targeted-regressions.json')):
        scope_input['verification'][name]=receipt(path)
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
                    check=aggregate(dict(cfg,audit_times=[t for t in ad['milestones'] if t<=T+1e-13]),audit);require(check['status']=='pass','aggregate_receipt_invalid:'+cohort)
                if 'spatial_audits' in audit:
                    audit['spatial_audits']={m:dict(value,checks=[compact_analysis(c) for c in value['checks']]) for m,value in audit['spatial_audits'].items()}
                row['audits'][str(T)]=audit
            final=row['audits']['0.5'];complete=all(x.get('status')=='pass' for x in row['comparisons'].values()) and final.get('status')=='pass'
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
        data=receipt(path);bodies={}
        for milestone in data.get('milestones',[]):
            target=Path(milestone['comparison']);bodies[str(milestone['time'])]=receipt(target)
        data['complete_comparisons']=bodies;data['audits']={};rows=[]
        for time,body in bodies.items():
            audit=receipt(r/'response-audits'/(path.parent.name+'-'+time)/'results.json')
            if 'spatial' in audit:audit['spatial']={m:dict(v,checks=[compact_analysis(c) for c in v['checks']]) for m,v in audit['spatial'].items()}
            data['audits'][time]=audit
            budgets={q.parent.name:receipt(q) for q in (r/'physical-budgets').glob(path.parent.name+'-*'+time+'*/results.json')}
            data.setdefault('physical_budgets',{})[time]=budgets
            rows.append(dict(time=float(time),status=body.get('status','not_completed'),audit_status=audit.get('status','not_completed'),physical_budget_status='reconstructed' if budgets and all(v['status']=='pass' for v in budgets.values()) else 'not_completed',companions={m:dict(last_time=float(time)) for m in MEMBERS} if len(body.get('provenance',[]))==8 else {}))
        result['responses'][path.parent.name]=data
        controller=load(path.parent/'request-manifest.json')['request']
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
    light={k:v for k,v in result.items() if k not in ('cases','cohorts','admissions','responses')};light.update(complete_results=dict(file=complete.name,sha256=sha(complete)),admissions={k:dict(status=v['status'],configuration=v['configuration'],source=v['source'],source_sha256=v['source_sha256']) for k,v in result['admissions'].items()},cohorts={k:dict(numerical_status=v['numerical_status'],comparison_statuses={t:x['status'] for t,x in v['comparisons'].items()},audit_statuses={t:x['status'] for t,x in v['audits'].items()}) for k,v in result['cohorts'].items()},response_statuses={k:v['status'] for k,v in result['responses'].items()})
    write_json(out/'summary.json',light);print(result['status']);return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
