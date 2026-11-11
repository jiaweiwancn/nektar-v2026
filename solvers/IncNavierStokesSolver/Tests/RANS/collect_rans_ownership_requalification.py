#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Supplement historical evidence with complete corrected-operator receipts."""
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
from qualify_rans_preservation import load, replay_request
from rans_request_manifest import sha, RequestManifest
from rans_ownership_contract import SCHEMA, OPERATOR, classify_scope
from rans_qualification import write_json
from rans_duration_contract import identity, prefix, require, required_times, MEMBERS
from collect_rans_balanced_spatial import audit_summaries


def short_matrix_coverage(root, short):
    """A pass label cannot replace any of the twelve required receipt bodies."""
    configurations=[dict(model=m,policy='Natural',wall=w)
                    for m in ('A93','A95','E94E95') for w in ('Accepted','EP')]
    errors=[];receipts=[];rows=(short or {}).get('configurations',{})
    if set(rows)!={prefix(**cfg) for cfg in configurations}:errors.append('short_configuration_coverage')
    for cfg in configurations:
        tag=prefix(**cfg)
        try:
            row=rows[tag];identity(row['configuration'],cfg)
            require(row['numerical_operator']==OPERATOR,'short_operator_mismatch')
            required_times(row['milestones'],(.0016,.016))
            require(len(row['milestones'])==2,'short_time_coverage')
            require(row['status']=='pass','short_configuration_failed')
            for milestone in row['milestones']:
                path=Path(milestone['receipt']).resolve();path.relative_to(root.resolve())
                require(sha(path)==milestone['sha256'],'short_receipt_hash_mismatch')
                body=load(path);identity(body['configuration'],cfg)
                require(body['numerical_operator']==OPERATOR,'short_receipt_operator_mismatch')
                require(abs(body['time']-milestone['time'])<=1e-13,'short_receipt_time_mismatch')
                require(body['status']==milestone['status']=='pass','short_receipt_failed')
                receipts.append(dict(configuration=cfg,time=body['time'],path=str(path),sha256=sha(path)))
        except (KeyError,ValueError,OSError) as error:errors.append(tag+':'+str(error))
    return dict(status='pass' if not errors and len(receipts)==12 and (short or {}).get('status')=='pass' else 'failed_unresolved_or_incomplete',
                required_times=[.0016,.016],required_receipts=12,verified_receipts=receipts,failures=errors)


def collect(root,output,previous_attempts):
    r=root.resolve();output.mkdir(parents=True,exist_ok=False);bodies={};index=[];costs=[];states=[]
    def include(path):
        if not path.exists():return None
        owner=path.parent
        if (path.parent/'completed-artifacts.json').exists():
            while owner!=r and not (owner/'request-manifest.json').exists():owner=owner.parent
            request=owner/'request-manifest.json'
            require(request.exists(),'completed_receipt_has_no_request_owner:'+str(path))
            if not RequestManifest(owner,load(request)['request'],True).completed(path.parent):
                raise ValueError('unsealed completed receipt: '+str(path))
        value=load(path);key=str(path.relative_to(r));bodies[key]=value
        index.append(dict(path=str(path),sha256=sha(path)));return value
    geometry=include(r/'geometric-regressions-v1/results.json')
    compat=include(r/'compatibility-v1/results.json');prescribed=include(r/'selected-prescribed-results.json')
    prescribed_supplement=include(r/'selected-prescribed-results-v2.json')
    if prescribed_supplement:prescribed=prescribed_supplement
    include(r/'build/a93-a95-e94e95-kernels.json')
    include(r/'ownership-short-qualification/progress.json')
    short=include(r/'ownership-short-qualification/results.json')
    supplement=include(r/'ownership-short-qualification-v2/results.json')
    if supplement:short=supplement
    for glob in ('short-audits/*/results.json','short-audits-v2/*/results.json','initial-gate-*.json','admissions/*.json',
        'duration-schedules/*/results.json','duration-audits/*/results.json','comparisons/*.json',
        'common-inner-qualification/*/results.json','response-cohorts/*/results.json',
        'response-cohorts/*/T-*/results.json','response-audits/*/results.json',
        'physical-budget-analyses/*/results.json','physical-budget-accuracy/*/results.json',
        'physical-budget-bindings/*/results.json','physical-budget-diagnosis/*/results.json',
        'physical-postprocessing/*.json','common-wall-response/results.json',
        'response-quadrature-*/results.json','duration-quadrature-*/results.json',
        'corrective-duration-response/*.json','corrective-duration-response-v2/*.json',
        'duration-residual-diagnosis/*/results.json',
        'selected-frozen-operator-comparison/*/results.json','actual-loaded-maps/*.json',
        'budget-axis-controls/*/results.json','budget-axis-controls/*/hypothesis.json',
        'budget-axis-controls/*/mesh.json','budget-axis-controls/*/unsupported-receipt.json',
        'budget-axis-controls/*/execution/results.json'):
        for path in sorted(r.glob(glob)):include(path)
    for path in sorted((r/'campaigns').glob('*/T-*/metrics.json')):
        m=load(path);index.append(dict(path=str(path),sha256=sha(path)))
        for name in ('case.xml','case.evolution.json','case.evolution.state',
                     'case.evolution.state.json','completed-artifacts.json','endpoint-samples.json'):
            f=path.parent/name
            if f.exists():index.append(dict(path=str(f),sha256=sha(f)))
        costs.append(dict(campaign=path.parents[1].name,time=m.get('time'),target=m.get('target'),
            status=m['status'],wall_seconds=m['wall_seconds'],process_cpu_seconds=m.get('process_cpu_seconds'),
            history=m.get('history_identity'),configuration_sha256=m['configuration_sha256'],
            initialization_receipt=m.get('initialization_receipt'),final_state_sha256=m.get('final_state_sha256')))
        native=load(path.parent/'case.evolution.json')
        domain_keys=('minimum_inner_k','minimum_outer_k','minimum_outer_epsilon','minimum_sampled_K',
            'minimum_projected_trace_epsilon','scalar_minimum_points_inner_k_outer_k_epsilon',
            'sampled_covariance_minima','covariance_minimum_points','covariance_minimum_sides',
            'covariance_minimum_geometries','minimum_A95_radicand','A95_radicand_minimum_point',
            'A95_radicand_side','A95_radicand_geometry','minimum_required_budget_ratio',
            'required_cutoff_count','required_cutoff_admissible','cutoff_failure_reason',
            'cutoff_limiting_point','cutoff_range','maximum_budget_relative_defect',
            'active_correction_profiles','inactive_baseline_profiles','fallback_profiles')
        domains=[dict(time=attempt['time'],stage=stage,**{k:d[k] for k in domain_keys if k in d})
                 for attempt in native.get('pimple_attempts',[]) for stage,d in attempt.items()
                 if stage.endswith('_domain') and isinstance(d,dict)]
        calls=sum(x.get('bound_calls',0) for x in native.get('pimple_attempts',[]))
        states.append(dict(campaign=path.parents[1].name,time=m.get('time'),status=m['status'],
            native_report_sha256=sha(path.parent/'case.evolution.json'),extrema=m.get('extrema'),
            final_state_audit=m.get('final_audit'),all_state_and_stage=audit_summaries(native),
            audited_stage_domains=domains,reference_error=m.get('reference_error'),
            actual_initial_drift=m.get('change_from_initial'),spatial_alpha_zero=m.get('spatial_alpha_zero'),
            source_and_constraint_maxima=m.get('source_and_constraint_maxima'),
            corrections=dict(calls=calls,signed_integral=0. if calls==0 else None,
                             absolute_integral=0. if calls==0 else None)))
    for path in sorted((r/'campaigns').glob('*/request-manifest.json')):include(path)
    for path in sorted(r.rglob('probe-process-cost.json')):
        # The source/build are symlinked to the owned isolated build area;
        # rglob does not follow those directories or any unrelated campaigns.
        costs.append(dict(diagnostic=str(path),**load(path)));index.append(dict(path=str(path),sha256=sha(path)))
    for name in ('runtime','runtime-resolution-v3','runtime-budget-v5'):
        include(r/name/'source-build.json')
    for name in ('targeted-native-executions.json','required-scope-cli-receipt.json','controller-launches.json',
                 'worker-guard-status.json','debugger-resource-incident.json','owned-orphan-MPI-cleanup.json',
                 'idle-watcher-resource-cleanup.json','independent-admission-launch.json',
                 'build-completion.json','source-identities.json','native-exclusive-costs.json',
                 'final-resource-census.json','historical-evidence-integrity.json',
                 'completion-wrapper-transition.json','legacy-probe-dependency-identity.json',
                 'short-coverage-negative-control.json','budget-cost-instrumentation.json',
                 'budget-postprocessing-v2-launch.json','budget-axis-controller-launch.json',
                 'budget-axis-control-results.json','A93-half-spatial-execution.json',
                 'A93-corrected-signed-load-execution.json','environment-and-cost-summary.json',
                 'completed-request-recheck.json'):
        include(r/name)
    for name in ('worker-bookings.jsonl','worker-guard-events.jsonl','resource-priority-events.jsonl'):
        path=r/name
        if path.exists():
            with path.open('rb') as source,gzip.open(output/(name+'.gz'),'wb',compresslevel=9) as dest:
                shutil.copyfileobj(source,dest)
            index.append(dict(path=str(path),sha256=sha(path)))
    schedules={}
    for model in ('A93','A95'):
        p=r/'corrective-duration-response-v2'/(model+'-N-Accepted-neutral.json')
        if not p.exists():p=r/'corrective-duration-response'/(model+'-N-Accepted-neutral.json')
        v=include(p)
        schedules[model]=dict(status=v['status'] if v else 'not_attempted',source=str(p),body=v)
    response=r/'corrective-duration-response/response-results.json';v=include(response)
    field_rows=[body for key,body in bodies.items() if key.startswith('response-cohorts/') and '/T-' in key]
    weak_rows=[body for key,body in bodies.items() if key.startswith('response-audits/')]
    budget_rows=[body for key,body in bodies.items() if key.startswith('physical-budget-bindings/')]
    simple=lambda v:dict(status=v.get('status','not_attempted') if v else 'not_attempted')
    response_coverage={(v.get('configuration',{}).get('wall'),v.get('amplitude')) for v in field_rows}
    response_complete=(len(field_rows)==3 and response_coverage=={('Accepted',2e-5),('Accepted',1e-5),('EP',2e-5)} and
        all(v.get('configuration',{}).get('model')=='A95' and v.get('configuration',{}).get('policy')=='Natural' and
            abs(v.get('time',-1)-.016)<=1e-13 and len(v.get('fields',{}))==17 and
            v.get('core_pressure_velocity_nonlinear_production_response_resolved') for v in field_rows))
    half=bodies.get('physical-postprocessing/half-amplitude.json',{})
    common_wall=bodies.get('common-wall-response/results.json',{})
    weak_complete=response_complete and len(weak_rows)==3
    for key,body in list(bodies.items()):
        if not key.startswith('response-cohorts/') or '/T-' not in key:continue
        name=key.split('/')[1];audit=bodies.get('response-audits/'+name+'-0.016/results.json',{})
        weak_complete &= (audit.get('status')=='pass' and audit.get('configuration')==body.get('configuration') and
            abs(audit.get('time',-1)-.016)<=1e-13 and audit.get('independent_coverage')=='complete' and
            set(audit.get('scalars',{}))==set(MEMBERS) and all(v['coverage']['status']=='pass' and
                v['coverage']['last_time']>=.016-1e-13 for v in audit.get('scalars',{}).values()))
    half_complete=(half.get('status')=='pass' and len(half.get('checks',[]))==68 and
                   all(math.isfinite(v['limit']) and v['limit']>=0 and
                       0<=v['departure_from_linear_scaling_RMS']<=v['limit'] for v in half.get('checks',[])))
    wall_complete=(common_wall.get('status')=='pass' and abs(common_wall.get('time',-1)-.016)<=1e-13 and
                   set(common_wall.get('fields',{}))=={'Accepted','EP'} and
                   all(len(v)==17 for v in common_wall.get('fields',{}).values()) and
                   len(common_wall.get('wall_differences',{}))==17 and all(
                       all(math.isfinite(v[k]) and v[k]>=0 for k in ('spatial_difference','temporal_difference','spatial_limit','temporal_limit')) and
                       v['spatial_difference']<=v['spatial_limit'] and v['temporal_difference']<=v['temporal_limit']
                       for fields in common_wall.get('fields',{}).values() for v in fields.values()))
    gates=dict(ownership_geometry=simple(geometry),history_cache_CLS=simple(compat),natural_short_matrix=short_matrix_coverage(r,short),
        selected_prescribed=simple(prescribed),A93_neutral=schedules['A93'],A95_neutral=schedules['A95'],
        A95_response_fields=dict(status='pass' if response_complete and all(x['status']=='pass' for x in field_rows) and
            half_complete and wall_complete else 'failed_unresolved_or_incomplete',
            half_amplitude_status='pass' if half_complete else 'incomplete_or_failed',
            common_wall_status='pass' if wall_complete else 'incomplete_or_failed'),
        A95_response_native_weak=dict(status='pass' if weak_complete else 'failed_unresolved_or_incomplete'),
        physical_budget_reconstruction=dict(status='pass' if len(budget_rows)==1 and budget_rows[0]['reconstruction_status']=='pass' else 'not_qualified'),
        physical_budget_accuracy=dict(status=budget_rows[0]['accuracy_status'] if len(budget_rows)==1 else 'not_qualified'))
    scope=dict(schema=SCHEMA,numerical_operator=OPERATOR,response_max_time=.016,out_of_scope=dict(A93_response=True),gates=gates)
    write_json(output/'required-scope-input.json',scope);decision=classify_scope(scope)
    write_json(output/'decision.json',decision);write_json(output/'costs-and-input-identities.json',dict(records=costs))
    with gzip.open(output/'all-state-stage-receipts.json.gz','wt',compresslevel=9) as stream:
        json.dump(states,stream,separators=(',',':'),allow_nan=False)
    with gzip.open(output/'complete-results.json.gz','wt',compresslevel=9) as stream:json.dump(bodies,stream,separators=(',',':'),allow_nan=False)
    prior={}
    for path in previous_attempts:
        entries=[]
        for name in ('retired-attempt.json','short-execution-stop.json','debugger-resource-incident.json',
                     'geometric-regressions-v1/results.json','compatibility-v1/results.json',
                     'ownership-short-qualification/results.json','ownership-short-qualification/progress.json'):
            p=path/name
            if p.exists():entries.append(dict(path=str(p),sha256=sha(p),body=load(p)))
        prior[str(path)]=entries
    write_json(output/'preserved-execution-failures.json',prior)
    witnesses=output/'failure-witnesses';witnesses.mkdir();failures=[]
    for key,body in bodies.items():
        if not key.startswith('duration-audits/'):continue
        failed=[v for v in body.get('evolved_reductions',[]) if v['status']!='pass']
        if not failed:continue
        first=min(failed,key=lambda v:v['time']);t=first['time'];cohort=key.split('/')[1].removesuffix('-'+str(t))
        witness=dict(configuration=body['configuration'],cohort=cohort,time=t,numerical_operator=OPERATOR,
            gate='evolved_spatial_reduction',reduction=first,aggregate=body['aggregate'],members={},
            meaning='post-step qualification failure in accepted physical states; no rejected timestep invented')
        for member in MEMBERS:
            folder=r/'campaigns'/(cohort+'-'+member);endpoint=folder/('T-'+format(t,'.4f'))
            q=load(folder/'request-manifest.json');require(RequestManifest(folder,q['request'],True).completed(endpoint),'failure_endpoint_unsealed')
            native=load(endpoint/'case.evolution.json');base=cohort+'-'+member
            shutil.copy2(endpoint/'case.xml',witnesses/(base+'.xml'))
            request=replay_request(load(folder/'initial-packet.json'),[],state=native['state'],time=native['time'])
            request.update(fixed_boundary=native['initialization_receipt']['live_fixed_boundary'],export_basis=True,export_resolution=True)
            with gzip.open(witnesses/(base+'-replay.json.gz'),'wt',compresslevel=9) as stream:json.dump(request,stream,separators=(',',':'),allow_nan=False)
            checkpoint=Path(native['checkpoint'])/'State.bin'
            if checkpoint.exists():
                with checkpoint.open('rb') as source,gzip.open(witnesses/(base+'-State.bin.gz'),'wb',compresslevel=9) as dest:shutil.copyfileobj(source,dest)
                index.append(dict(path=str(checkpoint),sha256=sha(checkpoint)))
            witness['members'][member]=dict(source=str(endpoint),configuration_sha256=sha(endpoint/'case.xml'),
                native_report_sha256=sha(endpoint/'case.evolution.json'),trajectory_request=q,
                source_fields={k:native[k] for k in ('signature','numerical_operator','time','dt','order','state','initial',
                    'fixed_reference','initialization_receipt','native_history','native_history_times','last_stage_history',
                    'last_stage_alpha','last_stage_time','failure','rejected_state') if k in native})
        dest=witnesses/(cohort+'-qualification-failure.json.gz')
        with gzip.open(dest,'wt',compresslevel=9) as stream:json.dump(witness,stream,separators=(',',':'),allow_nan=False)
        failures.append(dict(gate=witness['gate'],time=t,witness=dest.name))
    budget_requests=[*(r/'physical-budgets').glob('*-0.016-Q2-F32/request.json'),
                     *(r/'physical-budgets').glob('axis-*-Q2-F32/request.json')]
    for path in sorted(budget_requests):
        # The selected required cohort's eight primary/half/baseline requests
        # retain consecutive accepted and frozen-stage vectors for replay.
        if not path.parent.name.startswith(('A95-N-Accepted-response-2e-05-dt-','axis-')):continue
        base=path.parent.name;shutil.copy2(path.parent/'case.xml',witnesses/(base+'.xml'))
        with path.open('rb') as source,gzip.open(witnesses/(base+'-request.json.gz'),'wb',compresslevel=9) as dest:shutil.copyfileobj(source,dest)
        index.append(dict(path=str(path),sha256=sha(path)))
    fixtures=output/'geometry-and-guard-fixtures';fixtures.mkdir()
    for source in [*(r/'geometric-regressions-v1').glob('*/case.xml'),
                   *(r/'budget-axis-controls').glob('*/native-guard-witness/case.xml')]:
        base=source.parent.name
        if base=='native-guard-witness':base=source.parent.parent.name+'-'+base
        dest=fixtures/base;dest.mkdir()
        for name in ('case.xml','request.json','probe.log','packet.json','immutable-initialization-intent.json'):
            p=source.parent/name
            if not p.exists():continue
            if name=='request.json':
                with p.open('rb') as stream,gzip.open(dest/(name+'.gz'),'wb',compresslevel=9) as target:shutil.copyfileobj(stream,target)
            else:shutil.copy2(p,dest/name)
            index.append(dict(path=str(p),sha256=sha(p)))
    write_json(output/'qualification-failure-index.json',dict(failures=failures,
        physical_budget_accuracy='failed binding and replay requests retained separately from successful reconstruction',
        unsupported_control='StreamwiseP6V1 stopped by the existing equal-order tensor-basis guard before any evolution'))
    controllers=output/'linux-controllers';controllers.mkdir()
    for path in sorted(r.glob('*.py')):
        shutil.copy2(path,controllers/path.name);index.append(dict(path=str(path),sha256=sha(path)))
    write_json(output/'artifact-index.json',dict(authoritative_root=str(r),artifacts=index,
        note='complete retained heavy native histories/packets remain on Linux; compact bodies preserve all classifications and source hashes'))
    sums=[]
    for path in sorted(output.rglob('*')):
        if path.is_file():sums.append(sha(path)+'  '+str(path.relative_to(output)))
    (output/'SHA256SUMS').write_text('\n'.join(sums)+'\n');return decision


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--previous-attempt',type=Path,action='append',default=[])
    a=p.parse_args();result=collect(a.root,a.output,a.previous_attempt);print(result['status']);return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
