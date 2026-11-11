#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Supplement historical evidence with complete corrected-operator receipts."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
from qualify_rans_preservation import load
from rans_request_manifest import sha, RequestManifest
from rans_ownership_contract import SCHEMA, OPERATOR, classify_scope
from rans_qualification import write_json


def collect(root,output,previous_attempts):
    r=root.resolve();output.mkdir(parents=True,exist_ok=False);bodies={};index=[];costs=[]
    def include(path):
        if not path.exists():return None
        request=path.parent/'request-manifest.json'
        if request.exists() and (path.parent/'completed-artifacts.json').exists():
            if not RequestManifest(path.parent,load(request)['request'],True).completed(path.parent):
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
        'corrective-duration-response/*.json','duration-residual-diagnosis/*/results.json',
        'selected-frozen-operator-comparison/*/results.json'):
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
                 'final-resource-census.json','historical-evidence-integrity.json'):
        include(r/name)
    for name in ('worker-bookings.jsonl','worker-guard-events.jsonl','resource-priority-events.jsonl'):
        path=r/name
        if path.exists():
            with path.open('rb') as source,gzip.open(output/(name+'.gz'),'wb',compresslevel=9) as dest:
                shutil.copyfileobj(source,dest)
            index.append(dict(path=str(path),sha256=sha(path)))
    schedules={}
    for model in ('A93','A95'):
        p=r/'corrective-duration-response'/(model+'-N-Accepted-neutral.json');v=include(p)
        schedules[model]=dict(status=v['status'] if v else 'not_attempted',source=str(p),body=v)
    response=r/'corrective-duration-response/response-results.json';v=include(response)
    field_rows=[body for key,body in bodies.items() if key.startswith('response-cohorts/') and '/T-' in key]
    weak_rows=[body for key,body in bodies.items() if key.startswith('response-audits/')]
    budget_rows=[body for key,body in bodies.items() if key.startswith('physical-budget-bindings/')]
    simple=lambda v:dict(status=v.get('status','not_attempted') if v else 'not_attempted')
    response_coverage={(v.get('configuration',{}).get('wall'),v.get('amplitude')) for v in field_rows}
    response_complete=(len(field_rows)==3 and response_coverage=={('Accepted',2e-5),('Accepted',1e-5),('EP',2e-5)} and
        all(v.get('configuration',{}).get('model')=='A95' and v.get('configuration',{}).get('policy')=='Natural' and
            abs(v.get('time',-1)-.016)<=1e-13 and len(v.get('fields',{}))==17 for v in field_rows))
    half=bodies.get('physical-postprocessing/half-amplitude.json',{})
    gates=dict(ownership_geometry=simple(geometry),history_cache_CLS=simple(compat),natural_short_matrix=simple(short),
        selected_prescribed=simple(prescribed),A93_neutral=schedules['A93'],A95_neutral=schedules['A95'],
        A95_response_fields=dict(status='pass' if response_complete and all(x['status']=='pass' for x in field_rows) and
            half.get('status')=='pass' and len(half.get('checks',[]))==68 else 'failed_unresolved_or_incomplete'),
        A95_response_native_weak=dict(status='pass' if len(weak_rows)==3 and all(x['status']=='pass' for x in weak_rows) else 'failed_unresolved_or_incomplete'),
        physical_budget_reconstruction=dict(status='pass' if len(budget_rows)==1 and budget_rows[0]['reconstruction_status']=='pass' else 'not_qualified'),
        physical_budget_accuracy=dict(status=budget_rows[0]['accuracy_status'] if len(budget_rows)==1 else 'not_qualified'))
    scope=dict(schema=SCHEMA,numerical_operator=OPERATOR,response_max_time=.016,out_of_scope=dict(A93_response=True),gates=gates)
    write_json(output/'required-scope-input.json',scope);decision=classify_scope(scope)
    write_json(output/'decision.json',decision);write_json(output/'costs-and-input-identities.json',dict(records=costs))
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
