#!/usr/bin/env python3
"""Curate B1/B2 ledgers while retaining dense matrices/native fields on Linux."""
import argparse
from datetime import datetime, timezone
import hashlib
import itertools
import json
import os
from pathlib import Path
import shutil
import sys

os.sched_setaffinity(0,{48})
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root',required=True,type=Path)
parser.add_argument('--output',required=True,type=Path)
args=parser.parse_args()
args.output.mkdir(parents=True,exist_ok=False)
sys.path.insert(0,str(args.root))
import cls_epawf_m8_r1_nonlinear as analyzer


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def compressed(value):
    if isinstance(value,dict):
        result={}
        for key,item in value.items():
            if key.endswith('_branch_pattern') and isinstance(item,str):
                result[key]=dict(encoding='run_length',length=len(item),canonical_sha256=canonical(item),
                    runs=[[letter,len(list(run))] for letter,run in itertools.groupby(item)])
            else:
                result[key]=compressed(item)
        return result
    if isinstance(value,list):
        return [compressed(item) for item in value]
    return value


files=[]
def write(name,value,original,transform):
    path=args.output/name
    path.write_text(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n')
    files.append(dict(name=name,original_path=str(original),original_bytes=original.stat().st_size,
        original_sha256=sha(original),curated_bytes=path.stat().st_size,curated_sha256=sha(path),transform=transform))


def copy(name,original):
    shutil.copy2(original,args.output/name)
    files.append(dict(name=name,original_path=str(original),original_bytes=original.stat().st_size,
        original_sha256=sha(original),curated_bytes=(args.output/name).stat().st_size,
        curated_sha256=sha(args.output/name),transform='byte-for-byte copy'))


events=[]
for label,run,analysis in [('b1',args.root/'reproduction-001',args.root/'nonlinear-analysis-001'),
                           ('b2',args.root/'controlled/reproduction-001',args.root/'nonlinear-analysis-003')]:
    result_path=analysis/'validation.json';result=json.loads(result_path.read_text())
    trace_path=run/'case/case.nonlinear.json';trace=json.loads(trace_path.read_text())
    raw_path=run/'case/case.evolution.json';raw=json.loads(raw_path.read_text())
    source=Path(result['producing_runtime']['source_root'])
    audit=analyzer.Audit(raw,source)
    audit.gauge_variable_scale=audit.us/audit.H if label=='b2' else audit.us**2
    result['failed_stage_state']=trace['state']
    result['failed_stage_history']=trace['history']
    result['failed_stage_residual']=trace['residual']
    result['failed_stage_constraints']=trace['constraints']
    result['failed_stage_vectors_canonical_sha256']={key:canonical(trace[key]) for key in
        ('state','history','residual','constraints','cached_jacobian_state','cached_jacobian_history')}
    matrix_receipts={}
    for key in ('cached_jacobian','fresh_jacobian'):
        matrix=trace[key]
        matrix_receipts[key]=None if matrix is None else dict(rows=len(matrix),columns=len(matrix[0]),
            canonical_sha256=canonical(matrix),raw_trace_path=str(trace_path),raw_trace_sha256=sha(trace_path))
    result['dense_matrices_retained_on_linux']=matrix_receipts
    result['original_input_receipts']={name:dict(path=str(run/'case'/name),bytes=(run/'case'/name).stat().st_size,
        sha256=sha(run/'case'/name)) for name in ('case.xml','case.evolution.json','case.nonlinear.json')}
    # All finite-difference residual-block norms and complete branch strings
    # are retained (strings use lossless RLE). Full residual arrays stay in
    # the pinned raw trace; the failing state/history/residual remain above.
    if trace.get('directional_checks'):
        by_kind={item['kind']:item for item in trace['directional_checks']}
        for item in result['directional_checks']['families']:
            native=by_kind[item['kind']]
            item['direction']=native['direction']
            item['direction_canonical_sha256']=canonical(native['direction'])
            for level,original in zip(item['levels'],native['levels']):
                for side in ('plus','minus'):
                    value=original[side]
                    if value['evaluation_succeeded']:
                        constraints=audit.constraints(value['constraints'])
                        level['sides'][side]['full_residual_block_norms']=audit.norm_record(value['residual'],constraints)
                        level['sides'][side]['physical_branch_evidence']=value.get('physical')
                        level['sides'][side]['raw_residual_canonical_sha256']=canonical(value['residual'])
    result['curation']=dict(original_path=str(result_path),original_sha256=sha(result_path),
        raw_trace_path=str(trace_path),raw_trace_sha256=sha(trace_path),
        transform='Lossless RLE of all branch-pattern strings; adds full425-coefficient failed state/history/residual, constraints, '
            'direction vectors, every FD residual block norm, source/runtime/input pins and dense-matrix hashes. '
            'Dense matrices and complete native trajectories remain on Linux.')
    write('nonlinear-'+label+'-ledger.json',compressed(result),result_path,result['curation']['transform'])
    copy('nonlinear-'+label+'-runner.json',args.root/(analysis.name+'-runner.json'))
    copy('nonlinear-analyzer-'+label+'.py',analysis/'analyzer.py')
    for which in ('before','after'):
        path=run/('runtime-'+which+'.json');runtime=json.loads(path.read_text());omitted={}
        for key in ('source_inputs','documentation_tools','ldd'):
            value=runtime.pop(key,None)
            if value is not None:
                omitted[key]=dict(entries=len(value),canonical_sha256=canonical(value))
        runtime['curation']=dict(original_path=str(path),original_sha256=sha(path),original_bytes=path.stat().st_size,
            omitted=omitted,reason='Retain complete images/bindings/tools/environment/packages; omit only full source index, documentation tools and raw ldd with counts/hashes')
        write('nonlinear-'+label+'-runtime-'+which+'.json',runtime,path,'Explicit source/documentation/raw-ldd omissions only')
    before=json.loads((run/'runtime-before.json').read_text());after=json.loads((run/'runtime-after.json').read_text())
    for key in ('runtime_artifacts','runtime_bindings','runtime_environment','python_packages','tool_images'):
        if before[key]!=after[key]:
            raise ValueError('Runtime epoch changed: '+label+'/'+key)
    events.append(dict(epoch=label,source_commit=before['source_commit'],
        solver_sha256=before['runtime_artifacts']['build/IncNavierStokesSolver'],time=trace['time'],
        kind='first_failed_original_policy_stage',accepted_prefix_last_time=raw['time'],
        accepted_prefix_identical=result['accepted_prefix_equivalence'],
        original_trial_counts=result['original_policy']['classification_counts'],
        terminal_iteration=result['original_policy']['iterations'][-1],
        linear_or_assembly_failure=trace.get('assembly_failure'),fresh_retry=result.get('fresh_jacobian_same_stage_retry'),
        derivative_accuracy_failures=result.get('observed_Jacobian_accuracy_failures',[]),
        interpretation='Native original policy unchanged; diagnostic retries never advance history. Negative Newton iterates are not accepted trajectory states.'))

failed=args.root/'nonlinear-analysis-002/validation.json'
copy('nonlinear-analysis-002-failure.json',failed)
copy('nonlinear-analysis-002-runner.json',args.root/'nonlinear-analysis-002-runner.json')
copy('nonlinear-analyzer-failed-002.py',args.root/'nonlinear-analysis-002/analyzer.py')
write('nonlinear-events.json',compressed(dict(schema=1,events=events,
    inference_limits=['All terminal fresh fractions are classified from actual evaluations, not the broad terminal string.',
        'A fresh-J retry failing does not causally rule out stale-J effects.',
        'Residual-domain-admissible finite differences do not establish physical admissibility; base native-GLL k is negative.',
        'The pressure/gauge directional accuracy threshold failure is preserved without changing tolerances.'])),
    args.root/'controlled/reproduction-001/case/case.nonlinear.json','Chronological/event summary, lossless RLE branch evidence')
lines=['| Epoch | Producing source | Stage | Accepted / invalid / valid non-descent | Observation |',
       '|---|---|---:|---|---|']
for event in events:
    count=event['original_trial_counts']
    text='Original policy: '+str(count['accepted_descent'])+' / '+str(count['unevaluable_residual_or_norm'])+' / '+str(count['valid_non_descent'])
    lines.append('| '+event['epoch'].upper()+' | '+event['source_commit'][:7]+' | '+format(event['time'],'.6g')+' | '+text+' | last accepted .055; failed Newton iterate remains unaccepted |')
    if event['fresh_retry']:
        count=event['fresh_retry']['classification_counts']
        lines.append('| '+event['epoch'].upper()+' fresh retry | '+event['source_commit'][:7]+' | '+format(event['time'],'.6g')+' | '+str(count['accepted_descent'])+' / '+str(count['unevaluable_residual_or_norm'])+' / '+str(count['valid_non_descent'])+' | '+event['fresh_retry']['terminal']+' |')
(args.output/'nonlinear-events.md').write_text('\n'.join(lines)+'\n')
files.append(dict(name='nonlinear-events.md',original_path='derived from the two pinned ledgers by nonlinear-curate.py',
    original_sha256=sha(args.output/'nonlinear-events.md'),original_bytes=(args.output/'nonlinear-events.md').stat().st_size,
    curated_sha256=sha(args.output/'nonlinear-events.md'),curated_bytes=(args.output/'nonlinear-events.md').stat().st_size,
    transform='Derived table; exact code and both ledger inputs are retained'))
index=dict(schema=1,kind='M8-R1 B1/B2 independent nonlinear diagnosis',generated_utc=datetime.now(timezone.utc).isoformat(),
    curation_script_sha256=sha(__file__),files=files,
    dense_matrices_native_fields_accepted_histories_retained_on_linux=True,
    failed_state_history_and_residual_complete_425_entries_retained=True,
    source_epochs_separate={'b1':'426d48cd44f9f4af017244f6be4c1a6f62d36c01','b2':'360cbb7ac728b575cdf9af21e6d0963a1d88acb9'},
    preserved_analyzer_failure='Attempt002 detected the diagnostic gauge scale change .16 to .2; corrected epoch interpretation, no tolerance relaxation',
    policy='No solver/build/trajectory is launched by this curation; branch strings are losslessly compressed and matrix omissions have canonical hashes',
    canonical_encoding='json.dumps(sort_keys=True,separators=(comma,colon),allow_nan=False), UTF-8, ensure_ascii=True')
(args.output/'nonlinear-curation-index.json').write_text(json.dumps(index,indent=2,sort_keys=True)+'\n')
print(json.dumps(dict(output=str(args.output),files=len(files)+2,bytes=sum(p.stat().st_size for p in args.output.iterdir()))))
