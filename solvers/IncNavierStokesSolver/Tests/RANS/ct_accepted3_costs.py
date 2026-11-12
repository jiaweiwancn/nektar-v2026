#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Exclusive accepted evolution costs, separate from preparation and postmortem."""
import argparse,json,os
from pathlib import Path
from rans_qualification import THREADS,write_json
os.environ.update(THREADS)
from rans_request_manifest import RequestManifest,sha


def category(name):
    if name.startswith('geometry.'):return 'geometry_and_basis_rows'
    if name.startswith('initialization.'):return 'projection'
    if name.startswith(('domain.','constitutive.','pimple.check.')):return 'state_domain_and_constitutive'
    if name.endswith('.assembly') or name.startswith('assembly.'):return 'volume_interface_external_assembly'
    if name=='matching.constraints':return 'interface_and_constraints'
    if name=='dense.factorization' or name.endswith('.factorize'):return 'factorization'
    if name=='dense.triangular_solve' or name.endswith('.solve_validate'):return 'solves_and_linear_validation'
    if name.startswith(('pressure.','mass.')) or name=='flow.correction_validate':return 'mass_and_pressure_operations'
    if name.startswith(('audit.','native.acceptance_audit','native.physical_samples')):return 'physical_audits'
    if name.startswith(('io.','output.')):return 'snapshot_checkpoint_output_IO'
    return 'other_exclusive_orchestration_and_residual_work'


def scopes(path,after=None,before=None):
    names={};steps=set()
    for line in path.read_text().splitlines():
        q=json.loads(line)
        if after is not None and (q['step']<=0 or not after<q['time']<=before+1e-13):continue
        row=names.setdefault(q['name'],dict(exclusive_CPU=0.,exclusive_wall=0.,calls=0))
        row['exclusive_CPU']+=q['exclusive_process_cpu_seconds'];row['exclusive_wall']+=q['exclusive_wall_seconds'];row['calls']+=q.get('calls',1)
        if q['step']>0:steps.add(q['step'])
    groups={}
    for name,q in names.items():
        key=category(name);r=groups.setdefault(key,dict(exclusive_CPU=0.,exclusive_wall=0.))
        for k in r:r[k]+=q[k]
    total=sum(q['exclusive_CPU'] for q in names.values())
    return dict(by_scope=names,by_category=groups,total_instrumented_exclusive_CPU=total,accepted_steps=len(steps),
        CPU_per_accepted_step=total/len(steps) if steps else None,
        definition='exclusive intervals only; inclusive/nested parent totals never added')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    files=list(a.root.glob('**/results.json'))
    m=RequestManifest(a.output,dict(schema='CT3AcceptedCostSummaryV1',source_sha256=sha(Path(__file__)),
        receipts={str(f.relative_to(a.root)):sha(f) for f in files if f.name=='results.json'}))
    result=dict(status='pass',primary_segments=[],detailed_profiles=[],initializations=[]);reference_executions={}
    for f in files:
        q=json.loads(f.read_text())
        if q.get('schema')=='CT3InitialV1':
            for dt,row in q['members'].items():
                result['initializations'].append(dict(case=q['case'],dt=dt,slots=row['slots'],block_counts=row['block_counts'],
                    native_zero_step=row['executions'][0],observation_export=row['executions'][1],
                    timing=scopes(f.parent/dt/'timing.jsonl'),receipt_sha256=sha(f)))
        if q.get('schema')=='CT3SegmentV1':
            for dt,row in q['members'].items():
                reference_executions[row['native_report_sha256']]=row['execution']
                folder=Path(row['directory']);timing=scopes(folder/'timing.jsonl',q['start'],row['time'])
                result['primary_segments'].append(dict(case=q['case'],dt=dt,start=q['start'],target=q['target'],actual_time=row['time'],failure=row['failure'],
                    process=row['execution'],accepted_timing=timing,slots=len(json.loads((folder/'case.evolution.json').read_text())['state']),
                    rejected_attempt_timing=scopes(folder/'timing.jsonl',row['time'],q['target']) if row['failure'] else None,
                    timing_sha256=sha(folder/'timing.jsonl'),receipt_sha256=sha(f),
                    note='pinned primary runtime has broad scopes; geometry/law subcategories require separate equivalent detailed profile'))
        if q.get('bitwise_accepted_state_history_equivalence'):
            folder=f.parent/'run';row=dict(q,exclusive_grouping=scopes(folder/'timing.jsonl',.008,.02),
                full_program_grouping=scopes(folder/'timing.jsonl'),receipt_sha256=sha(f))
            row['reference_native_sha256']=json.loads((f.parent/'request-manifest.json').read_text())['request']['reference_sha256']
            row['whole_program_unattributed_CPU_including_MPI_helpers']=max(0.,q['execution']['process_cpu_seconds']-row['full_program_grouping']['total_instrumented_exclusive_CPU'])
            result['detailed_profiles'].append(row)
    for row in result['detailed_profiles']:
        original=reference_executions[row['reference_native_sha256']]
        row['ordinary_reference_execution']=original
        row['detailed_telemetry_CPU_overhead_ratio']=row['execution']['process_cpu_seconds']/original['process_cpu_seconds']
        row['interpretation']='exact exclusive measurements of the instrumented equivalent run; detailed clocks/maps add overhead. Do not interpret fine-grained percentages as ordinary-runtime fractions. Ordinary coarse scopes, costs and physical operation counts are separate.'
    result['profile_scope']='accepted .008-.02 segments where admitted; initialization, full-program output, independent audits and owned-launcher wait costs separate; no startup extrapolation'
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);print('cost summary',len(result['primary_segments']),len(result['detailed_profiles']));return 0
if __name__=='__main__':raise SystemExit(main())
