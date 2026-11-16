#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Publish the unresolved no-G4 path from sealed Linux metadata on the Mac.

This is receipt compaction and forecast arithmetic, not numerical validation.
If actual retention experiments exist, use the complete workflow collector.
"""
import argparse, json
from pathlib import Path
from ct_accepted15_collect import digest,read,write


def require(test,message):
    if not test:raise ValueError(message)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('snapshot','forecast','pair','output'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--equivalence',type=Path);a=p.parse_args()
    seal=read(a.snapshot.parent/'seal.json');require(seal['files']['snapshot.json']==digest(a.snapshot),'metadata snapshot changed')
    s=read(a.snapshot);f=read(a.forecast);pair=read(a.pair);art=s['artifacts']
    require(not s['source_working_tree'],'dirty private controller checkout')
    require(not s['unclosed_resources'],'unclosed Linux family; publication cannot close it')
    require('segments-v1/results.json' not in art,'actual G4 data require complete workflow collector')
    require(art['pair-v1/results.json']['sha256']==pair['external_results_sha256']==pair['seal']['files']['results.json'],'foreign short pair')
    require(pair['request']['sha256']==pair['seal']['request_sha256'],'pair request/seal mismatch')
    require(pair['pair']['qualification']['status']=='pass','short pair did not pass')
    for relative,suffix in [('training-v4/results.json','/training-v4.json'),('pair-v1/results.json','/pair-v1/results.json'),('static-v1/results.json','/static-v1/results.json'),('timestep-v1/results.json','/timestep-v1/results.json')]:
        values=[v for k,v in f['source_receipt_SHA256'].items() if k.endswith(suffix)]
        require(len(values)==1 and values[0]==art[relative]['sha256'],'forecast input differs: '+relative)
    require(f['whole_workflow_component_coverage'] is False and all(r['CPU_with_headroom'] is None and r['status']=='UNADMITTED' for r in f['horizons']),'partial forecast relabelled as complete')
    equivalence=None
    if 'equivalence-v1/results.json' in art:
        require(a.equivalence is not None and digest(a.equivalence)==art['equivalence-v1/results.json']['sha256'],'missing actual equivalence receipt')
        equivalence=read(a.equivalence)
    resources=s['closed_resources'];exception=s['exception_reserves'];CPU=sum(r['process_cpu_seconds'] for r in resources.values())+exception['CPU_seconds']
    live=s['directory_logical_bytes']['.'];peak=live;writes=0
    for folder in resources:
        watch=s['watches'].get(folder+'/package-budget-watch.json') or s['watches'].get(folder+'/live-budget.json') or {}
        sample=watch.get('last_sample',{});peak=max(peak,sample.get('peak_sampled_live_bytes',0));writes+=sample.get('sampled_cumulative_process_wchar',0)
    failed={k:v['resource'] for k,v in resources.items() if v['resource']['aggregate_status']!='pass' or v['resource']['project_status']!='pass'}
    reserve=dict(CPU=240,wall=240,bytes=256*1024**2)
    accounting=dict(schema='CT15ClosedPublicationAccountingV1',status='closed_receipts_with_explicit_publication_reserves',
        snapshot_sha256=digest(a.snapshot),source_at_snapshot=s['source_revision'],external_root=s['external_root'],
        closed_family_count=len(resources),closed_family_CPU_before_exception=sum(r['process_cpu_seconds'] for r in resources.values()),
        exception_reserves=exception,closed_CPU_including_exception=CPU,publication_reserve=reserve,
        CPU_upper_with_reserve=CPU+reserve['CPU'],elapsed_wall_at_snapshot=s['whole_elapsed_wall'],wall_upper_with_reserve=s['whole_elapsed_wall']+reserve['wall'],
        logical_live_bytes_at_snapshot=live,peak_observed_live_bytes=peak,new_live_upper_with_publication_reserve=peak+reserve['bytes'],
        sampled_cumulative_process_wchar_lower_bound=writes,write_accounting='observed per-family process wchar highwaters, including pipes; may miss short-lived processes and the unadmitted queue helpers; not exact physical disk writes',
        failed_worker_audits=failed,unclosed_resources=[],
        reserve_scope='fixed existing source/inspection/transport/publication allowance, including metadata exporters and Mac compaction; no second reserve added for those helpers',
        scope='all closed CT15 Linux families plus the documented exception and administrative reserves; historical and prospective campaign costs excluded')
    accounting['finite_CPU_wall_disk_checks']=dict(CPU=accounting['CPU_upper_with_reserve']<=24000,wall=accounting['wall_upper_with_reserve']<=28800,live_bytes=accounting['new_live_upper_with_publication_reserve']<=16*1024**3)
    gates=dict(G1='PASS' if equivalence and equivalence.get('status')=='pass' else 'UNATTEMPTED',G2='PASS',G3='UNATTEMPTED',G4='UNATTEMPTED',G5='FAIL',G6='PREPARED_FOR_FINAL_PUSH')
    decision=dict(schema='CT15UnresolvedClosureV1',status='UNRESOLVED',resource_admission_goal='NOT_ACHIEVED',gates=gates,
        gate_scope=dict(G1='small/static/changed-alpha controls passed; full trajectory comparison is required',G3='startup/step and partial workflow costs measured; unchanged heldout6-to8 not executed',G4='selected multi-cycle retirement/recovery remains unattempted',G5='no complete validated forecast; native-only extrapolation and measured probe AS with headroom exceed ceilings',G6='closed hash-bound metadata ready; final source/publication push and remote containment still required'),
        physical_T0_407='UNATTEMPTED',full_transit='UNQUALIFIED',snapshot_sha256=digest(a.snapshot),forecast_sha256=digest(a.forecast),pair_compact_sha256=digest(a.pair),
        historical_ledger=s['historical_ledger'],engineering_ledger=accounting,future_ledger=f,
        actual_equivalence=equivalence,source_and_runtime=pair['request']['request'],
        reproduction_scope='authoritative Linux numerical receipts; read-only Linux metadata export; hash checks and compaction on Mac. No new solve, physical classifier or operator validation is run by this publisher.')
    a.output.mkdir(exist_ok=False)
    write(a.output/'results.json',decision);write(a.output/'closed-publication-accounting.json',accounting)
    write(a.output/'seal.json',dict(files={n:digest(a.output/n) for n in ('results.json','closed-publication-accounting.json')}))
    print(json.dumps(dict(status='UNRESOLVED',gates=gates,CPU_upper=accounting['CPU_upper_with_reserve'],wall_upper=accounting['wall_upper_with_reserve'],worker_failures=list(failed))))


if __name__=='__main__':raise SystemExit(main())
