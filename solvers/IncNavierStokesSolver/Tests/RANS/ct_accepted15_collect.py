#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact CT15 evidence and conditional forecasts; missing gates never pass."""
import argparse,hashlib,json,math,os,time
from decimal import Decimal,ROUND_CEILING
from pathlib import Path

START=1794805033
GIB=1024**3
def read(path):return json.loads(path.read_text()) if path.is_file() else None
def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()
def write(path,value):
    with path.open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')
def count(T,dt):return int((Decimal(str(T))/Decimal(str(dt))).to_integral_value(rounding=ROUND_CEILING))
def envelope(v):return max(v)+max(v)-min(v) if v else 0.
def compact(value):
    """Keep decisions and witnesses; bind large external payloads by hash."""
    if isinstance(value,str) and len(value)>4096:
        return dict(external_string_characters=len(value),UTF8_sha256=hashlib.sha256(value.encode()).hexdigest())
    if isinstance(value,list) and len(value)>32:
        payload=json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
        return dict(external_array_entries=len(value),canonical_JSON_sha256=hashlib.sha256(payload).hexdigest())
    if isinstance(value,dict):return {k:compact(v) for k,v in value.items()}
    if isinstance(value,list):return [compact(v) for v in value]
    return value
def max_address_space(value):
    if isinstance(value,dict):return max([v*1024 for k,v in value.items() if k in ('peak_address_space_KiB','peak_virtual_KiB')]+[max_address_space(v) for v in value.values() if isinstance(v,(dict,list))]+[0])
    if isinstance(value,list):return max([max_address_space(v) for v in value]+[0])
    return 0
def logical_size(root):
    total=0
    for folder,dirs,files in os.walk(root):
        for name in files:
            try:total+=(Path(folder)/name).lstat().st_size
            except FileNotFoundError:pass
    return total

def historical_ledger(root,ct14):
    """Prior closed receipts are sunk research, never future campaign work."""
    result=dict(scope='historical CT11-CT14 only; no amount is added to CT15 engineering or prospective production',campaigns={})
    for n in (11,12,13,14):
        directory=ct14.parent/('ct-accepted-'+str(n)+'-20261115')
        receipts={str(f.relative_to(directory)):dict(sha256=digest(f),CPU=read(f)['process_cpu_seconds'],wall=read(f)['wall_seconds']) for f in directory.glob('resources-*/execution.json')}
        if not receipts:raise ValueError('missing authoritative historical receipts: '+str(directory))
        name='resource-accounting-v2.json' if n==11 else 'closed-publication-accounting.json'
        published=root/'source/docs/plans'/('cao-tamura-accepted-'+str(n)+'-evidence')/name
        q=read(published)
        result['campaigns'][str(n)]=dict(external_root=str(directory),closed_family_count=len(receipts),closed_Linux_family_CPU=sum(v['CPU'] for v in receipts.values()),summed_family_elapsed_wall=sum(v['wall'] for v in receipts.values()),receipt_hashes=receipts,published_accounting_sha256=digest(published),published_measured_CPU=q.get('process_CPU_seconds',q.get('closed',{}).get('closed_CPU_seconds')),published_upper_bounds=q.get('upper_bounds'),publication_scope=q.get('scope','prospective formula file; authoritative completed-family sums reported separately'))
    return result
def gates(root):
    q={k:read(root/v/'results.json') for k,v in dict(component='small-v7',static='static-v1',timestep='timestep-v1',pair='pair-v1',equivalence='equivalence-v1',segments='segments-v1').items()}
    status=lambda v:isinstance(v,dict) and v.get('status')=='pass'
    g1=all(status(q[k]) for k in ('component','static','timestep','equivalence'))
    g2=q['pair'] is not None and q['pair'].get('qualification',{}).get('status')=='pass'
    g3=q['segments'] is not None and q['segments'].get('cost_model_validation',{}).get('status')=='pass'
    retention_path=root/'segments-v1/retention-controls.json';retention=read(retention_path);retention_seal=read(root/'segments-v1/retention-controls-seal.json')
    retention_request=read(root/'segments-v1/request-manifest.json')
    g4=status(retention) and len(retention.get('cycles',[]))>=2 and retention_seal is not None and retention_request is not None and retention_seal['request_sha256']==retention_request['sha256'] and retention_seal['files']['retention-controls.json']==digest(retention_path)
    q['retention']=retention
    return q,{k:('PASS' if v else 'UNRESOLVED') for k,v in dict(G1=g1,G2=g2,G3=g3,G4=g4).items()}

def forecast(root,q):
    training=read(root/'training-v4/results.json')
    if training is None:return dict(status='UNATTEMPTED',reason='no complete selected timing-off training')
    cases=training['cases'];cost={}
    for key,field in (('CPU','process_cpu_seconds'),('wall','wall_seconds')):
        v=[cases[str(i)]['execution'][field] for i in (0,1,2)]
        cost[key]=dict(startup=max(0.,v[0]),step=max(0.,v[1]-v[0],(v[2]-v[0])/2,v[2]-v[1]))
    observer_rate=max(x['execution'].get('observer_self_CPU_seconds',0.)/x['execution']['wall_seconds'] for x in cases.values())
    parse=[max(0.,cases[str(i)]['execution'].get('controller_evolve_self_CPU',0.)-cases[str(i)]['execution'].get('observer_self_CPU_seconds',0.)) for i in (0,1,2)]
    parse_startup=parse[0];parse_step=max(0.,parse[1]-parse[0],(parse[2]-parse[0])/2,parse[2]-parse[1])
    segments=q['segments'];complete=segments is not None and (q.get('retention') or {}).get('status')=='pass' and q['pair'] is not None and q['pair'].get('qualification',{}).get('status')=='pass'
    validated=complete and segments.get('cost_model_validation',{}).get('status')=='pass'
    phase_CPU=phase_wall=retention_CPU=retention_wall=milestone_CPU=milestone_wall=0.
    durable_per_segment=0;transient_per_step=0;milestone_bytes=0
    measured_AS=max(max_address_space(q.get('pair')),max_address_space(q.get('static')),max_address_space(q.get('timestep')))
    if complete:
        audits=[read(root/'segments-v1'/('audit-replay-'+str(i)+'.json')) for i in (2,4)]
        phase_CPU=envelope([u['CPU'] for a in audits for u in a['phase_usage']['per_state']])
        phase_wall=envelope([u['wall'] for a in audits for u in a['phase_usage']['per_state']])
        full_audits=[a for a in audits if a.get('scheduled_full_audit',True)]
        milestone_CPU=envelope([a['phase_usage']['scalar_budgets']['CPU']+a['phase_usage']['constitutive']['CPU'] for a in full_audits])
        milestone_wall=envelope([a['phase_usage']['scalar_budgets']['wall']+a['phase_usage']['constitutive']['wall'] for a in full_audits])
        held=read(root/'segments-v1/heldout-audits.json')
        if held:
            milestone_CPU+=held['phase_usage']['actions']['CPU'];milestone_wall+=held['phase_usage']['actions']['wall']
        else:
            phases=q['pair']['workflow_phases'];names=('complete_endpoint_scalar_controls','endpoint_constitutive_probe')
            # An explicitly conservative per-member upper bound uses the whole
            # measured TWO-member endpoint wrappers and whole changed-alpha
            # family. Their nested native costs are not added again.
            milestone_CPU=sum(phases[k]['self_CPU']+phases[k]['children_CPU'] for k in names)+read(root/'resources-timestep-v1/execution.json')['process_cpu_seconds']
            milestone_wall=sum(phases[k]['wall'] for k in names)+read(root/'resources-timestep-v1/execution.json')['wall_seconds']
        retention_CPU=envelope([max(0.,x['cycle_usage']['CPU']-x['execution']['whole_segment_CPU']-x['audit_usage']['CPU']) for x in segments['cycles']])
        retention_wall=envelope([max(0.,x['cycle_usage']['wall']-x['execution']['whole_segment_wall']-x['audit_usage']['wall']) for x in segments['cycles']])
        durable_per_segment=envelope([x['durable_owned_bytes'] for x in segments['cycles']])
        transient_per_step=envelope([x['retired_bytes']/4 for x in segments['cycles']])
        measured_AS=max([x['peak_address_space_bytes'] for x in segments['cycles']]+[segments.get('held_out',{}).get('peak_address_space_bytes',0),measured_AS])
        milestone_bytes=sum(logical_size(root/'segments-v1/held-out-6-to8'/('actions-'+label)) for label in ('native','q48','q64')) if held else sum(logical_size(root/'timestep-v1'/k) for k in ('4000-native','4000-q48','4000-q64','8000-native','8000-q48','8000-q64'))
    # Selected short commands may become faster or slower later. These are
    # conditional extrapolations, never mathematical lower bounds on later PDEs.
    monitor_rate=segments.get('model',{}).get('external_monitor_CPU_per_wall_second',0.) if complete else 0.
    measured=dict(native=cost,observer_CPU_per_native_wall_second=observer_rate,controller_parse_startup_CPU=parse_startup,controller_parse_CPU_per_step=parse_step,external_monitor_CPU_per_wall_second=monitor_rate,independent_CPU_per_accepted_state=phase_CPU,independent_wall_per_accepted_state=phase_wall,retention_CPU_per_two_steps=retention_CPU,retention_wall_per_two_steps=retention_wall,full_milestone_CPU=milestone_CPU,full_milestone_wall=milestone_wall,short_two_step_durable_segment_bytes=durable_per_segment,transient_snapshot_bytes_per_actual_step=transient_per_step)
    if not complete:
        for key in ('external_monitor_CPU_per_wall_second','independent_CPU_per_accepted_state','independent_wall_per_accepted_state','retention_CPU_per_two_steps','retention_wall_per_two_steps','full_milestone_CPU','full_milestone_wall','short_two_step_durable_segment_bytes','transient_snapshot_bytes_per_actual_step'):
            measured[key]=None
    measured['milestone_bound_scope']=('actual heldout endpoint phases' if held else 'whole measured two-member endpoint wrappers plus whole changed-alpha family per future milestone; conservative scope bound, not exclusive timing or a validated heldout') if complete else 'unmeasured'
    paired=q['pair'].get('workflow_phases',{}) if q['pair'] else {};common_count=len(q['pair'].get('common_times',[])) if q['pair'] else 0
    paired_CPU=sum(paired.get(k,{}).get('self_CPU',0.)+paired.get(k,{}).get('children_CPU',0.) for k in ('paired_common_state_and_endpoint_comparison','paired_regional_time_series'))/max(1,common_count)
    paired_wall=sum(paired.get(k,{}).get('wall',0.) for k in ('paired_common_state_and_endpoint_comparison','paired_regional_time_series'))/max(1,common_count)
    measured['paired_CPU_per_common_state']=paired_CPU;measured['paired_wall_per_common_state']=paired_wall
    horizons=[('entry',Decimal('.407'))]+[(str(f)+'transit',Decimal(str(f))*Decimal('40.6216378576')) for f in (.1,1,3)]
    rows=[];qualification_passed=all(g=='PASS' for g in gates(root)[1].values())
    immutable=logical_size(root/'build-v10/runtime')
    packet=root/'static-v1/target/frozen.json'
    if packet.exists():immutable+=packet.stat().st_size
    for name,T in horizons:
        counts=[count(T,dt) for dt in (.00025,.000125)]
        for start in ('fresh_zero','own_qualified_T0.002'):
            additional=counts if start=='fresh_zero' else [counts[0]-8,counts[1]-16]
            segcounts=[math.ceil(v/32) for v in additional];N=sum(additional);S=sum(segcounts)
            # Required full audits: initial and .002 for fresh runs, each member
            # step512, and rounded final endpoint. Deduplicate coincident steps.
            audits=[]
            for i,n in enumerate(counts):
                full={0,(8,16)[i],n}|set(range(512,n+1,512))
                if start!='fresh_zero':full={j for j in full if j>(8,16)[i]}
                audits.append(len(full))
            A=sum(audits);initial=2 if start=='fresh_zero' else 0
            native_cpu=N*cost['CPU']['step']+(S+initial)*cost['CPU']['startup']
            native_wall=N*cost['wall']['step']+(S+initial)*cost['wall']['startup']
            # Retired payloads are verified by an actual second native replay.
            replay_cpu=N*cost['CPU']['step']+S*cost['CPU']['startup']
            replay_wall=N*cost['wall']['step']+S*cost['wall']['startup']
            controller=(native_wall+replay_wall)*observer_rate+(2*S+initial)*parse_startup+2*N*parse_step
            components=dict(native_CPU=native_cpu,replay_CPU=replay_cpu,controller_CPU=controller,independent_state_CPU=N*phase_CPU,paired_comparison_CPU=(additional[0]+1)*paired_CPU,milestone_CPU=A*milestone_CPU,retention_CPU=math.ceil(N/2)*retention_CPU,collection_transfer_publication_CPU=240.)
            wall=native_wall+replay_wall+controller+N*phase_wall+(additional[0]+1)*paired_wall+A*milestone_wall+math.ceil(N/2)*retention_wall+600.
            components['external_monitor_CPU']=wall*monitor_rate;cpu=sum(components.values())
            # A short segment's entire retained inventory per two steps is a
            # conservative finite-horizon growth bound, not a claim of constant
            # resident storage. It includes replay reports/anchors/receipts.
            durable=immutable+math.ceil(N/2)*durable_per_segment+A*milestone_bytes
            live=durable+2*32*transient_per_step+256*1024**2
            admitted=complete and qualification_passed and 1.25*cpu<=24000 and 1.25*wall<=28800 and 1.25*live<=24*GIB and 1.25*measured_AS<=16*GIB
            rows.append(dict(horizon=name,requested_T=str(T),actual_endpoints=[str(Decimal(n)*Decimal(str(dt))) for n,dt in zip(counts,(.00025,.000125))],from_state=start,start_is_qualified=(start=='fresh_zero' or q['pair'] is not None and q['pair'].get('qualification',{}).get('status')=='pass'),total_member_steps=counts,additional_member_steps=additional,segment_length=32,native_hard_segment_cap=64,segments=segcounts,full_audit_counts=audits,CPU_components=components,CPU_before_headroom=cpu,CPU_with_headroom=1.25*cpu,active_wall_with_headroom=1.25*wall,finite_retained_output_with_headroom=1.25*durable,peak_new_live_output_with_headroom=1.25*live,measured_short_AS_bytes=measured_AS,status='ADMITTED' if admitted else 'UNADMITTED',native_only_CPU_with_headroom=1.25*native_cpu,required_total_cost_reduction=max(1.,1.25*cpu/24000),scope_complete=complete))
            rows[-1]['native_only_active_wall_with_headroom']=1.25*native_wall
    for row in rows:
        row['peak_AS_with_headroom']=1.25*row['measured_short_AS_bytes']
        if complete:
            N=sum(row['additional_member_steps']);S=sum(row['segments']);initial=2 if row['from_state']=='fresh_zero' else 0;c=row['CPU_components']
            fixed=sum(c[k] for k in ('independent_state_CPU','paired_comparison_CPU','milestone_CPU','retention_CPU','collection_transfer_publication_CPU'))+(2*S+initial)*parse_startup+2*N*parse_step+(2*S+initial)*cost['CPU']['startup']
            row['other_fitted_checks_and_startups_CPU']=fixed
            row['optimistic_remaining_native_CPU_per_step_with_other_fitted_checks_fixed']=(24000/1.25-fixed)/(2*N)
            row['necessary_throughput_scope']='holding measured independent/retention/parse/startup costs fixed, and optimistically giving remaining native observers/monitors zero cost; replay doubles actual native step work. Negative allowance requires improvements beyond native evolution.'
        else:
            # A missing recovery experiment cannot turn unknown archive/audit
            # costs into zero-cost whole-workflow forecasts. Preserve the
            # measured partial extrapolation under an explicitly named scope.
            row['partial_observed_scope_CPU_with_headroom']=row['CPU_with_headroom']
            row['partial_observed_scope_wall_with_headroom']=row['active_wall_with_headroom']
            row['partial_observed_scope']='native evolution/replay, measured native observer/parse and paired comparisons, plus publication reserve only; missing independent audits, retention, full milestones and external monitors'
            row['required_total_cost_reduction']=None
            for key in ('independent_state_CPU','milestone_CPU','retention_CPU','external_monitor_CPU'):
                row['CPU_components'][key]=None
            for key in ('CPU_before_headroom','CPU_with_headroom','active_wall_with_headroom','finite_retained_output_with_headroom','peak_new_live_output_with_headroom'):
                row[key]=None
    return dict(status=('complete_conditional_forecast' if validated else 'complete_unvalidated_conditional_forecast') if complete else 'incomplete_native_replay_extrapolation',whole_workflow_component_coverage=complete,heldout_model_validated=bool(validated),headroom=1.25,headroom_applied='once to complete future CPU/wall/output/AS totals; not sunk research',measurements=measured,horizons=rows,unmeasured_additional_costs=[] if complete else ['complete heldout external-monitor validation','per-state independent audits','all full milestones/actions','retention and finite durable growth'],validity='short-history bounds are conditional extrapolations; all proposed campaigns require live spent-plus-remaining guards and segment fill/iteration reforecasting. No long-time physical qualification.',necessary_zero_overhead_CPU_per_step=24000/(1.25*4884),strict_necessary_bound_scope='arithmetic budget bound only; observed short native timing is not a proven lower bound on all future states')

def main():
    p=argparse.ArgumentParser()
    for key in ('root','ct14','output'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--active-control',type=Path);a=p.parse_args();a.output.mkdir(exist_ok=False)
    q,g=gates(a.root);future=forecast(a.root,q)
    g['G5']='PASS' if future.get('status')=='complete_conditional_forecast' and all(r['status']=='ADMITTED' for r in future['horizons'] if r['horizon']=='entry') else 'FAIL' if future.get('horizons') else 'UNRESOLVED'
    resources={};unclosed=[];worker_failures=[];sampled_writes=0;peak=0
    for folder in sorted(a.root.glob('resources-*')):
        f=folder/'execution.json'
        if f.is_file():
            r=read(f);resources[folder.name]=dict(receipt_sha256=digest(f),**r)
            if r.get('resource',{}).get('aggregate_status')!='pass' or r.get('resource',{}).get('project_status')!='pass':worker_failures.append(folder.name)
        elif folder!=a.active_control:unclosed.append(folder.name)
        watch=read(folder/'package-budget-watch.json') or read(folder/'live-budget.json')
        if watch:
            s=watch.get('last_sample',{});sampled_writes+=s.get('sampled_cumulative_process_wchar',0);peak=max(peak,s.get('peak_sampled_live_bytes',0))
    exception=read(a.root/'exception-reserves.json') or dict(CPU_seconds=0)
    CPU=sum(r['process_cpu_seconds'] for r in resources.values())+exception['CPU_seconds']
    inventory={}
    for folder in sorted(a.root.iterdir()):
        if folder.is_dir() and folder.name not in ('source','build'):
            for name in ('request-manifest.json','completed-artifacts.json','results.json','progress.json','execution.json'):
                path=folder/name
                if path.is_file():inventory[str(path.relative_to(a.root))]=dict(sha256=digest(path),bytes=path.stat().st_size)
    g['G6']='PENDING_FINAL_PUSH_AND_CLOSED_ACCOUNTING'
    result=dict(schema='CT15CompactPublicationV1',status='UNRESOLVED',gates=g,physical_T0_407='UNATTEMPTED',full_transit_scientific_and_operational_readiness='UNQUALIFIED',external_root=str(a.root),historical_ledger=historical_ledger(a.root,a.ct14),CT14_baseline_verification_sha256=digest(a.root/'baseline-v1/results.json'),engineering_ledger=dict(closed_CPU=CPU,publication_transport_CPU_reserve=240,total_CPU_upper_bound_before_current_collector=CPU+240,elapsed_wall_including_connection_outage=time.time()-START,logical_retained_bytes=logical_size(a.root),peak_sampled_live_bytes=peak,sampled_cumulative_process_wchar_lower_bound=sampled_writes,write_accounting_scope='sum of per-family /proc wchar highwaters; includes pipe writes and may miss short-lived processes; not exact physical disk bytes',exception_reserves=exception,unclosed_families=unclosed,failed_worker_audits=worker_failures,resource_receipts=resources),future_ledger=future,artifact_inventory=inventory)
    # The historical body contains stage paths/hashes but no numerical arrays.
    write(a.output/'results.json',result);write(a.output/'forecast.json',future);write(a.output/'artifact-inventory.json',inventory)
    write(a.output/'scientific-evidence.json',dict(schema='CT15CompactScientificEvidenceV1',scope='large arrays and signatures are external hash-bound payloads; this publication does not replace raw Linux evidence',gates=g,evidence={k:compact(v) for k,v in q.items()}))
    write(a.output/'seal.json',dict(files={f.name:digest(f) for f in a.output.iterdir() if f.is_file()}))
    print(json.dumps(dict(status=result['status'],gates=g,closed_CPU=CPU,forecast_status=future['status']),indent=2))
if __name__=='__main__':main()
