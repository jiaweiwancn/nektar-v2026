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
def logical_size(root):
    total=0
    for folder,dirs,files in os.walk(root):
        for name in files:
            try:total+=(Path(folder)/name).lstat().st_size
            except FileNotFoundError:pass
    return total
def gates(root):
    q={k:read(root/v/'results.json') for k,v in dict(component='small-v7',static='static-v1',timestep='timestep-v1',pair='pair-v1',equivalence='equivalence-v1',segments='segments-v1').items()}
    status=lambda v:isinstance(v,dict) and v.get('status')=='pass'
    g1=all(status(q[k]) for k in ('component','static','timestep','equivalence'))
    g2=q['pair'] is not None and q['pair'].get('qualification',{}).get('status')=='pass'
    g3=q['segments'] is not None and q['segments'].get('cost_model_validation',{}).get('status')=='pass'
    g4=status(q['segments']) and len(q['segments'].get('cycles',[]))>=2
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
    segments=q['segments'];complete=segments is not None and segments.get('status')=='pass'
    phase_CPU=phase_wall=retention_CPU=retention_wall=milestone_CPU=milestone_wall=0.
    durable_per_segment=0;transient_per_step=0;measured_AS=0;milestone_bytes=0
    if complete:
        audits=[read(root/'segments-v1'/('audit-replay-'+str(i)+'.json')) for i in (2,4)]
        phase_CPU=envelope([u['CPU'] for a in audits for u in a['phase_usage']['per_state']])
        phase_wall=envelope([u['wall'] for a in audits for u in a['phase_usage']['per_state']])
        full_audits=[a for a in audits if a.get('scheduled_full_audit',True)]
        milestone_CPU=envelope([a['phase_usage']['scalar_budgets']['CPU']+a['phase_usage']['constitutive']['CPU'] for a in full_audits])
        milestone_wall=envelope([a['phase_usage']['scalar_budgets']['wall']+a['phase_usage']['constitutive']['wall'] for a in full_audits])
        held=read(root/'segments-v1/heldout-audits.json');milestone_CPU+=held['phase_usage']['actions']['CPU'];milestone_wall+=held['phase_usage']['actions']['wall']
        retention_CPU=envelope([max(0.,x['cycle_usage']['CPU']-x['execution']['whole_segment_CPU']-x['audit_usage']['CPU']) for x in segments['cycles']])
        retention_wall=envelope([max(0.,x['cycle_usage']['wall']-x['execution']['whole_segment_wall']-x['audit_usage']['wall']) for x in segments['cycles']])
        durable_per_segment=envelope([x['durable_owned_bytes'] for x in segments['cycles']])
        transient_per_step=envelope([x['retired_bytes']/4 for x in segments['cycles']])
        measured_AS=max(x['peak_address_space_bytes'] for x in segments['cycles'])
        milestone_bytes=sum(logical_size(root/'segments-v1/held-out-6-to8'/('actions-'+label)) for label in ('native','q48','q64'))
    # Selected short commands may become faster or slower later. These are
    # conditional extrapolations, never mathematical lower bounds on later PDEs.
    monitor_rate=segments.get('model',{}).get('external_monitor_CPU_per_wall_second',0.) if complete else 0.
    measured=dict(native=cost,observer_CPU_per_native_wall_second=observer_rate,controller_parse_startup_CPU=parse_startup,controller_parse_CPU_per_step=parse_step,external_monitor_CPU_per_wall_second=monitor_rate,independent_CPU_per_accepted_state=phase_CPU,independent_wall_per_accepted_state=phase_wall,retention_CPU_per_two_steps=retention_CPU,retention_wall_per_two_steps=retention_wall,full_milestone_CPU=milestone_CPU,full_milestone_wall=milestone_wall,short_two_step_durable_segment_bytes=durable_per_segment,transient_snapshot_bytes_per_actual_step=transient_per_step)
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
            components=dict(native_CPU=native_cpu,replay_CPU=replay_cpu,controller_CPU=controller,independent_state_CPU=N*phase_CPU,milestone_CPU=A*milestone_CPU,retention_CPU=math.ceil(N/2)*retention_CPU,collection_transfer_publication_CPU=240.)
            wall=native_wall+replay_wall+controller+N*phase_wall+A*milestone_wall+math.ceil(N/2)*retention_wall+600.
            components['external_monitor_CPU']=wall*monitor_rate;cpu=sum(components.values())
            # A short segment's entire retained inventory per two steps is a
            # conservative finite-horizon growth bound, not a claim of constant
            # resident storage. It includes replay reports/anchors/receipts.
            durable=immutable+math.ceil(N/2)*durable_per_segment+A*milestone_bytes
            live=durable+2*32*transient_per_step+256*1024**2
            admitted=complete and qualification_passed and 1.25*cpu<=24000 and 1.25*wall<=28800 and 1.25*live<=24*GIB and measured_AS<=16*GIB
            rows.append(dict(horizon=name,requested_T=str(T),actual_endpoints=[str(Decimal(n)*Decimal(str(dt))) for n,dt in zip(counts,(.00025,.000125))],from_state=start,start_is_qualified=(start=='fresh_zero' or q['pair'] is not None and q['pair'].get('qualification',{}).get('status')=='pass'),total_member_steps=counts,additional_member_steps=additional,segment_length=32,native_hard_segment_cap=64,segments=segcounts,full_audit_counts=audits,CPU_components=components,CPU_before_headroom=cpu,CPU_with_headroom=1.25*cpu,active_wall_with_headroom=1.25*wall,finite_retained_output_with_headroom=1.25*durable,peak_new_live_output_with_headroom=1.25*live,measured_short_AS_bytes=measured_AS,status='ADMITTED' if admitted else 'UNADMITTED',native_only_CPU_with_headroom=1.25*native_cpu,required_total_cost_reduction=max(1.,1.25*cpu/24000),scope_complete=complete))
    return dict(status='complete_conditional_forecast' if complete else 'incomplete_native_replay_extrapolation',headroom=1.25,headroom_applied='once to complete future totals; not sunk research',measurements=measured,horizons=rows,validity='short-history bounds are conditional extrapolations; all proposed campaigns require live spent-plus-remaining guards and segment fill/iteration reforecasting. No long-time physical qualification.',necessary_zero_overhead_CPU_per_step=24000/(1.25*4884),strict_necessary_bound_scope='arithmetic budget bound only; observed short native timing is not a proven lower bound on all future states')

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
    result=dict(schema='CT15CompactPublicationV1',status='UNRESOLVED',gates=g,physical_T0_407='UNATTEMPTED',full_transit_scientific_and_operational_readiness='UNQUALIFIED',external_root=str(a.root),historical_ledger=read(a.root/'baseline-v1/results.json'),engineering_ledger=dict(closed_CPU=CPU,publication_transport_CPU_reserve=240,total_CPU_upper_bound_before_current_collector=CPU+240,elapsed_wall_including_connection_outage=time.time()-START,logical_retained_bytes=logical_size(a.root),peak_sampled_live_bytes=peak,sampled_cumulative_process_wchar_lower_bound=sampled_writes,write_accounting_scope='sum of per-family /proc wchar highwaters; includes pipe writes and may miss short-lived processes; not exact physical disk bytes',exception_reserves=exception,unclosed_families=unclosed,failed_worker_audits=worker_failures,resource_receipts=resources),future_ledger=future,artifact_inventory=inventory)
    # The historical body contains stage paths/hashes but no numerical arrays.
    write(a.output/'results.json',result);write(a.output/'forecast.json',future);write(a.output/'artifact-inventory.json',inventory)
    write(a.output/'seal.json',dict(files={f.name:digest(f) for f in a.output.iterdir() if f.is_file()}))
    print(json.dumps(dict(status=result['status'],gates=g,closed_CPU=CPU,forecast_status=future['status']),indent=2))
if __name__=='__main__':main()
