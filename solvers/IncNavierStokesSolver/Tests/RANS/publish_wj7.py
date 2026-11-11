#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact WJ-7 tables and semantic source/stage lineage; no governing solves."""
import argparse
import csv
import json
from pathlib import Path
import time
import subprocess
from wj_immutable import (Ledger,read,write_once,file_identity,file_hash,artifacts,verify_artifacts,digest)
from run_cls_epawf_picard import require,state_hash
from wj6_resources import volume,DEFAULT
from wj_qualification import component,aggregate
from wj7_qualification import expected_tuples,tuple_gate
from publish_wj6 import table,execution_stage_binding


def copy_exact(source,target):
    target.parent.mkdir(parents=True,exist_ok=True)
    with target.open('xb') as stream:stream.write(source.read_bytes())


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--parent',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    root=a.root.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    rows={k:[] for k in ('fields','interfaces','residuals','sensitivity','runs','counts','arithmetic','resources')}
    endpoints={};stages={};components={};heavy={};summaries={};coverage={}
    for policy in ('EP','Budget'):
        folder=root/'study-v2'/policy;s=read(folder/'summary.json')
        frozen=root/'snapshots'/(policy+'-'+file_hash(folder/'summary.json')+'.json');frozen.parent.mkdir(parents=True,exist_ok=True)
        raw=(folder/'summary.json').read_bytes()
        if frozen.exists():require(frozen.read_bytes()==raw,'frozen summary bytes differ')
        else:
            with frozen.open('xb') as stream:stream.write(raw)
        summaries[policy]=file_identity(frozen)
        quals=sorted(folder.glob('qualification-*.json'),key=lambda p:int(p.stem.split('-')[-1]))
        last=32;gates={}
        for path in quals:
            q=read(path);last=int(q['end']);gates[str(last)]=dict(gate=q['gate'],packet=file_identity(path))
            require(tuple_gate(q['sensitivity'],last)['status']=='PASS','publication tuple coverage missing/duplicated')
            for v in q['sensitivity']:
                order,group,step,scope,field,norm=v['key']
                rows['sensitivity'].append(dict(policy=policy,qualification_end=last,order=order,group=group,step=step,scope=scope,field=field,norm=norm,
                    status=v['status'],value=v['value'],uncertainty=v['uncertainty'],limit=v['limit'],margin=v['limit']-v['value']-v['uncertainty']))
        components[policy]=dict(imports=s['components'].get('imports',component('NOT_RUN')),prefixes=gates,
            final_target=component('PASS' if last==128 and gates['128']['gate']['status']=='PASS' else 'NOT_RUN','fixed T3 target'),
            completed_capacity_transaction_reused=s['components']['transactions'])
        actual={tuple(v['key']) for v in read(quals[-1])['sensitivity']}
        missing=expected_tuples(128)-actual
        coverage[policy]=dict(status='NOT_RUN' if missing else 'PASS',required_final=1008,observed=len(actual),
            missing=[list(v) for v in sorted(missing)],new_arithmetic_analyses_required=9,
            actual_new_arithmetic_analyses=sum(bool(c.get('arithmetic')) for t in s['trajectories'].values() for c in t['trajectory'] if c['steps']>(64 if t['dt']==1.6e-5 else 32)))
        endpoints[policy]={}
        for tag,t in s['trajectories'].items():
            for c8 in [t['initial']]+t['trajectory']:
                step=int(c8['steps']);q10=s['order10'][tag+'/'+str(step)]
                for order,c in ((8,c8),(10,q10)):
                    verify_artifacts({c['report']:dict(path=c['report'],bytes=Path(c['report']).stat().st_size,sha256=c['report_sha256'])})
                    prefix=dict(policy=policy,tag=tag,steps=step,time=c['time'],order=order,state_sha256=c['state_sha256'],
                        producer=c['report'],producer_sha256=c['report_sha256'],observation=c['path'],observation_sha256=c['sha256'])
                    for kind in ('continuum_error','initial_drift','initial_representation'):
                        for scope in ('full','central','constrained_boundary'):
                            for field,v in c['errors'][kind][scope]['fields'].items():
                                rows['fields'].append(dict(prefix,kind=kind,scope=scope,field=field,**{k:v[k] for k in ('L2','sampled_Linf','measure','normalized_L2','normalized_sampled_Linf')}))
                    rows['interfaces'].append(dict(prefix,**{k:v for k,v in c['errors']['interface'].items() if not isinstance(v,list)},
                        sampled_K_min=c['errors']['interface']['sampled_K'][0],sampled_K_max=c['errors']['interface']['sampled_K'][1],
                        cutoff_min=c['errors']['interface']['cutoff'][0],cutoff_max=c['errors']['interface']['cutoff'][1],
                        normalized_strong_divergence=c['errors']['strong_divergence']))
                    rows['residuals'].append(dict(prefix,uncancelled_spatial=json.dumps(c['spatial_residual']),actual_completed_stage=json.dumps(c['stage_defect'])))
                    if c.get('arithmetic'):rows['arithmetic'].append(dict(prefix,absolute_steps=c['arithmetic'].get('absolute_steps',step),
                        propagation=json.dumps(c['arithmetic']['accumulated_propagation_factors']),safety_factor=c['arithmetic']['safety_factor']))
            c=t['trajectory'][-1];d=read(c['report'])
            endpoints[policy][tag]=dict(absolute_steps=d['steps'],time=d['time'],dt=d['dt'],state_size=len(d['state']),
                report=file_identity(c['report']),checkpoint=artifacts(Path(c['checkpoint']).rglob('*')),signature=d['signature'],
                state_sha256=state_hash(d['state']),history_sha256=[state_hash(v) for v in d['native_history']],
                history_times=d['native_history_times'],fixed_reference_sha256=state_hash(d['fixed_reference']),
                completed_callback=dict(time=d['last_stage_time'],alpha=d['last_stage_alpha'],history_sha256=state_hash(d['last_stage_history'])))
        for name,e in s['execution'].items():
            if not e.get('report'):continue
            d,stage,allocations=execution_stage_binding(e['report']);attempts=d['pimple_attempts']
            fields=('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls',
                'bound_activations','additional_governing_equation_solves','additional_correction_or_restoration_solves')
            counts={key:sum(v.get(key,0) for v in attempts) for key in fields};equal=True
            for v in attempts:
                expected=[1790,187,1977,900,242] if len(d['state'])==4331 else [716,125,841,360,98]
                require([v[k] for k in ('velocity_size','pressure_gauge_size','assembled_flow_size','k_solve_size','epsilon_solve_size')]==expected,'native dimensions differ')
                require(v['scheduled_complete'] and all(v[k]==1 for k in fields[:5]),'scheduled solve/bound count differs')
                snap=read(v['k_snapshot']);equal=equal and state_hash(snap['raw'])==state_hash(snap['corrected'])
            rows['counts'].append(dict(policy=policy,case=name,attempts=len(attempts),raw_corrected_equal=equal,**counts))
            require(equal and counts['bound_activations']==0,'publication found active correction')
            rows['runs'].append(dict(policy=policy,case=name,absolute_steps=d['steps'],time=d['time'],actual_attempts=e['actual_attempts'],
                report=e['report']['path'],report_sha256=e['report']['sha256'],exit_code=e['exit_code'],capacity=d['dense_evolution_capacity']['effective']))
            stages[e['report']['path']]=dict(report=e['report'],steps=d['steps'],time=d['time'],
                last_stage_time=d['last_stage_time'],last_stage_alpha=d['last_stage_alpha'],last_stage_history_sha256=state_hash(d['last_stage_history']),
                callbacks=[{k:v[k] for k in ('history_sha256','incoming_sha256','time','alpha')} for v in attempts],
                retained_packets_sha256={k:digest(d[k]) for k in ('accepted_states','native_history','native_history_times','state_audits','wj_audit_packets')})
        for f in folder.rglob('*'):
            if f.is_file() and f.suffix!='.lock' and f.name!='summary.json':heavy[str(f)]=file_identity(f)
        write_once(out/(policy+'-phases.json'),components[policy])
        write_once(out/(policy+'-imports.json'),read(folder/'imports.json'))
        write_once(out/(policy+'-executions.json'),s['execution'])
        for path in quals:write_once(out/(policy+'-'+path.name),read(path))
        for version in ('retained-replay','retained-replay-v3','retained-replay-v5'):
            write_once(out/(policy+'-'+version+'.json'),read(root/version/(policy+'.json')))
    for c in read(root/'command-receipts.json'):
        rows['resources'].append(dict(label=c['label'],seconds=c['seconds'],exit_code=c['exit_code'],
            whole_process_peak_rss_bytes=c['whole_process_peak_rss_bytes'],parent_seconds=c['unaccounted_parent_seconds'],log_sha256=c['log']['sha256']))
    for name,data in rows.items():table(out/(name+'.csv'),data)
    for name in ('correction.json','native-equivalence.json','command-receipts.json','preflight.json','setup-cost-allowance.json',
                 'prefix-row-cache-index-EP.json','prefix-row-cache-index-Budget.json'):
        copy_exact(root/name,out/name)
    for f in (root/'study').rglob('*'):
        if f.is_file() and f.suffix!='.lock':heavy[str(f)]=file_identity(f)
    write_once(out/'endpoints.json',endpoints);write_once(out/'accepted-stage-lineage.json',stages)
    write_once(out/'final-target-coverage.json',coverage)
    for f in (root/'native-row-cache').rglob('*.json'):heavy[str(f)]=file_identity(f)
    write_once(out/'heavy-index.json',heavy)
    for folder in ('joint','resources'):
        for f in (root/folder).rglob('*.json'):
            write_once(out/folder/f.relative_to(root/folder),read(f))
    source=root/'source-v8'
    manifest=dict(schema='WJ7PublicationManifestV1',governing_source_revision=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip(),
        governing_source_tree=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD^{tree}'],text=True).strip(),
        publication_code=artifacts([Path(__file__)]),parent=file_identity(a.parent/'export-v1/transport.json'),
        native_equivalence=file_identity(root/'native-equivalence.json'),summaries=summaries,heavy_index=file_identity(out/'heavy-index.json'),
        components=components,scope='finite serial flat2D Dirichlet strip; T3*u*/H=.002048; not stationarity or global stability')
    write_once(out/'manifest.json',manifest)
    write_once(out/'receipt.json',dict(schema='WJ7CompactReceiptV1',manifest=file_identity(out/'manifest.json'),files=artifacts(f for f in out.rglob('*') if f.suffix!='.lock'),governing_attempts=0))
    print(out/'receipt.json')


if __name__=='__main__':main()
