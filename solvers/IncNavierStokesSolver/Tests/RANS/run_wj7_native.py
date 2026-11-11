#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Immutable imports, streamed replay and one joint WJ-7 continuation cadence."""
import argparse
import contextlib
import copy
import gc
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import uuid

from run_wj6_native import Study as WJ6Study
from wj_immutable import (Ledger,read,atomic,write_once,file_identity,file_hash,
    artifacts,verify_artifacts,digest,runtime,execution_intent)
from wj6_resources import DEFAULT,run_limited,volume,scientific_workers
from wj7_lineage import correct_stage_map,verify_import,bind_producer
from wj7_qualification import qualify,joint_ready,final_gate,observation_steps
from wj4_measure import AMP,BUDGET
from wj2_cases import normal_knots
from run_cls_epawf_picard import THREADS,require,state_hash
from wj_qualification import aggregate,component

SCRIPTS=Path(__file__).resolve().parent
PARENT_TRANSPORT='74649d94ba88f1ddb20e5f2ca4642a7a0bfb36da47290b195cf862cda79017df'


@contextlib.contextmanager
def session(ledger,label):
    key=label+'-'+uuid.uuid4().hex;ledger.begin_session(key)
    try:yield
    finally:ledger.finish_session(key)


def sealed_prefix(parent,policy):
    s=read(parent/'compact-v3'/(policy+'.json'))
    verify_artifacts({s['order10_row_supplement']['result']['path']:s['order10_row_supplement']['result']})
    supplement=read(s['order10_row_supplement']['result']['path'])
    verify_artifacts({supplement['records']['path']:supplement['records']})
    return s,read(supplement['records']['path']),supplement


def verify_parent(parent,source):
    transport=parent/'export-v1/transport.json'
    require(file_hash(transport)==PARENT_TRANSPORT,'sealed WJ-6 transport anchor differs')
    for item in read(transport)['files'].values():
        f=item['original'];verify_artifacts({f['path']:f})
    manifest=read(parent/'compact-v3/manifest.json')
    actual={name:runtime(Path(r['path'])) for name,r in manifest['runtime'].items()}
    require(actual==manifest['runtime'],'native executable/dependency closure differs')
    old=read(parent/'producer-source-content.json');checked={};registration=[]
    for name,wanted in old.items():
        # Match the sealed producer's declared 1800-input selection exactly.
        if not (name.endswith(('.cpp','.h','.hpp')) or name.endswith('CMakeLists.txt')):continue
        got=file_hash(source/name)
        if got!=wanted:
            require(name=='solvers/IncNavierStokesSolver/Tests/RANS/CMakeLists.txt','native source input changed: '+name)
            before=(parent/'source-v2'/name).read_text();after=(source/name).read_text()
            addition='''    add_test(NAME WJStageLineageContinuationMetadata COMMAND ${Python3_EXECUTABLE} -B
        ${CMAKE_CURRENT_LIST_DIR}/test_wj7_metadata.py)
    set_tests_properties(WJStageLineageContinuationMetadata PROPERTIES PROCESSORS 1 TIMEOUT 90
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
'''
            require(after.replace(addition,'')==before,'change extends beyond metadata-only CTest registration')
            registration.append(dict(path=name,before=wanted,after=got,scope='only add_test metadata; no native target/flags/ABI change'))
        checked[name]=got
    require(len(checked)==1800,'native source coverage differs')
    return dict(status='PASS',parent_manifest=file_identity(parent/'compact-v3/manifest.json'),
        sealed_transport=file_identity(transport),runtime=actual,native_inputs=checked,
        metadata_only_registration=registration,reused_regressions=file_identity(parent/'compact-v3/regression-verification.json'),
        transactions={p:file_identity(parent/'study'/p/'transaction-result.json') for p in ('EP','Budget')})


class Study(WJ6Study):
    def __init__(self,a,policy):
        self.a=copy.copy(a);self.a.policy=policy;self.a.output=a.root/'study'
        self.root=self.a.output/policy;self.root.mkdir(parents=True,exist_ok=True)
        self.summary=self.root/'summary.json'
        self.ledger=Ledger(a.root,policy,cap=800)
        self.s=read(self.summary) if self.summary.exists() else dict(schema='WJ7StudyV1',policy=policy,execution={},components={},trajectories={})
        self.runtime={str(p):runtime(p) for p in (a.solver,a.probe,Path(sys.executable))}
        self.code=dict(native_case_inputs=artifacts([SCRIPTS/n for n in ('wj2_cases.py','cls_epawf_m4_cases.py','cls_epawf_m3_cases.py','EquilibriumWallFlat.xml')]),
            execution_engine=artifacts([Path(__file__),SCRIPTS/'wj_immutable.py',SCRIPTS/'wj6_resources.py',SCRIPTS/'wj7_lineage.py']))
        self.measure_code=artifacts([SCRIPTS/n for n in ('run_wj7_native.py','run_wj6_native.py','run_wj5_native.py','wj3_measure.py','wj4_measure.py',
            'wj5_measure.py','wj6_qualification.py','wj6_arithmetic.py','wj7_qualification.py','wj7_cache.py','wj4_qualification.py','wj_qualification.py','test_wj2_native.py','validate_wj1.py')])
        self.measure_version=digest(self.measure_code)[:12]
        self.resource_file=self.root/'resource-design.json'
        self.design=dict(schema='WJ7FrozenRunIntentV1',policy=policy,runtime=self.runtime,source=artifacts([SCRIPTS.parents[3]/'docs/plans/wj7-contract.md']),
            code=self.code,parent=file_identity(a.root/'correction.json'),equivalence=file_identity(a.root/'native-equivalence.json'),
            retained_replay=file_identity(a.root/a.replay_name/(policy+'.json')),
            prefix_row_cache=file_identity(a.root/('prefix-row-cache-index-'+policy+'.json')),
            meshes=dict(L2=normal_knots(2)),xknots=[[0.,1.],[0.,.5,1.]],amplitudes=AMP,budgets=BUDGET,sensitivity_fraction=.1,
            start=dict(coarse=32,fine=64,baseline=32),targets=[64,96,128],dt=dict(coarse=3.2e-5,fine=1.6e-5,baseline=3.2e-5),
            final_time=.004096,absolute_steps=True,drift_origin='original WJ-6 native CP0/fixed_reference',observation_orders=[8,10],xorders=[4,6],
            maximum_native_steps_per_invocation=32,attempt_cap=800,seconds_cap=10800,closure_hold_seconds=600)
        write_once(self.root/'design.json',self.design)
        write_once(self.resource_file,dict(schema='WJ7FiniteResourcesV1',limits=DEFAULT,capacity=4331,native_block_hard_seconds=1500,
            design=file_identity(self.root/'design.json'),native_equivalence=file_identity(a.root/'native-equivalence.json'),headroom=1.25))
        self.s['components']['transactions']=component('PASS','reused unchanged-runtime WJ-6 actual-capacity transaction',
            evidence=file_identity(a.parent/'study'/policy/'transaction-result.json'))
        self.s['design_sha256']=file_hash(self.root/'design.json');self.save()

    def execute(self,folder,argv,inputs,settings,steps=0,outputs=(),measurement=False):
        argv=list(map(str,argv));folder.mkdir(parents=True,exist_ok=True)
        require(steps<=32,'native invocation exceeds predeclared32-step limit')
        kind='arithmetic' if settings.get('operation')=='passive_arithmetic' else 'observation' if measurement else 'native'
        intent=execution_intent(argv,[*inputs,self.resource_file],self.measure_code if measurement else self.code,self.runtime,dict(settings,resource_design=file_identity(self.resource_file)))
        write_once(folder/'execution-intent.json',intent);receipt=folder/'receipt.json'
        for value in self.runtime.values():
            require(file_hash(value['path'])==value['sha256'],'runtime changed')
            for f,h in value['libraries'].items():require(file_hash(f)==h,'runtime dependency changed')
        if receipt.exists():
            r=read(receipt);require(r['intent_sha256']==digest(intent),'cached requested intent differs')
            verify_artifacts(r['artifacts']);require(r['exit_code']==0,'failed/unknown launched execution cannot be reused');return r
        require(not (folder/'command.log').exists(),'unreconciled attempted execution; no automatic physical retry')
        seconds=DEFAULT[kind+'_seconds']
        if steps:
            projection=read(self.a.root/'resources'/self.a.policy/('before-'+str(self.a.end)+'.json'))
            seconds=max(120.,projection['native_rates'][settings['role']]*steps+120.)
            require(seconds<=1500,'RESOURCE_LIMIT: measured native timeout exceeds1500; never truncate')
        totals=self.ledger.totals()
        require(totals['seconds']+seconds+600<=10800,'RESOURCE_LIMIT: invocation would consume closure reserve')
        require(volume(self.root)+DEFAULT['invocation_output_bytes']<=DEFAULT['policy_output_bytes'],'RESOURCE_LIMIT: output reservation')
        label=str(folder.resolve());self.ledger.reserve(label,intent,steps,seconds)
        inspection=time.monotonic();not_launched=False
        try:metrics=run_limited(argv,folder,folder/'command.log',DEFAULT,kind,seconds)
        except RuntimeError as error:
            if 'RESOURCE_LIMIT:' not in str(error) or (folder/'command.log').exists() or (folder/'resource-time.txt').exists():raise
            not_launched=True;(folder/'command.log').write_text(str(error)+'\n')
            metrics=dict(exit_code=125,limit_reason=str(error),seconds=time.monotonic()-inspection,peak_rss_bytes=0,
                output_bytes=volume(folder),threads=THREADS,not_launched=True)
        report=folder/'case.evolution.json';d=read(report) if report.exists() else None
        actual=0 if not_launched else len(d['pimple_attempts']) if d else 0 if not steps else None
        files=[folder/'command.log',folder/'resource-time.txt',*outputs]
        if d and d.get('checkpoint'):files.extend(p for p in Path(d['checkpoint']).rglob('*') if p.is_file())
        r=dict(argv=argv,intent_sha256=digest(intent),exit_code=metrics['exit_code'],seconds=metrics['seconds'],resource=metrics,
            failure=metrics['limit_reason'],requested_steps=steps,actual_attempts=actual,artifacts=artifacts(files),threads=THREADS)
        write_once(receipt,r);self.ledger.finish(label,actual,r['seconds'],str(receipt));self.save();return r

    def observe(self,tag,end,d,folder,arithmetic):
        nx=1 if tag=='baseline' else 2;name=tag+'-'+str(end)+'-q8'
        m,c=self.measure(name,d,folder,nx=nx,assembly=True,diagnostics=True)
        c.update(tag=tag,dt=d['dt'],native_receipt=file_identity(folder/'receipt.json'))
        require(m['represented_initial_sha256']==self.s['origins'][tag]['represented_initial_sha256'],'original CP0 drift reference changed')
        if arithmetic:
            c['arithmetic']=self.arithmetic(name,d,m)
            c['arithmetic']['absolute_steps']=d['steps']
        self.cache_row(m,c,tag,end,8)
        del m;gc.collect()
        m,q=self.measure(tag+'-'+str(end)+'-q10x6',d,folder,nx=nx,order=10,xorder=6,diagnostics=True)
        q.update(tag=tag,dt=d['dt'],native_receipt=file_identity(folder/'receipt.json'))
        require(m['represented_initial_sha256']==self.s['origins'][tag]['represented_initial_sha256'],'order10 original CP0 changed')
        self.cache_row(m,q,tag,end,10)
        del m;gc.collect()
        changes=[abs(c['errors'][k][s]['fields'][f]['L2']-q['errors'][k][s]['fields'][f]['L2']) for k in ('continuum_error','initial_drift')
            for s in ('full','central','constrained_boundary') for f in AMP]
        c['quadrature']=component('PASS' if max(changes)<=1e-10 else 'UNRESOLVED',maximum=max(changes),limit=1e-10,independent=file_identity(q['path']))
        for rec,order in ((c,8),(q,10)):
            out=self.root/('measure-'+self.measure_version)/(tag+'-'+str(end)+'-qualified'+str(order)+'.json')
            write_once(out,rec);rec['qualified_record']=file_identity(out)
        self.s.setdefault('order10',{})[tag+'/'+str(end)]=q
        return c

    def cache_row(self,m,c,tag,step,order):
        from wj7_cache import build
        out=self.root/'native-row-cache'/(tag+'-'+str(step)+'-'+str(order)+'.json')
        write_once(out,build(m,c));c['row_cache']=file_identity(out)

    def cache_index(self,end):
        old=read(self.a.root/('prefix-row-cache-index-'+self.a.policy+'.json'))['records']
        return {key:c.get('row_cache') or old[key[0]+'/'+str(key[1])+'/'+str(key[2])] for key,c in self.index(end).items()}

    def import_endpoints(self):
        old,q10,supp=sealed_prefix(self.a.parent,self.a.policy);correction=read(self.a.root/'correction.json')
        require(old['components']['primary']['status']=='PASS' and supp['gate']['status']=='PASS','sealed Tx gates incomplete')
        self.s['trajectories']=copy.deepcopy(old['trajectories']);self.s['order10']=q10;self.s['origins']={};imports={}
        for tag,t in self.s['trajectories'].items():
            c=t['trajectory'][-1];report,identity=bind_producer(c,correction)
            initial,initial_identity=bind_producer(t['initial'],correction)
            origin_measure=read(t['initial']['path'])['measurements']
            self.s['origins'][tag]=dict(producer=initial_identity,checkpoint=artifacts(Path(t['initial']['checkpoint']).rglob('*')),
                fixed_reference_sha256=state_hash(initial['fixed_reference']),represented_initial_sha256=origin_measure['represented_initial_sha256'],
                native_time=initial['time'],absolute_steps=initial['steps'])
            del origin_measure;gc.collect()
            imported,r,folder=self.native('import-'+tag,cp=c['checkpoint'],dt=t['dt'],nx=t['nx'],role='import')
            self.valid(imported,r);check=verify_import(report,imported,initial)
            imports[tag]=dict(check=check,original_producer=identity,preceding_completed_stage=correction['entries'][identity['path']],
                checkpoint=artifacts(Path(c['checkpoint']).rglob('*')),readback=file_identity(folder/'case.evolution.json'),
                original_reference=self.s['origins'][tag],readback_startup_is_not_completed_callback=True)
        write_once(self.root/'imports.json',dict(status='PASS',endpoints=imports,original_summary=file_identity(self.a.parent/'compact-v3'/(self.a.policy+'.json'))))
        self.s['components']['imports']=component('PASS',evidence=file_identity(self.root/'imports.json'));self.save()

    def index(self,end):
        result={}
        for tag,t in self.s['trajectories'].items():
            for c in [t['initial']]+t['trajectory']:
                step=int(c['steps'])
                if step>(end*2 if tag=='fine' else end):continue
                require((tag,step,8) not in result,'duplicate q8 record')
                result[tag,step,8]=c
                q=self.s['order10'][tag+'/'+str(step)]
                require((tag,step,10) not in result,'duplicate q10 record');result[tag,step,10]=q
        return result

    def replay(self,end):
        out=self.root/('qualification-'+str(end)+'.json')
        require(not out.exists(),'immutable qualification already exists')
        previous=read(self.root/('qualification-'+str(end-32)+'.json')) if end>32 else None
        started=time.time();result=qualify(self.index(end),self.a.policy,self.design,end,self.cache_index(end),previous)
        if end==32:
            old=sealed_prefix(self.a.parent,self.a.policy)[2];differences=[]
            for row in result['legacy_sensitivity']:
                order,group,step,scope,field,norm=row['key'];order=int(order);step=int(step)
                before=old['gate']['components'][group+'_order'+str(order)+'_sensitivity']['components'][str(step)]['components'][scope+'/'+field+'/'+norm]
                require(before['status']==row['status'],'replay sensitivity decision differs')
                differences.extend([abs(before['value']-row['value']),abs(before['uncertainty']-row['uncertainty'])])
            result['retained_replay']=dict(leaves=504,maximum_absolute_value_or_uncertainty_difference=max(differences),
                decisions_exact=True,original=file_identity(Path(old['records']['path']).parent/'qualification.json'))
            require(max(differences)<=1e-15,'arithmetic replay differs beyond declared evaluation tolerance; no budget relaxation')
        result['seconds']=time.time()-started;write_once(out,result)
        self.s['components']['prefix_'+str(end)]=result['gate'];self.s['last_qualification']=file_identity(out);self.save();return result

    def row(self,end):
        start_time=time.monotonic();commands_before=self.ledger.totals()['command_seconds']
        for tag,target in (('coarse',end),('fine',2*end),('baseline',end)):
            t=self.s['trajectories'][tag];last=t['trajectory'][-1];start=int(last['steps']);cp=last['checkpoint']
            require(target>start,'already executed row cannot be relaunched')
            while start<target:
                stop=min(start+32,target);terminal=stop==target
                d,r,folder=self.native(tag+'-'+str(stop),stop-start,cp,t['dt'],t['nx'],arithmetic=terminal,role=tag)
                self.valid(d,r)
                self.s.setdefault('native_endpoints',{})[tag]=dict(steps=stop,checkpoint=d['checkpoint'],producer=file_identity(folder/'case.evolution.json'))
                self.save();cp=d['checkpoint'];start=stop
                if terminal:
                    c=self.observe(tag,stop,d,folder,True);t['trajectory'].append(c);self.save()
                    require(c['errors']['accuracy_gate']['status']=='PASS' and c['quadrature']['status']=='PASS',
                            'OBSTRUCTED: required cadence accuracy/quadrature failed; no further native advance')
        result=self.replay(end)
        elapsed=time.monotonic()-start_time;child=self.ledger.totals()['command_seconds']-commands_before
        parent=max(0,elapsed-child);cost=dict(end=end,whole_seconds=elapsed,child_seconds=child,parent_seconds=parent,
            prefix_qualification_seconds=result['seconds'],production_parent_seconds_per_pair=max(0,parent-result['seconds'])/3)
        write_once(self.root/('row-cost-'+str(end)+'.json'),cost);return result


def projection(a,policy,end):
    folder=a.root/'study'/policy;s=read(folder/'summary.json');ledger=Ledger(a.root,policy,cap=800)
    old=read(a.parent/'study'/policy/'T3-resource-projection.json')
    rates=dict(old['native_seconds_per_step_bounds'])
    current_receipts={}
    for name,e in s['execution'].items():
        verify_artifacts({e['receipt']['path']:e['receipt']});r=read(e['receipt']['path'])
        require(r['exit_code']==e['exit_code'] and r['actual_attempts']==e['actual_attempts'] and r['resource']==e['resource'],
                'projection progress differs from immutable execution receipt')
        current_receipts[name]=e['receipt']
    for tag in rates:
        current=[e['resource']['seconds']/e['actual_attempts']*1.25 for e in s['execution'].values()
            if e['role']==tag and e.get('actual_attempts') and e['exit_code']==0]
        if current:rates[tag]=max(rates[tag],max(current))
    remaining={tag:(256 if tag=='fine' else 128)-int(t['trajectory'][-1]['steps']) for tag,t in s['trajectories'].items()}
    rows=(128-(end-32))//32
    # Old native observation pair plus separately measured map export cost is
    # retained. We do NOT subtract all of the former q10 supplement.
    # Actual first-read native q8 and q10-with-maps costs, not the old extra
    # quadrature-only q10 read plus an entire duplicate repair campaign.
    prefix,q10,_=sealed_prefix(a.parent,policy);native_observation={}
    old_manifest=read(a.parent/'compact-v3/manifest.json');identity=old_manifest['artifact_index']
    verify_artifacts({identity['path']:identity});old_index=read(identity['path'])
    cost_inputs={}
    for tag,t in prefix['trajectories'].items():
        pairs=[]
        for c in [t['initial']]+t['trajectory']:
            q=q10[tag+'/'+str(int(c['steps']))]
            paths=[Path(v['path']).parent/'receipt.json' for v in (c,q)]
            for f in paths:
                require(str(f) in old_index,'missing sealed native row-export cost receipt')
                verify_artifacts({str(f):old_index[str(f)]});cost_inputs[str(f)]=old_index[str(f)]
            pairs.append(sum(read(f)['seconds'] for f in paths))
        native_observation[tag]=max(pairs)
    observation=sum(native_observation.values())*rows*1.25
    svd=sum(max(old['measured_svd_seconds'][tag]) for tag in rates)*rows*1.25
    replay=read(folder/('qualification-'+str(end-32)+'.json'))
    retained_count=len(replay['record_index'])/2
    # Existing prefix must be requalified each row under an enlarged envelope.
    expected_counts=[sum(len(observation_steps(tag,v)) for tag in rates) for v in range(end,129,32)]
    prefix_parent=sum(expected_counts)*replay['seconds']/retained_count*1.25
    # Full payload parsing/evaluation is directly timed by the replay. Retain
    # an explicit additional allowance for independent loads/diagnosis and
    # serialization, validated by the passive production replay below.
    profile=read(a.root/'production-profile'/ (policy+'.json'))
    observed_parent=[read(f)['production_parent_seconds_per_pair']+5 for f in folder.glob('row-cost-*.json')]
    production_parent=3*rows*max([profile['maximum_parent_pair_seconds'],*observed_parent])*1.25
    parent=prefix_parent+production_parent
    launches=sum(math.ceil(n/32) for n in remaining.values());launch_overhead=launches*5*1.25
    estimate=sum(remaining[t]*rates[t] for t in rates)+observation+svd+parent+launch_overhead+600
    totals=ledger.totals();fits=totals['seconds']+estimate<=10800 and totals['attempts']+sum(remaining.values())<=800
    archive=a.root/'resources'/policy/('before-'+str(end)+'-basis')
    write_once(archive/'ledger.json',read(ledger.file));write_once(archive/'native-export-receipts.json',cost_inputs)
    write_once(archive/'current-native-receipts.json',current_receipts)
    value=dict(schema='WJ7RemainingResourceProjectionV1',status='PASS' if fits else 'RESOURCE_LIMIT',before=end,
        native_rates=rates,remaining_steps=remaining,remaining_rows=rows,native_launches=launches,
        native_seconds=sum(remaining[t]*rates[t] for t in rates),observation_seconds=observation,svd_seconds=svd,
        streamed_prefix_parent_seconds=prefix_parent,production_parent_seconds=production_parent,additional_launch_overhead_seconds=launch_overhead,closure_seconds=600,
        headroom=1.25,totals=totals,projected_remaining_seconds=estimate,projected_total_seconds=totals['seconds']+estimate,
        measured_native_row_export_pairs=native_observation,production_profile=file_identity(a.root/'production-profile'/(policy+'.json')),
        basis=artifacts([a.parent/'study'/policy/'T3-resource-projection.json',folder/('qualification-'+str(end-32)+'.json'),
            archive/'ledger.json',archive/'native-export-receipts.json',archive/'current-native-receipts.json']),
        allocation=scientific_workers(),policy_bytes=volume(folder),hard_limits=DEFAULT)
    require(max(rates.values())*32+120<=1500,'RESOURCE_LIMIT: <=32-step timeout unsupported')
    write_once(a.root/'resources'/policy/('before-'+str(end)+'.json'),value);return value


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--parent',type=Path,required=True)
    p.add_argument('--solver',type=Path,required=True);p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--mode',choices=('prepare','prepare-caches','cached-replay','production-profile','import','replay','admit','cadence'),required=True)
    p.add_argument('--replay-name',default='retained-replay')
    p.add_argument('--end',type=int,choices=(32,64,96,128),default=32);a=p.parse_args()
    require(platform.system()=='Linux','Linux is authoritative');os.environ.update(THREADS)
    for name in ('root','parent','solver','probe'):setattr(a,name,getattr(a,name).resolve())
    a.root.mkdir(parents=True,exist_ok=True)
    if a.mode=='prepare':
        ls=[Ledger(a.root,policy,cap=800) for policy in ('EP','Budget')]
        with contextlib.ExitStack() as stack:
            for l in ls:stack.enter_context(session(l,'prepare'))
            eq=verify_parent(a.parent,SCRIPTS.parents[3]);write_once(a.root/'native-equivalence.json',eq)
            correct_stage_map(a.parent,a.root/'correction.json')
        print('semantic correction and native equivalence PASS');return
    if a.mode=='prepare-caches':
        from wj7_cache import build
        from wj7_qualification import load_measurement
        for policy in ('EP','Budget'):
            ledger=Ledger(a.root,policy,cap=800)
            with session(ledger,'prepare-caches'):
                old,q10,_=sealed_prefix(a.parent,policy);records={}
                for tag,t in old['trajectories'].items():
                    for c8 in [t['initial']]+t['trajectory']:
                        step=int(c8['steps'])
                        for order,c in ((8,c8),(10,q10[tag+'/'+str(step)])):
                            m=load_measurement(c);out=a.root/'native-row-cache'/policy/(tag+'-'+str(step)+'-'+str(order)+'.json')
                            write_once(out,build(m,c));records[tag+'/'+str(step)+'/'+str(order)]=file_identity(out)
                            del m;gc.collect()
                write_once(a.root/('prefix-row-cache-index-'+policy+'.json'),dict(records=records,
                    previous_full_replay=file_identity(a.root/'retained-replay'/(policy+'.json')),code=artifacts([SCRIPTS/'wj7_cache.py'])))
        return
    if a.mode in ('retained-replay','cached-replay'):
        for policy in ('EP','Budget'):
            ledger=Ledger(a.root,policy,cap=800)
            with session(ledger,'retained-replay'):
                old,q10,supp=sealed_prefix(a.parent,policy);index={}
                for tag,t in old['trajectories'].items():
                    for c in [t['initial']]+t['trajectory']:
                        step=int(c['steps']);index[tag,step,8]=c;index[tag,step,10]=q10[tag+'/'+str(step)]
                caches=None;previous=None
                if a.mode=='cached-replay':
                    info=read(a.root/('prefix-row-cache-index-'+policy+'.json'));verify_artifacts(info['code'])
                    caches={tuple([name.split('/')[0],int(name.split('/')[1]),int(name.split('/')[2])]):v for name,v in info['records'].items()}
                    verify_artifacts({info['previous_full_replay']['path']:info['previous_full_replay']});previous=read(info['previous_full_replay']['path'])
                start=time.time();result=qualify(index,policy,dict(amplitudes=AMP,budgets=BUDGET,sensitivity_fraction=.1),32,caches,previous)
                changes=[]
                for row in result['legacy_sensitivity']:
                    order,group,step,scope,field,norm=row['key'];order=int(order);step=int(step)
                    prior=supp['gate']['components'][group+'_order'+str(order)+'_sensitivity']['components'][str(step)]['components'][scope+'/'+field+'/'+norm]
                    require(prior['status']==row['status'],'retained replay decision differs')
                    changes.extend((abs(prior['value']-row['value']),abs(prior['uncertainty']-row['uncertainty'])))
                require(max(changes)<=1e-15,'retained arithmetic ordering differs; unchanged budgets required')
                result.update(seconds=time.time()-start,retained_replay=dict(leaves=504,decisions_exact=True,
                    maximum_absolute_value_or_uncertainty_difference=max(changes)),code=artifacts([SCRIPTS/'wj7_qualification.py',Path(__file__)]))
                write_once(a.root/a.replay_name/(policy+'.json'),result)
        print('both retained-data qualifications complete');return
    if a.mode=='production-profile':
        from run_wj5_native import diagnose,independent_loads,independent_weak,independent_near_h,flow_decomposition
        from wj7_qualification import complete_errors
        for policy in ('EP','Budget'):
            ledger=Ledger(a.root,policy,cap=800)
            with session(ledger,'production-profile'):
                old,q10,_=sealed_prefix(a.parent,policy);timings=[];checks=[]
                for tag,t in old['trajectories'].items():
                    c=t['trajectory'][-1];start=time.monotonic();raw=read(c['path']);m=raw['measurements'];d=read(c['report'])
                    errors=complete_errors(m,policy,c);require(errors==c['errors'],'production-profile accuracy changed')
                    dg=diagnose(m);require(dg==c['diagnosis'],'production-profile diagnosis changed')
                    ch=raw['assembly_checks'];ch['policy']=policy
                    independent_loads(ch,policy,d['offsets']);independent_weak(ch,m['weak_modes']);independent_near_h(ch,m);flow_decomposition(d,ch)
                    del raw,m,d,ch;gc.collect()
                    q=q10[tag+'/'+str(int(c['steps']))];raw=read(q['path']);m=raw['measurements']
                    require(complete_errors(m,policy,q)==q['errors'],'production-profile order10 changed');diagnose(m)
                    del raw,m;gc.collect();timings.append(time.monotonic()-start)
                    checks.append(dict(tag=tag,q8=file_identity(c['path']),q10=file_identity(q['path'])))
                value=dict(status='PASS',pair_seconds=timings,maximum_parent_pair_seconds=max(timings)+5.,
                    fixed_extra_serialization_verification_seconds=5.,inputs=checks,governing_attempts=0,
                    scope='measured complete passive parent work on actual retained fields; native export CPU accounted separately; plus25% in admission')
                write_once(a.root/'production-profile'/(policy+'.json'),value)
        return
    results={}
    if a.mode=='cadence':
        require(a.end in (64,96,128),'scientific cadence invalid')
        previous=read(a.root/'joint'/('barrier-'+str(a.end-32)+'.json'))
        admission=read(a.root/'joint'/('admission-'+str(a.end)+'.json'))
        require(previous['status']=='PASS' and admission['status']=='PASS','joint previous barrier/admission did not pass')
    for policy in ('EP','Budget'):
        s=Study(a,policy)
        with session(s.ledger,a.mode+'-'+str(a.end)):
            if a.mode=='import':s.import_endpoints();results[policy]=s.s['components']['imports']
            elif a.mode=='replay':
                retained=a.root/a.replay_name/(policy+'.json')
                if retained.exists():
                    value=read(retained);verify_artifacts(value['code'])
                    require(value['gate']['status']=='PASS' and len(s.index(32))==len(value['record_index']),'retained replay/import coverage differs')
                    for (tag,step,order),c in s.index(32).items():
                        old=value['record_index'][tag+'/'+str(step)+'/'+str(order)]
                        require(all(old[k]==c[k] for k in ('path','sha256','report','report_sha256')),'retained replay/import record differs')
                    write_once(s.root/'qualification-32.json',value);s.s['components']['prefix_32']=value['gate'];s.save();results[policy]=value
                else:results[policy]=s.replay(32)
            elif a.mode=='admit':results[policy]=projection(a,policy,a.end)
            elif a.mode=='cadence':
                try:results[policy]=s.row(a.end)
                except Exception as error:
                    status='RESOURCE_LIMIT' if 'RESOURCE_LIMIT' in str(error) else 'OBSTRUCTED'
                    results[policy]=dict(end=a.end,gate=component(status,str(error)))
                    write_once(a.root/'joint'/('failure-'+str(a.end)+'-'+policy+'.json'),dict(result=results[policy],execution=s.s['execution']))
                if results[policy]['gate']['status']!='PASS':break
            s.save()
    if a.mode in ('replay','cadence'):
        end=32 if a.mode=='replay' else a.end
        gate=component('PASS' if joint_ready(results,end) else aggregate({k:r['gate'] for k,r in results.items()})['status'],
            end=end,policies={k:dict(end=r['end'],gate=r['gate'],evidence=file_identity(a.root/'study'/k/('qualification-'+str(end)+'.json'))
                if (a.root/'study'/k/('qualification-'+str(end)+'.json')).exists() else None) for k,r in results.items()})
        if set(results)!={'EP','Budget'} and gate['status']=='PASS':gate['status']='NOT_RUN'
        write_once(a.root/'joint'/('barrier-'+str(end)+'.json'),gate)
    elif a.mode=='admit':
        write_once(a.root/'joint'/('admission-'+str(a.end)+'.json'),aggregate({p:component(v['status'],projection=file_identity(a.root/'resources'/p/('before-'+str(a.end)+'.json'))) for p,v in results.items()}))
    print(json.dumps({k:(r.get('gate',r)).get('status') for k,r in results.items()}))


if __name__=='__main__':main()
