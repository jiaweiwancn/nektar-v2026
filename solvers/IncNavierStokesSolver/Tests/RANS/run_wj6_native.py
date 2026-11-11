#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Resource-qualified uniform L2/Nx2 WJ pilot, Tx pair and conditional T3."""
import argparse
import copy
import json
import math
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time
import uuid
import xml.etree.ElementTree as ET

from run_wj5_native import Study as WJ5Study
from wj2_cases import case_xml,normal_knots
from wj3_measure import checkpoint_payload
from wj4_measure import comparison,AMP,BUDGET
from wj4_qualification import sensitivity_sequence
from wj5_measure import boundary_uncertainty
from wj6_qualification import paired_envelope,required_accuracy,common_extension_authorized,fresh_startup_gate
from wj6_resources import DEFAULT,run_limited,volume
from wj_immutable import (Ledger,read,write_once,atomic,artifacts,verify_artifacts,file_hash,
    file_identity,runtime,digest,xml_identity,execution_intent)
from wj_qualification import aggregate,component
from cls_epawf_m3_cases import set_info,write
from run_cls_epawf_picard import THREADS,require,state_hash
from compare_cls_epawf_pimple import same_history
from test_wj2_native import digest_object
from verify_wj4_transactions import verify as verify_transaction

SCRIPTS=Path(__file__).resolve().parent


def meta(path):return json.loads(Path(path).read_text())


class Study(WJ5Study):
    def __init__(self,a):
        require(platform.system()=='Linux','native Linux is authoritative')
        start=time.time();os.environ.update(THREADS);self.a=a
        self.root=a.output.resolve()/a.policy;self.root.mkdir(parents=True,exist_ok=True)
        self.summary=self.root/'summary.json';self.s=meta(self.summary) if self.summary.exists() else dict(schema='WJ6StudyV1',policy=a.policy,execution={},components={},trajectories={})
        self.runtime={str(p):runtime(p) for p in (a.solver,a.probe,Path(sys.executable))}
        self.code=dict(native_case_inputs=artifacts([SCRIPTS/n for n in ('wj2_cases.py','cls_epawf_m4_cases.py','cls_epawf_m3_cases.py','EquilibriumWallFlat.xml')]),
            execution_engine=artifacts([Path(__file__),SCRIPTS/'wj_immutable.py',SCRIPTS/'wj6_resources.py']))
        self.measure_code=artifacts([SCRIPTS/n for n in ('run_wj6_native.py','run_wj5_native.py','wj3_measure.py','wj4_measure.py','wj5_measure.py','wj6_qualification.py','wj6_arithmetic.py','wj4_qualification.py','wj_qualification.py','test_wj2_native.py','validate_wj1.py')])
        self.measure_version=digest(self.measure_code)[:12]
        self.resource=dict(schema='WJ6ResourceDesignV1',policy=a.policy,limits=DEFAULT,capacity=4331,
            package_attempts=800,package_seconds=10800,closure_hold_seconds=600,transaction_hold_attempts=9,
            measured_projection_headroom=1.25,native_block_hard_seconds=1500,
            native_code=self.code,runtime=self.runtime,contract=artifacts([SCRIPTS.parents[3]/'docs/plans/wj6-contract.md']))
        write_once(self.root/'resource-design.json',self.resource)
        self.resource_file=self.root/'resource-design.json'
        if a.mode in ('primary','transactions'):
            self.resource_file=self.root/'Tx-resource-projection.json'
            require(self.resource_file.exists(),'Tx needs a frozen measured resource projection')
        elif a.mode=='extend':
            self.resource_file=self.root/'T3-resource-projection.json'
            require(self.resource_file.exists(),'T3 needs a frozen measured resource projection')
        self.design=dict(schema='WJ6FrozenDesignV1',policy=a.policy,runtime=self.runtime,governing_code=self.code,
            resource_design=file_identity(self.root/'resource-design.json'),mesh=normal_knots(2),meshes=dict(L2=normal_knots(2)),
            xknots=[[0.,1.],[0.,.5,1.]],amplitudes=AMP,budgets=BUDGET,sensitivity_fraction=.1,
            attempt_cap_per_policy=800,wall_cap_seconds_per_policy=10800,source_order=16,basis_points=[8,8],
            pilot=dict(coarse=4,fine=8,baseline=4),primary=dict(T=.001024,coarse=32,fine=64,baseline=32),
            extension=dict(T=.004096,coarse=128,fine=256,baseline=128),
            scientific_cadence=dict(coarse=32,fine=64,baseline=32),arithmetic=dict(coarse=[4,16,32,64,128],fine=[8,32,64,128,256],baseline=[4,16,32,64,128]),
            boundary='reference differential Dirichlet both laterals/top on [0,1]x[0,2]; no-slip ground, modeled zero ground k flux',
            fresh_same_runtime_baseline=True,observations=dict(order=8,xorder=4,independent_order=10,independent_xorder=6,absolute_floor=1e-10,all_faces=True))
        write_once(self.root/'design.json',self.design)
        self.ledger=Ledger(Path(os.environ['WJ4_PACKAGE_ROOT']),a.policy,cap=800)
        self.session='wj6-'+uuid.uuid4().hex;self.ledger.begin_session(self.session,start)
        self.s['design_sha256']=file_hash(self.root/'design.json')
        self.s['source_revision']=subprocess.check_output(['git','-C',str(SCRIPTS.parents[3]),'rev-parse','HEAD'],text=True).strip()
        for key in ('pilot','resource_projection','primary','transactions','extension_resource','extension'):
            self.s['components'].setdefault(key,component('NOT_RUN','phase has not executed'))
        self.save()

    def execute(self,folder,argv,inputs,settings,steps=0,outputs=(),measurement=False):
        argv=list(map(str,argv));folder.mkdir(parents=True,exist_ok=True)
        kind=settings.get('resource_kind','observation' if measurement else 'native')
        intent=execution_intent(argv,[*inputs,self.resource_file],self.measure_code if measurement else self.code,self.runtime,
            dict(settings,resource_design=file_identity(self.resource_file)))
        write_once(folder/'execution-intent.json',intent);receipt=folder/'receipt.json'
        for value in self.runtime.values():
            require(file_hash(value['path'])==value['sha256'],'frozen runtime changed')
            for n,h in value['libraries'].items():require(file_hash(n)==h,'frozen runtime dependency changed')
        if receipt.exists():
            r=meta(receipt);require(r['intent_sha256']==digest(intent),'cached execution intent differs');verify_artifacts(r['artifacts']);return r
        require(not (folder/'command.log').exists(),'unreconciled command retained; choose a new explicitly versioned observer/attempt')
        limits=dict(DEFAULT)
        projection=meta(self.resource_file)
        seconds=limits[kind+'_seconds']
        if kind=='native' and steps and 'native_seconds_per_step_bound' in projection:
            rate=projection.get('native_seconds_per_step_bounds',{}).get(settings.get('role'),projection['native_seconds_per_step_bound'])
            seconds=min(1500.,max(120.,rate*steps+120.))
        totals=self.ledger.totals()
        hold=600. # Explicit closure headroom persists across every mode/resume.
        transaction_pending=self.s['components']['transactions']['status']!='PASS' and self.a.mode!='transactions'
        if transaction_pending:
            hold+=9*projection.get('native_seconds_per_step_bounds',{}).get('coarse',projection.get('native_seconds_per_step_bound',100.))
        require(totals['attempts']+steps+(9 if transaction_pending else 0)<=800,
                'RESOURCE_LIMIT: transaction unit reservation would be consumed')
        require(totals['seconds']+seconds+hold<=10800,'RESOURCE_LIMIT: invocation plus closure reserve exceeds package seconds')
        require(volume(self.root)<=limits['policy_output_bytes'],'RESOURCE_LIMIT: policy output ceiling')
        label=str(folder.resolve());self.ledger.reserve(label,intent,steps,seconds)
        not_launched=False
        inspection_start=time.monotonic()
        try:metrics=run_limited(argv,folder,folder/'command.log',limits,kind,seconds)
        except RuntimeError as error:
            # A denied worker reservation precedes Popen. Preserve and charge
            # its measured inspection time, but do not invent a native attempt.
            if 'RESOURCE_LIMIT:' not in str(error) or (folder/'command.log').exists() or (folder/'resource-time.txt').exists():raise
            not_launched=True
            with (folder/'command.log').open('x') as stream:stream.write(str(error)+'\n')
            metrics=dict(exit_code=125,limit_reason=str(error),seconds=time.monotonic()-inspection_start,peak_rss_bytes=0,
                peak_process_threads=0,output_bytes=volume(folder),timeout_seconds=seconds,
                limits=limits,allocation=None,threads=THREADS,not_launched=True)
        report=folder/'case.evolution.json';d=read(report) if report.exists() else None
        actual=0 if not_launched else len(d['pimple_attempts']) if d and 'pimple_attempts' in d else 0 if not steps else None
        files=[folder/'command.log',folder/'resource-time.txt',*outputs]
        if d and d.get('checkpoint') and Path(d['checkpoint']).exists():files.extend(p for p in Path(d['checkpoint']).rglob('*') if p.is_file())
        r=dict(argv=argv,intent_sha256=digest(intent),exit_code=metrics['exit_code'],seconds=metrics['seconds'],resource=metrics,
            failure='RESOURCE_LIMIT: '+str(metrics['limit_reason']) if metrics['limit_reason'] else None,
            requested_steps=steps,actual_attempts=actual,charged_attempts=actual if actual is not None else steps,
            artifacts=artifacts(files),threads=THREADS)
        write_once(receipt,r);self.ledger.finish(label,actual,r['seconds'],str(receipt));self.save();return r

    def native(self,name,steps=0,cp=None,dt=3.2e-5,nx=2,reject=0,arithmetic=False,role='primary'):
        folder=self.root/name;xml=folder/'case.xml';start=int(checkpoint_payload(cp)['steps']) if cp else 0
        options=dict(steps=int(steps),start=start,target=start+steps,restart=str(Path(cp).resolve()) if cp else None,
            dt=dt,nx=nx,zknots=self.design['meshes']['L2'],xknots=[i/nx for i in range(nx+1)],
            capacity=4331,resource_design_sha256=file_hash(self.resource_file),policy=self.a.policy,role=role,
            hard_reject=int(reject),arithmetic_time=(start+steps)*dt if arithmetic else None,source_order=16)
        if xml.exists():
            write_once(folder/'requested-intent.json',options)
            require(meta(folder/'xml-identity.json')==xml_identity(xml),'stale native XML')
        else:
            case_xml(folder,self.a.policy,steps,cp,dt=dt,level=2,xknots=options['xknots'],zknots=options['zknots'],hard_reject=reject)
            tree=ET.parse(xml).getroot()
            for key,value in dict(RANSWallDenseEvolutionCapacity='4331',RANSWallDenseEvolutionResourceDesign=options['resource_design_sha256'],
                    RANSWallPimpleStopOnBound='On',RANSWallPimpleSnapshotDirectory=str(folder/'pimple-states')).items():set_info(tree,key,value)
            if arithmetic:
                set_info(tree,'RANSWallPimpleArithmetic','On');set_info(tree,'RANSWallPimpleArithmeticTime',format(options['arithmetic_time'],'.17g'))
            write(tree,xml);write_once(folder/'requested-intent.json',options);write_once(folder/'xml-identity.json',xml_identity(xml))
        report=folder/'case.evolution.json'
        inputs=[xml]+([p for p in Path(cp).rglob('*') if p.is_file()] if cp else [])
        r=self.execute(folder,[self.a.solver,xml],inputs,options,steps,[report]);d=read(report) if report.exists() else None
        if d and r['exit_code']==0:
            require(d['steps']==start+steps and len(d['pimple_attempts'])==steps and abs(d['time']-(start+steps)*dt)<=2e-15,'native target mismatch')
            require(d['dense_evolution_capacity']['effective']==4331 and d['dense_evolution_capacity']['resource_design_sha256']==options['resource_design_sha256'],'native capacity/resource report differs')
        self.s['execution'][name]=dict(receipt=file_identity(folder/'receipt.json'),report=file_identity(report) if report.exists() else None,
            exit_code=r['exit_code'],actual_attempts=r['actual_attempts'],role=role,resource=r['resource'])
        self.save();return d,r,folder

    def measure(self,*args,order=8,xorder=4,**kwargs):
        return super().measure(*args,order=order,xorder=xorder,mesh='L2',**kwargs)

    def arithmetic(self,name,d,m):
        folder=self.root/('measure-'+self.measure_version)/(name+'-svd');out=folder/'arithmetic.json'
        producer=Path(d['checkpoint']).parent
        observation=self.root/('measure-'+self.measure_version)/name/'result.json'
        inputs=[producer/'case.evolution.json',observation]
        for field in ('momentum','pressure','k','epsilon'):
            identity=d['pimple_attempts'][-1][field+'_arithmetic'];require(file_hash(identity['path'])==identity['sha256'],'native arithmetic packet changed')
            p=read(identity['path']);inputs.extend([Path(identity['path']),Path(p['matrix']['path']),Path(p['mass_or_correction_map']['path'])])
        r=self.execute(folder,[sys.executable,'-B',SCRIPTS/'wj6_arithmetic.py','--report',producer/'case.evolution.json','--measurement',observation,'--output',out],
            inputs,dict(operation='passive_arithmetic',resource_kind='arithmetic',state_sha256=m['state_sha256']),outputs=[out],measurement=True)
        require(r['exit_code']==0,'RESOURCE_LIMIT: arithmetic/SVD failed' if r['failure'] else 'arithmetic/SVD process failed')
        return read(out)

    def observe(self,tag,end,d,folder,arithmetic):
        name=tag+'-'+str(end)+'-q8';nx=1 if tag=='baseline' else 2
        m,c=self.measure(name,d,folder,nx=nx,assembly=True)
        c['tag']=tag;c['dt']=d['dt'];c['native_receipt']=file_identity(folder/'receipt.json')
        if arithmetic:self.observe_arithmetic(name,d,m,c)
        qm,qc=self.measure(tag+'-'+str(end)+'-q10x6',d,folder,nx=nx,order=10,xorder=6,diagnostics=False)
        changes={}
        for kind in ('continuum_error','initial_drift'):
            for scope in ('full','central','constrained_boundary'):
                for field,row in c['errors'][kind][scope]['fields'].items():
                    changes[kind+'/'+scope+'/'+field]=abs(row['L2']-qc['errors'][kind][scope]['fields'][field]['L2'])
        c['quadrature']=component('PASS' if max(changes.values())<=1e-10 else 'UNRESOLVED',
            'independent actual-mesh normal8/10 and horizontal4/6, absolute floor1e-10',changes=changes,independent=file_identity(Path(qc['path'])))
        c['unknowns']=len(d['state'])
        qualified=self.root/('measure-'+self.measure_version)/(name+'-qualified.json')
        write_once(qualified,c);c['qualified_record']=file_identity(qualified)
        return c

    def advance(self,tag,end):
        trajectories=self.s['trajectories'];t=trajectories[tag]
        known=[t['initial']]+t['trajectory'];last=known[-1];start=int(last['steps'])
        require(end>start,'continuation must advance the retained trajectory')
        arithmetic=end in self.design['arithmetic'][tag]
        d,r,f=self.native(tag+'-'+str(end),end-start,last['checkpoint'],t['dt'],t['nx'],arithmetic=arithmetic,role=tag)
        self.valid(d,r)
        for row in d['pimple_attempts']:
            expected=[4331,1790,187,1977,900,242] if t['nx']==2 else [2493,716,125,841,360,98]
            actual=[len(d['state'])]+[row[k] for k in ('velocity_size','pressure_gauge_size','assembled_flow_size','k_solve_size','epsilon_solve_size')]
            require(actual==expected,'actual native solve dimensions differ')
        c=self.observe(tag,end,d,f,arithmetic);t['trajectory'].append(c)
        if c['errors']['accuracy_gate']['status']!='PASS':t.setdefault('accuracy_events',[]).append(dict(step=end,gate=c['errors']['accuracy_gate']))
        self.save();return c

    def pilot(self):
        for tag,dt,nx,end in (('coarse',3.2e-5,2,4),('fine',1.6e-5,2,8),('baseline',3.2e-5,1,4)):
            if tag not in self.s['trajectories']:
                d,r,f=self.native(tag+'-cp0',dt=dt,nx=nx,role=tag);self.valid(d,r)
                c=self.observe(tag,0,d,f,False)
                self.s['trajectories'][tag]=dict(initial=c,trajectory=[],dt=dt,nx=nx);self.save()
            if int(self.s['trajectories'][tag]['trajectory'][-1]['steps'] if self.s['trajectories'][tag]['trajectory'] else 0)<end:self.advance(tag,end)
        a,b=(checkpoint_payload(self.s['trajectories'][tag]['initial']['checkpoint']) for tag in ('coarse','fine'))
        require(state_hash(a['state'])==state_hash(b['state']) and [state_hash(v) for v in a['native_history']]==[state_hash(v) for v in b['native_history']],'fresh dt pair changed initial physical state/history')
        self.s['fresh_startup']=dict(exact_physical_primary_gauge_history=True,native_labels_differ=a['native_history_times']!=b['native_history_times'],
            checkpoints=artifacts([Path(self.s['trajectories'][tag]['initial']['checkpoint'])/'State.bin' for tag in ('coarse','fine')]))
        require(self.s['fresh_startup']['native_labels_differ'],'fresh dt pair lost native distinct time labels')
        self.s['components']['pilot']=self.qualify('pilot');self.save()

    def cost(self,phase='Tx'):
        target=self.design['primary' if phase=='Tx' else 'extension'];headroom=1.25
        rows=[r for r in self.s['execution'].values() if r.get('actual_attempts',0)]
        rates=[r['resource']['seconds']/r['actual_attempts'] for r in rows if r['exit_code']==0]
        require(rates,'missing actual capacity pilot cost')
        rate=max(rates)*headroom
        rates_by_tag={tag:max(r['resource']['seconds']/r['actual_attempts'] for r in rows if r['exit_code']==0 and r['role']==tag)*headroom
                      for tag in ('coarse','fine','baseline')}
        observations={};svd={}
        for tag,t in self.s['trajectories'].items():
            observations[tag]=[];svd[tag]=[]
            for c in [t['initial']]+t['trajectory']:
                q8=meta(Path(c['path']).parent/'receipt.json')['seconds']
                q10=meta(Path(c['quadrature']['independent']['path']).parent/'receipt.json')['seconds']
                observations[tag].append(q8+q10)
            for c in t['trajectory']:
                if c.get('arithmetic'):
                    folder=Path(c['path']).parent.with_name(Path(c['path']).parent.name+'-svd')
                    svd[tag].append(meta(folder/'receipt.json')['seconds'])
        remaining_by_tag={tag:max(0,int(target[tag])-int(t['trajectory'][-1]['steps'])) for tag,t in self.s['trajectories'].items()}
        remaining=sum(remaining_by_tag.values())
        remaining_observations=dict(coarse=3,fine=3,baseline=2) if phase=='Tx' else dict(coarse=3,fine=3,baseline=3)
        remaining_svd=dict(coarse=2,fine=2,baseline=2)
        tx=0 if self.s['components']['transactions']['status']=='PASS' else 9
        estimate=sum(remaining_by_tag[tag]*rates_by_tag[tag]+remaining_observations[tag]*max(observations[tag])*headroom+
                     remaining_svd[tag]*max(svd[tag])*headroom for tag in remaining_by_tag)
        estimate+=tx*rates_by_tag['coarse']+(60 if tx else 0)+600
        totals=self.ledger.totals();fits=totals['attempts']+remaining+tx<=800 and totals['seconds']+estimate<=10800
        value=dict(schema='WJ6MeasuredResourceProjectionV1',phase=phase,status='PASS' if fits else 'RESOURCE_LIMIT',
            parent_resource=file_identity(self.root/'resource-design.json'),native_seconds_per_step_bound=rate,
            native_seconds_per_step_bounds=rates_by_tag,remaining_by_tag=remaining_by_tag,
            measured_step_rates=rates,measured_observation_seconds=observations,measured_svd_seconds=svd,
            remaining_native_attempts=remaining,transaction_attempts_held=tx,closure_seconds_held=600,
            remaining_observation_allowance=remaining_observations,remaining_svd_allowance=remaining_svd,
            projected_remaining_seconds=estimate,headroom=headroom,work_totals=totals,
            actual_policy_bytes=volume(self.root),resource=DEFAULT,target=target)
        write_once(self.root/(phase+'-resource-projection.json'),value)
        self.s['components']['resource_projection' if phase=='Tx' else 'extension_resource']=component(value['status'],projection=file_identity(self.root/(phase+'-resource-projection.json')))
        self.save();return value

    def primary(self):
        require(meta(self.resource_file)['status']=='PASS','RESOURCE_LIMIT: frozen Tx resource gate did not pass')
        for tag,ends in (('coarse',(16,24,32)),('fine',(32,48,64)),('baseline',(16,32))):
            for end in ends:
                if int(self.s['trajectories'][tag]['trajectory'][-1]['steps'])<end:self.advance(tag,end)
        self.s['components']['primary']=self.qualify('primary');self.save()

    def qualify(self,phase,endcap=None):
        t=self.s['trajectories'];primary=phase=='primary';extension=phase=='extension'
        targets=dict(self.design['extension' if extension else 'primary' if primary else 'pilot'])
        if endcap is not None:targets.update(coarse=endcap,fine=2*endcap,baseline=endcap)
        parts={'fresh_startup':fresh_startup_gate(self.s.get('fresh_startup'))};observations={}
        if phase!='pilot':parts['pilot_gate']=self.s['components']['pilot']
        if extension:parts['primary_gate']=self.s['components']['primary']
        for tag,row in t.items():
            records={str(int(c['steps'])):c for c in [row['initial']]+row['trajectory']}
            expected=([0,4,16,24,32] if tag=='coarse' else [0,8,32,48,64] if tag=='fine' else [0,4,16,32]) if primary else [0,self.design['pilot'][tag]]
            if extension:expected=([0,4,16,24,32,64,96,128] if tag=='coarse' else [0,8,32,48,64,128,192,256] if tag=='fine' else [0,4,16,32,64,96,128])
            expected=[v for v in expected if v<=targets[tag]]
            parts[tag+'_accuracy']=required_accuracy(records,list(map(str,expected)))
            parts[tag+'_target']=component('PASS' if max(map(int,records))==targets[tag] else 'NOT_RUN')
            parts[tag+'_quadrature']=aggregate({k:c['quadrature'] for k,c in records.items()})
            parts[tag+'_arithmetic']=aggregate({str(step):component('PASS' if str(step) in records and records[str(step)].get('arithmetic') else 'NOT_RUN',
                'every declared arithmetic milestone required') for step in self.design['arithmetic'][tag] if step<=targets[tag]})
            for k,c in records.items():
                verify_artifacts({c['qualified_record']['path']:c['qualified_record']})
                require({name:value for name,value in c.items() if name!='qualified_record'}==meta(c['qualified_record']['path']),
                        'mutable progress record differs from immutable qualification packet')
                require(file_hash(c['path'])==c['sha256'] and file_hash(c['report'])==c['report_sha256'] and
                    file_hash(Path(c['checkpoint'])/'State.bin')==c['checkpoint_sha256'],'native observation/producer/checkpoint bytes changed')
                observations[tag,k]=read(c['path'])['measurements']
        if set(t)!={'coarse','fine','baseline'}:return aggregate(dict(parts,coverage=component('NOT_RUN')))
        for group,a,b in (('horizontal','baseline','coarse'),('timestep','coarse','fine')):
            times=[0,4,32] if primary else [0,4]
            if extension:times=[0,4,32,64,96,128]
            if endcap is not None:times=[v for v in times if v<=endcap]
            diffs={}
            for step in times:
                sa=step;sb=step*2 if b=='fine' else step
                if (a,str(sa)) not in observations or (b,str(sb)) not in observations:continue
                ma,mb=observations[a,str(sa)],observations[b,str(sb)]
                require(abs(ma['time']-mb['time'])<=2e-15,'paired physical times differ')
                diffs[str(step)]=comparison(ma,mb,self.a.policy,t[a]['initial']['rule_settings'])
            floors,boundary,coverage=paired_envelope(t[a],t[b])
            gate=sensitivity_sequence(diffs,floors,self.design,list(map(str,times)),scope_uncertainty={'constrained_boundary':boundary})
            parts[group+'_sensitivity']=gate
            self.s[group+'_comparisons']=diffs;self.s[group+'_uncertainty']=dict(floors=floors,boundary=boundary,coverage=coverage)
        self.s['comparisons_scope']='common physical composite observations; no unequal coefficient subtraction; one full-domain pressure gauge'
        return aggregate(parts)

    def transactions(self):
        require(meta(self.resource_file)['status']=='PASS','RESOURCE_LIMIT: transaction resource design failed')
        t=self.s['trajectories']['coarse'];c24=next(c for c in t['trajectory'] if int(c['steps'])==24);c32=next(c for c in t['trajectory'] if int(c['steps'])==32)
        whole=read(c32['report']);self.valid(whole,dict(exit_code=0))
        original=artifacts([p for p in Path(c24['checkpoint']).rglob('*') if p.is_file()])
        original.update(artifacts([Path(c24['report']),Path(c32['report'])]))
        copy_path=self.root/'transaction-private-CP24'
        if not copy_path.exists():shutil.copytree(c24['checkpoint'],copy_path)
        source={str(p.relative_to(c24['checkpoint'])):file_hash(p) for p in Path(c24['checkpoint']).rglob('*') if p.is_file()}
        require(source=={str(p.relative_to(copy_path)):file_hash(p) for p in copy_path.rglob('*') if p.is_file()},'private CP24 differs')
        frozen,r,folder=self.native('transaction-readback24',cp=copy_path,role='transaction_frozen_readback');self.valid(frozen,r)
        require(frozen['initialization_lift_calls_iterations']==[0,0],'frozen native readback called initialization Lift')
        rejected,r,rf=self.native('transaction-reject25',1,copy_path,reject=int(c24['steps'])+1,role='transaction_postsolve_reject')
        require(r['exit_code']!=0 and rejected and 'deliberate_PIMPLE_WJ_hard_audit_rejection' in rejected.get('failure',''),'actual-capacity end-of-step rejection reason differs')
        same_history(frozen,rejected)
        for key in ('state_audits','physical_samples','wj_audit_packets','accepted_states','audit_times','last_stage_history','last_stage_time','last_stage_alpha'):
            require(digest_object(frozen[key])==digest_object(rejected[key]),'rollback changed accepted callback/audit/sample packet: '+key)
        require(rejected['rejected_state'] and state_hash(rejected['rejected_state'])!=state_hash(rejected['state']),'failed candidate not separately retained')
        require(state_hash(rejected['rejected_state'])==state_hash(whole['accepted_states'][1]),'rejected step25 candidate differs from actual reference step25')
        attempt=rejected['pimple_attempts'][-1]
        require(rejected['implicit_callback_count']==1,'deliberate rejection must retain one actual callback')
        for key in ('history_sha256','incoming_sha256','time','alpha'):
            require(attempt[key]==whole['pimple_attempts'][0][key],'rejection incoming/callback context differs from reference')
        require(all(attempt[k]==1 for k in ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls')),'post-solve rejection lost scheduled counts')
        first,r,f=self.native('transaction-24-28',4,copy_path,role='transaction_resume');self.valid(first,r)
        second,r,f=self.native('transaction-28-32',4,first['checkpoint'],role='transaction_resume');self.valid(second,r)
        result=verify_transaction(whole,first,second,segment_length=8)
        joined=first['pimple_attempts']+second['pimple_attempts']
        for aa,bb in zip(whole['pimple_attempts'],joined):
            sa,sb=read(aa['k_snapshot']),read(bb['k_snapshot'])
            for key in ('committed','post_flow','pre_k','raw','corrected','history','fixed_reference','transport_residual','raw_increment','raw_linear_residual'):
                require(state_hash(sa[key])==state_hash(sb[key]),'retained transaction substage differs: '+key)
            for key in ('transport_rows','fixed_external_coefficients','reference_constraints'):
                require(digest_object(sa[key])==digest_object(sb[key]),'transaction functional/constraint differs: '+key)
            require(file_hash(sa['matrix_file'])==sa['matrix_sha256']==file_hash(sb['matrix_file'])==sb['matrix_sha256'],
                    'transaction frozen k matrix identity differs')
        verify_artifacts(original)
        require(source=={str(p.relative_to(copy_path)):file_hash(p) for p in copy_path.rglob('*') if p.is_file()},'transaction changed private input bytes')
        result.update(reference_producer=file_identity(c32['report']),reference_steps=[24,32],reference_accepted_indices=[0,8],
            separately_executed_reference=True,original_input_lineage=original,rollback_exact=True,postsolve_rejection_step=25,
            rejected_report=file_identity(rf/'case.evolution.json'),failed_candidate_sha256=state_hash(rejected['rejected_state']),
            failed_candidate_matches_reference_step25=True,all_retained_substages_and_frozen_k_matrices_equal=True,
            start_completed_stage_producer=file_identity(c24['report']),reference_start_scope='native restart-startup packet within actually executed24->32 report; completed CP24 stage uses its own producer')
        write_once(self.root/'transaction-result.json',result)
        self.s['components']['transactions']=component('PASS',result=file_identity(self.root/'transaction-result.json'));self.save()

    def extend(self):
        siblings={p:meta(self.a.output.resolve()/p/'summary.json')['components'] for p in ('EP','Budget')}
        require(common_extension_authorized(siblings),'T3 requires BOTH complete primary/transaction/resource gates')
        for end in (64,96,128):
            for tag,target in (('coarse',end),('fine',end*2),('baseline',end)):
                if int(self.s['trajectories'][tag]['trajectory'][-1]['steps'])<target:self.advance(tag,target)
                if any(c['errors']['accuracy_gate']['status']!='PASS' for c in self.s['trajectories'][tag]['trajectory']):
                    self.s['components']['extension']=component('OBSTRUCTED','accuracy crossing stops further extension; frozen T3 target remains');self.save();return
            gate=self.qualify('extension',endcap=end)
            self.s.setdefault('extension_cadence_gates',{})[str(end)]=gate;self.save()
            if gate['status']!='PASS':
                self.s['components']['extension']=component(gate['status'],'complete common-time cadence gate stopped extension',gate=gate);self.save();return
        self.s['components']['extension']=self.qualify('extension');self.save()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--solver',type=Path,required=True);p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--policy',choices=('EP','Budget'),required=True)
    p.add_argument('--mode',choices=('pilot','cost','primary','transactions','extension-cost','extend','report'),required=True);a=p.parse_args()
    a.solver=a.solver.resolve();a.probe=a.probe.resolve();s=Study(a)
    try:
        if a.mode=='cost':s.cost()
        elif a.mode=='extension-cost':s.cost('T3')
        elif a.mode!='report':getattr(s,a.mode)()
    except Exception as error:
        s.s['components'][a.mode]=component('RESOURCE_LIMIT' if 'RESOURCE_LIMIT' in str(error) else 'OBSTRUCTED',str(error))
        s.s.setdefault('failures',[]).append(dict(mode=a.mode,reason=str(error),code=s.code));s.save();raise
    finally:s.ledger.finish_session(s.session);s.save()
    print(s.summary)


if __name__=='__main__':main()
