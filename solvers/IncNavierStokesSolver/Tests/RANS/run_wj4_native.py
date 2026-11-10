#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded WJ-4 native finite-strip study. No missing component qualifies."""
import argparse
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import time
import uuid
import xml.etree.ElementTree as ET
from numpy.testing import assert_allclose
from wj2_cases import case_xml
from wj3_measure import checkpoint_payload,replay_witnesses,arithmetic,FIELDS
from wj4_measure import rule,boundary_rule,owned_rule,evaluate,comparison,sensitivity_gate,independent_weak,AMP,BUDGET
from wj_immutable import (Ledger,write_once,atomic,read,file_hash,artifacts,verify_artifacts,
                          runtime,execution_intent,digest,xml_identity)
from wj_qualification import component,aggregate,wj3_gates
from run_wj3_native import Study as WJ3Study
from run_cls_epawf_picard import THREADS,require,state_hash
from test_wj2_native import independent_loads,reference_comparison,digest_object
from compare_cls_epawf_pimple import same_history
from cls_epawf_m3_cases import set_info,write

SCRIPTS=Path(__file__).resolve().parent

class Study:
    def __init__(self,args):
        require(platform.system()=='Linux','compiled Linux is authoritative')
        started=time.time()
        self.a=args;self.root=args.output.resolve()/args.policy;self.root.mkdir(parents=True,exist_ok=True)
        self.summary=self.root/'summary.json';self.s=read(self.summary) if self.summary.exists() else dict(schema='WJ4StudyV1',policy=args.policy,execution={},components={})
        os.environ.update(THREADS)
        self.runtime={str(p):runtime(p) for p in (args.solver,args.probe)}
        native_code=artifacts([SCRIPTS/'wj2_cases.py',SCRIPTS/'cls_epawf_m4_cases.py',SCRIPTS/'cls_epawf_m3_cases.py',SCRIPTS/'EquilibriumWallFlat.xml'])
        self.code=dict(native_case_inputs=native_code,execution_engine=artifacts([SCRIPTS/'wj_immutable.py']))
        self.measure_code=artifacts([SCRIPTS/p for p in ('wj4_measure.py','wj3_measure.py','test_wj2_native.py','validate_wj1.py')])
        self.measure_version=digest(self.measure_code)[:12]
        design=dict(schema='WJ4FrozenIntentV1',policy=args.policy,runtime=self.runtime,governing_code=self.code,
            strip=[0.,1.,0.,2.],boundary='reference differential-field Dirichlet at both lateral boundaries and top',
            horizontal=dict(level=1,xknots=[[0.,1.],[0.,.5,1.]],T=.004096,dt=3.2e-5,steps=128),
            longer=dict(level=2,xknots=[0.,1.],T=.008192,dt=[3.2e-5,1.6e-5],steps=[256,512],cadence=[32,64],
                arithmetic_steps=[[32,128,256],[64,256,512]],window_endpoints=[[160,192,224,256],[320,384,448,512]]),
            amplitudes=AMP,budgets=BUDGET,sensitivity_fraction=.1,attempt_cap_per_policy=1500,wall_cap_seconds_per_policy=10800,
            source_order=16,basis_points=[8,8],observation_order_z=6,observation_order_x=4,independent_x_order=6,
            physics='WJHighRe Original DevStrainRawWorkV1; physical epsilon; uncapped calibrated wall; eta0; zero forces; one sweep U,p,k,epsilon; reuse/verify Off; StopOnBound On',
            observation_scope='common composite physical quadrature; explicit element/face traces; one full-domain pressure gauge')
        write_once(self.root/'design.json',design);self.design=design
        self.ledger=Ledger(Path(os.environ.get('WJ4_PACKAGE_ROOT',str(args.output.parent))),args.policy)
        self.session='driver-'+uuid.uuid4().hex
        self.ledger.begin_session(self.session,started)
        self.s.update(design_sha256=file_hash(self.root/'design.json'),runtime=self.runtime,measurement_code=self.measure_code,
            source_revision=subprocess.check_output(['git','-C',str(SCRIPTS.parents[3]),'rev-parse','HEAD'],text=True).strip())
        self.save()

    def save(self):
        self.s['ledger_totals']=self.ledger.totals() if hasattr(self,'ledger') else {}
        atomic(self.summary,self.s)

    def execute(self,folder,argv,inputs,settings,steps=0,outputs=(),measurement=False):
        argv=list(map(str,argv));folder.mkdir(parents=True,exist_ok=True)
        intent=execution_intent(argv,inputs,self.measure_code if measurement else self.code,self.runtime,settings)
        write_once(folder/'execution-intent.json',intent);receipt=folder/'receipt.json'
        if receipt.exists():
            for value in self.runtime.values():
                require(file_hash(value['path'])==value['sha256'],'runtime binary changed before reuse')
                for path,h in value['libraries'].items():require(file_hash(path)==h,'runtime dependency changed before reuse')
            r=read(receipt);require(r['intent_sha256']==digest(intent),'cached execution fingerprint mismatch')
            verify_artifacts(r['artifacts']);return r
        require(not (folder/'command.log').exists(),'unreconciled attempt retained; inspect ledger and choose a new attempt name')
        label=str(folder.relative_to(self.a.output));limit=max(120.,steps*16.+120.)
        self.ledger.reserve(label,intent,steps,limit)
        start=time.time();code=-1;failure=None
        try:
            with (folder/'command.log').open('x') as stream:
                code=subprocess.run(argv,cwd=folder,stdout=stream,stderr=subprocess.STDOUT,env=os.environ.copy(),timeout=limit).returncode
        except subprocess.TimeoutExpired:code=124;failure='RESOURCE_LIMIT: command timeout'
        except Exception as error:failure=str(error)
        report=folder/'case.evolution.json';d=read(report) if report.exists() else None
        actual=len(d['pimple_attempts']) if d else 0 if steps==0 else None
        files=[folder/'command.log',*outputs]
        if d and d.get('checkpoint'):files.append(Path(d['checkpoint'])/'State.bin')
        r=dict(argv=argv,intent_sha256=digest(intent),started_epoch=start,seconds=time.time()-start,exit_code=code,
            failure=failure,requested_steps=steps,actual_attempts=actual,charged_attempts=actual if actual is not None else steps,
            attempt_scope='native report' if d else 'unknown; reservation remains charged' if steps else 'no scheduled steps',artifacts=artifacts(files),threads=THREADS)
        write_once(receipt,r);self.ledger.finish(label,actual,r['seconds'],str(receipt));self.save();return r

    def native(self,name,steps=0,cp=None,dt=3.2e-5,level=1,nx=1,reject=0,arithmetic_time=None):
        folder=self.root/name;xml=folder/'case.xml';xknots=[i/nx for i in range(nx+1)]
        options=dict(steps=steps,restart=str(Path(cp).resolve()) if cp else None,dt=dt,level=level,xknots=xknots,
            policy=self.a.policy,reject=reject,arithmetic_time=arithmetic_time,stop_on_bound=True)
        if xml.exists():
            require((folder/'requested-intent.json').exists(),'existing XML has no immutable intent')
            write_once(folder/'requested-intent.json',options)
            require(read(folder/'xml-identity.json')==xml_identity(xml),'stale XML bytes/semantics')
        else:
            case_xml(folder,self.a.policy,steps,cp,dt=dt,level=level,xknots=xknots,hard_reject=reject)
            tree=ET.parse(xml).getroot();set_info(tree,'RANSWallPimpleStopOnBound','On')
            set_info(tree,'RANSWallPimpleSnapshotDirectory',str(folder/'pimple-states'))
            if arithmetic_time is not None:
                set_info(tree,'RANSWallPimpleArithmetic','On');set_info(tree,'RANSWallPimpleArithmeticTime',format(arithmetic_time,'.17g'))
            write(tree,xml);write_once(folder/'requested-intent.json',options);write_once(folder/'xml-identity.json',xml_identity(xml))
        report=folder/'case.evolution.json';inputs=[xml]+([Path(cp)/'State.bin'] if cp else [])
        r=self.execute(folder,[self.a.solver,xml],inputs,options,steps,[report])
        d=read(report) if report.exists() else None
        self.s['execution'][name]=dict(status='PASS' if r['exit_code']==0 and d and not d.get('failure') else 'OBSTRUCTED',
            receipt=str(folder/'receipt.json'),receipt_sha256=file_hash(folder/'receipt.json'),report=str(report),
            report_sha256=file_hash(report) if report.exists() else None,attempts=r['actual_attempts'])
        self.save();return d,r,folder

    def valid(self,d,r):WJ3Study.valid(self,d,dict(exit_code=r['exit_code']))

    def measure(self,name,d,producer,level=1,nx=1,assembly=False,xorder=4):
        # Each measurement version is write-once and explicitly bound to its
        # producing report. A zero-step restart report is not a stage producer.
        folder=self.root/('measure-'+self.measure_version)/name;folder.mkdir(parents=True,exist_ok=True)
        cp=Path(d['checkpoint']);report=producer/'case.evolution.json';xml=producer/'case.xml'
        before=artifacts([cp/'State.bin',report])
        physical=rule(xorder=xorder);owned=owned_rule(physical,level,[i/nx for i in range(nx+1)])
        packet=dict(operation='wj_measure',inner_stress_model='Linear',native_checkpoint=str(cp),
            checkpoint_sha256=file_hash(cp/'State.bin'),native_signature=d['signature'],state=d['state'],history=d['state'],time=d['time'],alpha=0.,
            native_differential_state=d['state'][:int(d['offsets'][4])],native_history=d['native_history'],native_history_times=d['native_history_times'],
            native_report=str(report),native_report_sha256=file_hash(report),observations=owned,export_assembly_checks=assembly)
        write_once(folder/'input.json',packet)
        output=folder/'result.json'
        r=self.execute(folder,[self.a.probe,xml,folder/'input.json',output],[xml,folder/'input.json',cp/'State.bin',report],
            dict(operation='wj_measure',assembly=assembly,observation_sha256=digest(physical),xorder=xorder),outputs=[output],measurement=True)
        require(r['exit_code']==0,'strict native measurement failed: '+str(folder/'command.log'));verify_artifacts(before)
        raw=read(output);m=raw['measurements']
        require(m['state_sha256']==m['state_after_sha256'] and m['stage_history_sha256']==m['stage_history_after_sha256'] and m['governing_solves']==0,'read-only mutation or solve')
        c=dict(path=str(output),sha256=file_hash(output),report_sha256=file_hash(report),checkpoint_sha256=file_hash(cp/'State.bin'),
            observation_sha256=digest(physical),boundary_quadrature_sha256=digest([[x,z,w] for (x,z),w in boundary_rule(xorder=xorder).items()]),time=m['time'],steps=d['steps'],state_sha256=m['state_sha256'],history_sha256=m['stage_history_sha256'],
            spatial_residual=m['spatial_mass_dual_norms'],stage_defect=m['completed_stage']['mass_dual_norms'] if m['completed_stage'] else None,
            errors=evaluate(m,self.a.policy,xorder),witness_replay=replay_witnesses(m['audit']))
        if assembly:
            checks=raw['assembly_checks'];checks['policy']=self.a.policy
            c['independent_loads']=independent_loads(checks,self.a.policy,d['offsets'])
            c['independent_weak_modes']=independent_weak(checks,m['weak_modes'])
            if d['steps']==0:c['continuum_hypotheses']=reference_comparison(checks,self.a.policy)['continuum']
        write_once(folder/'comparison.json',c)
        return m,c

    def inventory(self):
        rows=[]
        for level,nx,count in ((0,1,639),(0,2,1109),(0,4,2049),(1,1,1257),(1,2,2183),(2,1,2493),(2,2,4331)):
            name='mesh-L%d-Nx%d-cp0'%(level,nx)
            d,r,f=self.native(name,level=level,nx=nx);self.valid(d,r)
            require(len(d['state'])==count,'unexpected native startup mesh size')
            rows.append(dict(level=level,nx=nx,unknowns=len(d['state']),steps=d['steps'],report=str(f/'case.evolution.json'),
                evolution_supported=count<=3000,native_report_sha256=file_hash(f/'case.evolution.json')))
            self.s['mesh_inventory']=rows;self.save()
        self.s['components']['mesh_inventory']=component('PASS',meshes=rows);self.save()

    def horizontal(self):
        records={};measurements={}
        for nx in (1,2):
            z,r,f=self.native('horizontal-Nx%d-cp0'%nx,nx=nx);self.valid(z,r)
            require(len(z['state'])=={1:1257,2:2183}[nx],'unexpected native horizontal mesh size')
            m,c=self.measure('horizontal-Nx%d-cp0'%nx,z,f,nx=nx,assembly=True)
            records[str(nx)]=dict(initial=c,unknowns=len(z['state']),trajectory=[]);measurements[(nx,0)]=m
            cp=z['checkpoint'];last=z
            for end in (32,64,96,128):
                d,r,producer=self.native('horizontal-Nx%d-%d'%(nx,end),32,cp,nx=nx,arithmetic_time=end*3.2e-5 if end in (32,128) else None)
                try:self.valid(d,r)
                except ValueError as error:
                    records[str(nx)]['hard_obstruction']=str(error);break
                m,c=self.measure('horizontal-Nx%d-%d'%(nx,end),d,producer,nx=nx,assembly=end==128)
                if end in (32,128):
                    c['arithmetic']=self.arithmetic('horizontal-Nx%d-%d'%(nx,end),d,m)
                records[str(nx)]['trajectory'].append(c);measurements[(nx,end)]=m;cp=d['checkpoint'];last=d
                self.s['horizontal']=records;self.save()
                if c['errors']['accuracy_gate']['status']!='PASS':
                    records[str(nx)]['accuracy_obstruction']=dict(frozen_target=.004096,accepted_time=d['time'],gate=c['errors']['accuracy_gate'])
                    if end!=128:
                        _,loads=self.measure('horizontal-Nx%d-%d-obstruction-loads'%(nx,end),d,producer,nx=nx,assembly=True)
                        records[str(nx)]['independent_obstruction_loads']=loads
                    break
            records[str(nx)]['last_result']=str(producer/'case.evolution.json');self.s['horizontal']=records;self.save()
        common=sorted(set(t for n,t in measurements if n==1)&set(t for n,t in measurements if n==2))
        comparisons={str(t):comparison(measurements[(1,t)],measurements[(2,t)],self.a.policy) for t in common}
        floors={f:dict(L2=0.,sampled_Linf=0.) for f in FIELDS+('velocity',)}
        for row in records.values():
            for c in row['trajectory']:
                for f,v in c.get('arithmetic',{}).get('estimated_comparison_floor',{}).items():
                    floors[f]['L2']=max(floors[f]['L2'],v)
                    floors[f]['sampled_Linf']=max(floors[f]['sampled_Linf'],c['arithmetic']['estimated_sampled_comparison_floor'][f])
        gate=sensitivity_gate(comparisons[str(common[-1])],floors)
        required=all(row['trajectory'] and row['trajectory'][-1]['steps']==128 and not row.get('accuracy_obstruction') and not row.get('hard_obstruction') for row in records.values())
        self.s['horizontal_comparisons']=comparisons
        self.s['components']['horizontal']=aggregate(dict(target=component('PASS' if required else 'OBSTRUCTED','both frozen T3 targets and accuracy gates required'),sensitivity=gate))
        # Independent horizontal integration increase on both startup fields.
        q=[]
        for nx in (1,2):
            d=read(self.root/('horizontal-Nx%d-cp0'%nx)/'case.evolution.json');m,c=self.measure('horizontal-Nx%d-q6'%nx,d,self.root/('horizontal-Nx%d-cp0'%nx),nx=nx,xorder=6)
            base=records[str(nx)]['initial']['errors']['continuum_error']['full']['fields'];new=c['errors']['continuum_error']['full']['fields']
            delta={f:abs(base[f]['L2']-new[f]['L2']) for f in base};require(max(delta.values())<1e-10,'horizontal quadrature increase is material')
            q.append(dict(nx=nx,order=[4,6],L2_delta=delta,absolute_floor=1e-10))
        self.s['horizontal_quadrature']=q;self.save()

    def arithmetic(self,name,d,m):
        folder=self.root/('measure-'+self.measure_version)/name;path=folder/'arithmetic.json'
        if path.exists():return read(path)
        intent=dict(report_state=state_hash(d['state']),measurement_state=m['state_sha256'],code=self.measure_code)
        label=str(folder.relative_to(self.a.output))+'/svd';self.ledger.reserve(label,intent,0,300)
        start=time.time()
        try:value=arithmetic(d,m)
        except Exception:
            self.ledger.finish(label,0,time.time()-start,'failed offline arithmetic');raise
        value['scope']='finite frozen-system estimate at declared early/middle/end checkpoints; envelopes are not a nonlinear perturbation theorem'
        write_once(path,value);self.ledger.finish(label,0,time.time()-start,str(path));self.save();return value

    def cost(self):
        if (self.root/'cost-design.json').exists():
            self.s['cost_design']=read(self.root/'cost-design.json');self.save();return
        d,r,f=self.native('cost-fine-cp0',level=2);self.valid(d,r)
        s,rr,sf=self.native('cost-fine-4',4,d['checkpoint'],level=2,arithmetic_time=4*3.2e-5);self.valid(s,rr)
        m,c=self.measure('cost-fine-4',s,sf,level=2)
        arithmetic_packet=self.arithmetic('cost-fine-4',s,m)
        seconds=rr['seconds']/4
        measurement=read(Path(c['path']).parent/'receipt.json')['seconds']
        # Includes a process/startup allowance on every continuation block,
        # all32 remaining observation packets and six SVD envelope inspections.
        estimate=768*seconds+32*(r['seconds']+measurement)+6*60+600
        cost=dict(measured_four_step_seconds=rr['seconds'],measured_seconds_per_step=seconds,
            measured_startup_seconds=r['seconds'],measured_observation_seconds=measurement,
            conservative_remaining_pair_seconds=estimate,used_seconds=self.ledger.totals()['seconds'],
            practical_package_wall_cap_seconds=10800,attempt_cap=1500,frozen_before_longer=True,
            longer_T=.008192,longer_dt=[3.2e-5,1.6e-5],longer_steps=[256,512])
        write_once(self.root/'cost-design.json',cost);self.s['cost_design']=cost
        self.s['components']['resource_design']=component('PASS' if cost['used_seconds']+estimate<=10800 else 'RESOURCE_LIMIT',**cost);self.save()

    def longer(self):
        require((self.root/'cost-design.json').exists(),'cost design must be frozen before longer runs')
        cost=read(self.root/'cost-design.json')
        if cost['used_seconds']+cost['conservative_remaining_pair_seconds']>10800:
            self.s['components']['longer']=component('RESOURCE_LIMIT','predeclared full target does not fit measured budget');self.save();return
        trajectories={};measurements={};initial=None
        for j,(dt,n,block) in enumerate(((3.2e-5,256,32),(1.6e-5,512,64))):
            tag=('coarse','fine')[j];z,r,f=self.native('longer-'+tag+'-cp0',dt=dt,level=2);self.valid(z,r)
            if initial is not None:
                require(state_hash(z['state'])==state_hash(initial['state']) and digest_object(z['native_history'])==digest_object(initial['native_history']),'fresh dt startup physical fields/gauge/history differ')
            initial=z
            m,c=self.measure('longer-'+tag+'-cp0',z,f,level=2);measurements[(j,0)]=m
            row=dict(dt=dt,target_steps=n,target_time=.008192,initial=c,initial_primary_sha256=state_hash(z['state']),
                initial_history_sha256=digest_object(z['native_history']),native_time_labels=z['native_history_times'],trajectory=[])
            cp=z['checkpoint']
            for end in range(block,n+1,block):
                arithmetic_time=end*dt if end in (block,n//2,n) else None
                d,r,producer=self.native('longer-'+tag+'-'+str(end),block,cp,dt=dt,level=2,arithmetic_time=arithmetic_time)
                try:self.valid(d,r)
                except ValueError as error:row['hard_obstruction']=dict(reason=str(error),producer=str(producer),accepted_time=d['time'] if d else None);break
                m,c=self.measure('longer-'+tag+'-'+str(end),d,producer,level=2,assembly=end==n)
                if arithmetic_time is not None:c['arithmetic']=self.arithmetic('longer-'+tag+'-'+str(end),d,m)
                row['trajectory'].append(c);measurements[(j,int(end/block))]=m;cp=d['checkpoint']
                trajectories[tag]=row;self.s['longer']=trajectories;self.save()
                if c['errors']['accuracy_gate']['status']!='PASS':
                    row['accuracy_obstruction']=c['errors']['accuracy_gate']
                    if end!=n:
                        _,loads=self.measure('longer-'+tag+'-'+str(end)+'-obstruction-loads',d,producer,level=2,assembly=True)
                        row['independent_obstruction_loads']=loads
                    break
            row['last_producer']=str(producer);trajectories[tag]=row;self.s['longer']=trajectories;self.save()
        targets=all(v['trajectory'] and v['trajectory'][-1]['steps']==v['target_steps'] and not v.get('hard_obstruction') and not v.get('accuracy_obstruction') for v in trajectories.values())
        floors={f:dict(L2=0.,sampled_Linf=0.) for f in FIELDS+('velocity',)}
        for row in trajectories.values():
            for c in row['trajectory']:
                for f,v in c.get('arithmetic',{}).get('estimated_comparison_floor',{}).items():
                    floors[f]['L2']=max(floors[f]['L2'],v)
                    floors[f]['sampled_Linf']=max(floors[f]['sampled_Linf'],c['arithmetic']['estimated_sampled_comparison_floor'][f])
        pairs={}
        for window in sorted(set(t for j,t in measurements if j==0)&set(t for j,t in measurements if j==1)):
            d=comparison(measurements[(0,window)],measurements[(1,window)],self.a.policy)
            pairs[str(window)]=dict(time=measurements[(0,window)]['time'],fine_time=measurements[(1,window)]['time'],comparison=d,
                gate=sensitivity_gate(d,floors))
        self.s['timestep_comparisons']=pairs;self.s['arithmetic_envelope']=dict(estimated_comparison_floor=floors,scope='maximum early/middle/end finite frozen estimates; finite observation basis-row norms and native full-mass eigenvalues convert scalar/velocity uncertainties; pressure uses its coefficient estimate; not a nonlinear theorem')
        self.s['components']['longer']=aggregate(dict(target=component('PASS' if targets else 'OBSTRUCTED','full frozen T4 pair with inactive correction and accuracy required'),
            sensitivity=pairs[str(max(map(int,pairs)))]['gate'] if pairs else component('NOT_RUN')))
        self.stationarity(trajectories,measurements)
        self.transactions(trajectories)
        self.save()

    def stationarity(self,trajectories,measurements):
        rows=[]
        for j,tag in enumerate(('coarse','fine')):
            times=(160,192,224,256) if j==0 else (320,384,448,512);block=32 if j==0 else 64
            if not all((j,n//block) in measurements for n in times):continue
            for start,end in zip(times,times[1:]):
                a,b=measurements[(j,start//block)],measurements[(j,end//block)]
                diff=comparison(a,b,self.a.policy)['full']['fields'];dt=b['time']-a['time']
                V={f:2*diff[f]['normalized_L2']/dt for f in ('velocity','k','outer_epsilon')}
                residual=b['spatial_mass_dual_norms'];start_residual=a['spatial_mass_dual_norms']
                passed=max(V.values())<=1e-5 and max(abs(r[i]) for r in (start_residual,residual) for i in (0,1,2,3))<=1e-5 and max(abs(v) for r in (start_residual,residual) for v in r[4:])<=1e-10
                rows.append(dict(trajectory=tag,steps=[start,end],times=[a['time'],b['time']],V=V,spatial_residual=residual,start_spatial_residual=start_residual,status='PASS' if passed else 'NOT_REACHED'))
        self.s['stationarity']=dict(status='REACHED' if len(rows)==6 and all(r['status']=='PASS' for r in rows) else 'NOT_REACHED',windows=rows)

    def transactions(self,trajectories):
        coarse=trajectories.get('coarse',{});rows=coarse.get('trajectory',[])
        if len(rows)<2:
            self.s['components']['transactions']=component('NOT_RUN','no retained late32-step segment');return
        end=int(rows[-1]['steps']);start=end-32
        original=self.root/('longer-coarse-'+str(end));whole=read(original/'case.evolution.json')
        predecessor=self.root/('longer-coarse-'+str(start));first=read(predecessor/'case.evolution.json');cp=Path(first['checkpoint'])
        immutable=artifacts([cp/'State.bin'])
        a,rc,af=self.native('transaction-first16',16,cp,level=2);self.valid(a,rc)
        b,rc,bf=self.native('transaction-second16',16,a['checkpoint'],level=2);self.valid(b,rc)
        require(same_history(whole,b),'late16+16 primary/gauge/time/history mismatch')
        for key in ('state_audits','wj_audit_packets'):
            require(digest_object(whole[key][-1])==digest_object(b[key][-1]),'late accepted audit mismatch: '+key)
        require(digest_object(whole['physical_samples'])==digest_object(b['physical_samples']),'late physical samples mismatch')
        self.measure('transaction-actual-producer-readback',first,predecessor,level=2)
        frozen,rc,_=self.native('transaction-zero-readback',0,cp,level=2);self.valid(frozen,rc);require(same_history(first,frozen),'zero readback changed state/history/time')
        rejected,rc,rf=self.native('transaction-reject',1,cp,level=2,reject=start+1)
        require(rejected and 'deliberate_PIMPLE_WJ_hard_audit_rejection' in rejected.get('failure',''),'wrong deliberate rejection reason')
        require(same_history(first,rejected),'rejection changed accepted state/history/time')
        for key in ('state_audits','wj_audit_packets'):
            require(digest_object(frozen[key])==digest_object(rejected[key]),'rollback audit packet mismatch')
        require(digest_object(first['physical_samples'])==digest_object(rejected['physical_samples']),'rollback sample mismatch')
        require(state_hash(rejected['rejected_state'])!=state_hash(rejected['state']),'failed candidate not separate')
        for key in ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls'):
            require(rejected['pimple_attempts'][-1][key]==1,'lost rejected scheduled solve count')
        verify_artifacts(immutable)
        self.s['components']['transactions']=component('PASS',segment=[start,end],exact_primary_gauge_history_time=True,exact_accepted_audit_samples=True,
            immutable_checkpoint=immutable,producer=str(predecessor/'case.evolution.json'),failed_candidate_sha256=state_hash(rejected['rejected_state']),rejection_report=str(rf/'case.evolution.json'))

    def finish(self):
        required={k:self.s['components'].get(k,component('NOT_RUN','required component missing')) for k in ('mesh_inventory','horizontal','resource_design','longer','transactions')}
        self.s['qualification']=aggregate({**self.s['components'],**required})
        self.s['longer_interpretation']='CONDITIONAL finite reference-Dirichlet strip evidence' if required['horizontal']['status']!='PASS' else 'finite reference-Dirichlet strip evidence'
        self.save()

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--solver',type=Path,required=True);p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--policy',choices=('EP','Budget'),required=True)
    p.add_argument('--mode',choices=('inventory','horizontal','cost','longer','all','report'),required=True);a=p.parse_args()
    a.output=a.output.resolve();a.solver=a.solver.resolve();a.probe=a.probe.resolve();s=Study(a)
    try:
        for mode in (('inventory','horizontal','cost','longer') if a.mode=='all' else (a.mode,)):
            if mode!='report':getattr(s,mode)()
        if s.s.get('failure') and (s.s.get('failure_mode')==a.mode or a.mode=='all'):
            s.s.setdefault('resolved_harness_failures',[]).append(dict(reason=s.s.pop('failure'),mode=s.s.pop('failure_mode',None),
                resolution_source=s.s['source_revision'],measurement_code=s.measure_code))
            if a.mode in ('all','cost'):s.s['components'].pop(a.mode,None)
        s.finish()
    except Exception as error:
        s.s['components'][a.mode]=component('RESOURCE_LIMIT' if 'RESOURCE_LIMIT' in str(error) else 'FAIL',str(error));s.s['failure']=str(error);s.s['failure_mode']=a.mode;s.finish();raise
    finally:
        s.ledger.finish_session(s.session);s.save()
    print(s.summary)

if __name__=='__main__':main()
