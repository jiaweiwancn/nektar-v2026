#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded WJ-5 local normal refinement, retained-state diagnosis and checks."""
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
import numpy as np
from run_wj4_native import Study as WJ4Study
from wj2_cases import case_xml,local_normal_knots
from wj3_measure import checkpoint_payload,replay_witnesses,FIELDS
from wj4_measure import rule,boundary_rule,owned_rule,evaluate,comparison,independent_weak,AMP,BUDGET
from wj4_qualification import sensitivity_sequence,arithmetic_floors
from wj5_measure import diagnose,boundary_uncertainty
from wj_immutable import (Ledger,write_once,atomic,read,file_hash,artifacts,verify_artifacts,
                          runtime,execution_intent,digest,xml_identity)
from wj_qualification import component,aggregate
from run_cls_epawf_picard import THREADS,require,state_hash
from test_wj2_native import independent_loads,reference_comparison,digest_object
from compare_cls_epawf_pimple import same_history
from verify_wj4_transactions import verify as verify_transaction
from cls_epawf_m3_cases import set_info,write

SCRIPTS=Path(__file__).resolve().parent

class Study(WJ4Study):
    def __init__(self,a):
        require(platform.system()=='Linux','compiled Linux is authoritative')
        start=time.time();self.a=a;self.read_only_reprocess=False;self.native_parent=None
        os.environ.update(THREADS);os.environ['WJ_PACKAGE_ATTEMPT_CAP']='800'
        self.root=a.output.resolve()/a.policy;self.root.mkdir(parents=True,exist_ok=True)
        self.summary=self.root/'summary.json';self.s=read(self.summary) if self.summary.exists() else dict(schema='WJ5StudyV1',policy=a.policy,execution={},components={})
        self.runtime={str(p):runtime(p) for p in (a.solver,a.probe)}
        self.code=dict(native_case_inputs=artifacts([SCRIPTS/p for p in ('wj2_cases.py','cls_epawf_m4_cases.py','cls_epawf_m3_cases.py','EquilibriumWallFlat.xml')]),
            execution_engine=artifacts([SCRIPTS/'wj_immutable.py']))
        self.measure_code=artifacts([SCRIPTS/p for p in ('run_wj5_native.py','wj5_measure.py','wj4_measure.py','wj3_measure.py','wj4_qualification.py','wj_qualification.py','test_wj2_native.py','validate_wj1.py')]);self.measure_version=digest(self.measure_code)[:12]
        self.design=dict(schema='WJ5FrozenDesignV1',policy=a.policy,runtime=self.runtime,governing_code=self.code,
            retained_root=str(a.retained.resolve()),strip=[0.,1.,0.,2.],boundary='reference Dirichlet differential traces at both laterals/top; no-slip ground; zero modeled ground k flux',
            meshes={m:local_normal_knots(m) for m in ('L1','M1','M2')},xknots=[[0.,1.],[0.,.5,1.]],
            control=dict(mesh='L1',nx=2,dt=1.6e-5,steps=64,T=.001024,adverse_diagnostic=True,blocks=1),
            screen=dict(dt=3.2e-5,steps=32,T=.001024,escalation='M2 only if M1 fails/unresolved in either policy'),
            selected=dict(T=.004096,coarse_dt=3.2e-5,coarse_steps=128,fine_dt=1.6e-5,fine_steps=256,cadence=[32,64],
                arithmetic_steps=[[32,64,128],[64,128,256]],transactions='reuse late32 + extra16+16 + frozen readback + one post-solve rejection'),
            amplitudes=AMP,budgets=BUDGET,sensitivity_fraction=.1,attempt_cap_per_policy=800,wall_cap_seconds_per_policy=10800,
            observations=dict(order=6,xorder=4,independent_orders=[8,6],all_faces=True),source_order=16,basis_points=[8,8],
            physics='unchanged WJHighRe Original DevStrainRawWorkV1, direct epsilon, uncapped calibrated AWF, independent Cd, one U/p/k/epsilon sweep, eta0, zero forces, reuse/verify Off, WJNativeV1 MassProjection StopOnBound On')
        write_once(self.root/'design.json',self.design)
        self.s['design_sha256']=file_hash(self.root/'design.json');self.s['measurement_code']=self.measure_code
        self.s['source_revision']=subprocess.check_output(['git','-C',str(SCRIPTS.parents[3]),'rev-parse','HEAD'],text=True).strip()
        for k in ('retained_diagnosis','adverse_dt_control','local_screen','resource_design','selected_horizontal','selected_timestep','transactions'):
            self.s['components'].setdefault(k,component('NOT_RUN','required component not executed'))
        self.ledger=Ledger(Path(os.environ.get('WJ4_PACKAGE_ROOT',str(a.output.parent))),a.policy,cap=800)
        self.session='wj5-'+uuid.uuid4().hex;self.ledger.begin_session(self.session,start);self.save()

    def native(self,name,steps=0,cp=None,dt=3.2e-5,mesh='M1',nx=2,reject=0,arithmetic_time=None,role='local_screen',normal=16):
        folder=self.root/name;xml=folder/'case.xml';x=[i/nx for i in range(nx+1)];z=self.design['meshes'][mesh]
        options=dict(steps=steps,restart=str(Path(cp).resolve()) if cp else None,dt=dt,level=1,xknots=x,zknots=z,comparison_role=role,
            policy=self.a.policy,reject=reject,arithmetic_time=arithmetic_time,source_order=normal,stop_on_bound=True)
        if xml.exists():
            require((folder/'requested-intent.json').exists(),'existing XML has no immutable intent')
            write_once(folder/'requested-intent.json',options);require(read(folder/'xml-identity.json')==xml_identity(xml),'stale XML bytes/semantics')
        else:
            case_xml(folder,self.a.policy,steps,cp,dt=dt,level=1,xknots=x,zknots=z,hard_reject=reject,normal=normal)
            tree=ET.parse(xml).getroot();set_info(tree,'RANSWallPimpleStopOnBound','On');set_info(tree,'RANSWallPimpleSnapshotDirectory',str(folder/'pimple-states'))
            if arithmetic_time is not None:set_info(tree,'RANSWallPimpleArithmetic','On');set_info(tree,'RANSWallPimpleArithmeticTime',format(arithmetic_time,'.17g'))
            write(tree,xml);write_once(folder/'requested-intent.json',options);write_once(folder/'xml-identity.json',xml_identity(xml))
        report=folder/'case.evolution.json';inputs=[xml]+([Path(cp)/'State.bin'] if cp else [])
        r=self.execute(folder,[self.a.solver,xml],inputs,options,steps,[report]);d=read(report) if report.exists() else None
        self.s['execution'][name]=dict(status='PASS' if r['exit_code']==0 and d and not d.get('failure') else 'OBSTRUCTED',receipt=str(folder/'receipt.json'),
            receipt_sha256=file_hash(folder/'receipt.json'),report=str(report),report_sha256=file_hash(report) if report.exists() else None,attempts=r['actual_attempts'],role=role)
        self.save();return d,r,folder

    def retained(self,name,nx,step):
        folder=self.a.retained.resolve()/'study'/self.a.policy/('horizontal-Nx%d-%s'%(nx,'cp0' if step==0 else step))
        report=folder/'case.evolution.json';intent=read(folder/'execution-intent.json');receipt=read(folder/'receipt.json')
        require(receipt['intent_sha256']==digest(intent),'retained native intent changed');verify_artifacts(intent['inputs']);verify_artifacts(receipt['artifacts'])
        for group in intent['code'].values():verify_artifacts(group)
        for identity in intent['runtime'].values():
            require(file_hash(identity['path'])==identity['sha256'],'retained producer runtime changed')
            for p,h in identity['libraries'].items():require(file_hash(p)==h,'retained producer dependency changed')
        d=read(report);require(d['steps']==step,'retained step differs')
        m,c=self.measure(name,d,folder,mesh='L1',nx=nx,assembly=True)
        c['producer_lineage']=dict(inputs=artifacts([report,folder/'case.xml',folder/'receipt.json',folder/'execution-intent.json']),runtime=intent['runtime'])
        return d,m,c

    def measure(self,name,d,producer,mesh='M1',nx=2,assembly=False,order=6,xorder=4,diagnostics=True):
        folder=self.root/('measure-'+self.measure_version)/name;folder.mkdir(parents=True,exist_ok=True)
        cp=Path(d['checkpoint']);report=producer/'case.evolution.json';xml=producer/'case.xml';before=artifacts([cp/'State.bin',report])
        settings=dict(xmeshes=self.design['xknots'],zmeshes=[self.design['meshes'][mesh]],order=order,xorder=xorder,all_faces=True)
        physical=rule(**settings);owned=owned_rule(physical,1,[i/nx for i in range(nx+1)],self.design['meshes'][mesh])
        packet=dict(operation='wj_measure',inner_stress_model='Linear',native_checkpoint=str(cp),checkpoint_sha256=file_hash(cp/'State.bin'),
            native_signature=d['signature'],state=d['state'],history=d['state'],time=d['time'],alpha=0.,
            native_differential_state=d['state'][:int(d['offsets'][4])],native_history=d['native_history'],native_history_times=d['native_history_times'],
            native_report=str(report),native_report_sha256=file_hash(report),observations=owned,export_assembly_checks=assembly,export_wj5_diagnostics=diagnostics)
        write_once(folder/'input.json',packet);output=folder/'result.json'
        r=self.execute(folder,[self.a.probe,xml,folder/'input.json',output],[xml,folder/'input.json',cp/'State.bin',report],
            dict(operation='wj_measure',assembly=assembly,diagnostics=diagnostics,role='passive_readback',mesh=mesh,nx=nx,zknots=self.design['meshes'][mesh],
                observation_sha256=digest(physical),order=order,xorder=xorder),outputs=[output],measurement=True)
        require(r['exit_code']==0,'strict native measurement failed: '+str(folder/'command.log'));verify_artifacts(before)
        raw=read(output);m=raw['measurements'];require(m['state_sha256']==m['state_after_sha256'] and m['stage_history_sha256']==m['stage_history_after_sha256'] and m['governing_solves']==0,'read-only mutation or solve')
        errors=evaluate(m,self.a.policy,xorder,settings)
        # Every field and scope is required, including prescribed-boundary data.
        leaves={}
        for kind in ('continuum_error','initial_drift'):
            for scope in ('full','central','constrained_boundary'):
                for field,row in errors[kind][scope]['fields'].items():
                    for norm in ('L2','sampled_Linf'):
                        value=row['normalized_'+norm];limit=BUDGET[norm]
                        leaves[kind+'/'+scope+'/'+field+'/'+norm]=component('PASS' if math.isfinite(value) and value<=limit else 'OBSTRUCTED',value=value,limit=limit)
        errors['accuracy_gate']=aggregate(dict(original=errors['accuracy_gate'],all_scopes=aggregate(leaves)))
        c=dict(path=str(output),sha256=file_hash(output),report=str(report),report_sha256=file_hash(report),checkpoint=str(cp),checkpoint_sha256=file_hash(cp/'State.bin'),
            mesh=mesh,nx=nx,unknowns=len(d['state']),zknots=self.design['meshes'][mesh],rule_settings=settings,observation_sha256=digest(physical),
            time=m['time'],steps=d['steps'],state_sha256=m['state_sha256'],history_sha256=m['stage_history_sha256'],
            spatial_residual=m['spatial_mass_dual_norms'],stage_defect=m['completed_stage']['mass_dual_norms'] if m['completed_stage'] else None,
            errors=errors,gradient_replay=dict(witness=m['strong_divergence_witness'],max_defect=m['fast_independent_gradient_max_defect'],observations=m['independent_gradient_observations']),
            witness_replay=replay_witnesses(m['audit']))
        if diagnostics:c['diagnosis']=diagnose(m)
        if assembly:
            checks=raw['assembly_checks'];checks['policy']=self.a.policy;c['independent_loads']=independent_loads(checks,self.a.policy,d['offsets']);c['independent_weak_modes']=independent_weak(checks,m['weak_modes'])
            if d['steps']==0:c['continuum_hypotheses']=reference_comparison(checks,self.a.policy)['continuum']
        write_once(folder/'comparison.json',c);return m,c

    def observe_arithmetic(self,name,d,m,c):
        c['arithmetic']=self.arithmetic(name,d,m)
        c['boundary_uncertainty']=boundary_uncertainty(m,c['arithmetic'],c['rule_settings'])
        return c

    def diagnosis(self):
        records={}
        for nx in (1,2):
            for step in (0,32):
                name='retained-Nx%d-%d'%(nx,step);d,m,c=self.retained(name,nx,step)
                if step:self.observe_arithmetic(name,d,m,c)
                records[name]=c;self.s['retained']=records;self.save()
        self.s['components']['retained_diagnosis']=component('PASS',scope='immutable startup/Tx states, actual producer callbacks; no governing reconstruction')

    def control(self):
        z,r,f=self.native('control-L1-fine-cp0',dt=1.6e-5,mesh='L1',role='adverse_dt_diagnostic');self.valid(z,r)
        m,c0=self.measure('control-L1-fine-cp0',z,f,mesh='L1',assembly=True)
        d,r,f=self.native('control-L1-fine-64',64,z['checkpoint'],dt=1.6e-5,mesh='L1',role='adverse_dt_diagnostic',arithmetic_time=.001024);self.valid(d,r)
        m,c=self.measure('control-L1-fine-64',d,f,mesh='L1',assembly=True);self.observe_arithmetic('control-L1-fine-64',d,m,c)
        retained=self.s['retained']['retained-Nx2-32'];old=read(retained['path'])['measurements'];diff=comparison(old,m,self.a.policy,c['rule_settings'])
        self.s['adverse_dt_control']=dict(initial=c0,terminal=c,comparison=diff,engineering_result=c['errors']['accuracy_gate'],
            interpretation='one authorized diagnostic block; no accuracy-preservation claim or further adverse continuation')
        self.s['components']['adverse_dt_control']=component('PASS',diagnostic_completed=True,accuracy_outcome=c['errors']['accuracy_gate']['status']);self.save()

    def screen(self,mesh):
        if mesh=='M2':
            siblings=[read(self.a.output.resolve()/p/'summary.json') for p in ('EP','Budget')]
            require(any(s.get('screens',{}).get('M1',{}).get('terminal',{}).get('errors',{}).get('accuracy_gate',{}).get('status')!='PASS' for s in siblings),'M2 is not authorized when both M1 policies qualify')
        name='screen-'+mesh
        z,r,f=self.native(name+'-cp0',mesh=mesh);self.valid(z,r)
        require(len(z['state'])==dict(M1=2447,M2=2975)[mesh],'unexpected custom mesh size')
        m,c0=self.measure(name+'-cp0',z,f,mesh=mesh,assembly=True)
        d,r,f=self.native(name+'-32',32,z['checkpoint'],mesh=mesh,arithmetic_time=.001024);self.valid(d,r)
        m,c=self.measure(name+'-32',d,f,mesh=mesh,assembly=True);self.observe_arithmetic(name+'-32',d,m,c)
        # Independent normal AND horizontal observation quadrature increase.
        qm,qc=self.measure(name+'-32-q8x6',d,f,mesh=mesh,order=8,xorder=6,diagnostics=False)
        quadrature={}
        for kind in ('continuum_error','initial_drift'):
            for scope in ('full','central','constrained_boundary'):
                for field,row in c['errors'][kind][scope]['fields'].items():
                    error=abs(row['L2']-qc['errors'][kind][scope]['fields'][field]['L2']);require(error<=1e-10,'observation quadrature sensitivity exceeds declared absolute floor')
                    quadrature[kind+'/'+scope+'/'+field]=error
        entry=dict(initial=c0,terminal=c,unknowns=len(d['state']),run_seconds=r['seconds'],seconds_per_step=r['seconds']/32,
            quadrature=dict(order=[6,8],xorder=[4,6],absolute_floor=1e-10,changes=quadrature),status=c['errors']['accuracy_gate']['status'])
        self.s.setdefault('screens',{})[mesh]=entry;self.s['components']['local_screen']=component(entry['status'],mesh=mesh);self.save()

    def select(self):
        siblings={p:read(self.a.output.resolve()/p/'summary.json') for p in ('EP','Budget')};selected=None
        for mesh in ('M1','M2'):
            if all(s.get('screens',{}).get(mesh,{}).get('status')=='PASS' for s in siblings.values()):selected=mesh;break
        design=dict(schema='WJ5CommonSelectionV1',selected=selected,status='PASS' if selected else 'OBSTRUCTED',
            screen_identities={p:{m:dict(status=r['status'],comparison_sha256=r['terminal']['sha256'],zknots=r['terminal']['zknots']) for m,r in s.get('screens',{}).items()} for p,s in siblings.items()},
            target=self.design['selected'],reason='one common vector; worse outcome across both policies; unchanged budgets')
        write_once(self.a.output.resolve()/'selection.json',design);return design

    def resource_design(self,mesh):
        if (self.root/'remaining-design.json').exists():
            packet=read(self.root/'remaining-design.json')
            require(packet['mesh']==mesh and packet['target']==self.design['selected'],'frozen remaining design differs')
            totals=self.ledger.totals();require(totals['attempts']<=800 and totals['seconds']<=10800,'RESOURCE_LIMIT: resumed package budget exhausted')
            self.s['components']['resource_design']=component('PASS',**packet);return
        screen=self.s['screens'][mesh];cost=screen['seconds_per_step'];totals=self.ledger.totals()
        remaining_attempts=128+96+256+33
        # Measured Nx2 screened process cost; Nx1 is conservatively charged at
        # the same per-step cost. Fine dt uses the same supported dense size.
        estimate=1.2*cost*remaining_attempts+800
        packet=dict(mesh=mesh,measured_seconds_per_step=cost,remaining_attempts=remaining_attempts,
            remaining_seconds_estimate=estimate,observation_SVD_IO_reserve=800.,safety_factor=1.2,already_charged=totals,
            attempt_limit=800,wall_limit=10800,target=self.design['selected'])
        write_once(self.root/'remaining-design.json',packet)
        require(totals['attempts']+remaining_attempts<=800 and totals['seconds']+estimate<=10800,'RESOURCE_LIMIT: frozen selected comparison estimate exceeds work cap')
        self.s['components']['resource_design']=component('PASS',**packet);self.save()

    def full(self):
        selection=read(self.a.output.resolve()/'selection.json');mesh=selection['selected']
        if not mesh:
            self.s['components']['selected_horizontal']=component('OBSTRUCTED','neither supported local screen qualified');self.s['components']['selected_timestep']=component('NOT_RUN','screen obstruction');return
        self.resource_design(mesh);traces=self.s.setdefault('selected',{})
        for tag,nx,dt,block,target in (('Nx1',1,3.2e-5,32,128),('Nx2',2,3.2e-5,32,128),('Nx2half',2,1.6e-5,64,256)):
            trajectory=[]
            if tag=='Nx2':
                entry=self.s['screens'][mesh];initial=entry['initial'];trajectory=[entry['terminal']]
                d=read(entry['terminal']['report']);cp=d['checkpoint'];start=32
            else:
                z,r,f=self.native('selected-'+tag+'-cp0',dt=dt,mesh=mesh,nx=nx,role='selected_'+tag);self.valid(z,r)
                m,initial=self.measure('selected-'+tag+'-cp0',z,f,mesh=mesh,nx=nx,assembly=True);cp=z['checkpoint'];start=0
            cp0=checkpoint_payload(initial['checkpoint'])
            traces[tag]=dict(initial=initial,trajectory=trajectory,dt=dt,nx=nx,mesh=mesh,
                initial_primary_sha256=state_hash(cp0['state']),initial_history_sha256=digest([state_hash(v) for v in cp0['native_history']]),
                native_dt_signature_sha256=digest(cp0['signature']))
            self.save()
            for end in range(start+block,target+1,block):
                name='selected-'+tag+'-'+str(end);at=end*dt if end in (block,target//2,target) else None
                d,r,f=self.native(name,block,cp,dt=dt,mesh=mesh,nx=nx,role='selected_'+tag,arithmetic_time=at)
                if r['exit_code'] or not d or d.get('failure'):
                    traces[tag]['hard_obstruction']=dict(report=str(f/'case.evolution.json'),exit_code=r['exit_code'],reason=d.get('failure') if d else r.get('failure'));break
                self.valid(d,r);m,c=self.measure(name,d,f,mesh=mesh,nx=nx,assembly=True)
                if at is not None:self.observe_arithmetic(name,d,m,c)
                trajectory.append(c);traces[tag]['trajectory']=trajectory;self.save()
                if c['errors']['accuracy_gate']['status']!='PASS':traces[tag]['accuracy_obstruction']=dict(step=end,gate=c['errors']['accuracy_gate']);break
                cp=d['checkpoint']
        self.compare_selected();self.save()

    def compare_selected(self):
        traces=self.s['selected'];mesh=traces['Nx2']['mesh'];components={}
        # Complete WJ-4 gate extension: all startup/cadence times, fields/scopes.
        for key,a,b in (('selected_horizontal','Nx1','Nx2'),('selected_timestep','Nx2','Nx2half')):
            ta,tb=traces[a],traces[b];components[key]={};pairs={};scope_uncertainty={};global_floors=arithmetic_floors({a:ta,b:tb})
            for i in range(5):
                ca=ta['initial'] if i==0 else next((c for c in ta['trajectory'] if c['steps']==i*32),None)
                cb=tb['initial'] if i==0 else next((c for c in tb['trajectory'] if c['steps']==i*(64 if b=='Nx2half' else 32)),None)
                if not ca or not cb:continue
                require(abs(ca['time']-cb['time'])<=2e-15,'common physical times differ')
                ma=read(ca['path'])['measurements'];mb=read(cb['path'])['measurements'];pairs[str(i)]=comparison(ma,mb,self.a.policy,ca['rule_settings'])
            original=sensitivity_sequence(pairs,global_floors,self.design,[str(i) for i in range(5)])
            # Keep the original global estimator unchanged alongside the tighter
            # constraint-aware boundary estimate. Use the maximum over declared
            # early/middle/end packets and all measured boundary rows.
            boundary={f:dict(L2=0.,sampled_Linf=0.) for f in FIELDS+('velocity',)}
            for t in (ta,tb):
                for c in t['trajectory']:
                    if 'boundary_uncertainty' not in c:continue
                    for f,r in c['boundary_uncertainty']['fields'].items():
                        for norm in ('L2','sampled_Linf'):boundary[f][norm]=max(boundary[f][norm],r[norm])
            adjusted=sensitivity_sequence(pairs,global_floors,self.design,[str(i) for i in range(5)],scope_uncertainty={'constrained_boundary':boundary})
            targetok=all(len(t['trajectory'])==4 and not t.get('hard_obstruction') and not t.get('accuracy_obstruction') and t['initial']['errors']['accuracy_gate']['status']=='PASS' for t in (ta,tb))
            arithmetic_ok=all(all(any(c['steps']==step and c.get('arithmetic') for c in t['trajectory']) for step in ([64,128,256] if tag=='Nx2half' else [32,64,128])) for tag,t in ((a,ta),(b,tb)))
            gate=aggregate(dict(target=component('PASS' if targetok else 'OBSTRUCTED','complete frozen targets required'),arithmetic=component('PASS' if arithmetic_ok else 'NOT_RUN'),sensitivity=adjusted))
            self.s[key+'_comparisons']=pairs;self.s[key+'_uncertainty']=dict(global_floors=global_floors,constraint_aware_boundary=boundary,original_global_gate=original)
            self.s['components'][key]=gate
        a,b=traces['Nx2'],traces['Nx2half'];require(a['initial_primary_sha256']==b['initial_primary_sha256'] and a['initial_history_sha256']==b['initial_history_sha256'] and a['native_dt_signature_sha256']!=b['native_dt_signature_sha256'],'fresh dt pair startup physical fields/history/signatures differ')

    def transactions(self):
        t=self.s.get('selected',{}).get('Nx2',{});rows=t.get('trajectory',[])
        if not rows or rows[-1]['steps']!=128 or t.get('accuracy_obstruction') or t.get('hard_obstruction'):
            self.s['components']['transactions']=component('NOT_RUN','qualified selected Nx2 target unavailable');return
        mesh=t['mesh'];whole=read(rows[-1]['report']);pre=read(rows[-2]['report']);cp=pre['checkpoint'];immutable=artifacts([Path(cp)/'State.bin',Path(rows[-2]['report']),Path(rows[-1]['report'])])
        a,r,f=self.native('transaction-first16',16,cp,mesh=mesh,role='transaction');self.valid(a,r)
        b,r,g=self.native('transaction-second16',16,a['checkpoint'],mesh=mesh,role='transaction');self.valid(b,r)
        exact=verify_transaction(whole,a,b)
        self.measure('transaction-actual-producer',pre,Path(rows[-2]['report']).parent,mesh=mesh,diagnostics=False)
        frozen,r,h=self.native('transaction-frozen',0,cp,mesh=mesh,role='transaction');self.valid(frozen,r);require(same_history(pre,frozen),'frozen readback changed accepted state/history/time')
        rejected,r,j=self.native('transaction-reject',1,cp,mesh=mesh,reject=97,role='transaction')
        require(rejected and 'deliberate_PIMPLE_WJ_hard_audit_rejection' in rejected.get('failure',''),'wrong deliberate rejection reason');require(same_history(pre,rejected),'rejection changed accepted state/history/time')
        for key in ('state_audits','wj_audit_packets'):require(digest_object(frozen[key])==digest_object(rejected[key]),'rollback audit packet mismatch')
        require(digest_object(pre['physical_samples'])==digest_object(rejected['physical_samples']),'rollback physical sample mismatch')
        require(state_hash(rejected['state'])!=state_hash(rejected['rejected_state']),'failed candidate not separate')
        for key in ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls'):require(rejected['pimple_attempts'][-1][key]==1,'rejected scheduled solve count differs')
        verify_artifacts(immutable);self.s['components']['transactions']=component('PASS',segment=[96,128],exact_primary_gauge_history_time=True,exact_accepted_audit_samples=True,all_step=exact,
            immutable_inputs=immutable,producer=rows[-2]['report'],rejected_report=str(j/'case.evolution.json'),failed_candidate_sha256=state_hash(rejected['rejected_state']));self.save()

    def finish(self):
        self.s['qualification']=aggregate(self.s['components']);self.s['stationarity']='NOT_TARGETED; WJ-4 NOT_REACHED remains unchanged';self.save()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--solver',type=Path,required=True);p.add_argument('--probe',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--retained',type=Path,required=True);p.add_argument('--policy',choices=('EP','Budget'),required=True);p.add_argument('--mode',choices=('diagnosis','control','screen','screen2','select','full','transactions','report'),required=True);a=p.parse_args()
    a.solver=a.solver.resolve();a.probe=a.probe.resolve();s=Study(a)
    try:
        if a.mode=='screen':s.screen('M1')
        elif a.mode=='screen2':s.screen('M2')
        elif a.mode=='select':s.select()
        elif a.mode!='report':getattr(s,a.mode)()
        s.finish()
    except Exception as e:
        s.s['components'][a.mode]=component('RESOURCE_LIMIT' if 'RESOURCE_LIMIT' in str(e) else 'FAIL',str(e));s.s.setdefault('failures',[]).append(dict(mode=a.mode,reason=str(e),measurement_code=s.measure_code));s.finish();raise
    finally:s.ledger.finish_session(s.session);s.save()

if __name__=='__main__':main()
