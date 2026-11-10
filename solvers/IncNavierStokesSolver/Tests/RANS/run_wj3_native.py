#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded compiled-native WJ spatial/drift study, with resumable immutable runs.

Modes are explicit gates. 'scout' freezes the time target and practical budget
before 'main' is permitted. Existing command receipts are read, never rerun.
"""
import argparse
import math
import os
from pathlib import Path
import platform
import subprocess
import time
import uuid
import xml.etree.ElementTree as ET
import numpy as np
from numpy.testing import assert_allclose
from wj2_cases import case_xml
from wj3_measure import observations,errors,compare,replay_witnesses,arithmetic,FIELDS,wall,K0,checkpoint_payload
from test_wj2_native import digest_object,independent_loads,reference_comparison,resting_functions
from compare_cls_epawf_pimple import same_history
from cls_epawf_m3_cases import set_info,set_parameter,write
from run_cls_epawf_picard import THREADS,load,state_hash,sha,require,write_json

COUNTS=('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls')


class Study:
    def __init__(self,args):
        require(platform.system()=='Linux','Linux compiled native execution required')
        self.a=args;self.root=args.output.resolve()/args.policy;self.root.mkdir(parents=True,exist_ok=True)
        self.file=self.root/'summary.json';self.started=time.time();os.environ.update(THREADS)
        self.s=load(self.file) if self.file.exists() else dict(schema='WJ3BoundedStudyV1',policy=args.policy,
            source_revision=subprocess.check_output(['git','-C',str(Path(__file__).resolve().parents[4]),'rev-parse','HEAD'],text=True).strip(),
            solver_sha256=sha(args.solver),probe_sha256=sha(args.probe),receipts=[],scheduled_attempts=0,
            attempt_limit=500,ranks=1,threads=THREADS,arithmetic_policy='-ffp-contract=off -fno-fast-math',
            physical_duration_units='seconds; dimensionless duration T*u*/H = T/2')
        require(self.s['solver_sha256']==sha(args.solver) and self.s['probe_sha256']==sha(args.probe),'runtime changed during immutable study; use a new root')
        if self.s.get('failure'):
            self.s.setdefault('resolved_harness_failures',[]).append(dict(reason=self.s['failure'],prior_driver_revision=self.s.get('comparison_driver_revision')))
        self.s['comparison_driver_revision']=subprocess.check_output(['git','-C',str(Path(__file__).resolve().parents[4]),'rev-parse','HEAD'],text=True).strip()
        self.s['measurement_driver_sha256']=sha(Path(__file__))
        self.s['physical_comparison_sha256']=sha(Path(__file__).with_name('wj3_measure.py'))
        self.save()

    def save(self):write_json(self.file,self.s)

    def execute(self,command,folder):
        receipt=folder/'receipt.json'
        if receipt.exists():return load(receipt)
        start=time.time()
        with (folder/'command.log').open('w') as log:
            result=subprocess.run(list(map(str,command)),cwd=folder,stdout=log,stderr=subprocess.STDOUT,env=os.environ.copy(),timeout=3600)
        r=dict(command=list(map(str,command)),exit_code=result.returncode,seconds=time.time()-start,
            log=str(folder/'command.log'),log_sha256=sha(folder/'command.log'))
        write_json(receipt,r);self.s['receipts'].append(r);self.save();return r

    def native(self,name,steps=0,cp=None,dt=1e-6,level=1,normal=16,reject=0,arithmetic_time=None,floor_relative=None,resting=False):
        require(level in (0,1,2),'unsupported WJ normal refinement level')
        level=int(level)
        folder=self.root/name;xml=folder/'case.xml'
        if not xml.exists():
            if folder.exists() and not any(folder.iterdir()):folder.rmdir()
            require(self.s['scheduled_attempts']+steps<=500,'500-attempt WJ work cap would be exceeded')
            xml=case_xml(folder,self.a.policy,steps,cp,level=level,normal=normal,dt=dt,hard_reject=reject)
            tree=ET.parse(xml).getroot();set_info(tree,'RANSWallPimpleSnapshotDirectory',str(folder/'pimple-states'))
            set_info(tree,'RANSWallPimpleStopOnBound','On')
            if floor_relative is not None:set_parameter(tree,'RANSWallPimpleKFloorRelative',floor_relative)
            if resting:resting_functions(tree)
            if arithmetic_time is not None:
                set_info(tree,'RANSWallPimpleArithmetic','On');set_info(tree,'RANSWallPimpleArithmeticTime',format(arithmetic_time,'.17g'))
            write(tree,xml)
        r=self.execute([self.a.solver,xml],folder)
        report=load(folder/'case.evolution.json') if (folder/'case.evolution.json').exists() else None
        index=folder/'identity.json'
        if not index.exists():
            attempts=report.get('pimple_attempts',[]) if report else []
            self.s['scheduled_attempts']+=len(attempts)
            write_json(index,dict(xml_sha256=sha(xml),report_sha256=sha(folder/'case.evolution.json') if report else None,
                scheduled_attempts=len(attempts),requested_steps=steps,dt=dt,level=level,
                state_sha256=state_hash(report['state']) if report and 'state' in report else None,
                history_sha256=digest_object(report['native_history']) if report and 'native_history' in report else None,
                checkpoint=report.get('checkpoint') if report else None,
                checkpoint_sha256=sha(Path(report['checkpoint'])/'State.bin') if report and report.get('checkpoint') else None))
            self.save()
        return report,r,folder

    def valid(self,report,receipt):
        require(report is not None and receipt['exit_code']==0 and not report.get('failure'),
            'native obstruction: '+str(report.get('failure') if report else receipt))
        require(report['initial_k_bound_check']['activations']==0,'startup MassProjection changed reference')
        for packet in report['wj_audit_packets']:
            require(not packet['failure_categories'],'WJ hard finite-inventory audit failed')
        for attempt in report['pimple_attempts']:
            require(attempt['scheduled_complete'] and all(attempt[x]==1 for x in COUNTS),'incomplete or extra scheduled solves')
            require(attempt['additional_governing_equation_solves']==0,'governing diagnostic solve entered trajectory')
            require(attempt['bound_activations']==0,'material MassProjection: unmodified-equation baseline stopped')
            require(attempt['refreshed_trace_defect']<1e-11,'final endpoint-K epsilon trace mismatch')
            snapshot=load(Path(attempt['k_snapshot']))
            require(state_hash(snapshot['raw'])==state_hash(snapshot['corrected']),'inactive bound changed candidate bytes')

    def measure(self,name,report,casefolder,assembly=False,order=6):
        folder=self.root/name;folder.mkdir(exist_ok=True);cp=Path(report['checkpoint']);native=casefolder/'case.evolution.json'
        hashes={str(f.relative_to(cp)):sha(f) for f in cp.rglob('*') if f.is_file()}
        if not (folder/'input.json').exists():
            write_json(folder/'input.json',dict(operation='wj_measure',inner_stress_model='Linear',native_checkpoint=str(cp),
                checkpoint_sha256=sha(cp/'State.bin'),native_signature=report['signature'],state=report['state'],history=report['state'],
                time=report['time'],alpha=0.,native_differential_state=report['state'][:int(report['offsets'][4])],
                native_history=report['native_history'],native_history_times=report['native_history_times'],
                native_report=str(native),native_report_sha256=sha(native),observations=observations(order),export_assembly_checks=assembly))
        r=self.execute([self.a.probe,casefolder/'case.xml',folder/'input.json',folder/'result.json'],folder)
        require(r['exit_code']==0,'strict native readback failed: '+str(folder/'command.log'))
        result=load(folder/'result.json');m=result['measurements'];require(m['state_sha256']==m['state_after_sha256'] and
            m['stage_history_sha256']==m['stage_history_after_sha256'] and m['governing_solves']==0,'measurement mutated/solved')
        require(hashes=={str(f.relative_to(cp)):sha(f) for f in cp.rglob('*') if f.is_file()},'measurement changed saved bytes')
        compact=dict(path=str(folder/'result.json'),sha256=sha(folder/'result.json'),
            spatial_mass_dual_norms=m['spatial_mass_dual_norms'],completed_stage=m['completed_stage']['mass_dual_norms'] if m['completed_stage'] else None,
            state_sha256=m['state_sha256'],stage_history_sha256=m['stage_history_sha256'],
            field_errors=errors(m,self.a.policy),witness_replay=replay_witnesses(m['audit']))
        if assembly:
            compact['independent_loads']=independent_loads(result['assembly_checks'],self.a.policy,report['offsets'])
            compact['independent_profiles']=reference_comparison(result['assembly_checks'],self.a.policy)
        write_json(folder/'comparison.json',compact)
        return m,compact

    def intermediate(self,name,report,folder,step):
        candidates=[]
        for cp in folder.glob('checkpoint*'):
            if cp.is_dir() and (cp/'State.bin').exists():
                payload=checkpoint_payload(cp)
                if payload['steps']==step:candidates.append(payload)
        require(candidates,'missing immutable intermediate checkpoint')
        payload=candidates[0];view=dict(report);view.update(payload)
        return self.measure(name,view,folder)

    def smoke(self):
        z,r,f=self.native('smoke-cp0',level=0);self.valid(z,r)
        s,r,sf=self.native('smoke-step',1,z['checkpoint'],level=0,arithmetic_time=1e-6);self.valid(s,r)
        m,c=self.measure('smoke-measure',s,sf)
        require(m['completed_stage'] is not None and m['completed_stage']['alpha_mass']==1e6,'actual stage not measured')
        u,r,_=self.native('smoke-observer-off',1,z['checkpoint'],level=0);self.valid(u,r)
        require(same_history(s,u),'passive arithmetic observations changed trajectory')
        # A corrupted report binding fails its intended guard without altering
        # any saved report/checkpoint. It is a new input packet only.
        folder=self.root/'smoke-binding-negative';folder.mkdir(exist_ok=True)
        payload=load(self.root/'smoke-measure/input.json');payload['native_report_sha256']='0'*64
        write_json(folder/'input.json',payload)
        bad=self.execute([self.a.probe,sf/'case.xml',folder/'input.json',folder/'result.json'],folder)
        require(bad['exit_code']!=0 and 'WJ accepted report SHA256 differs' in (folder/'command.log').read_text(),'wrong stage binding rejection')
        near=1-1e-12;large=1.28e-4
        iz,ir,_=self.native('smoke-projection-cp0',level=0,dt=large,floor_relative=near,resting=True);self.valid(iz,ir)
        stopped,sr,stf=self.native('smoke-projection-stop',4,iz['checkpoint'],level=0,dt=large,floor_relative=near,resting=True)
        require(stopped and 'WJ_baseline_material_mass_projection_exit' in stopped.get('failure',''), 'wrong material-bound stopping reason')
        require(same_history(iz,stopped) and len(stopped['pimple_attempts'])==1,'material projection advanced history or ran further scheduled attempts')
        require(all(stopped['pimple_attempts'][0][c]==1 for c in COUNTS) and stopped['pimple_attempts'][0]['bound_activations']>0,
            'projection stop did not retain complete attempt/correction accounting')
        self.s['smoke']=dict(status='PASS',comparison=c,passive_observer_bitwise=True,arithmetic=arithmetic(s,m),
            material_projection_stop=dict(result=str(stf/'case.evolution.json'),accepted_unchanged=True,attempts=1))
        self.save()

    def spatial(self):
        rows=[]
        for level in (0,1,2):
            z,r,f=self.native('spatial-L'+str(level),level=level);self.valid(z,r)
            require(len(z['state'])==(639,1257,2493)[level],'unexpected tensor-cubic mesh size')
            m,c=self.measure('spatial-measure-L'+str(level),z,f,assembly=True)
            first,fr,ff=self.native('spatial-first-step-L'+str(level),1,z['checkpoint'],level=level)
            self.valid(first,fr);fm,fc=self.measure('spatial-first-measure-L'+str(level),first,ff)
            rows.append(dict(level=level,unknowns=len(z['state']),raw_reference_norms=z['raw_reference_norms'],comparison=c,
                first_step_pressure_adjustment=fc['field_errors']['initial_drift'],first_step_comparison=fc))
            self.s['spatial']=rows;self.save()
        expected_k=[.013013171987,.002727140552] if self.a.policy=='EP' else [.270481784118,.063982673543]
        for i in (0,1):
            assert_allclose(np.array(rows[i]['raw_reference_norms'])[[0,2,3]],
                [[.268457088879,.118646107517][i],expected_k[i],[.004010426499,.000605967544][i]],rtol=2e-10,atol=1e-12)
        z,r,f=self.native('spatial-L1-source32',level=1,normal=32);self.valid(z,r)
        m,c=self.measure('spatial-L1-source32-measure',z,f)
        base=np.array(rows[1]['raw_reference_norms']);delta=np.array(z['raw_reference_norms'])-base
        require(max(abs(delta))<1e-10,'source quadrature influence is not negligible')
        coarse=load(Path(rows[1]['comparison']['path']))['measurements']
        qm,qc=self.measure('spatial-L1-observation10',z,f,order=10)
        # Error-integral sensitivity is compared as scalar norms, since rules
        # differ; coefficients are never subtracted across meshes/rules.
        base_errors=errors(m,self.a.policy)['continuum_error']['fields'];new_errors=errors(qm,self.a.policy)['continuum_error']['fields']
        qdelta={name:new_errors[name]['L2']-base_errors[name]['L2'] for name in FIELDS}
        require(max(abs(v) for v in qdelta.values())<1e-10,'observation quadrature error is not negligible')
        residual=np.array([row['raw_reference_norms'] for row in rows]);rates=[]
        for j in (0,1):
            rates.append(dict(from_level=j,to_level=j+1,factors=(residual[j]/np.maximum(residual[j+1],1e-300)).tolist(),
                slopes_log2=[float(math.log2(a/b)) if a>1e-11 and b>1e-11 else None for a,b in zip(residual[j],residual[j+1])]))
        self.s['spatial_quadrature']=dict(source_order=[16,32],basis_points=[8,8],independent_source_order=32,
            independent_outer_order=24,observation_normal_orders=[6,10],source_norm_delta=delta.tolist(),
            L2_integral_delta=qdelta,absolute_floor=1e-10)
        self.s['spatial_rates']=rates;self.save()

    def scout(self):
        scouts=[];chosen=None;initial=None
        for dt in (1e-6,4e-6,1.6e-5,6.4e-5,1.28e-4):
            tag=format(dt,'.9g');z,r,f=self.native('scout-'+tag+'-cp0',dt=dt);self.valid(z,r)
            if initial is None:initial=z
            require(state_hash(z['state'])==state_hash(initial['state']) and digest_object(z['native_history'])==digest_object(initial['native_history']) and
                z['time']==initial['time']==0. and z['dt']==dt,'dt CP0 changed initial primary/gauge/physical history')
            s,rc,sf=self.native('scout-'+tag,4,z['checkpoint'],dt=dt,arithmetic_time=4*dt)
            status='admissible'
            try:self.valid(s,rc)
            except ValueError as error:status=str(error)
            row=dict(dt=dt,status=status,seconds=rc['seconds'],seconds_per_step=rc['seconds']/4,
                initial_sha256=state_hash(z['state']),initial_history_sha256=digest_object(z['native_history']),
                native_initial_history_labels=z['native_history_times'],native_dt_signature_sha256=digest_object(z['signature']),
                result=str(sf/'case.evolution.json'),accepted_steps=s['steps'] if s else None)
            if s:
                row['last_accepted_sha256']=state_hash(s['state']);row['failure']=s.get('failure')
                row['failed_candidate_sha256']=state_hash(s['rejected_state']) if s.get('rejected_state') else None
                if s.get('physical_range_events'):row['witness']=s['physical_range_events'][-1]
            scouts.append(row);self.s['scouts']=scouts;self.save()
            if status!='admissible':break
            chosen=dt
        require(chosen is not None,'no fully admissible four-step scout')
        cost=max(s['seconds_per_step'] for s in scouts if s['status']=='admissible')
        fine_cost=load(self.root/'spatial-first-step-L2/receipt.json')['seconds']
        fine_estimate=257*cost+128*fine_cost+120
        fine_fits=fine_estimate<3600
        # Declared BEFORE main trajectories. Cross-mesh level0 keeps this
        # bounded dense-backend pilot within a 60-minute per-policy remaining trajectory cap.
        # Level2 remains spatially qualified, temporally open.
        self.s['frozen_design']=dict(D=chosen,T=32*chosen,T_ustar_over_H=16*chosen,
            main_steps=[32,64,128],main_dt=[chosen,chosen/2,chosen/4],mesh_unknowns=1257,
            cross_mesh_level=2 if fine_fits else 0,
            cross_mesh_unknowns=2493 if fine_fits else 639,cross_mesh_steps=128,
            measured_max_seconds_per_step=cost,estimated_remaining_seconds=fine_estimate if fine_fits else 385*cost+120,
            measured_fine_one_step_process_seconds=fine_cost,
            practical_wall_seconds_per_policy=3600,safety_factor_over_arithmetic_floor=100,
            floor_policy='terminal frozen matrix SVD and long-double residual, arithmetic evaluation sensitivity, accumulated using frozen history-map mass-norm amplification and safety100; no spatial/stage residual used as noise',
            pressure_policy='common-time common overall gauge; algebraic first-step adjustment not required to scale with dt',
            frozen_before_main_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
            finer_temporal_mesh_open=not fine_fits,
            cross_mesh_cost_estimator='257 times measured 1257 scout step cost plus128 times measured 2493 full one-step process cost plus120 seconds observation reserve; choose2493 if below3600 seconds')
        self.save()

    def main_runs(self):
        design=self.s['frozen_design'];D,T=design['D'],design['T'];results=[];previous_initial=None
        for j,n in enumerate((32,64,128)):
            dt=D/2**j;z,r,f=self.native('main-'+str(n)+'-cp0',dt=dt);self.valid(z,r)
            if previous_initial is not None:
                require(state_hash(z['state'])==state_hash(previous_initial['state']) and digest_object(z['native_history'])==digest_object(previous_initial['native_history']),
                    'fresh dt startup primary/history mismatch')
            previous_initial=z
            s,rc,sf=self.native('main-'+str(n),n,z['checkpoint'],dt=dt,arithmetic_time=T)
            try:self.valid(s,rc)
            except ValueError as error:
                obstruction=dict(reason=str(error),frozen_T=T,result=str(sf/'case.evolution.json'),
                    accepted_steps=s['steps'] if s else None,dt=dt)
                self.s.setdefault('main_obstructions',[]).append(obstruction)
                results.append(dict(requested_steps=n,steps=s['steps'] if s else 0,dt=dt,result=str(sf/'case.evolution.json'),obstruction=obstruction))
                self.s['main']=results;self.save();continue
            require(abs(s['time']-T)<1e-15,'common-time trajectory missed frozen target')
            m,c=self.measure('main-'+str(n)+'-measure',s,sf)
            floor=arithmetic(s,m);results.append(dict(steps=n,dt=dt,comparison=c,arithmetic=floor,result=str(sf/'case.evolution.json')))
            self.s['main']=results;self.save()
        if any('obstruction' in row for row in results):
            common_steps=min(int(round(row['steps']*row['dt']/D)) for row in results)
            self.s['common_earlier_time']=common_steps*D
            if common_steps>0:
                earlier=[]
                for j,row in enumerate(results):
                    report=load(Path(row['result']));folder=Path(row['result']).parent
                    m,c=self.intermediate('main-earlier-'+str(j),report,folder,common_steps*2**j);earlier.append(m)
                self.s['earlier_comparisons']=dict(A=compare(earlier[0],earlier[1],self.a.policy),B=compare(earlier[1],earlier[2],self.a.policy),
                    interpretation='frozen target obstructed; common earlier time retained; no rate without arithmetic floor')
            self.save();return
        measurements=[load(Path(row['comparison']['path']))['measurements'] for row in results]
        A=compare(measurements[0],measurements[1],self.a.policy);B=compare(measurements[1],measurements[2],self.a.policy)
        rates={}
        for f in FIELDS+('velocity',):
            a,b=A['fields'][f]['L2'],B['fields'][f]['L2'];floor=max(row['arithmetic']['estimated_comparison_floor'][f] for row in results)
            rates[f]=dict(A=a,B=b,comparison_floor=floor,resolved=min(a,b)>floor,
                p=math.log2(a/b) if min(a,b)>floor else None)
        self.s['temporal']=dict(A=A,B=B,rates=rates)
        half=[]
        for j,row in enumerate(results):
            report=load(Path(row['result']));m,c=self.intermediate('main-half-'+str(j),report,Path(row['result']).parent,(16,32,64)[j]);half.append(m)
        self.s['common_half_time']=dict(time=T/2,A=compare(half[0],half[1],self.a.policy),B=compare(half[1],half[2],self.a.policy))
        self.save()
        # Same T, finest dt, independently projected other mesh.
        z,r,f=self.native('cross-cp0',dt=D/4,level=design['cross_mesh_level']);self.valid(z,r)
        s,rc,sf=self.native('cross-128',128,z['checkpoint'],dt=D/4,level=design['cross_mesh_level'],arithmetic_time=T)
        try:self.valid(s,rc)
        except ValueError as error:self.s['cross_obstruction']=dict(reason=str(error),result=str(sf/'case.evolution.json'));self.save();return
        m,c=self.measure('cross-128-measure',s,sf)
        self.s['cross_mesh']=dict(comparison=c,arithmetic=arithmetic(s,m),versus_main_finest=compare(m,measurements[-1],self.a.policy));self.save()
        self.transactions(D)

    def transactions(self,dt):
        whole=load(self.root/'main-32/case.evolution.json');cp0=load(self.root/'main-32-cp0/case.evolution.json')['checkpoint']
        first,r,f=self.native('restart-first16',16,cp0,dt=dt);self.valid(first,r)
        cp=Path(first['checkpoint']);before={str(p.relative_to(cp)):sha(p) for p in cp.rglob('*') if p.is_file()}
        second,r,sf=self.native('restart-second16',16,cp,dt=dt);self.valid(second,r)
        require(same_history(whole,second),'16+16 restart primary/gauge/history/time mismatch')
        for name in ('state_audits','wj_audit_packets'):
            require(digest_object(whole[name][-1])==digest_object(second[name][-1]),'restart accepted audit mismatch: '+name)
        require(digest_object(whole['physical_samples'])==digest_object(second['physical_samples']),'restart physical samples mismatch')
        self.measure('later-accepted-readback',first,f)
        readonly,rr,_=self.native('later-zero-restart',0,cp,dt=dt);self.valid(readonly,rr)
        rejected,r,rf=self.native('later-reject',1,cp,dt=dt,reject=17)
        require(rejected and 'deliberate_PIMPLE_WJ_hard_audit_rejection' in rejected.get('failure',''),'wrong end-of-step rejection')
        require(same_history(first,rejected),'hard rejection changed accepted primary/history/time')
        require(all(rejected['pimple_attempts'][-1][key]==1 for key in COUNTS),'rejection lost fully scheduled solve counts')
        require(state_hash(rejected['rejected_state'])!=state_hash(rejected['state']),'failed candidate not retained')
        require(digest_object(first['physical_samples'])==digest_object(rejected['physical_samples']),'rejection changed physical samples')
        require(digest_object(readonly['state_audits'])==digest_object(rejected['state_audits']) and digest_object(readonly['wj_audit_packets'])==digest_object(rejected['wj_audit_packets']), 'rejection changed accepted audit packets')
        require(before=={str(p.relative_to(cp)):sha(p) for p in cp.rglob('*') if p.is_file()},'transaction changed saved checkpoint bytes')
        self.s['transactions']=dict(exact_16_plus_16=True,later_frozen_readback=True,rollback_exact=True,
            checkpoint_hashes=before,failed_candidate_sha256=state_hash(rejected['rejected_state']),
            rejected_result=str(rf/'case.evolution.json'));self.save()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--solver',type=Path,required=True);p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--policy',choices=('EP','Budget'),required=True)
    p.add_argument('--mode',choices=('smoke','spatial','scout','main'),required=True);a=p.parse_args();a.solver=a.solver.resolve();a.probe=a.probe.resolve()
    if a.mode=='smoke':a.output=a.output/('run-'+uuid.uuid4().hex[:12])
    s=Study(a)
    try:getattr(s,{'main':'main_runs'}.get(a.mode,a.mode))();s.s[a.mode+'_gate']='OBSTRUCTED' if a.mode=='main' and s.s.get('main_obstructions') else 'PASS';s.s.pop('failure',None);s.save()
    except Exception as error:s.s[a.mode+'_gate']='FAIL';s.s['failure']=str(error);s.save();raise
    print(s.file)

if __name__=='__main__':main()
