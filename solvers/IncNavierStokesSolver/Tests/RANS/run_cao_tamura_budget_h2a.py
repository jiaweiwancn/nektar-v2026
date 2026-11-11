#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded H2a developing rough-flat qualification; never mapped production."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from copy import deepcopy
from datetime import datetime,timezone
import math
import os
from pathlib import Path
import platform
import subprocess
import time
import xml.etree.ElementTree as ET

from cao_tamura_budget_h1 import compatibility_case,HH,UH,USTAR
from compare_cls_epawf_pimple import same_history
from continue_cls_epawf_pimple import CampaignLock,atomic,freeze
from diagnose_cls_accepted_covariance import checkpoint,common_packet,run_probe
from run_cls_epawf_picard import THREADS,load,require,sha,state_hash,set_info,set_parameter,write


def prepare(folder,input,steps,dt=1e-5,restart=None,developing=True,optimized=True,packets='SelectedV1',interval=32):
    xml,definition=compatibility_case(folder,input/'input/finite-inflow.txt',input/'preflight.json',steps,dt=dt,restart=restart)
    tree=ET.parse(xml).getroot()
    for k,v in dict(RANSWallBoundaryPolicy='RoughFlatOpenV1' if developing else 'CompatibilityV1',
        RANSWallVelocityMassReuse='ExactV1' if optimized else 'Off',
        RANSWallFiniteReferenceReuse='ExactHeightV1' if optimized else 'Off',RANSWallAttemptPackets=packets).items():set_info(tree,k,v)
    set_parameter(tree,'RANSWallPacketInterval',interval)
    if developing:set_parameter(tree,'RANSWallOutletTotalPressure',0.)
    write(tree,xml);definition.update(version='CaoTamuraBudgetH2aV1',boundary='RoughFlatOpenV1' if developing else 'CompatibilityV1',
        case='rough-flat-finite-profile-developing-v1' if developing else 'rough-flat-finite-profile-compatibility-v1',
        boundary_scope='prescribed inlet, no-slip ground velocity with natural k, impermeable full-stress slip top, absolute-total-traction open outlet' if developing else
            'both streamwise boundaries and top fixed; H1 compatibility regression',xml_sha256=sha(xml),
        observation_x=0.,numerical_inlet_x=-.2,physical_hill_height_scale=HH,
        boundary_pressure='absolute selected total traction; no mean pressure gauge' if developing else 'compatibility common mean gauge',
        backflow='native fixed reference, upwind incoming turbulence and momentum; positive prescribed scalars' if developing else 'not applicable; prescribed streamwise sides',
        top='impermeable selected-stress slip, natural scalars, floor-only Yap' if developing else 'prescribed compatibility top',optimized=optimized,packets=packets)
    (xml.parent/'h1-definition.json').rename(xml.parent/'compatibility-template-definition.json')
    atomic(xml.parent/'definition.json',definition);atomic(xml.parent/'h2a-definition.json',definition);return xml


def physical(a,b):
    """Same mesh, native physical quadrature; pressure is ABSOLUTE, anchored."""
    aa=a['physical_samples'];bb=b['physical_samples'];out={}
    for i,name in enumerate(('U','W','k','epsilon','pi_inner','pi_outer')):
        w=aa['weights_'+str(i)];require(w==bb['weights_'+str(i)],'physical panels differ')
        out[name]=math.sqrt(max(0.,math.fsum(v*(x-y)**2 for v,x,y in zip(w,aa['values_'+str(i)],bb['values_'+str(i)]))))
    out['velocity']=math.hypot(out['U'],out['W']);out['pressure']=math.hypot(out['pi_inner'],out['pi_outer'])
    for name in ('wall_stress','wall_K'):
        w=aa['wall_weights'];require(w==bb['wall_weights'],'wall quadrature differs')
        out[name]=math.sqrt(max(0.,math.fsum(v*(x-y)**2 for v,x,y in zip(w,aa[name],bb[name]))))
    def mean(q):
        p=q['physical_samples'];return math.fsum(w*x for w,x in zip(p['wall_weights'],p['wall_stress']))/math.fsum(p['wall_weights'])
    out['mean_traction']=abs(mean(a)-mean(b));return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('solver','probe','h1','output'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--phase',choices=['focused','startup'],default='focused')
    p.add_argument('--primary-cap-seconds',type=float,default=1800)
    p.add_argument('--reuse-temporal-root',type=Path,help='reuse hash-validated completed temporal legs after the known missing-h0-samples reporting stop; never retry a numerical rejection')
    a=p.parse_args();require(platform.system()=='Linux','authoritative Linux native execution required')
    a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False);lock=CampaignLock(a.output)
    os.environ.update(THREADS);os.environ.update(VECLIB_MAXIMUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
    solver,image=freeze(a.solver.resolve(),a.output/'runtime');probe,probe_image=freeze(a.probe.resolve(),a.output/'probe-runtime')
    summary=dict(schema=1,phase=a.phase,status='started',started_UTC=datetime.now(timezone.utc).isoformat(),
        runtime=image,probe=probe_image,profile_sha256=sha(a.h1/'input/finite-inflow.txt'),
        h1=str(a.h1),threads={**THREADS,'VECLIB_MAXIMUM_THREADS':'1','NUMEXPR_NUM_THREADS':'1'},rank=1,tests={},observations={})
    def save(**q):summary.update(q);atomic(a.output/'progress.json',summary)
    save()
    def native(name,steps,dt=1e-5,cp=None,developing=True,optimized=True,packets='SelectedV1',interval=32,mutate=None,expect=True):
        folder=a.output/name;xml=prepare(folder,a.h1,steps,dt,cp,developing,optimized,packets,interval)
        if mutate:
            t=ET.parse(xml).getroot();mutate(t);write(t,xml)
        env=dict(os.environ,LD_LIBRARY_PATH=str(solver.parent/'lib'),TMPDIR=str(folder),NEKTAR_RANS_TIMING_FILE=str(folder/'timing.jsonl'))
        cmd=['/usr/bin/time','-f','%e %U %S %M','-o',str(folder/'usage.txt'),str(solver),str(xml)]
        before=time.monotonic()
        with (folder/'solver.log').open('w') as log:
            child=subprocess.Popen(cmd,cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT,pass_fds=(lock.stream.fileno(),))
            atomic(folder/'launch.json',dict(pid=child.pid,command=cmd,started_UTC=datetime.now(timezone.utc).isoformat()))
            code=child.wait(timeout=1800)
        file=folder/'case.evolution.json';q=load(file) if file.exists() else None
        receipt=dict(code=code,seconds=time.monotonic()-before,usage=(folder/'usage.txt').read_text().strip(),
            report=str(file),report_sha256=sha(file) if file.exists() else None,xml_sha256=sha(xml),failure=q.get('failure') if q else None,local_steps_requested=steps)
        atomic(folder/'receipt.json',receipt)
        if expect:
            require(code==0 and q and not q['failure'],'native rejection '+name+': '+str(receipt))
            require(q['evaluation_reuse']['mode']=='ExactStepV1' and q['audit_revision']=='NativeElementSideV2','hard audit/reuse policy differs')
            for row in q['pimple_attempts']:
                for k in ('outer_sweeps','momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves','bound_calls'):
                    require(row[k]==1,'wrong governing schedule '+k)
                require(row['additional_governing_equation_solves']==0 and row['solved_epsilon_retained_bitwise'],'extra PDE or epsilon replacement')
            if developing:require(abs(q['developing_boundary']['net_volume_flux'])<1e-10 and q['boundary_flux_correction']==0,'open weak mass balance/top correction failed')
        else:require(code!=0 or (q and q.get('failure')),'expected native rejection absent')
        return q,receipt
    def observe(name,q,fixture=False,x=0.,source=None):
        source=source or a.output
        cp=Path(q['checkpoint']);saved=checkpoint(cp);packet=dict(common_packet(cp),operation='finite_inflow_h1',
            state=saved['coupled'],history=saved['coupled'],time=saved['time'],alpha=0.,inner_stress_model='Linear',
            native_report_file=str(source/name/'case.evolution.json'),native_report_sha256=sha(source/name/'case.evolution.json'),
            permitted_queries_file=str(a.h1/'input/permitted-queries.txt'),permitted_queries_sha256=sha(a.h1/'input/permitted-queries.txt'),
            observation_x=x,boundary_fixture=fixture)
        result,receipt=run_probe(probe,source/name/'case.xml',packet,a.output/(name+'-readback'))
        require(result['governing_solves']==result['physical_commits']==0 and result['initialization_work']==[0,0],'readback advanced or lifted')
        atomic(a.output/(name+'-observation.json'),dict(receipt=receipt,diagnostic=result))
        return dict(path=str(a.output/(name+'-observation.json')),sha256=sha(a.output/(name+'-observation.json')))
    try:
        if a.phase=='focused':
            old=load(a.h1/'whole8/case.evolution.json');h0=load(a.h1/'initial/case.evolution.json');parent=Path(h0['checkpoint'])
            base,rb=native('H1-off',8,cp=parent,developing=False,optimized=False,packets='EveryAttempt')
            opt,ro=native('H1-on',8,cp=parent,developing=False,optimized=True,interval=8)
            require(same_history(old,base) and same_history(base,opt),'H1 compatibility or exact performance reuse changed native state/history')
            for before,after in zip(base['pimple_attempts'],opt['pimple_attempts']):
                for key in ('history_sha256','incoming_sha256','post_flow_sha256','pre_k_sha256','raw_k_sha256','corrected_k_sha256'):
                    require(before[key]==after[key],'cache/output changed native stage '+key)
            summary['tests']['H1_regression_and_optimization']=dict(status='passed',baseline=rb,optimized=ro,
                bitwise_state_history_stages=True,speedup=rb['seconds']/ro['seconds']);save()
            origin,_=native('normal-origin-rejected',0,developing=False,mutate=lambda t:set_parameter(t,'RANSWallReferenceOriginZ',-1e-5),expect=False)
            require(not origin or not origin.get('native_history'),'unsupported fresh origin formed histories')
            require('zero normal origin' in (a.output/'normal-origin-rejected/solver.log').read_text(),'normal-origin rejected for wrong reason')
            tangent,_=native('tangent-origin',0,developing=False,mutate=lambda t:set_parameter(t,'RANSWallReferenceOriginX',-.3))
            require(state_hash(tangent['state'])==state_hash(h0['state']),'tangential origin altered incoming height/seed')
            zero,_=native('open-h0',0);cp=Path(zero['checkpoint']);parent_hash=sha(cp/'State.bin')
            summary['observations']['open-h0']=observe('open-h0',zero,True);save()
            whole,rw=native('open8-on',8,cp=cp,interval=8);uncached,ru=native('open8-off',8,cp=cp,optimized=False,packets='EveryAttempt')
            require(same_history(whole,uncached),'developing exact-cache/output path changed native state/history')
            for before,after in zip(uncached['pimple_attempts'],whole['pimple_attempts']):
                for key in ('history_sha256','incoming_sha256','post_flow_sha256','pre_k_sha256','raw_k_sha256','corrected_k_sha256'):
                    require(before[key]==after[key],'open cache changed stage '+key)
            first,_=native('split4',4,cp=cp);second,_=native('restart4',4,cp=Path(first['checkpoint']))
            require(same_history(whole,second),'open split/uninterrupted restart differs')
            reject,_=native('rejection1',1,cp=cp,interval=1000,mutate=lambda t:set_parameter(t,'RANSWallPimpleHardAuditRejectStep',1),expect=False)
            require(same_history(zero,reject),'open hard-rejection rollback changed accepted state/history')
            require(reject['attempt_packet_status']['last_written_packet'] and sha(cp/'State.bin')==parent_hash,'rejection packet or parent preservation failed')
            # Reject at an UNSELECTED step after a real prefix, so the deferred
            # packet is flushed by native acceptance rollback rather than cadence.
            late=Path(whole['checkpoint']);startup,_=native('rollback-parent',0,cp=late)
            rejected,_=native('rejection9-unselected',1,cp=late,interval=1000,mutate=lambda t:set_parameter(t,'RANSWallPimpleHardAuditRejectStep',9),expect=False)
            require(same_history(startup,rejected),'late deferred-packet rollback changed histories')
            packet=Path(rejected['attempt_packet_status']['last_written_packet']);require(packet.exists() and load(packet)['time']>whole['time'],'unselected failed attempt not persisted')
            for label,mutation in [('policy',lambda t:set_info(t,'RANSWallBoundaryPolicy','CompatibilityV1')),
                ('outlet-datum',lambda t:set_parameter(t,'RANSWallOutletTotalPressure',.01)),
                ('dt',lambda t:set_parameter(t,'TimeStep',2e-5))]:
                q,r=native('reject-'+label,0,cp=cp,mutate=mutation,expect=False)
                require(q and 'signature mismatch' in q['failure'],'incompatible boundary/dt did not reach checkpoint guard')
            summary['observations']['open8-on']=observe('open8-on',whole)
            summary['tests'].update(fresh_origin=True,backflow_outflow_fixture=True,open_cache_bitwise=True,open_baseline=ru,open_optimized=rw,
                restart_bitwise=True,rollback_bitwise=True,deferred_rejected_packet=str(packet),boundary_interpretation_restart_rejected=True)
            save(status='passed',completed_UTC=datetime.now(timezone.utc).isoformat())
        else:
            if a.reuse_temporal_root:
                previous=load(a.reuse_temporal_root/'progress.json')
                require(previous['status']=='validation_or_numerical_stop' and previous.get('failure')=="'weights_0'",
                    'only the explicit missing-h0-sample reporting stop can reuse completed legs')
                require(previous['runtime']['sha256']==image['sha256'] and previous['profile_sha256']==summary['profile_sha256'] and
                    Path(previous['h1']).resolve()==a.h1.resolve(),'reused temporal physics/runtime differs')
                summary['reused_temporal_root']=str(a.reuse_temporal_root.resolve())
                summary['reused_temporal_progress_sha256']=sha(a.reuse_temporal_root/'progress.json')
            policy=dict(version='H2aStartupTemporalV1',steps='IMEX1 sequential candidate, expected smooth order1',
                dt=[4e-5,2e-5,1e-5],conditional_extra_dt=5e-6,conditional_extra_reason='only if finest primary level fails measured accuracy/sensitivity screen',common_times=[.0002,.0004],milestones_tstar=[.02,.1,1.],
                normalization='physical area RMS: U_H=1, K=completed K(Hw), epsilon=E_C(Hw); absolute pressure datum, no mean subtraction',
                scaled_error_targets=dict(velocity=1e-4,k=1e-4,epsilon=1e-3,pressure=1e-3,mean_traction=1e-4),
                response_relative_target=.01,numerical_floor=1e-12,coherent_order_band=[.8,1.2],
                primary_wall_cap_seconds=a.primary_cap_seconds,rule='choose largest dt satisfying absolute and response-relative targets at both common times; coherent Richardson estimates, otherwise conservative 4x adjacent sensitivity for coarse/middle and 2x for finest; no stability-only selection')
            atomic(a.output/'temporal-policy.json',policy);summary['policy']=policy;save()
            def one_dt(dt):
                label={4e-5:'coarse',2e-5:'middle',1e-5:'fine',5e-6:'extra-fine'}[dt]
                reuse=a.reuse_temporal_root and (a.reuse_temporal_root/(label+'-h0/receipt.json')).exists()
                def retained(name,target):
                    folder=a.reuse_temporal_root/name;receipt=load(folder/'receipt.json');q=load(folder/'case.evolution.json')
                    require(receipt['code']==0 and not q['failure'] and receipt['report_sha256']==sha(folder/'case.evolution.json') and
                        receipt['xml_sha256']==sha(folder/'case.xml'),'retained temporal receipt/decision differs')
                    require(q['dt']==dt and q['steps']==target and abs(q['time']-target*dt)<1e-14,'retained temporal time map differs')
                    saved=checkpoint(Path(q['checkpoint']))
                    require(state_hash(saved['coupled'])==state_hash(q['state']) and saved['time']==q['time'],
                        'retained native checkpoint differs from report')
                    return q,receipt
                seed,rs=retained(label+'-h0',0) if reuse else native(label+'-h0',0,dt)
                # Zero-step native reports intentionally omit final sample arrays.
                # Strict readback evaluates the ACTUAL saved h0, without Lift or
                # history advance, using precisely the same physical norm rows.
                observation=observe(label+'-h0',seed,source=a.reuse_temporal_root if reuse else None)
                seed['physical_samples']=load(Path(observation['path']))['diagnostic']['physical_samples']
                cp=Path(seed['checkpoint']);out={0.:seed};receipts=[]
                for T in policy['common_times']:
                    current=out[max(out)];target=round(T/dt)
                    if reuse:
                        current,r=retained(label+'-step'+str(target),target);receipts.append(r);out[T]=current;continue
                    while current['steps']<target:
                        end=min(target,round(current['steps'])+32)
                        current,r=native(label+'-step'+str(end),end-round(current['steps']),dt,cp,interval=end if end==target else 32)
                        cp=Path(current['checkpoint']);receipts.append(r)
                    out[T]=current
                return label,dt,out,receipts
            with ThreadPoolExecutor(max_workers=3) as pool:levels=list(pool.map(one_dt,policy['dt']))
            require(len({state_hash(v[2][0.]['state']) for v in levels})==1,'time levels changed physical native h0')
            K0=load(a.h1/'preflight.json')['candidates'][1]['K'];E0=load(a.h1/'preflight.json')['candidates'][1]['scan'][7]['epsilon_H']
            area=.4*.8;ea=.4*(.8-.015);normal=dict(U=UH*math.sqrt(area),W=UH*math.sqrt(area),velocity=UH*math.sqrt(area),
                k=K0*math.sqrt(area),epsilon=E0*math.sqrt(ea),pressure=UH**2*math.sqrt(area),
                mean_traction=USTAR**2,wall_stress=USTAR**2*math.sqrt(.4),wall_K=K0*math.sqrt(.4))
            def compare(triplet):
                comparisons={}
                for T in policy['common_times']:
                    qc,qm,qf=[v[2][T] for v in triplet];dc=physical(qc,qm);df=physical(qm,qf)
                    drift=[physical(v[2][T],v[2][0.]) for v in triplet];metrics={}
                    for k in normal:
                        c,f=dc[k],df[k];floor=normal[k]*policy['numerical_floor'];order=math.log2(c/f) if c>floor and f>floor else None
                        coherent=order is not None and .8<=order<=1.2
                        fine=f/(2**order-1) if coherent else None
                        metrics[k]=dict(coarse_middle=c,middle_fine=f,order=order,
                            fine_Richardson=fine,middle_Richardson=2**order*fine if coherent else None,
                            normalized_difference=f/normal[k],response_relative=f/max(drift[1][k],drift[2][k],floor))
                    comparisons[str(T)]=dict(fields=metrics,drift_from_native_h0=drift)
                return comparisons
            def choose(triplet,comparison):
                for idx in (0,1,2):
                    good=True
                    for T in policy['common_times']:
                        for k,target in policy['scaled_error_targets'].items():
                            m=comparison[str(T)]['fields'][k]
                            if m['fine_Richardson'] is not None:err=m['fine_Richardson']*2**(m['order']*(2-idx))
                            else:err=(4*m['coarse_middle'] if idx==0 else (4 if idx==1 else 2)*m['middle_fine'])
                            response=physical(triplet[idx][2][T],triplet[idx][2][0.])[k]
                            if err/normal[k]>target or (response>100*normal[k]*1e-12 and err/response>.01):good=False
                    if good:return idx
                return None
            comparisons=compare(levels);selected=choose(levels,comparisons)
            summary['temporal']=comparisons;summary['normalizations']=normal;save()
            if selected is None:
                extra=one_dt(policy['conditional_extra_dt']);levels.append(extra)
                require(state_hash(extra[2][0.]['state'])==state_hash(levels[0][2][0.]['state']),'extra level changed physical h0')
                refined=compare(levels[1:]);local=choose(levels[1:],refined)
                summary['conditional_temporal']=refined;save()
                if local is None:
                    q=extra[2][.0004];name=extra[0]+'-step'+str(round(.0004/extra[1]))
                    summary['observations']['finest-bounded']=observe(name,q)
                    save(status='temporal_accuracy_gate_not_met',accepted_steps=q['steps'],accepted_time=q['time'],checkpoint=q['checkpoint'],
                        interpretation='bounded startup triplets executed; no longer trajectory selected on survival alone',completed_UTC=datetime.now(timezone.utc).isoformat())
                    return
                selected=local+1
            summary['selected_level']=levels[selected][0];dt=levels[selected][1];q=levels[selected][2][.0004]
            summary['selected_dt']=dt;summary['completed_tstar']=q['time']*UH/HH
            source_name=levels[selected][0]+'-step'+str(round(.0004/dt))
            source=a.reuse_temporal_root if a.reuse_temporal_root and (a.reuse_temporal_root/source_name).exists() else None
            summary['observations']['temporal-selected']=observe(source_name,q,source=source)
            save(status='primary_startup_running',accepted_steps=q['steps'],accepted_time=q['time'],checkpoint=q['checkpoint'])
            start=time.monotonic();mean_cost=sum(r['seconds'] for r in levels[selected][3])/sum(r['local_steps_requested'] for r in levels[selected][3])
            cumulative=0.;positive=0.;negative=0.;absolute=0.;primary=[]
            for receipt in levels[selected][3]:
                for row in load(Path(receipt['report']))['pimple_attempts']:
                    bound=row['k_bound'];terms=bound['integrals_inner_outer_near_positive_negative']
                    cumulative+=bound['integral_total'];positive+=terms[3];negative+=terms[4];absolute+=abs(terms[3])+abs(terms[4])
            for tstar in (.02,.1,1.):
                target=round(tstar*HH/UH/dt);remaining=target-round(q['steps']);forecast=remaining*mean_cost
                if time.monotonic()-start+forecast>a.primary_cap_seconds:
                    save(status='bounded_cost_stop',next_tstar=tstar,next_target_steps=target,forecast_seconds=forecast,
                        elapsed_primary_seconds=time.monotonic()-start);break
                while q['steps']<target:
                    end=min(target,round(q['steps'])+32);name='primary-step'+str(end)
                    if time.monotonic()-start+(end-q['steps'])*mean_cost>a.primary_cap_seconds:
                        save(status='bounded_cost_stop',next_tstar=tstar,next_target_steps=target,accepted_steps=q['steps'],accepted_time=q['time'],
                            checkpoint=q['checkpoint'],elapsed_primary_seconds=time.monotonic()-start,completed_UTC=datetime.now(timezone.utc).isoformat());return
                    q,r=native(name,end-round(q['steps']),dt,Path(q['checkpoint']),interval=end if end==target else 32)
                    for z in q['pimple_attempts']:
                        bound=z['k_bound'];terms=bound['integrals_inner_outer_near_positive_negative']
                        cumulative+=bound['integral_total'];positive+=terms[3];negative+=terms[4];absolute+=abs(terms[3])+abs(terms[4])
                    primary.append(r)
                    mean_cost=.5*mean_cost+.5*r['seconds']/len(q['pimple_attempts'])
                    save(accepted_steps=q['steps'],accepted_time=q['time'],checkpoint=q['checkpoint'],checkpoint_sha256=sha(Path(q['checkpoint'])/'State.bin'),
                        cumulative_bound_signed_integral=cumulative,cumulative_bound_positive_integral=positive,cumulative_bound_negative_integral=negative,cumulative_bound_absolute_integral=absolute,primary_receipts=primary,measured_seconds_per_step=mean_cost)
                summary['observations']['tstar-'+str(tstar)]=observe(name,q)
                save(completed_tstar=tstar,accepted_steps=q['steps'],accepted_time=q['time'])
            else:save(status='requested_startup_targets_completed')
            save(completed_UTC=datetime.now(timezone.utc).isoformat(),elapsed_primary_seconds=time.monotonic()-start)
    except subprocess.TimeoutExpired as error:save(status='infrastructure_timeout',failure=str(error));raise
    except Exception as error:save(status='validation_or_numerical_stop',failure=str(error));raise
    print(a.output/'progress.json')


if __name__=='__main__':main()
