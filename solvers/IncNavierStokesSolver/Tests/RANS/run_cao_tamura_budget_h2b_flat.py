#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H2b transaction faults and bounded finite/developing post-startup IVPs."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import math,os,platform,subprocess,time
from pathlib import Path
from datetime import datetime,timezone
import xml.etree.ElementTree as ET
from continue_cls_epawf_pimple import atomic,freeze,CampaignLock
from run_cao_tamura_budget_h2a import prepare,physical
from run_cls_epawf_picard import THREADS,load,require,sha,state_hash,set_info,set_parameter,write
from diagnose_cls_accepted_covariance import checkpoint,common_packet,run_probe


def identical(a,b):
    return all(state_hash(a[k])==state_hash(b[k]) for k in ('state','fixed_reference','fixed_defect','native_history_times')) and \
        [state_hash(v) for v in a['native_history']]==[state_hash(v) for v in b['native_history']]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('solver','probe','h1','parent','output'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--validation-only',action='store_true',help='focused final-runtime import/transaction checks; reuse the completed temporal triplet')
    p.add_argument('--cadence-only',action='store_true',help='three-step nonzero-epoch first/periodic/unselected packet check')
    a=p.parse_args();require(platform.system()=='Linux','Linux native qualification required')
    require(sha(a.parent/'State.bin')=='7485f1b78243f6cb09d7406f3c685caa8779cf868c5826f9d8f9aca025dd6605','wrong retained H2a parent')
    a.output.mkdir(parents=True,exist_ok=False);lock=CampaignLock(a.output);os.environ.update(THREADS)
    solver,image=freeze(a.solver.resolve(),a.output/'runtime');probe,probe_image=freeze(a.probe.resolve(),a.output/'probe-runtime')
    parent=checkpoint(a.parent);summary=dict(status='started',source='H2b-flat-v1',runtime=image,probe=probe_image,
        parent=str(a.parent),parent_sha256=sha(a.parent/'State.bin'),parent_dt=parent['dt'],physical_t0=parent['time'],tests={},temporal={})
    def save(**kw):summary.update(kw);atomic(a.output/'progress.json',summary)
    def run(name,steps,dt=5e-6,restart=None,imported=False,mutate=None,success=True):
        xml=prepare(a.output/name,a.h1,steps,dt,restart,packets='SelectedV1',interval=1000)
        t=ET.parse(xml).getroot()
        if imported:
            for k,v in dict(RANSWallDiagnosticIVP='FiniteDevelopingIMEX1V1',RANSWallIVPParent=str(a.parent),
                RANSWallIVPParentSHA256=sha(a.parent/'State.bin'),RANSWallIVPParentTimeStep=format(parent['dt'],'.17g')).items():set_info(t,k,v)
        if mutate:mutate(t)
        write(t,xml)
        definition=load(xml.parent/'definition.json')
        start=checkpoint(restart)['time'] if restart else (parent['time'] if imported else 0.)
        definition.update(version='CaoTamuraBudgetH2bFiniteLocalIVPV1',xml_sha256=sha(xml),
            physical_t0=start,requested_local_elapsed=steps*dt,endpoint_time=start+steps*dt,
            diagnostic_import=imported,parent_checkpoint=str(a.parent),parent_sha256=summary['parent_sha256'],
            diagnostic_mutation=bool(mutate),scope='finite/developing compatibility or local post-startup IVP, not a developed precursor')
        atomic(xml.parent/'definition.json',definition)
        env=dict(os.environ,LD_LIBRARY_PATH=str(solver.parent/'lib'),TMPDIR=str(xml.parent),NEKTAR_RANS_TIMING_FILE=str(xml.parent/'timing.jsonl'))
        cmd=['/usr/bin/time','-f','%e %U %S %M','-o',str(xml.parent/'usage.txt'),str(solver),str(xml)]
        start=time.monotonic()
        with (xml.parent/'solver.log').open('w') as log:
            child=subprocess.Popen(cmd,cwd=xml.parent,env=env,stdout=log,stderr=subprocess.STDOUT,pass_fds=(lock.stream.fileno(),));code=child.wait(timeout=1800)
        report=xml.with_suffix('.evolution.json');q=load(report);receipt=dict(command=cmd,code=code,seconds=time.monotonic()-start,
            report=str(report),report_sha256=sha(report),xml_sha256=sha(xml),usage=(xml.parent/'usage.txt').read_text())
        atomic(xml.parent/'receipt.json',receipt)
        require((code==0 and not q['failure'])==success,'unexpected native outcome '+name+': '+q['failure'])
        return q,receipt
    def observe(name,q):
        cp=Path(q['checkpoint']);saved=checkpoint(cp);file=a.output/name/'case.evolution.json'
        packet=dict(common_packet(cp),operation='finite_inflow_h1',state=saved['coupled'],history=saved['coupled'],time=saved['time'],alpha=0.,
            inner_stress_model='Linear',native_report_file=str(file),native_report_sha256=sha(file),observation_x=0.)
        d,rr=run_probe(probe,file.parent/'case.xml',packet,a.output/(name+'-readback'));require(d['initialization_work']==[0,0] and d['governing_solves']==d['physical_commits']==0,'readback lifted/advanced')
        return d,rr
    try:
        if a.cadence_only:
            q,rr=run('epoch-periodic3',3,imported=True,mutate=lambda t:set_parameter(t,'RANSWallPacketInterval',2))
            selected=[v['attempt_packet_selected'] for v in q['pimple_attempts']]
            require(selected==[True,True,False] and q['attempt_packet_status']['native_attempt']==3,
                'packet cadence used absolute time instead of native epoch count')
            require(all(v['bound_calls']==v['momentum_predictor_solves']==v['pressure_system_solves']==v['epsilon_transport_solves']==v['k_transport_solves']==1
                for v in q['pimple_attempts']),'cadence fixture changed one-sweep schedule')
            save(status='cadence_passed',tests=dict(nonzero_epoch_first_periodic_unselected=True),
                selected=selected,receipt=rr,physical_steps=3)
            return
        zero,_=run('ordinary-zero',0,restart=a.parent);seed,_=run('import-zero',0,imported=True)
        require(identical(zero,seed) and seed['initialization_lift_calls_iterations']==[0,0],'finite import changed seed/reference/history bits')
        control,_=run('ordinary-one',1,restart=a.parent);ivp,_=run('import-one',1,imported=True)
        require(identical(control,ivp),'same-dt imported first step differs from native restart')
        for key in ('history_sha256','incoming_sha256','post_flow_sha256','pre_k_sha256','raw_k_sha256','corrected_k_sha256'):
            require(control['pimple_attempts'][0][key]==ivp['pimple_attempts'][0][key],'same-dt imported stage differs '+key)
        require(ivp['attempt_packet_status']['native_attempt']==1 and ivp['attempt_packet_status']['last_written_packet'], 'nonzero-time new epoch missed its first packet')
        summary['tests']['same_dt_import_state_reference_history_stages_bitwise']=True;save()
        def fault(t):set_info(t,'RANSWallPacketFault','OnFlush')
        def hard(t):fault(t);set_parameter(t,'RANSWallPimpleHardAuditRejectStep',parent['steps']+1)
        bad,_=run('reject-unselected-output-fault',1,restart=a.parent,mutate=hard,success=False)
        require(identical(bad,zero) and 'deliberate_PIMPLE_hard_audit_rejection' in bad['failure'] and bad['attempt_packet_status']['output_failure'] and not bad['attempt_packet_status']['pending_writer'],'packet failure masked rejection/rollback')
        def assigned(t):fault(t);set_parameter(t,'RANSWallPimplePostAssignRejectStep',parent['steps']+2)
        late,_=run('reject-assigned-output-fault',2,restart=a.parent,mutate=assigned,success=False)
        require(identical(late,control) and 'post-assignment explicit rejection' in late['failure'],'assigned candidate/history not rolled back')
        for k in ('last_stage_history','last_stage_time','last_stage_alpha','physical_samples','state_audits','stage_audits','explicit_audits','accepted_states','audit_times'):
            require(late[k]==control[k],'post-assignment rollback failed '+k)
        require(late['rollback_receipt']['assigned_candidate_sha256']!=late['rollback_receipt']['restored_state_sha256'] and
            late['attempt_packet_status']['output_failure'] and not late['attempt_packet_status']['pending_writer'],'post-assignment test did not exercise changed candidate/failed writer')
        for rejected in (bad,late):
            receipt=rejected['rollback_receipt'];attempt=rejected['pimple_attempts'][-1]
            require(state_hash(receipt['attempted_stage_history'])==receipt['attempted_stage_history_sha256']==attempt['history_sha256'] and
                receipt['attempted_stage_time']==attempt['time'] and receipt['attempted_stage_alpha']==attempt['alpha'],
                'packet-output failure lost actual rejected callback history/time/alpha')
        changed,_=run('ordinary-changed-dt-rejected',0,dt=1e-5,restart=a.parent,success=False)
        require('signature mismatch' in changed['failure'],'ordinary dt mismatch bypassed')
        summary['tests'].update(packet_fault_rollback=True,post_assignment_all_coupled_diagnostics_rollback=True,ordinary_changed_dt_rejected=True)
        def ordinary_reject(t):set_parameter(t,'RANSWallPimpleHardAuditRejectStep',parent['steps']+1)
        normal,_=run('normal-unselected-rejection',1,restart=a.parent,mutate=ordinary_reject,success=False)
        require(identical(normal,zero) and Path(normal['attempt_packet_status']['last_written_packet']).exists(),
            'normal unselected rejection lost packet or rollback')
        def selected(t):set_info(t,'RANSWallAttemptPackets','EveryAttempt')
        retained,_=run('normal-selected-acceptance',1,restart=a.parent,mutate=selected)
        require(identical(retained,control) and Path(retained['attempt_packet_status']['last_written_packet']).exists(),
            'selected packet changed accepted state/history or was lost')
        summary['tests']['normal_selected_and_unselected_packets']=True
        summary['faults']=dict(unselected=dict(reason=bad['failure'],packet=bad['attempt_packet_status'],rollback=bad['rollback_receipt']),
            assigned=dict(reason=late['failure'],packet=late['attempt_packet_status'],rollback=late['rollback_receipt']))
        if a.validation_only:
            require(sha(a.parent/'State.bin')==summary['parent_sha256'],'retained parent changed')
            save(status='validation_passed',completed_UTC=datetime.now(timezone.utc).isoformat())
            return
        policy=dict(dt=[4e-5,2e-5,1e-5],elapsed=[.0002,.0004],expected_order=1.,relative_response_target=.01,
            normalized_targets=dict(velocity=1e-4,k=1e-4,epsilon=1e-3,pressure=1e-3,mean_traction=1e-4),
            floor=1e-12,scope='same retained complete physical state; own IMEX1 histories, no accumulated T0 error claim, no automatic long continuation')
        summary['policy']=policy;save(status='temporal_running')
        def level(dt):
            label={4e-5:'coarse',2e-5:'middle',1e-5:'fine'}[dt];q,_=run(label+'-seed',0,dt,imported=True)
            d,rr=observe(label+'-seed',q);q['physical_samples']=d['physical_samples'];out={0.:q};receipts=[]
            current=None
            for elapsed in policy['elapsed']:
                count=round(elapsed/dt)-(round((out[max(out)]['time']-parent['time'])/dt) if elapsed else 0)
                name=label+'-elapsed'+str(elapsed);q,r=run(name,count,dt,restart=Path(current['checkpoint']) if current else None,imported=current is None)
                d,rr=observe(name,q);out[elapsed]=q;receipts.append(r)
                atomic(a.output/(name+'-observation.json'),dict(receipt=rr,diagnostic=d));current=q
            return label,dt,out,receipts
        with ThreadPoolExecutor(max_workers=3) as pool:levels=list(pool.map(level,policy['dt']))
        H=.015;K=.0250145796399802;E=.1050450587819008;area=.32;normal=dict(U=math.sqrt(area),W=math.sqrt(area),velocity=math.sqrt(area),
            k=K*math.sqrt(area),epsilon=E*math.sqrt(.4*(.8-H)),pressure=math.sqrt(area),mean_traction=.005515726603512427,
            wall_stress=.005515726603512427*math.sqrt(.4),wall_K=K*math.sqrt(.4))
        for elapsed in policy['elapsed']:
            dc=physical(levels[0][2][elapsed],levels[1][2][elapsed]);df=physical(levels[1][2][elapsed],levels[2][2][elapsed]);metrics={}
            for k in normal:
                c,f=dc[k],df[k];order=math.log2(c/f) if min(c,f)>normal[k]*1e-12 else None
                fine=f/(2**order-1) if order and .8<=order<=1.2 else None
                metrics[k]=dict(coarse_middle=c,middle_fine=f,order=order,fine_Richardson=fine,
                    normalized_difference=f/normal[k],drift=[physical(v[2][elapsed],v[2][0.])[k] for v in levels])
            summary['temporal'][str(elapsed)]=metrics
        qualifiers=[]
        for i,(label,dt,out,rr) in enumerate(levels):
            passes=True
            for metrics in summary['temporal'].values():
                for k,target in policy['normalized_targets'].items():
                    m=metrics[k];error=m['fine_Richardson']*2**(m['order']*(2-i)) if m['fine_Richardson'] is not None else (4*m['coarse_middle'] if i==0 else (4 if i==1 else 2)*m['middle_fine'])
                    response=m['drift'][i]
                    if error/normal[k]>target or (response>100e-12*normal[k] and error/response>.01):passes=False
            qualifiers.append(dict(level=label,dt=dt,passes=passes,receipts=rr,final_checkpoint=out[.0004]['checkpoint']))
        summary['qualification']=qualifiers;summary['normalizations']=normal
        require(sha(a.parent/'State.bin')==summary['parent_sha256'],'retained parent changed')
        save(status='completed',completed_UTC=datetime.now(timezone.utc).isoformat(),governing_steps_in_temporal_triplet=70)
    except Exception as e:save(status='failed',failure=str(e));raise
    print(a.output/'progress.json')

if __name__=='__main__':main()
