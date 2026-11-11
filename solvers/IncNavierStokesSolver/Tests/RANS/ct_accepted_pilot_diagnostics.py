#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only failed-pilot diagnosis plus an identical minimal timed reproducer."""
import argparse,json,os,shutil,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS,write_json,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from ct_accepted_open import probe,run,require,row_value
from ct_accepted_pilot import replay,grid,physical_observations,constitutive_comparison,compare,budget_classification
from run_cls_epawf_picard import set_parameter,set_info,write
from compare_cls_epawf_pimple import same_history
from rans_request_manifest import RequestManifest,runtime_identity,sha


def radicand_witness(initial,report):
    o=initial['observation'];s=np.array(report['rejected_state']);off=np.array(o['offsets'],int);minimum=1.;witness=None;count=0
    for q in o['volume']:
        if not q['region']:continue
        v=q['velocity_k'];G=np.zeros((3,3));G[0]=row_value(v,s,off[0])[1];G[2]=row_value(v,s,off[1])[1]
        k=row_value(v,s,off[2])[0];epsilon=row_value(q['epsilon'],s,off[3])[0]
        S=(G+G.T)/2-np.eye(3)*np.trace(G)/3;O=(G-G.T)/2;sn=np.linalg.norm(S);on=np.linalg.norm(O);tau=k/epsilon
        require(k>0 and epsilon>0,'rejected_volume_scalar_invalid_before_radicand')
        if sn==0:
            require(on==0,'rejected_pure_rotation_domain');C=1/6.5
        else:C=1/(6.5+3*np.linalg.eigvalsh(S)[-1]/sn*np.hypot(sn,on)*tau)
        rad=float(1-9*C*C*sn*sn*tau*tau);count+=1
        if rad<minimum:
            minimum=rad;witness=dict(coordinate=q['coordinate'],geometry=v['geometry'],side=0,region='outer',
                generating_geometry=q['generating_geometry'],raw_gradient=G.tolist(),k=k,epsilon=epsilon,c_mu=float(C),
                physical_radicand=rad,raw_divergence=float(np.trace(G)),strain_norm=float(sn),rotation_norm=float(on),
                isotropic_production_work=float(-2*k*np.trace(G)/3),production='undefined: actual A95 state outside radicand domain',
                native_velocity_row=v,native_epsilon_row=q['epsilon'])
    require(minimum<0 and witness is not None,'native_failure_not_reproduced_on_actual_outer_volume_rows')
    return dict(status='confirmed_constitutive_failure',minimum=minimum,points=count,witness=witness,
        stage='after momentum predictor and consistent-mass pressure correction, before epsilon transport',
        attempted_time=report['pimple_attempts'][-1]['time'],last_accepted_time=report['time'],
        basis='sealed native initial packet; outer quadrature rows fixed by this mesh/order, unchanged during evolution',
        rejected_state_sha256=exact_vector_hash(s),native_failure=report['failure'])


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--runtime',type=Path,required=True);p.add_argument('--pilot',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
    manifest=RequestManifest(a.output,dict(schema='CTAcceptedPilotPostmortemV1',pilot_results_sha256=sha(a.pilot/'results.json'),
        runtime={n:runtime_identity(a.runtime/n) for n in ('IncNavierStokesSolver','RANSPreservationProbe')},
        reports={n:sha(a.pilot/n/'case.evolution.json') for n in ('A95-dt-0.002','A95-dt-0.001')},
        common_prefix_time=.002,reproducer='fresh unchanged dt=.002, two attempts only; no continuation or source/field modification',
        required_pilot_status='failed remains failed'))
    result=dict(status='fail',pilot_qualification='failed',reconstructions={},witnesses={},executions=[]);endpoints={}
    try:
        for name in ('A95-dt-0.002','A95-dt-0.001'):
            source=a.pilot/name;report=json.loads((source/'case.evolution.json').read_text());initial=json.loads((a.pilot/(name+'-zero')/'initial.json').read_text())
            require('negative square-root radicand' in report['failure'],'foreign_pilot_failure')
            result['witnesses'][name]=radicand_witness(initial,report)
            c=json.loads((source/'construction.json').read_text());points,weights=grid(c)
            times=report['audit_times'];states=report['accepted_states'];require(len(times)==len(states),'accepted_time_state_coverage')
            index=next(i for i,t in enumerate(times) if abs(t-.002)<1e-13);state=states[index]
            row=dict(time=times[index],state_sha256=exact_vector_hash(state),qualification='diagnostic common prefix, intended endpoint not reached',budgets={});result['reconstructions'][name]=row
            folder=a.output/name;folder.mkdir();shutil.copy2(source/'case.xml',folder/'case.xml')
            packet,e=probe(a.runtime,folder,replay(initial,state,times[index],points,export_basis=True,export_constitutive=True),'prefix');result['executions'].append(e)
            if e['exit_code']==0:
                row['observations']=physical_observations(packet,weights,c);row['constitutive']=constitutive_comparison(packet,c);endpoints[name]=packet
            else:row['observation_failure']=packet
            attempt=next(x for x in report['pimple_attempts'] if abs(x['time']-.002)<1e-13 and x['scheduled_complete'])
            snapshot=json.loads(Path(attempt['k_snapshot']).read_text())
            for label,orders,face in [('native',[32,8,10,10,10],16),('q48f16',[48,12,12,12,12],16),('q64f16',[64,16,16,16,16],16),('q64f32',[64,16,16,16,16],32)]:
                sub=folder/label;sub.mkdir();root=ET.parse(source/'case.xml').getroot()
                for k,v in zip(('RANSWallNormalOrder','RANSWallTangentialOrder','RANSWallMassOrder','RANSWallVolumeOrder','RANSWallInterfaceOrder'),orders,strict=True):set_parameter(root,k,v)
                write(root,sub/'case.xml')
                request=replay(initial,state,times[index],[],export_basis=False,export_physical_scalar_budget=True,
                    previous_state=states[index-1],stage_history=snapshot['history'],export_frozen_scalar_reconciliation=True,physical_face_order=face)
                b,e=probe(a.runtime,sub,request,'budget');result['executions'].append(e)
                row['budgets'][label]=budget_classification(b,c) if e['exit_code']==0 else dict(status='failed_forward_or_audit',packet=b)
            write_json(a.output/'progress.json',result)
        if len(endpoints)==2:result['common_prefix_timestep_comparison']=compare(endpoints['A95-dt-0.002'],endpoints['A95-dt-0.001'],weights,c)
        folder=a.output/'minimal-timed-reproducer';folder.mkdir();source=a.pilot/'A95-dt-0.002';root=ET.parse(source/'case.xml').getroot()
        set_parameter(root,'NumSteps',2);set_info(root,'RANSWallRestart','None');set_info(root,'RANSWallCheckpoint',str(folder/'checkpoint'));write(root,folder/'case.xml')
        os.environ['NEKTAR_RANS_TIMING_FILE']=str(folder/'timing.jsonl')
        e=run([a.runtime/'IncNavierStokesSolver','case.xml'],folder,'solver.log');result['executions'].append(e);os.environ.pop('NEKTAR_RANS_TIMING_FILE')
        r=json.loads((folder/'case.evolution.json').read_text());old=json.loads((source/'case.evolution.json').read_text())
        require(e['exit_code']!=0 and r['failure']==old['failure'] and same_history(r,old) and r['rejected_state']==old['rejected_state'],'minimal_or_telemetry_reproducer_changed_failure')
        timing={}
        for line in (folder/'timing.jsonl').read_text().splitlines():
            t=json.loads(line);row=timing.setdefault(t['name'],dict(exclusive_process_cpu_seconds=0.,exclusive_wall_seconds=0.,calls=0))
            row['exclusive_process_cpu_seconds']+=t['exclusive_process_cpu_seconds'];row['exclusive_wall_seconds']+=t['exclusive_wall_seconds'];row['calls']+=1
        result['minimal_reproducer']=dict(status='pass',actual_failure_unchanged=True,state_and_native_history_bitwise_equal=True,
            rejected_state_bitwise_equal=True,configuration_sha256=sha(folder/'case.xml'),timing_sha256=sha(folder/'timing.jsonl'),timing=timing,
            source='fresh native history; unchanged physical functions, mesh, coefficients and dt; two attempts; telemetry only')
        result['status']='pass'
    except Exception as error:result['diagnostic_failure']=repr(error)
    write_json(a.output/'results.json',result);manifest.seal(a.output,['results.json']+[str(p.relative_to(a.output)) for p in a.output.rglob('case.xml')])
    print(result['status'],result.get('diagnostic_failure',''));return 0 if result['status']=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
