#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""CT7's separately contracted finer CT6 comparison and stage diagnosis.

Absolute pressure remains the admission observable. Decompositions and
first-order extrapolation are read-only and never replace the failed CT6 gate.
"""
import argparse,json,math,os
from functools import lru_cache
from pathlib import Path
from rans_qualification import THREADS,write_json,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from ct_accepted6_compare import observations,norms,NAMES
from ct_accepted6_audit import common_grid
from ct_accepted6_operator import geometry
from ct_accepted6_metrics import wall_nut
from ct_accepted2_audit import values
from ct_accepted_feasibility import parameters
from ct_accepted_open import require
from ct_accepted_open import row_value
from ct_accepted6_continuous import Continuous
from validate_a93_a95_e94e95 import oracle
from rans_request_manifest import RequestManifest,sha


def pressure_parts(record):
    _,_,_,v,truth,w,region,_=record;e=v[:,4]-truth[:,4];mean=float(w@e/w.sum())
    regions={str(r):float(w[region==r]@e[region==r]/w[region==r].sum()) for r in (0,1)}
    return dict(absolute_RMS=float(np.sqrt(w@(e*e)/w.sum())),common_offset=mean,
        variation_RMS=float(np.sqrt(w@((e-mean)**2)/w.sum())),region_mean_error=regions,
        difference_of_region_mean_errors=regions['1']-regions['0'],
        pressure_scale=.1,units='native nondimensional Pi; admission uses absolute RMS divided by .1',
        scope='diagnosis only; no gauge removal or retiming in pressure admission')


def stage_metrics(packet,state,n):
    o=values(packet,state,len(packet['test_basis']))['observation'];v=np.asarray(o['grid_values']);g=np.asarray(o['grid_gradients']);p=parameters('A95')
    tractions={};lag=[];mechanical={};pressure={};minimum=1e300
    for i in range(n,len(v)):
        q=packet['observation']['grid'][i][0];side=packet['observation']['grid_sides'][i];region=side>0
        G=np.zeros((3,3));G[0]=g[i,0];G[2]=g[i,1];k=v[i,2];pi=v[i,4]
        pressure[str((q,side))]=float(pi);mechanical[str((q,side))]=float(pi-2*k/3)
        try:
            if region:
                law=oracle('A95',k,v[i,3],G);R=law['R'];T=-R+np.eye(3)*(2*k/3-pi)
                lag.append(float(v[i,3]-(k/p['q'])**1.5/(p['kappa']*.5)))
            else:
                nt=wall_nut('A95',k,.5);R=2*k/3*np.eye(3)-nt*(G+G.T)
                T=(nt+1.5e-5)*(G+G.T)-pi*np.eye(3)
            *_,normal=geometry(q,.5);tractions[str((q,side))]=(T@normal).tolist()
            minimum=min(minimum,float(np.linalg.eigvalsh(R)[0]))
        except (ValueError,FloatingPointError) as error:tractions[str((q,side))]=dict(unavailable=str(error),diagnostic_only=True)
    return dict(pressure=pressure,mechanical_pressure=mechanical,full_interface_traction=tractions,
        refreshed_dependent_epsilon_sample_maximum=max(map(abs,lag),default=None),sampled_covariance_minimum=minimum,
        interface_inventory='13 fixed chart stations with explicit sides, sample diagnostics; not a surface integral')


def diagnose(member,audit):
    receipt=json.loads((member/'results.json').read_text());report=json.loads((member/'trajectory/case.evolution.json').read_text())
    a=json.loads((audit/'results.json').read_text());path=Path(a['common_observations']['path'])
    require(sha(path)==a['common_observations']['sha256'],'changed stage diagnostic basis packet')
    packet=json.loads(path.read_text());_,_,_,n=common_grid();off=packet['observation']['offsets'];w=np.asarray(common_grid()[2])
    initial=values(packet,report['initial'],n)['observation'];g=np.asarray(initial['grid_gradients']);div=g[:,0,0]+g[:,1,2]
    initial_path=Path(receipt['initial_packet']['path'])
    require(sha(initial_path)==receipt['initial_packet']['sha256'],'changed initial native rows')
    native=json.loads(initial_path.read_text())['observation'];initial_load=np.zeros(off[6]-off[4])
    for point in native['volume']:
        velocity=point['velocity_k'];state=np.asarray(report['initial']);gradient=np.zeros((3,3))
        gradient[0]=row_value(velocity,state,off[0])[1];gradient[2]=row_value(velocity,state,off[1])[1]
        pressure=point['pressure'];block=4+point['region']
        initial_load[np.asarray(pressure['ids'],int)+off[block]-off[4]]+=point['weight']*np.asarray(pressure['phi'])*np.trace(gradient)
    continuum=Continuous('A95')
    @lru_cache(None)
    def jump(q,time):return float(continuum.ray(q,time)['jump'])
    def time_errors(state,time,dt):
        represented=np.asarray(values(packet,state,n)['observation']['grid_values'])[:,4];out={}
        for label,t in [('target',time),('previous',time-dt),('midpoint',time-dt/2)]:
            truth=np.array([.03*(point[0]+.2*point[2])+.02*t+.01*t*t+(jump(point[0],t) if region else 0.)
                for point,region in zip(packet['observation']['grid'][:n],packet['observation']['grid_regions'][:n],strict=True)])
            error=represented-truth;mean=float(w@error/w.sum())
            out[label]=dict(physical_time=t,absolute_scaled_RMS=float(np.sqrt(w@(error*error)/w.sum())/.1),
                common_offset=mean,spatial_variation_RMS=float(np.sqrt(w@((error-mean)**2)/w.sum())))
        return out
    result=dict(initial_divergence_RMS=float(np.sqrt(w@(div*div)/w.sum())),initial_pressure_error=receipt['initial_metrics']['fields']['Pi'],
        initial_native_weak_continuity=dict(maximum_absolute_load=float(max(abs(initial_load))),
            inner_maximum=float(max(abs(initial_load[:off[5]-off[4]]))),outer_maximum=float(max(abs(initial_load[off[5]-off[4]:]))),
            rows=len(initial_load),scope='all actual free pressure tests, physical gradients and native volume quadrature; distinct from pointwise divergence'),stages=[])
    attempts=report['pimple_attempts']
    for index in sorted(set([0,len(attempts)-1])):
        attempt=attempts[index];fp=Path(attempt['flow_snapshot']);kp=Path(attempt['k_snapshot'])
        require(sha(fp)==attempt['flow_snapshot_sha256'] and sha(kp)==attempt['k_snapshot_sha256'],'changed computed stage snapshots')
        f=json.loads(fp.read_text());k=json.loads(kp.read_text());previous=report['accepted_states'][index]
        phases={'previous_accepted':previous,'current_time_lifted_input':f['input'],'momentum_predictor':f['momentum_predictor'],
                'pressure_corrected_flow':f['corrected_flow'],'post_epsilon':k['pre_k'],'post_k':k['corrected']}
        row=dict(time=f['time'],dt=1/f['alpha'],flow_snapshot_sha256=sha(fp),k_snapshot_sha256=sha(kp),
            phase_state_sha256={name:exact_vector_hash(s) for name,s in phases.items()},metrics={name:stage_metrics(packet,s,n) for name,s in phases.items()},
            scalar_stage_pressure_coefficient_change=float(max(abs(np.asarray(k['corrected'])[off[4]:off[6]]-np.asarray(f['corrected_flow'])[off[4]:off[6]]))),
            current_lift_minus_previous_by_block=[float(max(abs(np.asarray(f['input'])[off[b]:off[b+1]]-np.asarray(previous)[off[b]:off[b+1]]))) for b in range(6)],
            predictor_linear_residual=attempt['predictor_linear_residual'],pressure_linear_residual=attempt['pressure_linear_residual'],
            momentum_splitting_identity_defect=attempt['momentum_splitting_identity_defect'],native_full_pressure_action_defect=attempt['native_full_G_action_defect'],
            frozen_momentum_defect=attempt['frozen_momentum_defect'],refreshed_trace_defect=attempt['refreshed_trace_defect'])
        row['pressure_time_comparisons']={name:time_errors(phases[name],f['time'],1/f['alpha'])
            for name in ('previous_accepted','current_time_lifted_input','pressure_corrected_flow','post_epsilon','post_k')}
        result['stages'].append(row)
    result['time_semantics']='stored Pi is updated by dz in the flow solve at declared target time with target-time forcing/essential lift; epsilon and k then change constitutive normal stress without updating Pi. This is frozen-stage splitting, not permission to relabel pressure time or repair it.'
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('old','fine','audit','output'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--fine-p3',type=Path);p.add_argument('--audit-p3',type=Path)
    a=p.parse_args();old=a.old
    members={.002:old/'evolution-v2/A95-p4-dt002',.001:old/'evolution-v3/A95-p4-dt001',.0005:a.fine}
    audits={.002:old/'physical-audits-v3/A95-p4-dt002',.001:old/'physical-audits-v3/A95-p4-dt001',.0005:a.audit}
    m=RequestManifest(a.output,dict(schema='CT7FinerTemporalDiagnosisV1',source=sha(Path(__file__)),
        member_receipts={str(k):sha(v/'results.json') for k,v in members.items()},audit_receipts={str(k):sha(v/'results.json') for k,v in audits.items()},
        old_CT6_comparison=sha(old/'primary-comparisons-v2/results.json'),fine_p3=sha(a.fine_p3/'results.json') if a.fine_p3 else None,
        audit_p3=sha(a.audit_p3/'results.json') if a.audit_p3 else None,absolute_pressure=True))
    data={dt:observations(member,audits[dt]) for dt,member in members.items()}
    require(len({r[0]['zero']['state_sha256'] for r in data.values()})==1,'fresh same-space initial coefficients differ')
    require(all(r[7]==data[.001][7] for r in data.values()),'foreign continuous inputs')
    D0=norms(data[.002][3]-data[.001][3],data[.001][5],data[.001][6],'A95')
    D1=norms(data[.001][3]-data[.0005][3],data[.001][5],data[.001][6],'A95');orders={}
    for name in NAMES:
        x=D0[name]['scaled_weighted_RMS'];y=D1[name]['scaled_weighted_RMS'];resolved=min(x,y)>1e-6
        order=math.log2(x/y) if resolved else None
        orders[name]=dict(D002_001=x,D001_0005=y,order=order,signal_floor=1e-6,band=[.8,1.2],status=('pass' if .8<=order<=1.2 else 'fail') if resolved else 'unresolved_floor')
    temporal_pass=all(r['status'] in ('pass','unresolved_floor') for r in orders.values()) and orders['Pi']['status']=='pass'
    result=dict(schema='CT7FinerTemporalDiagnosisV1',orders=orders,temporal_status='pass' if temporal_pass else 'failed_or_unresolved',
        pressure_decomposition={str(dt):pressure_parts(r) for dt,r in data.items()},diagnosis={str(dt):diagnose(members[dt],audits[dt]) for dt in members},
        actual_initial_coefficients_identical=True,history_policy='independently initialized zero-step native history per dt, unchanged method',
        matched_p3_eligible=temporal_pass,old_CT6_decision='unchanged failed original triplet/degree comparison',curved_model_transfer_eligible=False)
    if a.fine_p3:
        require(a.audit_p3 is not None and temporal_pass,'unadmitted matched P3 comparison')
        coarse=observations(a.fine_p3,a.audit_p3);fine=data[.0005]
        ce=norms(coarse[3]-coarse[4],coarse[5],coarse[6],'A95');fe=norms(fine[3]-fine[4],fine[5],fine[6],'A95');delta=norms(coarse[3]-fine[3],fine[5],fine[6],'A95')
        oldp3=observations(old/'evolution-v3/A95-p3-dt001',old/'physical-audits-v3/A95-p3-dt001')
        require(oldp3[0]['zero']['state_sha256']==coarse[0]['zero']['state_sha256'],'matched P3 initialization differs')
        spread3=norms(oldp3[3]-coarse[3],coarse[5],coarse[6],'A95');rows={}
        for name in NAMES:
            x=ce[name]['scaled_weighted_RMS'];y=fe[name]['scaled_weighted_RMS'];h=delta[name]['scaled_weighted_RMS'];spread=max(spread3[name]['scaled_weighted_RMS'],D1[name]['scaled_weighted_RMS'])
            margin=.85*x-y;status='unresolved_floor' if max(x,y)<=1e-7 else 'pass' if margin>2e-8 else 'fail' if margin < -2e-8 else 'unresolved_uncertainty'
            rows[name]=dict(coarse_error=x,fine_error=y,ratio=y/x if x else None,reduction_status=status,spatial_difference=h,
                timestep_spread=spread,timestep_limit=max(1e-7,.1*h),uncertainty_status='pass' if spread<=max(1e-7,.1*h) else 'unresolved_temporal_contamination')
        result['matched_degree']=rows
        result['spatial_temporal_candidate_checks_pass']=all(r['reduction_status'] in ('pass','unresolved_floor') and r['uncertainty_status']=='pass' for r in rows.values())
        result['transfer_requirement']='independent complete A95 primitive admission, including actual-stage budgets, frozen quadrature, all shared prerequisites and closed resources; comparison success alone never admits transfer'
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);print('CT7 new temporal gate',result['temporal_status'],flush=True)
    return 0 if temporal_pass else 1

if __name__=='__main__':raise SystemExit(main())
