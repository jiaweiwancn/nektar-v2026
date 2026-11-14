#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only CT8 absolute-pressure uncertainty and stored-stage diagnosis.

No solve, field modification, pressure regauging or admission change is made.
Run under ct_accepted3_command.py with numerical threads set before import.
"""
import argparse,json,os
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from ct_accepted_feasibility import write_json,require
from ct_accepted_open import row_value
from ct_accepted8_native import col,reference
from ct_accepted5_traction import profiles
from rans_request_manifest import RequestManifest,sha

def read(p):return json.loads(Path(p).read_text())
def verified(p,h):
    p=Path(p);require(sha(p)==h,'changed source artifact:'+str(p));return read(p)

def decomposition(v,w,regions,scale):
    out={}
    for name,mask in [('whole',np.ones(len(w),dtype=bool)),('inner',regions==0),('outer',regions==1)]:
        a=v[mask];weights=w[mask];mean=float(weights@a/weights.sum());rms=float(np.sqrt(weights@(a*a)/weights.sum()));variation=float(np.sqrt(weights@((a-mean)**2)/weights.sum()))
        out[name]=dict(signed_mean=mean,weighted_RMS=rms,spatial_variation_RMS=variation,scaled_RMS=rms/scale,scale=scale)
    out['outer_minus_inner_mean']=out['outer']['signed_mean']-out['inner']['signed_mean']
    return out

def compact_traction(p):
    return {k:p[k] for k in ('state_sha256','mean_wall_shear','normalized_wall_traction_error','normalized_interface_traction_RMS','wall_spatial_standard_deviation','minimum_wall','maximum_wall')}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('preflight','pair','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();pre=read(a.preflight/'results.json');pair=read(a.pair/'results.json');c=pre['case'];n=pre['global_count'];w=np.asarray(pre['weights']);scale=c['ustar']**2
    require(pair['target']==.008 and set(pair['members'])=={'0.002','0.001'},'wrong paired scope')
    m=RequestManifest(a.output,dict(schema='CT8ReadOnlyPressureV1',preflight=sha(a.preflight/'results.json'),pair=sha(a.pair/'results.json'),source=sha(Path(__file__)),scope='unchanged absolute admission; offsets and variation diagnostic only; no solves'))
    result=dict(schema='CT8ReadOnlyPressureV1',members={},equal_time_spread=[],status='pass',qualification='original temporal:Pi failure remains; diagnostic completion is not preservation admission')
    seen={};times={}
    for name,row in pair['members'].items():
        initial=verified(Path(pre['initial'][name]['directory'])/'initial.json',pre['initial'][name]['packet_sha256']);o=initial['observation'];off=o['offsets'];regions=np.asarray(o['grid_regions'][:n]);initial_state=initial['observation']['state']
        native=verified(Path(row['directory'])/'case.evolution.json',row['native_report_sha256'])
        require(native['state']==native['accepted_states'][-1] and native['time']==.008 and not native['failure'],'unequal or rejected terminal state')
        def at(s):
            s=np.asarray(s)
            return np.asarray([row_value(q['pressure'],s,off[4+int(r)])[0] for q,r in zip(initial['test_basis'][:n],regions,strict=True)])
        # Parallel equilibrium has piecewise constant augmented pressure. This
        # avoids recomputing unrelated continuous k/velocity quadratures.
        start=at(initial_state);truth=np.where(regions==1,reference(col(c),c['H'],True)[3],0.)
        states=[];seen[name]={};times[name]={}
        for t,s in zip(native['audit_times'],native['accepted_states'],strict=True):
            now=at(s);seen[name][round(t,12)]=now;times[name][round(t,12)]=exact_vector_hash(s)
            states.append(dict(time=t,state_sha256=exact_vector_hash(s),drift=decomposition(now-start,w,regions,scale),reference_error=decomposition(now-truth,w,regions,scale)))
        stages=[];attempts=[q for q in native['pimple_attempts'] if q['scheduled_complete']]
        for attempt in [attempts[0],attempts[-1]]:
            flow=verified(attempt['flow_snapshot'],attempt['flow_snapshot_sha256']);ks=verified(attempt['k_snapshot'],attempt['k_snapshot_sha256'])
            snapshots={k:flow[k] for k in ('input','momentum_predictor','corrected_flow')};snapshots.update(post_epsilon=ks['pre_k'],post_k=ks['corrected'])
            require(flow['corrected_flow']==ks['post_flow'],'foreign flow/k snapshot')
            phase={}
            for label,s in snapshots.items():
                q=dict(state_sha256=exact_vector_hash(s),pressure_drift=decomposition(at(s)-start,w,regions,scale))
                # Predictor stress admissibility is diagnostic, not a new method gate.
                try:q['traction']=compact_traction(profiles(initial,s,c))
                except (ValueError,RuntimeError,AssertionError) as error:q['traction']=dict(status='diagnostic_unavailable',reason=str(error))
                phase[label]=q
            phase['post_epsilon']['pressure_coefficients_equal_corrected_flow']=ks['pre_k'][off[4]:off[6]]==flow['corrected_flow'][off[4]:off[6]]
            phase['post_k']['pressure_coefficients_equal_corrected_flow']=ks['corrected'][off[4]:off[6]]==flow['corrected_flow'][off[4]:off[6]]
            identity={key:attempt[key] for key in ('predictor_linear_residual','pressure_linear_residual','momentum_splitting_identity_defect','native_full_G_action_defect','frozen_momentum_defect','boundary_defect') if key in attempt}
            for key in ('refreshed_trace_defect','epsilon_frozen_linear_residual','k_frozen_linear_residual'):
                if key in attempt:identity[key]=attempt[key]
            identity['algebraic_continuity_maximum']=max(map(abs,attempt.get('algebraic_corrected_continuity',[])),default=0.)
            stages.append(dict(time=flow['time'],flow_snapshot=attempt['flow_snapshot'],flow_snapshot_sha256=attempt['flow_snapshot_sha256'],k_snapshot=attempt['k_snapshot'],k_snapshot_sha256=attempt['k_snapshot_sha256'],phases=phase,discrete_identities=identity))
        result['members'][name]=dict(native_report_sha256=row['native_report_sha256'],initial_packet_sha256=pre['initial'][name]['packet_sha256'],actual_initial_sha256=exact_vector_hash(initial_state),history_sha256=row['history_sha256'],initialized_projection=decomposition(start-truth,w,regions,scale),accepted_states=states,first_and_last_stages=stages)
    for t in sorted(set(seen['0.002'])&set(seen['0.001'])):
        result['equal_time_spread'].append(dict(time=t,states={k:times[k][t] for k in times},pressure=decomposition(seen['0.002'][t]-seen['0.001'][t],w,regions,scale)))
    result['interpretation']=dict(pressure_time='stored absolute Pi after flow update at native step time; scalar stages retain pressure coefficients and change current constitutive/trace data',uncertainty='offset, variation and interface-mean diagnostics do not change absolute pressure metric or remove temporal:Pi',limits='a single two-member comparison cannot prove an asymptotic order or isolate spatial/splitting effects; no further timestep or mesh search in CT8')
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);print('CT8 read-only pressure diagnosis complete; original qualification unchanged',flush=True);return 0
if __name__=='__main__':raise SystemExit(main())
