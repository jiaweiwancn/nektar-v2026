#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only final return/window checks for one completed recovery leg."""
import argparse
from pathlib import Path
import numpy as np
from analyze_wp2_approach import assess,native_states,tensor_differences,K0
from wp2_all_state import Rows,US
from wp2_common_quadrature import resample
from wp2_independent_reference import evaluate
from audit_wp2_dense_profiles import gates
from wp2_qualification import quarter_budget_differences
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import require,sha,state_hash


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--model',choices=['CLS','standard'],required=True);p.add_argument('--trajectory',type=Path)
    a=p.parse_args();a.output.mkdir(exist_ok=False);progress_path=a.trajectory or a.root/f'preservation-recovery-{a.model}/progress.json'
    level=load(progress_path)['case_definition']['level'];require(level in (4,5),'unsupported recovery row catalogue')
    rows=Rows(a.root/f'geometry-rows-{"confirmation" if level==5 else "selected"}/result.json')
    assessment,end,_=assess(a.root,a.model,'recovery',rows,a.output/'all-states.json',progress_path)
    states,context=native_states(a.root,a.model,'recovery',progress_path);context.update(ustar=end['ustar'],developing_boundary=end['developing_boundary'])
    dense=[];reference=None;final_start=end['time']-.3
    times=[];acc0=None;acc1=None;tensor0={};tensor1={};first=None
    for _,(t,s) in sorted(states.items()):
        v=rows.tuple(s,context);sample=dict(context['wp2_common_samples'],values=v['values'])
        actual=resample(dict(context,wp2_common_samples=sample,physical_samples=dict(wall_K=v['endpoint_K'])),28)
        if reference is None:reference=evaluate(context,actual['physical'])
        dense.append(dict(actual_time=t,state_sha256=state_hash(s),**gates(actual,reference)))
        if t>=final_start-2e-14:
            if first is None:
                require(abs(t-final_start)<2e-14,'dense final window has no exact retained starting time')
                first=actual
                acc0=np.zeros_like(actual['values']);acc1=np.zeros_like(acc0)
                tensor0={n:np.zeros_like(actual[n]) for n in ['R','T']};tensor1={n:np.zeros_like(actual[n]) for n in ['R','T']}
            tau=t-final_start;times.append(tau);acc0+=actual['values'];acc1+=tau*actual['values']
            for name in tensor0:tensor0[name]+=actual[name];tensor1[name]+=tau*actual[name]
    atomic(a.output/'dense-profiles.json',dict(history=dense,final_passed=dense[-1]['passed']))
    progress=load(progress_path);start=end['time']-.3
    old=next(x for x in progress['observations'] if abs(x['time']-start)<2e-14);left=load(Path(old['report']))
    window=quarter_budget_differences(left['wp2_common_samples'],end['wp2_common_samples'],US,K0[a.model],.0065)
    tensors=tensor_differences(left,end,rows,rows)
    tt=np.asarray(times);den=float(np.sum((tt-np.mean(tt))**2));require(den>0,'dense recovery trend is unresolved')
    slope=(acc1-np.mean(tt)*acc0)/den
    base=dict(columns=context['wp2_common_samples']['columns'],weights=actual['weights'].tolist(),regions=actual['regions'].tolist(),reference_values=reference['values'].tolist())
    dense_window=quarter_budget_differences(dict(base,values=first['values']),dict(base,values=actual['values']),US,K0[a.model],.0065)
    from analyze_wp2_approach import TREQ
    dense_trend=quarter_budget_differences(dict(base,values=reference['values']+TREQ*slope),dict(base,values=reference['values']),US,K0[a.model],.0065)
    tensor_window={};tensor_trend={}
    for side in [0,1]:
        mask=actual['regions']==side;w=actual['weights'][mask]
        for name in tensor0:
            d=actual[name][mask]-first[name][mask]
            ds=TREQ*(tensor1[name][mask]-np.mean(tt)*tensor0[name][mask])/den
            tensor_window[str(side)+'_'+name]=float(np.sqrt(np.sum(w[:,None,None]*d*d)/sum(w))/US**2)
            tensor_trend[str(side)+'_'+name]=float(np.sqrt(np.sum(w[:,None,None]*ds*ds)/sum(w))/US**2)
    dense_window_pass=dense_window['passed'] and dense_trend['passed'] and max(tensor_window.values())<=.0025 and max(tensor_trend.values())<=.0025
    result=dict(version='WP2RecoveryReturnV3',model=a.model,actual_time=end['time'],state_sha256=state_hash(end['state']),
        final_physical_passed=assessment['final_physical_and_full_tensor_passed'],final_GL28_profiles_passed=dense[-1]['passed'],
        final_window=window,final_window_full_tensors=tensors,final_window_trend=assessment['final_window']['transit_scaled_linear_trend'],
        final_report=progress['observations'][-1]['report'],final_report_sha256=sha(Path(progress['observations'][-1]['report'])),
        native_progress=str(progress_path),native_progress_sha256=sha(progress_path),
        GL28_final_window=dense_window,GL28_transit_scaled_trend=dense_trend,GL28_full_tensor_window=tensor_window,GL28_full_tensor_trend=tensor_trend,
        GL28_window_passed=dense_window_pass,
        source_sha256={n:sha(Path(__file__).with_name(n)) for n in [Path(__file__).name,'analyze_wp2_approach.py','wp2_all_state.py','wp2_common_quadrature.py','wp2_independent_reference.py','audit_wp2_dense_profiles.py']},
        accuracy_failure_reporting='failed recovery accuracy checks remain failed and are retained; scalar/covariance requirements are never bypassed',
        passed=assessment['final_physical_and_full_tensor_passed'] and dense[-1]['passed'] and window['passed'] and tensors['passed'] and assessment['final_window']['transit_scaled_linear_trend']['passed'] and dense_window_pass,
        scope='same complete recovery trajectory, no field or reference alteration; finite-time return and measured final-window sensitivity; zero PDE solves in this assessment')
    atomic(a.output/'result.json',result)


if __name__=='__main__':main()
