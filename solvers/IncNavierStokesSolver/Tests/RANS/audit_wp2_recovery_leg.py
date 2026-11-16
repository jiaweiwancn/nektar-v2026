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
    p.add_argument('--model',choices=['CLS','standard'],required=True);a=p.parse_args();a.output.mkdir(exist_ok=False)
    rows=Rows(a.root/'geometry-rows-selected/result.json')
    assessment,end,_=assess(a.root,a.model,'recovery',rows,a.output/'all-states.json')
    states,context=native_states(a.root,a.model,'recovery');context.update(ustar=end['ustar'],developing_boundary=end['developing_boundary'])
    dense=[];reference=None
    for _,(t,s) in sorted(states.items()):
        v=rows.tuple(s,context);sample=dict(context['wp2_common_samples'],values=v['values'])
        actual=resample(dict(context,wp2_common_samples=sample,physical_samples=dict(wall_K=v['endpoint_K'])),28)
        if reference is None:reference=evaluate(context,actual['physical'])
        dense.append(dict(actual_time=t,state_sha256=state_hash(s),**gates(actual,reference)))
    atomic(a.output/'dense-profiles.json',dict(history=dense,final_passed=dense[-1]['passed']))
    progress=load(a.root/f'preservation-recovery-{a.model}/progress.json');start=end['time']-.3
    old=next(x for x in progress['observations'] if abs(x['time']-start)<2e-14);left=load(Path(old['report']))
    window=quarter_budget_differences(left['wp2_common_samples'],end['wp2_common_samples'],US,K0[a.model],.0065)
    tensors=tensor_differences(left,end,rows,rows)
    result=dict(version='WP2RecoveryReturnV1',model=a.model,actual_time=end['time'],state_sha256=state_hash(end['state']),
        final_physical_passed=assessment['final_physical_and_full_tensor_passed'],final_GL28_profiles_passed=dense[-1]['passed'],
        final_window=window,final_window_full_tensors=tensors,final_window_trend=assessment['final_window']['transit_scaled_linear_trend'],
        final_report=progress['observations'][-1]['report'],final_report_sha256=sha(Path(progress['observations'][-1]['report'])),
        native_progress_sha256=sha(a.root/f'preservation-recovery-{a.model}/progress.json'),
        passed=assessment['final_physical_and_full_tensor_passed'] and dense[-1]['passed'] and window['passed'] and tensors['passed'] and assessment['final_window']['transit_scaled_linear_trend']['passed'],
        scope='same complete recovery trajectory, no field or reference alteration; finite-time return and measured final-window sensitivity; zero PDE solves in this assessment')
    atomic(a.output/'result.json',result)


if __name__=='__main__':main()
