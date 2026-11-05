#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compare saved native single-sweep trajectories in physical quadrature norms."""
import argparse
import math
from pathlib import Path
from run_cls_epawf_picard import load, require, write_json, state_hash


def norm(weights, values):
    return math.sqrt(max(0., math.fsum(w*v*v for w,v in zip(weights, values))))


def differences(a, b):
    require(abs(a['time']-b['time']) < 1e-13 and a['initial'] == b['initial'],
            'comparison changed time or initial state')
    aa, bb = a['physical_samples'], b['physical_samples']
    require(aa and bb, 'comparison needs a completed numerical trajectory')
    result = {}
    pressure_diff, pressure_weights = [], []
    for i in range(6):
        w = aa['weights_'+str(i)]
        require(w == bb['weights_'+str(i)], 'changed quadrature')
        delta = [x-y for x,y in zip(aa['values_'+str(i)], bb['values_'+str(i)])]
        if i < 4:
            result[('u','v','k','epsilon')[i]] = norm(w, delta)
        else:
            pressure_weights.extend(w);pressure_diff.extend(delta)
    result['velocity'] = math.hypot(result.pop('u'), result.pop('v'))
    mean = math.fsum(w*d for w,d in zip(pressure_weights, pressure_diff))/math.fsum(pressure_weights)
    result['pressure_common_gauge_shift'] = mean
    result['pressure'] = norm(pressure_weights, [d-mean for d in pressure_diff])
    w = aa['wall_weights']
    require(w == bb['wall_weights'], 'changed wall quadrature')
    for name in ('wall_stress','wall_B','wall_lambda','wall_K'):
        result[name] = norm(w, [x-y for x,y in zip(aa[name],bb[name])])
    return result


def describe(report):
    audit = report['state_audits'][-1]
    row = report['pimple_attempts'][-1]
    domain = audit['r2_scalar_domain']
    return dict(dt=report['dt'], accepted_steps=report['steps'], final_time=report['time'],
        continuity_gauge=row['corrected_continuity_gauge'],
        strong_continuity=audit['strong_continuity'],
        minimum_inner_k=domain[2], minimum_outer_k=domain[3], minimum_outer_epsilon=domain[4],
        minimum_sampled_K=domain[5], minimum_projected_epsilon=domain[6],
        minimum_B=domain[7], minimum_lambda=domain[8],
        imposed_trace_defect=row['imposed_trace_defect'],refreshed_trace_defect=row['refreshed_trace_defect'],
        original_residual=max(row['refreshed_original_block_norms']),
        qualified=audit['r2_qualified'][0] == 1,
        first_qualification_exit=report.get('physical_range_events', [{}])[0].get('time') if report.get('physical_range_events') else None,
        seconds=sum(r['seconds'] for r in report['pimple_attempts']))


def same_history(a, b):
    for key in ('state','native_history','native_history_times','time','steps'):
        # Hash the raw binary64, including signed zeros, for vector payloads.
        x,y = a[key],b[key]
        if key == 'native_history':
            require([state_hash(v) for v in x] == [state_hash(v) for v in y], key+' changed')
        elif isinstance(x,list):
            require(state_hash(x) == state_hash(y), key+' changed')
        else:
            require(x == y, key+' changed')
    return True


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('reports',type=Path,nargs='+')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--identical-history',action='store_true')
    a=p.parse_args();reports=[load(f) for f in a.reports]
    result=dict(reports=[str(f.resolve()) for f in a.reports])
    if a.identical_history:
        require(len(reports)==2,'history comparison requires two reports')
        result['bitwise_state_time_and_all_native_history_equal']=same_history(*reports)
    else:
        result['runs']=[describe(r) for r in reports]
        diffs=[differences(x,y) for x,y in zip(reports,reports[1:])]
        result['successive_physical_L2_differences']=diffs
        if len(diffs)==2:
            require(abs(reports[0]['dt']/reports[1]['dt']-2)<1e-12 and
                    abs(reports[1]['dt']/reports[2]['dt']-2)<1e-12,'order estimate requires dt halving')
            result['observed_successive_difference_orders']={k:math.log2(diffs[0][k]/diffs[1][k])
                for k in diffs[0] if k!='pressure_common_gauge_shift' and diffs[0][k]>0 and diffs[1][k]>0}
    write_json(a.output,result)
    print(a.output)


if __name__=='__main__':
    main()
