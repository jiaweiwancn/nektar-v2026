#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Quarter-budget comparisons on every executed common time at volume GL28."""
import argparse,time,resource
from pathlib import Path
import numpy as np
from analyze_wp2_approach import native_states,K0,TREQ
from wp2_all_state import Rows,US
from wp2_common_quadrature import resample
from wp2_independent_reference import evaluate
from wp2_qualification import quarter_budget_differences
from run_cls_epawf_picard import require,sha,state_hash
from continue_cls_epawf_pimple import atomic


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(exist_ok=False);started=time.monotonic();summary={}
    rows={name:Rows(a.root/f'geometry-rows-{name}/result.json') for name in ['selected','confirmation']}
    for model in ['CLS','standard']:
        states={};contexts={}
        for phase in ['coarse','fine','confirmation']:
            states[phase],contexts[phase]=native_states(a.root,model,phase)
            from cao_tamura_h2d_state import load
            end=load(Path(load(a.root/f'preservation-{phase}-{model}/progress.json')['observations'][-1]['report']))
            contexts[phase].update(ustar=end['ustar'],developing_boundary=end['developing_boundary'])
        times=sorted(set(states['coarse'])&set(states['fine'])&set(states['confirmation']));history=[];base=None
        require(states['coarse'][times[-1]][0]>=TREQ,'dense comparison misses physical horizon')
        for key in times:
            actual={};hashes={};physical_times={}
            for phase in states:
                t,s=states[phase][key];context=contexts[phase];catalog=rows['confirmation' if phase=='confirmation' else 'selected']
                v=catalog.tuple(s,context)
                actual[phase]=resample(dict(context,wp2_common_samples=dict(context['wp2_common_samples'],values=v['values']),
                                            physical_samples=dict(wall_K=v['endpoint_K'])),28)
                hashes[phase]=state_hash(s);physical_times[phase]=t
            require(max(physical_times.values())-min(physical_times.values())<2e-14,'dense comparison relabels physical times')
            if base is None:
                reference=evaluate(contexts['coarse'],actual['coarse']['physical'])
                base=dict(columns=contexts['coarse']['wp2_common_samples']['columns'],weights=actual['coarse']['weights'].tolist(),
                          regions=actual['coarse']['regions'].tolist(),reference_values=reference['values'].tolist())
            comparisons={}
            for name,right in [('timestep','fine'),('mesh','confirmation')]:
                left=actual['coarse'];other=actual[right]
                require(np.array_equal(left['physical'],other['physical']) and np.array_equal(left['weights'],other['weights']),'dense physical grids differ')
                check=quarter_budget_differences(dict(base,values=left['values'].tolist()),dict(base,values=other['values'].tolist()),US,K0[model],.0065)
                tensors={}
                for side in [0,1]:
                    mask=left['regions']==side;w=left['weights'][mask]
                    for tensor in ['R','T']:
                        tensors[str(side)+'_'+tensor]=float(np.sqrt(np.sum(w[:,None,None]*(left[tensor][mask]-other[tensor][mask])**2)/sum(w))/US**2)
                comparisons[name]=dict(principal=check,full_tensor_RMS=tensors,passed=check['passed'] and max(tensors.values())<=.0025)
                require(comparisons[name]['passed'],'dense quarter-budget comparison failed')
            history.append(dict(actual_times=physical_times,state_sha256=hashes,comparisons=comparisons))
        result=dict(version='WP2GL28CommonTimeDifferencesV1',model=model,history=history,common_native_times=len(times),
                    points_per_state=len(actual['coarse']['weights']),passed=True,
                    scope='fixed seen GL28 union panels, actual native primitive polynomials, fresh selected stress evaluation, absolute pressure, no time interpolation')
        atomic(a.output/(model+'.json'),result);summary[model]=dict(path=str(a.output/(model+'.json')),sha256=sha(a.output/(model+'.json')),passed=True)
    atomic(a.output/'summary.json',dict(models=summary,seconds=time.monotonic()-started,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                                      source_sha256=sha(Path(__file__)),governing_solves=0))


if __name__=='__main__':main()
