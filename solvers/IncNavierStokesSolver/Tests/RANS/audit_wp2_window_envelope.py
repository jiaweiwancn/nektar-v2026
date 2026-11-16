#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Conservative all-pairs final-window change bounds on the fixed GL28 grid."""
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
from cao_tamura_h2d_state import load


def range_checks(minimum,maximum,reference,columns,weights,regions,K):
    base=dict(columns=columns,weights=weights.tolist(),regions=regions.tolist(),reference_values=reference['values'].tolist())
    principal=quarter_budget_differences(dict(base,values=maximum['values']),dict(base,values=minimum['values']),US,K,.0065)
    tensors={}
    for side in [0,1]:
        mask=regions==side;w=weights[mask]
        for name in ['R','T']:
            d=maximum[name][mask]-minimum[name][mask]
            tensors[str(side)+'_'+name]=float(np.sqrt(np.sum(w[:,None,None]*d*d)/sum(w))/US**2)
    return dict(principal=principal,full_tensor_RMS_range_bound=tensors,
                passed=principal['passed'] and max(tensors.values())<=.0025,
                interpretation='componentwise range gives an upper bound for every pair of sampled times; the envelope is not a constructed physical state')


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--models',nargs='+',choices=['CLS','standard'],default=['CLS','standard'])
    p.add_argument('--phases',nargs='+',choices=['coarse','fine','confirmation','recovery'],default=['coarse','fine','confirmation'])
    p.add_argument('--trajectory',type=Path);a=p.parse_args();a.output.mkdir(exist_ok=False);started=time.monotonic();results={}
    require(not a.trajectory or (len(a.models)==1 and a.phases==['recovery']),'an explicit recovery trajectory requires exactly one model/phase')
    for model in a.models:
        results[model]={}
        for phase in a.phases:
            progress_path=a.trajectory or a.root/f'preservation-{phase}-{model}/progress.json';progress_hash=sha(progress_path);progress=load(progress_path)
            require(progress['status']=='completed_bounded_trajectory','window trajectory incomplete')
            states,context=native_states(a.root,model,phase,progress_path);end=load(Path(progress['observations'][-1]['report']))
            context.update(ustar=end['ustar'],developing_boundary=end['developing_boundary'])
            rows=Rows(a.root/f'geometry-rows-{"confirmation" if phase=="confirmation" else "selected"}/result.json')
            times=[];hashes=[];minimum=None;maximum=None;reference=None
            for _,(t,s) in sorted(states.items()):
                if t<end['time']-.3-2e-14:continue
                v=rows.tuple(s,context);actual=resample(dict(context,wp2_common_samples=dict(context['wp2_common_samples'],values=v['values']),
                                                          physical_samples=dict(wall_K=v['endpoint_K'])),28)
                if minimum is None:
                    minimum={n:actual[n].copy() for n in ['values','R','T']};maximum={n:actual[n].copy() for n in minimum}
                    reference=evaluate(context,actual['physical'])
                else:
                    for n in minimum:np.minimum(minimum[n],actual[n],out=minimum[n]);np.maximum(maximum[n],actual[n],out=maximum[n])
                times.append(t);hashes.append(state_hash(s))
            require(times[-1]-times[0]>=TREQ/4-2e-14,'all-pairs window is too short')
            checks=range_checks(minimum,maximum,reference,context['wp2_common_samples']['columns'],actual['weights'],actual['regions'],K0[model])
            require(sha(progress_path)==progress_hash,'trajectory changed while assessing its window')
            result=dict(version='WP2AllPairWindowEnvelopeV1',model=model,phase=phase,actual_times=times,state_sha256=hashes,
                        duration=times[-1]-times[0],required=TREQ/4,trajectory=str(progress_path),trajectory_sha256=progress_hash,
                        points_per_state=len(actual['weights']),checks=checks,source_sha256=sha(Path(__file__)),governing_solves=0)
            target=a.output/f'{model}-{phase}.json';atomic(target,result);results[model][phase]=dict(path=str(target),sha256=sha(target),passed=checks['passed'])
    atomic(a.output/'summary.json',dict(results=results,seconds=time.monotonic()-started,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))


if __name__=='__main__':main()
