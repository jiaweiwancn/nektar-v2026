#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only exact-time pressure-control preservation and sensitivity checks."""
import argparse,time,resource
from pathlib import Path
import numpy as np
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import require,sha,state_hash
from analyze_wp2_approach import assess,native_states,K0
from wp2_all_state import Rows,US
from wp2_common_quadrature import resample
from wp2_independent_reference import evaluate
from wp2_qualification import quarter_budget_differences
from audit_wp2_dense_profiles import gates
from audit_wp2_window_envelope import range_checks


def dense(rows,s,context):
    v=rows.tuple(s,context)
    return resample(dict(context,wp2_common_samples=dict(context['wp2_common_samples'],values=v['values']),
                        physical_samples=dict(wall_K=v['endpoint_K'])),36)


def main():
    p=argparse.ArgumentParser()
    for n in ['root','output','selected-rows','confirmation-rows']:p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--model',choices=['CLS','standard'],required=True);a=p.parse_args();a.output.mkdir(exist_ok=False)
    started=time.monotonic();r=a.root;catalogues=[Rows(a.selected_rows),Rows(a.confirmation_rows)];summary={};all_states={};contexts={}
    for phase in ['coarse','fine','confirmation','recovery']:
        path=r/('pressure-'+phase+'-'+a.model)/'progress.json';progress=load(path)
        require(progress['case_definition']['pressure_space']=='inner-p2-v1','pressure case identity')
        rows=catalogues[int(phase=='confirmation')]
        result,end,_=assess(r,a.model,phase,rows,a.output/(phase+'-all-states.json'),path,independent_order=36)
        states,context=native_states(r,a.model,phase,path);context.update(ustar=end['ustar'],developing_boundary=end['developing_boundary'])
        reference=None;profiles=[];minimum=None;maximum=None
        for key,(t,s) in sorted(states.items()):
            actual=dense(rows,s,context)
            if reference is None:reference=evaluate(context,actual['physical'])
            profiles.append(dict(time=t,state_sha256=state_hash(s),**gates(actual,reference)))
            if t>=end['time']-.3-2e-14:
                if minimum is None:minimum={n:actual[n].copy() for n in ['values','R','T']};maximum={n:actual[n].copy() for n in minimum}
                else:
                    for n in minimum:np.minimum(minimum[n],actual[n],out=minimum[n]);np.maximum(maximum[n],actual[n],out=maximum[n])
        window=range_checks(minimum,maximum,reference,context['wp2_common_samples']['columns'],actual['weights'],actual['regions'],K0[a.model])
        atomic(a.output/(phase+'-GL36-profiles.json'),dict(profiles=profiles,final_window_envelope=window))
        passed=result['final_physical_and_full_tensor_passed'] and profiles[-1]['passed'] and window['passed'] and result['final_window']['transit_scaled_linear_trend']['passed']
        if phase!='recovery':passed=passed and all(v['passed'] for v in profiles) and all(v['physical_checks']['passed'] for v in result['history'])
        summary[phase]=dict(passed=passed,actual_time=end['time'],final_profile_passed=profiles[-1]['passed'],window_passed=window['passed'],
            trend_passed=result['final_window']['transit_scaled_linear_trend']['passed'],progress_sha256=sha(path),
            final_native_checks=result['history'][-1]['physical_checks'],final_state_sha256=state_hash(end['state']))
        all_states[phase]=states;contexts[phase]=context
    comparisons={}
    for left,right in [('coarse','fine'),('coarse','confirmation')]:
        pairs=[]
        for key in sorted(all_states[left].keys() & all_states[right].keys()):
            tl,sl=all_states[left][key];tr,sr=all_states[right][key]
            require(abs(tl-tr)<3e-14,'actual common times differ beyond addition tolerance')
            l=dense(catalogues[0],sl,contexts[left]);rr=dense(catalogues[int(right=='confirmation')],sr,contexts[right])
            ref=evaluate(contexts[left],l['physical']);base=dict(columns=contexts[left]['wp2_common_samples']['columns'],weights=l['weights'].tolist(),regions=l['regions'].tolist(),reference_values=ref['values'].tolist())
            check=quarter_budget_differences(dict(base,values=l['values']),dict(base,values=rr['values']),US,K0[a.model],.0065)
            tensors={}
            for region in [0,1]:
                mask=l['regions']==region;w=l['weights'][mask]
                for name in ['R','T']:
                    tensors[str(region)+'_'+name]=float(np.sqrt(np.sum(w[:,None,None]*(l[name][mask]-rr[name][mask])**2)/sum(w))/US**2)
            pairs.append(dict(actual_left_time=tl,actual_right_time=tr,principal=check,full_tensor_RMS_over_ustar2=tensors,
                passed=check['passed'] and max(tensors.values())<=.0025))
        require(len(pairs)>=61,'missing exact-time sensitivity coverage')
        comparison=dict(passed=all(v['passed'] for v in pairs),pairs=pairs)
        atomic(a.output/(left+'-'+right+'-GL36.json'),comparison);comparisons[left+'_'+right]=dict(passed=comparison['passed'],common_times=len(pairs))
    atomic(a.output/'result.json',dict(version='WP2InnerP2PhysicalAssessmentV1',model=a.model,phases=summary,comparisons=comparisons,
        passed=all(x['passed'] for x in summary.values()) and all(x['passed'] for x in comparisons.values()),
        seconds=time.monotonic()-started,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        source_sha256=sha(Path(__file__)),catalogue_sha256=[sha(x.path) for x in catalogues],governing_solves=0,
        scope='GL36 declared before pressure-control tuning; actual retained polynomials, native hard audits and unchanged absolute regional pressure; no interpolation in time'))


if __name__=='__main__':main()
