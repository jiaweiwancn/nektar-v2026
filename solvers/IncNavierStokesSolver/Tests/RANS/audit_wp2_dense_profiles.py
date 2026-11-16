#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""All-state GL28 profile/gradient/full-stress accuracy, without PDE solves.

Complements the native exact sided/essential/flux inventories. Neither the
selected state nor the fixed analytical reference is projected or modified.
"""
import argparse,time,resource
from pathlib import Path
import numpy as np
from analyze_wp2_approach import native_states,K0
from wp2_all_state import Rows,US
from wp2_common_quadrature import resample
from wp2_independent_reference import evaluate
from run_cls_epawf_picard import require,sha,state_hash
from continue_cls_epawf_pimple import atomic
from cao_tamura_h2d_state import load


def gates(actual,reference):
    v=actual['values'];r=reference['values'];w=actual['weights'];region=actual['regions'];out={}
    def add(name,value,limit):
        finite=bool(np.isfinite(value));out[name]=dict(value=float(value) if finite else None,limit=limit,passed=finite and bool(value<=limit))
    for side in [0,1]:
        m=region==side;weight=w[m];volume=sum(weight);d=v[m]-r[m];prefix=str(side)+'_'
        add(prefix+'velocity_L2',np.sqrt(np.sum(weight[:,None]*d[:,:2]**2)/np.dot(weight,r[m,0]**2)),.001)
        add(prefix+'W_RMS',np.sqrt(np.dot(weight,d[:,1]**2)/volume)/US,.001)
        for j,name in [(0,'U'),(1,'W'),(2,'k'),(3,'epsilon'),(6,'nut')]:
            if j==3 and side==0:continue
            scale=np.maximum(abs(r[m,j]),US) if j==0 else US if j==1 else abs(r[m,j])
            add(prefix+name+'_point',np.max(abs(d[:,j])/scale),.01)
            if j in [2,3,6]:add(prefix+name+'_local_relative_RMS',np.sqrt(np.dot(weight,(d[:,j]/scale)**2)/volume),.01)
        dg=actual['G'][m]-reference['G'][m]
        add(prefix+'full_gradient_L2',np.sqrt(np.sum(weight[:,None,None]*dg**2)/np.sum(weight[:,None,None]*reference['G'][m]**2)),.01)
        for j,name in [(4,'pi'),(5,'mechanical_pressure'),(7,'shear'),(8,'normal_stress')]:add(prefix+name+'_RMS',np.sqrt(np.dot(weight,d[:,j]**2)/volume)/US**2,.01)
        for name in ['R','T']:add(prefix+name+'_full_tensor_RMS',np.sqrt(np.sum(weight[:,None,None]*(actual[name][m]-reference[name][m])**2)/volume)/US**2,.01)
    z=actual['physical'][:,2];edges=[0.,.0065,.012,.025,.05,.1,.2,.4,.8]
    for lo,hi in zip(edges,edges[1:]):
        m=(z>lo)&(z<hi);weight=w[m];volume=sum(weight)
        for j,name in [(0,'U'),(1,'W'),(2,'k'),(3,'epsilon'),(6,'nut')]:
            if j==3 and hi<=.0065:continue
            scale=np.maximum(abs(r[m,j]),US) if j==0 else US if j==1 else abs(r[m,j])
            add('band_'+str(lo)+'_'+str(hi)+'_'+name,np.sqrt(np.dot(weight,((v[m,j]-r[m,j])/scale)**2)/volume),.001 if j in [0,1] else .01)
    return dict(passed=all(bool(q['passed']) for q in out.values()),gates=out)


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--phases',nargs='+',default=['coarse','fine','confirmation','recovery'])
    a=p.parse_args();require(not a.output.exists(),'dense profile audit output already exists');a.output.mkdir()
    rows={n:Rows(a.root/f'geometry-rows-{n}/result.json') for n in ['selected','confirmation']};result={}
    for model in ['CLS','standard']:
        result[model]={}
        for phase in a.phases:
            start=time.monotonic();catalog=rows['confirmation' if phase=='confirmation' else 'selected']
            states,context=native_states(a.root,model,phase);endpoint=load(Path(load(a.root/f'preservation-{phase}-{model}/progress.json')['observations'][-1]['report']))
            # Boundary datum is immutable, and is used only by the separate
            # analytical reference. No actual state values enter its definition.
            context['developing_boundary']=endpoint['developing_boundary'];context['ustar']=endpoint['ustar']
            reference=None;history=[]
            for _,(t,s) in sorted(states.items()):
                tuple_data=catalog.tuple(s,context);sample=dict(context['wp2_common_samples'],values=tuple_data['values'])
                actual=resample(dict(context,wp2_common_samples=sample,physical_samples=dict(wall_K=tuple_data['endpoint_K'])),28)
                if reference is None:reference=evaluate(context,actual['physical'])
                checks=gates(actual,reference)
                history.append(dict(actual_time=t,state_sha256=state_hash(s),**checks))
                if phase!='recovery' and not checks['passed']:
                    atomic(a.output/f'failed-{model}-{phase}.json',dict(history=history,model=model,phase=phase,scope='read-only failed dense-profile audit; no state change'))
                    raise RuntimeError('GL28 full profile accuracy failed')
            d=dict(version='WP2GL28PhysicalProfilesV1',model=model,phase=phase,history=history,
                final_passed=history[-1]['passed'],all_states_passed=all(q['passed'] for q in history),
                point_count=len(actual['weights']),seconds=time.monotonic()-start,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                source_identity={name:sha(Path(__file__).with_name(name)) for name in [Path(__file__).name,'wp2_independent_reference.py','wp2_common_quadrature.py','wp2_all_state.py']},
                context_sha256=sha(catalog.path),scope='same frozen union-panel GL28; actual native polynomials and selected tensors; independent analytical log/AWF primitive including nonconstant ReferenceK; exact native sided checks remain complementary; finite inventory')
            atomic(a.output/f'{model}-{phase}.json',d);result[model][phase]=dict(path=str(a.output/f'{model}-{phase}.json'),sha256=sha(a.output/f'{model}-{phase}.json'),
                final_passed=d['final_passed'],all_states_passed=d['all_states_passed'],seconds=d['seconds'])
    atomic(a.output/'summary.json',result)


if __name__=='__main__':main()
