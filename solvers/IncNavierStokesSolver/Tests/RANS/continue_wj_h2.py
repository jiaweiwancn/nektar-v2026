#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""One revised, finite H2 design using the admitted refined half-dt history.

Original study-v1 and diagnosis-v1 remain immutable. No unchanged rejected
candidate is retried and no further timestep design is authorized here.
"""
import argparse
import math
import os
from pathlib import Path
import time
import numpy as np
from run_wj_h2 import execute,load,save,sha
from study_wj_h2 import accepted,observation_rule,comparison
from wj_h2_cases import case_xml
from wj_h2_profiles import Completion,QE,KAPPA


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--solver',type=Path,required=True)
    p.add_argument('--observer',type=Path,required=True);a=p.parse_args();root=a.root.resolve();out=root/'study-v2';out.mkdir(exist_ok=False)
    old=load(root/'study-v1/progress.json');diagnosis=load(root/'diagnosis-v1/result.json');measure=load(root/'observations-v3/summary.json')
    assert old['status']=='OBSTRUCTED' and diagnosis['qualification']=='PASS' and diagnosis['accepted_steps']==8
    for key in ('time08-q8','time08-q10','normal-fine08-q8','normal-fine08-q10'):
        assert measure['comparisons'][key]['status']=='PASS'
    index=load(root/'observations-v3/index.json')
    for key in ('base08-nativeweak','refined-fine08-nativeweak'):
        assert index[key]['independent_weak']['status']=='PASS'
    last={key:dict(old['runs'][name]) for key,name in [('base','base-pilot-on'),('fine','fine-pilot')]}
    last['refined']=dict(diagnosis,checkpoint=load(diagnosis['report'])['checkpoint'],dt=.01,refined=True,
        steps=diagnosis['accepted_steps'],status='PASS')
    for key,r in last.items():
        d=load(r['report']);assert sha(r['report'])==r['report_sha256'] and accepted(d) and d['time']==.08
        assert sha(a.solver)==load(root/'study-v1/intent.json')['solver']['sha256']
    immutable={r['report']:sha(r['report']) for r in last.values()}
    immutable[str(root/'study-v1/progress.json')]=sha(root/'study-v1/progress.json')
    save(out/'intent.json',dict(schema='WJH2OneRevisedDesignV1',driver_sha256=sha(__file__),
        solver=str(a.solver),solver_sha256=sha(a.solver),observer=str(a.observer),observer_sha256=sha(a.observer),
        imported=last,immutable=immutable,target=8.,time_pair={'base':.02,'fine':.01},
        pure_normal_pair={'fine':.01,'refined':.01},maximum_native_block=32,
        common_cadence=.5,no_further_halving=True,original_failed_refined_dt=.02,
        additional_missing_input_allowance=32,role='new bounded numerical design; original hard exit is never erased',
        weak_quadrature='same native nodes verified; increased-order residual difference remains separately UNRESOLVED'))
    state=dict(status='STARTED',target=8.,imported=last,runs={},observations={},comparisons={},first_limiting_gate=None)
    status=out/'progress.json';save(status,state)
    env=dict(LD_LIBRARY_PATH=os.environ['LD_LIBRARY_PATH']);c=Completion();eps=(c.K/QE)**1.5/(KAPPA*c.H)
    designs=[load(Path(r['xml']).with_name('case-design.json')) for r in last.values()]
    rules={o:observation_rule(designs,o,x) for o,x in ((8,4),(10,6))}
    rates={key:load(root/'operations'/name/'receipt.json')['seconds']/n for key,name,n in
        [('base','study-base-pilot-on',4),('fine','study-fine-pilot',8),('refined','diagnosis-refined-half-dt',8)]}
    for cadence in range(1,17):
        book=load(root/'accounting.json');spent=sum(v.get('seconds',v['reserved_seconds']) for v in book['entries'].values())+book['prior_profile_setup_seconds']
        charged=sum(v['units'] for v in book['entries'].values())
        remaining={key:round(8/(.02 if key=='base' else .01))-last[key]['steps'] for key in last}
        raw=sum(remaining[k]*rates[k] for k in last)+60*(17-cadence)+300
        projection=dict(cadence=cadence,measured_seconds_per_step_including_short_startup=rates,
            spent=spent,charged=charged,remaining_steps=remaining,remaining_raw_seconds=raw,
            headroom=1.25,closure=600,projected_seconds=spent+1.25*raw+600,
            projected_units=charged+sum(remaining.values())+32,seconds_ceiling=21600,units_ceiling=2200)
        projection['status']='PASS' if projection['projected_seconds']<=21600 and projection['projected_units']<=2200 else 'RESOURCE_LIMIT'
        save(out/f'admission-{cadence:02d}.json',projection)
        if projection['status']!='PASS':
            state.update(status='RESOURCE_LIMIT',first_limiting_gate='remaining numerical design admission');save(status,state);break
        target=cadence*.5
        # The tighter refined principal margin is evaluated first. No other
        # trajectory advances after a failed barrier/attempt.
        failed=False
        for key in ('refined','base','fine'):
            dt=.02 if key=='base' else .01;goal=round(target/dt)
            while last[key]['steps']<goal:
                start=last[key]['steps'];count=min(32,goal-start);label=f'{key}-{start:04d}-{start+count:04d}'
                xml=case_xml(out/label,'physical',count,dt,key=='refined','ExactStepV1',restart=last[key]['checkpoint'])
                okay=execute(root,'revised-'+label,[a.solver,xml],xml.parent,count,900,extra_env=env)
                rp=xml.with_suffix('.evolution.json');d=load(rp) if rp.exists() else None
                record=dict(xml=str(xml),xml_sha256=sha(xml),report=str(rp),report_sha256=sha(rp) if rp.exists() else None,
                    status='PASS' if okay and d is not None and accepted(d) else 'OBSTRUCTED',
                    checkpoint=d['checkpoint'] if d else None,steps=d['steps'] if d else start,time=d['time'] if d else None,
                    failure=d['failure'] if d else 'missing report',requested_steps=count,dt=dt,refined=key=='refined')
                state['runs'][label]=record;save(status,state)
                if record['status']!='PASS':
                    state.update(status='OBSTRUCTED',first_limiting_gate=label);save(status,state);failed=True;break
                assert d['steps']==start+count and d['initialization_lift_calls_iterations']==[0,0]
                last[key]=record
            if failed:break
        if failed:break
        for order,(points,weights) in rules.items():
            current={}
            for key in ('base','fine','refined'):
                record=last[key];d=load(record['report']);label=f'{cadence:02d}-{key}-q{order}';folder=out/label;folder.mkdir()
                packet=dict(operation='open_checks',inner_stress_model='Linear',native_report=record['report'],
                    native_report_sha256=record['report_sha256'],native_checkpoint=d['checkpoint'],
                    checkpoint_sha256=sha(Path(d['checkpoint'])/'State.bin'),physical_observations=points,
                    assembly=False,reuse=False,weak=False)
                save(folder/'packet.json',packet)
                assert execute(root,'revised-observe-'+label,[a.observer,record['xml'],folder/'packet.json',folder/'native.json'],folder,0,600,extra_env=env)
                current[key]=load(folder/'native.json');state['observations'][label]=dict(path=str(folder/'native.json'),sha256=sha(folder/'native.json'))
            for name,x,y in [('time','base','fine'),('normal','fine','refined')]:
                row=comparison(current[x],current[y],weights,c.K,eps);label=f'{cadence:02d}-{name}-q{order}';state['comparisons'][label]=row;save(status,state)
                if row['status']!='PASS':
                    state.update(status='UNRESOLVED',first_limiting_gate=label);failed=True;break
            del current
            if failed:break
        if failed:break
        assert all(sha(path)==value for path,value in immutable.items()),'revised design changed prior evidence'
        state['last_common_time']=target;save(status,state)
    else:
        state.update(status='PASS',completed_time=8.,development='apply separate residual/station-profile window criteria; completion alone is not developed')
        save(status,state)
    print(state['status'],state['first_limiting_gate'])
    raise SystemExit(0 if state['status']=='PASS' else 2)


if __name__=='__main__':main()
