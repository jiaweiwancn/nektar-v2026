#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Executed common-time recovery comparison; no inferred future qualification."""
import argparse
from pathlib import Path
import numpy as np
from cao_tamura_h2d_state import load
from wp2_all_state import Rows,US
from wp2_qualification import quarter_budget_differences
from run_cls_epawf_picard import require,sha,state_hash
from continue_cls_epawf_pimple import atomic


def states(progress):
    p=load(progress);out={};context=None
    for observation in p['observations']:
        d=load(Path(observation['report']));context=d
        for s,m in zip(d['accepted_states'],d['wp2_per_state_monitoring']):
            key=round(m['time'],12)
            if key in out:require(state_hash(out[key][1])==state_hash(s),'recovery restart changed state')
            else:out[key]=(m['time'],s)
    require(out,'no executed observations');return p,out,context


def main():
    p=argparse.ArgumentParser()
    for n in ['root','left','right','output']:p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();left_hash,right_hash=sha(a.left),sha(a.right)
    left,ls,lc=states(a.left);right,rs,rc=states(a.right);rows=Rows(a.root/'geometry-rows-selected/result.json')
    require(left['model']==right['model'],'recovery models differ')
    require(state_hash(ls[0.][1])==state_hash(rs[0.][1]),'recovery timesteps did not start from the exact same pulse')
    for key in ['fixed_reference','fixed_defect']:require(state_hash(lc[key])==state_hash(rc[key]),'recovery '+key+' differs')
    result=[]
    for key in sorted(set(ls)&set(rs)):
        lt,l=ls[key];rt,r=rs[key];require(abs(lt-rt)<2e-14,'common recovery time exceeds arithmetic tolerance')
        lv=rows.tuple(l,lc);rv=rows.tuple(r,rc);cold=rows.tuple(lc['fixed_reference'],lc)
        check=quarter_budget_differences(dict(lc['wp2_common_samples'],values=lv['values']),dict(rc['wp2_common_samples'],values=rv['values']),US,lc['parameters'][5],.0065)
        response={}
        for col,name in [(0,'U'),(1,'W'),(2,'k'),(3,'outer_epsilon'),(4,'pi'),(5,'mechanical_pressure')]:
            mask=~rows.inner if col==3 else np.ones(len(rows.inner),dtype=bool);w=rows.weight[mask]
            den=float(np.sqrt(np.dot(w,(rv['values'][mask,col]-cold['values'][mask,col])**2)))
            error=float(np.sqrt(np.dot(w,(lv['values'][mask,col]-rv['values'][mask,col])**2)))
            response[name]=dict(dimensional_difference_L2=error,fine_departure_from_native_cold_L2=den,
                                fraction_of_fine_response=error/den if den>1e-13 else None,comparison_floor=1e-13)
        result.append(dict(left_actual_time=lt,right_actual_time=rt,left_state_sha256=state_hash(l),right_state_sha256=state_hash(r),principal=check,response=response))
    require(sha(a.left)==left_hash and sha(a.right)==right_hash,'live progress changed during read-only comparison; retry a stable completed boundary')
    atomic(a.output,dict(version='WP2RecoveryCommonTimeControlV1',left=dict(path=str(a.left),sha256=left_hash,dt=left['dt']),
        right=dict(path=str(a.right),sha256=right_hash,dt=right['dt']),same_seed_sha256=state_hash(ls[0.][1]),results=result,
        quarter_budgets_passed_on_executed_common_times=all(x['principal']['passed'] for x in result),
        scope='actual native common times only; same exact perturbation; absolute pressures and complete outer epsilon; no temporal order or recovery claim'))


if __name__=='__main__':main()
