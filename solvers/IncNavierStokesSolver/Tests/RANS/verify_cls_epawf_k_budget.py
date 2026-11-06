#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Close the actual frozen discrete k increment budget, not a continuum surrogate.

Two offline mass-dual functionals convert free-row loads into physical global
and first-cell integrals. Production equations and histories are untouched.
"""
import argparse
from pathlib import Path
import json
import numpy as np
from run_cls_epawf_picard import load,require,write_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--campaign',type=Path,required=True);p.add_argument('--prior',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.campaign=a.campaign.resolve();a.prior=a.prior.resolve()
    spec=load(a.campaign/'analysis-final/contact0025/result.json')
    mass=np.array(spec['free_k_mass']);rows=np.array(spec['free_k_rows'],dtype=int);offset=int(spec['offsets'][2])
    functionals=np.array([spec['k_integral_load'],spec['near_k_integral_load']])[:,rows-offset]
    dual=np.linalg.solve(mass.T,functionals.T).T
    prior=load(a.prior/'fixed-horizon/bound-0/case.evolution.json')
    attempts={round(v['time']/.000625):v for v in prior['pimple_attempts'] if v['time']<=.0025+1e-13}
    for line in (a.campaign/'main/steps.jsonl').read_text().splitlines():
        v=json.loads(line);attempts[v['global_step']]=v['attempt']
    progress=load(a.campaign/'main/progress.json');report=load(Path(progress['segment'])/'case.evolution.json')
    attempts[int(report['steps'])+1]=report['pimple_attempts'][-1]
    result=dict(schema=1,offline_mass_inverse_applications=2,production_governing_solves=0,steps=[],windows=[])
    for step,row in sorted(attempts.items()):
        packet=load(Path(row['k_snapshot']));n=len(rows)
        matrix=np.frombuffer(Path(packet['matrix_file']).read_bytes(),dtype='<f8').reshape(n,n)
        raw=np.array(packet['raw_increment']);r=np.array(packet['transport_residual']);h=1/packet['alpha']
        candidate=np.array(packet['corrected']);before=np.array(packet['pre_k']);unbounded=np.array(packet['raw'])
        correction=(candidate-unbounded)[rows];net=(candidate-before)[rows]
        history=np.array(packet['history']);k_before=before[offset:offset+len(spec['k_integral_load'])]
        k_history=history[offset:offset+len(spec['k_integral_load'])]
        lraw=matrix@raw-(mass@raw)/h
        input_load=-h*(dual@r);implicit_increment=-h*(dual@lraw)
        repair=functionals@correction;actual=functionals@net
        defect=actual-input_load-implicit_increment-repair
        require(np.max(np.abs(defect))<1e-12,'physical frozen increment budget did not close')
        result['steps'].append(dict(step=step,accepted=step<=report['steps'],time=packet['time'],
            physical_global_near_increment=actual.tolist(),frozen_input_load_contribution=input_load.tolist(),
            frozen_L_increment_contribution=implicit_increment.tolist(),bound_contribution=repair.tolist(),
            closure_defect=defect.tolist(),k_pre_minus_native_history_max=float(np.max(np.abs(k_before-k_history)))))
    for lo,hi in ((0,4),(4,8),(8,12)):
        rows_window=[v for v in result['steps'] if lo<v['step']<=hi and v['accepted']]
        result['windows'].append(dict(start=lo*.000625,end=hi*.000625,
            **{k:np.sum([v[k] for v in rows_window],axis=0).tolist() for k in
                ('physical_global_near_increment','frozen_input_load_contribution','frozen_L_increment_contribution','bound_contribution','closure_defect')}))
    result['status']='passed';write_json(a.output,result);print(a.output)


if __name__=='__main__':main()
