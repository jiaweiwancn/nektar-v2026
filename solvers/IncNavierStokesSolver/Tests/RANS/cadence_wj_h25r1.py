#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Exact dense/sparse/uninterrupted control parity, with fresh audit provenance."""
import argparse
from pathlib import Path
from run_wj_h4 import load
from wj_h15_observer import identity,verify,atomic,vector_hash
from wj_h17_observer import native_load,equal,independent_parity


def run(root):
    base=root/'resumption-r1';wide=base/'science/C-wide-v1/monitor';records=[]
    labels=('C-sparse-restart-v1','C-dense-restart-v1')
    for label in labels:
        d=base/'science'/label;q=load(d/'qualification.json')
        assert q['status']=='PASS_UNPERTURBED_RESTART_PILOT'
        log=q['bounded_stdout'];assert log['status']=='PASS_BOUNDED_STDOUT'and log['rotations']>0
        for pin in log['retained']:verify(pin)
        for step in range(9,13):
            p=d/f'monitor/step-{step}.report.json';r=native_load(p);old=native_load(wide/f'step-{step}.report.json')
            for k in ('state','fixed_reference','fixed_defect','native_history','native_history_times','time','steps','dt','signature','last_stage_history','last_stage_alpha','last_stage_time'):
                equal(r[k],old[k],f'{label} {step} {k}')
            both_full=r['native_full_audit_status']==old['native_full_audit_status']=='PASS'
            if both_full:
                equal(r['native_hard_audit'],old['native_hard_audit'],'fresh native complete physical map')
                independent_parity(load(d/f'monitor/step-{step}.result.json'),load(d/f'monitor/step-{step}.inventory.json'),
                    load(wide/f'step-{step}.result.json'),load(wide/f'step-{step}.inventory.json'))
            records.append(dict(label=label,step=step,time=r['time'],report=identity(p),oracle=identity(wide/f'step-{step}.report.json'),
                state=vector_hash(r['state']),fresh_native_and_independent_map_exact=both_full))
    dense=base/'science/C-dense-restart-v1/monitor';moving=[]
    for step in range(9,13):
        result=load(dense/f'step-{step}.result.json');assert result['status']=='PASS'
        interface=next(v for v in result['summary']['interface_budget']['rows']if v['x']==0.)
        r=native_load(dense/f'step-{step}.report.json');ready=load(dense/f'step-{step}.ready.json')
        structure=load(verify(ready['structure']));offsets=structure['offsets']
        stations=sorted({float(e[end][0])for e in structure['fields'][0]['mesh']for end in ('lower','upper')})[1:-1]
        traces=[next(v for v in result['summary']['interface_budget']['rows']if v['x']==x)for x in stations]
        moving.append(dict(step=step,physical_interface_x=0.,K=interface['K'],cutoff=interface['rc'],
            fixed_interior_stations=stations,interior_K=[v['K']for v in traces],interior_cutoff=[v['rc']for v in traces],
            k_coefficients=vector_hash(r['state'][offsets[2]:offsets[3]]),
            epsilon_coefficients=vector_hash(r['state'][offsets[3]:offsets[4]])))
    for field in ('k_coefficients','epsilon_coefficients'):
        assert len({r[field]for r in moving})>1,'control did not change '+field
    for field in ('interior_K','interior_cutoff'):
        assert len({tuple(r[field])for r in moving})>1,'control did not change '+field
    out=base/'science/cadence-comparison-v1.json';assert not out.exists()
    atomic(out,dict(status='PASS_UNPERTURBED_CADENCE_AND_RESTART_ONLY',records=records,moving_fields=moving,
        scientific_strides=dict(wide=2,sparse=4,dense=1),recovery_stride=1,
        repaired_diagnostic='v1 wrongly demanded changing K at the prescribed, fixed inlet x=-10; v2 retains the change requirement on fixed interior mesh-vertex traces, not a moving witness',
        limitation='finite prescribed seed rejected before evolution; no production cadence or full-duration admission'))
    print('exact cadence/restart/maps PASS; changing K/epsilon/cutoff verified',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);run(p.parse_args().root)
