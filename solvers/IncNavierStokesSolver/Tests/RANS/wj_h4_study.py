#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded H4 import, native controls and retrospective physical inventories."""
import argparse
import os
from pathlib import Path
import subprocess
import time
import numpy as np
from numpy.testing import assert_allclose
from run_wj_h4 import load,save,sha,THREADS
from wj_h3_analysis import actual,identity,state_context
from wj_h4_points import inventory,scan,owners,evaluate,owner_point,vector_hash,aligned
from validate_wj1 import Probe

PARENT=actual('/media/jay/新加卷/nektar-v2026-runs/wj-h3-20261112')
RUNTIME=actual('/media/jay/新加卷/nektar-v2026-runs/wj-h2-20261111/runtime-v4/IncNavierStokesSolver')
OBSERVER=actual('/media/jay/新加卷/nektar-v2026-runs/wj-h2-20261111/observer-v6/RANSWallResolutionProbe')
KERNEL=actual('/media/jay/新加卷/nektar-v2026-runs/wj-h2-20261111/build-v1/WJKernelProbe')
SOURCE=Path(__file__).resolve().parents[4]
CASES={'base-x1-half':dict(refined=False,nx=1,dt=.005,steps=16,size=1062),
       'refined-x1':dict(refined=True,nx=1,dt=.01,steps=8,size=1966),
       'refined-x1-half':dict(refined=True,nx=1,dt=.005,steps=16,size=1966),
       'base-x3-half':dict(refined=False,nx=3,dt=.005,steps=16,size=2630),
       'local-x3':dict(refined=False,nx=3,dt=.01,steps=8,size=2878,extra_physical_z=[.8818157002137355,1.1947077409162703]),
       'local-x3-half':dict(refined=False,nx=3,dt=.005,steps=16,size=2878,extra_physical_z=[.8818157002137355,1.1947077409162703])}


def runtime_check(root):
    inputs=load(SOURCE/'docs/plans/wj-h2-linux/native-source-inputs.json')
    for f in inputs:assert sha(SOURCE/f['path'])==f['sha256'],f['path']
    libs=load(SOURCE/'docs/plans/wj-h2-linux/loaded-runtime.json')
    for f in libs['loaded_files']:assert sha(actual(f['path']))==f['sha256'],f['path']
    assert sha(RUNTIME)=='9ae37a1aba62755b10ad99c9484c498cef9d927a75fb950f00344f4ffc5d06fe'
    assert sha(OBSERVER)=='6f35b837f0dd13eb155e0fe34f0dc670c609ee6cbe739609a5de288dd49ed06d'
    save(root/'receipts/runtime-equivalence.json',dict(status='PASS',native_inputs=len(inputs),closure_files=len(libs['loaded_files']),
        solver=identity(RUNTIME),observer=identity(OBSERVER),kernel=identity(KERNEL),source_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=SOURCE,text=True).strip()))


def bind(report,step=None,rejected=False):
    state,history,alpha,t,p=state_context(report,step,rejected)
    assert vector_hash(state)==p['state_sha256']
    if p['callback_history']:assert vector_hash(p['callback_history'])==p['stage_history_sha256']
    if not rejected:aligned(report['steps'] if step is None else step,report['dt'],t)
    return state,dict(time=t,alpha=alpha,packet=p,state_sha256=p['state_sha256'],callback_history_sha256=p['stage_history_sha256'],
        native_status='rejected candidate' if rejected else 'accepted under original native inventory',step=step)


def assess_report(root,label,rp,inv,first_step=None):
    """Every accepted state, including odd half steps; no callback substitution."""
    folder=root/'scans'/label;folder.mkdir(parents=True,exist_ok=False);report=load(rp);result=[]
    start=report['steps']+1-len(report['accepted_states']);earliest=None
    selections=[(step,False) for step in range(start,report['steps']+1) if first_step is None or step>=first_step]
    if report['failure'] and report['physical_range_events']:selections.append((None,True))
    probe=Probe(str(KERNEL))
    for step,rejected in selections:
        state,context=bind(report,step,rejected);packet=scan(report,state,inv)
        # Independently replay every state's limiting outer action, including
        # accepted states newly found inadmissible by the expanded inventory.
        w=packet['minimum']['row'];G=np.array(w['G']);h=np.outer(w['principal']['polarization'],w['principal']['xi'])
        c=probe.closure(w['k'],w['epsilon'],G);d=probe.closure(w['k'],w['epsilon'],G,(0,0,h))
        replay=float(-np.sum(h*d['R'])/c['scalar'][1]);assert_allclose(replay,w['principal']['in_plane'],rtol=2e-9,atol=2e-11)
        packet['compiled_in_plane_replay']=replay
        if packet['status']!='PASS' and not rejected and earliest is None:earliest=context['time']
        packet.update(context,producer=identity(rp),prefix_qualified=earliest is None,
            qualification='PASS' if packet['status']=='PASS' and earliest is None else 'FAIL' if packet['status']!='PASS' else 'UNQUALIFIED_AFTER_EARLIER_EXPANDED_FAILURE')
        path=folder/(f'step-{step}.json' if not rejected else 'rejected.json');save(path,packet)
        result.append(dict(path=str(path),sha256=sha(path),step=step,time=context['time'],native_status=context['native_status'],
            expanded_status=packet['status'],prefix_qualified=packet['prefix_qualified'],minimum=packet['minimum'],failures=len(packet['failures']),
            state_sha256=context['state_sha256'],callback_history_sha256=context['callback_history_sha256']))
    probe.close();summary=dict(producer=identity(rp),inventory_sha256=inv['sha256'],states=result,earliest_extra_failure=earliest,
        retrospective=True,scope='no claim to an earlier native precommit stop')
    save(folder/'index.json',summary);return summary


def preflight(root):
    runtime_check(root)
    published=load(SOURCE/'docs/plans/wj-h3-linux/runs.json');wanted=('B-zero-x1','B-zero-x3','A-dt02','A-dt01')
    selected={k:next(v['producer'] for v in published if v['case']==k) for k in wanted}
    inv_reports=[]
    for name in (*wanted,'B-frozen-x1','B-frozen-x3'):
        record=next(v['producer'] for v in published if v['case']==name);rp=actual(record['path']);assert sha(rp)==record['sha256'];inv_reports.append(load(rp))
    inv=inventory(inv_reports);save(root/'inventory-v0.json',inv)
    result={}
    for name,record in selected.items():
        rp=actual(record['path']);assert sha(rp)==record['sha256'];result[name]=assess_report(root,'retained-'+name,rp,inv)
    # Check representative existing native samples on their exact declared owner.
    entry=load(PARENT/'native-observations-v1/index.json')['B-zero-x3'];p=actual(entry['native']['path']);assert sha(p)==entry['native']['sha256']
    native=load(p)['weak_checks'];report=load(selected['B-zero-x3']['path'])
    errors=[]
    for row in native['samples'][::191]:
        el=next(v for v in report['fields'][0]['mesh'] if v['geometry']==row['geometry']);o=owner_point(report,el,row['physical'])
        v=evaluate(report,report['state'],row['physical'],o);assert_allclose(v['G'],row['G'],rtol=3e-10,atol=2e-12)
        assert_allclose(v['k'],row['k'],rtol=3e-12,atol=2e-13);errors.append(float(np.max(abs(np.array(v['G'])-row['G']))))
    missed=any(result[k]['earliest_extra_failure'] is not None for k in ('B-zero-x1','B-zero-x3'))
    # A .01 is also an affected native-accepted trajectory if expanded points
    # reveal an omitted negative; A's already-rejected candidate is not a gap.
    missed=missed or result['A-dt01']['earliest_extra_failure'] is not None or result['A-dt02']['earliest_extra_failure'] is not None
    save(root/'preflight.json',dict(status='EXPANDED_GAP' if missed else 'PASS',strict_path_required=missed,
        cases=result,inventory=identity(root/'inventory-v0.json'),native_sample_max_error=max(errors),native_sample_count=len(errors)))
    # The two fixed splits are admitted only as a retained-data normal question.
    # A sign change at its normal face and a small native B margin provide that
    # question; resource admission remains separate and is updated after pilots.
    question=any(s['minimum']['margin']<.01 for k in wanted for s in result[k]['states'])
    save(root/'optional-activation.json',dict(activated=question,resource_admission='PENDING_MEASURED_CORE_COST',
        physical_splits=[.8818157002137355,1.1947077409162703],
        evidence='retained limiting B neighborhood and A normal-face sign/timestep sensitivity',
        preflight_sha256=sha(root/'preflight.json'),no_fresh_evolution_before_choice=True))
    print('Retained expanded preflight',missed,'points',len(inv['points']),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('preflight',),required=True)
    a=p.parse_args();preflight(a.root)
