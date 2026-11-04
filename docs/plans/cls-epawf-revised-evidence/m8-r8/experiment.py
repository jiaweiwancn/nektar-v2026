#!/usr/bin/env python3
"""Bind the R7 terminal ray and run the focused J25 check on Linux."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def bind(r7, out, h, r7_tool):
    ledger = h.load(r7/'source/docs/plans/cls-epawf-revised-evidence/m8-r7/decision-ledger.json') if (
        r7/'source/docs/plans/cls-epawf-revised-evidence/m8-r7/decision-ledger.json').exists() else None
    # The producing R7 checkout predates its compact decision ledger.
    if ledger is None:
        ledger = h.load(out.parent/'reviewed-r7-ledger.json')
    p = ledger['pilot']
    receipt = h.load(h.bound(p['receipt']))
    collection = h.load(h.bound(receipt['durable_collection']['file']))
    events = r7_tool.events(collection)
    base = collection['last_snapshot']
    matrix = p['latest_complete_matrix']['matrix']
    j, = [x for x in collection['jacobians'] if x['iteration']==25 and x['status']=='complete']
    assert base['iteration']==25 and base['snapshot']==p['terminal_snapshot']['snapshot']
    assert j['event']['matrix']==matrix and not j['event']['details']['loaded_bound_matrix']
    assert base['history']==ledger['target']['original_h0'] and not base['history_advanced']
    d, = [x for x in events if x['kind']=='r5_direction' and x['iteration']==25]
    cut, = [x for x in events if x['kind']=='r6_cut_audit' and x['name']=='iterate25.r7_full']
    assert cut['status']=='passed_original_cut_scan' and max(cut['worst_native_violation'],cut['worst_extended_violation'])<=5e-12
    history, state, direction = (h.vector(pin) for pin in (base['history'],base['state'],d['step']))
    old = [x for x in events if x['kind']=='r5_trial' and x['iteration']==25]
    assert [x['fraction'] for x in old]==[1.,.5,.25,.125]
    for trial in old:
        assert h.vector(trial['state'])==[a+trial['fraction']*b for a,b in zip(state,direction)]
        h.bound(trial['residual'])
        assert trial['classification']!='accepted_comparison'
    for pin in (base['snapshot'],base['residual'],base['constraints'],matrix,d['action'],d['step']):
        h.bound(pin)
    out.mkdir()
    h.save(out/'binding.json',dict(base=base,matrix=matrix,producer=p['receipt'],direction=d,
        cut= {k:cut[k] for k in ('status','worst_native_violation','worst_extended_violation')},
        trial_records=old,assembly=j['event']['details']['assembly'],manifest=receipt['manifest']))
    h.save(out/'packet.json',dict(state=state,history=history,alpha=400,time=.0025,only_external_direction=True,
        r5_matrix_binding=dict(state_sha256=base['state']['sha256'],history_sha256=base['history']['sha256'],
            residual_sha256=base['residual']['sha256'],matrix=matrix,
            assembly_json=[json.dumps(a) for a in j['event']['details']['assembly']])))
    source_probe=r7/'cases-001/interior/probe.xml'
    (out/'probe.xml').write_bytes(source_probe.read_bytes())


def action_check(r7, binding, h):
    tests = r7/'source/solvers/IncNavierStokesSolver/Tests/RANS'
    sys.path.insert(0,str(tests))
    import numpy as np
    from validate_cls_epawf_m8_r3_actions import validate_case
    pin=h.load(binding/'binding.json')
    runtime=h.load(r7/'runtime-fingerprint.json')
    exe=h.pin(r7/'build/RANSNonpressureActionProbe')
    assert exe['sha256']==runtime['runtime_artifacts']['build/RANSNonpressureActionProbe']
    raw=binding/'action.json'
    cmd=['taskset','-c','49',exe['path'],str(binding/'probe.xml'),str(binding/'packet.json'),
         str(raw),str(h.bound(pin['direction']['step']))]
    env=dict(os.environ,LD_LIBRARY_PATH=str(r7/'build/geometry-lib'))
    with (binding/'action.log').open('x') as log:
        subprocess.run(cmd,cwd=binding,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=120,check=True)
    data=h.load(raw)
    assert h.pin(data['matrix_file'])==pin['matrix'] and np.asarray(data['directions'][0]['stage_action'],
        dtype='<f8').tobytes()==h.bound(pin['direction']['action']).read_bytes()
    result=validate_case(data,r7/'source')
    h.save(binding/'action-check.json',dict(success=result['success'],binding=h.pin(binding/'binding.json'),
        probe=exe,action=h.pin(raw),command=cmd,result=result))
    if not result['success']:
        raise RuntimeError('J25 independent action gate failed')


def prepare(epoch, r7, binding, out, h):
    check=h.load(binding/'action-check.json')
    assert check['success'] and check['binding']==h.pin(binding/'binding.json')
    data=h.load(binding/'binding.json')
    base=data['base']
    previous=h.load(h.bound(data['manifest']))
    old=h.bound(previous['current_target_xml'])
    exe=h.pin(epoch/'build/IncNavierStokesSolver')
    runtime_path=epoch/'runtime-fingerprint.json'
    runtime=h.load(runtime_path)
    assert runtime['success'] and runtime['runtime_artifacts']['build/IncNavierStokesSolver']==exe['sha256']
    tree=ET.parse(old)
    info=tree.find('.//SOLVERINFO')
    selectors={row.get('PROPERTY'):row for row in info}
    for key in ('RANSWallR7Policy','RANSWallR7Radius','RANSWallR7InteriorTarget','RANSWallR5Updates'):
        if key in selectors:
            info.remove(selectors[key])
    selectors['RANSWallStageTest'].set('VALUE','DirectionExperimentR5')
    selectors['RANSWallR3PredictorFile'].set('VALUE',str(out/'predictor.json'))
    selectors['RANSWallR3ProgressFile'].set('VALUE',str(out/'case.r3-progress.jsonl'))
    ET.SubElement(info,'I',PROPERTY='RANSWallR8Experiment',VALUE='On')
    out.mkdir()
    xml=out/'case.xml'
    tree.write(xml,encoding='unicode')
    assert h.predictor.protected_xml(xml)==previous['protected_physical_xml_sha256']
    d=data['direction']
    h.save(out/'predictor.json',dict(previous,role='retained_r7_terminal_iterate_only',
        state=base['state'],history=base['history'],producing_source_commit=runtime['source_commit'],
        current_executable=exe,current_target_xml=h.pin(xml),current_runtime=h.pin(runtime_path),
        origin=dict(alpha=400,time=.0025,converged=False,inner_stress_model='SaturatedR2',
            iteration=25,snapshot=base['snapshot'],residual=base['residual'],constraints=base['constraints']),
        frozen_action=dict(matrix=data['matrix'],state_sha256=base['state']['sha256'],
            history_sha256=base['history']['sha256'],residual_sha256=base['residual']['sha256'],
            assembly_json=json.dumps(data['assembly'])),
        r8_ray=dict(step=d['step'],action=d['action'],trials=[
            {k:x[k] for k in ('fraction','state','residual')} for x in data['trial_records']])))


def summary(epoch, out, h):
    receipt=h.load(epoch/'local-run/receipt.json')
    manifest=h.load(h.bound(receipt['manifest']))
    collection=h.load(h.bound(receipt['durable_collection']['file']))
    initial,=collection['completed_iterates']
    matrix,=collection['jacobians']
    assert all(initial[k]['sha256']==(manifest[k] if k in ('state','history') else manifest['origin'][k])['sha256']
               for k in ('state','history','residual','constraints'))
    assert matrix['event']['matrix']['sha256']==manifest['frozen_action']['matrix']['sha256']
    assert matrix['event']['details']['loaded_bound_matrix']
    prefix=collection['progress']
    import hashlib
    payload=Path(prefix['path']).read_bytes()[:prefix['complete_prefix_bytes']]
    assert hashlib.sha256(payload).hexdigest()==prefix['complete_prefix_sha256']
    events=[json.loads(l) for l in payload.splitlines()]
    trials=[x for x in events if x['kind']=='r8_trial']
    work=next(x for x in events if x['kind']=='r8_work')
    assert work['native_trial_states']==len(trials)<=20
    assert collection['terminal']['terminal']=='r8_fixed_comparison_completed'
    report=dict(receipt=h.pin(epoch/'local-run/receipt.json'),initial=initial['snapshot'],
        J25=matrix['event']['matrix'],results=[x for x in events if x['kind'].startswith('r8_')],
        cuts=[x for x in events if x['kind']=='r6_cut_audit'])
    h.save(out/'comparison.json',report)
    for t in report['results']:
        k=t['kind']
        if k=='r8_trial':
            print(k,t['name'],t['fraction'],t['classification'],max(t.get('actual',[0])),
                  t.get('wall',{}).get('minimum'),t.get('linear_agreement'),t.get('quadratic_agreement'))
        if k in ('r8_radius_selection','r8_curvature_outcome','r8_work','r8_direction_failure','r8_components'):
            print({x:y for x,y in t.items() if x not in ('schema','policy','alpha','time','action_policy','sequence','solve_number')})

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('bind','action','prepare','summary'))
    p.add_argument('--r7',type=Path,required=True)
    p.add_argument('--binding',type=Path,required=True)
    p.add_argument('--epoch',type=Path)
    p.add_argument('--out',type=Path)
    a=p.parse_args()
    if platform.system()!='Linux':
        raise ValueError('Use the designated Linux host')
    r7=a.r7.resolve()
    evidence=r7/'source/docs/plans/cls-epawf-revised-evidence'
    h=module(evidence/'m8-r4/prepare.py','r4_prepare')
    if a.mode=='bind':
        # The checked-in R7 utility may postdate the producing solver.
        r7_tool=module(a.binding.parent/'r7-experiment.py','r7_experiment')
        bind(r7,a.binding,h,r7_tool)
    elif a.mode=='action':
        action_check(r7,a.binding,h)
    else:
        h=module(a.epoch/'source/docs/plans/cls-epawf-revised-evidence/m8-r4/prepare.py','r8_prepare')
        if a.mode=='prepare':
            prepare(a.epoch,r7,a.binding,a.out,h)
        else:
            summary(a.epoch,a.out,h)


if __name__=='__main__':
    main()
