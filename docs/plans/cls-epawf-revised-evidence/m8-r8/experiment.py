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

def pilot_prepare(epoch, binding, local, out, h):
    comparison=h.load(local/'comparison.json')
    cost=h.load(h.bound(comparison['receipt']))
    records=comparison['results']
    selected,=[x for x in records if x['kind']=='r8_selected']
    trials=[x for x in records if x['kind']=='r8_trial']
    ordinary,=[x for x in trials if x['name']=='r8_ordinary0']
    chosen,=[x for x in trials if x['name']==selected['name']]
    assert chosen['name']=='r8_curvature' and chosen['classification']=='accepted'
    assert max(chosen['actual'])<max(ordinary['actual']) and chosen['wall']['minimum']>ordinary['wall']['minimum']
    assert chosen['box']['maximum']==0. and ordinary['box']['maximum']==0.
    action=h.load(binding/'action-check.json')
    assert action['success'] and action['binding']==h.pin(binding/'binding.json')
    old=h.load(h.bound(cost['manifest']))
    assert comparison['J25']['sha256']==old['frozen_action']['matrix']['sha256']
    runtime_path=epoch/'runtime-fingerprint.json'
    runtime=h.load(runtime_path)
    exe=h.pin(epoch/'build/IncNavierStokesSolver')
    assert runtime['success'] and runtime['runtime_artifacts']['build/IncNavierStokesSolver']==exe['sha256']
    out.mkdir()
    tree=ET.parse(h.bound(old['current_target_xml']))
    info=tree.find('.//SOLVERINFO')
    selectors={row.get('PROPERTY'):row for row in info}
    selectors['RANSWallR8Experiment'].set('VALUE','Pilot')
    selectors['RANSWallStageTest'].set('VALUE','ExactTargetR3')
    selectors['RANSWallR3PredictorFile'].set('VALUE',str(out/'predictor.json'))
    selectors['RANSWallR3ProgressFile'].set('VALUE',str(out/'case.r3-progress.jsonl'))
    ET.SubElement(info,'I',PROPERTY='RANSWallR5Updates',VALUE='10')
    xml=out/'case.xml'
    tree.write(xml,encoding='unicode')
    assert h.predictor.protected_xml(xml)==old['protected_physical_xml_sha256']
    replay=old['r8_ray']
    fit,=[t for t in replay['trials'] if t['fraction']==.125]
    holdout={k:ordinary[k] for k in ('fraction','state','residual')}
    h.save(out/'predictor.json',dict(old,current_executable=exe,
        current_target_xml=h.pin(xml),current_runtime=h.pin(runtime_path),producing_source_commit=runtime['source_commit'],
        r8_ray=dict(replay,trials=[fit,holdout])))
    h.save(out/'gate.json',dict(schema=1,source_commit=runtime['source_commit'],
        executable_sha256=exe['sha256'],absolute_action_tolerance=1e-7,required_state_action_gates_met=True,
        active_dependencies_covered=True,useful_joint_native_descent=True,base_state=old['state'],
        maximum_updates=10,local_evidence=comparison['receipt'],local_source_commit=cost['source_commit'],
        local_executable=cost['executable'],current_first_matrix=comparison['J25'],
        action=h.pin(binding/'action-check.json'),comparison=h.pin(local/'comparison.json'),
        selected_fraction=chosen['fraction'],selected_state=chosen['state'],selected_residual=chosen['residual']))
    (out/'budget.txt').write_text(
        "One pilot with ten accepted updates, nine reachable native trials and one current direction per base. "
        "The first base reuses its verified current J25, direction and two samples. A fresh matrix and live samples "
        "are required after acceptance. Cooperative stop at a completed state after 300 solver seconds. "
        "The external cap is 450 native seconds. No restart.\n")


def verify_pilot(epoch, out, h):
    receipt=h.load(epoch/'pilot-run/receipt.json')
    collection=h.load(h.bound(receipt['durable_collection']['file']))
    gate=h.load(out/'gate.json')
    saved=collection['completed_iterates']
    assert receipt['maximum_seconds']==450 and collection['protocol']['maximum_Newton_updates']==10
    assert 2<=len(saved)<=11
    assert saved[0]['state']['sha256']==gate['base_state']['sha256']
    assert saved[1]['state']['sha256']==gate['selected_state']['sha256']
    assert saved[1]['residual']['sha256']==gate['selected_residual']['sha256']
    for row in collection['jacobians']:
        if row['status']=='complete':
            bound=row['event']['details']['loaded_bound_matrix']
            assert bound==(row['iteration']==0)
            if bound:
                assert row['event']['matrix']['sha256']==gate['current_first_matrix']['sha256']
    import hashlib
    prefix=collection['progress']
    payload=Path(prefix['path']).read_bytes()[:prefix['complete_prefix_bytes']]
    assert hashlib.sha256(payload).hexdigest()==prefix['complete_prefix_sha256']
    progress=[json.loads(l) for l in payload.splitlines()]
    trials=[x for x in progress if x['kind']=='r8_trial']
    updates={x['details']['iteration']+1:x['details'] for x in progress if x['kind']=='newton_update'}
    maximum=max((sum(t['iteration']==i for t in trials) for i in range(len(saved))),default=0)
    assert maximum<=9
    for left,right in zip(saved,saved[1:]):
        update=updates[right['iteration']]
        assert max(right['scaled_block_norms'])<(1.-1e-4*update['accepted_fraction'])*max(left['scaled_block_norms'])
        assert update['box']['maximum']==0. and right['scalar_domain']['admissible']
        assert right['scalar_domain']['r4_wall_trace']['minimum']>=0. and not right['history_advanced']
    cut_rows=[x for x in progress if x['kind']=='r6_cut_audit' and x['status']=='passed_original_cut_scan']
    largest_cut=max((max(x['worst_native_violation'],x['worst_extended_violation']) for x in cut_rows),default=0.)
    assert largest_cut<=5e-12
    h.save(out/'verification.json',dict(success=True,receipt=h.pin(epoch/'pilot-run/receipt.json'),
        first_step_reproduced=True,maximum_new_trials=maximum,trials=len(trials),accepted_updates=len(saved)-1,
        terminal=collection['terminal'],largest_original_row_violation=largest_cut))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('bind','action','prepare','summary','pilot','verify-pilot'))
    p.add_argument('--r7',type=Path,required=True)
    p.add_argument('--binding',type=Path,required=True)
    p.add_argument('--epoch',type=Path)
    p.add_argument('--out',type=Path)
    p.add_argument('--local',type=Path)
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
        elif a.mode=='summary':
            summary(a.epoch,a.out,h)
        elif a.mode=='pilot':
            pilot_prepare(a.epoch,a.binding,a.local,a.out,h)
        else:
            verify_pilot(a.epoch,a.out,h)


if __name__=='__main__':
    main()
