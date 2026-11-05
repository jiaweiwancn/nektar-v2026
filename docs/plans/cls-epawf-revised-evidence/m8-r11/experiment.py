#!/usr/bin/env python3
"""Bind and collect the fixed four-policy R11 experiment on authoritative Linux."""
import argparse
import importlib.util
import json
import hashlib
import struct
from pathlib import Path
import platform
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode=True


def module(path, name):
    spec=importlib.util.spec_from_file_location(name,path)
    obj=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def selected_trial(policy, trials):
    """Resolve the selector's exact chosen payload, including multiple accepted trials."""
    if not policy['accepted']:
        return None
    chosen=[t for t in trials if t['comparison']==policy['name'] and
            t['fraction']==policy['fraction'] and
            t['state']['sha256']==policy['selected_state']['sha256'] and
            t['residual']['sha256']==policy['selected_residual']['sha256']]
    if not chosen or any(t['classification']!='accepted' for t in chosen):
        raise ValueError('Chosen native trial is missing or rejected')
    return chosen[0]  # Exact aliases have identical fraction, state, residual and acceptance.


def prepare(e, out, h, evidence):
    previous=module(evidence/'m8-r10/experiment.py','r10_prepare')
    previous.prepare(e,out,h,evidence)
    r10=h.load(evidence/'m8-r10/decision-ledger.json')
    manifest=h.load(out/'predictor.json')
    tree=ET.parse(out/'case.xml')
    info=tree.find('.//SOLVERINFO')
    selectors={r.get('PROPERTY'):r for r in info}
    info.remove(selectors['RANSWallR10Experiment'])
    ET.SubElement(info,'I',PROPERTY='RANSWallR11Experiment',VALUE='On')
    tree.write(out/'case.xml',encoding='unicode')
    assert h.predictor.protected_xml(out/'case.xml')==manifest['protected_physical_xml_sha256']
    anchors=r10['accepted_model_anchors']
    assert len(anchors)==5
    replay=[dict(name='R10_anchor'+str(a['index']),direction=a['direction'],fraction=.125,anchor=True,
                 Phi_full=a['Phi']) for a in anchors]
    for name,fraction in (('linear',.125),('source',.125),('source',.0625),('source',.04541015625)):
        trial,=[t for t in r10['native_trials'] if t['comparison']=='r10_'+name and t['fraction']==fraction]
        replay.append(dict(name='R10_'+name+'_'+str(fraction),direction=anchors[0 if name=='linear' else -1]['direction'],
                           fraction=fraction,anchor=False,state=trial['state'],actual=trial['residual']))
    for row in replay:
        h.bound(row['direction'])
        if 'actual' in row:
            h.bound(row['state'])
            h.bound(row['actual'])
    manifest.update(current_target_xml=h.pin(out/'case.xml'),r11_replay=replay)
    (out/'predictor.json').write_text(json.dumps(manifest,indent=2)+'\n')
    budget=dict(affine_solves_per_policy=1,tangents_per_nonlinear_policy=3,
        armijo_trials_per_tangent=10,source_values_per_nonlinear_policy=50,native_trials_per_policy=9,
        cooperative_comparison_seconds=600,native_supervisor_seconds=650,f_star=.125,
        slope_multiplier=-2e-4,
        comparators={'Lf':['L1'],'N1':['L1'],'Nf':['L1','Lf if qualified']},
        seed_cost='N1 charged L1 construction; Nf charged Lf construction and any required L1 fallback',
        duplicate_cost='Exact native states reused; each policy charged the original evaluation cost',
        selection='Highest gain/second promoted policy. Otherwise retain lowest-E qualified trial.',
        replay_separate=True,physical_value_and_source_compatibility=h.load(out/'budget.json')['physical_value_and_source_compatibility'])
    (out/'budget.json').write_text(json.dumps(budget,indent=2)+'\n')


def summary(e, out, h, evidence):
    receipt=h.load(e/'local-run/receipt.json')
    collection=h.load(h.bound(receipt['durable_collection']['file']))
    events=module(evidence/'m8-r7/experiment.py','r7_events').events(collection)
    base=collection['completed_iterates'][0]
    assert base['state']['sha256']=='eb063e15f16b2610f1c5a832280d87894f567e8dab321560202e6dd9af95194c'
    assert base['history']['sha256']=='b8ba7917424057a9f6f16a026d82d53ede43985871537facea170ef6e1841b71'
    matrix,=[j for j in collection['jacobians'] if j['status']=='complete']
    assert matrix['event']['matrix']['sha256']=='2194965dee9c5960ba4e1c1604f76e8e85372cabef2d5ced7c8a394f593c8dc6'
    assert matrix['event']['details']['loaded_bound_matrix']
    assert collection['terminal']['terminal']=='r11_fixed_comparison_completed'
    trials=[r for r in events if r['kind']=='r8_trial']
    policies=[r for r in events if r['kind']=='r11_policy']
    chosen={}
    for p in policies:
        chosen[p['name']]=selected_trial(p,trials)
        assert sum(t['comparison']==p['name'] for t in trials)<=9
    for arm in ('N1','Nf'):
        anchors=[r for r in events if r['kind']=='r11_model_anchor' and r['name']==arm]
        for a,b in zip(anchors,anchors[1:]):
            assert b['Psi']<a['Psi'] and b['native_slope']<=-2e-4*max(base['scaled_block_norms'])
            assert b['fractional_C_identity_error']<=5e-12
        work=[r for r in events if r['kind']=='r11_model_work' and r['name']==arm]
        if work:assert work[0]['source_values']<=50 and work[0]['tangent_solves']<=3
    decision,=[r for r in events if r['kind']=='r11_decision']
    last=collection['last_snapshot']
    if decision['retained_policy'] is not None:
        retained=chosen[decision['retained_policy']]
        assert retained is not None and last['state']['sha256']==retained['state']['sha256']
        assert last['residual']['sha256']==retained['residual']['sha256']
    assert all(not s['history_advanced'] for s in collection['completed_iterates'])
    report=h.load(h.bound(receipt['native_report']))
    def digest(values):
        return hashlib.sha256(struct.pack('<'+str(len(values))+'d',*values)).hexdigest()
    h0=base['history']['sha256']
    start=report['stage_test_initial_native_history']
    assert report['time']==report['steps']==0 and digest(report['state'])==h0
    assert [digest(v) for v in start['history']]==[digest(v) for v in report['native_history']]
    assert start['times']==report['native_history_times'] and digest(start['state'])==h0
    level,=[v for v in report['stage_test_levels'] if v.get('attempted') and
            digest(v['state'])==last['state']['sha256'] and digest(v['history'])==h0]
    assert digest(level['residual'])==last['residual']['sha256'] and not level['native_history_advanced']
    h.save(out,dict(receipt=h.pin(e/'local-run/receipt.json'),collection=receipt['durable_collection']['file'],
        native_history_verified=True,native_report=receipt['native_report'],
        initial=base,selected=last,matrix=matrix['event']['matrix'],events=events,
        policies=policies,chosen=chosen,decision=decision))
    # Reuse a complete independent audit only after exact state/residual/history binding.
    r10=h.load(evidence/'m8-r10/decision-ledger.json')
    old_snapshot=h.load(h.bound(r10['final']['snapshot']))
    if all(last[k]['sha256']==old_snapshot[k]['sha256'] for k in ('state','residual','history')):
        h.bound(r10['final']['audit'])
        reuse=dict(audit=r10['final']['audit'],state=last['state'],
             residual=last['residual'],history=last['history'],reason='Exact audited payload bytes reused')
        reuse_path=e/'retained-audit-reuse.json'
        if reuse_path.exists():
            assert h.load(reuse_path)==reuse
        else:
            h.save(reuse_path,reuse)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('prepare','summary'))
    p.add_argument('--epoch',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if platform.system()!='Linux':raise ValueError('Use the designated Linux server')
    e=a.epoch.resolve()
    evidence=e/'source/docs/plans/cls-epawf-revised-evidence'
    h=module(evidence/'m8-r4/prepare.py','r11_prepare')
    (prepare if a.mode=='prepare' else summary)(e,a.out.resolve(),h,evidence)


if __name__=='__main__':
    main()
