#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check the matched two-update 96-response Krylov experiment, not M8 convergence.

This intentionally specific evidence check compares the 16- and 64-Krylov
experiments. Both completed searches must have qualified nonroot receipts;
passing does not establish a preferred basis size or general performance.
"""
from pathlib import Path
import argparse
import json, hashlib, datetime, sys, xml.etree.ElementTree as ET

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('run_root', type=Path)
parser.add_argument('output', type=Path)
args = parser.parse_args()
R = args.run_root
names = ['response-krylov16-space96', 'response-krylov64-space96']
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def rows(name):
    return [json.loads(l[15:]) for l in (R/name/'target/solver.log').open()
            if l.startswith('EP-AWF Reduced ')]
def session(name):
    root = ET.parse(R/name/'target/case.xml').getroot()
    selected = []
    for p in root.findall('.//PARAMETERS/P'):
        if (p.text or '').strip().split('=')[0].strip() == 'RANSWallReducedKrylov':
            selected.append(p.text.strip())
            p.text = 'RANSWallReducedKrylov = MATCHED'
    assert len(selected) == 1, selected
    return ET.tostring(root), selected[0]
def omit_time(value):
    if isinstance(value, dict): return {k: omit_time(v) for k,v in value.items() if k != 'seconds'}
    if isinstance(value, list): return [omit_time(v) for v in value]
    return value
def selected_action_checks(records):
    names = {'conditional_newton','response_krylov_1','response_krylov_16',
             'response_krylov_17','response_krylov_64','ground_lift',
             'opposite_sample_-2.600000','opposite_sample_1.600000',
             'normalization_-2.600000','normalization_1.600000',
             'epsilon_lower_sample','epsilon_lower_trace',
             'contact_value_0','contact_tangent_0'}
    selected = []
    for r in records:
        if r['event'] != 'response_check' or r['name'] not in names:
            continue
        r = dict(r)
        for key in ['shift_domain','half_domain']:
            r[key] = {k:v for k,v in r[key].items()
                      if k not in ['profile_branch_pattern','inventory_counts']}
        selected.append(r)
    return selected

reviews = [json.loads((R/(n+'-completed-review.json')).read_text()) for n in names]
logs = [rows(n) for n in names]
sessions = [session(n) for n in names]
assert sessions[0][0] == sessions[1][0], 'physical/solver inputs differ beyond the declared Krylov parameter'
assert reviews[0]['manifest']['IncNavierStokesSolver']['sha256'] == reviews[1]['manifest']['IncNavierStokesSolver']['sha256']
assert reviews[0]['manifest']['IncNavierStokesSolver']['libraries'] == reviews[1]['manifest']['IncNavierStokesSolver']['libraries']
assert reviews[0]['manifest']['expected_h0_sha256'] == reviews[1]['manifest']['expected_h0_sha256']
initial = [[r for r in rs if r['event']=='initial'][0] for rs in logs]
assert initial[0]['state_sha256'] == initial[1]['state_sha256']
assert initial[0]['norms'] == initial[1]['norms']
checks = [[r for r in rs if r['event']=='response_check' and r['iteration']==0] for rs in logs]
assert len(checks[0]) >= 17 and len(checks[1]) >= 65
assert all(omit_time(a)==omit_time(b) for a,b in zip(checks[0][:17], checks[1][:17])), 'common response prefix is not bitwise identical'
out = dict(scope='Matched original-target saved-guess searches: different response-space composition only; qualified nonroots, not physical trajectory success',
    utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    command=sys.argv, script_sha256=sha(Path(__file__)),
    producing_code_commit='90dd80152ef5a8b69e92e05b7048c60f782cc2b3',
    only_xml_difference=[x[1] for x in sessions], same_executable=True,
    same_original_h0=True, same_nonlinear_start=True,
    common_first_model_response_prefix=17, common_prefix_bitwise_except_elapsed_seconds=True,
    original_initial_norms=initial[0]['norms'], cases={})
for name,x,rs in zip(names,reviews,logs):
    assert x['operator_checks']=='passed_operator_checks'
    assert x['native']['steps']==0 and x['native']['time']==0
    assert x['native']['accepted_state_and_history_preserved'] and x['native']['readback_norm_difference']==0
    a=x['independent_audit']; assert a['unknowns']==2165 and a['alpha']==400 and a['time']==.0025
    assert a['qualification']['r2_qualified'][0]==1
    assert max(a['original_norms']) > 2e-12, 'this comparison is scoped to completed nonroots'
    directions=[r for r in rs if r['event']=='direction']
    updates=[r for r in rs if r['event']=='accepted']
    assert len(directions)==len(updates)==2
    iterations=[]
    for d,u in zip(directions,updates):
        assert d['response_rank']==d['basis_dimension']==96
        trials=[r for r in rs if r['event']=='trial' and r['iteration']==d['iteration'] and r['accepted']]
        assert len(trials)==1
        t=trials[0]; assert t['displacement']<=d['radius']*(1+1e-8)
        assert u['qualified'] and u['ground']['admissible']
        iterations.append(dict(numerical_iteration=d['iteration'],rank=d['response_rank'],
            krylov_directions=sum(n.startswith('response_krylov_') for n in d['basis_names']),
            active_constraints=d['active_constraints'],
            active=[{k:v for k,v in cut.items() if k!='row'} for cut in d['active']],
            predicted_norm=t['predicted_norm'],actual_norm=t['norm'],
            predicted_violation=t['predicted_violation'],actual_violation=t['violation'],
            model_error=t['model_error'],H_error=t['non_k_error'],
            actual_displacement=t['displacement'],radius=d['radius'],
            ground=u['ground'],outer_range=u['outer_range'],
            state_sha256=u['state_sha256'],work=u['work'],seconds=u['seconds']))
    out['cases'][name]=dict(review=str(R/(name+'-completed-review.json')),
        review_sha256=sha(R/(name+'-completed-review.json')),
        native=x['native'],non_k=x['non_k'],response_summary=x['response_summary'],
        selected_two_size_signed_action_checks=selected_action_checks(rs),
        source_predictions=[r for r in rs if r['event']=='source_ray_prediction'],
        rejected_trial_witnesses=[dict(iteration=r['iteration'],fraction=r['fraction'],
            norm=r['norm'],predicted_norm=r['predicted_norm'],
            violation=r['violation'],predicted_violation=r['predicted_violation'],
            H_error=r['non_k_error'],
            witnesses=r['qualification_inventory']['witnesses'])
            for r in rs if r['event']=='trial' and not r['accepted']],
        iterations=iterations,original_final_norms=a['original_norms'],
        independent_audit_receipt=a['receipt'],fresh_qualification=a['qualification'],
        elapsed_reduced_seconds=rs[-1]['seconds'])
out['conclusion']='64 Krylov directions improve the first update, but 16 reach a lower residual after two updates; no uniform ordering or stationarity follows. Both searches continue separately.'
p=args.output; p.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n'); print(p)
