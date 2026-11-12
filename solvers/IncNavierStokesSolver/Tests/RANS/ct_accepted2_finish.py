#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Final read-only CT2 primitives, negative controls and admission decision.
A completed diagnosis does not imply a completed physical pilot.
"""
import argparse,json,os,copy
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from ct_accepted_feasibility import write_json
from ct_accepted_open import require
from ct_accepted2_diagnosis import inventory,invariant,WITNESS
from ct_accepted2_storage import artifact,canonical_location
from ct_accepted_pilot import budget_classification
from rans_request_manifest import RequestManifest,sha,runtime_identity

def physical_gate(case,T=.008):
    # Recompute from physical primitives, ignoring any enclosing pass label.
    reasons=[]
    if case['failure']:reasons.append('native_physical_rejection')
    if abs(case['time']-T)>1e-13:reasons.append('required_horizon_not_reached')
    if case.get('fixed_domain_failures',0):reasons.append('independent_fixed_grid_constitutive_failure')
    return dict(status='pass' if not reasons else 'fail',reasons=reasons)
def gate_tests():
    good=dict(failure='',time=.008,fixed_domain_failures=0,status='pass')
    require(physical_gate(good)['status']=='pass','good_gate_control')
    cases={'native_rejection':dict(good,failure='A95 negative square-root radicand',status='pass'),
           'short_horizon':dict(good,time=.004,status='pass'),
           'dense_domain':dict(good,fixed_domain_failures=1,status='pass')}
    out={}
    for name,c in cases.items():
        q=physical_gate(c);require(q['status']=='fail','pass_label_overrode_primitive:'+name);out[name]=dict(status='pass',expected_rejection=q)
    G=np.diag([-6.,0.,6.]);v=invariant(G,1.,1.)
    require(np.trace(G)==0 and v['radicand']<0,'tracefree_A95_negative_control')
    out['exactly_incompressible_A95']=dict(status='pass',expected='negative radicand despite exact zero divergence',witness=v,threshold_a_tau=6.5/(np.sqrt(18)-3))
    return out

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    r=a.root.resolve();rt=r/'runtime';os.environ['LD_LIBRARY_PATH']=str(rt/'lib')
    manifest=RequestManifest(a.output,dict(schema='CT2FinalPrimitivesV1',runtime=runtime_identity(rt/'IncNavierStokesSolver'),reports={str(f.relative_to(r)):sha(f) for f in r.glob('**/case.evolution.json')},source_sha256=sha(Path(__file__)),required_layouts=['B','X','Zplus','XZplus','Zminus'],dt=[.002,.001],required_T=.008,conditional_T=[.02,.2]))
    result=dict(schema='CT2ScientificDecisionV1',diagnostic_execution='pass',negative_controls=gate_tests(),cases={},conditional_admission={},budget_normalization_test={})
    for label,prep,tra in [('primary','initial-v2','matrix-v1'),('companion','companion-initial-v1','companion-v1')]:
        folders=sorted((r/tra/'work').glob('*'))
        for folder in folders:
            if not folder.is_dir() or not (folder/'case.evolution.json').exists():continue
            name=folder.name;report=json.loads((folder/'case.evolution.json').read_text());ini=json.loads((r/prep/'work'/name/'initial.json').read_text());o=ini['observation'];off=np.array(o['offsets'],int);s=np.array(report['state'])
            native=inventory(o['volume'],[q['coordinate'] for q in o['volume']],s,off);fixed=inventory(ini['test_basis'],o['grid'],s,off)
            row=dict(time=report['time'],failure=report['failure'],steps=report['steps'],state_sha256=exact_vector_hash(s),native_report_sha256=sha(folder/'case.evolution.json'),fixed_domain_failures=fixed['domain_failure_count'],actual_endpoint_native=native,actual_endpoint_fixed=fixed,physical_gate=physical_gate(dict(time=report['time'],failure=report['failure'],fixed_domain_failures=fixed['domain_failure_count'])),last_attempt_time=report['pimple_attempts'][-1]['time'])
            result['cases'][name]=row
    expected={x+'-'+str(dt) for x in ('B','X','Zplus','XZplus','Zminus') for dt in (.002,.001)}
    require(set(result['cases'])==expected,'missing_required_native_trajectory')
    for layout in ('B','X','Zplus','XZplus','Zminus'):
        members=[result['cases'][layout+'-'+str(dt)] for dt in (.002,.001)]
        result['conditional_admission'][layout]=dict(status='not_admitted' if any(x['physical_gate']['status']!='pass' for x in members) else 'requires_full_temporal_scalar_representation_gate',required_T=.008,extension_attempted=False)
    # Geometry normalization control is diagnostic only; it does not relabel
    # any physical case using a fictitiously changed volume.
    path=r/'primary-audits-v1/B-0.002/q64f32/budget.json'
    if path.exists():
        packet=json.loads(path.read_text());c=json.loads((r/'matrix-v1/work/B-0.002/construction.json').read_text())
        one=budget_classification(packet,c);changed=copy.deepcopy(c);changed['x']=[2*x for x in c['x']];two=budget_classification(packet,changed)
        require(two['physical_volume']==2*one['physical_volume'] and all(abs(two['equations'][k]['scaled_frozen_defect']-.5*one['equations'][k]['scaled_frozen_defect'])<1e-15 for k in one['equations']),'hardcoded_budget_length_remaining')
        result['budget_normalization_test']=dict(status='pass',actual_volume=one['physical_volume'],synthetic_doubled_volume=two['physical_volume'],scope='unit check only; actual pilot budgets retain actual geometry')
    for name,file in [('boundary_pressure_implementation','boundary-verification-v1/results.json'),('history_rollback_closed_compatibility','history-controls-v1/results.json'),('wall_reference_study','wall-heights-v2/results.json'),('shared_model_smokes','shared-smokes-v1/results.json'),('independent_projection_and_flow_replay','extras-v1/results.json')]:
        f=r/file
        if f.exists():
            data=json.loads(f.read_text());result[name]=dict(status=data.get('shared_path_smoke_status',data.get('status','unresolved')),receipt_sha256=sha(f))
        else:result[name]=dict(status='missing_required_receipt')
    controls=[]
    for f in r.glob('launch-*/worker-guard-highwater.jsonl'):
        for line in f.read_text().splitlines():
            q=json.loads(line);controls.append(dict(control=f.parent.name,**q))
    maximum=max(q['active_before_actions'] for q in controls);ownedmax=max(sum(p['threads'] for p in q['processes'] if p['owned']) for q in controls)
    # The launcher parent and its monitor are outside the leased child family.
    # Count both; otherwise a 16-thread Git/build family is falsely accepted.
    projectmax=ownedmax+2
    result['resource_compliance']=dict(status='pass' if maximum<=72 and projectmax<=16 else 'failed_project_allocation' if maximum<=72 else 'failed_aggregate_cap',aggregate_status='pass' if maximum<=72 else 'fail',project_status='pass' if projectmax<=16 else 'fail',maximum_observed_active=maximum,maximum_observed_owned_family=ownedmax,control_helpers_outside_family=2,maximum_observed_project=projectmax,aggregate_limit=72,project_limit=16,highwater_samples=controls,historical_CT1_79='unchanged failed historical gate')
    result['constitutive_admissibility']=dict(primary_matrix='failed and unchanged',selected_Zminus_short='pass' if all(result['cases']['Zminus-'+str(dt)]['physical_gate']['status']=='pass' for dt in (.002,.001)) else 'fail')
    result['selected_continuations']={};latest=None
    for name in ('continue02-v1','continue2-v1'):
        file=r/name/'results.json'
        if file.exists():
            d=json.loads(file.read_text());result['selected_continuations'][name]=d;latest=d
    result['physical_pilot_completion']='not_completed_T0p2'
    result['temporal_and_representation_uncertainty']='short/milestone receipts classified at their actual times; no order inferred from dt pairs'
    result['scalar_budget_accuracy']={}
    for name in ('primary-audits-v1','companion-audits-v1','continue02-audits-v1','continue2-audits-v1'):
        file=r/name/'results.json'
        if file.exists():
            d=json.loads(file.read_text());result['scalar_budget_accuracy'][name]=dict(status=d['scalar_accuracy_status'],mass_status=d['mass_status'],receipt_sha256=sha(file))
    final=r/'continue2-audits-v1/results.json'
    if latest and latest.get('physical_completion')=='pass' and latest.get('cases') and all(abs(c['time']-.2)<1e-13 and not c['failure'] for c in latest['cases'].values()) and final.exists():
        d=json.loads(final.read_text());comparisons=latest['comparisons']['Zminus'];end=next((x for x in comparisons if abs(x['time']-.2)<1e-13),None)
        from ct_accepted2_continue import physical_budgets
        primitive={}
        for name,row in latest['cases'].items():
            file=Path(latest['trajectory_directory'])/name/'case.evolution.json'
            require(sha(file)==row['native_report_sha256'],'changed_final_native_report')
            actual=json.loads(file.read_text());dt=float(name.split('-')[-1])
            require(not actual['failure'] and abs(actual['time']-.2)<1e-13 and actual['steps']==round(.2/dt),'wrong_final_native_horizon_or_steps')
            primitive[name]=physical_budgets(d['cases'][name])
        result['final_budget_primitive_reclassification']=primitive
        if d['mass_status']=='pass' and d['scalar_accuracy_status']=='pass' and end and end['status']=='pass' and all(q['status']=='pass' for q in primitive.values()):result['physical_pilot_completion']='pass_T0p2_bounded_open_floor_only'
    result['precursor_admission']=False;result['benchmark_admission']=False
    result['next_gate']='mechanism-based bounded decision; no automatic longer continuation and no input retuning'
    write_json(a.output/'results.json',result);manifest.seal(a.output,['results.json'])
    print('diagnosis completed;',result['physical_pilot_completion'])
    # Scientific gate failure is visible at the command line.
    return 0 if result['physical_pilot_completion']=='pass_T0p2_bounded_open_floor_only' else 1
if __name__=='__main__':raise SystemExit(main())
