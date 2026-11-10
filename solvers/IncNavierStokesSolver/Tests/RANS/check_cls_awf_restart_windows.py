#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Confirm steady windows across the T=5 restart, including all outer epsilon.

Use the existing strict checkpoint probe unchanged. Interior GL4 points on each
native flat cubic outer element integrate squared epsilon differences exactly.
No PDE, Lift, projection, history advance or production configuration change.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone
import json
import math
import os
from pathlib import Path
import shutil

from analyse_cls_epawf_height_comparison import differences
from analyse_cls_accepted_awf import normalized
from continue_cls_epawf_pimple import atomic, header, load_progress
from diagnose_cls_accepted_covariance import checkpoint, common_packet, run_probe, value
from run_cls_epawf_picard import THREADS, load, require, sha, state_hash


def lineage(campaign, first):
    sources=[];seen=set()
    while True:
        require(campaign not in seen, 'cyclic campaign prefix');seen.add(campaign)
        p=load_progress(campaign/'progress.json');sources.append((campaign,p))
        if p.get('prefix_steps',0)<=first:break
        require(p.get('prefix_campaign'), 'missing native prefix lineage')
        prefix=Path(p['prefix_campaign']);parent=load_progress(prefix/'progress.json')
        require(parent['accepted_steps']==p['prefix_steps'], 'prefix accepted frontier changed')
        report=load(Path(p['prefix_report']));saved=checkpoint(Path(parent['checkpoint']))
        require(sha(Path(p['prefix_report']))==p['prefix_report_sha256'] and
            sha(Path(parent['checkpoint'])/'State.bin')==parent['checkpoint_sha256'], 'prefix report/checkpoint hash differs')
        require(state_hash(saved['coupled'])==state_hash(report['state']) and
            [state_hash(v) for v in saved['history']]==[state_hash(v) for v in report['native_history']] and
            state_hash(saved['history_times'])==state_hash(report['native_history_times']), 'prefix native state/history differs')
        campaign=prefix
    records={}
    for path,p in reversed(sources):
        for line in (path/'steps.jsonl').read_text().splitlines():
            try:row=json.loads(line)
            except json.JSONDecodeError:break
            n=int(row['global_step'])
            if n in records:require(records[n]['checkpoint_sha256']==row['checkpoint_sha256'], 'conflicting native step identities')
            records[n]=row
    return sources,records


def endpoint(record):
    folder=Path(record['segment']);matches=[]
    for file in folder.glob('checkpoint*/State.bin'):
        if header(file.parent)['steps']==record['global_step'] and sha(file)==record['checkpoint_sha256']:matches.append(file.parent)
    require(matches, 'retained exact native checkpoint unavailable: '+str(record['global_step']))
    return sorted(matches)[0],folder/'case.xml'


def outer_rule(definition,order):
    import numpy as np
    require(definition['velocity_scalar_degree']==3 and definition['H']==.5, 'full epsilon rule is the flat cubic H=.5 pilot')
    x=definition['tangential_knots'];z=[v for v in definition['radial_knots'] if v>=definition['H']]
    nodes,weights=np.polynomial.legendre.leggauss(order);points=[];physical=[]
    for xa,xb in zip(x,x[1:]):
        for za,zb in zip(z,z[1:]):
            for a,wx in zip(nodes,weights):
                for b,wz in zip(nodes,weights):
                    points.append([float(xa+(a+1)*(xb-xa)/2),0.,float(za+(b+1)*(zb-za)/2)])
                    physical.append(float(wx*wz*(xb-xa)*(zb-za)/4))
    require(abs(math.fsum(physical)-(x[-1]-x[0])*(z[-1]-z[0]))<1e-12, 'physical outer measure differs')
    return points,physical


def epsilon_norm(weights,left,right):
    require(len(weights)==len(left)==len(right), 'epsilon quadrature shape differs')
    return math.sqrt(math.fsum(w*(a-b)**2 for w,a,b in zip(weights,left,right)))


def readback(probe,record,definition,folder,order=4):
    folder.mkdir(parents=True,exist_ok=False);cp,xml=endpoint(record);copy=folder/'checkpoint';copy.mkdir()
    shutil.copy2(cp/'State.bin',copy/'State.bin');shutil.copy2(xml,folder/'case.xml')
    require(sha(copy/'State.bin')==record['checkpoint_sha256'], 'immutable native copy differs')
    saved=checkpoint(copy);points,weights=outer_rule(definition,order)
    packet=dict(common_packet(copy),operation='selected_original_operator',state=saved['coupled'],history=saved['coupled'],
        alpha=0.,time=saved['time'],inner_stress_model=definition['model_choices']['inner_stress'],
        witnesses=points,reporting_breakpoints=definition['radial_knots'])
    result,receipt=run_probe(probe,folder/'case.xml',packet,folder/'readback')
    require(result['governing_solves']==result['physical_commits']==0 and not result['native_history_advanced'] and
        result['initialization_work']==[0,0], 'native observation advanced or lifted the state')
    selected=result['diagnostics']['selected_state'];physical=result['diagnostics']['physical'];audit=result['complete_audits']['selected_state']
    require(physical['state_sha256']==selected['state_sha256']==state_hash(saved['coupled']) and
        physical['state_unchanged'] and physical['common_reporting_queries_valid'] and
        audit['accepted_awf_hard_admissible'][0]==1, 'native selected-state observation is inadmissible')
    witnesses=selected['one_sided_witnesses'];require(len(witnesses)==len(points), 'epsilon interior element membership changed')
    values=[];error=0.;geometry=[]
    for p,row in zip(points,witnesses):
        require(not row['inner'] and row['gradient_region_matches_law'] and
            max(abs(a-b) for a,b in zip(p,row['coordinate']))<2e-13 and
            row['epsilon_semantics'].startswith('actual retained field'), 'epsilon point has incompatible provenance')
        reconstructed=value(row['epsilon_evaluation_row'],saved['coupled'],int(selected['offsets'][3]))
        error=max(error,abs(reconstructed-row['epsilon']));values.append(row['epsilon']);geometry.append(row['geometry'])
    require(error<2e-12 and all(math.isfinite(v) and v>0 for v in values), 'independent retained-epsilon reconstruction differs')
    return dict(time=saved['time'],step=saved['steps'],physical=physical,audit=audit,epsilon_values=values,
        epsilon_weights=weights,epsilon_geometry=geometry,receipt=receipt,checkpoint_sha256=sha(copy/'State.bin'),
        state_sha256=state_hash(saved['coupled']),history_sha256=[state_hash(v) for v in saved['history']],
        history_times_sha256=state_hash(saved['history_times']),epsilon_reconstruction_error=error)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--campaign',action='append',required=True)
    parser.add_argument('--probe',type=Path,required=True);parser.add_argument('--policy',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);a=parser.parse_args();os.environ.update(THREADS)
    a.output.mkdir(parents=True,exist_ok=False);policy=load(a.policy);steps=[7200,8000,8800,9600]
    require(policy['primary_physical_window']==.5 and policy['required_consecutive_windows']==3, 'scientific policy differs')
    out=dict(captured_Asia_Shanghai=datetime.now(timezone(timedelta(hours=8))).isoformat(),policy=policy,policy_sha256=sha(a.policy),
        probe_sha256=sha(a.probe),script_sha256=sha(Path(__file__)),steps=steps,cases={},status='running',
        governing_solves=0,physical_commits=0,native_histories_advanced=False,
        epsilon_metric='Full transported H<=z<=Z; interior GL4x4 on each native flat cubic element; fixed common E0*sqrt(12) normalization; prescribed and matching trace contributions retained.')
    atomic(a.output/'summary.json',out)
    def case(item):
        label,path=item.split('=',1);require(Path(label).name==label, 'invalid label')
        folder=a.output/label;folder.mkdir();sources,records=lineage(Path(path),steps[0]);definition=sources[0][1]['case_definition']
        require(sources[0][1]['accepted_steps']>=steps[-1] and all(n in records for n in range(steps[0],steps[-1]+1)), 'three final windows are incomplete')
        observations={n:readback(a.probe,records[n],definition,folder/f'step-{n:06d}') for n in steps}
        summaries={str(n):{k:v for k,v in q.items() if k not in ('physical','epsilon_values','epsilon_weights','epsilon_geometry')} for n,q in observations.items()}
        for n,q in observations.items():
            summaries[str(n)].update(common_normalized_spatial_blocks=q['physical']['common_normalized_spatial_first4'],
                dimensional_spatial_blocks=q['physical']['dimensional_mass_dual_spatial_first4'],mean_wall_stress=q['physical']['mean_wall_stress'],
                matching_defect=q['physical']['matching_trace_defect'])
        windows=[]
        for first,last in zip(steps,steps[1:]):
            left,right=observations[first],observations[last];duration=right['time']-left['time']
            require(abs(duration-.5)<1e-10 and left['epsilon_geometry']==right['epsilon_geometry'], 'native window geometry/time differs')
            change=normalized(differences(right['physical'],left['physical']))
            eps=epsilon_norm(left['epsilon_weights'],left['epsilon_values'],right['epsilon_values'])
            eps_normal=eps/(policy['common_normalization']['epsilon']*math.sqrt(policy['common_normalization']['outer_epsilon_volume']))
            bound=[records[n]['attempt']['k_bound'] for n in range(first+1,last+1)]
            positive=math.fsum(v['integrals_inner_outer_near_positive_negative'][3] for v in bound)
            negative=math.fsum(v['integrals_inner_outer_near_positive_negative'][4] for v in bound)
            scale=policy['common_normalization']['K']*24/5
            row=dict(start=left['time'],end=right['time'],normalized_field_change=change,full_outer_epsilon_L2_change=eps,
                full_outer_epsilon_normalized_change=eps_normal,relative_traction_change=abs(right['physical']['mean_wall_stress']-left['physical']['mean_wall_stress'])/.16,
                signed_correction_integral=math.fsum(v['integral_total'] for v in bound),absolute_correction_integral=positive-negative,
                positive_correction_integral=positive,negative_correction_integral=negative,bound_activations=sum(v['activations'] for v in bound),
                normalized_absolute_supply_rate=(positive-negative)/duration/scale)
            row.update(fields_pass=all(change[k]<=limit for k,limit in policy['normalized_window_L2_change_maximum'].items()) and eps_normal<=policy['normalized_window_L2_change_maximum']['epsilon'],
                residuals_pass=max(right['physical']['common_normalized_spatial_first4'])<=policy['normalized_selected_spatial_residual_maximum_each_block'],
                matching_pass=right['physical']['matching_trace_defect']/policy['common_normalization']['epsilon']<=policy['relative_matching_defect_maximum'])
            row['traction_pass']=row['relative_traction_change']<=policy['relative_wall_stress_window_change_maximum']
            row['bound_supply_pass']=row['normalized_absolute_supply_rate']<=policy['bound_supply_per_time_relative_to_K0_volume_per_wall_time_maximum']
            row['screen_pass']=all(row[k] for k in ('fields_pass','residuals_pass','matching_pass','traction_pass','bound_supply_pass'));windows.append(row)
        # Independently overintegrate the last full-outer epsilon difference.
        check=[readback(a.probe,records[n],definition,folder/f'GL5-step-{n:06d}',5) for n in steps[-2:]]
        fine=epsilon_norm(check[0]['epsilon_weights'],check[0]['epsilon_values'],check[1]['epsilon_values'])
        require(abs(fine-windows[-1]['full_outer_epsilon_L2_change'])<2e-12, 'GL4/GL5 physical epsilon norms disagree')
        passed=sum(w['screen_pass'] for w in windows)
        return label,dict(campaign=path,production_status=sources[0][1]['status'],frontier_steps=sources[0][1]['accepted_steps'],frontier_time=sources[0][1]['accepted_time'],
            lineage=[str(p) for p,_ in sources],observations=summaries,windows=windows,passing_windows=passed,
            full_outer_epsilon_GL4_GL5_absolute_difference=abs(fine-windows[-1]['full_outer_epsilon_L2_change']),
            classification='candidate_steady_selected_operator_at_declared_screen' if passed==3 else 'not_steady_by_declared_screen',
            caveat='Exploratory selected-model screen; not continuum/spatial accuracy, exact steady root, or neutral-reference preservation.')
    with ThreadPoolExecutor(max_workers=2) as pool:
        for label,result in pool.map(case,a.campaign):out['cases'][label]=result;atomic(a.output/'summary.json',out)
    out.update(status='passed',completed_Asia_Shanghai=datetime.now(timezone(timedelta(hours=8))).isoformat());atomic(a.output/'summary.json',out)
    print(a.output/'summary.json')


if __name__=='__main__':main()
