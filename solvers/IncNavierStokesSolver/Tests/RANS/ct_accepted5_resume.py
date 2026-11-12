#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Explicit CT5 postprocessing recovery with unchanged completed trajectories."""
import argparse,ast,json,os,subprocess
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
from ct_accepted5_run import Campaign,basic_gate,traction_gate,physical,read
from ct_accepted_pilot import replay,artifact_hash
from ct_accepted_open import require
from ct_accepted4_reference import write_json
from rans_request_manifest import RequestManifest,sha,runtime_identity,canonical

def helper_equivalence(initial):
    g=['git','-c','index.threads=1','-c','core.preloadIndex=false','show','b6a5c3d:solvers/IncNavierStokesSolver/Tests/RANS/ct_accepted_pilot.py']
    source=subprocess.check_output(g,cwd=Path(__file__).resolve().parents[4],text=True);node=next(x for x in ast.parse(source).body if isinstance(x,ast.FunctionDef) and x.name=='replay')
    ns={};exec(compile(ast.Module(body=[node],type_ignores=[]),'old_replay','exec'),ns)
    for extra in ({'export_constitutive':True},{'export_physical_scalar_budget':True,'previous_state':[],'stage_history':[]}):
        old=ns['replay'](initial,initial['observation']['state'],.008,[],**extra);new=replay(initial,initial['observation']['state'],.008,[],**extra)
        require(canonical(old)==canonical(new),'changed_default_physical_replay_request')
    points=[[0.,0.,.15],[0.,0.,.15]];request=replay(initial,initial['observation']['state'],0.,points,grid_sides=[-1,1]);require(request['grid_sides']==[-1,1],'explicit_interface_sides_not_retained')
    try:replay(initial,[],0.,points,grid_sides=[-1])
    except (ValueError,RuntimeError) as err:require('replay_grid_side_length_mismatch' in str(err),'wrong_side_length_rejection')
    else:raise RuntimeError('invalid_side_inventory_accepted')
    return dict(status='pass',default_constitutive_and_budget_requests='byte-identical canonical dictionaries',new_explicit_sides=request['grid_sides'],bad_side_length='expected rejection')

def attach(camp,model):
    ip=camp.root/'initial'/model;prep=read(ip/'results.json');require(prep['status']=='pass','completed_initialization_failed')
    seal=read(ip/'completed-artifacts.json');require(all(sha(ip/f)==h for f,h in seal['files'].items()),'changed_initialization')
    camp.model=model;camp.c=prep['case'];camp.initial_hash=prep['members']['0.001']['state_sha256']
    for dt,row in prep['members'].items():
        require(row['state_sha256']==camp.initial_hash and row['identity']['runtime']==runtime_identity(camp.rt/'IncNavierStokesSolver'),'foreign_initial_state_or_runtime')
        camp.inits[dt]=(ip,prep)
    q=read(camp.root/'pairs'/model/'T0.008-0.001/results.json');seal=read(camp.root/'pairs'/model/'T0.008-0.001/completed-artifacts.json')
    require(all(sha(camp.root/'pairs'/model/'T0.008-0.001'/f)==h for f,h in seal['files'].items()),'changed_completed_pair')
    require(q['basic']==basic_gate(q['members'],q['temporal'],.008,(.002,.001),camp.initial_hash) and q['traction_uncertainty']==traction_gate(q['members'],camp.c,.008),'forged_completed_milestone')
    for dt,row in q['members'].items():
        f=Path(row['directory']);n=read(f/'case.evolution.json');require(sha(f/'case.evolution.json')==row['native_report_sha256'] and artifact_hash(n['checkpoint'])==row['history_sha256'],'changed_completed_trajectory_or_history')
        camp.members[(dt,.008)]=row
    old=read(camp.old/'initial-v3/A95-6.0-P4/results.json')['case'];require(physical(old)==physical(camp.c),'continuous_input_changed')
    return prep,q

def maybe_half(camp,row,q,target,a):
    """One conditional primary fresh history, including valid early gates."""
    spatial_uncertain=camp.model=='A95' and row.get('spatial',{}).get('status')!='resolved_improvement'
    allowed=all(r.endswith(':traction') or r=='temporal_fields_or_coverage' for r in q['basic']['reasons'])
    if camp.model!='A95' or camp.quarter or not allowed or q['basic']['physical_completion']!='pass' or not (q['traction_uncertainty']['status']=='ambiguous' or spatial_uncertain):return q
    camp.quarter=True;row['quarter_initial']=camp.prepare(camp.model,(.0005,),True)
    for t in (.008,.02) if target==.02 else (.008,):
        fine=camp.pair(t,(.001,.0005));row['milestones'][str(t)+'-quarter']=fine
        if t==.008:
            original=camp.root;camp.root=a.root/'recovery-analysis-quarter';row['spatial']=camp.spatial(fine);camp.root=original
        if fine['qualification']['status']!='pass' or row['spatial']['status']!='resolved_improvement':return fine if t==target else q
    return fine

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','ct4','runtime','diagnosis','contract','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');qfile=a.root/'pairs/A95/T0.008-0.001/results.json'
    m=RequestManifest(a.output,dict(schema='CT5ExactTrajectoryRecoveryV1',completed_pair_sha256=sha(qfile),original_controller_request_sha256=sha(a.root/'controller-v1/request-manifest.json'),
        contract_sha256=sha(a.contract),runtime=runtime_identity(a.runtime/'IncNavierStokesSolver'),sources={f:sha(Path(__file__).with_name(f)) for f in ('ct_accepted5_resume.py','ct_accepted5_run.py','ct_accepted_pilot.py')}))
    result=dict(status='incomplete',models={},reused_completed_native_members=2,new_governing_runtime=False)
    for model in ('A95','A93','E94E95'):
        if model!='A95' and result['models']['A95']['status']!='pass':result['models'][model]=dict(status='not_admitted_primary_failed');continue
        camp=Campaign(a);row=dict(milestones={});result['models'][model]=row
        if model=='A95':
            prep,q=attach(camp,model);row['initial']=dict(initialization=prep,physical_input_before=physical(camp.c),physical_input_after=physical(camp.c),source='validated existing CT5 initialization')
            result['helper_equivalence']=helper_equivalence(read(a.root/'initial/A95/0.001/initial.json'))
            row['milestones']['0.008']=q
            # The failed analysis destination is immutable. New diagnostic
            # destination, same common physical inventory and sealed histories.
            original=camp.root;camp.root=a.root/'recovery-analysis-v1';row['spatial']=camp.spatial(q);camp.root=original
        else:
            row['initial']=camp.prepare(model);q=camp.pair(.008);row['milestones']['0.008']=q
        q=maybe_half(camp,row,q,.008,a);row['milestones']['0.008']=q
        if q['qualification']['status']=='pass' and (model!='A95' or row['spatial']['status']=='resolved_improvement'):
            q=camp.pair(.02,tuple(q['dts']));row['milestones']['0.02']=q
            q=maybe_half(camp,row,q,.02,a);row['milestones']['0.02']=q
        row['status']='pass' if all(t in row['milestones'] and row['milestones'][t]['qualification']['status']=='pass' for t in ('0.008','0.02')) and (model!='A95' or row['spatial']['status']=='resolved_improvement') else 'failed_or_unresolved'
        row['half_step_control']='not_needed' if all(q['traction_uncertainty']['status']!='ambiguous' for q in row['milestones'].values()) else 'required_not_yet_executed'
        write_json(a.output/'progress.json',result)
    result['status']='pass' if all(q['status']=='pass' for q in result['models'].values()) else 'failed_or_unresolved';write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);print('CT5 recovered campaign', {k:v['status'] for k,v in result['models'].items()});return 0 if result['status']=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
