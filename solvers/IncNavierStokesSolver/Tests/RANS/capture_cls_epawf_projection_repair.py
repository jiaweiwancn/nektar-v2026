#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Archive live campaign publications without modifying their processes/files."""
import argparse
import datetime
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

from run_cls_epawf_picard import load,require,sha,write_json,THREADS
from continue_cls_epawf_pimple import header,process_identity,load_progress


def snapshot(campaign,output):
    """A rotating pathname is evidence only after a hash-validated copy."""
    for attempt in range(20):
        progress=load_progress(campaign/'trajectory/progress.json')
        temporary=output.with_name(output.name+'.partial')
        if temporary.exists():shutil.rmtree(temporary)
        try:
            shutil.copytree(Path(progress['checkpoint']),temporary)
            receipt=header(temporary)
            if receipt['sha256']!=progress['checkpoint_sha256'] or receipt['steps']!=progress['accepted_steps'] or receipt['time']!=progress['accepted_time']:
                time.sleep(.05);continue
            os.replace(temporary,output)
            identities={}
            status=load_progress(campaign/'status.json')
            for name,expected in (('driver',dict(pid=progress['pid'],start=progress['process_start'])),
                                  ('native',progress.get('native_identity')),('supervisor',status.get('supervisor_identity'))):
                if expected:
                    actual=process_identity(expected['pid'])
                    identities[name]=dict(recorded=expected,actual=actual,
                        matches=actual is not None and all(str(actual.get(k))==str(v) for k,v in expected.items()))
            return dict(progress=progress,status=status,processes=identities,
                archived_checkpoint=dict(path=str(output),**receipt))
        except FileNotFoundError:
            if attempt==19:raise
            time.sleep(.05)
    raise RuntimeError('no stable native publication; no checkpoint fabricated')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--uniform',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False);r=a.root.resolve();u=a.uniform.resolve()
    out=dict(schema=1,captured_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),root=str(r),
        source_revision=subprocess.check_output(['git','-C',str(r/'source'),'rev-parse','HEAD'],text=True).strip(),
        source_status=subprocess.check_output(['git','-C',str(r/'source'),'status','--porcelain'],text=True).strip(),
        threads=THREADS,ranks=1,repair_verification=load(r/'repair-verification.json'),
        repaired=snapshot(r/'repaired-campaign',a.output/'repaired-checkpoint'),
        uniform=snapshot(u/'uniform-campaign',a.output/'uniform-checkpoint'),uniform_milestones=[],
        numerical_source_sha256={},orchestration_source_sha256={})
    core=r/'source/solvers/IncNavierStokesSolver/EquationSystems'
    for name in ('RANSWallKBoundProjection.h','RANSWallEvolutionKBound.cpp','RANSWallPhysicsSignature.h',
                 'RANSWallEvolutionPimple.cpp','RANSWallEvolutionPatch.h','RANSWallEvolutionPatch.cpp',
                 'RANSEquilibriumWallEvolution.cpp','RANSWallEvolutionFeasible.cpp','RANSWallEvolutionContact.cpp',
                 'RANSWallTimeIntegration.cpp','RANSHistory.cpp'):
        out['numerical_source_sha256'][name]=sha(core/name)
    for name in ('run_cls_epawf_pimple.py','continue_cls_epawf_pimple.py','supervise_cls_epawf_pimple.py',
                 'monitor_cls_epawf_constitutive_v2.py','replay_cls_epawf_projection_limit.py',
                 'analyse_cls_epawf_projection_repair.py','capture_cls_epawf_projection_repair.py'):
        out['orchestration_source_sha256'][name]=sha(Path(__file__).with_name(name))
    out['runtimes']={}
    for name,folder in (('repaired',r/'repaired-campaign/trajectory/runtime'),('uniform',u/'uniform-campaign/trajectory/runtime'),
                         ('probe',r/'probe-budget-runtime')):
        pinned=load(folder/'identity.json')['frozen']
        require(sha(Path(pinned['path']))==pinned['sha256'] and
                all(sha(Path(path))==digest for path,digest in pinned['libraries'].items()),'frozen runtime changed: '+name)
        out['runtimes'][name]=pinned
    original=r.parent/'cls-epawf-pimple-constitutive-v2-20261106/campaign/trajectory/progress.json'
    old=load(original);out['original_terminal']={k:old.get(k) for k in
        ('status','accepted_steps','accepted_time','failure','checkpoint','checkpoint_sha256')}
    out['original_terminal'].update(path=str(original),sha256=sha(original))
    completed=[int(f.parent.name.split('-')[-1]) for f in (u/'uniform-monitor').glob('step-*/summary.json')]
    for step in sorted({160,800,max(completed)}):
        folder=u/f'uniform-monitor/step-{step:06d}';d=load(folder/'summary.json');s=d['states']['step'+str(step)]
        h0=d['states']['initial']
        out['uniform_milestones'].append(dict(step=step,time=s['time'],source=str(folder/'summary.json'),sha256=sha(folder/'summary.json'),
            actual_initial_sha256=d['actual_initial_sha256'],governing_solves=d['governing_solves'],
            initialized_reference={k:h0[k] for k in ('global_k_integral','near_k_integral','mean_wall_stress','spatial_block_norms')},
            qualification_events=[e for segment in d['completed_native_segments'] for e in segment['qualification_events']][:1],
            state={k:s[k] for k in ('spatial_block_norms','epsilon_trace_defect','global_k_integral','near_k_integral',
                'mean_wall_stress','rms_wall_stress','field_L2_norms_u_v_k_epsilon_pinner_pouter','flow_through_time','H_over_ustar',
                'profiles','wall_B','wall_lambda')},window=d['windows'][-1],correction=d['cumulative_correction'],
            elapsed_wall_times=s['time']/s['H_over_ustar'],elapsed_initial_flow_times=s['time']/s['flow_through_time']))
    records=[json.loads(x) for x in (r/'repaired-campaign/trajectory/steps.jsonl').read_text().splitlines()]
    through=out['repaired']['progress']['accepted_steps'];records=[v for v in records if v['global_step']<=through]
    out['repaired_steps']=[]
    for row in records:
        v=row['attempt'];b=v['k_bound'];require(v['scheduled_complete'] and v['governing_solves_complete'],'accepted incomplete schedule')
        require(all(v[k]==1 for k in ('momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves','bound_calls')),'changed schedule')
        require(v['additional_governing_equation_solves']==0,'extra production PDE solve')
        out['repaired_steps'].append(dict(step=row['global_step'],time=row['time'],checkpoint_sha256=row['checkpoint_sha256'],
            projection={k:b[k] for k in ('separation_passes','inventory_rebuilds','active_set_iterations','kkt_solves','mass_factor_applications',
                'seconds_total','seconds_algebra','seconds_inventory','seconds_separation','primal_floor_defect','corrected_minimum','corrected_ground',
                'mass_norm_correction','mass_norm_raw_step','correction_to_step_ratio','integral_total','integrals_inner_outer_near_positive_negative',
                'discovered_constraint_families','maximum_work_one_add','minimum_observed_qr_residual')},
            original_refreshed_stage_blocks=v['refreshed_original_block_norms'],continuity_gauge=v['corrected_continuity_gauge'],
            matching=v['refreshed_trace_defect'],B_min=v['k_domain']['minimum_b'],lambda_min=v['k_domain']['minimum_lambda'],
            added_A_delta=v['k_added_frozen_residual'],storage_alpha_M_delta=v['k_added_storage_residual']))
    for name in ('readback73','readback80'):
        f=r/name/'summary.json'
        if f.exists():
            d=load(f);out[name]=dict(path=str(f),sha256=sha(f),states=d['states'],windows=d['windows'],
                cumulative_correction=d['cumulative_correction'],window_budgets=d['window_integrated_budgets'])
    artifacts=[]
    def index(path,role):
        require(path.exists(),'missing indexed artifact: '+str(path))
        artifacts.append(dict(path=str(path),role=role,bytes=path.stat().st_size,sha256=sha(path)))
    for entry in out['repair_verification']['artifacts']:index(Path(entry['path']),entry['role'])
    for folder in ('replay-old-controls','replay-total20000','replay-outer-budget1','earlier-active6','rollback-active73'):
        index(r/folder/'manifest.json',folder+' manifest')
    for f in ('launched-jobs.json','repair-verification.json','qp-tests2.log','driver-tests2.log','native-corrector-checks-h0/summary.json',
              'solver-runtime/identity.json','probe-budget-runtime/identity.json','repaired-campaign/trajectory/runtime/identity.json'):
        index(r/f,'runtime/validation provenance')
    for who in ('repaired','uniform'):index(Path(out[who]['archived_checkpoint']['path'])/'State.bin',who+' immutable latest accepted native checkpoint')
    for folder in (a.output,r/'readback73',r/'readback80',r/'figures'):
        if folder.exists():
            for f in sorted(folder.glob('*.png')):index(f,'scientific figure')
    write_json(a.output/'artifacts.json',artifacts)
    out['artifact_index']=dict(path=str(a.output/'artifacts.json'),sha256=sha(a.output/'artifacts.json'))
    write_json(a.output/'native.json',out);print(a.output/'native.json')


if __name__=='__main__':main()
