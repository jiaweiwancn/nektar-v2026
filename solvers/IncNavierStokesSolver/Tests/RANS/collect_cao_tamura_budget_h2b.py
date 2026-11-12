#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Reduce executed H2b receipts; retain heavy native evidence by path and hash."""
import argparse
from collections import defaultdict
from datetime import datetime, timezone, timedelta
import hashlib
import json
import math
from pathlib import Path
import statistics
import numpy as np
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import load, require, sha, state_hash
from diagnose_cls_accepted_covariance import value

LARGE = {'physical_mass_action','native_pressure_action','native_spatial_residual',
         'scalar_frozen_load','epsilon_frozen_load','scalar_free_rows','epsilon_free_rows','projection_only'}
QP = ('success','failure','inventory_policy','floor','arithmetic_reserve','raw_minimum','corrected_minimum',
      'calls','activations','separation_passes','inventory_rebuilds','active_set_iterations','kkt_solves',
      'mass_factorizations','mass_factor_applications','primal_floor_defect','ground_passed','ground_floor_defect',
      'current_inventory_passed','completed_primal_dual_coherent','finite_pool_complete','stationarity_mass_dual',
      'complementarity','dual_feasibility_defect','primal_objective','primal_dual_gap_working_set','seconds_total',
      'seconds_algebra','seconds_inventory','seconds_separation','raw_sha256','corrected_sha256')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--h1',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--source-revision',default='working-tree')
    a=p.parse_args();r=a.root.resolve();a.output.mkdir(parents=True,exist_ok=False)
    def read(path):return load(r/path)
    def save(name,q):atomic(a.output/(name+'.json'),q)
    def result(root):return read(root+'/composite-hill-probe/result.json')
    on=result('curved-lazy-final');before=result('curved-final-on');off=result('curved-final-off');datum=result('curved-nonzero-datum')
    for q in (before,off):
        require(on['state_sha256']==q['state_sha256'] and on['audit']==q['audit'] and on['scalar_inventory']==q['scalar_inventory'],
                'exact optimization changed state/audit/inventory')
        for name in LARGE-{'projection_only'}:require(on['checks'][name]==q['checks'][name],'optimization changed native action '+name)
        require(on['checks']['projection_only']['corrected_sha256']==q['checks']['projection_only']['corrected_sha256'],
                'optimization changed active correction')
    require(on['checks']['status']=='passed' and on['audit']['accepted_awf_hard_admissible']==[0],
            'expected operator-pass/seed-rejection distinction lost')
    require(datum['checks']['Pout']==.013 and abs(datum['checks']['outlet_datum_traction_integral']-.0104)<1e-13,
            'nonzero mapped outlet datum failed')
    native=read('mapped-h0-v1/hill-h0/case.evolution.json')
    require(native['steps']==native['time']==native['implicit_callback_count']==0 and
            state_hash(native['state'])==on['state_sha256'],'native rejected seed differs from diagnostic fixture')
    spectra=[]
    for name in ('minimum_inner','minimum_outer'):
        w=on['seed_constitutive'][name];g=np.array([value(w['velocity_evaluation_row'],on['state'],int(on['seed_constitutive']['offsets'][b]),True) for b in (0,1)])
        require(np.max(abs(g-np.array(w['velocity_gradient'])[[0,2],:]))<1e-10,'independent coefficient-derived gradient differs')
        tensor=np.array(w['selected']['tensor']);eig=np.linalg.eigvalsh(tensor)
        require(abs(eig[0]-w['selected']['eigenvalues'][0])<1e-13,'independent tensor eigenvalue differs')
        spectra.append(dict(name=name,eigenvalues=eig.tolist(),gradient_reconstruction_max=float(np.max(abs(g-np.array(w['velocity_gradient'])[[0,2],:])))))
    # The pre-2D lift uses the NATIVE 1D radial projection. Separately check
    # the exact incoming AWF radial gradient at the decisive inner foot.
    preflight=load(a.h1/'preflight.json');c=preflight['candidates'][1]
    w=on['seed_constitutive']['minimum_inner'];K=w['k'];nu=native['parameters'][0]
    nt=w['selected']['eddy_viscosity'];H=.015;curvature=.04*math.pi**2/(2*.1**2)
    upstream_traction=native['ustar']**2
    exact_radial_gradient=upstream_traction/(nu+nt)
    turning=c['finite_UH']*curvature/(1-H*curvature)
    exact_inner=dict(scope='same tangential-following streamfunction seed with exact flat incoming AWF radial derivative; H1 calibrated UH=u_star^2 J0; not a general hill solution',
        radial_gradient=exact_radial_gradient,turning_dW_dx=turning,
        minimum_linear_covariance=2*K/3-nt*(exact_radial_gradient+turning))
    require(exact_inner['minimum_linear_covariance']<0,'analytical inner seed interpretation changed')
    qp=on['checks']['projection_only'];objective=on['checks']['projection_physical_mass_objective']
    require(abs(objective-qp['primal_objective'])<=1e-10*max(1e-30,objective),'physical mass objective differs from raw-centered QP')
    geometry=dict(checks={k:v for k,v in on['checks'].items() if k not in LARGE},
        projection={k:on['checks']['projection_only'].get(k) for k in QP},
        scalar_inventory=on['scalar_inventory'],setup=on['seed_setup'],
        nonzero_outlet=dict(Pout=.013,integral=datum['checks']['outlet_datum_traction_integral']),
        flat_limit=read('operators-v5/progress.json')['tests']['composite_flat_limit'],
        native_source_budget=on['native_budget'])
    save('geometry',geometry)
    seed=dict(status='initialization_rejected',hill_physical_steps=0,attempted_time=0.,native_failure=native['failure'],
        state_sha256=on['state_sha256'],native_report_sha256=sha(r/'mapped-h0-v1/hill-h0/case.evolution.json'),
        reference_sha256=state_hash(native['fixed_reference']),defect_sha256=state_hash(native['fixed_defect']),
        history_sha256=[state_hash(s) for s in native['native_history']],history_times=native['native_history_times'],
        complete_audit=on['audit'],one_sided_pre2D_lift_witnesses=on['seed_target_witnesses'],
        minimum_inner=on['seed_constitutive']['minimum_inner'],minimum_outer=on['seed_constitutive']['minimum_outer'],
        independent_spectra=spectra,exact_radial_inner_control=exact_inner,
        classification='incompatible chosen velocity seed; 1D/2D representation worsens it; no proof of universal hill-model failure')
    for item in ('minimum_inner','minimum_outer'):
        for key in ('velocity_evaluation_row','epsilon_evaluation_row','endpoint_evaluation_row'):seed[item].pop(key,None)
    save('seed',seed)
    temporal=read('flat-local-v2/progress.json')
    seed_observation=read('flat-local-v2/coarse-seed-readback/result.json')
    face=next(v for v in seed_observation['developing_boundary']['faces'] if v[0]==0 and v[1]==-1)
    bulk=-face[3]/face[2]
    flow_scales=dict(inlet_bulk_velocity=bulk,domain_flowthrough=.4/bulk,inlet_to_observation=.2/bulk,
        Hw_over_ustar=.015/native['ustar'],elapsed_over_domain_flowthrough=.0004/(.4/bulk),
        physical_tstar_start=.0008/.04,physical_tstar_end=.0012/.04)
    save('temporal',dict(status=temporal['status'],parent=temporal['parent'],parent_sha256=temporal['parent_sha256'],
        physical_t0=temporal['physical_t0'],steps=temporal['governing_steps_in_temporal_triplet'],
        policy=temporal['policy'],normalizations=temporal['normalizations'],comparisons=temporal['temporal'],
        qualification=temporal['qualification'],runtime=temporal['runtime'],probe=temporal['probe'],flow_scales=flow_scales))
    representation=read('representation-v1/progress.json')['comparison'];rep={k:representation[k] for k in ('scales','before_projection','after_projection','steps','selection')}
    for key,name in [('before','retained-mesh-h0'),('after','split-mesh-h0')]:
        d=read('representation-v1/'+name+'-readback/result.json');s=d['physical_samples'];b=d['native_turbulence_budget']
        rep[key]=dict(traction=math.fsum(x*y for x,y in zip(s['wall_weights'],s['wall_stress']))/math.fsum(s['wall_weights']),
            reaction=b['ground_reaction'],inner_rates=b['inner_rates'],outer_rates=b['outer_rates'],
            independent_source_integration=b['independent_source_integration'],weak_components=b['weak_components'],
            complete_audit=d['complete_audit'],slots=len(read('representation-v1/'+name+'/case.evolution.json')['state']))
    save('representation',rep)
    timing=read('timing-v1/progress.json');pairs=[]
    for pair in range(3):
        aa=next(q for q in timing['timings'] if q['pair']==pair and not q['on']);bb=next(q for q in timing['timings'] if q['pair']==pair and q['on'])
        pairs.append(dict(pair=pair,off=aa['receipt'],on=bb['receipt'],wall_speedup=aa['receipt']['seconds']/bb['receipt']['seconds']))
    exclusive=defaultdict(float);inclusive=defaultdict(float)
    for line in (r/'timing-v1/continued32/timing.jsonl').read_text().splitlines():
        q=json.loads(line)
        if q['step']>0:exclusive[q['name']]+=q['exclusive_process_cpu_seconds'];inclusive[q['name']]+=q['process_cpu_seconds']
    costs=dict(pairs=pairs,median_wall_speedup=statistics.median(q['wall_speedup'] for q in pairs),
        continued32=timing['continued32'],exclusive_CPU_32_steps=dict(sorted(exclusive.items(),key=lambda kv:-kv[1])),
        inclusive_CPU_32_steps=dict(inclusive),do_not_add_inclusive_to_exclusive=True,
        curved_before=read('curved-final-on/progress.json')['observations']['composite-hill']['receipt'],
        curved_after=read('curved-lazy-final/progress.json')['observations']['composite-hill']['receipt'],
        curved_off=read('curved-final-off/progress.json')['observations']['composite-hill']['receipt'],
        curved_cache=on['evaluation_reuse'],slow_reference=read('operators-v4/progress.json')['failure'],
        distance_kernel=read('build-mapped/cao-tamura-h1-kernels.json')['geometry']['physical_distance_bound_comparison'])
    save('performance',costs)
    final=read('composite-flat-lazy-final/progress.json');previous=read('composite-flat-steps/progress.json')
    require(final['status']=='completed_flat_limit_steps','final composite step check incomplete')
    for key in ('state','fixed_reference','fixed_defect','native_history_times'):
        require(state_hash(read('composite-flat-lazy-final/composite-flat2/case.evolution.json')[key])==
                state_hash(read('composite-flat-steps/composite-flat2/case.evolution.json')[key]),'lazy guard changed continued '+key)
    old=read('composite-flat-steps/composite-flat2/case.evolution.json');new=read('composite-flat-lazy-final/composite-flat2/case.evolution.json')
    require([state_hash(s) for s in old['native_history']]==[state_hash(s) for s in new['native_history']],
            'lazy guard changed native history slots')
    for x,y in zip(old['pimple_attempts'],new['pimple_attempts']):
        for key in ('history_sha256','incoming_sha256','post_flow_sha256','pre_k_sha256','raw_k_sha256','corrected_k_sha256'):
            require(x[key]==y[key],'lazy guard changed governing stage '+key)
    validation_root='rollback-context-final'
    va=read(validation_root+'/progress.json');require(va['status']=='validation_passed','final import/transaction checks incomplete')
    cadence=read('epoch-cadence-final/progress.json');require(cadence['status']=='cadence_passed','epoch cadence test incomplete')
    faults=va['faults']
    for q in faults.values():
        history=q['rollback'].pop('attempted_stage_history',None)
        if history is not None:q['rollback']['attempted_stage_history_size']=len(history)
    verify=dict(flat_import_and_packet_rollback=va['tests'],faults=faults,epoch_packet_cadence=cadence,
        final_composite_steps=final,prior_composite_steps=previous,
        curved_cache_state_audit_inventory_correction_bitwise=True,lazy_guard_state_audit_operators_bitwise=True,
        lazy_guard_native_stages_and_history_bitwise=True,
        independent_spectral_and_gradient_checks=spectra,kernel_checks=read('focused-final-kernels.json'),
        native_seed_rejected_without_PDE=True,curved_split_restart_and_dt_pair='not executed: initial covariance gate failed',
        optional_diagnostics_preserve_native_outcome=True)
    save('verification',verify)
    initial=read('initial-status.json');p2=initial['observer'];latest=p2['diagnostics'][str(max(map(int,p2['diagnostics'])))]
    status=dict(capture_Asia_Shanghai=datetime.now(timezone(timedelta(hours=8))).isoformat(),
        modified_P2=dict(status=initial['campaign']['status'],accepted_steps=initial['campaign']['accepted_steps'],
            accepted_time=initial['campaign']['accepted_time'],checkpoint_sha256=initial['campaign']['checkpoint_sha256'],
            assessed_steps=p2['assessed_steps'],last_assessment=latest,scope='read once at H2b start; runtime/cap/policy unchanged'),
        hill=dict(status='initialization_rejected',accepted_physical_steps=0,attempted_time=0.,reason='inner and outer covariance'),
        flat_temporal=dict(status='completed',physical_start=.0008,physical_end=.0012,steps=70),
        developed_precursor='not qualified',experimental_hill_validation='not attempted; data held out')
    save('native',status)
    save('source-revisions',dict(reviewed_head='5cbaed815ff54d46312af657663aa4f346a5fe40',tested_source_revision=a.source_revision,
        final_solver=va['runtime'],final_probe=read('curved-lazy-final/progress.json')['probe'],
        source_hashes=read('h2b-implementation-hashes.json'),profile_sha256=sha(a.h1/'input/finite-inflow.txt'),
        full_native_source_verification=read('final-full-native-verification.json'),
        storage=dict(actual_root=str(r),historical_mount='/media/jay/新加卷',actual_mount='/media/jay/新加卷1',UUID='5AB8A669B8A642FD',
                     reason='same NTFS volume mounted at a different name after server reboot; no checkpoint bytes/signatures changed'),
        resources=dict(worker_cap=72,build_parallelism=4,native_rank=1,computational_threads=1,maximum_concurrent_native_jobs=4),
        compilation=dict(type='Release',flags='-fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS',fresh_build=str(r/'build-mapped'))))
    # Do not hash the complete source/build tree or copy large native payloads
    # into Git. Every retained run/packet/checkpoint/closure remains indexed.
    roots=['flat-local-v2','flat-final-validation','flat-manifest-final','rollback-receipt-final','rollback-context-final','epoch-cadence-final','timing-v1','representation-v1','mapped-h0-v1','operators-v3','operators-v4','operators-v5',
           'curved-final-on','curved-final-off','curved-nonzero-datum','curved-lazy-final','composite-flat-steps','composite-flat-lazy-final']
    files=[]
    for root in roots:
        for f in sorted((r/root).rglob('*')):
            if f.is_file() and not f.is_symlink():files.append(dict(path=str(f),bytes=f.stat().st_size,sha256=sha(f)))
    for f in sorted(r.glob('*.log')):files.append(dict(path=str(f),bytes=f.stat().st_size,sha256=sha(f)))
    p2path=Path(initial['campaign']['checkpoint'].replace('/media/jay/新加卷/','/media/jay/新加卷1/'))/'State.bin'
    require(sha(p2path)==initial['campaign']['checkpoint_sha256'],'P2 checkpoint differs from captured terminal ledger')
    save('artifact-index',dict(root=str(r),files=files,total_bytes=sum(f['bytes'] for f in files),
        external_immutable_inputs=[dict(path=str(p2path),sha256=sha(p2path)),
            dict(path=str(a.h1/'input/finite-inflow.txt'),sha256=sha(a.h1/'input/finite-inflow.txt'))],
        source_revision=a.source_revision,scope='heavy native checkpoints, stage vectors, matrices, logs and frozen library closures; prior experiments preserved'))
    print(a.output)


if __name__=='__main__':main()
