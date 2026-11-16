#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H25 finite transient initial-data producer, using the compiled strict helper.

This creates new initial data/history, never an evolved-state restart bridge.
The unchanged internal reference/defect and the new physical qh0 stay distinct.
"""
import argparse, copy, hashlib, os, subprocess
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from run_wj_h4 import THREADS, load
os.environ.update(THREADS)
from wj_h15_observer import identity, verify, atomic, vector_hash, assess, BoundedIndex, pkey, add_witnesses
from wj_h16_observer import compatibility
from wj_h17_observer import native_load, load_compatible, equal
from wj_h20_observer import common_points
from wj3_measure import checkpoint_payload
from cls_epawf_m3_cases import set_info, set_parameter, write
from wj_h7_analysis import Rows

DT = (.0025, .005, .01, .02, .04, .08)


def checked_seed(seed, cp, source, perturbed, dt):
    target = seed.get('target_discretization', False)
    assert seed['schema'] == ('WJH25TargetNativeSeedV2' if target else 'WJH25NativeSeedV2') and seed['status'] == 'PASS'
    assert seed['steps'] == cp['steps'] == 0 and seed['time'] == cp['time'] == 0
    assert seed['dt'] == cp['dt'] == dt and seed['signature'] == cp['signature']
    for key in ('state', 'fixed_reference', 'fixed_defect', 'eta', 'native_history', 'native_history_times'):
        equal(seed[key], cp[key], 'new own-dt seed ' + key)
    equal(cp['eta'], source['eta'], 'immutable eta')
    if not target:
        for key in ('fixed_reference', 'fixed_defect'):
            equal(cp[key], source[key], 'immutable original ' + key)
    else:
        assert len(cp['state']) in (12445,24688) and len(source['state']) == 24688
        assert seed['source_checkpoint_sha256'] == 'd27ee0e2f9b9bf8639d47449c5e8e338212363417faca7bf014d1008f44f69a3'
    if not perturbed:
        equal(cp['state'], seed['base_state'] if target else source['state'], 'unperturbed spatial identity')
    else:
        assert vector_hash(cp['state']) != vector_hash(seed['base_state'] if target else source['state'])
        assert seed['scalar_fraction'] == .01 and seed['velocity_scale'] == .001
    assert seed['governing_solves'] == 0
    assert seed['initial_velocity_applications'] == int(target)
    assert seed['native_history_initializations'] == 1
    if target:
        assert seed['spatial_work'][0] == 1 + int(perturbed)
        assert seed['spatial_work'][1] >= seed['spatial_work'][0]
    else:
        assert seed['spatial_work'] == ([1, 1] if perturbed else [0, 0])
    assert seed['perturbed'] is perturbed
    assert seed['last_stage_alpha'] == 0 and seed['last_stage_history'] == []


def rebound_cache(root, seed, directory, runtime):
    """Keep historical pins intact; reuse only proven invariant FE support."""
    prior = load(root/'science/cache-compatible-v2.json')
    c = compatibility(root, verify(prior['historical_manifest']), directory/'cache-base.json')
    c.update(new_native_identity=identity(runtime),
        permitted_difference='H25 own-dt initial history and explicitly new transient qh0. Immutable native IDs/signs/geometry/basis/quadrature/analytical reference only; all initial/current fields recomputed.')
    atomic(directory/'cache-compatible.json', c)
    cache, manifest = load_compatible(directory/'cache-compatible.json')
    for key in ('fields', 'offsets', 'map_parameters', 'geometry_signature', 'parameters', 'quadrature_orders'):
        equal(seed[key], cache.report[key], 'dt-independent row identity ' + key)
    old_reference = copy.deepcopy(cache.reference_identity)
    cache.qh0 = list(seed['state'])
    cache.initial = cache.base.evaluate(cache.qh0, principal=False)
    cache.variable.clear()
    cache.reference_identity = dict(old_reference, qh0=vector_hash(seed['state']),
        physical_signature_sha256=hashlib.sha256(seed['signature'].encode()).hexdigest())
    atomic(directory/'cache-rebind.json', dict(status='PASS', historical=old_reference,
        actual_initial=cache.reference_identity, seed=identity(directory/'native.json'),
        immutable_rebuilds=0, state_reference_recomputed=True,
        rule='new qh0 computed once on all base rows and each required extension; old manifest and pickle never edited'))
    return cache


def run(root, label, dt, perturbed, target_xml=None, runtime=None):
    assert dt in DT and label.replace('-', '').isalnum()
    v = load(root/'inputs/baseline.json')['own_cp0']
    for pin in v.values(): verify(pin)
    d = root/'resumption-r1/science'/label; d.mkdir()
    x = ET.parse(target_xml or v['xml']['path']).getroot()
    runtime = runtime or root/'resumption-r1/science/runtime-build-seed-v1.json'
    record = load(runtime)
    probe = verify(next(e['wrapper'] for e in record['executables'] if Path(e['wrapper']['path']).name=='RANSWallResolutionProbe'))
    for key, value in dict(Time=0, TimeStep=dt, NumSteps=0).items(): set_parameter(x, key, value)
    for key, value in dict(RANSWallNativeLookup='TensorExactV1', RANSWallScalarEvaluation='Off',
        RANSWallRestart='None', RANSWallRetention='Off', RANSWallMonitoringPolicy='Off',
        RANSWallPimpleSnapshots='Off', RANSWallPimpleFlowVectors='Off', RANSWallPimpleArithmetic='Off',
        RANSWallPimpleVerify='Off').items(): set_info(x, key, value)
    write(x, d/'case.xml')
    request = dict(operation='h25_target_seed' if target_xml else 'h25_native_seed', inner_stress_model='Linear', source_dt=.0025,
        source_xml=v['xml'], native_report=v['report']['path'], native_report_sha256=v['report']['sha256'],
        native_checkpoint=str(Path(v['checkpoint']['path']).parent), checkpoint_sha256=v['checkpoint']['sha256'],
        checkpoint_base=str(d/'CP'), perturbed=perturbed)
    atomic(d/'request.json', request)
    atomic(d/'execution-intent.json', dict(xml=identity(d/'case.xml'), request=identity(d/'request.json'),
        runtime=identity(runtime), source=v, target_xml=identity(target_xml) if target_xml else None,
        new_initial_data=True, spatial_projection_calls=int(bool(target_xml))+int(perturbed),
        spatial_lift_calls=int(bool(target_xml))+int(perturbed),initial_velocity_applications=int(bool(target_xml)),
        native_history_initializations=1, governing_attempts=0, later_error_window=[0, .08]))
    with (d/'native.log').open('x') as log:
        code = subprocess.run([str(probe), str(d/'case.xml'),
            str(d/'request.json'), str(d/'native.json')], stdout=log, stderr=subprocess.STDOUT).returncode
    atomic(d/'launch-result.json', dict(exit=code))
    assert code == 0, 'initial-data admission failed; preserve first native failure'
    seed = native_load(d/'native.json')
    cp = checkpoint_payload(seed['checkpoint'], include_reference=True)
    source = checkpoint_payload(verify(v['checkpoint']).parent, include_reference=True)
    checked_seed(seed, cp, source, perturbed, dt)
    fixed = list(load(verify(load(root/'inputs/baseline.json')['support']))['points'])
    if perturbed:
        xt,zt=seed['streamfunction_x_knots'],seed['streamfunction_z_knots']
        normalization=[xt[2],0.,(zt[1]+zt[2])/2]
        scalar_peak=[xt[2],0.,zt[2]]
        fixed.extend((normalization,scalar_peak))
    if target_xml:
        # Rebuild target geometry/basis/owners; never relabel the fine cache.
        prior=load(root/'science/cache-compatible-v2.json')
        c=compatibility(root,verify(prior['historical_manifest']),d/'fine-config-compatibility-base.json')
        c.update(new_native_identity=identity(runtime),permitted_difference='fine context for physical model constants only; all target FE rows freshly built')
        atomic(d/'fine-config-compatibility.json',c)
        fine,_=load_compatible(d/'fine-config-compatibility.json')
        points={pkey(p):p for p in fixed}
        for packet in seed['wj_audit_packets']: add_witnesses(points,packet,fine.report,fine.config)
        index=BoundedIndex(seed)
        rows=Rows(seed,list(points.values()),fine.config,owner_lookup=index,geometry_cache=index)
        initial=rows.evaluate(seed['state'],principal=False)
        stats=dict(target_immutable_rows_rebuilt=len(rows.rows),historical_rows_relabelled=False)
    else:
        cache = rebound_cache(root, seed, d, runtime)
        points = common_points(cache, fixed, {'seed': seed}, {'seed': seed})
        rows, initial, stats = cache.view(points)
    physical, data = assess(rows, seed['state'], initial, 0)
    perturbation=None
    if perturbed:
        base=rows.evaluate(seed['base_state'],principal=False)
        changes={name:float(np.max(np.abs(data[name]-base[name]))) for name in ('U','G','k','epsilon','K','rc')}
        assert all(np.isfinite(v) and v>0 for v in changes.values()), 'transient failed to change required physical fields/current K/cutoff'
        norm_ids=[i for i,(point,owner) in enumerate(rows.rows[:rows.point_count]) if pkey(point)==pkey(normalization)]
        peak_ids=[i for i,(point,owner) in enumerate(rows.rows[:rows.point_count]) if pkey(point)==pkey(scalar_peak)]
        assert norm_ids and peak_ids
        velocity_error=max(abs(float(data['U'][i,0]-base['U'][i,0])-.001) for i in norm_ids)
        k_error=max(abs(float(data['k'][i]-base['k'][i])-.01*rows.config.K)/rows.config.K for i in peak_ids)
        assert velocity_error<=4e-11 and k_error<=4e-11, 'physical finite perturbation normalization differs'
        perturbation=dict(changes=changes,velocity_normalization_point=normalization,
            scalar_peak=scalar_peak,velocity_error=velocity_error,k_relative_error=k_error,
            coordinates='physical nondimensional coordinates of the retained full-height model',
            epsilon='1% local reference scale before required current-K dependent trace lift; trace change reported separately')
    atomic(d/'independent.json', dict(status=physical['summary']['status'], physical=physical, cache=stats,
        seed=identity(d/'native.json'), actual_qh0=vector_hash(seed['state']), new_support_points=len(points),
        finite_perturbation=perturbation,
        distinction='analytical q_eq, new transient qh0, immutable original spatial anchor and fixed internal reference are separate'))
    assert physical['summary']['status'] == 'PASS', 'new initial-data physical/quality obstruction'
    atomic(d/'qualification.json', dict(status='PASS_INITIAL_DATA_ONLY', checkpoint=identity(Path(seed['checkpoint'])/'State.bin'),
        report=identity(d/'native.json'), xml=identity(d/'case.xml'), independent=identity(d/'independent.json'),
        perturbed=perturbed, dt=dt, spatial_original=v['checkpoint'], actual_qh0=vector_hash(seed['state']),
        governing_attempts=0, temporal_accuracy='NOT_RUN'))
    print(label, 'PASS_INITIAL_DATA_ONLY', flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--root', type=Path, required=True)
    p.add_argument('--label', required=True); p.add_argument('--dt', type=float, required=True)
    p.add_argument('--perturbed', action='store_true')
    p.add_argument('--target-xml',type=Path);p.add_argument('--runtime',type=Path)
    a = p.parse_args(); run(a.root, a.label, a.dt, a.perturbed,a.target_xml,a.runtime)
