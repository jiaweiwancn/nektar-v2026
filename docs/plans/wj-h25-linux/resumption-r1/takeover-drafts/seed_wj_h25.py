#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H25 finite transient initial-data producer, using the compiled strict helper.

This creates new initial data/history, never an evolved-state restart bridge.
The unchanged internal reference/defect and the new physical qh0 stay distinct.
"""
import argparse, copy, hashlib, os, subprocess
from pathlib import Path
import xml.etree.ElementTree as ET
from run_wj_h4 import THREADS, load
os.environ.update(THREADS)
from wj_h15_observer import identity, verify, atomic, vector_hash, assess
from wj_h16_observer import compatibility
from wj_h17_observer import native_load, load_compatible, equal
from wj_h20_observer import common_points
from wj3_measure import checkpoint_payload
from cls_epawf_m3_cases import set_info, set_parameter, write

DT = (.0025, .005, .01, .02, .04, .08)


def checked_seed(seed, cp, source, perturbed, dt):
    assert seed['schema'] == 'WJH25NativeSeedV1' and seed['status'] == 'PASS'
    assert seed['steps'] == cp['steps'] == 0 and seed['time'] == cp['time'] == 0
    assert seed['dt'] == cp['dt'] == dt and seed['signature'] == cp['signature']
    for key in ('state', 'fixed_reference', 'fixed_defect', 'eta', 'native_history', 'native_history_times'):
        equal(seed[key], cp[key], 'new own-dt seed ' + key)
    for key in ('fixed_reference', 'fixed_defect', 'eta'):
        equal(cp[key], source[key], 'immutable original ' + key)
    if not perturbed:
        equal(cp['state'], source['state'], 'unperturbed spatial identity')
    else:
        assert vector_hash(cp['state']) != vector_hash(source['state'])
        assert seed['scalar_fraction'] == .01 and seed['velocity_scale'] == .001
    assert seed['governing_solves'] == seed['initial_velocity_applications'] == 0
    assert seed['native_history_initializations'] == 1
    assert seed['spatial_work'] == ([1, 1] if perturbed else [0, 0])
    assert seed['perturbed'] is perturbed
    assert seed['last_stage_alpha'] == 0 and seed['last_stage_history'] == []


def rebound_cache(root, seed, directory):
    """Keep historical pins intact; reuse only proven invariant FE support."""
    prior = load(root/'science/cache-compatible-v2.json')
    c = compatibility(root, verify(prior['historical_manifest']), directory/'cache-base.json')
    c.update(new_native_identity=identity(root/'science/runtime-build-scalar-v1.json'),
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


def run(root, label, dt, perturbed):
    assert dt in DT and label.replace('-', '').isalnum()
    v = load(root/'inputs/baseline.json')['own_cp0']
    for pin in v.values(): verify(pin)
    d = root/'science'/label; d.mkdir()
    x = ET.parse(v['xml']['path']).getroot()
    for key, value in dict(Time=0, TimeStep=dt, NumSteps=0).items(): set_parameter(x, key, value)
    for key, value in dict(RANSWallNativeLookup='TensorExactV1', RANSWallScalarEvaluation='Off',
        RANSWallRestart='None', RANSWallRetention='Off', RANSWallMonitoringPolicy='Off',
        RANSWallPimpleSnapshots='Off', RANSWallPimpleFlowVectors='Off', RANSWallPimpleArithmetic='Off',
        RANSWallPimpleVerify='Off').items(): set_info(x, key, value)
    write(x, d/'case.xml')
    request = dict(operation='h25_native_seed', inner_stress_model='Linear', source_dt=.0025,
        source_xml=v['xml'], native_report=v['report']['path'], native_report_sha256=v['report']['sha256'],
        native_checkpoint=str(Path(v['checkpoint']['path']).parent), checkpoint_sha256=v['checkpoint']['sha256'],
        checkpoint_base=str(d/'CP'), perturbed=perturbed)
    atomic(d/'request.json', request)
    atomic(d/'execution-intent.json', dict(xml=identity(d/'case.xml'), request=identity(d/'request.json'),
        runtime=identity(root/'science/runtime-build-scalar-v1.json'), source=v,
        new_initial_data=True, spatial_projection_calls=int(perturbed), spatial_lift_calls=int(perturbed),
        native_history_initializations=1, governing_attempts=0, later_error_window=[0, .08]))
    with (d/'native.log').open('x') as log:
        code = subprocess.run([str(root/'runtime-scalar-v1/RANSWallResolutionProbe'), str(d/'case.xml'),
            str(d/'request.json'), str(d/'native.json')], stdout=log, stderr=subprocess.STDOUT).returncode
    atomic(d/'launch-result.json', dict(exit=code))
    assert code == 0, 'initial-data admission failed; preserve first native failure'
    seed = native_load(d/'native.json')
    cp = checkpoint_payload(seed['checkpoint'], include_reference=True)
    source = checkpoint_payload(verify(v['checkpoint']).parent, include_reference=True)
    checked_seed(seed, cp, source, perturbed, dt)
    cache = rebound_cache(root, seed, d)
    fixed = load(verify(load(root/'inputs/baseline.json')['support']))['points']
    points = common_points(cache, fixed, {'seed': seed}, {'seed': seed})
    rows, initial, stats = cache.view(points)
    physical, data = assess(rows, seed['state'], initial, 0)
    atomic(d/'independent.json', dict(status=physical['summary']['status'], physical=physical, cache=stats,
        seed=identity(d/'native.json'), actual_qh0=vector_hash(seed['state']), new_support_points=len(points),
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
    a = p.parse_args(); run(a.root, a.label, a.dt, a.perturbed)
