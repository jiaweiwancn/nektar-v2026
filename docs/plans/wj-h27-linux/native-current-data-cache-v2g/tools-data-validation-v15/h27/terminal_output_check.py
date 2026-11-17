#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent Linux terminal bracketing/reduction FOUNDATION, not full H27-D.

Manifest WJH27TerminalFoundationV1 binds scheduled,dense,prefix,cold physical
runs, aligned,initial MMS runs, a failed MMS index, and the pre-launch cold
checkpoint. Original signed owners are retained. No native stress/basis action
or time-reduction routine is called.
"""
import argparse
import bisect
import json
import math
import os
from pathlib import Path
import resource
import sys
import time
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[variable] = '1'
import numpy as np
from numpy.polynomial.legendre import leggauss
import history_recovery_check as hr
import implicit_algorithm_check as fields
import physical_cache_history_check as history
import physical_oracle as po
import physical_transient_oracle as transient

require, read = po.require, history.read
BITS = dict(periodic=1, begin_left=2, begin_right=4, begin_exact=8,
            end_left=16, end_right=32, end_exact=64, final=128)
HASHES = ('state_sha256', 'previous_sha256', 'advection_sha256', 'previous_advection_sha256')
LIMITS = dict(U=1e-4, k=1e-3, Pi=1e-3, R=1e-3, G=1e-2,
              epsilon=1e-2, Ew=1e-2, K=1e-3, Gamma=1e-3, cutoff=1e-3, traction=1e-3)
WINDOW = (.0173, .0537)
SCIENTIFIC_STRIDE = 2


def clock(deadline):
    require(time.monotonic() < deadline, 'terminal checker exhausted cooperative analysis cap')


def input_file(spec, root):
    return history.input_path(spec, root)


def session(spec, root):
    path, evidence = input_file(spec, root)
    xml = ET.parse(path).getroot(); parameters = {}; solver = {}
    for item in xml.findall('.//PARAMETERS/P'):
        key, value = ''.join(item.itertext()).split('=', 1); key = key.strip()
        require(key not in parameters, 'duplicate session parameter'); parameters[key] = value.strip()
    for item in xml.findall('.//SOLVERINFO/I'):
        key = item.attrib['PROPERTY']; require(key not in solver, 'duplicate solver option')
        solver[key] = item.attrib['VALUE']
    return path, evidence, parameters, solver


def expected_schedule(start, count, dt, stride, begin, end):
    """Endpoint location in the complete clock grid, independently of step loop."""
    result = {}; grid = [(start+i)*dt for i in range(count+1)]
    def add(step, mask): result[step] = result.get(step, 0) | mask
    for step in range(start+1, start+count+1):
        if stride > 0 and step % stride == 0: add(step, BITS['periodic'])
    for label, boundary in [('begin', begin), ('end', end)]:
        if boundary < grid[0] or boundary > grid[-1]: continue
        at = bisect.bisect_left(grid, boundary)
        if grid[at] == boundary:
            if count: add(start+at, BITS[label+'_exact'])
        else:
            require(0 < at < len(grid), 'missing interior endpoint brackets')
            add(start+at-1, BITS[label+'_left']); add(start+at, BITS[label+'_right'])
    add(start+count, BITS['final'])
    return result


def record_hash(record):
    return transient.state_hash(transient.state(record['coefficients'], record['sizes']))


def load_run(name, spec, root, control, ref, deadline):
    clock(deadline); path, evidence = input_file(spec['result'], root); packet = read(path); data = packet['state']
    index_path, index_id = input_file(spec['index'], root); index = read(index_path)
    xml, xml_id, parameters, solver = session(spec['session'], root)
    start = 6 if name == 'cold' else 0
    count = 6 if name in ('prefix', 'cold') else 12 if name in ('scheduled', 'dense') else 4
    physical = name in ('scheduled', 'dense', 'prefix', 'cold')
    dt = .005 if physical else .02
    window = WINDOW if physical else ((.02, .06) if name == 'aligned' else (0., .06))
    stride = 1 if name == 'dense' else SCIENTIFIC_STRIDE if physical else 0
    require(index['schema'] == 'WJH27NativeScientificIndexV1' and index['run_complete'] is True,
            'missing or incomplete native scientific index')
    require(index['reason_bits'] == BITS and index['terminal_output'] is True and index['coefficients_included'] is True,
            'scientific output policy differs')
    require((index['terminal_begin'], index['terminal_end']) == window, 'requested terminal interval changed')
    require(index['dt'] == dt and index['initial_step'] == start and index['initial_time'] == start*dt and
            index['final_step'] == start+count and index['final_time'] == (start+count)*dt, 'scientific index clock differs')
    for key in ('source_identity', 'reference_identity', 'signature'):
        require(index[key] == data[key], 'index/native identity differs: '+key)
    require(data['source_identity'] == control['source_identity'] and data['reference_identity'] == po.TABLE_SHA,
            'run source/reference identity differs')
    require(data['manufactured'] is (not physical) and data['dt'] == dt and data['step'] == start+count and data['time'] == (start+count)*dt,
            'run physical mode or endpoint differs')
    for key, value in dict(SolverType='WJNativeCoupledV2', WJNativeMode='Physical' if physical else 'Manufactured',
            WJImplicitSolver='CoupledJFNKV1', WJTerminalOutput='True', WJOutputCoefficients='True').items():
        require(solver[key] == value, 'session solver policy differs: '+key)
    for key, value in dict(NumSteps=count, TimeStep=dt, WJTerminalBegin=window[0], WJTerminalEnd=window[1],
                           WJScientificSteps=stride, IO_InfoSteps=1).items():
        require(float(parameters[key]) == value, 'session parameter differs: '+key)
    prefix = Path(str(path)[:-len('.result.json')])
    require(path.name.endswith('.result.json') and hr.resolve(root, solver['WJOutputPrefix']) == prefix and
            hr.resolve(root, packet['scientific_index']) == index_path, 'output/index prefix binding differs')
    require(hr.identity(hr.resolve(root, solver['WJReferenceTable']))['sha256'] == po.TABLE_SHA, 'session table differs')
    require(bool(solver.get('WJRestart')) is (name == 'cold'), 'session restart mode differs')
    steps = packet['complete_step_seconds']
    require(len(steps) == count and all(math.isfinite(v) and v > 0 for v in steps), 'incomplete accepted-step timing')
    reserve = 128*np.finfo(float).eps
    require(abs(math.fsum(steps)-packet['complete_advancing_seconds']) <= reserve*packet['complete_advancing_seconds'],
            'complete step timing does not close')
    log_path = Path(str(prefix)+'.steps.jsonl'); logs = [fields.native_json(line) for line in log_path.read_text().splitlines()]
    require([row['step'] for row in logs] == list(range(start+1, start+count+1)), 'accepted log is not complete/consecutive')
    previous = None
    for row in logs:
        require(row['time'] == row['step']*dt and all(row[key] == data[key] for key in ('source_identity', 'reference_identity', 'signature')),
                'accepted log identity/clock differs')
        require(all(math.isfinite(v) and 0 <= v <= data['nonlinear_tolerance'] for v in row['residuals']) and
                all(math.isfinite(row[key]) and row[key] > 0 for key in history.GUARDS), 'accepted log residual/admissibility gate failed')
        if previous is not None:
            require(row['previous_sha256'] == previous['state_sha256'] and
                    row['previous_advection_sha256'] == previous['advection_sha256'], 'accepted log history chain differs')
        if physical: history.accepted(row, control, False)
        previous = row
    require(all(logs[-1][key] == data[key] for key in HASHES), 'endpoint and final accepted log differ')
    if physical:
        history.configuration(data, control)
        mesh = history.session_geometry(spec['session'], root, ref)
        owners = history.inventory(data, control, ref, mesh)
    else: mesh, owners = None, None
    expected = expected_schedule(start, count, dt, stride, *window)
    require([frame['step'] for frame in index['frames']] == sorted(expected), 'native endpoint/periodic schedule differs')
    records = {}; frame_ids = []; by_step = {row['step']: row for row in logs}
    for frame in index['frames']:
        clock(deadline); step = frame['step']; frame_path = hr.resolve(root, frame['path']); identity = hr.identity(frame_path)
        require(identity['sha256'] == frame['sha256'] and frame_path.stat().st_size == frame['bytes'] and
                math.isfinite(frame['write_seconds']) and frame['write_seconds'] > 0, 'scientific file identity/timing differs')
        require(frame_path == Path(str(prefix)+'.science-'+str(step)+'.json'), 'scientific owner prefix differs')
        require(frame['reason_mask'] == expected[step] and frame['time'] == step*dt, 'scientific endpoint reason/clock differs')
        record = read(frame_path)
        require(record_hash(record) == record['state_sha256'], 'scientific coefficient SHA differs')
        require(record['step'] == step and record['time'] == frame['time'] and
                all(record[key] == data[key] for key in ('source_identity', 'reference_identity', 'signature', 'sizes', 'dt')),
                'scientific frame metadata differs')
        origin = ('imported_accepted_history' if start else 'fresh_initial_observation') if step == start else 'current_process_accepted_step'
        require(frame['origin'] == origin, 'scientific frame provenance differs')
        if step > start:
            require(all(record[key] == by_step[step][key] for key in (*HASHES, 'residuals', *history.GUARDS)),
                    'scientific frame does not bind its accepted log')
        else:
            require(origin != 'current_process_accepted_step', 'initial/imported observation falsely marked accepted')
        if physical:
            for key in ('geometry', 'elements', 'free_velocity_ids', 'k_free_ids', 'epsilon_free_ids'):
                require(record[key] == control[key], 'scientific original-owner metadata differs: '+key)
        records[step] = record; frame_ids.append(identity)
    actual = {p.resolve() for p in prefix.parent.glob(prefix.name+'.science-[0-9]*.json')}
    require(actual == {hr.resolve(root, f['path']) for f in index['frames']}, 'actual scientific retention differs from index')
    frame_total = math.fsum(f['write_seconds'] for f in index['frames'])
    require(abs(frame_total-packet['scientific_frame_seconds']) <= reserve*frame_total and
            0 < packet['scientific_index_seconds'] <= packet['final_export_seconds'], 'scientific I/O timing accounting differs')
    states, binary = history.binary_history(Path(str(prefix)+'.final.bin'), data, dt, deadline)
    require(binary['step'] == data['step'] and all(binary['hashes'][role] == data[key] for role, key in
        zip(history.HISTORY_NAMES, HASHES)), 'complete endpoint binary/history differs')
    return dict(packet=packet, data=data, index=index, records=records, logs=logs, mesh=mesh, owners=owners,
                histories=history.exported(data, states) if physical else states, session_solver=solver,
                evidence=dict(result=evidence, index=index_id, session=xml_id, log=hr.identity(log_path),
                    frames=frame_ids, checkpoint=binary, expected_schedule=expected,
                    complete_advancing_seconds=packet['complete_advancing_seconds'],
                    complete_step_seconds=steps, scientific_frame_seconds=frame_total,
                    scientific_index_seconds=packet['scientific_index_seconds']))


class Sampler:
    """Use only the independent owner's basis/stress sampler; never snap bounds."""
    trace = po.PhysicalOracle.trace
    sample = po.PhysicalOracle.sample

    def __init__(self, record, owners, ref, deadline):
        self.ref, self.owner, self.deadline = ref, owners, deadline
        self.n = record['sizes'][0]; self.s = transient.state(record['coefficients'], record['sizes'])

    def clock(self): clock(self.deadline)


def observations(record, owners, mesh, ref, deadline):
    sampler = Sampler(record, owners, ref, deadline); groups = {}
    a, gauss_weight = leggauss(5); nodes = np.r_[a, -1., 1.]; weights = np.r_[gauss_weight, 0., 0.]
    r_tolerance = mesh['r_tolerance']
    cutoff0 = float(ref.cutoff(np.array([ref.Kin]))[0])
    def add(region, field, value, scale, e, x, r, weight):
        key = region+'_'+field; v = np.asarray(value, dtype=float).reshape(len(x), -1)
        scale = np.broadcast_to(np.asarray(scale, dtype=float).reshape(-1, 1), (len(x), v.shape[1]))
        require(np.all(np.isfinite(v)) and np.all(np.isfinite(scale)) and np.all(scale > 0), 'invalid physical observation/scale')
        g = groups.setdefault(key, dict(region=region, field=field, value=[], scale=[], x=[], r=[], weight=[], owner=[]))
        for name, block in dict(value=v/scale, scale=scale, x=x, r=r, weight=weight,
                                owner=np.full(len(x), e['geometry'], dtype=np.int64)).items(): g[name].append(block)
    for e in sorted(record['elements'], key=lambda v: v['geometry']):
        clock(deadline); xl, xr, rl, rr = e['bounds']
        x, r = np.meshgrid((xl+xr)/2+(xr-xl)*nodes/2, (rl+rr)/2+(rr-rl)*nodes/2)
        x, r = x.ravel(), r.ravel(); weight = (weights[:, None]*weights[None, :]*(xr-xl)*(rr-rl)/4).ravel()
        value = sampler.sample(e, x, r); region = 'outer' if e['epsilon'] else 'inner'
        require(np.all(value['covariance'] > 0), 'independently sampled scientific covariance is nonpositive')
        gradient_scale = ref.utau/(ref.kappa*r) if e['epsilon'] else ref.tau/(ref.nu+ref.aux(np.full(len(r), ref.Kin), r))
        for field, observed, scale in [('U', np.column_stack((value['U'], value['W'])), 1.),
                ('k', value['k'], ref.Kin), ('Pi', value['Pi'], ref.tau),
                ('R', value['R'], ref.tau), ('G', value['G'], gradient_scale)]:
            add(region, field, observed, scale, e, x, r, weight)
        add(region, 'epsilon' if e['epsilon'] else 'Ew', value['E'],
            ref.A/r if e['epsilon'] else ref.A/np.maximum(r, cutoff0), e, x, r, weight)
        at_h = abs(rl-ref.H) <= r_tolerance and bool(e['epsilon'])
        at_ground = abs(rl) <= r_tolerance
        if at_h or at_ground:
            xx = (xl+xr)/2+(xr-xl)*nodes/2; radii = np.full(len(xx), rl)
            face = sampler.sample(e, xx, radii); face_weight = weights*(xr-xl)/2
            if at_h:
                for field, observed, scale in [('K', face['K'], ref.Kin),
                        ('Gamma', face['K']**1.5/ref.cl/ref.H, ref.A/ref.H),
                        ('cutoff', face['cutoff'], ref.H)]:
                    add('interface', field, observed, scale, e, xx, radii, face_weight)
            if at_ground:
                traction = np.column_stack((face['V'][:, 0, 2], face['V'][:, 2, 2]-face['Pi']))
                add('ground', 'traction', traction, ref.tau, e, xx, radii, face_weight)
    for group in groups.values():
        for key in ('value', 'scale', 'x', 'r', 'weight', 'owner'): group[key] = np.concatenate(group[key], axis=0)
        require(float(np.sum(group['weight'])) > 0, 'zero spatial observation measure')
    require(len(groups) == 16, 'incomplete observation inventory')
    return groups


def clipped_segments(times, begin, end):
    require(len(times) >= 2 and all(b > a for a, b in zip(times[:-1], times[1:])) and
            times[0] <= begin < end <= times[-1], 'terminal window has incomplete/gapped clock coverage')
    segments = []
    for i, (left, right) in enumerate(zip(times[:-1], times[1:])):
        lo, hi = max(left, begin), min(right, end)
        if hi > lo: segments.append(dict(i=i, left=left, right=right, lo=lo, hi=hi,
            weight=hi-lo, lower_fraction=(lo-left)/(right-left), upper_fraction=(hi-left)/(right-left)))
    require(segments[0]['lo'] == begin and segments[-1]['hi'] == end and
            all(a['hi'] == b['lo'] for a, b in zip(segments[:-1], segments[1:])), 'terminal coverage is discontinuous')
    require(abs(math.fsum(s['weight'] for s in segments)-(end-begin)) <= 64*np.finfo(float).eps*(end-begin),
            'terminal integration weights shorten/change the requested interval')
    return segments


def reduce(records, owners, mesh, ref, deadline, path):
    ordered = sorted(records); times = [records[step]['time'] for step in ordered]
    segments = clipped_segments(times, *WINDOW); length = WINDOW[1]-WINDOW[0]
    frames = [observations(records[step], owners, mesh, ref, deadline) for step in ordered]
    first = frames[0]
    for frame in frames[1:]:
        require(frame.keys() == first.keys(), 'observation fields changed')
        for key in first:
            require(all(np.array_equal(first[key][name], frame[key][name]) for name in ('scale', 'x', 'r', 'weight', 'owner')),
                    'observation owners/locations/scales changed')
    arrays = {}; result = {}; checks = {}; gauss, gw = leggauss(3)
    for key, group in first.items():
        clock(deadline); mean_integral = np.zeros_like(group['value']); second_integral = mean_integral.copy()
        gauss_mean = mean_integral.copy(); gauss_second = mean_integral.copy(); abs_mean = mean_integral.copy(); abs_second = mean_integral.copy()
        endpoints = []
        for s in segments:
            a, b = frames[s['i']][key]['value'], frames[s['i']+1][key]['value']
            lo = a+s['lower_fraction']*(b-a); hi = a+s['upper_fraction']*(b-a); w = s['weight']
            mean_integral += w*(lo+hi)/2
            second_integral += w*((lo+hi)**2+lo**2+hi**2)/6
            abs_mean += w*(abs(lo)+abs(hi))/2; abs_second += w*(lo**2+abs(lo*hi)+hi**2)/3
            for node, weight in zip(gauss, gw):
                value = a+(((s['lo']+s['hi'])/2+w*node/2-s['left'])/(s['right']-s['left']))*(b-a)
                gauss_mean += weight*w*value/2; gauss_second += weight*w*value*value/2
            endpoints.append((w, lo, hi))
        rounding = 256*np.finfo(float).eps
        require(np.all(abs(mean_integral-gauss_mean) <= rounding*abs_mean) and
                np.all(abs(second_integral-gauss_second) <= rounding*abs_second), 'analytic/time-Gauss first/second moment disagreement')
        mean = mean_integral/length; variance = np.zeros_like(mean)
        for w, lo, hi in endpoints:
            a, b = lo-mean, hi-mean; variance += w*((a+b)**2+a*a+b*b)/6/length
        require(np.all(variance >= 0) and np.all(np.isfinite(variance)), 'nonpositive/nonfinite integrated centered variance')
        deviation = np.sqrt(variance)
        result[key] = dict(mean=mean, sd=deviation, group=group)
        for name, array in dict(mean_scaled=mean, standard_deviation_scaled=deviation,
                second_moment_scaled=second_integral/length, scale=group['scale'], x=group['x'], r=group['r'],
                weight=group['weight'], owner=group['owner']).items(): arrays[key+'__'+name] = array
        checks[key] = dict(points=len(mean), components=mean.shape[1], spatial_measure=float(np.sum(group['weight'])),
            first_moment_time_Gauss_MAX=float(np.max(abs(mean_integral-gauss_mean))),
            second_moment_time_Gauss_MAX=float(np.max(abs(second_integral-gauss_second))),
            centered_second_identity_MAX=float(np.max(abs(second_integral/length-mean*mean-variance))))
    require(not path.exists(), 'reduced array output is immutable')
    with path.open('xb') as stream: np.savez_compressed(stream, **arrays)
    with np.load(path, allow_pickle=False) as loaded:
        require(set(loaded.files) == set(arrays) and all(np.array_equal(loaded[key], value) for key, value in arrays.items()),
                'scientific reduction changed an exported array')
    return result, dict(artifact=hr.identity(path), bytes=path.stat().st_size, interval=list(WINDOW),
        window_measure=length, frames=[dict(step=step, time=records[step]['time']) for step in ordered],
        clipped_segments=segments, checks=checks, physical_arrays_recoverable='mean_scaled*fixed_scale; sd_scaled*fixed_scale',
        raw_scientific_frames_retained=True, continuum_supremum_claim=False, spatial_quadrature_qualified=False)


def compare(a, b, limits):
    norms = po.Norms()
    for key in a:
        g = a[key]['group']
        for moment in ('mean', 'sd'):
            norms.add(moment+'/'+g['region'], g['field'], a[key][moment]-b[key][moment], g['weight'])
    result = norms.report(limits)
    return result


def run(manifest, folder, table, deadline, out, output):
    require(manifest['schema'] == 'WJH27TerminalFoundationV1' and tuple(manifest['window']) == WINDOW and
            manifest['scientific_stride'] == SCIENTIFIC_STRIDE, 'wrong terminal manifest/window/cadence')
    root = hr.resolve(folder, manifest.get('root', '.')); ref = po.Reference(read(table))
    require(hr.identity(table)['sha256'] == po.TABLE_SHA, 'unqualified reference table')
    names = {'scheduled', 'dense', 'prefix', 'cold', 'aligned', 'initial'}
    require(set(manifest['runs']) == names, 'incomplete terminal run matrix')
    scheduled_path, _ = input_file(manifest['runs']['scheduled']['result'], root); control = read(scheduled_path)['state']
    require(control['source_identity'] == manifest['source_identity'] and control['manufactured'] is False, 'wrong source or physical control')
    runs = {name: load_run(name, manifest['runs'][name], root, control, ref, deadline) for name in sorted(names)}
    out['runs'] = {name: run['evidence'] for name, run in runs.items()}
    source_keys = ('geometry', 'elements', 'free_velocity_ids', 'k_free_ids', 'epsilon_free_ids')
    for name in ('dense', 'prefix', 'cold'):
        require(all(runs[name]['data'][key] == control[key] for key in source_keys), 'cross-run original ownership differs')
        require(runs[name]['mesh']['geometry_sha256'] == runs['scheduled']['mesh']['geometry_sha256'], 'session meshes differ')
    for step, record in runs['scheduled']['records'].items():
        dense = runs['dense']['records'][step]
        require(all(record[key] == dense[key] for key in HASHES), 'output schedule changed the uninterrupted trajectory')
    for name in ('aligned', 'initial'):
        require(all(runs[name]['data'][key] == runs['aligned']['data'][key] for key in HASHES), 'MMS output interval changed trajectory')
    checkpoint_path, checkpoint_id = input_file(manifest['cold_input'], root)
    states, binary = history.binary_history(checkpoint_path, control, .005, deadline, checkpoint_id['sha256'])
    prefix = runs['prefix']; cold = runs['cold']
    require(binary['step'] == 6 and checkpoint_path == Path(str(prefix['evidence']['result']['path'])[:-len('.result.json')]+'.final.bin') and
            hr.resolve(root, cold['session_solver']['WJRestart']) == checkpoint_path, 'cold process did not import declared prefix6 checkpoint')
    require(all(binary['hashes'][role] == prefix['data'][key] for role, key in zip(history.HISTORY_NAMES, HASHES)), 'restart boundary histories differ')
    require(cold['logs'][0]['previous_sha256'] == binary['hashes']['current'] and
            cold['logs'][0]['previous_advection_sha256'] == binary['hashes']['advection'], 'cold first accepted step lost imported history')
    endpoint = {}
    for role, values in runs['scheduled']['histories'].items():
        require(role in cold['histories'], 'cold endpoint lost an accepted/explicit history role')
        endpoint[role] = fields.coefficient_delta(values, cold['histories'][role], control['sizes'], ref)
        require(all(v['scaled_MAX'] <= 1e-8 for v in endpoint[role].values()), 'cold endpoint exceeds fixed history coefficient tolerance: '+role)
    ownership = manifest['segment_ownership']
    require(ownership == [dict(run='prefix', first_step=0, last_step=6, interval='[0,6]'),
                          dict(run='cold', first_step=6, last_step=12, interval='(6,12]')],
            'restart segment ownership differs')
    stitched = dict(prefix['records']); duplicates = []
    for step, record in cold['records'].items():
        require(6 <= step <= 12, 'cold index contains an abandoned/past branch')
        if step in stitched:
            require(step == 6 and all(record[key] == stitched[step][key] for key in HASHES), 'conflicting duplicate restart frame')
            duplicates.append(step)
        else: stitched[step] = record
    prefix_weight, cold_weight = 6*.005-WINDOW[0], WINDOW[1]-6*.005
    require(prefix_weight > 0 and cold_weight > 0, 'restart does not split the requested terminal interval')
    out['restart_boundary'] = dict(checkpoint=binary, endpoint_scaled_coefficients=endpoint,
        duplicate_boundary_frames=duplicates, ownership=ownership,
        prefix_window_weight=prefix_weight, cold_window_weight=cold_weight,
        both_processes_contribute=True, complete_history_link=True)
    equation = history.equation_report(manifest['current_equation_report'], root,
        Path(cold['evidence']['result']['path']), cold['data'], cold['histories'], po.TABLE_SHA)
    out['independent_current_equation_report'] = equation
    # A real failed step must publish an incomplete index of accepted history.
    failed_path, failed_id = input_file(manifest['failed_index'], root); failed = read(failed_path)
    failure_path, failure_id = input_file(manifest['failed_run'], root); failure = read(failure_path)
    require(failed['schema'] == 'WJH27NativeScientificIndexV1' and failed['run_complete'] is False and
            failed['initial_step'] == failed['final_step'] == 0 and failed['initial_time'] == failed['final_time'] == 0 and
            failure['status'] == 'FAILED_WITH_ACCEPTED_HISTORY_PRESERVED' and failure['state']['step'] == 0 and
            failed['source_identity'] == control['source_identity'] == failure['state']['source_identity'],
            'failed step index falsely claims accepted progress')
    require(len(failed['frames']) == 1 and failed['frames'][0]['step'] == 0 and
            failed['frames'][0]['origin'] == 'fresh_initial_observation', 'failed initial observation is falsely accepted')
    ff = failed['frames'][0]; fp = hr.resolve(root, ff['path']); fr = read(fp)
    require(hr.identity(fp)['sha256'] == ff['sha256'] and record_hash(fr) == failure['state']['state_sha256'] and
            fr['state_sha256'] == failure['state']['state_sha256'], 'failed preserved initial frame identity differs')
    out['failed_step'] = dict(index=failed_id, failure=failure_id, frame=hr.identity(fp), accepted_steps=0)
    try: clipped_segments([.015, .03], *WINDOW)
    except ValueError: out['incomplete_window_control'] = 'REJECTED'
    else: raise ValueError('incomplete terminal coverage was accepted')
    reductions = {}; data = {}
    for name, records in [('scheduled', runs['scheduled']['records']), ('dense', runs['dense']['records']), ('stitched', stitched)]:
        data[name], reductions[name] = reduce(records, runs['scheduled']['owners'], runs['scheduled']['mesh'], ref,
            deadline, output.with_name(output.stem+'-'+name+'.npz'))
    out['reductions'] = reductions
    out['comparisons'] = dict(scheduled_vs_dense=compare(data['scheduled'], data['dense'], LIMITS),
        cold_segments_vs_uninterrupted=compare(data['stitched'], data['scheduled'], {key: 1e-8 for key in LIMITS}))
    require(all(v['status'] == 'MET' for rows in out['comparisons'].values() for v in rows.values()),
            'scheduled/dense or restart terminal statistic exceeds prospective limit; numeric comparisons retained')
    out.update(status='PASS_TERMINAL_OUTPUT_FOUNDATION_ONLY', foundation_gate_pass=True,
        exact_uninterrupted_hashes_at_common_clocks=True, original_owner_bounds_preserved=True,
        unapproved_v4_implemented=False, full_H27_qualification=False, H27_D_pass=False,
        long_pilots_qualified=False, production_storage_qualified=False, monitoring_escalation_qualified=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('manifest', 'table', 'out'): parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--seconds', type=float, required=True); args = parser.parse_args(); began = time.monotonic()
    if args.out.exists(): parser.error('immutable output exists')
    out = dict(schema='WJH27TerminalOutputFoundationCheckV1', status='FAIL', foundation_gate_pass=False,
        full_H27_qualification=False, H27_D_pass=False, TASK3_READY=False, workers=1,
        prospective_limits=LIMITS, window=list(WINDOW), scientific_stride=SCIENTIFIC_STRIDE,
        manifest=hr.identity(args.manifest), table=hr.identity(args.table),
        sources={key: hr.identity(path) for key, path in dict(checker=__file__, reference_stress_basis=po.__file__,
            histories=history.__file__, fields=fields.__file__, transient=transient.__file__, JSON_reader=hr.__file__).items()})
    try:
        require(sys.platform.startswith('linux') and sys.byteorder == 'little', 'authoritative little-endian Linux analysis required')
        require(math.isfinite(args.seconds) and 0 < args.seconds <= 590, 'explicit analysis cap must be in (0,590]')
        soft, hard = resource.getrlimit(resource.RLIMIT_AS)
        cap = min([8*1024**3]+[v for v in (soft, hard) if v != resource.RLIM_INFINITY]); resource.setrlimit(resource.RLIMIT_AS, (cap, hard))
        out['memory_cap_bytes'] = cap
        run(read(args.manifest), args.manifest.resolve().parent, args.table.resolve(), began+args.seconds, out, args.out)
    except Exception as error: out.update(exception=type(error).__name__, reason=str(error))
    out['elapsed_seconds'] = time.monotonic()-began; out['maximum_resident_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
    with args.out.open('x') as stream: json.dump(out, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps({key: out[key] for key in ('status', 'reason', 'elapsed_seconds') if key in out}))
    if not out['foundation_gate_pass']: raise SystemExit(2)


if __name__ == '__main__': main()
