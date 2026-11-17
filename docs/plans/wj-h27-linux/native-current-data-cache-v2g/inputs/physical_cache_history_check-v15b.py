#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent Linux-only physical inverse-cache history FOUNDATION check.

WJH27PhysicalCacheHistoryV1: runs maps control,prefix,cold,recovery9,reload12
and optional retry to .result.json paths; input_checkpoints maps cold,
recovery9,reload12 to {path,sha256}; current_equation_reports maps each of
control,cold,recovery9 and optional retry to an independent transient report
path or {path,sha256}. sessions maps every run to its exact XML path or
{path,sha256}. root defaults to the manifest directory. No launches.
The fresh binary reader below takes dt explicitly; no imported DT is mutated.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import struct
import sys
import time
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[name] = '1'
import numpy as np
import history_recovery_check as hr
import implicit_algorithm_check as fields
import physical_oracle as po
import physical_transient_oracle as transient

DT = .005
FINAL_STEP = 12
LIMIT = 1e-8
STARTS = dict(control=0, prefix=0, cold=3, recovery9=9, reload12=12, retry=0)
HISTORY_NAMES = ('current', 'previous', 'advection', 'previous_advection')
HASH_NAMES = ('state_sha256', 'previous_sha256', 'advection_sha256', 'previous_advection_sha256')
GUARDS = ('minimum_k', 'minimum_epsilon', 'minimum_covariance_over_k', 'minimum_principal_over_mu')
TIMERS = ('constitutive', 'wall', 'assembly', 'flow', 'scalar', 'rebuild',
          'guards', 'forcing', 'residual', 'nonlinear', 'io', 'other')
require = hr.require


def clock(deadline):
    if time.monotonic() >= deadline: raise TimeoutError('physical cache/history check reached its admitted cap')


def read(path):
    return hr.native_json(Path(path).read_text())


def replace_suffix(path, old, new):
    require(path.name.endswith(old), 'invalid native artifact suffix: '+str(path))
    return path.with_name(path.name[:-len(old)]+new)


def input_path(spec, root):
    path = hr.resolve(root, spec if isinstance(spec, str) else spec['path'])
    evidence = hr.identity(path)
    if isinstance(spec, dict): require(evidence['sha256'] == spec['sha256'], 'input SHA256 differs: '+str(path))
    return path, evidence


def scientific_time(step, value, dt, label):
    require(type(step) is int and step >= 0 and math.isfinite(value) and
            abs(value-step*dt) <= 1e-13, 'scientific clock differs: '+label)


def session(spec, root, name, result, checkpoint, table_sha, cache_policy='OperatorDriftV1', terminal_zero=False):
    path, evidence = input_path(spec, root)
    xml = ET.parse(path).getroot(); parameters = {}; solver = {}
    for item in xml.findall('.//PARAMETERS/P'):
        key, value = ''.join(item.itertext()).strip().split('=', 1); key = key.strip()
        require(key not in parameters, 'duplicate session parameter'); parameters[key] = value.strip()
    for item in xml.findall('.//SOLVERINFO/I'):
        key = item.attrib['PROPERTY']; require(key not in solver, 'duplicate session solver option')
        solver[key] = item.attrib['VALUE']
    steps = (3 if name == 'prefix' else FINAL_STEP)-STARTS[name]
    require(float(parameters['NumSteps']) == steps and float(parameters['TimeStep']) == DT and
            float(parameters['IO_InfoSteps']) == 1, 'session executed steps/dt/info stride differs')
    for key, value in dict(SolverType='WJNativeCoupledV2', WJNativeMode='Physical',
            WJImplicitSolver='CoupledJFNKV1', WJStepPreconditioner=cache_policy,
            WJOutputCoefficients='True', WJNonlinearTrace='True').items():
        require(solver[key] == value, 'session numerical identity differs: '+key)
    require(float(parameters.get('WJStepPCRelativeDrift', '.02')) == .02 and
            float(parameters.get('WJStepPCMaximumUses', '8')) == 8, 'session inverse-cache thresholds differ')
    require(hr.resolve(root, solver['WJOutputPrefix']) == replace_suffix(result, '.result.json', ''), 'session output prefix differs')
    require(hr.identity(hr.resolve(root, solver['WJReferenceTable']))['sha256'] == table_sha, 'session reference table differs')
    if checkpoint is None: require(not solver.get('WJRestart', ''), 'fresh session unexpectedly restarts')
    else: require(hr.resolve(root, solver['WJRestart']) == Path(checkpoint['path']), 'session restart input differs')
    require(float(parameters.get('WJInjectRetryStep', '-1')) == (4 if name == 'retry' else -1), 'session retry injection differs')
    if name == 'control':
        require(float(parameters['IO_CheckSteps']) == 3 and float(parameters['WJScientificSteps']) == 3,
                'control recovery/science stride differs')
    if terminal_zero:
        require(name == 'reload12' and solver.get('WJTerminalOutput') == 'True' and
                float(parameters['WJTerminalBegin']) == .0173 and float(parameters['WJTerminalEnd']) == .0537,
                'zero-step terminal export session differs')
    else: require(solver.get('WJTerminalOutput', 'False') == 'False', 'unexpected terminal output changes control retention')
    return dict(file=evidence, NumSteps=steps, restart=None if checkpoint is None else checkpoint['path'],
                output_prefix=str(replace_suffix(result, '.result.json', '')))


def axis_index(axis, value, tolerance, label):
    matches = [i for i, coordinate in enumerate(axis) if abs(coordinate-value) <= tolerance]
    require(len(matches) == 1, 'missing/ambiguous XML interval endpoint: '+label)
    return matches[0]


def session_geometry(spec, root, ref):
    """Authoritative flat tensor inventory from XML; never changes native bounds."""
    path, evidence = input_path(spec, root); geometry = ET.parse(path).getroot().find('GEOMETRY')
    require(geometry is not None and geometry.get('DIM') == geometry.get('SPACE') == '2', 'expected planar session geometry')
    vertices = geometry.find('VERTEX'); edges = geometry.find('EDGE'); elements = geometry.find('ELEMENT')
    require(vertices is not None and edges is not None and elements is not None and not vertices.attrib,
            'expected explicit untransformed session vertices/edges/elements')
    curved = geometry.find('CURVED'); require(curved is None or len(curved) == 0, 'curved session geometry is unsupported')
    raw = {}
    for item in vertices:
        identifier = int(item.attrib['ID']); coordinates = tuple(map(float, ''.join(item.itertext()).split()))
        require(item.tag == 'V' and identifier not in raw and len(coordinates) == 3 and
                all(math.isfinite(v) for v in coordinates), 'invalid/duplicate XML vertex')
        raw[identifier] = coordinates
    require(raw, 'empty XML vertex inventory')
    xvalues = sorted({v[0] for v in raw.values()}); rvalues = sorted({v[1]-ref.d0 for v in raw.values()})
    xspan, rspan = xvalues[-1]-xvalues[0], rvalues[-1]-rvalues[0]
    require(xspan > 0 and rspan > 0, 'degenerate XML geometry')
    xtol, rtol = 64*np.finfo(float).eps*xspan, 64*np.finfo(float).eps*rspan
    def classify_axis(values, tolerance):
        axis = []
        for value in values:
            if not axis or value-axis[-1] > tolerance: axis.append(value)
        require(all(b-a > 2*tolerance for a, b in zip(axis[:-1], axis[1:])), 'ambiguous XML axis intervals')
        return axis
    xaxis, raxis = classify_axis(xvalues, xtol), classify_axis(rvalues, rtol)
    nx, nr = len(xaxis)-1, len(raxis)-1
    require(nx*nr == 700 and abs(xaxis[0]+10) <= xtol and abs(xaxis[-1]-15) <= xtol and
            abs(raxis[0]) <= rtol and abs(raxis[-1]-ref.top) <= rtol, 'XML is not the retained 700-cell physical domain')
    hindex = axis_index(raxis, ref.H, rtol, 'H'); require(0 < hindex < nr, 'H must be an interior native edge')
    positions = {}
    for identifier, (x, z, transverse) in raw.items():
        require(abs(transverse) <= max(xtol, rtol), 'nonplanar XML vertex')
        positions[identifier] = (axis_index(xaxis, x, xtol, 'vertex x'), axis_index(raxis, z-ref.d0, rtol, 'vertex r'))
    require(len(set(positions.values())) == len(positions) == (nx+1)*(nr+1) and
            set(positions.values()) == {(i, j) for i in range(nx+1) for j in range(nr+1)}, 'XML vertices do not cover tensor intersections uniquely')
    edge_map = {}
    for item in edges:
        identifier = int(item.attrib['ID']); ids = tuple(map(int, ''.join(item.itertext()).split()))
        require(item.tag == 'E' and identifier not in edge_map and len(ids) == 2 and len(set(ids)) == 2 and
                all(v in positions for v in ids), 'invalid/duplicate XML edge')
        a, b = (positions[v] for v in ids)
        require(abs(a[0]-b[0])+abs(a[1]-b[1]) == 1, 'XML edge is not one adjacent axis interval')
        edge_map[identifier] = ids
    cells = {}
    for item in elements:
        identifier = int(item.attrib['ID']); ids = tuple(map(int, ''.join(item.itertext()).split()))
        require(item.tag == 'Q' and identifier not in cells and len(ids) == len(set(ids)) == 4 and
                all(v in edge_map for v in ids), 'invalid/duplicate XML quadrilateral')
        points = {positions[v] for edge in ids for v in edge_map[edge]}
        require(len(points) == 4, 'XML quadrilateral has invalid endpoint inventory')
        i, j = min(v[0] for v in points), min(v[1] for v in points)
        corners = {(i, j), (i+1, j), (i, j+1), (i+1, j+1)}
        require(points == corners and len({frozenset(edge_map[v]) for v in ids}) == 4,
                'XML quadrilateral does not cover exactly one tensor cell')
        cells[identifier] = (i, j)
    require(len(cells) == len(set(cells.values())) == 700 and
            set(cells.values()) == {(i, j) for i in range(nx) for j in range(nr)}, 'XML cells do not cover all 700 intervals uniquely')
    signature = hashlib.sha256(json.dumps(dict(vertices=raw, edges=edge_map, cells=cells),
        sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
    return dict(session=evidence, geometry_sha256=signature, x=xaxis, r=raxis, nx=nx, nr=nr,
                x_tolerance=xtol, r_tolerance=rtol, H_index=hindex, cells=cells)


def binary_history(path, reference, dt, deadline, expected_sha=None):
    """Independent fixed-layout reader: unsigned64, IEEE64, four full states."""
    clock(deadline); data = path.read_bytes(); evidence = hr.identity(path)
    require(64 < len(data) < 128*1024**2, 'invalid binary checkpoint size')
    if expected_sha is not None: require(evidence['sha256'] == expected_sha, 'pre-launch checkpoint SHA differs')
    payload, trailer = data[:-64], data[-64:]
    require(hashlib.sha256(payload).hexdigest().encode('ascii') == trailer, 'independent binary checksum mismatch')
    magic, signature, body = payload.split(b'\n', 2)
    require(magic == b'WJ_H27_NATIVE_COUPLED_V2' and signature.decode('ascii') == reference['signature'],
            'checkpoint model/algorithm/source signature differs')
    cursor = 0
    def unpack(fmt):
        nonlocal cursor
        length = struct.calcsize(fmt); require(cursor+length <= len(body), 'truncated checkpoint metadata')
        result = struct.unpack_from(fmt, body, cursor); cursor += length; return result
    step, when, stored_dt = unpack('<Qdd')
    require(stored_dt == dt, 'checkpoint dt differs'); scientific_time(step, when, dt, 'checkpoint')
    sizes = reference['sizes']; histories = {}
    for name in HISTORY_NAMES:
        value = {}
        for key, expected in zip(('flow', 'k', 'epsilon'), (2*sizes[0]+sizes[1], sizes[2], sizes[3])):
            length, = unpack('<Q'); require(length == expected and cursor+8*length <= len(body), 'checkpoint field inventory differs')
            value[key] = np.frombuffer(body, dtype='<f8', count=length, offset=cursor).copy(); cursor += 8*length
            require(np.all(np.isfinite(value[key])), 'nonfinite checkpoint coefficients')
        histories[name] = value
    require(cursor == len(body), 'checkpoint has trailing data')
    return histories, dict(**evidence, payload_sha256=trailer.decode('ascii'), signature=signature.decode('ascii'),
        step=step, time=when, dt=dt, hashes={key: transient.state_hash(value) for key, value in histories.items()})


def configuration(record, control):
    for key in ('schema', 'model', 'source_identity', 'reference_identity', 'signature', 'sizes', 'gauge', 'manufactured',
                'basis', 'quadratureX', 'quadratureR', 'linear_tolerance', 'nonlinear_tolerance', 'maximum_iterations',
                'anderson_depth', 'stress_split', 'constitutive_backend', 'steady_velocity_preconditioner',
                'mms_variant', 'mms_body_amplitude', 'implicit_solver', 'dt'):
        require(record[key] == control[key], 'native identity/configuration differs: '+key)
    require(record['schema'] == 'WJH27NativeCoupledStateV2' and record['model'] == 'FiniteNuIsotropicWallFullWJCollarV2' and
            record['manufactured'] is False and record['gauge'] is False and record['basis'] == 'Modified_A' and
            record['implicit_solver'] == 'CoupledJFNKV1', 'not active physical v2 with open traction outlet')
    require(record['dt'] == DT and record['nonlinear_tolerance'] == 3e-12 and record['linear_tolerance'] == 1e-12 and
            record['anderson_depth'] == 0 and record['mms_body_amplitude'] == 0., 'prospective physical tolerance/time policy differs')
    for key in ('source_identity', 'reference_identity', 'signature'):
        require(len(record[key]) == 64 and all(c in '0123456789abcdef' for c in record[key]), 'missing full '+key)
    cap = record['coupled_newton']['step_seconds_cap']
    require(math.isfinite(cap) and cap > 0 and cap == control['coupled_newton']['step_seconds_cap'], 'step cap binding differs')
    cache = record['coupled_newton']
    expected_policy = control['coupled_newton']['preconditioner_policy']
    require(expected_policy in ('OperatorDriftV1', 'CurrentDataDriftV1'), 'unsupported physical inverse-cache policy')
    fields.cache_policy_statistics(record, expected_policy)
    require(cache['preconditioner_policy'] == expected_policy and cache['relative_drift_limit'] == .02 and
            cache['maximum_uses'] == 8 and cache['preconditioner_events_recorded'] is True,
            'fixed numerical cache policy/event recording differs')
    require(cache['uncommitted_guard'] == 'positivity/root/covariance' and
            cache['physical_acceptance_guard'] == 'fresh full planar principal and complete current residual',
            'current scientific guard provenance differs')


def inventory(record, control, ref, mesh):
    # Exact cross-run ownership is checked BEFORE any geometric grouping.
    for key in ('geometry', 'elements', 'free_velocity_ids', 'k_free_ids', 'epsilon_free_ids'):
        require(record[key] == control[key], 'physical signed ownership differs: '+key)
    require(len(record['elements']) == 700 and len({e['geometry'] for e in record['elements']}) == 700, 'expected retained 700 unique elements')
    g = record['geometry']
    for key, value in dict(left=-10., right=15., top=ref.top, d0=ref.d0, H=ref.H, h=ref.h).items():
        require(abs(g[key]-value) <= 64*np.finfo(float).eps*max(abs(value), ref.H), 'retained physical geometry differs: '+key)
    xaxis, raxis = mesh['x'], mesh['r']; xtol, rtol = mesh['x_tolerance'], mesh['r_tolerance']
    columns = {i: {} for i in range(mesh['nx'])}; classified = {}
    seen = {key: set() for key in ('u', 'p', 'k', 'epsilon')}; owners = {}; sizes = record['sizes']
    for e in record['elements']:
        xl, xr, rl, rr = e['bounds']
        require(all(math.isfinite(v) for v in e['bounds']) and xr > xl and rr > rl, 'invalid original physical bounds')
        left, right = axis_index(xaxis, xl, xtol, 'native left'), axis_index(xaxis, xr, xtol, 'native right')
        bottom, top = axis_index(raxis, rl, rtol, 'native bottom'), axis_index(raxis, rr, rtol, 'native top')
        require(right == left+1 and top == bottom+1 and bottom not in columns[left], 'native cell overlaps/skips XML intervals')
        require(e['geometry'] in mesh['cells'] and mesh['cells'][e['geometry']] == (left, bottom), 'native geometry ID/cell differs from exact XML owner')
        columns[left][bottom] = e; classified[e['geometry']] = (left, bottom)
        outer = bottom >= mesh['H_index']
        require(bool(e['epsilon']) == outer, 'epsilon is not exclusively outer')
        for key, size in zip(('u', 'p', 'k', 'epsilon'), sizes):
            ids, signs = e[key], e[po.SIGNS[key]]
            if key == 'epsilon' and not outer:
                require(not ids and not signs, 'inner epsilon map present'); continue
            count = 9 if key == 'p' else 16
            require(len(ids) == len(signs) == len(set(ids)) == count and all(type(v) is int and 0 <= v < size for v in ids) and
                    all(v in (-1., 1.) for v in signs), 'invalid native signed map: '+key)
            seen[key].update(ids)
    for column in columns.values():
        require(set(column) == set(range(mesh['nr'])), 'native column does not cover every XML radial interval')
        howner = column[mesh['H_index']]; require(howner['epsilon'], 'missing outer H owner')
        for e in column.values(): owners[e['geometry']] = howner
    for key, size in zip(('u', 'p', 'k', 'epsilon'), sizes):
        require(seen[key] == set(range(size)), 'incomplete global inventory: '+key)
    for ids, size in ((record['free_velocity_ids'][0], sizes[0]), (record['free_velocity_ids'][1], sizes[0]),
                      (record['k_free_ids'], sizes[2]), (record['epsilon_free_ids'], sizes[3])):
        require(len(ids) == len(set(ids)) and all(type(v) is int and 0 <= v < size for v in ids), 'invalid free mask')
    for field, size, free, sides in (('u', sizes[0], record['free_velocity_ids'][0], (0, 2)),
            ('u', sizes[0], record['free_velocity_ids'][1], (0, 2, 3)),
            ('k', sizes[2], record['k_free_ids'], (0,)), ('epsilon', sizes[3], record['epsilon_free_ids'], (0, 2))):
        fixed = set(); bottom = mesh['H_index'] if field == 'epsilon' else 0
        for e in record['elements']:
            if not e[field]: continue
            ix, ir = classified[e['geometry']]
            exterior = (ix == 0, ix+1 == mesh['nx'], ir == bottom, ir+1 == mesh['nr'])
            for side in sides:
                if exterior[side]: fixed.update(e[field][j] for j in po.EDGES[side])
        require(set(free) == set(range(size))-fixed, 'physical component/scalar essential mask differs')
    return owners


def exported(record, binary, reload=False):
    sizes = record['sizes']; values = dict(current=transient.state(record['coefficients'], sizes),
        previous=binary['previous'], advection=transient.state(record['current_advection'], sizes),
        previous_advection=transient.state(record['previous_advection'], sizes))
    scientific_time(record['step'], record['time'], DT, 'native endpoint')
    for key, header in zip(HISTORY_NAMES, HASH_NAMES):
        require(transient.state_hash(values[key]) == transient.state_hash(binary[key]) == record[header], 'binary/exported history binding differs: '+key)
    if reload:
        require('input_history' not in record, 'zero-step reload unexpectedly claims a final advancing input history')
    else:
        h = record['input_history']; require(h['step'] == record['step']-1 and h['dt'] == DT, 'pre-step clock/step label differs')
        scientific_time(h['step'], h['time'], DT, 'input history')
        require(transient.state_hash(transient.state(h['current'], sizes)) == record['previous_sha256'] and
                transient.state_hash(transient.state(h['advection'], sizes)) == record['previous_advection_sha256'], 'stored pre-step/current cache shift differs')
        values['older'] = transient.state(h['previous'], sizes)
        values['older_advection'] = transient.state(h['previous_advection'], sizes)
        if record['step'] == 1:
            require(transient.state_hash(values['previous']) == transient.state_hash(values['older']) and
                    transient.state_hash(values['previous_advection']) == transient.state_hash(values['older_advection']), 'BE initial duplicate history differs')
    for key in ('advection', 'previous_advection', 'older_advection'):
        if key in values: require(np.all(values[key]['flow'][2*sizes[0]:] == 0.), 'explicit history contains pressure loads')
    return values


def timings(record):
    seconds = record['seconds']; values = record['disjoint_seconds']
    require(math.isfinite(seconds) and seconds > 0 and set(values) == set(TIMERS), 'invalid complete step timing inventory')
    reserve = 128*np.finfo(float).eps*seconds
    require(all(math.isfinite(v) and v >= -reserve for v in values.values()) and
            abs(math.fsum(values.values())-seconds) <= reserve, 'disjoint timing does not close')
    return dict(seconds=seconds, disjoint=values, attributed_fraction=1-values['other']/seconds,
                warm_cost_gate='NOT_EVALUATED')


def accepted(record, control, injected):
    configuration(record, control)
    residual = record['residuals']
    require(len(residual) == 6 and all(math.isfinite(v) and 0 <= v <= 3e-12 for v in residual), 'accepted current residual exceeds unchanged native limit')
    for key in GUARDS: require(math.isfinite(record[key]) and record[key] > 0, 'nonpositive accepted guard: '+key)
    require(record['retries'] == int(injected), 'unexpected retry step/count')
    stats = record['coupled_newton']
    require(type(stats['counts_complete']) is bool and stats['counts_complete'] is (not injected), 'interrupted statistics completeness differs')
    for key in ('krylov_iterations', 'jacobian_actions', 'residual_evaluations', 'preconditioner_applications', 'rejected_trials'):
        require(type(stats[key]) is int and stats[key] >= 0, 'invalid Newton count: '+key)
    require(stats['residual_evaluations'] > 0, 'accepted step lacks live residual evaluation')
    if injected:
        require('qualification injected recoverable Krylov failure' in record['retry_reason'] and
                stats['preconditioner_applications'] > 0, 'retry lost reason/live callback work')
    else: require(not record['retry_reason'], 'unexpected retry reason')
    return dict(step=record['step'], time=record['time'], native_residuals=residual,
                guards={key: record[key] for key in GUARDS}, newton=stats, rebuilds=record['rebuilds'], timing=timings(record))


def cache_events(record, previous_event, cold, injected, cache_policy='OperatorDriftV1'):
    stats = record['coupled_newton']; reasons = ('cold', 'alpha', 'drift', 'age', 'control', 'retry', 'injected')
    fields.cache_policy_statistics(record, cache_policy)
    require(stats['preconditioner_policy'] == cache_policy and stats['relative_drift_limit'] == .02 and
            stats['maximum_uses'] == 8 and stats['preconditioner_events_recorded'] is True,
            'prospective inverse-cache policy/events differ')
    events = stats['preconditioner_events']; require(isinstance(events, list) and events, 'no recorded Newton preparation events')
    for key in ('refreshes', 'reuses', 'maximum_used_age'):
        require(type(stats[key]) is int and stats[key] >= 0, 'invalid preconditioner count: '+key)
    counts = dict.fromkeys(reasons, 0); reuses = 0; current = None if cold else previous_event
    require(cold or current is not None, 'missing previous prepared cache lineage')
    attempts = []; next_base = 0; last_attempt = 0
    for e in events:
        attempt, base = e['attempt'], e['base']
        require(type(attempt) is int and type(base) is int and attempt in (0, 1) and base >= 0, 'invalid preparation attempt/base')
        if attempt != last_attempt:
            require(injected and last_attempt == 0 and attempt == 1 and len(attempts) == 1,
                    'unexpected cache retry sequence')
            current = None; next_base = 0; last_attempt = 1
        require(base == next_base, 'preparation base indices skipped/repeated'); next_base += 1; attempts.append(attempt)
        alpha = (1. if record['step'] == 1 else 1.5)/DT
        require(e['alpha'] == alpha and math.isfinite(e['drift']) and e['drift'] >= 0,
                'cache uses wrong BE/BDF2 alpha or invalid drift')
        require(type(e['age_before']) is int and type(e['age_used']) is int and type(e['refreshed']) is bool,
                'invalid prepared cache age/refresh types')
        before = 0 if current is None else current['age_used']
        require(e['age_before'] == before, 'cache age lineage differs across bases/steps/restart')
        if injected and attempt == 0: expected = 'injected'
        elif attempt == 1: expected = 'retry'
        elif current is None: expected = 'cold'
        elif current['alpha'] != alpha: expected = 'alpha'
        elif e['drift'] > .02: expected = 'drift'
        elif before >= 8: expected = 'age'
        else: expected = 'reuse'
        require(e['reason'] == expected and e['refreshed'] is (expected != 'reuse'), 'refresh reason/precedence differs from declared policy')
        require(e['age_used'] == (before+1 if expected == 'reuse' else 1) and 1 <= e['age_used'] <= 8,
                'prepared cache use count is not one per Newton base')
        if current is None: require(e['drift'] == 0., 'cold cache has a spurious saved-cache drift')
        if expected == 'reuse': reuses += 1
        else: counts[expected] += 1
        current = e
    require(attempts == ([0]+[1]*(len(events)-1) if injected else [0]*len(events)) and
            (not injected or len(events) >= 2), 'retry did not perform one injected refresh then cold all-Newton refreshes')
    require(stats['refresh_reasons'] == counts and stats['refreshes'] == sum(counts.values()) and stats['reuses'] == reuses and
            stats['maximum_used_age'] == max(e['age_used'] for e in events), 'cache summary counters disagree with events')
    require(record['rebuilds'] == 3*stats['refreshes'], 'every declared refresh must perform exactly three block-factor rebuilds')
    if record['step'] == 2:
        require(not cold and events[0]['reason'] == 'alpha', 'BE to BDF2 did not refresh the changed alpha')
    return dict(policy=cache_policy, refreshes=stats['refreshes'], reuses=reuses, refresh_reasons=counts,
                maximum_used_age=stats['maximum_used_age'], events=events,
                drift_measure=stats.get('drift_measure', 'CurrentInverseOperatorFrobeniusV1'),
                drift_scope='Reported current-versus-prepared metric validated for policy consistency; not independently reconstructed here.'), current


def compare(a, b, data, owners, ref, deadline):
    result = dict(coefficients={}, owner_fields={}, dependent={}, exact_hashes={})
    common = [key for key in a if key in b]
    for key in common:
        clock(deadline); delta = fields.coefficient_delta(a[key], b[key], data['sizes'], ref)
        require(all(v['scaled_MAX'] <= LIMIT for v in delta.values()), 'history coefficient difference exceeds 1e-8: '+key)
        result['coefficients'][key] = delta
        result['exact_hashes'][key] = dict(first=transient.state_hash(a[key]), second=transient.state_hash(b[key]),
                                          equal=transient.state_hash(a[key]) == transient.state_hash(b[key]))
    for key in ('current', 'previous', 'older'):
        if key not in common: continue
        result['owner_fields'][key] = [fields.physical_delta(a[key], b[key], data, ref, q, grid, deadline)
                                      for q, grid in ((8, 9), (12, 17))]
        require(all(value.get('status', 'MET') == 'MET' for row in result['owner_fields'][key] for value in row['norms'].values()), 'owner field difference exceeds 1e-8')
        result['dependent'][key] = fields.dependent_delta(a[key], b[key], data, owners, ref, deadline)
        require(result['dependent'][key]['gate_pass'], 'current K/Gamma/cutoff/Ew differences exceed 1e-8')
    result['unavailable_history_roles'] = sorted(set(a)-set(b))
    return result


def equation_report(spec, root, native_path, data, histories, table_sha):
    path, evidence = input_path(spec, root); report = read(path)
    require(report['schema'] == transient.SCHEMA and report['native']['sha256'] == hr.identity(native_path)['sha256'] and
            report['table']['sha256'] == table_sha and report['source_identity'] == data['source_identity'] and
            report['state_sha256'] == data['state_sha256'], 'independent equation report identity differs')
    for key in ('accepted_time_discrete_equations_status', 'independent_quadrature_status', 'accepted_history_boundary_status'):
        require(report[key] == 'MET', 'independent current-equation subgate failed: '+key)
    core = report['core']
    require(core['transient_residual_status'] == 'MET' and core['history_cache_checks']['status'] == 'MET' and
            core['accounting_identity_status'] == 'MET' and all(v['status'] == 'MET' for v in report['boundary_checks'].values()),
            'independent residual/cache/accounting/boundary subgate failed')
    require(all(math.isfinite(core['residuals'][key]) and 0 <= core['residuals'][key] <= limit for key, limit in po.RESIDUAL_LIMITS.items()),
            'independent residual numbers exceed frozen limits')
    scheme = report['scheme']; require(scheme['dt'] == DT and scheme['startup'] is False and scheme['alpha'] == 1.5/DT and
        scheme['mass_coefficients'] == [1.5, -2., .5] and scheme['explicit_coefficients'] == [2., -1.], 'independent final scheme differs')
    old = report['input_history_identity']
    require(old['step'] == data['step']-1 and abs(old['time']-(data['step']-1)*DT) <= 1e-13, 'independent history clock differs')
    for header, key in (('current_sha256', 'previous'), ('previous_sha256', 'older'),
                        ('advection_sha256', 'previous_advection'), ('previous_advection_sha256', 'older_advection')):
        require(old[header] == transient.state_hash(histories[key]), 'independent history/cache binding differs: '+key)
    return dict(file=evidence, overall_status=report['status'], current_equations='MET',
                residuals=core['residuals'], overall_NOT_MET_allowed=True,
                scope='Full supplied transient oracle verifies these accepted equations; other physical/mesh/long-time gates remain open.')


def zero_terminal_export(packet, path, control, root):
    """Imported history serialization is charged to export, not advancement."""
    prefix = replace_suffix(path, '.result.json', '')
    index_path = Path(str(prefix)+'.science-index.json')
    require(hr.resolve(root, packet['scientific_index']) == index_path, 'zero-step scientific index prefix differs')
    index = read(index_path); record = packet['state']
    require(index['schema'] == 'WJH27NativeScientificIndexV1' and index['run_complete'] is True and
            index['terminal_output'] is True and index['coefficients_included'] is True and
            index['terminal_begin'] == .0173 and index['terminal_end'] == .0537,
            'zero-step terminal index configuration differs')
    for key in ('source_identity', 'reference_identity', 'signature'):
        require(index[key] == record[key], 'zero-step index identity differs')
    require(index['initial_step'] == index['final_step'] == 12 and
            index['initial_time'] == index['final_time'] == 12*DT and index['dt'] == DT,
            'zero-step index claims a different scientific clock')
    require(len(index['frames']) == 1, 'zero-step terminal export must own exactly one final frame')
    frame = index['frames'][0]; frame_path = Path(str(prefix)+'.science-12.json'); actual = read(frame_path)
    require(hr.resolve(root, frame['path']) == frame_path and frame['step'] == 12 and frame['time'] == 12*DT and
            frame['origin'] == 'imported_accepted_history' and frame['reason_mask'] == 128,
            'zero-step frame is falsely marked as an advancing observation')
    require(hr.identity(frame_path)['sha256'] == frame['sha256'] and frame_path.stat().st_size == frame['bytes'],
            'zero-step scientific frame byte identity differs')
    configuration(actual, control)
    require(actual['step'] == 12 and actual['time'] == 12*DT and 'input_history' not in actual and
            all(actual[key] == record[key] for key in HASH_NAMES), 'zero-step exported frame history differs')
    require(transient.state_hash(transient.state(actual['coefficients'], actual['sizes'])) == actual['state_sha256'],
            'zero-step terminal current coefficients differ')
    for key, header in (('current_advection', 'advection_sha256'), ('previous_advection', 'previous_advection_sha256')):
        require(transient.state_hash(transient.state(actual[key], actual['sizes'])) == actual[header],
                'zero-step terminal explicit history differs: '+key)
    require({p.resolve() for p in prefix.parent.glob(prefix.name+'.science-[0-9]*.json')} == {frame_path},
            'zero-step terminal scientific retention differs')
    frame_seconds, index_seconds, export = packet['scientific_frame_seconds'], packet['scientific_index_seconds'], packet['final_export_seconds']
    require(all(math.isfinite(v) and v > 0 for v in (frame_seconds, index_seconds, export)) and
            frame_seconds == frame['write_seconds'] and frame_seconds+index_seconds <= export+128*np.finfo(float).eps*export,
            'zero-step terminal frame/index work is outside final export')
    return dict(index=hr.identity(index_path), frame=hr.identity(frame_path), provenance=frame['origin'],
                scientific_frame_seconds=frame_seconds, scientific_index_seconds=index_seconds,
                final_export_seconds=export, contained_in_final_export=True,
                advanced_residual_evidence=False, full_terminal_window_assessment=False)


def run(manifest, folder, table_path, deadline, out):
    root = hr.resolve(folder, manifest.get('root', '.')); runs = manifest['runs']
    needed = {'control', 'prefix', 'cold', 'recovery9', 'reload12'}
    require(needed <= set(runs) <= needed|{'retry'}, 'wrong physical history run inventory')
    require(set(manifest['sessions']) == set(runs), 'every native operation needs its exact session XML')
    paths = {name: hr.resolve(root, value) for name, value in runs.items()}
    packets = {name: read(path) for name, path in paths.items()}; control = packets['control']['state']
    cache_policy = manifest.get('preconditioner_policy', 'OperatorDriftV1')
    require(cache_policy in ('OperatorDriftV1', 'CurrentDataDriftV1') and
            control['coupled_newton']['preconditioner_policy'] == cache_policy, 'manifest/cache policy differs')
    terminal_zero = manifest.get('terminal_zero_export', False)
    require(type(terminal_zero) is bool, 'terminal-zero-export flag must be boolean')
    require(hr.identity(table_path)['sha256'] == po.TABLE_SHA == control['reference_identity'], 'table is not the qualified active v2 table')
    ref = po.Reference(read(table_path)); mesh = session_geometry(manifest['sessions']['control'], root, ref)
    owners = inventory(control, control, ref, mesh)
    if 'source_identity' in manifest: require(manifest['source_identity'] == control['source_identity'], 'manifest source differs')
    out.update(root=str(root), source_identity=control['source_identity'], signature=control['signature'],
               preconditioner_policy=cache_policy, terminal_zero_export_requested=terminal_zero, runs={}, comparisons={},
               input_checkpoints={}, equation_reports={}, scales=dict(U=1., Pi=ref.tau, k=ref.Kin, epsilon_coefficient=ref.A/ref.H,
               epsilon_field='Ain/r', cutoff=ref.H), original_owner_bounds_preserved=True,
               authoritative_session_geometry=mesh)
    inputs = {}; input_states = {}
    require(set(manifest['input_checkpoints']) == {'cold', 'recovery9', 'reload12'}, 'missing pre-launch checkpoint inventory')
    for name, spec in manifest['input_checkpoints'].items():
        require(isinstance(spec, dict) and set(spec) >= {'path', 'sha256'}, 'input checkpoint needs pre-launch SHA256')
        path, evidence = input_path(spec, root)
        states, binary = binary_history(path, control, DT, deadline, evidence['sha256'])
        require(binary['step'] == STARTS[name], 'input checkpoint wrong step')
        inputs[name] = binary; input_states[name] = states; out['input_checkpoints'][name] = binary
    values = {}; logs = {}; final_states = {}
    for name, packet in packets.items():
        clock(deadline); path = paths[name]; reload = name == 'reload12'; end = 3 if name == 'prefix' else FINAL_STEP
        require(packet['schema'] == 'WJH27RegisteredCoupledRunV2' and packet['factory'] == 'WJNativeCoupledV2' and
                packet['normal_main'] is True and packet['steady_mode'] is False, 'not a completed registered advancing/import operation')
        require(packet['restart'] is (name in inputs), 'restart flag differs from declared role')
        record = packet['state']; configuration(record, control)
        run_mesh = session_geometry(manifest['sessions'][name], root, ref)
        require(run_mesh['geometry_sha256'] == mesh['geometry_sha256'], 'session geometry/owners differ from authoritative control XML')
        inventory(record, control, ref, mesh)
        session_evidence = session(manifest['sessions'][name], root, name, path, inputs.get(name), po.TABLE_SHA,
                                   cache_policy, terminal_zero and reload)
        require(session_evidence['file'] == run_mesh['session'], 'session XML changed during audit')
        session_evidence['geometry_sha256'] = run_mesh['geometry_sha256']
        require(record['step'] == end, 'wrong endpoint step'); scientific_time(end, record['time'], DT, name)
        require(len(packet['complete_step_seconds']) == end-STARTS[name], 'executed-step inventory differs')
        require(all(math.isfinite(value) and value > 0 for value in packet['complete_step_seconds']), 'invalid complete step timing')
        complete = packet['complete_advancing_seconds']; total = math.fsum(packet['complete_step_seconds'])
        require(math.isfinite(complete) and complete >= 0 and abs(complete-total) <= 128*np.finfo(float).eps*total,
                'complete advancing time does not sum to actual step calls')
        final_path = replace_suffix(path, '.result.json', '.final.bin')
        binary, evidence = binary_history(final_path, control, DT, deadline)
        require(evidence['step'] == end, 'final binary clock differs')
        values[name] = exported(record, binary, reload); final_states[name] = binary
        log_path = replace_suffix(path, '.result.json', '.steps.jsonl')
        rows = [hr.native_json(line) for line in log_path.read_text().splitlines() if line.strip()]
        require([row['step'] for row in rows] == list(range(STARTS[name]+1, end+1)), 'missing/extra accepted step rows')
        previous = None; previous_event = None; diagnostics = []
        for index, row in enumerate(rows):
            clock(deadline); injected = name == 'retry' and row['step'] == 4
            scientific_time(row['step'], row['time'], DT, 'step log')
            diagnosis = accepted(row, control, injected)
            cache, previous_event = cache_events(row, previous_event, index == 0, injected, cache_policy)
            diagnosis['cache'] = cache; diagnostics.append(diagnosis)
            if previous is not None:
                require(row['previous_sha256'] == previous['state_sha256'] and row['previous_advection_sha256'] == previous['advection_sha256'], 'accepted history/cache chain differs')
            elif name in inputs:
                require(row['previous_sha256'] == inputs[name]['hashes']['current'] and row['previous_advection_sha256'] == inputs[name]['hashes']['advection'], 'cold continuation did not consume declared checkpoint')
            previous = row
        logs[name] = rows
        if reload:
            require(not packet['complete_step_seconds'] and not rows and packet['complete_advancing_seconds'] == 0., 'zero-step reload performed advancement')
            require(evidence['sha256'] == inputs[name]['sha256'], 'zero-step reload changed complete binary history')
            cache = record['coupled_newton']
            require(cache['preconditioner_events'] == [] and cache['refreshes'] == cache['reuses'] == record['rebuilds'] == 0,
                    'zero-step reload claims advancing inverse preparations')
            require(cache['maximum_used_age'] == 0 and all(v == 0 for v in cache['refresh_reasons'].values()),
                    'zero-step reload claims advancing cache ages/reasons')
            if terminal_zero: out['zero_step_terminal_export'] = zero_terminal_export(packet, path, control, root)
        else:
            accepted(record, control, False)
            require(all(record[key] == rows[-1][key] for key in HASH_NAMES), 'final native/log identity differs')
            require(record['coupled_newton'] == rows[-1]['coupled_newton'] and record['rebuilds'] == rows[-1]['rebuilds'],
                    'final/log inverse-cache diagnostics differ')
        out['runs'][name] = dict(result=hr.identity(path), session=session_evidence, log=hr.identity(log_path), final_checkpoint=evidence,
            executed_steps=end-STARTS[name], accepted_endpoint_step=end, accepted_endpoint_time=end*DT,
            history_hashes={key: transient.state_hash(v) for key, v in values[name].items()}, steps=diagnostics,
            advanced_step_residual_evidence=not reload,
            reload_note='Zero executed steps: native default residual/statistic zeros are not convergence evidence.' if reload else None)
    refreshes = sum(row['coupled_newton']['refreshes'] for row in logs['control'])
    reuses = sum(row['coupled_newton']['reuses'] for row in logs['control'])
    require(refreshes > 0 and reuses > 0, 'control did not exercise both actual cache refresh and reuse')
    out['cache_exercised'] = dict(control_refreshes=refreshes, control_reuses=reuses,
        BE_to_BDF2_alpha_refresh=True, cold_restarts_start_unprepared=True,
        observed_control_refresh_reasons=sorted({event['reason'] for row in logs['control']
            for event in row['coupled_newton']['preconditioner_events'] if event['refreshed']}),
        cold_all_Newton_retry_exercised='retry' in runs)
    require(inputs['cold']['path'] == str(replace_suffix(paths['prefix'], '.result.json', '.final.bin')), 'cold input is not prefix final file')
    for key in HISTORY_NAMES: require(transient.state_hash(input_states['cold'][key]) == transient.state_hash(final_states['prefix'][key]), 'cold input differs from full prefix history')
    for name in ('cold', 'recovery9', 'reload12', 'retry'):
        if name not in values: continue
        out['comparisons'][name] = compare(values['control'], values[name], control, owners, ref, deadline)
        lifts = {}
        for key, scale in (('essential_flow_lift', 1.), ('essential_k_lift', ref.Kin), ('essential_epsilon_lift', ref.A/ref.H)):
            a, b = np.asarray(control[key]), np.asarray(packets[name]['state'][key])
            require(a.shape == b.shape and a.ndim == 1 and np.all(np.isfinite(a)) and np.all(np.isfinite(b)), 'essential lift inventory differs')
            error = float(np.max(abs(a-b), initial=0.))/scale
            require(error <= LIMIT, 'essential lift disagreement exceeds 1e-8'); lifts[key] = error
        out['comparisons'][name]['essential_lifts'] = lifts
    prefix = replace_suffix(paths['control'], '.result.json', '')
    retained = {Path(str(prefix)+'.recovery-'+str(step)+'.bin') for step in (9, 12)}
    removed = [Path(str(prefix)+'.recovery-'+str(step)+'.bin') for step in (3, 6)]
    require(all(not p.exists() for p in removed), 'recovery3/6 not actually rotated')
    require({p.resolve() for p in prefix.parent.glob(prefix.name+'.recovery-*.bin')} == retained ==
            {hr.resolve(root, p) for p in packets['control']['retained_recovery']}, 'actual/reported recovery9/12 inventory differs')
    require(len(packets['control']['recovery_seconds']) == 4 and packets['control']['rotated_recovery_count'] == 2,
            'control did not report four writes/two rotations')
    require(not list(prefix.parent.glob(prefix.name+'.recovery-*.bin.tmp')), 'unfinished recovery temporary')
    science_paths = {Path(str(prefix)+'.science-'+str(step)+'.json') for step in (3, 6, 9, 12)}
    index_path = Path(str(prefix)+'.science-index.json'); has_index = 'scientific_index' in packets['control']
    require(index_path.exists() is has_index, 'scientific descriptor presence differs from native result')
    expected_paths = science_paths|({index_path} if has_index else set())
    require({p.resolve() for p in prefix.parent.glob(prefix.name+'.science-*.json')} == expected_paths, 'science retention differs')
    out['rotation'] = dict(removed_absent=[str(p) for p in removed], retained=[hr.identity(p) for p in sorted(retained)], science=[])
    index_frames = {}
    if has_index:
        require(hr.resolve(root, packets['control']['scientific_index']) == index_path, 'scientific descriptor prefix differs')
        descriptor = read(index_path)
        require(descriptor['schema'] == 'WJH27NativeScientificIndexV1' and descriptor['run_complete'] is True and
                descriptor['terminal_output'] is False and descriptor['coefficients_included'] is True,
                'unexpected control scientific descriptor policy')
        require(all(descriptor[key] == control[key] for key in ('source_identity', 'reference_identity', 'signature', 'dt')) and
                descriptor['initial_step'] == 0 and descriptor['initial_time'] == 0. and descriptor['final_step'] == 12 and
                descriptor['final_time'] == 12*DT, 'scientific descriptor source/clock differs')
        require([f['step'] for f in descriptor['frames']] == [3, 6, 9, 12], 'scientific descriptor numbered-frame inventory differs')
        index_frames = {f['step']: f for f in descriptor['frames']}
        total = math.fsum(f['write_seconds'] for f in descriptor['frames'])
        require(all(math.isfinite(f['write_seconds']) and f['write_seconds'] > 0 for f in descriptor['frames']) and
                abs(total-packets['control']['scientific_frame_seconds']) <= 128*np.finfo(float).eps*total and
                0 < packets['control']['scientific_index_seconds'] <= packets['control']['final_export_seconds'],
                'scientific descriptor/frame timing differs')
        out['rotation']['scientific_index'] = hr.identity(index_path)
    for step in (3, 6, 9, 12):
        clock(deadline); path = Path(str(prefix)+'.science-'+str(step)+'.json'); record = read(path)
        if has_index:
            frame = index_frames[step]
            require(hr.resolve(root, frame['path']) == path and hr.identity(path)['sha256'] == frame['sha256'] and
                    path.stat().st_size == frame['bytes'] and frame['time'] == step*DT and
                    frame['reason_mask'] == 1 and frame['origin'] == 'current_process_accepted_step',
                    'numbered scientific frame/descriptor binding differs')
        configuration(record, control); inventory(record, control, ref, mesh); accepted(record, control, False)
        require(record['step'] == step and all(record[key] == logs['control'][step-1][key] for key in HASH_NAMES), 'science/log binding differs')
        require(record['coupled_newton'] == logs['control'][step-1]['coupled_newton'], 'science/log cache event binding differs')
        h = record['input_history']; pseudo = dict(current=transient.state(record['coefficients'], record['sizes']),
            previous=transient.state(h['current'], record['sizes']), advection=transient.state(record['current_advection'], record['sizes']),
            previous_advection=transient.state(record['previous_advection'], record['sizes']))
        scientific_values = exported(record, pseudo)
        item = dict(file=hr.identity(path), hashes={key: transient.state_hash(v) for key, v in scientific_values.items()})
        if step == 3:
            out['comparisons']['prefix_at_step3'] = compare(scientific_values, values['prefix'], control, owners, ref, deadline)
        if step in (9, 12):
            name = 'recovery9' if step == 9 else 'reload12'
            require(inputs[name]['path'] == str(Path(str(prefix)+'.recovery-'+str(step)+'.bin')), 'restart did not import actual retained recovery file')
            for key in HISTORY_NAMES: require(transient.state_hash(input_states[name][key]) == transient.state_hash(scientific_values[key]), 'retained recovery/science history differs')
            item['checkpoint'] = inputs[name]
        out['rotation']['science'].append(item)
    required_reports = {'control', 'cold', 'recovery9'}|({'retry'} if 'retry' in values else set())
    require(set(manifest['current_equation_reports']) == required_reports, 'full transient reports missing for advancing comparisons')
    for name, spec in manifest['current_equation_reports'].items():
        clock(deadline); out['equation_reports'][name] = equation_report(spec, root, paths[name], packets[name]['state'], values[name], po.TABLE_SHA)
    out.update(status='PASS_PHYSICAL_CACHE_HISTORY_FOUNDATION_ONLY', foundation_gate_pass=True,
               retry_exercised='retry' in values, cache_equations_verified_by_supplied_independent_reports=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('manifest', 'table', 'out'): parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--seconds', type=float, required=True); args = parser.parse_args(); began = time.monotonic()
    if args.out.exists(): parser.error('immutable output already exists')
    out = dict(schema='WJH27PhysicalCacheHistoryCheckV1', status='FAIL', foundation_gate_pass=False,
        full_H27_qualification=False, H27_D_pass=False, physical_long_restart_qualified=False,
        temporal_order_qualified=False, terminal_handling_qualified=False, long_pilots_qualified=False,
        scope='Retained 700-cell active-v2 short physical history/cache FOUNDATION only.',
        sources={key: hr.identity(path) for key, path in dict(checker=__file__, JSON_reader=hr.__file__,
            field_comparison=fields.__file__, reference_and_basis=po.__file__, Jacobi_basis=po.fixtures.__file__,
            history_arrays=transient.__file__).items()}, manifest=hr.identity(args.manifest), table=hr.identity(args.table),
        prospective=dict(dt=DT, final_step=FINAL_STEP, final_time=.06, comparison_limit=LIMIT,
            native_residual_limit=3e-12, retained_recovery=[9, 12], rotated_recovery=[3, 6]), workers=1)
    try:
        require(sys.platform.startswith('linux') and sys.byteorder == 'little', 'authoritative little-endian Linux analysis required')
        require(math.isfinite(args.seconds) and 0 < args.seconds <= 300, 'explicit analysis cap must be in (0,300] seconds')
        soft, hard = resource.getrlimit(resource.RLIMIT_AS)
        cap = min([8*1024**3]+[v for v in (soft, hard) if v != resource.RLIM_INFINITY])
        resource.setrlimit(resource.RLIMIT_AS, (cap, hard)); out['memory_cap_bytes'] = cap
        manifest = read(args.manifest); require(manifest['schema'] == 'WJH27PhysicalCacheHistoryV1', 'wrong manifest schema')
        run(manifest, args.manifest.resolve().parent, args.table.resolve(), began+args.seconds, out)
    except Exception as error: out.update(exception=type(error).__name__, reason=str(error))
    out['elapsed_seconds'] = time.monotonic()-began
    out['maximum_resident_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x') as stream: json.dump(out, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps({key: out[key] for key in ('status', 'reason', 'elapsed_seconds') if key in out}))
    if not out['foundation_gate_pass']: raise SystemExit(2)


if __name__ == '__main__': main()
