#!/usr/bin/env python3
"""Read-only finite native-state coverage of M8 constitutive certificates.

The certificate bounds a declared continuum constitutive box.  This audit
asks whether actual accepted native evaluations fall inside it; it never
turns finite point coverage into a polynomial-everywhere bound, invariant
state set, outer-turbulence positivity proof or terrain qualification.
"""
import argparse
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import platform
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
THREADS = {name: '1' for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                                'MKL_NUM_THREADS', 'BLIS_NUM_THREADS')}
os.environ.update(THREADS)
import numpy as np

from validate_cls_epawf_m4 import Space
from validate_cls_epawf_m5_geometry import geometry


def sha(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            result.update(chunk)
    return result.hexdigest()


def number(value):
    return float(Fraction(str(value)))


def reconstruction_rules(raw, session):
    """Read actual native GLL counts from the source-bound input XML.

    Current reports carry degrees/maps but do not carry each expansion's
    native GLL point count.  Missing or ambiguous input makes native-GLL
    coverage inconclusive; this audit does not guess a default.
    """
    if session is None:
        return {}, 'Native GLL rule unavailable: no matching session supplied'
    root = ET.parse(session).getroot()
    entries = []
    for entry in root.findall('EXPANSIONS/E'):
        variables = {value.strip() for value in entry.get('FIELDS', '').split(',')}
        if 'u' not in variables:
            continue
        counts = [int(value.strip()) for value in entry.get('NUMPOINTS', '').split(',') if value.strip()]
        modes = [int(value.strip()) for value in entry.get('NUMMODES', '').split(',') if value.strip()]
        kinds = [value.strip() for value in entry.get('POINTSTYPE', '').split(',')]
        if len(counts) == 1:
            counts *= raw['dimension']
        if len(modes) == 1:
            modes *= raw['dimension']
        if len(kinds) == 1:
            kinds *= raw['dimension']
        if len(counts) != raw['dimension'] or len(modes) != raw['dimension'] or not all(
                value == 'GaussLobattoLegendre' for value in kinds):
            return {}, 'Native velocity expansion is not an explicit tensor GLL rule'
        entries.append((modes, counts))
    result = {}
    for element in raw['fields'][0]['mesh']:
        matches = {tuple(counts) for modes, counts in entries if all(
                   value == element['degree']+1 for value in modes)}
        if len(matches) != 1:
            return {}, 'Ambiguous or absent native GLL rule for geometry '+str(element['geometry'])
        result[element['geometry']] = list(next(iter(matches)))
    return result, None


def tensor_points(element, axes, nodes):
    result = []
    for indices in itertools.product(*(range(len(value)) for value in nodes)):
        point = [0., 0., 0.]
        for axis, values, index in zip(axes, nodes, indices):
            point[axis] = element['lower'][axis]+(values[index]+1)*(element['upper'][axis]-element['lower'][axis])/2
        result.append(point)
    return np.asarray(result)


def gll(count):
    if count < 2:
        raise ValueError('Invalid GLL count')
    return np.r_[-1., np.polynomial.legendre.Legendre.basis(count-1).deriv().roots(), 1.]


def point_sets(raw, fields, native_rules):
    """Yield actual outer evaluation sets; repeated traces retain their side."""
    d, H = raw['dimension'], raw['parameters'][1]
    source_order, normal_order = map(int, raw['quadrature_orders'])
    axes = fields[0].axes
    outer = [element for element in fields[0].elements if element['lower'][2] >= H-1e-12]
    for element in outer:
        q = np.polynomial.legendre.leggauss(source_order)[0]
        yield 'source_volume_outer', element['geometry'], tensor_points(element, axes, [q]*d)
        mass_order = max(normal_order, source_order) if raw['mapped'] else max(8, element['degree']+3)
        q = np.polynomial.legendre.leggauss(mass_order)[0]
        yield 'mass_volume_outer', element['geometry'], tensor_points(element, axes, [q]*d)
        if element['geometry'] in native_rules:
            yield 'reconstruction_native_outer', element['geometry'], tensor_points(
                element, axes, [gll(count) for count in native_rules[element['geometry']]])
        q = np.linspace(-1., 1., 2*element['degree']+3)
        yield 'reconstruction_dense_outer', element['geometry'], tensor_points(element, axes, [q]*d)
    horizontal_axes = (0,) if d == 2 else (0, 1)
    for element in fields[d+1].elements:
        if abs(element['lower'][2]-H) > 1e-12:
            continue
        for name, nodes in [('interface_outer', np.polynomial.legendre.leggauss(8)[0]),
                            ('trace_transfer_outer', gll(max(8, element['degree']+3)))]:
            points = []
            for indices in itertools.product(range(len(nodes)), repeat=d-1):
                point = [0., 0., H]
                for axis, index in zip(horizontal_axes, indices):
                    point[axis] = element['lower'][axis]+(nodes[index]+1)*(element['upper'][axis]-element['lower'][axis])/2
                points.append(point)
            points = np.asarray(points)
            yield name, element['geometry'], points
            if raw['upstream_transition'] and name == 'trace_transfer_outer':
                for index, distance in enumerate(raw['transition_controls'][2:4]):
                    moved = points.copy()
                    moved[:, 2] = distance/geometry(raw['map_parameters'], points)['scale']
                    # A physical sensor can lie in another normal element;
                    # grouping by actual outer geometry preserves its trace.
                    groups = {}
                    for point in moved:
                        selected = fields[d+1].element(point)
                        groups.setdefault(selected['geometry'], []).append(point)
                    for identity, values in groups.items():
                        yield 'trace_transfer_sensor_'+str(index), identity, np.asarray(values)
        if raw['upstream_transition']:
            nodes = np.polynomial.legendre.leggauss(source_order)[0]
            points = []
            for indices in itertools.product(range(source_order), repeat=d-1):
                point = [0., 0., H]
                for axis, index in zip(horizontal_axes, indices):
                    point[axis] = element['lower'][axis]+(nodes[index]+1)*(element['upper'][axis]-element['lower'][axis])/2
                points.append(point)
            points = np.asarray(points)
            for index, distance in enumerate(raw['transition_controls'][2:4]):
                moved = points.copy()
                moved[:, 2] = distance/geometry(raw['map_parameters'], points)['scale']
                groups = {}
                for point in moved:
                    identity = fields[d+1].element(point)['geometry']
                    groups.setdefault(identity, []).append(point)
                for identity, values in groups.items():
                    yield 'physical_outer_sensor_'+str(index), identity, np.asarray(values)


def local_values(raw, fields, state, identity, points):
    """Independent Cartesian gradient, including C1-foot one-sided metrics."""
    d, offsets = raw['dimension'], raw['offsets']
    element = fields[0].by_id[identity]
    indices, phi, chart_rows = fields[0].rows(element, points)
    foot_side = np.zeros(len(points), dtype=int)
    if int(raw['map_parameters'][0]) == 4:
        length = raw['map_parameters'][12]
        midpoint = (element['lower'][0]+element['upper'][0])/2
        for foot in (-length, length):
            at = np.abs(points[:, 0]-foot) <= 32*np.finfo(float).eps*np.maximum(
                np.maximum(1., abs(foot)), np.abs(points[:, 0]))
            foot_side[at] = -1 if midpoint < foot else 1
    metric = geometry(raw['map_parameters'], points, foot_side)
    rows = np.einsum('qbi,qba->qai', chart_rows, metric['inverse'])
    gradient = np.zeros((len(points), 3, 3))
    for block, axis in enumerate((0, 2) if d == 2 else (0, 1, 2)):
        gradient[:, axis] = np.einsum('qai,i->qa', rows, state[offsets[block]+indices])
    values = []
    for block in (d, d+1):
        element = fields[block].by_id[identity]
        ids, basis, _ = fields[block].rows(element, points)
        values.append(basis@state[offsets[block]+ids])
    return gradient, values[0], values[1], metric['distance']


def verify(raw, certificate, session=None, maximum_examples=20):
    if certificate.get('status') == 'failed':
        raise ValueError('A failed certificate cannot support coverage qualification')
    d = raw['dimension']
    fields = [Space(record, d) for record in raw['fields']]
    native_rules, native_issue = reconstruction_rules(raw, session)
    flow, normal = np.asarray(raw['robust_reference_flow']), np.asarray(raw['robust_reference_normal'])
    if not (np.all(np.isfinite(flow)) and np.all(np.isfinite(normal)) and
            abs(flow@normal) < 1e-12 and abs(flow@flow-1) < 1e-12 and abs(normal@normal-1) < 1e-12):
        raise ValueError('Invalid incoming orthonormal reference frame')
    frame = np.column_stack([flow, np.cross(normal, flow), normal])
    beta, uniform, _, _, _, scale = raw['robust_parameters'][:6]
    boxes = []
    for index, box in enumerate(certificate.get('finite_boxes', [])):
        if number(box['beta']) != beta:
            continue
        center = np.asarray(box['gradient_center'], float)
        radius = number(box['gradient_entry_radius'])
        krange = list(map(number, box['k_over_k0']))
        erange = list(map(number, box['epsilon_over_epsilon0']))
        candidates = []
        for candidate in box.get('fixed_gamma_candidates', []):
            Gamma = number(candidate['gamma_over_nut0'])
            candidates.append(dict(gamma_over_nut0=Gamma,
                certificate_joint_gate=bool(candidate['joint_constitutive_gate']),
                certificate_momentum_target_met=bool(candidate['momentum_target_met']),
                certificate_covariance_positive=bool(candidate['certified_positive']),
                actual_gamma_exact_match=uniform == 0 and scale == Gamma,
                note='Covariance proof is for this fixed Gamma; a larger Gamma is not automatically covered'))
        boxes.append(dict(certificate_box_index=index, beta=beta, center=center.tolist(),
            radius=radius, k_range=krange, epsilon_range=erange,
            required_gamma_over_nut0=number(box['required_gamma_over_nut0']),
            target_over_nut0=number(box['target_over_nut0']),
            fixed_gamma_candidates=candidates, sample_count=0, violation_count=0,
            violations_by_set={}, examples=[]))
    if not boxes:
        raise ValueError('No finite certificate box for the actual common cubic factor')
    sets = list(point_sets(raw, fields, native_rules))
    set_counts = {}
    extrema = dict(gradient_entry_deviation=0., k_over_K0=[math.inf, -math.inf],
                   epsilon_over_epsilon0=[math.inf, -math.inf], gamma_over_nut0=[math.inf, -math.inf])
    states = raw.get('accepted_states', [])
    times = raw.get('audit_times', [])
    if not states or len(states) != len(times):
        raise ValueError('Missing accepted native state/time history')
    us, qstar, kappa = raw['ustar'], raw['parameters'][6], raw['parameters'][7]
    for state_index, (coefficients, epoch) in enumerate(zip(states, times)):
        state = np.asarray(coefficients, float)
        for name, identity, points in sets:
            gradient, k, epsilon, r = local_values(raw, fields, state, identity, points)
            if not np.all(np.isfinite(gradient)) or not np.all(np.isfinite(k+epsilon+r)) or np.any(r <= 0):
                raise ValueError('Invalid outer coefficient reconstruction')
            G0 = us/(kappa*r)
            normalized = np.einsum('ai,qab,bj->qij', frame, gradient, frame)/G0[:, None, None]
            k_ratio = k/(qstar*us*us)
            epsilon_ratio = epsilon/(us**3/(kappa*r))
            actual_gamma = scale+uniform/(kappa*us*r)
            count = len(points)
            set_counts[name] = set_counts.get(name, 0)+count
            center = np.zeros((3, 3)); center[0, 2] = 1.
            extrema['gradient_entry_deviation'] = max(extrema['gradient_entry_deviation'], float(np.max(np.abs(normalized-center))))
            for key, values in [('k_over_K0', k_ratio), ('epsilon_over_epsilon0', epsilon_ratio),
                                ('gamma_over_nut0', actual_gamma)]:
                extrema[key][0] = min(extrema[key][0], float(np.min(values)))
                extrema[key][1] = max(extrema[key][1], float(np.max(values)))
            for box in boxes:
                deviation = np.max(np.abs(normalized-np.asarray(box['center'])), axis=(1, 2))
                good = ((deviation <= box['radius']) & (k_ratio >= box['k_range'][0]) &
                        (k_ratio <= box['k_range'][1]) & (epsilon_ratio >= box['epsilon_range'][0]) &
                        (epsilon_ratio <= box['epsilon_range'][1]))
                bad = np.flatnonzero(~good)
                box['sample_count'] += count
                box['violation_count'] += len(bad)
                box['violations_by_set'][name] = box['violations_by_set'].get(name, 0)+len(bad)
                for j in bad[:max(0, maximum_examples-len(box['examples']))]:
                    box['examples'].append(dict(state=state_index, time=epoch, locations=name,
                        geometry=identity, chart=points[j].tolist(),
                        gradient_entry_deviation=float(deviation[j]), normalized_gradient=normalized[j].tolist(),
                        k_over_K0=float(k_ratio[j]), epsilon_over_epsilon0=float(epsilon_ratio[j])))
    complete = native_issue is None
    native_count_checks = []
    for index, audit in enumerate(raw.get('state_audits', [])):
        expected_source = sum(len(points) for name, _, points in sets if name == 'source_volume_outer')
        expected_interface = sum(len(points) for name, _, points in sets if name == 'interface_outer')
        for name, expected in [('source_volume_outer', expected_source), ('interface_outer', expected_interface)]:
            native_count_checks.append(dict(state=index, locations=name, expected=expected,
                                            native=audit[name][0], equal=expected == audit[name][0]))
            complete = complete and expected == audit[name][0]
        expected_dense = sum((2*element['degree']+3)**d for element in fields[0].elements)
        native_count_checks.append(dict(state=index, locations='reconstruction_dense_all_regions',
                                        expected=expected_dense, native=audit['reconstruction_dense'][0],
                                        equal=expected_dense == audit['reconstruction_dense'][0]))
        complete = complete and expected_dense == audit['reconstruction_dense'][0]
        if native_rules:
            expected_native = sum(math.prod(native_rules[element['geometry']]) for element in fields[0].elements)
            native_count_checks.append(dict(state=index, locations='reconstruction_native_all_regions',
                expected=expected_native, native=audit['reconstruction_native'][0],
                equal=expected_native == audit['reconstruction_native'][0]))
            complete = complete and expected_native == audit['reconstruction_native'][0]
    for box in boxes:
        box['finite_state_box_covered'] = box['violation_count'] == 0
        box['finite_fraction_inside'] = 1-box['violation_count']/box['sample_count']
        box['momentum_gamma_condition_met'] = extrema['gamma_over_nut0'][0] >= box['required_gamma_over_nut0']
        box['matching_joint_certificate'] = any(
            candidate['certificate_joint_gate'] and candidate['actual_gamma_exact_match']
            for candidate in box['fixed_gamma_candidates'])
        box['all_declared_finite_sets_covered_by_joint_certificate'] = bool(
            complete and box['finite_state_box_covered'] and box['momentum_gamma_condition_met'] and
            box['matching_joint_certificate'] and raw['model_variant'] == 'Regularized')
    accepted_complete = not raw.get('failure') and not raw.get('rejected_state')
    applies = any(box['all_declared_finite_sets_covered_by_joint_certificate'] for box in boxes)
    return dict(schema=1, status='finite_coverage_verified' if applies and accepted_complete else
                'outside_certificate' if any(box['violation_count'] for box in boxes) else 'inconclusive',
                finite_coverage_pass=applies and accepted_complete,
                model_variant=raw['model_variant'], robust_parameters=raw['robust_parameters'],
                accepted_records=len(states), accepted_trajectory_complete=accepted_complete,
                state_scope='Every supplied accepted record, including its initial state; a failed trajectory remains incomplete',
                normalization='A*=Qincoming^T grad(U) Qincoming/G0; k/K0; epsilon/epsilon0; G0=us/(kappa*r),K0=Qstar*us^2,epsilon0=us^3/(kappa*r)',
                incoming_basis=frame.tolist(), per_set_evaluations=set_counts, extrema=extrema,
                native_sampling_complete=complete, native_GLL_metadata_issue=native_issue,
                native_count_checks=native_count_checks, certificate_boxes=boxes,
                inner_scope='Inner stress and source admissibility require separate gates; this is an outer constitutive certificate',
                mms_scope='Manufactured states may intentionally lie outside the neutral reference box; successful equation verification does not extend this certificate',
                excluded_claims=['polynomial-everywhere state coverage', 'PDE state-envelope invariance',
                                 'outer k/epsilon positivity theorem', 'coupled-mode/long-time stability',
                                 'physical terrain qualification'],
                limitations='Double-precision coefficient reconstruction is a finite numerical audit. The outward-rounded interval certificate supplies the separate declared mathematical bound.')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--native', type=Path, required=True)
    p.add_argument('--certificate', type=Path, required=True)
    p.add_argument('--session', type=Path)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--maximum-examples', type=int, default=20)
    a = p.parse_args()
    if platform.system() != 'Linux':
        raise RuntimeError('Scientific numerical coverage reconstruction must run on Linux')
    session = a.session
    if session is None:
        # An adjacent case.xml is a locator only; its actual bytes are hashed.
        guessed = a.native.with_name(a.native.name.removesuffix('.evolution.json')+'.xml')
        if guessed.is_file():
            session = guessed
    report = dict(schema=1, started_utc=datetime.now(timezone.utc).isoformat(),
        inputs=dict(native=dict(path=str(a.native), sha256=sha(a.native)),
                    certificate=dict(path=str(a.certificate), sha256=sha(a.certificate)),
                    session=dict(path=str(session), sha256=sha(session)) if session else None),
        driver_sha256=sha(__file__), numerical_threads=THREADS)
    try:
        report.update(verify(json.loads(a.native.read_text()), json.loads(a.certificate.read_text()),
                             session, a.maximum_examples))
        report['audit_completed'] = True
    except Exception as error:
        report.update(status='audit_failed', finite_coverage_pass=False, audit_completed=False, error=str(error))
    for name, value in report['inputs'].items():
        if value and sha(value['path']) != value['sha256']:
            report.update(status='input_changed', finite_coverage_pass=False, audit_completed=False)
    report['finished_utc'] = datetime.now(timezone.utc).isoformat()
    a.output.parent.mkdir(parents=True, exist_ok=True)
    if a.output.exists():
        raise RuntimeError('Refusing to overwrite an existing coverage receipt')
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({key: report[key] for key in ('status', 'audit_completed', 'finite_coverage_pass')}, allow_nan=False))
    # Noncoverage is a scientific finding, not a corrupted audit.  This exit
    # code attests only successful reconstruction; callers must inspect the
    # explicit finite_coverage_pass field before using a local certificate.
    return 0 if report['audit_completed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
