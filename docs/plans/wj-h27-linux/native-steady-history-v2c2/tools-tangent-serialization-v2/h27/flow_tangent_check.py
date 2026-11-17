#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent Linux check of native frozen-scalar planar stress Jacobians.

Native owned rule coordinates/weights/maps are inputs. The Jacobi basis and
full 3D WJ value/root are independently evaluated by the existing physical
oracle. Finite differences change only G or native velocity coefficients;
k, epsilon, auxiliary, radius, chi, quadrature, and scalar ownership stay fixed.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import time

import physical_oracle as po
import numpy as np

AXES = (0, 2)


def finite(value, shape=None):
    result = np.asarray(value, dtype=float)
    po.require(np.all(np.isfinite(result)) and (shape is None or result.shape == shape), 'nonfinite or malformed probe array')
    return result


def relative(actual, expected):
    return float(np.max(abs(actual-expected), initial=0.))/max(float(np.max(abs(expected), initial=0.)), 1e-30)


def field(component, x, r, direction=False):
    if direction:
        return (.02+.01*x-.013*r+.007*x*x*r-.004*x*r*r if component == 0 else
                -.017+.009*x+.011*r+.005*x*r+.002*x*x*r*r)
    return (.03+.07*x+.11*r+.03*x*r+.017*x*x*r if component == 0 else
            .01-.015*x+.023*r+.006*x*r+.009*x*r*r)


def wall(k, epsilon, G, nu, aux, r, H):
    k, epsilon, aux, r = np.broadcast_arrays(k, epsilon, aux, r)
    M = G+G.transpose(0, 2, 1)
    s = np.clip((r/H-.4)/.35, 0., 1.)
    chi = s**3*(10-15*s+6*s*s)
    R = 2*k[:, None, None]*np.eye(3)/3
    R[:, 0, 2] = R[:, 2, 0] = -aux*M[:, 0, 2]
    R[:, 1, 2] = R[:, 2, 1] = -aux*M[:, 1, 2]
    active = chi > 0
    if np.any(active):
        full = po.full_wj(k[active], epsilon[active], G[active])[0]
        R[active] = (1-chi[active, None, None])*R[active]+chi[active, None, None]*full
    return nu*M-R+2*k[:, None, None]*np.eye(3)/3


def finite_columns(k, epsilon, G, nu, aux, r, H, step):
    columns = []
    for row in AXES:
        for column in AXES:
            plus, minus = G.copy(), G.copy()
            plus[:, row, column] += step; minus[:, row, column] -= step
            columns.append((wall(k, epsilon, plus, nu, aux, r, H)-
                            wall(k, epsilon, minus, nu, aux, r, H))/(2*step))
    return np.stack(columns, axis=1)


def sparse_input(packet, size):
    po.require(packet['rows'] == packet['cols'] == size, 'wrong full tangent matrix dimensions')
    matrix = np.zeros((size, size)); seen = set()
    for row, column, value in packet['entries']:
        po.require(type(row) is int and type(column) is int and 0 <= row < size and
                   0 <= column < size and (row, column) not in seen and math.isfinite(value),
                   'invalid or duplicate native sparse tangent entry')
        seen.add((row, column)); matrix[row, column] = value
    return matrix


def verify(data, deadline, diagnostics=None):
    po.require(data['schema'] == 'WJH27FrozenPlanarFlowTangentProbeV1' and data['status'] == 'RECORDED', 'wrong tangent probe schema/status')
    n = data['velocity_size']; po.require(type(n) is int and 0 < n <= 512, 'bounded native probe velocity inventory differs')
    velocity, direction = finite(data['velocity'], (2*n,)), finite(data['direction'], (2*n,))
    matrix, fast = (sparse_input(data[key], 2*n) for key in ('matrix', 'fast_matrix'))
    action, combined = finite(data['assembled_action'], (2*n,)), finite(data['combined_direction_action'], (2*n,))
    errors = dict(native_matrix_action=relative(matrix@direction, action),
                  native_combined_direction=relative(action, combined), backend_matrix=relative(fast, matrix))
    sign = finite(data['basis_sign_transform'], (n,))
    po.require(np.array_equal(sign, np.where(np.arange(n) % 3 == 1, -1., 1.)), 'signed-basis control is missing or changed')
    errors['signed_basis_covariance'] = relative(finite(data['sign_transformed_action'], (2*n,)),
                                                np.tile(sign, 2)*action)
    po.require(errors['native_matrix_action'] < 2e-12 and errors['native_combined_direction'] < 2e-12 and
               errors['backend_matrix'] < 2e-11 and errors['signed_basis_covariance'] < 2e-12,
               'native four-column/scatter/backend/signed-basis identity failed')
    po.require(data.get('finite_difference_serialization') == 'MaxDigits10V2',
               'missing full-precision FD serialization and round-trip diagnostics')
    serialization = []
    if diagnostics is not None:
        diagnostics['finite_difference_serialization'] = serialization
    for row in data['finite_differences']:
        value = finite(row['action'], (2*n,))
        raw = finite(row['native_scalar_serialized_action'], (2*n,))
        native_error = float(row['native_in_memory_relative_error'])
        native_difference = float(row['native_in_memory_max_abs_error'])
        difference = float(np.max(abs(value-action)))
        scale = max(float(np.max(abs(action))), 1e-30)
        po.require(math.isfinite(native_error) and math.isfinite(native_difference) and
                   native_error >= 0 and native_difference >= 0 and np.array_equal(value, raw),
                   'FD vector did not round-trip identically to native scalar serialization')
        po.require(abs(difference-native_difference) <= 16*np.finfo(float).eps*scale and
                   abs(relative(value, action)-native_error) <= 16*np.finfo(float).eps,
                   'serialized FD error differs from its native in-memory measurement')
        legacy_error = relative(finite(row['legacy_precision6_action'], (2*n,)), action)
        serialization.append(dict(step=row['step'], exact_scalar_vector_round_trip=True,
            native_in_memory_max_abs_error=native_difference,
            native_in_memory_relative_error=native_error,
            decoded_relative_error=relative(value, action),
            round_trip_relative_error_discrepancy=abs(relative(value, action)-native_error),
            legacy_precision6_relative_error=legacy_error,
            legacy_precision6_gate='PASS' if legacy_error <= 2e-8 else 'FAIL',
            unchanged_FD_gate=2e-8))
        errors['native_velocity_FD_'+str(row['step'])] = relative(value, action)
        po.require(errors['native_velocity_FD_'+str(row['step'])] <= 2e-8, 'native full-row velocity FD disagrees')
    po.require([row['step'] for row in data['finite_differences']] == [2e-5, 1e-5], 'changed prescribed velocity-FD ladder')

    point_rows = []
    po.require(len(data['points']) == 9, 'point derivative coverage changed')
    for point in data['points']:
        po.require(time.monotonic() < deadline, 'tangent checker cooperative wall cap reached')
        G = finite(point['G'], (3, 3))[None, :, :]
        actual = finite(point['dVdG'], (4, 3, 3))
        row = dict(name=point['name'], errors=[])
        for step in (2e-5, 1e-5):
            expected = finite_columns(np.array([point['k']]), np.array([point['epsilon']]), G,
                point['nu'], np.array([point['auxiliary']]), np.array([point['r']]), point['H'], step)[0]
            error = float(np.max(abs(actual-expected)))/point['mu']
            row['errors'].append(error)
            po.require(error <= 2e-8, 'independent full-tensor point FD failed: '+point['name'])
        if point['r'] <= .4*point['H']:
            exact = np.zeros((4, 3, 3)); exact[0, 0, 0] = exact[3, 2, 2] = 2*point['nu']
            for column in (1, 2): exact[column, 0, 2] = exact[column, 2, 0] = point['nu']+point['auxiliary']
            row['exact_lower_wall_error_over_mu'] = float(np.max(abs(actual-exact)))/point['mu']
            po.require(row['exact_lower_wall_error_over_mu'] <= 2e-12, 'lower-wall exact tangent differs')
        if point['name'] == 'negative_shear_unit_control':
            po.require(actual[1, 0, 2] < 0, 'negative Newton shear tangent was replaced by positive diffusion')
            row['negative_shear_retained'] = float(actual[1, 0, 2])
        point_rows.append(row)

    independent = [np.zeros_like(matrix), np.zeros_like(matrix)]
    stress_actions = [np.zeros(2*n), np.zeros(2*n)]
    native_fd_values = [finite(row['action'], (2*n,)) for row in data['finite_differences']]
    representation, quadrature, total_points, negative_signs = 0., 0., 0, 0
    coverage = set(); owners = set()
    for rule in data['rules']:
        po.require(time.monotonic() < deadline, 'tangent checker cooperative wall cap reached')
        owner = rule['geometry']; po.require(owner not in owners, 'duplicate probe rule owner'); owners.add(owner)
        ids = np.asarray(rule['u'], dtype=int); signs = finite(rule['us'], (16,))
        po.require(ids.shape == (16,) and len(set(ids)) == 16 and np.all((ids >= 0) & (ids < n)) and
                   np.all(np.isin(signs, (-1., 1.))), 'invalid signed probe map')
        coverage.update(ids); negative_signs += int(np.sum(signs < 0))
        x, z, weight = (finite(rule[key]) for key in ('x', 'z', 'w'))
        po.require(x.ndim == 1 and x.shape == z.shape == weight.shape and np.all(weight > 0), 'invalid physical quadrature data')
        r = z-data['d0']; total_points += len(x)
        xl, xr, rl, rr = rule['bounds']; area = (xr-xl)*(rr-rl)
        po.require(area > 0 and np.all((x >= xl-1e-13) & (x <= xr+1e-13) &
                                      (r >= rl-1e-13) & (r <= rr+1e-13)), 'quadrature escaped its owner')
        quadrature = max(quadrature, abs(weight.sum()-area)/area,
                         abs(weight@x-area*(xl+xr)/2)/area, abs(weight@r-area*(rl+rr)/2)/area)
        b, dx, dr = po.basis(rule, 'u', x, r)
        G, dG = np.zeros((len(x), 3, 3)), np.zeros((len(x), 3, 3))
        for component, axis in enumerate(AXES):
            v, d = velocity[component*n+ids], direction[component*n+ids]
            G[:, axis, 0], G[:, axis, 2] = dx@v, dr@v
            dG[:, axis, 0], dG[:, axis, 2] = dx@d, dr@d
            representation = max(representation, float(np.max(abs(b@v-field(component, x, r)))),
                                 float(np.max(abs(b@d-field(component, x, r, True)))))
        k, epsilon = data['Kin']*(1+.02*x+.03*r), data['Ein']*(1+.01*x+.07*r)
        aux = np.where(r < data['H'], .0003+.001*r+.00001*x, 0.)
        derivatives = (dx, dr)
        for level, step in enumerate((2e-5, 1e-5)):
            tangent = finite_columns(k, epsilon, G, data['nu'], aux, r, data['H'], step)
            for row, row_axis in enumerate(AXES):
                for column in range(2):
                    block = np.zeros((16, 16))
                    for test, test_axis in enumerate(AXES):
                        for trial in range(2):
                            factor = weight*tangent[:, 2*column+trial, row_axis, test_axis]
                            block += derivatives[test].T@(factor[:, None]*derivatives[trial])
                    independent[level][np.ix_(row*n+ids, column*n+ids)] += block
            # Independent value-only FD in the full coefficient direction;
            # this route never reads native dVdG or the native Jacobian matrix.
            dV = (wall(k, epsilon, G+step*dG, data['nu'], aux, r, data['H'])-
                  wall(k, epsilon, G-step*dG, data['nu'], aux, r, data['H']))/(2*step)
            for row, axis in enumerate(AXES):
                stress_actions[level][row*n+ids] += dx.T@(weight*dV[:, axis, 0])+dr.T@(weight*dV[:, axis, 2])
    po.require(coverage == set(range(n)) and owners, 'probe lacks full signed global velocity coverage')
    po.require(representation <= 2e-12 and quadrature <= 2e-12, 'independent basis/profile or physical quadrature moments failed')
    for level, step in enumerate((2e-5, 1e-5)):
        errors['independent_full_matrix_'+str(step)] = relative(independent[level], matrix)
        errors['independent_coefficient_direction_'+str(step)] = relative(stress_actions[level], action)
        errors['independent_vs_native_value_FD_'+str(step)] = relative(stress_actions[level], native_fd_values[level])
        po.require(max(errors[key] for key in errors if key.endswith(str(step))) <= 2e-8,
                   'independent signed global tangent matrix/directional action failed')
    errors['independent_FD_halving_matrix_change'] = relative(independent[0], independent[1])
    po.require(errors['independent_FD_halving_matrix_change'] <= 2e-8, 'independent tangent FD halving is unresolved')
    rejections = {row['name']: row['reason'] for row in data['rejections']}
    po.require(set(rejections) == {'owner_mismatch', 'point_count', 'nonfinite_tangent', 'negative_weight',
                                  'duplicate_global_id', 'spanwise_gradient', 'nonpositive_k'} and
               len(rejections) == len(data['rejections']) and all(rejections.values()), 'invalid-input controls are incomplete')
    return dict(errors=errors, finite_difference_serialization=serialization,
        point_checks=point_rows, basis_representation_MAX=representation,
        physical_quadrature_moment_MAX=quadrature, native_rule_owners=len(owners), native_points=total_points,
        negative_signed_map_entries=negative_signs,
        independent_signed_basis_transformation='EXERCISED_FROM_ACTUAL_NATIVE_RULES',
        negative_sign_coverage='EXERCISED' if negative_signs else 'NOT_EXERCISED_BY_THIS_NATIVE_FIXTURE',
        matrix_asymmetry_over_MAX=float(np.max(abs(matrix-matrix.T)))/max(float(np.max(abs(matrix))), 1e-30),
        invalid_input_controls=rejections,
        scope='Frozen-scalar full-row stress Jacobian only; no pressure/transport/time/Gamma block, nonlinear convergence, or speed qualification')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=90.)
    args = parser.parse_args()
    po.require(sys.platform.startswith('linux'), 'numerical tangent checking is Linux-only')
    po.require(math.isfinite(args.seconds) and args.seconds > 0 and not args.output.exists(), 'finite cap and immutable output required')
    began = time.monotonic()
    report = dict(schema='WJH27FrozenPlanarFlowTangentCheckV1', status='FAIL', input=po.identity(args.input),
        checker=po.identity(__file__), independent_basis=po.identity(po.fixtures.__file__),
        independent_full_WJ=po.identity(po.__file__), workers=1, wall_cap_seconds=args.seconds)
    try:
        report.update(verify(json.loads(args.input.read_text()), began+args.seconds, report)); report['status'] = 'PASS'
    except Exception as error:
        report.update(exception=type(error).__name__, reason=str(error))
    report['elapsed_seconds'] = time.monotonic()-began
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps({key: report[key] for key in ('status', 'reason', 'elapsed_seconds') if key in report}))
    if report['status'] != 'PASS': raise SystemExit(2)


if __name__ == '__main__':
    main()
