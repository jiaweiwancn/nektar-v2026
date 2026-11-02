#!/usr/bin/env python3
"""Independent NumPy geometry audit for the M5 native probes.

No production geometry code or exported metric is used to construct the oracle.
Finite differences and finite quadrature have explicitly finite tolerances;
these checks do not establish temporal or constitutive stability.
"""
import argparse
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
import numpy as np

PARAMETER_NAMES = ('kind', 'dimension', 'height', 'heightSlopeX', 'heightSlopeY',
                   'heightAmplitude', 'heightWaveX', 'heightWaveY', 'slopeX',
                   'slopeY', 'radius', 'amplitude', 'lengthX', 'lengthY')
KINDS = ('Flat', 'Plane', 'Cylinder', 'Gaussian', 'CosineHill')


def parameter_array(parameters):
    if isinstance(parameters, dict):
        record = dict(parameters)
        if isinstance(record['kind'], str):
            record['kind'] = KINDS.index(record['kind'])
        return np.array([record[name] for name in PARAMETER_NAMES], dtype=float)
    array = np.asarray(parameters, dtype=float)
    if array.shape != (14,):
        raise ValueError('normal-map parameters must have length 14')
    return array


def geometry(parameters14, chart, foot_side=0):
    """Analytic map/metrics for chart arrays (...,3), in Cartesian components.

    Parameters follow the native report's 14-array. ``inverse[...,b,a]`` is
    dchart_b/dx_a. ``matching_area`` is always at rho=H0 for the same q1/q2.
    ``foot_side`` means the q1 left/right limit at an exact cosine foot.
    """
    p = parameter_array(parameters14)
    if not np.all(np.isfinite(p)):
        raise ValueError('nonfinite geometry parameter')
    kind, dimension = int(p[0]), int(p[1])
    if p[0] != kind or p[1] != dimension or kind not in range(5) or dimension not in (2, 3):
        raise ValueError('invalid analytic geometry kind/dimension')
    H0, hx, hy, ha, hkx, hky, sx, sy, radius, amp, lx, ly = p[2:]
    q = np.asarray(chart, dtype=float)
    if q.shape[-1:] != (3,) or not np.all(np.isfinite(q)) or H0 <= 0:
        raise ValueError('invalid chart or reference height')
    if dimension == 2 and (hy != 0 or sy != 0 or np.any(q[..., 1] != 0)):
        raise ValueError('nonzero unused 2D direction')
    if np.any(q[..., 2] < 0):
        raise ValueError('negative normal distance')
    x, y, rho = np.moveaxis(q, -1, 0)
    side = np.broadcast_to(np.asarray(foot_side), q.shape[:-1])
    if np.any(~np.isin(side, (-1, 0, 1))):
        raise ValueError('invalid foot side')
    H = H0+hx*x+hy*y+ha*np.sin(hkx*x)*np.cos(hky*y)
    H1 = hx+ha*hkx*np.cos(hkx*x)*np.cos(hky*y)
    H2 = hy-ha*hky*np.sin(hkx*x)*np.sin(hky*y)
    if np.any(H <= 0):
        raise ValueError('nonpositive matching height')
    scale, r = H/H0, rho*H/H0
    zeros, ones = np.zeros_like(x), np.ones_like(x)
    if kind == 2:
        if radius <= 0 or np.any(side != 0):
            raise ValueError('invalid cylinder or foot side')
        angle = x/radius
        normal = np.stack((np.sin(angle), zeros, np.cos(angle)), axis=-1)
        t1 = np.stack((np.cos(angle), zeros, -np.sin(angle)), axis=-1)
        t2 = np.stack((zeros, ones, zeros), axis=-1)
        base = np.stack((radius*np.sin(angle), y, radius*(np.cos(angle)-1)), axis=-1)
        dn1, dn2 = t1/radius, np.zeros_like(t1)
        B = np.einsum('...i,...j->...ij', t1, t1)/radius
        old_curvature = -ones/radius
    else:
        f, fx, fy, fxx, fxy, fyy = [np.zeros_like(x) for _ in range(6)]
        if kind in (0, 1):
            f, fx, fy = sx*x+sy*y, zeros+sx, zeros+sy
        elif kind == 3:
            if lx <= 0 or (dimension == 3 and ly <= 0):
                raise ValueError('invalid Gaussian lengths')
            ix, iy = 1/lx**2, 1/ly**2 if dimension == 3 else 0.
            f = amp*np.exp(-x*x*ix-y*y*iy)
            fx, fy = -2*x*ix*f, -2*y*iy*f
            fxx, fxy, fyy = (4*x*x*ix*ix-2*ix)*f, 4*x*y*ix*iy*f, (4*y*y*iy*iy-2*iy)*f
        elif kind == 4:
            if lx <= 0:
                raise ValueError('invalid cosine half length')
            at_foot = np.abs(np.abs(x)-lx) <= 32*np.finfo(float).eps*np.maximum(np.maximum(1., lx), np.abs(x))
            if np.any((side != 0) & ~at_foot):
                raise ValueError('one-sided metric away from foot')
            interior = (np.abs(x) < lx) & ~at_foot
            a = np.pi/lx
            f = np.where(interior, amp*(1+np.cos(a*x))/2, 0.)
            fx = np.where(interior, -amp*a*np.sin(a*x)/2, 0.)
            fxx = np.where(interior, -amp*a*a*np.cos(a*x)/2, 0.)
            inner_limit = at_foot & (((x < 0) & (side > 0)) | ((x > 0) & (side < 0)))
            fxx = np.where(inner_limit, amp*a*a/2, fxx)
        if kind != 4 and np.any(side != 0):
            raise ValueError('one-sided metric on smooth surface')
        base = np.stack((x, y, f), axis=-1)
        t1, t2 = np.stack((ones, zeros, fx), axis=-1), np.stack((zeros, ones, fy), axis=-1)
        norm = np.sqrt(1+fx*fx+fy*fy)
        normal = np.stack((-fx, -fy, ones), axis=-1)/norm[..., None]
        dn1 = -np.stack((fxx, fxy, zeros), axis=-1)/norm[..., None]
        dn2 = -np.stack((fxy, fyy, zeros), axis=-1)/norm[..., None]
        dn1 -= normal*((fx*fxx+fy*fxy)/norm**2)[..., None]
        dn2 -= normal*((fx*fxy+fy*fyy)/norm**2)[..., None]
        tangents = np.stack((t1, t2), axis=-1)
        gram = np.einsum('...ia,...ib->...ab', tangents, tangents)
        reciprocal = tangents@np.linalg.inv(gram)
        hessian = np.stack((np.stack((fxx, fxy), axis=-1), np.stack((fxy, fyy), axis=-1)), axis=-2)
        B = -(reciprocal@hessian@np.swapaxes(reciprocal, -1, -2))/norm[..., None, None]
        old_curvature = fxx/(1+fx*fx)**1.5
    base_area = np.linalg.norm(np.cross(t1, t2), axis=-1)
    tangent1 = t1/np.linalg.norm(t1, axis=-1)[..., None]
    tangent2 = np.cross(normal, tangent1)
    frame = np.stack((tangent1, tangent2), axis=-1)
    principal = np.linalg.eigvalsh(np.swapaxes(frame, -1, -2)@B@frame)
    if np.any(1+r[..., None]*principal <= 0):
        raise ValueError('folded analytic normal offset')
    A = np.eye(3)+r[..., None, None]*B
    Br = B@np.linalg.inv(A)
    jacobian = np.stack((t1+r[..., None]*dn1+normal*(rho*H1/H0)[..., None],
                         t2+r[..., None]*dn2+normal*(rho*H2/H0)[..., None],
                         normal*scale[..., None]), axis=-1)
    matching_area = np.cross(t1+H[..., None]*dn1+normal*H1[..., None],
                             t2+H[..., None]*dn2+normal*H2[..., None])
    result = dict(point=base+r[..., None]*normal, jacobian=jacobian, inverse=np.linalg.inv(jacobian),
                  determinant=np.linalg.det(jacobian), normal=normal, tangent1=tangent1, tangent2=tangent2,
                  B=B, Br=Br, J=np.linalg.det(A), baseArea=base_area, distance=r, height=H, scale=scale,
                  kappa1=principal[..., 1], kappa2=principal[..., 0], matching_area=matching_area,
                  height_gradient=np.stack((H1, H2), axis=-1), base_tangents=np.stack((t1, t2), axis=-1),
                  old_graph_curvature=old_curvature)
    if any(not np.all(np.isfinite(value)) for value in result.values()):
        raise ValueError('nonfinite analytic metric')
    return result


class Audit:
    def __init__(self):
        self.checks, self.errors = 0, {}

    def require(self, good, message):
        self.checks += 1
        if not good:
            raise RuntimeError(message)

    def near(self, actual, expected, name, tolerance=4e-11):
        a, b = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
        self.require(a.shape == b.shape and np.all(np.isfinite(a)) and np.all(np.isfinite(b)), 'invalid '+name)
        error = float(np.max(np.abs(a-b), initial=0.))/max(1., float(np.max(np.abs(b), initial=0.)))
        self.errors[name] = max(error, self.errors.get(name, 0.))
        self.require(error <= tolerance, f'{name}: {error} > {tolerance}')
        return error


def _load(report):
    return report if isinstance(report, dict) else json.loads(Path(report).read_text())


def _sample(audit, exported, expected):
    aliases = {'r': 'distance', 'H': 'height'}
    for key in ('point', 'normal', 'tangent1', 'tangent2', 'B', 'Br', 'jacobian', 'inverse',
                'determinant', 'baseArea', 'r', 'H', 'scale', 'J'):
        audit.near(exported[key], expected[aliases.get(key, key)], 'analytic_'+key)
    audit.near(exported['kappas'], [expected['kappa1'], expected['kappa2']], 'principal_curvatures')


def _certificate(audit, parameters, lower, upper, minimum_stretch):
    p, lo, hi = parameter_array(parameters), np.asarray(lower), np.asarray(upper)
    kind, dimension, H0, hx, hy, ha = p[:6]
    audit.require(np.all(np.isfinite(lo)) and np.all(np.isfinite(hi)) and hi[0] > lo[0]
                  and hi[2] > lo[2] >= 0, 'invalid certified chart bounds')
    audit.require((dimension == 2 and lo[1] == hi[1] == 0)
                  or (dimension == 3 and hi[1] > lo[1]), 'invalid certified span bounds')
    hmin = H0+min(hx*lo[0], hx*hi[0])+min(hy*lo[1], hy*hi[1])-abs(ha)
    hmax = H0+max(hx*lo[0], hx*hi[0])+max(hy*lo[1], hy*hi[1])+abs(ha)
    audit.require(hmin > 0, 'whole-domain positive height certificate')
    rmax, bound = hi[2]*hmax/H0, 0.
    if kind == 2:
        audit.require(p[10] > 0 and (hi[0]-lo[0])/p[10] < 2*np.pi, 'nonwrapping outward cylinder certificate')
    else:
        if kind == 3:
            # e^-|v|^2 ||2 vv^T-I|| <= 1; this spectral bound is
            # independent of the native implementation's row-sum bound.
            bound = 2*abs(p[11])/min(p[12], p[13] if dimension == 3 else p[12])**2
        elif kind == 4:
            bound = abs(p[11])*np.pi**2/(2*p[12]**2)
        audit.require(rmax*bound < 1-minimum_stretch, 'global graph normal-reach certificate')
    return dict(minimum_height=float(hmin), maximum_normal_distance=float(rmax),
                graph_hessian_bound=float(bound), normal_stretch_lower_bound=float(1-rmax*bound))


def circle_normal_balance(radius=20., height=2.):
    """Independent Cartesian FD of U=c*r*t and its centrifugal pressure law."""
    c, s, r = .17, .73, 1.1
    theta = s/radius
    point = np.array([(radius+r)*np.sin(theta), 0., (radius+r)*np.cos(theta)-radius])
    normal = np.array([np.sin(theta), 0., np.cos(theta)])
    def velocity(p):
        angle = np.arctan2(p[0], p[2]+radius)
        distance = np.hypot(p[0], p[2]+radius)-radius
        return c*distance*np.array([np.cos(angle), 0., -np.sin(angle)])
    exact = -(c*r)**2/(radius+r)
    errors = []
    for h in (2e-4, 1e-4, 5e-5):
        gradient = np.zeros((3, 3))
        for axis in (0, 2):
            shift = np.eye(3)[axis]*h
            gradient[:, axis] = (velocity(point+shift)-velocity(point-shift))/(2*h)
        errors.append(abs(normal@(gradient@velocity(point))-exact))
    def pressure(distance):
        return c*c*(distance*distance/2-radius*distance+radius*radius*np.log1p(distance/radius))
    h = 1e-4
    dp = (pressure(height+h)-pressure(height-h))/(2*h)
    expected_dp = (c*height)**2/(radius+height)
    return dict(normal_acceleration=exact, fd_errors=errors, pressure_gradient=dp,
                expected_pressure_gradient=expected_dp, pressure_error=abs(dp-expected_dp))


def validate_geometry(report):
    data, audit = _load(report), Audit()
    audit.require(data['schema'] == 1 and data['success'], 'native geometry probe failed')
    cases = {case['name']: case for case in data['cases']}
    required = {'flat_2', 'flat_3', 'tilted_plane_2', 'tilted_plane_3', 'outward_circle_2',
                'outward_circle_3', 'variable_height_circle_2', 'variable_height_circle_3',
                'gaussian_2', 'gaussian_3', 'variable_height_gaussian_3',
                'cosine_explicit_feet_2', 'cosine_explicit_feet_3'}
    audit.require(required <= cases.keys(), 'missing positive geometry coverage')
    certificates = {}
    for case in cases.values():
        parameters = parameter_array(case['parameters'])
        certificates[case['name']] = _certificate(audit, parameters, case['lower'], case['upper'], case['minimum_stretch'])
        for item in case['samples']:
            native = item['map']; q = np.asarray(native['chart']); expected = geometry(parameters, q)
            _sample(audit, native, expected)
            audit.near(item['inverse_point'], q, 'inverse_roundtrip', 4e-10)
            for axis, cofactor in enumerate(item['cofactors']):
                if cofactor is not None:
                    audit.near(cofactor, expected['determinant']*expected['inverse'][axis], 'face_cofactor')
            fd = item['finite_differences']
            audit.require(len(fd) >= 3, 'missing finite-difference refinement sequence')
            for level in fd:
                audit.near(level['map_jacobian'], expected['jacobian'], 'fd_map_jacobian', 2e-7)
                audit.near(level['inverse_jacobian'], expected['inverse'], 'fd_inverse_jacobian', 2e-7)
                audit.near(level['normal_gradient'], expected['Br'], 'fd_physical_normal_gradient', 2e-7)
                h, piola = level['chart_h'], np.zeros(3)
                for axis in ((0, 2) if parameters[1] == 2 else (0, 1, 2)):
                    shift = np.eye(3)[axis]*h
                    gp, gm = geometry(parameters, q+shift), geometry(parameters, q-shift)
                    piola += (gp['determinant']*gp['inverse'][axis]-gm['determinant']*gm['inverse'][axis])/(2*h)
                audit.near(level['cofactor_divergence'], piola, 'independent_fd_piola', 3e-8)
                audit.near(level['cofactor_divergence'], np.zeros(3), 'fd_piola_identity', 2e-7)
            for column, floor in ((0, 2e-8), (1, 1e-7), (2, 2e-8)):
                audit.require(fd[-1]['errors'][column] <= max(floor, 1.25*fd[0]['errors'][column]), 'FD refinement failed')
            matching = item['matching']; qm = q.copy(); qm[2] = parameters[2]
            gm = geometry(parameters, qm); _sample(audit, matching['sample'], gm)
            area = gm['matching_area']; normal = area/np.linalg.norm(area)
            audit.near(matching['area'], area, 'matching_area')
            audit.near(matching['N_H'], normal, 'matching_normal')
            audit.near(matching['m_from_area'], np.dot(matching['U'], area)/gm['baseArea'], 'matching_mass_flux')
            audit.near(matching['m_from_formula'], matching['m_from_area'], 'variable_height_mass_identity')
            if parameters[0] in (3, 4) and parameters[1] == 2:
                audit.near(np.trace(native['B']), -expected['old_graph_curvature'], 'old_graph_curvature_sign')
    for foot in data['foot_interfaces']:
        parameters = parameter_array(cases[foot['case']]['parameters']); q = np.asarray(foot['chart'])
        for label, side in (('left', -1), ('right', 1)):
            expected = geometry(parameters, q, side); _sample(audit, foot[label], expected)
            audit.near(foot[label+'_lateral_cofactor'], expected['determinant']*expected['inverse'][0], 'C1_lateral_cofactor')
        audit.near(foot['left']['point'], foot['right']['point'], 'C1_position')
        audit.near(foot['left']['normal'], foot['right']['normal'], 'C1_normal')
        audit.near(foot['left_lateral_cofactor'], foot['right_lateral_cofactor'], 'C1_shared_area')
        audit.require(foot['B_jump'] > 1e-3, 'C1 curvature jump was hidden')
    rejections = data['negative_controls']; names = {record['name'] for record in rejections}
    audit.require(data['rejections'] == len(rejections) >= 26 and all(record['error'] for record in rejections), 'invalid rejection receipts')
    audit.require({'graph_tube_beyond_certified_reach', 'actual_folded_graph_offset', 'wrapping_cylinder_arc',
                   'nonpositive_height_between_positive_corners', '2D_spanwise_plane_slope',
                   '2D_spanwise_height_slope', 'nonpositive_cylinder_radius', 'nonfinite_parameter',
                   'negative_rho_domain', 'empty_normal_domain', 'off_foot_one_sided_metric'} <= names,
                  'missing negative geometry coverage')
    audit.require(data['height_slope_omission_control'] > 1e-4, 'missing tilted-interface control')
    circle = circle_normal_balance()
    audit.require(max(circle['fd_errors']) < 2e-8 and circle['pressure_error'] < 2e-7 and circle['pressure_gradient'] > 0,
                  'independent circular normal-momentum identity')
    return dict(success=True, checks=audit.checks, errors=audit.errors, cases=len(cases), rejections=len(rejections), certificates=certificates,
                circle_normal_balance=circle, qualification='analytic geometry and finite-difference diagnostics; no CFD stability claim')


def validate_wall(report):
    data, audit = _load(report), Audit()
    audit.require(data['schema'] == 1, 'unexpected wall probe schema')
    raw = np.asarray(data['raw_residual'])
    regional = np.sum(np.asarray(data['region_residuals']), axis=0)
    audit.near(regional, raw, 'regional_raw_residual_sum', 2e-10)
    levels = data['weak_metric']
    audit.require(len(levels) >= 3, 'missing weak metric quadrature sequence')
    for level in levels:
        for field in level['fields']:
            weak, seam = np.asarray(field['weak_residual']), np.asarray(field['seam_load'])
            audit.require(weak.shape == seam.shape == (field['size'], 3), 'weak metric row shape')
            audit.near(field['maximum'], np.max(np.abs(weak), initial=0.), 'reported_all_mode_metric_max', 1e-12)
            audit.near(field['seam_maximum'], np.max(np.abs(seam), initial=0.), 'reported_seam_max', 1e-12)
    for first, last in zip(levels[0]['fields'], levels[-1]['fields']):
        audit.require(last['maximum'] < 1e-10 and last['maximum'] <= max(1e-10, 1.25*first['maximum']), 'weak Piola quadrature did not converge')
        audit.near(last['seam_load'], np.zeros_like(last['seam_load']), 'equal_opposite_seam_load', 1e-10)
        audit.near(last['constant_flux_boundary'], np.zeros(3), 'constant_flux_divergence_theorem', 1e-10)
    ranks = data['ranks']; owners = data['cell_ownership']
    audit.require(len({entry[0] for entry in owners}) == len(owners) and all(0 <= entry[1] < ranks for entry in owners), 'invalid cell ownership')
    shared = data['interface']['shared_advective_seam_load']
    audit.near(shared, np.zeros_like(shared), 'shared_advective_interface_loads', 1e-10)
    return dict(success=True, checks=audit.checks, errors=audit.errors,
                qualification='actual native weak loads at reported finite quadrature orders')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report', nargs='?')
    parser.add_argument('--geometry-report')
    parser.add_argument('--wall-report', action='append', default=[])
    parser.add_argument('--output')
    args = parser.parse_args()
    filename = args.geometry_report or args.report
    if not filename:
        parser.error('a geometry report is required')
    result = dict(geometry=validate_geometry(filename), walls=[validate_wall(path) for path in args.wall_report])
    text = json.dumps(result, indent=2, allow_nan=False)+'\n'
    if args.output:
        Path(args.output).write_text(text)
    else:
        print(text, end='')


if __name__ == '__main__':
    main()
