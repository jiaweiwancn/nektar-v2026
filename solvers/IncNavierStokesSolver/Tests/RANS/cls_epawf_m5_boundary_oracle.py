"""Independent physical surface-L2 audit of the M5 boundary flux projection.

Only native topology/basis metadata and original/lifted coefficients are read.
Exterior quadrature, cofactor measures and the constrained minimum are rebuilt
independently; exported native quadrature weights and matrices are not inputs.
"""
import sys

sys.dont_write_bytecode = True
import numpy as np


def boundary_projection(raw, spaces, geometry, near, require):
    """Check the actual generic-state boundary correction and return a receipt.

    ``spaces`` is the root oracle's list of MappedSpace objects. Its rows method
    supplies independently evaluated native tensor bases, not native exports.
    The generic velocity polynomials are exactly in the native space, so their
    original projection also supplies the uncorrected prescribed boundary data.
    This audit concerns the boundary trace; it does not assert pointwise zero
    divergence of the generic, deliberately unconverged volume field.
    """
    field = spaces[0]
    dimension, size = int(raw['dimension']), int(field.size)
    axes = (0, 2) if dimension == 2 else (0, 1, 2)
    lower, upper = np.asarray(field.lower, float), np.asarray(field.upper, float)
    parameters = np.asarray(raw['map_parameters'], float)
    state, lifted = np.asarray(raw['state'], float), np.asarray(raw['lifted_state'], float)
    offsets = np.asarray(raw['offsets'], int)
    require(dimension in (2, 3) and parameters.shape == (14,), 'boundary projection geometry shape')
    require(state.shape == lifted.shape and np.all(np.isfinite(state)) and np.all(np.isfinite(lifted)),
            'boundary projection state shape/finiteness')
    require(all(offsets[b+1]-offsets[b] == size for b in range(dimension)),
            'boundary projection requires common native velocity space')

    def close(a, b):
        return abs(a-b) < 1e-11

    top = field.boundary(lambda axis, side, element:
                         axis == 2 and side == 1 and close(element['upper'][2], upper[2]))
    fixed = field.boundary(lambda axis, side, element:
                           (axis != 2 and close(element['upper'][axis] if side else element['lower'][axis],
                                                upper[axis] if side else lower[axis]))
                           or (axis == 2 and side == 0 and close(element['lower'][2], lower[2])))
    free = np.asarray(sorted(top-fixed), dtype=int)
    require(len(free) > 0, 'boundary projection has no top interior trace modes')
    require(np.all((free >= 0) & (free < size)), 'boundary projection trace mode range')
    other = np.ones(size, dtype=bool)
    other[free] = False

    # L[a,i] = integral_boundary phi_i n_a dS. M is the physical scalar
    # top-trace mass. All Cartesian components share this scalar trace space.
    functional = np.zeros((3, size))
    mass = np.zeros((size, size))
    area, face_count, quadrature_points = 0., 0, 0
    order = 20
    z, w = np.polynomial.legendre.leggauss(order)
    for element in field.elements:
        lo, hi = np.asarray(element['lower'], float), np.asarray(element['upper'], float)
        for axis in axes:
            free_axes = tuple(a for a in axes if a != axis)
            indices = np.indices((order,)*len(free_axes)).reshape(len(free_axes), -1).T
            for side in (-1, 1):
                position = lo[axis] if side < 0 else hi[axis]
                if not close(position, lower[axis] if side < 0 else upper[axis]):
                    continue
                points = np.zeros((len(indices), 3))
                points[:, axis] = position
                weights = np.ones(len(indices))
                for column, coordinate in enumerate(free_axes):
                    half = (hi[coordinate]-lo[coordinate])/2
                    points[:, coordinate] = lo[coordinate]+half*(z[indices[:, column]]+1)
                    weights *= half*w[indices[:, column]]
                metric = geometry(parameters, points)
                cofactor = side*metric['determinant'][:, None]*metric['inverse'][:, axis, :]
                ids, phi, _ = field.rows(element, points)
                ids = np.asarray(ids, dtype=int)
                require(phi.shape == (len(points), len(ids)) and np.all(np.isfinite(phi)),
                        'boundary projection independent basis shape')
                for physical_axis in axes:
                    np.add.at(functional[physical_axis], ids, phi.T@(weights*cofactor[:, physical_axis]))
                if axis == 2 and side == 1:
                    surface_weights = weights*np.linalg.norm(cofactor, axis=1)
                    local = phi.T@(surface_weights[:, None]*phi)
                    np.add.at(mass, (ids[:, None], ids[None, :]), local)
                    area += float(np.sum(surface_weights))
                face_count += 1
                quadrature_points += len(points)
    require(area > 0 and np.isfinite(area), 'boundary projection positive physical top area')
    top_mass = mass[np.ix_(free, free)]
    require(np.all(np.isfinite(top_mass)) and np.all(np.linalg.eigvalsh(top_mass) > 0),
            'boundary projection positive physical trace mass')
    before = np.stack([state[offsets[b]:offsets[b]+size] for b in range(dimension)])
    after = np.stack([lifted[offsets[b]:offsets[b]+size] for b in range(dimension)])
    delta = after-before
    loads = functional[np.asarray(axes)]
    flux_before = float(np.sum(loads*before))
    flux_after = float(np.sum(loads*after))
    direction = np.zeros_like(delta)
    for block in range(dimension):
        direction[block, free] = np.linalg.solve(top_mass, loads[block, free])
    denominator = float(np.sum(loads*direction))
    require(np.isfinite(denominator) and denominator > 0, 'boundary projection nondegenerate flux constraint')
    expected = -flux_before*direction/denominator
    near(delta[:, free], expected[:, free], 'physical_surface_L2_top_correction', tolerance=1e-10)
    near(delta[:, other], np.zeros_like(delta[:, other]), 'boundary_projection_unchanged_other_velocity_modes', tolerance=1e-10)
    near(flux_after, 0., 'boundary_projection_zero_net_flux', tolerance=1e-10)
    near(float(np.sum(loads*delta)), -flux_before, 'boundary_projection_accounts_for_flux_change', tolerance=1e-10)
    multiplier = flux_before/denominator
    kkt = np.stack([top_mass@delta[b, free]+multiplier*loads[b, free] for b in range(dimension)])
    near(kkt, np.zeros_like(kkt), 'boundary_projection_surface_L2_stationarity', tolerance=1e-10)
    square = sum(float(delta[b]@mass@delta[b]) for b in range(dimension))/area
    require(square >= -1e-24 and np.isfinite(square), 'boundary projection nonnegative physical correction norm')
    rms = float(np.sqrt(max(0., square)))
    predicted_rms = float(abs(flux_before)/np.sqrt(denominator*area))
    near(rms, predicted_rms, 'boundary_projection_minimum_surface_RMS', tolerance=1e-10)
    near(raw['boundary_flux_correction'], rms, 'boundary_projection_reported_surface_RMS', tolerance=1e-10)
    return dict(quadrature_order=order, exterior_faces=face_count, exterior_quadrature_points=quadrature_points,
                top_interior_modes=len(free), cartesian_components=list(axes), top_area=area,
                flux_before=flux_before, flux_after=flux_after, flux_constraint_denominator=denominator,
                physical_surface_rms=rms, predicted_surface_rms=predicted_rms,
                maximum_top_coefficient_error=float(np.max(np.abs(delta[:, free]-expected[:, free]), initial=0.)),
                maximum_other_coefficient_change=float(np.max(np.abs(delta[:, other]), initial=0.)),
                maximum_stationarity_residual=float(np.max(np.abs(kkt), initial=0.)),
                mechanism='prescribed top velocity changed by a physical surface-L2 minimum-norm projection',
                changed_only_top_interior_velocity=True, flux_change_accounted_for_by_boundary_coefficients=True,
                scope='boundary compatibility audit; not a claim of pointwise incompressibility or a manufactured volume mass source')
