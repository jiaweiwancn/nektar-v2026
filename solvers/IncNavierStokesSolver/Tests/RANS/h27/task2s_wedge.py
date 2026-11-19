#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Frozen wedge modes and angular fits to the retained coupled solution."""

import argparse
import math
from pathlib import Path
import sys
import time

import numpy as np
from scipy.optimize import brentq

from physical_oracle import NUMERICAL, PhysicalOracle, load_native, load_reference, require
from resources import atomic, identity


FIELDS = ("U", "W", "Ux", "Ur", "Wx", "Wr", "Pi")
BUDGETS = np.array([NUMERICAL["U"], NUMERICAL["W"], *([NUMERICAL["G"]]*4), NUMERICAL["Pi"]])
RADII = (.0035, .005, .0071, .01, .0141, .02, .0283, .04)
ANGLES = (14, 24, 34, 44, 54, 64, 74)
HOLD_RADII, HOLD_ANGLES = (.0071, .02), (24, 54)
ANNULI = ((.0035, .02), (.0071, .04))
NEAR = dict(radii_code=(.0025,.003,.0035,.004,.0045,.005,.0055,.006),
            angles_deg=(12,18,24,30,36,42,48),heldout_radii=(.0035,.005),
            heldout_angles_deg=(18,36),annuli_code=((.0025,.005),(.0035,.006)))
ANCHOR_RADIUS = .01
BACKGROUND_MODELS = (("exact_reference", 3), ("exact_reference", 4),
                     ("local_auxiliary", 3), ("local_auxiliary", 4))


def real_columns(values):
    return np.stack((values.real, values.imag), axis=-1).reshape((*values.shape[:-1], 4))


class Wedge:
    """Use s=15-x, u_s=-U, and positive r and W. The outlet normal is -e_s."""

    def __init__(self, nu, mu):
        require(0 < nu < mu, "the distinct-root pencil requires positive auxiliary viscosity")
        self.nu, self.mu = nu, mu
        self.sin_theta, self.cos_theta = math.sqrt(nu/mu), math.sqrt(1-nu/mu)
        self.theta = math.asin(self.sin_theta)
        self.a = np.exp(1j*np.array([self.theta, math.pi-self.theta]))
        self.lambda_ = brentq(lambda exponent: self.closed_determinants(exponent)[0],
                              0., 1., xtol=1e-14)
        matrix, _ = self.boundary_matrices(self.lambda_)
        _, sigma, vh = np.linalg.svd(matrix)
        self.null_residual = float(sigma[-1]/sigma[0])
        self.vector = vh[-1]
        angles = np.deg2rad([0., *ANGLES, 90.])
        columns = self.monomial_fields(ANCHOR_RADIUS*np.cos(angles), ANCHOR_RADIUS*np.sin(angles))
        values = (columns[:, 3]+columns[:, 4]) @ self.vector
        anchor = int(np.argmax(np.abs(values)))
        require(abs(values[anchor]) > 1e-12, "wedge eigenmode lacks a normalizing shear")
        self.vector = self.vector / values[anchor]
        self.anchor_angle_deg = float(np.rad2deg(angles[anchor]))

    def monomial(self, s, r, power, ds=0, dr=0):
        z = np.asarray(s)[:, None]+np.asarray(r)[:, None]*self.a
        factor = math.prod(power-j for j in range(ds+dr))
        return factor*self.a**dr*z**(power-ds-dr)

    def pressure(self, s, r, ds=0, dr=0):
        # The s-momentum equation determines Pi_s. For lambda<1 the
        # homogeneous pressure of degree lambda-1 has no additive constant.
        p = 1+self.lambda_
        multiplier = p*(p-1)*((2*self.nu-self.mu)*self.a+self.mu*self.a**3)
        return multiplier*self.monomial(s, r, p-2, ds, dr)

    def monomial_fields(self, s, r):
        # u_s=psi_r and W=-psi_s. Reflect the physical x component and x
        # derivative when returning [U,W,Ux,Ur,Wx,Wr,Pi] columns.
        p = 1+self.lambda_
        return real_columns(np.stack((
            -self.monomial(s,r,p,0,1), -self.monomial(s,r,p,1,0),
            self.monomial(s,r,p,1,1), -self.monomial(s,r,p,0,2),
            self.monomial(s,r,p,2,0), -self.monomial(s,r,p,1,1),
            self.pressure(s,r)), axis=1))

    def fields(self, s, r):
        return self.monomial_fields(s, r) @ self.vector

    def boundary_matrices(self, exponent):
        # At r=0,s=1 the scaled no-slip streamfunction rows are 1,a.
        # At s=0,r=1, t_x=2*nu*(u_s)_s-Pi and
        # t_r=-mu*((u_s)_r+W_s). The natural rows are scaled by mu*p*(p-1).
        psi = np.ones_like(self.a)
        az = (1j*np.array([self.theta, math.pi-self.theta]))
        shear = (1-self.a**2)*np.exp((exponent-1)*az)
        normal = self.a*shear
        minus_w = np.exp(exponent*az)
        return (real_columns(np.array([psi,self.a,shear,normal])),
                real_columns(np.array([psi,self.a,minus_w,normal])))

    def closed_determinants(self, exponent):
        d, c = self.sin_theta, self.cos_theta
        orig = 8*d*d*(d*d-np.cos((math.pi-2*self.theta)*exponent)-c*c*np.cos(math.pi*exponent))
        transverse_essential = 4*d*c*c*np.sin(math.pi*exponent)
        return orig, transverse_essential

    def audit(self):
        grid = np.linspace(.01, .99, 99)
        determinants = np.array([[np.linalg.det(m) for m in self.boundary_matrices(x)] for x in grid])
        formulas = np.column_stack(self.closed_determinants(grid))
        det_error = np.max(np.abs(determinants-formulas), axis=0)
        nu, mu, p = self.nu, self.mu, 1+self.lambda_
        polynomial_error = float(np.max(np.abs(self.a**4+(4*nu/mu-2)*self.a**2+1)))
        angle = np.deg2rad(ANGLES)
        s, r = ANCHOR_RADIUS*np.cos(angle), ANCHOR_RADIUS*np.sin(angle)
        d30 = self.monomial(s,r,p,3,0)
        d21 = self.monomial(s,r,p,2,1)
        d12 = self.monomial(s,r,p,1,2)
        d03 = self.monomial(s,r,p,0,3)
        pi_s, pi_r = self.pressure(s,r,1,0), self.pressure(s,r,0,1)
        residual = np.stack(((2*nu-mu)*d21+mu*d03-pi_s, (mu-2*nu)*d12-mu*d30-pi_r), axis=1)
        norm = max(np.max(np.abs(pi_s)), np.max(np.abs(pi_r)))
        direct_residual = float(np.max(np.abs(residual))/norm)
        mode_residual = float(np.max(np.abs(real_columns(residual)@self.vector))/norm)
        radius = np.array([.004, .01, .032])
        ground = self.fields(radius, np.zeros_like(radius))
        outlet = self.fields(np.zeros_like(radius), radius)
        traction = np.column_stack((2*nu*outlet[:, 2]-outlet[:, 6],
                                    mu*(outlet[:, 3]+outlet[:, 4])))
        boundary_error = max(float(np.max(np.abs(ground[:, :2]))/ANCHOR_RADIUS),
                             float(np.max(np.abs(traction))/mu))
        div_error = float(np.max(np.abs(self.fields(s, r)[:, [2, 5]].sum(axis=1))))
        require(max(*det_error, self.null_residual, direct_residual, mode_residual, boundary_error,
                    div_error, polynomial_error) < 5e-12, "independent frozen wedge algebra check failed")
        return dict(
            status="PASS",
            coefficients=dict(nu=nu, mu=mu, theta_rad=self.theta, characteristic_polynomial_error=polynomial_error),
            original=dict(velocity_exponent=self.lambda_, gradient_pressure_exponent=self.lambda_-1,
                          normalized_boundary_null_residual=self.null_residual, real_roots_in_open_unit_interval=1,
                          determinant_formula="8 d^2 [d^2-cos((pi-2 theta) lambda)-c^2 cos(pi lambda)]",
                          uniqueness_basis="Strictly increasing between opposite endpoint signs"),
            diagnostic_amendment=dict(real_roots_in_open_unit_interval=0,
                determinant_formula="4 d c^2 sin(pi lambda)",
                outward_condition="W=0, t_x=Vxx-Pi prescribed, t_r reaction",
                lower_pressure_compatibility="For smooth incompressible fields, Wr=Ux=0 and Pi=-t_x. In the collar Pi=Vxx-t_x."),
            determinant_error=det_error.tolist(), momentum_monomial_relative_max=direct_residual,
            momentum_combined_relative_max=mode_residual, divergence_max=div_error,
            original_boundary_max=boundary_error, signed_real_column_coefficients=self.vector.tolist(),
            anchor=dict(radius_code=ANCHOR_RADIUS, angle_degrees=self.anchor_angle_deg,
                        unit_response="physical Ur+Wx equals +1 at the anchor"))


def measure(path, ref, deadline, points):
    data, archive = load_native(path)
    require(not data["manufactured"] and archive["steady_mode"] and archive["steady_solver"]["converged"]
            and data["step"] == 0 and data["time"] == 0., "wedge requires a converged zero-time flat state")
    oracle = PhysicalOracle(data, ref, 20, deadline)
    right = oracle.geometry["right"]
    tips = [e for e in oracle.elements if e["bounds"][1] == right and e["bounds"][2] == 0.]
    require(len(tips) == 1, "multiple ground–outlet vertex owners in one mesh")
    vertex = tips[0]
    xl, xr, rl, rr = vertex["bounds"]
    width = np.array([xr-xl, rr-rl])
    K = float(oracle.trace(vertex, "k", np.array([right]))[0])
    aux0 = float(ref.aux(np.array([K]), np.zeros(1))[0])
    mode = Wedge(ref.nu, ref.nu+aux0)
    lower_elements = [e for e in oracle.elements if e["bounds"][1] >= right-max(RADII)-1e-12
                      and e["bounds"][2] <= max(RADII)+1e-12]
    fields = np.full((len(points), len(FIELDS)), np.nan)
    owners = []
    observations = []
    for index, (s, r) in enumerate(points):
        oracle.clock()
        x = right-s
        require(r < .4*ref.H, "corner annulus reached the WJ collar")
        at_point = []
        for e in lower_elements:
            lo, hi, bottom, top = e["bounds"]
            if lo-1e-12 <= x <= hi+1e-12 and bottom-1e-12 <= r <= top+1e-12:
                a = oracle.sample(e, np.array([x]), np.array([r]))
                grad = a["G"][0]
                values = np.array([a["U"][0], a["W"][0], grad[0,0], grad[0,2],
                                   grad[2,0], grad[2,2], a["Pi"][0]/ref.tau]).tolist()
                k, interface_k = float(a["k"][0]), float(a["K"][0])
                local_aux = float(ref.aux(np.array([interface_k]), np.array([r]))[0])
                covariance_margin = 2/3-local_aux*abs(grad[0,2]+grad[2,0])/k
                require(abs(covariance_margin-a["covariance"][0]) < 5e-13,
                        "auxiliary-only covariance reconstruction disagrees")
                at_point.append(dict(geometry=e["geometry"], bounds=e["bounds"], fields=values,
                                     D=float(a["div"][0]), k=k, epsilon=float(a["E"][0]),
                                     interface_K=interface_k, covariance_over_k=float(covariance_margin),
                                     transverse_principal_minimum_over_mu=ref.nu/(ref.nu+local_aux)))
        require(at_point, "annulus sample without an original physical owner")
        single = at_point[0] if len(at_point) == 1 else None
        if single:
            fields[index] = single["fields"]
        owners.append(single)
        observations.append(dict(index=index, x=x, r=r, beyond_terminal_widths=bool(s>width[0] and r>width[1]),
                                 one_sided=at_point))
    return dict(state=identity(path), source=data["source_identity"], K=K, width=width, fields=fields, owners=owners,
                report=dict(vertex_owner=vertex["geometry"], terminal_dx_dr=width.tolist(),
                            local_K=K, pencil=mode.audit(), samples=observations))


def smooth_background(ref, K, s, r, base, degree):
    # Both backgrounds solve no slip and start Pi at the smooth-part value
    # from the fixed normal traction. Velocity corrections use solenoidal
    # streamfunction monomials r^j*s^i with j>=2.
    values = np.zeros((len(s), len(FIELDS)))
    values[:, 6] = -ref.outlet(np.zeros(1))[0,0]/ref.tau
    if base == "exact_reference":
        inlet = ref.inlet(r)
        values[:, 0], values[:, 3] = inlet["U"], inlet["G"]
    else:
        a = ref.kappa*math.sqrt(K/ref.qe)
        mu0 = ref.nu+float(ref.aux(np.array([K]), np.zeros(1))[0])
        values[:, 0] = ref.reference_shear/a*np.log1p(a*r/mu0)
        values[:, 3] = ref.reference_shear/(mu0+a*r)

    L = ANCHOR_RADIUS
    columns = []
    for total in range(2, degree+1):
        for j in range(2, total+1):
            i = total-j
            derivative = {}
            for ds, dr in ((1,0),(0,1),(1,1),(0,2),(2,0)):
                derivative[ds,dr] = np.zeros_like(s)
                if i >= ds and j >= dr:
                    coefficient = math.prod(i-t for t in range(ds))*math.prod(j-t for t in range(dr))
                    derivative[ds,dr] = coefficient*L**(2-ds-dr)*(s/L)**(i-ds)*(r/L)**(j-dr)
            drs = derivative[1,1]
            columns.append(np.column_stack((-derivative[0,1],-derivative[1,0],
                                            drs,-derivative[0,2],derivative[2,0],-drs,np.zeros_like(s))))
    for total in range(1, degree-1):
        for i in range(total+1):
            pressure_column = np.zeros_like(values)
            pressure_column[:,6] = (s/L)**i*(r/L)**(total-i)
            columns.append(pressure_column)
    return values, np.stack(columns, axis=2)


def owner_weights(owners, mask):
    ids = [e["geometry"] for e, keep in zip(owners, mask) if keep]
    unique, membership, counts = np.unique(ids, return_inverse=True, return_counts=True)
    return 1/np.sqrt(counts[membership]), set(unique.tolist())


def norms(values, weights):
    return dict(RMS=np.sqrt(np.average(values**2, axis=0, weights=weights**2)).tolist(),
                MAX=np.max(np.abs(values), axis=0).tolist(),
                budget_RMS=float(np.sqrt(np.mean(np.average((values/BUDGETS)**2, axis=0, weights=weights**2)))))


def fit_window(fields, envelope, owners, mask, parts, background, columns, singular):
    training = mask & parts["training"]
    groups = {name: mask & flag for name, flag in parts.items() if name != "training"}
    if not np.any(training) or not all(np.any(group) for group in groups.values()):
        return dict(status="INSUFFICIENT_SAMPLES")
    weights, train_owners = owner_weights(owners, training)
    columns_x = {tuple(e["bounds"][:2]) for e, keep in zip(owners, training) if keep}
    columns_r = {tuple(e["bounds"][2:]) for e, keep in zip(owners, training) if keep}
    weight = weights[:,None]/BUDGETS
    A = (columns[training]*weight[:,:,None]).reshape(-1, columns.shape[2])
    b = (singular[training]*weight).ravel()
    y = ((fields[training]-background[training])*weight).ravel()
    coeff, _, rank, _ = np.linalg.lstsq(A, y, rcond=1e-12)
    q = b-A@np.linalg.lstsq(A, b, rcond=1e-12)[0]
    energy = float(q@q)
    independence = float(np.linalg.norm(q)/np.linalg.norm(b))
    if rank != A.shape[1] or independence <= 1e-6:
        return dict(status="ILL_CONDITIONED", rank=int(rank), mode_independence=independence)
    amp = float(q@y/energy)
    nuisance = np.linalg.lstsq(A, y-amp*b, rcond=1e-12)[0]
    regular = background+columns@coeff
    fitted = background+columns@nuisance+amp*singular
    delta = float(np.abs(q)@(envelope[training]*weight).ravel()/energy) if envelope is not None else None
    observed_divergence = fields[:,2]+fields[:,5]
    result = dict(status="FIT", amplitude=amp, sampled_mesh_amplitude_bound=delta,
                  mode_independence=independence, training_x_columns=len(columns_x), training_r_rows=len(columns_r),
                  training_owner_count=len(train_owners), training_points=int(np.count_nonzero(training)), heldout={})
    for name, group in dict(training=training, **groups).items():
        test_weights, test_owners = owner_weights(owners, group)
        metrics = dict(points=int(np.count_nonzero(group)), owners=len(test_owners),
                       owners_unseen_in_training=len(test_owners-train_owners),
                       regular_only=norms(fields[group]-regular[group], test_weights),
                       singular_and_regular=norms(fields[group]-fitted[group], test_weights),
                       D_RMS=float(np.sqrt(np.average(observed_divergence[group]**2, weights=test_weights**2))),
                       D_MAX=float(np.max(np.abs(observed_divergence[group]))))
        if envelope is not None:
            metrics["observed_intermesh_difference"] = norms(envelope[group], test_weights)
        result["heldout"][name] = metrics
    return result


def examine(cases, ref, mode_fields, points, radii, angles, resolved, sampling):
    s, r = points.T
    mask = np.all([np.isfinite(case["fields"]).all(axis=1) for case in cases], axis=0)
    if resolved:
        mask &= np.all([np.all(points > case["width"], axis=1) for case in cases], axis=0)
    hold_r = np.isin(radii,sampling["heldout_radii"])
    hold_angle = np.isin(angles,sampling["heldout_angles_deg"])
    parts = dict(training=~(hold_r|hold_angle), ray=hold_angle & ~hold_r,
                 radius=hold_r & ~hold_angle, ray_and_radius=hold_angle & hold_r)
    result = dict(resolved_sampling=resolved, common_points=int(np.count_nonzero(mask)), fits={})
    if not np.any(mask):
        result["status"] = "INSUFFICIENT_COMMON_ANNULI"
        return result
    fit_fields = np.array([case["fields"] for case in cases])
    envelope = np.max(np.abs(np.diff(fit_fields, axis=0)), axis=0) if len(cases) > 1 else None
    owners = cases[0]["owners"]  # Same coarse-owner weights at each compared level.
    for base, degree in BACKGROUND_MODELS:
        background, columns = smooth_background(ref, cases[-1]["K"], s, r, base, degree)
        for lo, hi in sampling["annuli_code"]:
            window = mask & (radii >= lo) & (radii <= hi)
            record = dict(background=base, degree=degree, annulus_code=[lo, hi],
                          levels=[fit_window(values, envelope, owners, window, parts,
                                             background, columns, mode_fields) for values in fit_fields])
            result["fits"][f"{base}-{degree}/{lo:g}-{hi:g}"] = record
    result["status"] = "DIAGNOSTIC_FITS"
    return result


def classify(last_three, cases):
    if len(cases) < 3 or not last_three["fits"]:
        return "INCONCLUSIVE", "Three common resolved corner grids are unavailable."
    widths = np.array([case["width"] for case in cases[-3:]])
    changes = np.diff(widths, axis=0)
    contracted = np.all(changes <= 1e-12) and np.all(np.any(changes < -1e-12, axis=0))
    if not contracted or not np.all(np.any(changes < -1e-12, axis=1)):
        return "INCONCLUSIVE", "The last three grids lack nested decreases in both terminal directions across the sequence."
    fits = [fit for window in last_three["fits"].values() for fit in window["levels"]]
    if not all(fit["status"] == "FIT" and fit["training_x_columns"] >= 2 and fit["training_r_rows"] >= 2
               and min(fit["heldout"][group]["owners_unseen_in_training"] for group in ("ray","radius")) >= 2
               for fit in fits):
        return "INCONCLUSIVE", "Some annuli lack independently held-out original owners or a resolved fit."
    amplitudes = np.array([fit["amplitude"] for fit in fits])
    mesh_floor = max(fit["sampled_mesh_amplitude_bound"] for fit in fits)
    if (np.min(np.abs(amplitudes)) <= 3*mesh_floor or np.min(amplitudes)*np.max(amplitudes) <= 0
            or np.ptp(amplitudes) > .5*np.max(np.abs(amplitudes))):
        return "INCONCLUSIVE", "Mesh uncertainty or annulus/background sensitivity leaves the amplitude unresolved."
    if not all(fit["heldout"][group]["singular_and_regular"]["budget_RMS"] <
               .8*fit["heldout"][group]["regular_only"]["budget_RMS"]
               for fit in fits for group in ("ray","radius")):
        return "INCONCLUSIVE", "The angular singular model lacks a repeatable held-out improvement."
    return ("REPEATABLE_NONZERO_MODE_ON_RESOLVED_SAMPLES",
            "Common resolved fits exceed the observed numerical envelope and improve independent held-out predictions.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", required=True, type=Path)
    parser.add_argument("--native", required=True, action="append", type=Path,
                        help="Registered steady state, with coarser grids first.")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--seconds", required=True, type=float)
    parser.add_argument("--near",action="store_true",help="Supplement with shorter annuli to limit coefficient variation.")
    args = parser.parse_args()
    require(sys.platform.startswith("linux"), "coupled numerical evidence is collected on Linux")
    require(math.isfinite(args.seconds) and args.seconds > 0., "explicit finite wall cap required")
    deadline = time.monotonic()+args.seconds
    ref = load_reference(args.table)
    require(ref.exact_inner, "corner diagnosis requires the immutable exact-inner reference")
    sampling = NEAR if args.near else dict(radii_code=RADII,angles_deg=ANGLES,
        heldout_radii=HOLD_RADII,heldout_angles_deg=HOLD_ANGLES,annuli_code=ANNULI)
    radii, angles = np.meshgrid(sampling["radii_code"],sampling["angles_deg"],indexing="ij")
    radii, angles = radii.ravel(), angles.ravel()
    phase = np.deg2rad(angles)
    points = radii[:,None]*np.column_stack((np.cos(phase), np.sin(phase)))
    cases = [measure(path, ref, deadline, points) for path in args.native]
    require(len({case["source"] for case in cases}) == 1, "coupled angular fits require one source identity")
    require(len({case["state"]["sha256"] for case in cases}) == len(cases), "duplicate raw coupled state")
    K = np.array([cases[-1]["K"]])
    mode = Wedge(ref.nu, ref.nu+float(ref.aux(K, np.zeros(1))[0]))
    shape = mode.fields(*points.T)
    shape[:,6] /= ref.tau
    studies = dict(all_levels_unresolved_control=examine(cases,ref,shape,points,radii,angles,False,sampling),
                   latest_individual=examine(cases[-1:],ref,shape,points,radii,angles,True,sampling),
                   latest_pair=examine(cases[-2:],ref,shape,points,radii,angles,True,sampling),
                   latest_three=examine(cases[-3:],ref,shape,points,radii,angles,True,sampling))
    status, reason = classify(studies["latest_three"], cases)
    for case in cases:
        case["report"].update(state=case["state"], source_identity=case["source"])
    require(time.monotonic() < deadline, "wedge fit exceeded its cooperative cap")
    atomic(args.out, dict(schema="WJH27Task2RSAngularWedgeV1", table=identity(args.table),
        status=status, reason=reason, fitted_shape_pencil=mode.audit(), field_order=list(FIELDS),
        field_scales=[1.]*6+[ref.tau], fit_budget_weights=BUDGETS.tolist(),
        pressure_treatment="Absolute physical pressure in tau0 units, with fixed smooth-part Pi(0)=-t_x(0). No fitted constant.",
        sampling=dict(**sampling, supplemental_near_profile=args.near,
                      largest_relative_shear_coefficient_change=float(
                          np.max(abs(ref.nu+ref.aux(K,points[:,1])-mode.mu))/mode.mu),
                      owner_policy="Record all ties. Fit unambiguous original owners only and weight each coarse owner equally.",
                      resolved_policy="Both s and r exceed their terminal owner widths on every compared grid."),
        amplitude_uncertainty="Observed common-point intermesh envelope propagated through the nuisance-orthogonal mode. Window and background spread are separate sensitivities.",
        divergence_scope="The profile is solenoidal. Raw Ux and Wr retain the measured D, which is orthogonal to model columns with equal gradient weights.",
        admissibility_scope="Sampled auxiliary-only covariance and principal values on the rays. Divergent shear with finite k would threaten covariance. Continuous behavior is unresolved.",
        levels=[case["report"] for case in cases], studies=studies))
    print(status, reason, flush=True)


if __name__ == "__main__":
    main()
