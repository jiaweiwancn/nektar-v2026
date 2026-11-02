#!/usr/bin/env python3
"""Linux verification of the revised native flat multiregion wall stage.

The residual oracle uses dimensional indexed CLS, Jacobi bases and a separate
physical quadrature. Native maps are structural inputs, not a second stress or
source oracle. Raw matrices/fields stay in the Linux run directory.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile
import time
import re
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
THREADS = {k: "1" for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS")}
os.environ.update(THREADS)
import numpy as np
import scipy
from scipy.special import eval_jacobi
from validate_cls_epawf_m2 import Profiles
from validate_cls_epawf_m1 import indexed
from cls_epawf_m3_cases import parallel_session, evolving_session, set_info, set_parameter

checks = 0
errors = {}


def require(condition, message):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(message)


def near(actual, expected, name, scale=1., tolerance=1e-10):
    a, b = np.asarray(actual), np.asarray(expected)
    require(a.shape == b.shape and np.all(np.isfinite(a)) and np.all(np.isfinite(b)), "Invalid "+name)
    error = float(np.max(np.abs(a-b)))/max(scale, float(np.max(np.abs(b))))
    errors[name] = max(error, errors.get(name, 0.))
    require(error <= tolerance, f"{name}: scaled error {error} > {tolerance}")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def basis(degree, x):
    x = np.atleast_1d(np.asarray(x, dtype=float))
    p = [(1-x)/2, (1+x)/2]
    d = [np.full_like(x, -.5), np.full_like(x, .5)]
    for j in range(2, degree+1):
        n = j-2
        jac = eval_jacobi(n, 1, 1, x)
        dp = 0 if n == 0 else (n+3)/2*eval_jacobi(n-1, 2, 2, x)
        p.append((1-x*x)/4*jac)
        d.append(-x/2*jac+(1-x*x)/4*dp)
    return np.array(p).T, np.array(d).T


class Space:
    def __init__(self, record):
        self.elements = record["mesh"]
        self.size = record["global_coefficients"]
        self.by_id = {e["geometry"]: e for e in self.elements}

    def element(self, x, r, side=0):
        possible = []
        for e in self.elements:
            a, b, c, d = e["bounds"]
            if a-1e-12 <= x <= b+1e-12 and c-1e-12 <= r <= d+1e-12:
                if side > 0 and r >= d-1e-12 or side < 0 and r <= c+1e-12:
                    continue
                possible.append(e)
        require(bool(possible), "Independent physical sampler missed a point")
        return min(possible, key=lambda e: e["geometry"])

    def rows(self, e, x, r):
        a, b, c, d = e["bounds"]
        px, dx = basis(e["degree"], 2*(np.asarray(x)-a)/(b-a)-1)
        pr, dr = basis(e["degree"], 2*(np.asarray(r)-c)/(d-c)-1)
        values = (pr[:, :, None]*px[:, None, :]).reshape(len(px), -1)
        gx = (pr[:, :, None]*dx[:, None, :]).reshape(len(px), -1)*2/(b-a)
        gr = (dr[:, :, None]*px[:, None, :]).reshape(len(px), -1)*2/(d-c)
        signs = np.array(e["sign"])
        return np.array(e["ids"], dtype=int), values*signs, gx*signs, gr*signs

    def evaluate(self, state, offset, x, r, side=0):
        e = self.element(x, r, side)
        ids, p, dx, dr = self.rows(e, [x], [r])
        coefficients = state[offset+ids]
        return float(p@coefficients), np.array([float(dx@coefficients), float(dr@coefficients)])


def check_layout(raw):
    summaries = []
    for record in raw["layout"]:
        space = Space(record)
        x, r, value, dx, dr = record["probe"]
        near([value, dx, dr], [1+.2*x+.03*r*r, .2, .06*r], "native_polynomial_evaluation")
        expected = np.zeros(space.size)
        z, w = np.polynomial.legendre.leggauss(12)
        for e in space.elements:
            a, b, c, d = e["bounds"]
            X, R = np.meshgrid(a+(z+1)*(b-a)/2, c+(z+1)*(d-c)/2)
            X, R = X.ravel(), R.ravel()
            W = np.outer(w, w).ravel()*(b-a)*(d-c)/4
            ids, p, _, _ = space.rows(e, X, R)
            np.add.at(expected, ids, p.T@(W*(1+.2*X+.03*R*R)))
        near(record["native_load"], expected, "native_subdomain_CG_assembly")
        defect = np.zeros(space.size)
        coefficients = np.array(record["nonpolynomial_projection"])
        z, w = np.polynomial.legendre.leggauss(24)
        for e in space.elements:
            a, b, c, d = e["bounds"]
            X, R = np.meshgrid(a+(z+1)*(b-a)/2, c+(z+1)*(d-c)/2)
            X, R = X.ravel(), R.ravel()
            W = np.outer(w, w).ravel()*(b-a)*(d-c)/4
            ids, p, _, _ = space.rows(e, X, R)
            np.add.at(defect, ids, p.T@(W*(p@coefficients[ids]-np.sin(.4*X+.2*R)-.01*R**3)))
        near(defect, np.zeros(space.size), "physical_L2_projection_all_mode_defect")
        summaries.append({k: record[k] for k in ("elements", "global_coefficients", "local_coefficients")})
    require(summaries[0]["elements"] == summaries[3]["elements"], "Full-domain k space differs")
    require(summaries[4]["elements"] == summaries[2]["elements"], "Epsilon is not outer-only")
    return summaries


def check_pressure(raw):
    pressure = raw["pressure"]
    B, C = np.array(pressure["divergence_matrix"]), np.array(raw["trace_constraint"])
    singular = np.linalg.svd(B, compute_uv=False)
    rank = int(np.sum(singular > 1e-11*singular[0]))
    require(rank == B.shape[0]-1, "Mixed pressure has extra null modes")
    require(np.linalg.matrix_rank(C, tol=1e-11*np.linalg.norm(C, 2)) == C.shape[0], "Interface constraints are dependent")
    A, rhs, solution = map(np.array, (pressure["system"], pressure["rhs"], pressure["solution"]))
    near(A@solution, rhs, "augmented_native_pressure_solve")
    require(np.linalg.matrix_rank(A, tol=1e-11*np.linalg.norm(A, 2)) == len(rhs), "Global pressure gauge did not remove exactly one null mode")
    k = raw["parameters"][5]
    near(pressure["metrics"][2], -11*k/75, "mechanical_reference_pressure_jump", k)
    near(pressure["metrics"][:2], [0, 0], "pressure_probe_velocity_and_mean")
    near(pressure["metrics"][3], 0, "full_normal_traction_matching", k)
    z, w = np.polynomial.legendre.leggauss(16)
    H = raw["parameters"][1]
    for side, record in enumerate((raw["layout"][1], raw["layout"][2])):
        space = Space(record)
        expected = np.zeros(space.size)
        for e in space.elements:
            a, b, c, d = e["bounds"]
            if not (abs((c if side else d)-H) < 1e-12):
                continue
            X, R = a+(z+1)*(b-a)/2, np.full_like(z, H)
            ids, p, _, _ = space.rows(e, X, R)
            np.add.at(expected, ids, p.T@(w*(b-a)/2*(.17+.03*X)*(1 if side == 0 else -1)))
        near(raw["shared_loads"][side], expected, "physical_common_flux_load_"+str(side))
    near(raw["merged_shared_load"], np.zeros(raw["layout"][0]["global_coefficients"]), "opposite_common_flux_dual_loads")
    near(raw["strong_Dirichlet_ground_load_control"], np.zeros(raw["layout"][0]["global_coefficients"]), "strong_no_slip_ground_load_control")
    Z, UZ, KZ, EZ, nt, dU, dE, gradient_flux = raw["top_reference"]
    us = math.sqrt(k/raw["parameters"][6])
    near([nt*dU, dE, gradient_flux], [us**2, -us**3/(raw["parameters"][7]*Z**2), -us**4/(1.3*Z)], "compatible_top_reference_flux_factors")
    for sampled, expected in zip(raw["top_reference_trace_samples"], (UZ, KZ, EZ)):
        near(sampled, [expected]*3, "native_compatible_top_Dirichlet_values")
    return {"pressure_rows": B.shape[0], "pressure_rank": rank,
            "pressure_singular_values": singular.tolist(), "interface_rows": C.shape[0],
            "interface_rank": C.shape[0], "augmented_rows": len(rhs),
            "metrics": pressure["metrics"]}


def check_condensation(raw):
    case = next(c for c in raw["generic_cases"] if "condensation" in c)
    c = case["condensation"]
    J, r = np.array(c["jacobian"]), np.array(case["residual"])
    I, R = np.array(c["inner"], dtype=int), np.array(c["retained"], dtype=int)
    require(len(set(I) & set(R)) == 0 and sorted(np.r_[I, R]) == list(range(len(r))), "Invalid condensation partition")
    require(np.all(I < raw["offsets"][2]), "Pressure or inner k was silently condensed")
    require(all(i in R for i in range(raw["offsets"][4], len(r))), "Relative pressure/gauge modes were discarded")
    ii, ir, ri, rr = J[np.ix_(I, I)], J[np.ix_(I, R)], J[np.ix_(R, I)], J[np.ix_(R, R)]
    S = rr-ri@np.linalg.solve(ii, ir)
    rhs = r[R]-ri@np.linalg.solve(ii, r[I])
    near(c["schur"], S, "independent_Schur_matrix")
    near(c["reduced_residual"], rhs, "independent_condensed_residual")
    expected = np.linalg.solve(J, -r)
    near(c["full_step"], expected, "independent_full_Newton_step")
    near(c["condensed_step"], expected, "condensed_recovered_Newton_step")
    direction, lifted = np.array(c["retained_direction"]), np.array(c["lifted_direction"])
    near(lifted[R], direction, "retained_action_injection")
    near(lifted[I], -np.linalg.solve(ii, ir@direction), "implicit_inner_velocity_action")
    action = J@lifted
    near(c["full_action"], action, "uncondensed_directional_action")
    near(action[I], np.zeros(len(I)), "eliminated_directional_residual")
    near(c["reduced_action"], action[R], "condensed_directional_action")
    differences = []
    for sample in c["finite_difference_actions"]:
        error = float(np.max(np.abs(np.array(sample["action"])-action)))/max(1., float(np.max(np.abs(action))))
        differences.append({"step": sample["step"], "scaled_error": error})
    require(differences[-1]["scaled_error"] < 1e-6 and
            differences[-1]["scaled_error"] < .1*differences[0]["scaled_error"]+1e-8,
            "Physical residual action lacks finite-difference convergence")
    return {"full_rows": len(r), "eliminated_velocity_rows": len(I), "retained_rows": len(R),
            "all_pressure_modes_retained": True, "physical_action_difference_sequence": differences,
            "condition_number": float(np.linalg.cond(J))}


def wall_integrals(profile, r):
    r = np.asarray(r)
    lo = max(0., float(profile.zv))
    z = np.minimum(r, lo)
    a, d, nu, h = map(float, (profile.a, profile.d, profile.nu, profile.h))
    j0 = z/nu+np.where(r > z, np.log1p(a*(r-z)/(d+a*z))/a, 0.)
    lower = np.maximum(z, h)
    span = np.maximum(0, r-lower)
    j1 = np.maximum(z-h, 0)**2/(2*nu)+span/a-(d+a*h)/a**2*np.log1p(a*span/(d+a*lower))
    return j0, j1


def check_parallel(raw):
    field = Space(raw["layout"][1])
    nu, H, h, rK, rm, K, q, kappa = raw["parameters"]
    profile = Profiles({"geometry": [nu, H, h, rK, rm], "metric": [0, 0], "k": K, "q": q, "kappa": kappa, "chi": 1})
    intervals = sorted({(e["bounds"][2], e["bounds"][3], e["degree"]) for e in field.elements})
    # Independent one-dimensional finite-element assembly. No native restriction
    # or native radial matrix is reused in this solve.
    mapping, local = {}, []
    for lower, upper, degree in intervals:
        ids = []
        for j in range(degree+1):
            key = (lower, lower, 0) if j == 0 else (upper, upper, 0) if j == 1 else (lower, upper, j)
            if key not in mapping:
                mapping[key] = len(mapping)
            ids.append(mapping[key])
        local.append(ids)
    n = len(mapping)
    A, unit_force = np.zeros((n, n)), np.zeros(n)
    z, w = np.polynomial.legendre.leggauss(48)
    for (lower, upper, degree), ids in zip(intervals, local):
        cuts = sorted({lower, upper, *[b for b in (h, float(profile.zv)) if lower < b < upper]})
        for a, b in zip(cuts[:-1], cuts[1:]):
            r = a+(z+1)*(b-a)/2
            p, dp = basis(degree, 2*(r-lower)/(upper-lower)-1)
            dp *= 2/(upper-lower)
            W = w*(b-a)/2
            viscosity = nu+float(profile.a)*np.maximum(r-float(profile.zv), 0)
            A[np.ix_(ids, ids)] += dp.T@((W*viscosity)[:, None]*dp)
            unit_force[ids] += p.T@(W*(r > h))
    ground, top = mapping[(0., 0., 0)], mapping[(H, H, 0)]
    free = [i for i in range(n) if i not in (ground, top)]
    summaries = []
    for case in raw["parallel_cases"]:
        metric = case["metrics"]
        UH, F = metric[:2]
        coefficients = np.zeros(n)
        coefficients[top] = UH
        coefficients[free] = np.linalg.solve(A[np.ix_(free, free)], -F*unit_force[free]-A[free, top]*UH)
        reactions = A@coefficients+F*unit_force
        near(metric[2:4], [-reactions[ground], reactions[top]], "independent_parallel_weak_traction", K/q)
        traction = UH/float(profile.J0)-F*float(profile.J1/profile.J0)
        near(metric[4:6], [traction, traction+F*(H-h)], "parallel_analytical_forcing_response", K/q)
        samples = np.array(case["samples"])
        expected, derivative = [], []
        for r in samples[:, 0]:
            pos = next(i for i, (a, b, _) in enumerate(intervals) if a-1e-13 <= r <= b+1e-13)
            a, b, degree = intervals[pos]
            p, dp = basis(degree, [2*(r-a)/(b-a)-1])
            expected.append(float(p@coefficients[local[pos]]))
            derivative.append(float(dp@coefficients[local[pos]])*2/(b-a))
        near(samples[:, 1], expected, "native_parallel_field", abs(UH))
        near(samples[:, 2], derivative, "native_parallel_physical_derivative", abs(UH)/H)
        j0, j1 = wall_integrals(profile, samples[:, 0])
        near(samples[:, 3], traction*j0+F*j1, "parallel_original_profile", abs(UH))
        native = np.array(case["coefficients"])
        for x in (.13, .47, .91):
            for r in samples[::10, 0]:
                value, gradient = field.evaluate(native, 0, x, r)
                near(value, expected[int(round(r/H*100))], "native_constant_tangential_restriction", abs(UH))
                near(gradient[0], 0, "native_parallel_tangential_derivative", abs(UH))
        near(metric[8:10], [0, 0], "parallel_ground_and_matching_constraints", abs(UH))
        near(metric[12], 0, "parallel_weak_energy_budget", abs(UH*traction)+1)
        summaries.append({"forcing": F, "radial_unknowns": case["radial_unknowns"],
                          "velocity_RMS_error": metric[6], "gradient_RMS_error": metric[7],
                          "weak_traction_error": max(abs(metric[2]-metric[4]), abs(metric[3]-metric[5])),
                          "weak_ground_traction": metric[2], "weak_matching_traction": metric[3],
                          "point_ground_traction": metric[10], "point_matching_traction": metric[11],
                          "energy_residual": metric[12]})
    return summaries


def residual_oracle(raw, case):
    # Layout records 0,3,4,5,6 are full velocity, full k, outer epsilon,
    # inner pressure, outer pressure. Maps come from native AssemblyMapCG.
    spaces = [Space(raw["layout"][i]) for i in (0, 3, 4, 5, 6)]
    v, ks, es, pw, po = spaces
    offsets = raw["offsets"]
    s, history = np.array(case["lifted_state"]), np.array(case["history"])
    result = np.zeros(raw["stage_unknowns"])
    nu, H, h, rK, rm, K0, qstar, kappa = raw["parameters"]
    chi, alpha = case["chi"], case["alpha"]
    columns = sorted({tuple(e["bounds"][:2]) for e in pw.elements})
    model = lambda K: Profiles({"geometry": [nu, H, h, rK, rm], "metric": [0, 0],
                               "k": K, "q": qstar, "kappa": kappa, "chi": chi})
    for x, sampled, forcing in case["forcing_samples"]:
        uH = (.1+.03*x)*H
        near(sampled, K0+.02*x, "fixed_physical_K_sampler", K0)
        near(forcing, 2*uH*.03*H+.04-2*.02/3-(.04*x+.01*H), "one_sided_outer_conservative_forcing")
    for _, value in case["potential_cancellation"]:
        near(value, 0, "mechanical_pressure_and_body_potential_cancellation")
    def point_fields(space, e, x, r, block):
        ids, p, dx, dr = space.rows(e, x, r)
        coefficients = s[offsets[block]+ids]
        return ids, p, dx, dr, p@coefficients, dx@coefficients, dr@coefficients
    def add(block, ids, p, dx, dr, W, source, fx, fr):
        np.add.at(result, offsets[block]+ids, p.T@(W*source)+dx.T@(W*fx)+dr.T@(W*fr))
    def integrate(geometry, x, r, W, inner, profile=None, forcing=0):
        ve = v.by_id[geometry]
        ids, p, dx, dr, U, ux, ur = point_fields(v, ve, x, r, 0)
        _, _, _, _, V, vx, vr = point_fields(v, ve, x, r, 1)
        ki, kp, kx, kr, K, Kx, Kr = point_fields(ks, ks.by_id[geometry], x, r, 2)
        ps = pw if inner else po
        pi, pp, px, pr, pressure, _, _ = point_fields(ps, ps.by_id[geometry], x, r, 4 if inner else 5)
        gradient = np.zeros((len(x), 3, 3))
        gradient[:, 0, 0], gradient[:, 0, 2] = ux, ur
        gradient[:, 2, 0], gradient[:, 2, 2] = vx, vr
        strain = gradient+gradient.transpose(0, 2, 1)
        if inner:
            nt = float(profile.a)*np.maximum(r-float(profile.zv), 0)
            reynolds = 2*K[:, None, None]/3*np.eye(3)-nt[:, None, None]*strain
            production = -np.sum(reynolds*gradient, axis=(1, 2))
            stress = (nu+nt)[:, None, None]*strain-pressure[:, None, None]*np.eye(3)
            reference = np.asarray(profile.fields(r)[1], dtype=float)
            original = K0*0+float(profile.k)**1.5/(2.55*np.maximum(r, float(profile.ze)))
            sk = production-(1-chi)*original-chi*reference+.03+.01*x+.02*r
            dk = nu+nt
            bx, br = .01+.02*x+(r <= h)*forcing, -.02+.03*r
        else:
            ei, ep, ex, er, E, Ex, Er = point_fields(es, es.by_id[geometry], x, r, 3)
            closures = [indexed(g, k, e, "physical") for g, k, e in zip(gradient, K, E)]
            reynolds = np.array([c["stress"].real for c in closures])
            nt = np.array([c["nut"].real for c in closures])
            production = -np.sum(reynolds*gradient, axis=(1, 2))
            stress = -reynolds+(2*K/3-pressure)[:, None, None]*np.eye(3)
            dk, de = nt, nt/1.3
            sk = production-E-.02+.01*x+.03*r
            ratio = K**1.5/(2.5*E*r)
            se = 1.44*E/K*production-1.92*E*E/K+.83*E*E/K*np.maximum(ratio-1, 0)*ratio**2+.01*(1+x+r)
            bx, br = .04*x+.01*r, -.03+.02*r
        store_u = alpha*(U-p@history[offsets[0]+ids])
        store_v = alpha*(V-p@history[offsets[1]+ids])
        store_k = alpha*(K-kp@history[offsets[2]+ki])
        add(0, ids, p, dx, dr, W, store_u-bx, stress[:, 0, 0]-U*U, stress[:, 0, 2]-U*V)
        add(1, ids, p, dx, dr, W, store_v-br, stress[:, 2, 0]-V*U, stress[:, 2, 2]-V*V)
        add(2, ki, kp, kx, kr, W, store_k-sk, dk*Kx-U*K, dk*Kr-V*K)
        add(4 if inner else 5, pi, pp, px, pr, W, ux+vr, np.zeros_like(x), np.zeros_like(x))
        if not inner:
            store_e = alpha*(E-ep@history[offsets[3]+ei])
            add(3, ei, ep, ex, er, W, store_e-se, de*Ex-U*E, de*Er-V*E)
    z, w = np.polynomial.legendre.leggauss(24)
    for left, right in columns:
        for x, wx in zip(left+(z+1)*(right-left)/2, w*(right-left)/2):
            K, _ = ks.evaluate(s, offsets[2], x, rK, -1 if rK == H else 0)
            profile = model(K)
            u, gu = v.evaluate(s, offsets[0], x, H, 1)
            _, gp = po.evaluate(s, offsets[5], x, H, 1)
            _, gk = ks.evaluate(s, offsets[2], x, H, 1)
            forcing = 2*u*gu[0]+gp[0]-2*gk[0]/3-(.04*x+.01*H)
            breaks = {0., H, h, float(profile.ze), float(profile.rc)}
            if profile.zv > 0:
                breaks.add(float(profile.zv))
            start = min(float(profile.rc), nu/float(profile.a), float(profile.ze))
            while start < H:
                breaks.add(start)
                start *= 3
            breaks.update(e["bounds"][2] for e in pw.elements)
            breaks.update(e["bounds"][3] for e in pw.elements)
            points = sorted(r for r in breaks if 0 <= r <= H)
            for lower, upper in zip(points[:-1], points[1:]):
                R = lower+(z+1)*(upper-lower)/2
                X = np.full_like(R, x)
                geometry = pw.element(x, (lower+upper)/2)["geometry"]
                integrate(geometry, X, R, wx*w*(upper-lower)/2, True, profile, forcing)
    for e in po.elements:
        left, right, lower, upper = e["bounds"]
        X, R = np.meshgrid(left+(z+1)*(right-left)/2, lower+(z+1)*(upper-lower)/2)
        W = np.outer(w, w).ravel()*(right-left)*(upper-lower)/4
        integrate(e["geometry"], X.ravel(), R.ravel(), W, False)
    volume = (max(e["bounds"][1] for e in v.elements)-min(e["bounds"][0] for e in v.elements))*max(e["bounds"][3] for e in v.elements)
    for block, space in ((4, pw), (5, po)):
        load = np.zeros(space.size)
        for e in space.elements:
            a, b, c, d = e["bounds"]
            X, R = np.meshgrid(a+(z+1)*(b-a)/2, c+(z+1)*(d-c)/2)
            ids, p, _, _ = space.rows(e, X.ravel(), R.ravel())
            np.add.at(load, ids, p.T@(np.outer(w, w).ravel()*(b-a)*(d-c)/4/volume))
        result[offsets[block]:offsets[block]+space.size] += s[offsets[6]]*load
        result[offsets[6]] += s[offsets[block]:offsets[block]+space.size]@load
    for row, value in case["constraints"]:
        result[int(row)] = case["state"][int(row)]-value
    for block in range(6):
        end = offsets[block+1]
        near(case["residual"][offsets[block]:end], result[offsets[block]:end], "independent_residual_block_"+str(block))
    near(case["residual"][offsets[6]], result[offsets[6]], "one_global_pressure_gauge")
    return {"chi": chi, "unknowns": len(result), "max_residual": float(np.max(np.abs(result)))}


def manufactured_values(x, r, time, parameters, derivative=False):
    _, H, _, rK, _, K, q, kappa = parameters
    Z = 4.
    E = (K/q)**1.5/(kappa*H)
    U, W = .2+.03*math.sin(time), .01+.005*math.cos(time)
    Ut, Wt = .03*math.cos(time), -.005*math.sin(time)
    s, so = r/H, (r-H)/(Z-H)
    gate, dg = 6*s*s-8*s**3+3*s**4, (12*s-24*s*s+12*s**3)/H
    a = 2*s-s*s
    g = 4*x*(1-x)*r*r*(r-rK)**2*(r-H)**2/H**6
    if derivative:
        ui, wi = -Wt*x*dg+Ut*a, Wt*gate
        ki, ko, eo = .01*math.cos(time)*g, -.005*math.sin(time)*so*so, .001*math.cos(time)*so*so
        return np.array([ui if r <= H else Ut, wi if r <= H else Wt,
                         ki if r <= H else ko, eo if r >= H else 0., 0., 0.])
    ui, wi = -W*x*dg+U*a, W*gate
    ki, ko = K+(.04+.01*math.sin(time))*g, K+(.02+.005*math.cos(time))*so*so
    eo = E+(.003+.001*math.sin(time))*so*so
    px, pzi, pzo = .02*(1+.1*math.sin(time)), .03*math.cos(time), .04*math.sin(time)
    mean = (pzi*H**3+pzo*(Z-H)**3)/(3*Z)
    pi, po = px*(x-.5)+pzi*(r-H)**2-mean, px*(x-.5)+pzo*(r-H)**2-mean
    return np.array([ui if r <= H else U, wi if r <= H else W,
                     ki if r <= H else ko, eo if r >= H else 0., pi, po])


def scaled_residual_norms(report, metadata):
    parameters = metadata["parameters"]
    H, K, q, kappa = parameters[1], parameters[5], parameters[6], parameters[7]
    us, scales = math.sqrt(K/q), []
    time_scale = H/us
    for i, field_scale in zip((0, 0, 3, 4), (us, us, K, us**3/(kappa*H))):
        volume = sum((e["bounds"][1]-e["bounds"][0])*(e["bounds"][3]-e["bounds"][2])
                     for e in metadata["layout"][i]["mesh"])
        scales.append(field_scale/time_scale*math.sqrt(volume))
    for i in (5, 6):
        volume = sum((e["bounds"][1]-e["bounds"][0])*(e["bounds"][3]-e["bounds"][2])
                     for e in metadata["layout"][i]["mesh"])
        scales.append(us/H*math.sqrt(volume))
    return (np.array(report["residual_norms"][:6])/scales).tolist()


def check_evolving(report, metadata, expected_solved):
    require(report["solved"] == expected_solved and not report["native_time_map"], "Stage solve/audit was confused with a native timestep")
    spaces = [Space(metadata["layout"][i]) for i in (0, 0, 3, 4, 5, 6)]
    offsets = metadata["offsets"]
    state, history = np.array(report["state"]), np.array(report["history"])
    parameters, time = metadata["parameters"], report["time"]
    H, K = parameters[1], parameters[5]
    degree = metadata["layout"][0]["mesh"][0]["degree"]
    z, weights = np.polynomial.legendre.leggauss(14)
    L2, rateL2, inner_rate_squared = np.zeros(6), np.zeros(4), np.zeros(3)
    source_rate_error = 0.
    for block, space in enumerate(spaces):
        for e in space.elements:
            a, b, c, d = e["bounds"]
            X, R = np.meshgrid(a+(z+1)*(b-a)/2, c+(z+1)*(d-c)/2)
            X, R = X.ravel(), R.ravel()
            W = np.outer(weights, weights).ravel()*(b-a)*(d-c)/4
            ids, p, _, _ = space.rows(e, X, R)
            numerical = p@state[offsets[block]+ids]
            exact = np.array([manufactured_values(x, r, time, parameters)[block] for x, r in zip(X, R)])
            L2[block] += np.sum(W*(numerical-exact)**2)
            if block < 4:
                rate = report["alpha"]*(p@(state[offsets[block]+ids]-history[offsets[block]+ids]))
                expected = np.array([manufactured_values(x, r, time, parameters, True)[block] for x, r in zip(X, R)])
                rateL2[block] += np.sum(W*expected**2)
                if block < 3 and d <= H:
                    inner_rate_squared[block] += np.sum(W*expected**2)
                source_rate_error = max(source_rate_error, float(np.sqrt(np.sum(W*(rate-expected)**2))))
    L2, rateL2 = np.sqrt(L2), np.sqrt(rateL2)
    scaled = scaled_residual_norms(report, metadata)
    if degree >= 6:
        require(max(scaled) < 1e-9, "Resolved evolving manufactured stage exceeds the frozen scaled residual tolerance")
        require(max(L2) < 1e-8 and source_rate_error < 1e-8, "Manufactured state/storage failed to converge")
    inner_rate = np.sqrt(inner_rate_squared)
    require(inner_rate[0]+inner_rate[1] > 1e-4 and inner_rate[2] > 1e-7, "Manufactured inner storage vanished")
    if degree >= 6:
        samples = np.array(report["forcing_samples"])
        near(samples[:, 1], np.full(3, K), "evolving_physical_sampled_K")
        near(samples[:, 2], np.full(3, -.03*math.cos(time)), "evolving_outer_one_sided_forcing")
    v, ks, es, pw, po = spaces[0], spaces[2], spaces[3], spaces[4], spaces[5]
    active_yap_points, inactive_yap_points = 0, 0
    normal_points = np.polynomial.legendre.leggauss(max(8, degree+3))[0]
    for e in es.elements:
        _, _, lower, upper = e["bounds"]
        for r in lower+(normal_points+1)*(upper-lower)/2:
            exact = manufactured_values(.5, r, time, parameters)
            active = exact[2]**1.5/(2.5*exact[3]*r) > 1
            active_yap_points += int(active)
            inactive_yap_points += int(not active)
    require(active_yap_points > 0 and inactive_yap_points > 0, "Manufactured case missed a Yap branch")
    trace, no_slip, k_ground, traction, k_flux, epsilon = [], [], [], [], [], []
    for x in (.13, .47, .91):
        uw, guw = v.evaluate(state, offsets[0], x, H, -1)
        ww, gww = v.evaluate(state, offsets[1], x, H, -1)
        uo, guo = v.evaluate(state, offsets[0], x, H, 1)
        wo, gwo = v.evaluate(state, offsets[1], x, H, 1)
        kw, gkw = ks.evaluate(state, offsets[2], x, H, -1)
        ko, gko = ks.evaluate(state, offsets[2], x, H, 1)
        eo, _ = es.evaluate(state, offsets[3], x, H, 1)
        pi, _ = pw.evaluate(state, offsets[4], x, H, -1)
        pressure_o, _ = po.evaluate(state, offsets[5], x, H, 1)
        sampled, _ = ks.evaluate(state, offsets[2], x, parameters[3])
        profile = Profiles({"geometry": parameters[:5], "metric": [0, 0], "k": sampled,
                            "q": parameters[6], "kappa": parameters[7], "chi": 1})
        d = parameters[0]+float(profile.a)*max(H-float(profile.zv), 0)
        G = np.array([[guo[0], 0, guo[1]], [0, 0, 0], [gwo[0], 0, gwo[1]]])
        closure = indexed(G, ko, eo, "physical")
        R = closure["stress"].real
        tw = np.array([d*(guw[1]+gww[0]), -pi+2*d*gww[1]])
        to = np.array([-R[0, 2], -pressure_o+2*ko/3-R[2, 2]])
        trace.append(max(abs(uw-uo), abs(ww-wo), abs(kw-ko)))
        traction.append(float(np.max(np.abs(tw-to))))
        k_flux.append(abs(d*gkw[1]-closure["nut"].real*gko[1]))
        epsilon.append(abs(eo-float(profile.A/H)))
        u0, _ = v.evaluate(state, offsets[0], x, 0, 1)
        w0, _ = v.evaluate(state, offsets[1], x, 0, 1)
        _, g0 = ks.evaluate(state, offsets[2], x, 0, 1)
        no_slip.append(max(abs(u0), abs(w0)))
        k_ground.append(abs((parameters[0]+float(profile.a)*max(-float(profile.zv), 0))*g0[1]))
        near(ww, .01+.005*math.cos(time), "manufactured_shared_nonzero_normal_velocity")
    if degree >= 6:
        for name, value in (("trace", trace), ("traction", traction), ("k_flux", k_flux),
                            ("epsilon", epsilon), ("ground_velocity", no_slip), ("ground_k_flux", k_ground)):
            near(value, np.zeros(3), "evolving_interface_"+name, tolerance=1e-8)
    near(trace, np.zeros(3), "one_shared_velocity_k_trace")
    # Independently integrate regional mass budgets using physical traces and
    # volume derivatives. The nonzero shared normal flow crosses H only once.
    mass_budgets = []
    for lower, upper in ((0., H), (H, 4.)):
        divergence = 0.
        for e in v.elements:
            a, b, c, d = e["bounds"]
            if c < lower or d > upper:
                continue
            X, R = np.meshgrid(a+(z+1)*(b-a)/2, c+(z+1)*(d-c)/2)
            ids, _, dx, dr = v.rows(e, X.ravel(), R.ravel())
            div = dx@state[offsets[0]+ids]+dr@state[offsets[1]+ids]
            divergence += float(np.dot(np.outer(weights, weights).ravel(), div))*(b-a)*(d-c)/4
        radial = lower+(z+1)*(upper-lower)/2
        side = sum(wt*(v.evaluate(state, offsets[0], 1., r)[0]
                       -v.evaluate(state, offsets[0], 0., r)[0])
                   for r, wt in zip(radial, weights))*(upper-lower)/2
        normal = sum(wt*(v.evaluate(state, offsets[1], (x+1)/2, upper, -1)[0]
                         -v.evaluate(state, offsets[1], (x+1)/2, lower, 1)[0])
                     for x, wt in zip(z, weights))/2
        near([divergence, side+normal], np.zeros(2), "evolving_regional_mass_budget")
        mass_budgets.append({"normal_interval": [lower, upper], "volume_divergence": divergence,
                             "side_mass_flux": side, "normal_mass_flux": normal,
                             "net_boundary_mass_flux": side+normal})
    return {"degree": degree, "time": time, "alpha": report["alpha"], "solved": report["solved"],
            "iterations": report["iterations"], "physical_block_residual_norms": report["residual_norms"],
            "scaled_block_residual_norms": scaled, "active_yap_normal_points": active_yap_points,
            "inactive_yap_normal_points": inactive_yap_points,
            "field_L2_errors": L2.tolist(), "storage_rate_L2": rateL2.tolist(),
            "inner_storage_rate_L2": inner_rate.tolist(), "regional_mass_budgets": mass_budgets,
            "storage_rate_L2_error": source_rate_error, "shared_trace_max": max(trace),
            "traction_jump_max": max(traction), "k_diffusive_flux_jump_max": max(k_flux),
            "epsilon_trace_error_max": max(epsilon), "ground_velocity_max": max(no_slip),
            "ground_k_flux_max": max(k_ground)}


def qualify_cases(build, output, template, report, base):
    def run(command, directory, name, expected_failure=False, timeout=900):
        directory.mkdir(parents=True, exist_ok=True)
        start = time.monotonic()
        result = subprocess.run(list(map(str, command)), cwd=directory, capture_output=True, text=True, timeout=timeout)
        log = directory/(name+".log")
        log.write_text(result.stdout+result.stderr)
        report["commands"].append({"command": list(map(str, command)), "exit_code": result.returncode,
                                   "expected_failure": expected_failure, "seconds": time.monotonic()-start,
                                   "log_sha256": digest(log), "binary_sha256": digest(command[0])})
        require((result.returncode != 0) == expected_failure, f"Unexpected {name} exit; see {log}")
        return result
    probe, solver = build/"RANSEquilibriumPatchProbe", build/"IncNavierStokesSolver"
    def solve_session(session, name, metadata, evolving=False):
        print("M3 running "+name, flush=True)
        result = run([solver, session], session.parent, name)
        data = json.loads(session.with_suffix(".patch.json").read_text())
        if evolving:
            summary = check_evolving(data, metadata, True)
        else:
            require(data["solved"] and data["iterations"] > 0 and not data["native_time_map"], "No actual nonlinear stage solve")
            scaled = scaled_residual_norms(data, metadata)
            require(max(scaled) < 1e-9, "Stage scaled physical residual exceeds target")
            near(data["history_as_state_residual"], np.zeros(len(data["state"])), "stationary_manufactured_target_residual")
            parsed = [float(v) for v in re.findall(r"L (?:2|inf) error \(variable [^)]+\) : ([^\n]+)", result.stdout)]
            require(len(parsed) == 10 and max(parsed) < 1e-9, "Manufactured native output fields missed the exact state")
            summary = {"iterations": data["iterations"], "physical_block_residual_norms": data["residual_norms"],
                       "scaled_block_residual_norms": scaled, "field_errors": parsed}
        return data, summary
    stationary = []
    for condensed in (False, True):
        name = "stationary-condensed" if condensed else "stationary-full"
        session = output/name/"case.xml"
        root = ET.parse(template).getroot()
        set_info(root, "RANSWallCondense", "True" if condensed else "False")
        session.parent.mkdir()
        ET.ElementTree(root).write(session)
        data, summary = solve_session(session, name, base)
        stationary.append(data)
        report[name] = summary
    near(stationary[0]["state"], stationary[1]["state"], "stationary_full_condensed_state")
    # The native source readback must reproduce the independently specified
    # target. The old max(a,0) expression is a deliberately failing API control.
    K, q, kappa, H = base["parameters"][5], base["parameters"][6], base["parameters"][7], base["parameters"][1]
    E0 = (K/q)**1.5/(kappa*H)
    near(np.array(base["fixture_source_probes"])[:, 2:], [[E0, 1.92*E0*E0/K]]*3, "native_source_expression_readback")
    bad = output/"literal-function-control"/"case.xml"
    bad.parent.mkdir()
    root = ET.parse(template).getroot()
    root.find("CONDITIONS/FUNCTION[@NAME='ScalarSourceOuter']/E[@VAR='eps_tilde']").set(
        "VALUE", "1.92*E0^2/K0-0.83*E0^2/K0*max(YapScale/y-1,0)*(YapScale/y)^2")
    ET.ElementTree(root).write(bad)
    run([probe, bad, bad.parent/"sources.json", "sources"], bad.parent, "literal-control")
    bad_values = json.loads((bad.parent/"sources.json").read_text())["fixture_source_probes"]
    require(np.max(np.abs(np.array(bad_values)[:, 3]-1.92*E0*E0/K)) > 1e-4, "Expected 5.9 literal-argument incompatibility was not reproduced")
    report["source_expression_control"] = {"original_max_literal": "failed_as_observed_on_5.9",
                                            "equivalent_unary_positive_part": "passed", "legacy_readback": bad_values}
    invalid = output/"invalid-matching-height"/"case.xml"
    invalid.parent.mkdir()
    root = ET.parse(template).getroot()
    set_parameter(root, "RANSWallHeight", 1.9)
    ET.ElementTree(root).write(invalid)
    result = run([probe, invalid, invalid.parent/"invalid.json", "layout"], invalid.parent, "invalid-H", True)
    require("element crosses ground or matching height" in result.stderr, "Wrong invalid-geometry rejection")
    model = Profiles({"geometry": base["parameters"][:5], "metric": [0, 0], "k": K, "q": q, "kappa": kappa, "chi": 1})
    ell = float(model.nu/model.a-model.zv)
    def knots(factor):
        result = {0., H, base["parameters"][2]}
        r = 0.
        while True:
            r = (r+ell)*factor-ell
            if r >= H:
                break
            result.add(r)
        return sorted(result)
    p_results = []
    for degree in (4, 6, 8):
        name = f"parallel-p{degree}"
        session = parallel_session(template, output/name/"case.xml", degree, knots(1.4))
        print("M3 running "+name, flush=True)
        run([probe, session, session.parent/"parallel.json", "parallel"], session.parent, name)
        raw = json.loads((session.parent/"parallel.json").read_text())
        p_results.append({"degree": degree, "inner_elements": len(knots(1.4))-1, "cases": check_parallel(raw)})
    for j in range(3):
        e4, e6, e8 = (p["cases"][j]["velocity_RMS_error"] for p in p_results)
        require(e6 < .1*e4+1e-12 and e8 < .1*e6+1e-12 and e8 < 1e-9, "Native parallel p refinement failed")
        require(p_results[-1]["cases"][j]["weak_traction_error"] < 1e-10, "Refined AWF traction mismatch")
    report["parallel_p_refinement"] = p_results
    h_results = []
    coarse = knots(2.5)
    for level in (1, 2, 4):
        grid = sorted({a+(b-a)*j/level for a, b in zip(coarse[:-1], coarse[1:]) for j in range(level+1)})
        name = f"parallel-h{level}"
        session = parallel_session(template, output/name/"case.xml", 4, grid)
        run([probe, session, session.parent/"parallel.json", "parallel"], session.parent, name)
        raw = json.loads((session.parent/"parallel.json").read_text())
        h_results.append({"subdivision": level, "inner_elements": len(grid)-1, "cases": check_parallel(raw)})
    for j in range(3):
        e = [r["cases"][j]["velocity_RMS_error"] for r in h_results]
        rates = [math.log(e[i]/e[i+1], 2) for i in range(2)]
        require(min(rates) > 3.5, "Native parallel h refinement has insufficient consistency")
        h_results[-1]["cases"][j]["h_rates"] = rates
    report["parallel_h_refinement"] = h_results
    metadata_by_degree, audits = {}, []
    specifications = [(4, .3, 5.), (6, .3, 5.), (8, .3, 5.), (6, 0., 5.), (6, .9, 5.), (6, .3, 2.)]
    for degree, at_time, alpha in specifications:
        name = f"evolving-audit-p{degree}-t{at_time}-a{alpha}"
        session = evolving_session(template, output/name/"case.xml", degree, at_time, "Residual", alpha=alpha, perturb=False)
        if degree not in metadata_by_degree:
            run([probe, session, session.parent/"layout.json", "layout"], session.parent, "layout")
            metadata_by_degree[degree] = json.loads((session.parent/"layout.json").read_text())
        print("M3 running "+name, flush=True)
        run([solver, session], session.parent, name)
        data = json.loads(session.with_suffix(".patch.json").read_text())
        audits.append(check_evolving(data, metadata_by_degree[degree], False))
    coarse_error = max(audits[0]["scaled_block_residual_norms"])
    require(coarse_error > 1e-8 and max(audits[1]["scaled_block_residual_norms"]) < .01*coarse_error,
            "Manufactured finite-space defect did not refine away")
    report["evolving_manufactured_audits"] = audits
    solutions = []
    for condensed in (False, True):
        name = "evolving-condensed" if condensed else "evolving-full"
        session = evolving_session(template, output/name/"case.xml", 6, .3, "Solve", condensed)
        data, summary = solve_session(session, name, metadata_by_degree[6], True)
        solutions.append(data)
        report[name] = summary
    near(solutions[0]["state"], solutions[1]["state"], "evolving_full_condensed_state")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--suite", choices=["residuals", "all"], default="all")
    args = parser.parse_args()
    require(sys.platform.startswith("linux"), "Authoritative numerical validation requires Linux")
    root = Path(__file__).resolve().parents[4]
    output = args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix="cls-epawf-m3-", dir=args.build.resolve()))
    if args.output:
        output.mkdir(parents=True, exist_ok=False)
    report = {"schema": 1, "milestone": "cls_epawf_revised_m3", "suite": args.suite,
              "utc": datetime.now(timezone.utc).isoformat(), "host": platform.node(),
              "python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
              "cpu_affinity": sorted(os.sched_getaffinity(0)), "mpi_ranks": 1,
              "thread_environment": THREADS, "commands": [], "verification_status": "running"}
    try:
        session = output/"EquilibriumWallFlat.xml"
        shutil.copyfile(root/"solvers/IncNavierStokesSolver/Tests/RANS/EquilibriumWallFlat.xml", session)
        command = [str(args.build.resolve()/"RANSEquilibriumPatchProbe"), str(session), str(output/"patch.json")]
        start = time.monotonic()
        result = subprocess.run(command, cwd=output, capture_output=True, text=True, timeout=600)
        (output/"native.log").write_text(result.stdout+result.stderr)
        report["commands"].append({"command": command, "exit_code": result.returncode,
                                   "seconds": time.monotonic()-start, "binary_sha256": digest(command[0])})
        require(result.returncode == 0, "Native patch probe failed: "+str(output/"native.log"))
        raw = json.loads((output/"patch.json").read_text())
        report["layout"] = check_layout(raw)
        report["pressure_and_interface"] = check_pressure(raw)
        report["residual_cases"] = [residual_oracle(raw, c) for c in raw["generic_cases"]]
        report["condensation"] = check_condensation(raw)
        report["coarse_parallel_reduction"] = check_parallel(raw)
        if args.suite == "all":
            qualify_cases(args.build.resolve(), output,
                          root/"solvers/IncNavierStokesSolver/Tests/RANS/EquilibriumWallFlat.xml", report, raw)
        report["verification_status"] = "passed"
        report["milestone_closure"] = args.suite == "all"
    except Exception as error:
        report["verification_status"] = "failed"
        report["error"] = str(error)
        raise
    finally:
        report["checks"], report["scaled_errors"] = checks, errors
        report["source_commit"] = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
        report["source_status"] = subprocess.check_output(["git", "-C", str(root), "status", "--porcelain"], text=True).splitlines()
        names = ["RANSWallFieldLayout.h", "RANSWallFieldLayout.cpp", "RANSWallPatch.h", "RANSWallPatch.cpp",
                 "RANSWallInterface.h", "RANSWallInterface.cpp", "RANSEquilibriumWallSystem.h", "RANSEquilibriumWallSystem.cpp"]
        sources = ["solvers/IncNavierStokesSolver/EquationSystems/"+name for name in names]
        sources += ["library/UnitTests/SolverUtils/RANSEquilibriumPatchProbe.cpp",
                    "solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_epawf_m3.py",
                    "solvers/IncNavierStokesSolver/Tests/RANS/cls_epawf_m3_cases.py",
                    "solvers/IncNavierStokesSolver/Tests/RANS/EquilibriumWallFlat.xml",
                    "solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_epawf_m1.py",
                    "solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_epawf_m2.py",
                    "solvers/IncNavierStokesSolver/Tests/RANS/CMakeLists.txt",
                    "solvers/IncNavierStokesSolver/CMakeLists.txt",
                    "solvers/IncNavierStokesSolver/EquationSystems/IncNavierStokes.h",
                    "docs/literature/cls_neutral_abl_equilibrium_awf_note_revised.tex",
                    "docs/plans/cls-epawf-revised-contract-v1.json"]
        sources += ["library/SolverUtils/Turbulence/"+name+extension
                    for name in ("RANSModels", "AnalyticalWallFunction", "EquilibriumAWFProfiles")
                    for extension in (".h", ".cpp")]
        report["source_sha256"] = {name: digest(root/name) for name in sources}
        report["artifacts_sha256"] = {p.name: digest(p) for p in output.iterdir() if p.is_file()}
        (output/"validation.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"verification_status": report["verification_status"], "checks": checks,
                      "milestone_closure": report["milestone_closure"], "report": str(output/"validation.json")}, indent=2))


if __name__ == "__main__":
    main()
