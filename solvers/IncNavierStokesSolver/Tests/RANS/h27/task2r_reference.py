#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Construct and independently qualify the H27 finite-nu exact-inner table."""

import argparse
from dataclasses import dataclass
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import platform

import numpy as np
from numpy.polynomial import chebyshev as cheb
from numpy.polynomial.legendre import leggauss
import scipy
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import PPoly
from scipy.optimize import brentq

from inner_k_columns import planar_response, stress
from resources import atomic, identity

KIN, H, D0, ZT, L0 = .008809015145515192, .15, .0875, 20., .04
NU, HS, Z0, UTAU = 6.423173803526447e-5, .11541061001735275, .005676270321216798, .051
QE, NE, CD, KAPPA = 81 / (2 * math.sqrt(143)), 81 / 20, 572 / 6561, .3819178635467331
B, C, TAU0 = .4 * H, .75 * H, UTAU ** 2
RT = ZT - D0
AIN = KIN ** 1.5 / (KAPPA * QE ** 1.5)
UTARGET = UTAU / KAPPA * math.log(H / Z0)
Q_LIMIT = 3.6 / .65
MIN_K, MAX_K = .9 * KIN, 1.1 * KIN
ARCHIVE = "1d4fff1e0efe09820bf3189d2c670ac4b55303536f7af45664760fad4421a619"
SCHEMA, POLICY = "WJH27ExactInnerFiniteNuTablesV1", "ExactInnerFiniteNuV1"


def require(ok, reason):
    if not ok:
        raise ValueError(reason)


def blend(r):
    y = np.clip((np.asarray(r) - B) / (C - B), 0., 1.)
    return y ** 3 * (10. - 15. * y + 6. * y * y)


@dataclass(frozen=True)
class Column:
    K: float
    h: float = HS

    @property
    def A(self):
        return (self.K / QE) ** 1.5 / KAPPA

    @property
    def slope(self):
        return KAPPA * math.sqrt(self.K / QE)

    @property
    def shear(self):
        return self.K / QE + NU * QE * self.A / (H * self.K)

    @property
    def intercept(self):
        x = self.h * math.sqrt(self.K) / (70 * NU)
        m0 = .5 - .4 * x ** .7
        m1 = 1 - .79 * x ** (-.28)
        require(x > 1 and m1 > max(0., m0), "roughness branch changed")
        return -10.7 * NU * math.expm1(m1 * math.log(x)) / math.sqrt(self.K)

    def auxiliary(self, r):
        return self.slope * (r - self.intercept)

    def lower_production_integral(self, r):
        d = NU + self.auxiliary(0.)
        x = self.slope * r / d
        return self.shear ** 2 / self.slope * (math.log1p(x) - (NU / d) * x / (1 + x))

    def lower_velocity(self, r):
        return self.shear / self.slope * math.log1p(self.slope * r / (NU + self.auxiliary(0.)))


def reduced_shear(k, e, alpha, chi, total):
    """Increasing shear with s=q/N. The auxiliary-only law has no q limit."""
    def residual(s):
        q = 1.8 * s * (1 + s * s) / (1 - .35 * s * s)
        return alpha * e * q / k + chi * k * .6 * s / (1 + s * s) - total

    require(residual(1.) > 0, "no increasing upper-column branch")
    s = brentq(residual, 0., 1., xtol=9e-16)
    n = 1.8 * (1 + s * s) / (1 - .35 * s * s)
    return s * n * e / k


def q_root_shear(q, tangent=False):
    """Direct cubic root in q coordinates for the independent upper IVP."""
    q2 = q * q
    b, p = .216 + 1.005 * q2, .36 + (.35 / 3) * q2
    # The selected shear interval has a nonnegative discriminant.
    discriminant = q2 * ((.432 * 1.005 - .3888 * (.35 / 3)) +
                          q2 * (1.005 ** 2 - 1.08 * (.35 / 3) ** 2 -
                                q2 * (.35 / 3) ** 3))
    w = (b + math.sqrt(discriminant)) ** (1 / 3)
    n = .6 + w + p / w
    denominator = n * n + q2
    if tangent:
        dn = 2 * q * (.35 * n + 1.8) / (3 * n * n - 3.6 * n - .35 * q2)
        return .6 * ((n + q * dn) * denominator - n * q * (2 * n * dn + 2 * q)) / denominator ** 2
    return .6 * n * q / denominator


def independent_shear(k, e, alpha, chi, total):
    def residual(q):
        return alpha * e * q / k + chi * k * q_root_shear(q) - total

    require(residual(Q_LIMIT) > 0., "independent shear branch capacity")
    return brentq(residual, 0., Q_LIMIT, xtol=1.8e-15) * e / k


def integrate(rhs, edges, initial, method, tol):
    segments = []
    value = initial
    for hi, lo in zip(edges[:-1], edges[1:]):
        solution = solve_ivp(rhs, (hi, lo), value, method=method, dense_output=True,
                             rtol=tol, atol=.02 * tol, max_step=.018 * H)
        require(solution.success and solution.t[-1] == lo, "incomplete column IVP")
        value = solution.y[:, -1]
        segments.append((lo, hi, solution.sol))
    return value, segments


def state(segments, r):
    for lo, hi, sol in segments:
        if lo <= r <= hi:
            return sol(r)
    raise ValueError("radius leaves solved segments")


def upper(column, independent):
    a, zv, A, K, T = column.slope, column.intercept, column.A, column.K, column.shear
    cho = independent_shear if independent else reduced_shear
    method, tol = ("RK45", 2.4e-13) if independent else ("DOP853", 3e-14)

    def rhs(r, z):
        k = K * (math.exp(z[0]) if independent else z[0])
        require(k > 0., "nonpositive upper k")
        e = A / r
        x = float(blend(r))
        alpha = NU + (1. - x) * a * (r - zv)
        g = cho(k, e, alpha, x, T) if x > 0. else T / alpha
        d = NU + CD * k * k / e
        p = (T - NU * g) * g
        return [A * z[1] / (d * (k if independent else K)), (e - p) / A, g]

    full, segments = integrate(rhs, sorted({H, C, column.h, B}, reverse=True),
                               [0. if independent else 1., 0., 0.], method, tol)
    if independent:
        # Integrate the auxiliary shear budget independently of the reduced formula.
        pi = quad(lambda r: (T - NU * T / (NU + a * (r - zv))) *
                  T / (NU + a * (r - zv)), 0., B, epsabs=1e-16, epsrel=1e-12)[0]
    else:
        pi = column.lower_production_integral(B)
    logc = math.log(B / H) + 1 - full[1] - pi / A
    return logc, column.lower_velocity(B) - full[2], segments, full


@lru_cache(maxsize=None)
def logc(factor, independent=True):
    return upper(Column(KIN * factor), independent)[0]


def complete(column, c, solved):
    _, uh, above, end = solved
    A, K, a, zv, T = column.A, column.K, column.slope, column.intercept, column.shear
    edges = sorted({0., c, B, C, column.h, H})
    require(edges == [0., c, B, C, column.h, H], "material ordering changed")

    def rhs(r, y):
        k = K * math.exp(y[0])
        e = A / max(r, c)
        g = T / (NU + a * (r - zv))
        return [A * y[1] / ((NU + CD * k * k / e) * k), (e - (T - NU * g) * g) / A]

    _, below = integrate(rhs, [B, c, 0.], end[:2], "RK45", 2.4e-13)

    def evaluate(r):
        r = np.asarray(r)
        inside = r < B
        # Lower and upper IVPs have different state dimensions.
        y = np.empty((r.size, 3))
        if inside.any():
            y[inside] = [np.r_[state(below, s), column.lower_velocity(s) - uh] for s in r[inside]]
        if (~inside).any():
            y[~inside] = [state(above, s) for s in r[~inside]]
        k, Q, u = K * np.exp(y[:, 0]), A * y[:, 1], uh + y[:, 2]
        x, aux, e = blend(r), column.auxiliary(r), A / np.maximum(r, c)
        g = np.array([independent_shear(kk, ee, NU + (1 - xx) * aa, xx, T)
                      if xx > 0. else T / (NU + aa) for kk, ee, xx, aa in zip(k, e, x, aux)])
        R, p = tensor(k, e, g, aux, x)
        return dict(U=u, k=k, Q=Q, E=e, G=g, P=p, chi=x, aux=aux, R=R,
                    kr=Q / (NU + CD * k * k / e), Qr=e - p, Pi=2 * k / 3 - R[:, 2, 2] - K / NE)

    return edges, evaluate


def tensor(k, e, g, aux, x):
    grad = np.zeros((len(k), 3, 3))
    grad[:, 0, 2] = g
    R = stress(k, e, grad, aux, x)
    return R, -R[:, 0, 2] * g


def quadrature(edges, order):
    a, b = np.array(edges[:-1]), np.array(edges[1:])
    x, w = leggauss(order)
    return (.5 * (a[:, None] + b[:, None]) + .5 * (b - a)[:, None] * x).ravel(), \
           (.5 * (b - a)[:, None] * w).ravel()


def fit_cutoff():
    training = np.cos(np.pi * (np.arange(18) + .5) / 18)
    coeffs = cheb.chebfit(training, [logc(1. + .1 * x, False) for x in training], 12)
    derivative = cheb.chebder(coeffs) * 2 / (MAX_K - MIN_K)
    radius = lambda K: H * math.exp(cheb.chebval((2 * K - MIN_K - MAX_K) / (MAX_K - MIN_K), coeffs))
    tail = sum(abs(derivative[1:]))
    reserve = 32 * len(coeffs) * np.finfo(float).eps * (abs(derivative[0]) + tail)
    slope_range = [float(derivative[0] - tail - reserve), float(derivative[0] + tail + reserve)]
    lower, upper_radius = radius(MAX_K), radius(MIN_K)
    require(slope_range[1] < 0 and 0 < lower < upper_radius < B, "represented cutoff enclosure")

    def inverse(r):
        if r == lower:
            return MAX_K
        if r == upper_radius:
            return MIN_K
        require(lower < r < upper_radius, "inverse cuts outside qualification")
        return brentq(lambda K: radius(K) - r, MIN_K, MAX_K, xtol=5e-18)

    controls = np.linspace(.9, 1.1, 35)
    # Interleave near both ends, where the slope and extrapolation reserves matter most.
    heldout = np.unique(np.r_[controls, (.9 + .1 * (controls[1:7] - .9)),
                                  (1.1 - .1 * (1.1 - controls[-7:-1]))])
    require(min(abs((1. + .1 * training)[:, None] - heldout).ravel()) > 1e-10,
            "training and held-out K overlap")
    rows = []
    for factor in heldout:
        K = KIN * factor
        column = Column(K)
        c = radius(K)
        true_log, uh, pieces, end = upper(column, True)
        direct = H * math.exp(true_log)
        edges, evaluate = complete(column, c, (true_log, uh, pieces, end))
        endpoints = evaluate(np.array([0., H]))
        r, w = quadrature(edges, 48)
        values = evaluate(r)
        r = np.unique(np.r_[np.linspace(0., H, 257), edges])
        v = evaluate(r)
        guard, _ = planar_response(v["k"], v["E"], v["G"], v["aux"], v["chi"], NU)
        require(guard["minimum_covariance_over_k"] > .25 and guard["minimum_principal_over_mu"] > .01
                and guard["tangent_refinement_difference"] <= 1e-8, f"sampled full-tensor gate at K/Kin={factor}")
        shear_error = np.max(abs(NU * v["G"] - v["R"][:, 0, 2] - column.shear)) / TAU0

        def fd(step, calculation=logc):
            if .9 <= factor - 2 * step and factor + 2 * step <= 1.1:
                return sum(wt * calculation(factor + i * step) for i, wt in ((-2, 1), (-1, -8), (1, 8), (2, -1))) / (12 * step)
            sign = 1. if factor < 1. else -1.
            return sign * sum(wt * calculation(factor + sign * i * step)
                              for i, wt in enumerate((-25, 48, -36, 16, -3))) / (12 * step)

        slope = cheb.chebval((2 * K - MIN_K - MAX_K) / (MAX_K - MIN_K), derivative)
        fd1, fd2 = fd(.001), fd(.0005)
        slope_error = max(abs(slope * KIN - fd1), abs(slope * KIN - fd2))
        shear_tangent = NU + (1 - v["chi"]) * v["aux"]
        active = v["chi"] > 0
        q = (v["k"] * v["G"] / v["E"])[active]
        shear_tangent[active] += (v["chi"] * v["k"] ** 2 / v["E"])[active] * [
            q_root_shear(qi, tangent=True) for qi in q]
        require(np.min(shear_tangent) > 0, "nonpositive shear slope")
        epsilon_error = []
        for r in (.8 * c, 1.2 * c, H):
            e = column.A / max(r, c)
            analytic = e * (1.5 / K - (slope if r < c else 0.))
            direct_e = lambda f: Column(KIN * f).A / max(r, H * math.exp(logc(f)))
            epsilon_error.append(abs(analytic * KIN - fd(.0005, direct_e)) / e)
        recovered = inverse(c)
        # For true values just outside the represented endpoint range, bound the
        # inverse bias using the certified minimum slope instead of clamping.
        direct_inverse = abs(inverse(direct) - K) / KIN if lower < direct < upper_radius \
            else abs(math.log(direct / c) / slope_range[1]) / KIN
        require(np.min(v["k"]) > 0. and np.min(v["E"]) > 0. and np.min(v["G"]) > 0.,
                "nonpositive represented validation column")
        rows.append(dict(K_factor=float(factor), K=float(K), c=float(c), independent_c=direct,
                         relative_cutoff_error=abs(c / direct - 1),
                         K_in_times_derivative_error=float(slope_error),
                         K_times_log_derivative_error=float(factor * slope_error),
                         coarse_fd=float(fd1), fine_fd=float(fd2), K_c_prime_over_c=float(K * slope),
                         c_prime=float(c * slope), inverse_K_error_over_Kin=max(direct_inverse, abs(recovered - K) / KIN),
                         represented_inverse_residual=float(abs(radius(recovered) - c)),
                         ground_flux_over_A=float(endpoints["Q"][0] / column.A),
                         interface_flux_over_A=float(endpoints["Q"][1] / column.A),
                         interface_K_error_over_Kin=float((endpoints["k"][1] - K) / KIN),
                         independent_budget_over_A=float(w @ (values["P"] - values["E"]) / column.A),
                         total_shear_error_over_tau0=float(shear_error),
                         epsilon_K_derivative_fixed_radius_error=float(max(epsilon_error)),
                         minimum_shear_tangent=float(np.min(shear_tangent)),
                         min_k_over_K=float(np.min(v["k"]) / K), minimum_shear=float(np.min(v["G"])),
                         UH=uh, guard=guard))

    gates = {key: max(abs(row[key]) for row in rows) for key in (
        "relative_cutoff_error", "K_in_times_derivative_error", "K_times_log_derivative_error",
        "inverse_K_error_over_Kin", "represented_inverse_residual", "ground_flux_over_A",
        "interface_flux_over_A", "interface_K_error_over_Kin", "independent_budget_over_A",
        "total_shear_error_over_tau0", "epsilon_K_derivative_fixed_radius_error")}
    budget = max(gates[key] for key in ("ground_flux_over_A", "interface_flux_over_A",
                                       "interface_K_error_over_Kin", "independent_budget_over_A",
                                       "total_shear_error_over_tau0"))
    require(budget <= 1e-10 and gates["relative_cutoff_error"] <= 1e-10
            and gates["K_in_times_derivative_error"] <= 1e-8
            and gates["K_times_log_derivative_error"] <= 1e-8
            and gates["epsilon_K_derivative_fixed_radius_error"] <= 1e-8
            and gates["inverse_K_error_over_Kin"] <= 1e-9
            and gates["represented_inverse_residual"] <= 256 * np.finfo(float).eps * H,
            f"held-out cutoff gates {gates}")
    return dict(argument="K", value="log(c/H)", domain=[MIN_K, MAX_K],
                degree=len(coeffs) - 1, coefficients=coeffs.tolist(),
                logcutoff_dK_coefficients=derivative.tolist(), heldout_K_count=len(rows),
                training_K_factors=(1. + .1 * training).tolist(), validations=rows,
                log_derivative_bounds=slope_range, range=[lower, upper_radius], maxima=gates,
                maximum_budget_or_ground_flux_error=budget,
                maximum_K_in_times_derivative_error=gates["K_in_times_derivative_error"])


def independent_bvp(root, r):
    with np.load(root / "independent/collocation_solution.npz") as data:
        poly = PPoly(data["coefficients"], data["breaks"])
        c = B / (1 + math.exp(-data["parameters"][1]))
        hs = math.exp(data["parameters"][0])
    lower = r <= c
    stride = np.where(lower, c, H - c)
    actual = np.where(lower, r, r - c) / stride
    v, dv = poly(actual), poly.derivative()(actual) / stride[:, None]
    k = KIN * np.exp(np.where(lower, v[:, 1], v[:, 4]))
    u, G = np.where(lower, v[:, 0], v[:, 3]), np.where(lower, dv[:, 0], dv[:, 3])
    Q = AIN * np.where(lower, v[:, 2], v[:, 5])
    kr = k * np.where(lower, dv[:, 1], dv[:, 4])
    e = AIN / np.maximum(r, c)
    x = blend(r)
    R, P = tensor(k, e, G, Column(KIN, hs).auxiliary(r), x)
    return dict(U=u, G=G, k=k, kr=kr, Q=Q, Qr=AIN * np.where(lower, dv[:, 2], dv[:, 5]),
                Pi=2 * k / 3 - R[:, 2, 2] - KIN / NE, P=P)


def make_inlet(repro, cutoff):
    c = H * math.exp(cheb.chebval((2 * KIN - MIN_K - MAX_K) / (MAX_K - MIN_K), cutoff["coefficients"]))
    column = Column(KIN)
    material, evaluate = complete(column, c, upper(column, True))
    edges = sorted(set(material + [.001, .0025, .005, .01, .02, .045, .085, .13]))
    degree = 26
    nodes = np.cos(np.pi * np.arange(degree + 1) / degree)
    error = {field: 0. for field in ("U", "G", "k", "kr", "Q", "Qr", "Pi", "P")}
    scales = dict(U=1., G=1 / H, k=KIN, kr=KIN / H, Q=AIN, Qr=AIN / H, Pi=TAU0, P=AIN / H)
    panels = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        midpoint, half = (lo + hi) / 2, (hi - lo) / 2
        sample = evaluate(np.clip(midpoint + half * nodes, lo, hi))
        panels.append(dict(domain=[lo, hi]))
        for field, grad in (("U", "G"), ("k", "kr"), ("Q", "Qr")):
            values = cheb.chebint(cheb.chebfit(nodes, sample[grad], degree),
                                  k=[sample[field][-1]], lbnd=-1., scl=half)
            values[:2] += (sample[field][0] - cheb.chebval(1., values)) / 2
            panels[-1][field] = values.tolist()
        for field in ("Pi", "P"):
            panels[-1][field] = cheb.chebfit(nodes, sample[field], degree).tolist()
        check = np.unique(np.r_[-1., 1., leggauss(71)[0], 2 * (np.arange(89) + .27) / 89 - 1])
        r = np.clip(midpoint + half * check, lo, hi)
        witness = independent_bvp(repro, r)
        values = {key: cheb.chebval(check, panels[-1][key]) for key in ("U", "k", "Q", "Pi", "P")}
        values.update({key: cheb.chebval(check, cheb.chebder(panels[-1][field]) / half)
                       for key, field in (("G", "U"), ("kr", "k"), ("Qr", "Q"))})
        for key in values:
            error[key] = max(error[key], float(np.max(abs(values[key] - witness[key])) / scales[key]))
        # Differential defects of the frozen profiles, using independent full stress.
        e = AIN / np.maximum(r, c)
        d = NU + CD * values["k"] ** 2 / e
        R, P = tensor(values["k"], e, values["G"], column.auxiliary(r), blend(r))
        defects = [np.max(abs(d * values["kr"] - values["Q"])) / AIN,
                   np.max(abs(values["Qr"] - e + P)) * H / AIN,
                   np.max(abs(NU * values["G"] - R[:, 0, 2] - column.shear)) / TAU0]
        require(max(defects) <= 1e-8, f"differential defect on [{lo}, {hi}]: {defects}")
        panels[-1]["differential_defects"] = list(map(float, defects))
    require(max(error.values()) <= 1e-8, f"independent fixed-scale profile gate {error}")
    return dict(cutoff=c, panels=panels, material_edges=material, scaled_errors=error, scales=scales)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode()).hexdigest()


def qualify(repro, out):
    bvp = json.loads((repro / "independent/collocation_refined.json").read_text())
    audit = json.loads((repro / "audit/audit_summary.json").read_text())
    archive_members = dict(line.split()[::-1] for line in (repro / "MANIFEST.sha256").read_text().splitlines())
    validate = lambda h: upper(Column(KIN, h), False)[1] - UTARGET
    calibrated = brentq(validate, .10, .13, xtol=6e-16)
    require(bvp["validated_branch"] and abs(calibrated / HS - 1) <= 1e-8
            and abs(calibrated / bvp["hs"] - 1) <= 1e-8,
            "independent fixed roughness calibration gate")
    cutoff = fit_cutoff()
    inlet = make_inlet(repro, cutoff)
    require(abs(inlet["cutoff"] / bvp["cutoff"] - 1) <= 1e-8 and
            abs(inlet["cutoff"] / .02439603217823862 - 1) <= 1e-8,
            "independent inlet cutoff gate")
    candidate = Column(KIN)
    _, evaluate = complete(candidate, inlet["cutoff"], upper(candidate, True))
    ends = evaluate(np.array([0., H]))
    g_out = UTAU / (KAPPA * H)
    top_shear = TAU0 + NU * UTAU / (KAPPA * RT)
    edge_defects = [ends["U"][0], ends["U"][1] - UTARGET, (ends["k"][1] - KIN) / KIN,
                    (ends["G"][1] - g_out) * H, ends["Pi"][1] / TAU0,
                    (candidate.shear - TAU0 - NU * g_out) / TAU0,
                    (top_shear - .002601430748812393) / TAU0]
    require(max(map(abs, edge_defects)) <= 1e-10 and
            max(map(abs, bvp["offgrid_collocation_scaled"])) <= 1e-8 and
            abs(bvp["budget_error_over_A"]) <= 1e-10 and bvp["max_bc_error"] <= 1e-10,
            "independent endpoints, budget or off-grid collocation gate")

    calibration = dict(h=HS, reproduced_h=calibrated, independent_h=bvp["hs"], independent_c=bvp["cutoff"],
                       roughness_intercept=candidate.intercept, top_velocity=UTARGET,
                       reference_shear=candidate.shear, historical_shear=TAU0, ground_k=float(ends["k"][0]),
                       interface_gradient=float(ends["G"][1]), endpoint_defects=list(map(float, edge_defects)),
                       inner_production=bvp["integrals"]["P"], inner_destruction=bvp["integrals"]["E"])
    loads = dict(outlet_inner="(-Kin/Ne-2*Pi, Tref(Kin))",
                 outlet_outer="(-Kin/Ne, tau0+nu*utau/(kappa*r))",
                 top_x=top_shear, interface_x=-KIN / NE, interface_r=candidate.shear,
                 top_e_out=NU * AIN / RT ** 2 + CD * KIN ** 2 / (1.3 * RT),
                 outlet_normal=[1, 0], top_normal=[0, 1])
    units = dict(L0_m=L0, H_code=H, H_m=H * L0, d0_code=D0, d0_m=D0 * L0,
                 interface_elevation_code=D0 + H, interface_elevation_m=(D0 + H) * L0,
                 top_elevation_code=ZT, top_elevation_m=ZT * L0, top_radius_code=RT, top_radius_m=RT * L0,
                 roughness_m=HS * L0, cutoff_m=inlet["cutoff"] * L0, z0e_code=Z0, z0e_m=Z0 * L0,
                 x_code=[-10, 15], x_m=[-.4, .6], velocity_scale="U0", time_scale="L0/U0",
                 resource_time="elapsed computational wall seconds")
    import inner_k_columns
    import physical_oracle

    source = dict(generator=identity(__file__)["sha256"], zip=ARCHIVE, archive_members=archive_members,
                  independent_tensor=identity(physical_oracle.__file__)["sha256"],
                  independent_symbol=identity(inner_k_columns.__file__)["sha256"])
    geometry = dict(H=H, d0=D0, rT=RT, force_support=HS, collar=[B, C], material_edges=inlet["material_edges"])
    table = dict(schema=SCHEMA, policy=POLICY, status="PASS",
                 model="IsotropicLowerWallFullWJCollarV2", transport="ScalarCdV2",
                 coefficient_convention="Ascending Chebyshev T_j coefficients; x=(2*argument-lo-hi)/(hi-lo)",
                 parameters=dict(K=KIN, H=H, h=HS, nu=NU, cd=CD, sigma_k=1., inner_diffusion="scalar",
                                 endpoint_flux_over_A=0.),
                 constants=dict(q_e=QE, N_e=NE, kappa=KAPPA, c_length=KAPPA * QE ** 1.5),
                 inlet_binding=dict(d0=D0, top_model_radius=RT, u_tau=UTAU, K_in=KIN, roughness=HS, A_in=AIN,
                                    epsilon_in_H=AIN / H, z0e=Z0),
                 guard=dict(sigma_e=1.3, reference_endpoint_flux_over_A=0.),
                 K_interval=[MIN_K, MAX_K], K_factor_interval=[.9, 1.1], cutoff_table=cutoff,
                 calibration=calibration, inlet=inlet, loads=loads, geometry=geometry, units=units, sources=source)
    table["bindings"] = {name: digest(table[name]) for name in
                         ("calibration", "cutoff_table", "inlet", "loads", "geometry", "units", "sources")}
    diagnostic = dict(inner_parallel_exact=True,
                      outer_rhs_integrated=dict(momentum=NU * UTAU / KAPPA * (1 / RT - 1 / H), k=0.,
                                                epsilon=NU * AIN * (1 / H ** 2 - 1 / RT ** 2)),
                      independently_integrated=dict(momentum=quad(lambda r: -NU * UTAU / (KAPPA * r * r), H, RT)[0],
                                                    k=0., epsilon=quad(lambda r: 2 * NU * AIN / r ** 3, H, RT)[0]))
    out.mkdir()
    atomic(out / "tables-exact-inner-v1.json", table)
    atomic(out / "qualification.json", dict(R1="PASS", archive=ARCHIVE,
        table=identity(out / "tables-exact-inner-v1.json"), calibration=calibration, geometry=geometry,
        fixed_scale_inlet_errors=inlet["scaled_errors"], cutoff_maxima=cutoff["maxima"], loads=loads,
        native_dependencies=dict(geometry=table["bindings"]["geometry"], policy=POLICY,
                                 legacy_table_sha256="b1bac2202d5991865a19c193528182954211f11f954fb5ca9b3f23998262250d"),
        heldout_K_count=cutoff["heldout_K_count"], sampled_full_tensor=audit, continuum_diagnostic=diagnostic,
        versions=dict(python=platform.python_version(), scipy=scipy.__version__, numpy=np.__version__)))
    print("R1 PASS", identity(out / "tables-exact-inner-v1.json")["sha256"], cutoff["maxima"], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reproduction", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    qualify(args.reproduction, args.out)


if __name__ == "__main__":
    main()
