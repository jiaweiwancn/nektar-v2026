#!/usr/bin/env python3
"""Independent Linux oracle for the revised EP-AWF spectral ray profiles.

Production budgets and redistribution fluxes use analytic metric primitives;
dual loads use a separate 64-point rule with factor-three physical panels.
Native segment/quad bases are checked against Jacobi polynomials. No coupled
patch, global curved mesh, or native solver timestep is claimed by this test.
"""
import argparse
import copy
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
import time

sys.dont_write_bytecode = True
THREADS = {key: "1" for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                               "MKL_NUM_THREADS", "BLIS_NUM_THREADS",
                               "VECLIB_MAXIMUM_THREADS")}
os.environ.update(THREADS)
import numpy as np
import scipy
from scipy.optimize import brentq
from scipy.special import eval_jacobi, eval_legendre

NOTE = "docs/literature/cls_neutral_abl_equilibrium_awf_note_revised.tex"
NOTE_SHA = "d35e437c2f2f8b25a46f7a3f285c6f1ad8b48d7b554aa1d161315f59eeb6de0f"
CONTRACT = "docs/plans/cls-epawf-revised-contract-v1.json"
CONTRACT_SHA = "4e7bd75022409b8dc2ef31cc73425779b646a97535760b31eeea4ec49087de8a"
LD = np.longdouble
errors = {}
checks = 0


def require(condition, message):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(message)


def near(actual, expected, name, scale, tolerance=1e-10):
    a, b = np.asarray(actual, dtype=LD), np.asarray(expected, dtype=LD)
    require(a.shape == b.shape, "Shape mismatch: " + name)
    require(np.all(np.isfinite(a)) and np.all(np.isfinite(b)), "Nonfinite: " + name)
    error = float(np.max(np.abs(a-b))/max(LD(scale), np.max(np.abs(b))))
    errors[name] = max(error, errors.get(name, 0.0))
    require(error <= tolerance, f"{name}: scaled error {error} > {tolerance}")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class Profiles:
    """Dimensional original law and analytic J-weighted budget primitives."""
    def __init__(self, case):
        self.nu, self.H, self.h, self.rK, self.rm = map(LD, case["geometry"])
        self.k, self.q, self.kappa = (LD(case[x]) for x in ("k", "q", "kappa"))
        self.chi = LD(case["chi"])
        k1, k2 = map(LD, case["metric"])
        self.s1, self.s2 = k1+k2, k1*k2
        root = np.sqrt(self.k)
        if self.h == 0:
            intercept = LD("10.7")
        else:
            x = self.h*root/(70*self.nu)
            m = max(LD(".5")-LD(".4")*x**LD(".7"),
                    1-LD(".79")*x**LD("-.28"))
            intercept = LD("10.7")*(1-x**m)
        self.zv = self.nu*intercept/root
        self.a = LD(".2295")*root
        self.d = self.nu-self.a*self.zv
        self.ze = LD("5.1")*self.nu/root
        self.t0 = self.k/self.q
        self.us = np.sqrt(self.t0)
        self.A = self.us**3/self.kappa
        self.lo = max(LD(0), self.zv)
        self.w0 = self.d+self.a*self.lo
        self.c0 = 1-self.s1*self.d/self.a+self.s2*(self.d/self.a)**2
        self.c1 = self.s1/self.a-2*self.s2*self.d/self.a**2
        self.c2 = self.s2/self.a**2
        z = min(self.H, self.lo)
        self.J0 = z/self.nu
        if self.H > z:
            self.J0 += np.log1p(self.a*(self.H-z)/(self.d+self.a*z))/self.a
        self.J1 = max(z-self.h, 0)**2/(2*self.nu)
        lower = max(z, self.h)
        if self.H > lower:
            log = np.log1p(self.a*(self.H-lower)/(self.d+self.a*lower))
            self.J1 += (self.H-lower)/self.a-(self.d+self.a*self.h)/self.a**2*log
        self.UH = self.t0*self.J0
        self.Pj = self.production_primitive(self.H)
        self.Bj = self.cap_primitive(self.H, self.A, self.rm)
        require(self.Pj >= self.Bj, "Oracle fixture is inadmissible")
        high = math.log(float(self.rm))
        low = high-1
        def budget(logrc):
            return float(self.cap_primitive(self.H, self.A, np.exp(LD(logrc)))-self.Pj)
        while budget(low) < 0:
            low -= 4
        self.rc = np.exp(LD(brentq(budget, low, high, xtol=5e-14, rtol=1e-14)))

    def jacobian(self, r):
        r = np.asarray(r, dtype=LD)
        return 1+self.s1*r+self.s2*r*r

    def volume(self, r):
        return r+self.s1*r*r/2+self.s2*r*r*r/3

    def production_primitive(self, r):
        r = np.maximum(np.asarray(r, dtype=LD), self.lo)
        delta = self.a*(r-self.lo)
        w = self.w0+delta
        # Integrate J(r)*(w-nu)/w^2 after w=d+a*r. log1p and
        # difference products avoid cancellation near the viscous breakpoint.
        return self.t0**2/self.a*(
            (self.c0-self.nu*self.c1)*np.log1p(delta/self.w0)
            +(self.c1-self.nu*self.c2)*delta
            +self.c2*delta*(w+self.w0)/2
            -self.nu*self.c0*delta/(w*self.w0))

    def cap_primitive(self, r, amplitude, cutoff):
        r = np.asarray(r, dtype=LD)
        below, above = np.minimum(r, cutoff), np.maximum(r, cutoff)
        return amplitude*(self.volume(below)/cutoff+np.log(above/cutoff)
                          +self.s1*(above-cutoff)+self.s2*(above-cutoff)*(above+cutoff)/2)

    def fields(self, r, velocity=None, forcing=(0, 0), supplied=None):
        r = np.asarray(r, dtype=LD)
        nt = self.a*np.maximum(r-self.zv, 0)
        p0 = self.t0**2*nt/(self.nu+nt)**2
        u = np.array([self.UH, 0], dtype=LD) if velocity is None else np.asarray(velocity, dtype=LD)
        f = np.asarray(forcing, dtype=LD)
        grad = (u[:, None]/self.J0+(np.maximum(r-self.h, 0)-self.J1/self.J0)*f[:, None])/(self.nu+nt)
        production = nt*np.sum(grad*grad, axis=0) if supplied is None else np.full_like(r, supplied)
        ea = self.k**LD("1.5")/(LD("2.55")*np.maximum(r, self.ze))
        ec = self.A/np.maximum(r, self.rc)
        dissipation = (1-self.chi)*ea+self.chi*ec
        redistribution = self.chi*(ec-p0)
        source = production-(1-self.chi)*ea-self.chi*p0
        return np.array([production, p0, dissipation, redistribution, source])

    def flux(self, r):
        return self.chi*(self.cap_primitive(r, self.A, self.rc)
                         -self.production_primitive(r))/self.jacobian(r)

    def quadrature(self, order=64):
        # This rule intentionally differs from the native Newton GL rule and
        # factor-two panels; none of the exported native weights are reused.
        breaks = {LD(0), self.H}
        for r in (self.h, self.zv, self.ze, self.rm, self.rc):
            if 0 < r < self.H:
                breaks.add(r)
        r = min(self.nu/self.a, self.ze, self.rm, self.rc)
        while r < self.H:
            breaks.add(r)
            r *= 3
        x, w = (np.asarray(v, dtype=LD) for v in np.polynomial.legendre.leggauss(order))
        points, weights = [], []
        intervals = sorted(breaks)
        for lower, upper in zip(intervals[:-1], intervals[1:]):
            half = (upper-lower)/2
            points.extend(lower+half*(1+x))
            weights.extend(half*w)
        return np.asarray(points), np.asarray(weights)


def basis(degree, r, height, derivative=False):
    x = np.asarray(2*r/height-1, dtype=float)
    result = [np.full_like(x, -.5) if derivative else (1-x)/2,
              np.full_like(x, .5) if derivative else (1+x)/2]
    for mode in range(2, degree+1):
        n = mode-2
        p = eval_jacobi(n, 1, 1, x)
        if derivative:
            dp = 0 if n == 0 else (n+3)/2*eval_jacobi(n-1, 2, 2, x)
            result.append(-x*p/2+(1-x*x)*dp/4)
        else:
            result.append((1-x*x)*p/4)
    return np.asarray(result, dtype=LD)*(2/height if derivative else 1)


def check_case(case):
    p = Profiles(case)
    scale = float(p.Pj)
    ref = case["reference"]
    expected = [p.us, p.t0, p.UH, p.A, p.Pj, p.Bj, (p.Pj-p.Bj)/p.Bj,
                p.H*np.sqrt(p.k)/p.nu, p.rc if p.chi else 0, p.chi]
    # Componentwise scaling avoids hiding a small cutoff error behind Re_H.
    for i, target in enumerate(expected):
        near(ref[i], target, "reference_"+str(i), max(abs(float(target)), 1e-16))
    near(ref[10], 0, "metric_budget", scale)
    require((ref[11] > 0) == (p.chi > 0), "Inactive branch solved a cutoff")
    require(p.cap_primitive(p.H, p.A, p.rc/2) > p.Pj,
            "Metric cutoff budget did not decrease from its left bracket")
    require(p.cap_primitive(p.H, p.A, min(2*p.rc, p.rm)) < p.Pj,
            "Metric cutoff budget did not decrease toward its admissible tail")
    q = np.asarray(case["quadrature"], dtype=LD)
    require(q.shape[1] == 3 and np.all(q[:, 1:] > 0), "Bad physical quadrature")
    require(np.all((q[:, 0] > 0) & (q[:, 0] < p.H)), "Quadrature outside physical ray")
    near(q[:, 2], p.jacobian(q[:, 0]), "physical_jacobian", 1)
    near(np.sum(q[:, 1]*q[:, 2]), p.volume(p.H), "physical_volume", float(p.H))
    u = (ref[2], 0)
    for sample in case["samples"]:
        r = np.array([sample["r"]], dtype=LD)
        for key, f, supplied in (("zero", (0, 0), None), ("forced", (.03, .01), None),
                                 ("supplied_negative_production", (0, 0), -.05)):
            target = p.fields(r, u, f, supplied)[:, 0]
            near(sample[key], target, "profile_"+key, max(scale/float(p.H), float(np.max(np.abs(target)))))
        near(sample["flux"], p.flux(r)[0], "analytic_flux", scale)
    near(case["flux_endpoints"], [0, 0], "zero_flux_endpoints", scale)
    r, w = p.quadrature()
    W = w*p.jacobian(r)
    values = p.fields(r, u)
    fa = p.fields(r, u, (-.02, 0))
    fb = p.fields(r, u, (.03, .01))
    fc = p.fields(r, u, (.03, 0))
    # Independent primitive and quadrature agreement is itself checked.
    near(np.sum(W*values[1]), p.Pj, "analytic_production_budget", scale)
    near(np.sum(W*p.A/np.maximum(r, p.rc)), p.Pj, "analytic_cutoff_budget", scale)
    checked_modes = 0
    maximum_source = 0.0
    for native in case["bases"]:
        degree = native["degree"]
        B, dB = basis(degree, r, p.H), basis(degree, r, p.H, True)
        probes_r = (np.array([-.81, .15, .77], dtype=LD)+1)*p.H/2
        bp, dp = basis(degree, probes_r, p.H), basis(degree, probes_r, p.H, True)
        probes = np.empty((degree+1, 6), dtype=LD)
        probes[:, ::2], probes[:, 1::2] = bp, dp
        near(native["basis_probes"], probes, "native_basis_value_derivative", 1)
        target = {"source": B@(W*values[4]), "forced_a": B@(W*fa[4]),
                  "forced_b": B@(W*fb[4]), "forced_scalar": B@(W*fc[4]),
                  "mean_only": B@(W*(values[0]-values[2])),
                  "redistribution": B@(W*values[3])}
        endpoint = basis(degree, np.array([p.H]), p.H)[:, 0]*p.jacobian(p.H)*p.flux(p.H)
        target["flux_action"] = -dB@(W*p.flux(r))+endpoint
        for key, value in target.items():
            near(native[key], value, "native_dual_load_"+key, scale)
        near(native["redistribution"], native["flux_action"], "source_flux_integration_by_parts", scale)
        near(np.asarray(native["forced_b"])-native["forced_a"], B@(W*(fb[0]-fa[0])),
             "all_mode_vector_forcing_increment", scale)
        near(np.asarray(native["forced_scalar"])-native["forced_a"], B@(W*(fc[0]-fa[0])),
             "all_mode_scalar_forcing_increment", scale)
        near(native["sampled_polynomial_k"], p.k*(1+LD(".2")*p.rK/p.H+LD(".1")*(p.rK/p.H)**2),
             "physical_K_sampler", float(p.k))
        if p.chi == 1:
            require(np.all(np.asarray(native["source"]) == 0), "Reference did not cancel pointwise")
            near(sum(native["mean_only"][:2]), 0, "mean_only_zero_mean", scale)
            require(max(abs(x) for x in native["mean_only"]) > 1e-5*scale,
                    "Mean-only failing higher-mode control disappeared")
        checked_modes += degree+1
        maximum_source = max(maximum_source, max(abs(x) for x in native["source"]))
    if "quad_ray_loads" in case:
        normal = basis(4, r, p.H)@(W*fb[4])
        tangent = basis(4, np.array([(.27+1)/2]), LD(1))[:, 0]
        near(case["quad_ray_loads"], [tangent[i % 5]*normal[i // 5] for i in range(25)],
             "all_native_quad_ray_rows", scale)
        checked_modes += 25
    moments = None
    if "legendre_mean_only" in case:
        x = np.asarray(2*r/p.H-1, dtype=float)
        moments = [float(np.sum(W*eval_legendre(i, x)*(values[0]-values[2]))) for i in range(7)]
        near(case["orthonormal_mean_only"], [moments[i]*math.sqrt((2*i+1)/2) for i in range(7)],
             "native_orthonormal_control", scale)
        near(case["legendre_mean_only"], moments, "native_legendre_control", scale)
        checked_modes += 7
        note = [0, .0173551439079, -.0220336637626, .0220527999939,
                -.0197059614435, .0160477076949, -.0116765463336]
        near(moments, note, "published_mean_only_control", 1, 2e-12)
    if case["name"] == "note":
        # The note's stronger dimensional tolerances, beyond the general scaled target.
        for native in case["bases"]:
            near(native["redistribution"], native["flux_action"], "note_dimensional_redistribution", 1, 2e-12)
            B = basis(native["degree"], r, p.H)
            near(np.asarray(native["forced_scalar"])-native["forced_a"], B@(W*(fc[0]-fa[0])),
                 "note_dimensional_forcing_increment", 1, 3e-12)
    if case["name"] == "curved_note" and p.chi:
        near(ref[8], .009458614009, "published_curved_cutoff", 1, 5e-13)
        flat_rc = p.H*np.exp(1-p.Pj/p.A)
        require(abs(flat_rc-p.rc) > LD(".001")*p.rc, "Flat-cutoff rejection control lacks separation")
    return {"name": case["name"], "order": case["order"], "chi": case["chi"],
            "native_modes": checked_modes, "physical_points": len(q), "oracle_points": len(r),
            "cutoff": ref[8], "production_budget": ref[4], "minimum_budget": ref[5],
            "budget_residual": ref[10], "max_abs_reference_source_load": maximum_source,
            "legendre_mean_only": moments}


def refinement_and_independence(cases):
    lookup = {(c["name"], c["order"], c["chi"]): c for c in cases}
    pairs = 0
    for c in cases:
        if c["order"] != 16:
            continue
        d = lookup[c["name"], 32, c["chi"]]
        scale = c["reference"][4]
        near(c["reference"][8], d["reference"][8], "quadrature_refinement_cutoff",
             max(c["reference"][8], 1e-16))
        for a, b in zip(c["bases"], d["bases"]):
            for key in ("source", "forced_a", "forced_b", "forced_scalar", "mean_only", "redistribution", "flux_action"):
                near(a[key], b[key], "quadrature_refinement_"+key, scale)
        pairs += 1
    for order in (16, 32):
        for chi in (0, .37, 1):
            c = lookup["note", order, chi]
            for name in ("sample_at_H", "tail_at_H"):
                other = lookup[name, order, chi]
                near(c["reference"][8], other["reference"][8], "independent_rK_rm_cutoff",
                     max(c["reference"][8], 1e-16))
                for a, b in zip(c["bases"], other["bases"]):
                    near(a["forced_b"], b["forced_b"], "independent_rK_rm_load", c["reference"][4])
    return {"order_16_32_pairs": pairs, "native_segment_degrees": [2, 4, 6, 8],
            "weights": [0, .37, 1], "fixed_physical_sampler_including_H": "passed",
            "independent_tail_start_including_H": "passed"}


def controls(raw):
    require(raw["fallbacks"] == [2, 0, 0, 3], "Wrong inactive/inadmissible branch")
    def gate(x, lo, hi):
        s = min(1, max(0, (x-lo)/(hi-lo)))
        return 1-10*s**3+15*s**4-6*s**5
    for spatial, sensor, margin, reynolds, actual in raw["switch_samples"]:
        expected = spatial*gate(sensor*sensor, .02**2, .1**2)*(1-gate(margin, .05, .1))*(1-gate(reynolds, 100, 200))
        near(actual, expected, "independent_C2_gate", 1, 1e-13)
    require(raw["cancelled_forcing_weight"] == 0, "Independent forcing sensors cancelled")
    range_control = raw["range_controls"]
    extreme = copy.deepcopy(raw["cases"][0])
    extreme["chi"] = 1
    extreme["geometry"][0] = range_control["subnormal_nu"]
    extreme["geometry"][2] = 0
    p = Profiles(extreme)
    near(range_control["cutoff"], p.rc, "representable_subnormal_cutoff", float(p.rc))
    near(range_control["cap"], p.A/p.rc, "representable_subnormal_cap", float(p.A/p.rc))
    near(range_control["production_budget"], p.Pj, "subnormal_production_budget", float(p.Pj))
    near(range_control["budget_residual"], 0, "subnormal_metric_budget", float(p.Pj))
    require(0 < range_control["cutoff"] < sys.float_info.min, "No subnormal-range control")
    extreme["geometry"][0] = range_control["overflow_nu"]
    overflowing = Profiles(extreme)
    require(overflowing.A/overflowing.rc > sys.float_info.max, "Overflow control is representable")
    require(range_control["overflow_weight"] == 0 and range_control["overflow_reason"] == 4,
            "Unrepresentable cap did not explicitly fall back")
    c = raw["positivity_counterexample"]
    require(c["sampled_k"] > 0 and c["r0"] != c["rK"], "No distinct positive K sampler")
    require([c[k] for k in ("analytic_value", "analytic_first", "analytic_second", "wall_derivative")] == [0, 0, 0, 0],
            "Counterexample is not a flat nonnegative zero")
    for key in ("polynomial_value", "polynomial_first", "polynomial_second", "native_wall_derivative"):
        near(c[key], 0, "native_counterexample_"+key, c["sampled_k"], 2e-10)
    require(c["production"] == 0 and c["dissipation"] > 0, "Wrong local viscous-sublayer balance")
    expected_ea = c["sampled_k"]**2/(2.55*5.1*1.5e-5)
    near(c["dissipation"], expected_ea, "independent_counterexample_sink", expected_ea)
    near(c["source"], -expected_ea, "counterexample_source", expected_ea)
    near(c["rhs"], c["source"]+1.5e-5*c["polynomial_second"], "counterexample_local_rhs", expected_ea)
    near(c["local_euler_value"], c["polynomial_value"]+c["dt"]*c["rhs"], "counterexample_local_step", expected_ea*c["dt"])
    require(c["rhs"] < 0 and c["local_euler_value"] < 0 and not c["full_timestep"],
            "Expected local positivity failure was hidden or overclaimed")
    return {"native_rejected_inputs": raw["summary"]["rejected_inputs"],
            "floating_point_range": range_control,
            "inactive_cutoff_iterations": raw["fallbacks"][2],
            "cancelling_pressure_convection_weight": raw["cancelled_forcing_weight"],
            "physical_inner_k_positivity_gate": "failed_as_specified",
            "counterexample": c, "native_full_solver_timestep": "not_tested_M4"}


def rejection_control(case):
    mutated = copy.deepcopy(case)
    mutated["bases"][0]["source"][0] += .001*case["reference"][4]
    saved = errors.copy()
    try:
        try:
            check_case(mutated)
        except RuntimeError as error:
            require("native_dual_load_source" in str(error), "Wrong oracle rejection")
            return "corrupted_native_source_load_rejected"
        raise RuntimeError("Independent oracle accepted a corrupted source load")
    finally:
        errors.clear()
        errors.update(saved)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", required=True, type=Path)
    parser.add_argument("--output", type=Path, help="Fresh output directory; default beside probe")
    args = parser.parse_args()
    require(sys.platform.startswith("linux"), "Authoritative numerical validation requires Linux")
    root = Path(__file__).resolve().parents[4]
    require(digest(root/NOTE) == NOTE_SHA and digest(root/CONTRACT) == CONTRACT_SHA, "Frozen specification changed")
    output = args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix="cls-epawf-m2-", dir=args.probe.resolve().parent))
    if args.output:
        output.mkdir(parents=True, exist_ok=False)
    report = {"schema": 1, "milestone": "cls_epawf_revised_m2", "verification_status": "running",
              "utc": datetime.now(timezone.utc).isoformat(), "host": platform.node(),
              "python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
              "cpu_affinity": sorted(os.sched_getaffinity(0)), "mpi_ranks": 1,
              "thread_environment": THREADS, "commands": []}
    try:
        command = [str(args.probe.resolve()), "--output", str(output/"profiles.json")]
        start = time.monotonic()
        result = subprocess.run(command, cwd=output, capture_output=True, text=True, timeout=180)
        (output/"native.log").write_text(result.stdout+result.stderr)
        report["commands"].append({"command": command, "exit_code": result.returncode,
                                   "seconds": time.monotonic()-start, "binary_sha256": digest(args.probe)})
        require(result.returncode == 0, "Native profile probe failed: "+str(output/"native.log"))
        raw = json.loads((output/"profiles.json").read_text())
        require(len(raw["cases"]) == 62, "Incomplete native matrix")
        report["cases"] = [check_case(c) for c in raw["cases"]]
        report["refinement"] = refinement_and_independence(raw["cases"])
        report["controls"] = controls(raw)
        report["oracle_rejection_control"] = rejection_control(raw["cases"][0])
        report["native_summary"] = raw["summary"]
        report["verification_status"] = "passed"
        report["scope"] = "local spectral ray profiles and dual loads; full-strain production is caller-supplied; global field layout, patch, curved mesh, time-map and robust-model qualification remain later milestones"
    except Exception as error:
        report["verification_status"] = "failed"
        report["error"] = str(error)
        raise
    finally:
        report["checks"] = checks
        report["scaled_errors"] = errors
        report["source_commit"] = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
        report["source_status"] = subprocess.check_output(["git", "-C", str(root), "status", "--porcelain"], text=True).splitlines()
        sources = [NOTE, CONTRACT, "library/UnitTests/SolverUtils/RANSEquilibriumProfilesProbe.cpp",
                   "library/SolverUtils/CMakeLists.txt", "solvers/IncNavierStokesSolver/CMakeLists.txt"]
        sources += ["library/SolverUtils/Turbulence/"+p for p in (
            "RANSModels.h", "RANSModels.cpp", "AnalyticalWallFunction.h", "AnalyticalWallFunction.cpp",
            "EquilibriumAWFProfiles.h", "EquilibriumAWFProfiles.cpp")]
        sources += ["solvers/IncNavierStokesSolver/Tests/RANS/"+p for p in ("validate_cls_epawf_m2.py", "CMakeLists.txt")]
        report["source_sha256"] = {p: digest(root/p) for p in sources}
        report["artifacts_sha256"] = {p.name: digest(p) for p in output.iterdir() if p.is_file()}
        (output/"validation.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"verification_status": report["verification_status"], "checks": checks,
                      "physical_inner_k_positivity_gate": report["controls"]["physical_inner_k_positivity_gate"],
                      "report": str(output/"validation.json")}, indent=2))


if __name__ == "__main__":
    main()
