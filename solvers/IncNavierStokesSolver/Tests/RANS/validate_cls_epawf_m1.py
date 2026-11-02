#!/usr/bin/env python3
"""Linux M1 audit of the original CLS and scalar AWF, with independent oracles.

The indexed/quadrature checks reuse the earlier independent M2 oracle.
New complex-step checks differentiate a dimensional indexed tensor, whereas
the native probe samples the production dimensionless matrix implementation.
Passing the counterexample tests records a FAILED physical stability gate.
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
from scipy.integrate import quad
import validate_m2 as oracle

NOTE = "docs/literature/cls_neutral_abl_equilibrium_awf_note_revised.tex"
NOTE_SHA256 = "d35e437c2f2f8b25a46f7a3f285c6f1ad8b48d7b554aa1d161315f59eeb6de0f"
errors = {}
checks = 0


def require(condition, message):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(message)


def near(actual, expected, name, tolerance=1e-11, scale=1.0):
    a, b = np.asarray(actual), np.asarray(expected)
    require(np.all(np.isfinite(a)) and np.all(np.isfinite(b)), "Nonfinite " + name)
    error = float(np.max(np.abs(a-b))) / max(scale, float(np.max(np.abs(b))))
    errors[name] = max(errors.get(name, 0.0), error)
    require(error <= tolerance, f"{name}: scaled error {error} > {tolerance}")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def indexed(gradient, k, epsilon, branch):
    """Dimensional index form; no production routines or finite differences."""
    g = np.asarray(gradient, dtype=complex)
    s, o, eye = g+g.T, g-g.T, np.eye(3)
    ss, oo = np.sum(s*s), np.sum(o*o)
    qs, qo = k/epsilon*np.sqrt(ss/2), k/epsilon*np.sqrt(oo/2)
    q = qs if branch == "strain" else qo if branch == "rotation" else (
        qs if qs.real >= qo.real else qo)
    with np.errstate(over="ignore", under="ignore"):
        c = .3/(1+.35*q**1.5)*(-np.expm1(-.36*np.exp(.75*q)))
    nt = c*k*k/epsilon
    t1 = np.einsum("ik,jk->ij", s, s)-ss*eye/3
    t2 = np.einsum("ik,jk->ij", o, s)+np.einsum("jk,ik->ij", o, s)
    t3 = np.einsum("ik,jk->ij", o, o)-oo*eye/3
    t4 = np.einsum("ki,lj,kl->ij", s, o, s)+np.einsum("kj,li,kl->ij", s, o, s)
    quadratic = k*nt/epsilon*(-.1*t1+.1*t2+.26*t3)
    cubic = nt*k*k/epsilon**2*(-10*c*c*t4+5*c*c*(oo-ss)*s)
    stress = 2*k*eye/3-nt*s+quadratic+cubic
    return {"stress": stress, "quadratic": quadratic, "cubic": cubic,
            "cMu": c, "nut": nt, "strain": qs, "rotation": qo,
            "production": -np.sum(stress*g)}


def check_state(case, state):
    k, ep, branch = case["k"], case["epsilon"], case["branch"]
    expected = indexed(state["gradient"], k, ep, branch)
    g = np.asarray(state["gradient"])
    for key, value in expected.items():
        near(state[key], value.real, "native_indexed_"+key, scale=k)
    near(state["flux"], case["nu"]*(g+g.T)-expected["stress"].real,
         "native_total_flux", scale=k)
    ratio = k**1.5/(2.5*ep*case["distance"])
    yap = ep*ep/k*.83*max(ratio-1, 0)*ratio*ratio
    near(state["yap"], yap, "native_geometric_yap", scale=ep)
    near(state["k_source"], expected["production"].real-ep,
         "native_k_source", scale=ep)
    near(state["epsilon_source"], 1.44*ep/k*expected["production"].real
         -1.92*ep*ep/k+yap, "native_epsilon_source", scale=ep*ep/k)


def tangents(report):
    results = []
    for case in report["cases"]:
        check_state(case, case["base"])
        g = np.asarray(case["base"]["gradient"])
        direction = np.asarray(case["direction"])
        near(np.sum(direction*direction), 1, "unit_rank_one_direction")
        near(np.trace(direction), 0, "incompressible_direction")
        hcs = 1e-30
        analytic = indexed(g+1j*hcs*direction, case["k"], case["epsilon"],
                           case["branch"])
        dr = analytic["stress"].imag/hcs
        dt = case["nu"]*(direction+direction.T)-dr
        dp = analytic["production"].imag/hcs
        forms, action_errors = [], []
        for sample in case["samples"]:
            step = sample["step"]
            for sign, side in ((1, "plus"), (-1, "minus")):
                check_state(case, sample[side])
                near(sample[side]["gradient"], g+sign*step*direction,
                     "probe_direction", tolerance=2e-15)
            # Actual physical smooth-branch differences must stay on that
            # branch. Formal tie derivatives are labelled separately.
            gap = case["base"]["strain"]-case["base"]["rotation"]
            if case["branch"] == "physical" and abs(gap) > 1e-12:
                for side in ("plus", "minus"):
                    child_gap = sample[side]["strain"]-sample[side]["rotation"]
                    require(gap*child_gap > 0, "Finite difference crossed invariant tie")
            actual = (np.asarray(sample["plus"]["flux"])-
                      np.asarray(sample["minus"]["flux"]))/(2*step)
            # Central-difference truncation plus the amplification of rounded
            # native stress samples. The complex-step oracle has no subtractive
            # cancellation; the native near-tie samples necessarily do.
            stress_scale = max(1, float(np.max(np.abs(case["base"]["flux"]))))
            bound = (1e-8*max(1, float(np.max(np.abs(dt))))+50*step*step
                     +16*np.finfo(float).eps*stress_scale/step)
            error = float(np.max(np.abs(actual-dt)))
            require(error <= bound, f"Tangent action mismatch in {case['name']}: {error}")
            forms.append(float(np.sum(direction*actual)))
            action_errors.append(error)
        require(min(action_errors) < 1e-6, "No resolved tangent difference")
        results.append({"name": case["name"], "branch": case["branch"],
                        "gradient": g.tolist(), "principal_form": float(np.sum(direction*dt)),
                        "native_forms": forms, "steps": [s["step"] for s in case["samples"]],
                        "minimum_action_error": min(action_errors),
                        "production_derivative": float(dp),
                        "reference_covariance_minimum": float(np.linalg.eigvalsh(
                            np.asarray(case["base"]["stress"])).min())})
    physical = {c["name"]: c for c in report["cases"] if c["branch"] == "physical"}
    block = np.zeros((2, 2))
    for j, name in enumerate(("transverse_x", "transverse_z")):
        sample = physical[name]["samples"][-1]
        derivative = (np.asarray(sample["plus"]["flux"])-
                      np.asarray(sample["minus"]["flux"]))/(2*sample["step"])
        block[:, j] = derivative[1, [0, 2]]
    q, c = report["q"], report["cMu"]
    expected = np.array([[1-20*c, -.36*q], [.36*q, 1.]])
    near(block, expected, "transverse_flux_block", tolerance=3e-9)
    require(block[0, 0] < -.96, "Adverse original CLS diffusion was removed")
    for row in results:
        if row["name"] == "transverse_x":
            require(row["reference_covariance_minimum"] > 1,
                    "Counterexample covariance lost realizability")
            near(row["production_derivative"], 0, "transverse_production_derivative")
        if row["name"] == "nearby_smooth_inplane":
            require(all(value < -.02 for value in row["native_forms"]),
                    "Nearby smooth 2D adverse direction was not reproduced")
            near(row["principal_form"], -.02596757790,
                 "nearby_smooth_inplane_form", tolerance=1e-10)
    # Both limiting tangents are independently exercised and are different.
    formal = {r["branch"]: r["principal_form"] for r in results
              if r["name"] == "reference_inplane"}
    require(abs(formal["strain"]-formal["rotation"]) > .01,
            "The two invariant-branch tangents were not distinguished")
    for name, branch in (("strain_limit", "strain"), ("rotation_limit", "rotation")):
        sequence = [r for r in results if r["name"] == name]
        require(len(sequence) == 2, "Missing physical limiting sequence")
        require(abs(sequence[-1]["principal_form"]-formal[branch]) <
                abs(sequence[0]["principal_form"]-formal[branch]),
                "Physical branch sequence does not approach its formal tangent")
    molecular = next(r for r in results if r["name"] == "transverse_molecular_control")
    near(molecular["principal_form"], 1-20*c+.025, "molecular_tangent_shift")
    return {"physical_stability_gate": "failed_as_specified",
            "transverse_flux_block": block.tolist(), "cases": results,
            "scope": "frozen positive k/epsilon; specified local directions, not an all-state certificate"}


def profiles(report):
    q, kap, cmu = report["q"], report["kappa"], report["cMu"]
    near(q, 3.190675081624288, "q_star")
    near(cmu, .0982278953282113, "Cmu_star")
    near(kap, .440553905186659, "kappa_c")
    covariance = np.array([[.92, 0, -1/q], [0, .56, 0], [-1/q, 0, .52]])
    expected_eigenvalues = [0.3482098773, .56, 1.0917901227]
    eigenvalues = np.linalg.eigvalsh(covariance)
    near(eigenvalues, expected_eigenvalues, "covariance_eigenvalues", tolerance=5e-11)
    ry = kap*q**1.5/2.5
    hy = .83*max(ry-1, 0)*ry*ry
    near(q*kap*kap/1.3+hy, .48, "kappa_yap_balance")
    require(hy > .0036, "Weakly active reference Yap correction was omitted")
    for p in report["reference_profiles"]:
        z, k, ep, u = p["z"], p["k"], p["epsilon"], .4
        case = {"k": k, "epsilon": ep, "branch": "physical", "distance": z, "nu": 0}
        check_state(case, p["state"])
        near(k, q*u*u, "constant_reference_k")
        near(ep, u**3/(kap*z), "reference_epsilon")
        near(p["velocity"], 5.7+u/kap*math.log(z/2), "reference_velocity")
        near(p["state"]["stress"], k*covariance, "reference_full_stress", scale=k)
        near(p["state"]["production"], ep, "reference_production", scale=ep)
        near(p["state"]["k_source"], 0, "reference_k_balance", scale=ep)
        near(p["state"]["epsilon_source"]+u**4/(1.3*z*z), 0,
             "reference_epsilon_balance", scale=ep*ep/k)
        near(p["epsilon_gradient"], -ep/z, "reference_epsilon_gradient", scale=ep/z)
        near(p["epsilon_gradient_flux"], -u**4/(1.3*z),
             "reference_epsilon_flux", scale=u**4/z)
    original = next(c for c in report["cases"] if c["name"] == "transverse_x" and c["branch"] == "physical")
    changed = next(c for c in report["cases"] if c["name"] == "same_state_double_wall_distance")
    near(original["base"]["stress"], changed["base"]["stress"], "distance_independent_stress")
    require(original["base"]["yap"] > changed["base"]["yap"] > 0,
            "Yap did not use the changed geometric wall distance")
    return {"q": q, "Cmu": cmu, "kappa": kap, "yap_ratio": ry,
            "yap_function": hy, "covariance_eigenvalues": eigenvalues.tolist(),
            "heights": [p["z"] for p in report["reference_profiles"]]}


def wall_balances(kernels):
    cases = kernels["walls"]+kernels["standard_walls"]
    energy_errors, reversed_cases = [], 0
    for c in cases:
        r, g = c["reference"], c["geometry"]
        # Each native integral/member has already been checked by independent
        # quadrature. This checks their combined reference energy bookkeeping.
        power = r["traction"]*r["U_H"]
        balance = power-g[1]*r["P0"]-g[0]*r["traction"]**2*c["integrals"][2]
        near(balance, 0, "reference_wall_energy", scale=power)
        energy_errors.append(abs(balance)/power)
        reversed_cases += any(v < 0 for v in c["original"]["traction"])
    require(reversed_cases > 0, "Signed/reversed wall traction was not exercised")
    note = [c for c in kernels["walls"] if c["name"] == "note"]
    require(len(note) == 4, "Missing illustrative forcing cases")
    by_forcing = {c["forcing"][0]: c for c in note}
    a, b = by_forcing[-.02], by_forcing[.03]
    j0, j1 = a["integrals"][:2]
    slope = (b["original"]["traction"][0]-a["original"]["traction"][0])/.05
    near(slope, -j1/j0, "finite_wall_forcing_slope")
    H, rough = a["geometry"][1:3]
    interface_increment = ((b["original"]["traction"][0]+.03*(H-rough))
                          -(a["original"]["traction"][0]-.02*(H-rough)))
    near(interface_increment/.05, H-rough-j1/j0, "finite_interface_forcing_slope")
    # Independent production polynomial for the published dimensional example.
    k, nu = a["k"], a["geometry"][0]
    alpha = .09*2.55
    zv = a["intercept"]*nu/math.sqrt(k)  # intercept independently audited above
    def nt(z):
        return alpha*math.sqrt(k)*max(0, z-zv)
    def integrate(f):
        return quad(f, 0, rough, epsabs=1e-13, epsrel=3e-13)[0]+quad(
            f, rough, H, epsabs=1e-13, epsrel=3e-13)[0]
    B = j1/j0
    moments = [integrate(lambda z, n=n: nt(z)/(nu+nt(z))**2*
                         (max(z-rough, 0)-B)**n)/H for n in range(3)]
    for c in note:
        u, f = np.asarray(c["velocity"])/j0, np.asarray(c["forcing"])
        expected = moments[0]*np.dot(u, u)+2*moments[1]*np.dot(u, f)+moments[2]*np.dot(f, f)
        near(c["original"]["production"], expected, "finite_production_polynomial", scale=expected)
    return {"wall_cases": len(cases), "signed_reversal_cases": reversed_cases,
            "maximum_reference_energy_relative_error": max(energy_errors),
            "wall_forcing_slope_m": slope,
            "interface_forcing_slope_m": interface_increment/.05,
            "production_polynomial": moments,
            "negative_virtual_intercept_verified": a["zv"] < 0,
            "note_reference": by_forcing[0]["reference"]}


def rejection_controls(kernels):
    """Ensure the independent checker rejects tensor and wall-sign defects."""
    saved = dict(oracle.errors)
    detected = []
    try:
        altered = copy.deepcopy(next(c for c in kernels["bulk"] if c["model"] == "cls"))
        altered["stress"][0][0] += .01*altered["k"]
        try:
            oracle.check_bulk(altered)
        except RuntimeError:
            detected.append("incorrect_native_stress")
        wall = copy.deepcopy(kernels["walls"][0])
        wall["original"]["traction"][0] *= -1
        try:
            oracle.check_wall(wall, wall["q"], wall["kappa"])
        except RuntimeError:
            detected.append("incorrect_wall_traction_sign")
    finally:
        oracle.errors.clear()
        oracle.errors.update(saved)
    require(len(detected) == 2, "An independent oracle accepted a deliberately corrupted output")
    return detected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kernel-probe", required=True, type=Path)
    parser.add_argument("--constitutive-probe", required=True, type=Path)
    parser.add_argument("--output", type=Path, help="Fresh directory; default keeps a new directory beside the probe")
    args = parser.parse_args()
    require(sys.platform.startswith("linux"), "Authoritative numerical validation requires Linux")
    root = Path(__file__).resolve().parents[4]
    require(digest(root/NOTE) == NOTE_SHA256, "The governing revised note changed")
    output = args.output.resolve() if args.output else Path(tempfile.mkdtemp(
        prefix="cls-epawf-m1-", dir=args.kernel_probe.resolve().parent))
    if args.output:
        output.mkdir(parents=True, exist_ok=False)
    report = {"schema": 1, "milestone": "cls_epawf_revised_m1",
              "verification_status": "running", "utc": datetime.now(timezone.utc).isoformat(),
              "host": platform.node(), "python": platform.python_version(),
              "numpy": np.__version__, "scipy": scipy.__version__,
              "cpu_affinity": sorted(os.sched_getaffinity(0)), "mpi_ranks": 1,
              "thread_environment": THREADS, "commands": []}
    try:
        for name, binary in (("kernels", args.kernel_probe), ("constitutive", args.constitutive_probe)):
            command = [str(binary.resolve()), "--output", str(output/(name+".json"))]
            print("M1 running "+name, flush=True)
            start = time.monotonic()
            result = subprocess.run(command, cwd=output, capture_output=True, text=True, timeout=60)
            (output/(name+".log")).write_text(result.stdout+result.stderr)
            report["commands"].append({"command": command, "exit_code": result.returncode,
                                       "seconds": time.monotonic()-start,
                                       "binary_sha256": digest(binary)})
            require(result.returncode == 0, f"{name} failed; see {output}/{name}.log")
        kernels = json.loads((output/"kernels.json").read_text())
        constitutive = json.loads((output/"constitutive.json").read_text())
        oracle.errors.clear()
        for case in kernels["references"]:
            oracle.reference(case)
        for case in kernels["bulk"]:
            oracle.check_bulk(case)
        for case in kernels["walls"]+kernels["standard_walls"]:
            oracle.check_wall(case, case["q"], case["kappa"])
        require(max(oracle.errors.values()) <= 1e-11, "Indexed/quadrature error exceeds the approved 1e-11 target")
        report["reference"] = profiles(constitutive)
        report["constitutive"] = tangents(constitutive)
        report["awf"] = wall_balances(kernels)
        report["oracle_rejection_controls"] = rejection_controls(kernels)
        report["native_kernel_summary"] = kernels["summary"]
        report["independent_oracle_errors"] = {k: float(v) for k, v in oracle.errors.items()}
        report["additional_scaled_errors"] = dict(errors)
        report["checks"] = checks
        report["verification_status"] = "passed"
        report["scope"] = "original high-Re CLS and scalar AWF only; no wall-patch assembly, native time-map, positivity or terrain qualification"
    except Exception as error:
        report["verification_status"] = "failed"
        report["error"] = str(error)
        raise
    finally:
        report["source_commit"] = subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
        report["source_status"] = subprocess.check_output(
            ["git", "-C", str(root), "status", "--porcelain"], text=True).splitlines()
        sources = [NOTE, "library/UnitTests/SolverUtils/RANSKernelProbe.cpp",
                   "library/UnitTests/SolverUtils/RANSConstitutiveProbe.cpp",
                   "solvers/IncNavierStokesSolver/CMakeLists.txt"]
        sources += ["library/SolverUtils/Turbulence/"+p for p in (
            "RANSModels.h", "RANSModels.cpp", "AnalyticalWallFunction.h",
            "AnalyticalWallFunction.cpp", "RANSWallStrip.h")]
        sources += ["solvers/IncNavierStokesSolver/Tests/RANS/"+p for p in (
            "validate_cls_epawf_m1.py", "validate_m2.py", "validate_m1.py", "CMakeLists.txt")]
        report["source_sha256"] = {p: digest(root/p) for p in sources}
        report["artifacts_sha256"] = {p.name: digest(p) for p in output.iterdir() if p.is_file()}
        (output/"validation.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"verification_status": report["verification_status"],
                      "physical_stability_gate": report["constitutive"]["physical_stability_gate"],
                      "checks": checks, "report": str(output/"validation.json")}, indent=2))


if __name__ == "__main__":
    main()
