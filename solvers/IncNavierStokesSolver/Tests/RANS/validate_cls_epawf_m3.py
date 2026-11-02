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
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile
import time

sys.dont_write_bytecode = True
THREADS = {k: "1" for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS")}
os.environ.update(THREADS)
import numpy as np
import scipy
from scipy.special import eval_jacobi
from validate_cls_epawf_m2 import Profiles
from validate_cls_epawf_m1 import indexed

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
    near(raw["merged_shared_load"], np.zeros(raw["layout"][0]["global_coefficients"]), "opposite_common_flux_dual_loads")
    return {"pressure_rows": B.shape[0], "pressure_rank": rank,
            "pressure_singular_values": singular.tolist(), "interface_rows": C.shape[0],
            "interface_rank": C.shape[0], "augmented_rows": len(rhs),
            "metrics": pressure["metrics"]}


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
        require(args.suite == "residuals", "Full M3 qualification is not implemented yet; residuals is a development pilot")
        report["verification_status"] = "passed"
        report["milestone_closure"] = False
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
                    "solvers/IncNavierStokesSolver/Tests/RANS/EquilibriumWallFlat.xml"]
        report["source_sha256"] = {name: digest(root/name) for name in sources}
        report["artifacts_sha256"] = {p.name: digest(p) for p in output.iterdir() if p.is_file()}
        (output/"validation.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"verification_status": report["verification_status"], "checks": checks,
                      "milestone_closure": report["milestone_closure"], "report": str(output/"validation.json")}, indent=2))


if __name__ == "__main__":
    main()
