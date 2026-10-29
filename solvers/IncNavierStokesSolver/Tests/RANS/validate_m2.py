#!/usr/bin/env python3
"""Independent index-form and quadrature audit of the production RANS kernels.

Run on Linux. NumPy/SciPy implement the independent oracle; the C++ probe uses
closed analytic wall moments and native point kernels. No external services.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys

sys.dont_write_bytecode = True
# Set these before importing NumPy/SciPy: BLAS may create its pool on import.
THREAD_ENVIRONMENT = {name: "1" for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                      "MKL_NUM_THREADS", "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")}
os.environ.update(THREAD_ENVIRONMENT)
if sys.platform.startswith("linux"):
    os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[-16:])
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from validate_m1 import validate as validate_column

errors = {}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def near(actual, expected, name, rtol=3e-11, floor=1e-100):
    a, b = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
    require(np.all(np.isfinite(a)) and np.all(np.isfinite(b)), f"Nonfinite {name}")
    scale = max(floor, float(np.max(np.abs(a))), float(np.max(np.abs(b))))
    relative = float(np.max(np.abs(a-b))) / scale
    errors[name] = max(errors.get(name, 0.0), relative)
    require(relative <= rtol, f"{name}: scaled error {relative} exceeds {rtol}")


def cmu(q):
    # Direct published expression; overflow in the nested exponential has the
    # well-defined limiting damping factor 1, unlike a NaN or a clipped q.
    q = np.longdouble(q)
    with np.errstate(over="ignore"):
        return np.longdouble(.3)/(1+np.longdouble(.35)*q**np.longdouble(1.5)) * (
            1-np.exp(-np.longdouble(.36)*np.exp(np.longdouble(.75)*q)))


def check_bulk(case):
    k, epsilon = np.longdouble(case["k"]), np.longdouble(case["epsilon"])
    p = case["parameters"]
    g = np.asarray(case["gradient"], dtype=np.longdouble)
    s, o, eye = g+g.T, g-g.T, np.eye(3, dtype=np.longdouble)
    ss, oo = np.sum(s*s), np.sum(o*o)
    strain, rotation = k/epsilon*np.sqrt(ss/2), k/epsilon*np.sqrt(oo/2)
    coefficient = cmu(max(strain, rotation)) if case["model"] == "cls" else np.longdouble(p[0])
    nut = coefficient*k*k/epsilon
    q2 = np.zeros((3, 3), dtype=np.longdouble)
    q3 = q2.copy()
    if case["model"] == "cls":
        c1, c2, c3 = np.longdouble("-.1"), np.longdouble(".1"), np.longdouble(".26")
        c4, c5, c6, c7 = -10*coefficient**2, 0, -5*coefficient**2, 5*coefficient**2
        # Dimensional, uncontracted index form from the original relation.
        # This deliberately differs from the kernel's dimensionless matrix law.
        t1 = np.einsum("ik,jk->ij", s, s)-ss*eye/3
        t2 = np.einsum("ik,jk->ij", o, s)+np.einsum("jk,ik->ij", o, s)
        t3 = np.einsum("ik,jk->ij", o, o)-oo*eye/3
        t4 = np.einsum("ki,lj,kl->ij", s, o, s)+np.einsum("kj,li,kl->ij", s, o, s)
        t5 = o@o@s+s@o@o-2*np.trace(s@o@o)*eye/3
        q2 = k*nut/epsilon*(c1*t1+c2*t2+c3*t3)
        q3 = nut*k*k/epsilon**2*(c4*t4+c5*t5+(c6*ss+c7*oo)*s)
    stress = 2*k*eye/3-nut*s+q2+q3
    production = -np.sum(stress*g)
    ratio = k**np.longdouble(1.5)/(p[6]*epsilon*case["distance"])
    yap = epsilon**2/k*p[5]*max(ratio-1, 0)*ratio**2 if case["model"] == "cls" else 0
    ep = p[1]*epsilon/k*production
    ed = p[2]*epsilon**2/k
    near(case["stress"], stress, "tensor_stress", floor=float(k))
    near(case["quadratic"], q2, "tensor_quadratic", floor=float(k)*1e-12)
    near(case["cubic"], q3, "tensor_cubic", floor=float(k)*1e-12)
    for key, expected in (("cMu", coefficient), ("nut", nut), ("strain", strain),
                          ("rotation", rotation), ("divergence", np.trace(g)),
                          ("production", production), ("k_source", production-epsilon),
                          ("epsilon_production", ep), ("epsilon_destruction", ed),
                          ("yap", yap), ("epsilon_source", ep-ed+yap)):
        near(case[key], expected, "bulk_"+key, floor=1e-13 if key == "divergence" else 1e-100)
    near(case["diffusion"], [case["nu"]+nut, case["nu"]+nut/p[3], case["nu"]+nut/p[4]],
         "diffusion")


def reference(case):
    p = case["parameters"]
    if case["model"] == "cls":
        q = brentq(lambda x: float(cmu(x)*x*x)-1, 1, 8, xtol=5e-15)
        def function(kappa):
            r = kappa*q**1.5/p[6]
            return q*kappa*kappa/p[4]+p[5]*max(r-1, 0)*r*r-(p[2]-p[1])
        kappa = brentq(function, 0, 2, xtol=5e-15)
        coefficient = cmu(q)
    else:
        q = 1/math.sqrt(p[0])
        kappa = math.sqrt((p[2]-p[1])*p[4]/q)
        coefficient = p[0]
    near(case["q"], q, "reference_q")
    near(case["kappa"], kappa, "reference_kappa")
    near(case["cMu"], coefficient, "reference_cMu")
    return q, kappa


def check_wall(case, qstar, kappa):
    nu, height, roughness, point_p = case["geometry"]
    k = case["k"]
    c_mu, c_length, smooth = case["parameters"]
    root_k = np.sqrt(np.longdouble(k))
    alpha = np.longdouble(c_mu)*c_length
    if roughness == 0:
        intercept = smooth
    else:
        x = roughness*root_k/(70*np.longdouble(nu))
        exponent = max(np.longdouble(.5)-np.longdouble(.4)*x**np.longdouble(.7),
                       1-np.longdouble(.79)*x**np.longdouble(-.28))
        intercept = smooth*(1-x**exponent)
    zv, ze = nu*intercept/root_k, 2*c_length*nu/root_k
    def nut(z):
        return nu*max(0, alpha*(z*root_k/nu-intercept))
    def d(z):
        return nu+nut(z)
    def epsilon_a(z):
        return 2*nu*k/ze**2 if z < ze else k**1.5/(c_length*z)
    breaks = [np.longdouble(0), np.longdouble(height)]
    for z in (roughness, zv, ze):
        if 0 < z < height:
            breaks.append(z)
    z = min(ze, nu/(alpha*root_k))
    while z < height:
        breaks.append(z)
        z *= 2
    breaks = sorted(set(breaks))
    def integrate(function, upper=height, extra=()):
        # SciPy integrates a local [0,1] coordinate on each physical segment.
        # Extended-precision interpolation avoids loss of the tiny z-h in the
        # thin forcing-layer case. No analytic primitive is used here.
        pts = [np.longdouble(0), np.longdouble(upper)]
        pts += [p for p in breaks if 0 < p < upper]
        pts += [np.longdouble(z) for z in extra if 0 < z < upper]
        pts = sorted(set(pts))
        samples = pts + [(a+b)/2 for a, b in zip(pts[:-1], pts[1:])]
        scale = max(abs(function(s)) for s in samples)
        if not scale:
            return 0.
        result = sum((b-a)*quad(lambda s: function(a+(b-a)*np.longdouble(s))/scale,
                          0, 1, epsabs=1e-14, epsrel=3e-13, limit=200)[0]
                     for a, b in zip(pts[:-1], pts[1:]))
        return scale*result
    j0 = integrate(lambda z: 1/d(z))
    j1 = integrate(lambda z: max(z-roughness, 0)/d(z))
    j2 = integrate(lambda z: 1/d(z)**2)
    weight = integrate(lambda z: nut(z)/d(z)**2)
    # b passes through zero at h*=70. Scale that comparison by the smooth
    # intercept (and its dimensional height), rather than a vanishing b.
    near(case["intercept"], intercept, "wall_intercept", floor=smooth)
    near(case["zv"], zv, "wall_viscosity_intercept", floor=nu*smooth/root_k)
    near(case["ze"], ze, "wall_dissipation_cutoff")
    for actual, expected, name in zip(case["integrals"], (j0, j1, j2, weight),
                                       ("j0", "j1", "j2", "positive_weight")):
        near(actual, expected, "wall_"+name, rtol=2e-10)
    u, forcing = np.array(case["velocity"]), np.array(case["forcing"])
    traction = u/j0-j1/j0*forcing
    def grad(z):
        return (traction+forcing*max(z-roughness, 0))/d(z)
    production = integrate(lambda z: nut(z)*np.dot(grad(z), grad(z)))/height
    epsilon_mean = integrate(epsilon_a)/height
    near(case["original"]["traction"], traction, "wall_traction")
    near(case["original"]["production"], production, "wall_production", rtol=2e-10)
    near(case["original"]["epsilon_mean"], epsilon_mean, "wall_original_mean")
    near(case["original"]["epsilon_P"], epsilon_a(point_p), "wall_original_P")
    near(case["original"]["epsilon_H"], epsilon_a(height), "wall_original_H")
    for sample in case["samples"]:
        z = sample["z"]
        near(sample["nu_t"], nut(z), "wall_nut", rtol=2e-10, floor=nu)
        near(sample["D"], d(z), "wall_total_viscosity")
        near(sample["gradient"], grad(z), "wall_gradient", rtol=2e-10)
        velocity = traction*integrate(lambda t: 1/d(t), upper=z)+forcing*integrate(
            lambda t: max(t-roughness, 0)/d(t), upper=z) if z else [0., 0.]
        near(sample["U"], velocity, "wall_reconstruction", floor=1e-13)
        near(sample["epsilon_A"], epsilon_a(z), "wall_original_profile")
    star, t0 = math.sqrt(k/qstar), k/qstar
    numerator = star**3/kappa
    p0 = t0*t0/height*weight
    lam = height*p0/numerator
    cutoff = height*math.exp(1-lam)
    ref = case["reference"]
    for key, expected in (("ustar", star), ("traction", t0), ("U_H", t0*j0),
                          ("A", numerator), ("P0", p0), ("lambda", lam),
                          ("minimum_lambda", 1+math.log(height/point_p)), ("cutoff", cutoff),
                          ("log_z0", math.log(height)-kappa*t0*j0/star)):
        near(ref[key], expected, "equilibrium_"+key, rtol=2e-10)
    require(ref["admissible"] == (cutoff <= point_p), "Cutoff admissibility disagrees")
    for blend in case["blends"]:
        chi = blend["requested"] if ref["admissible"] else 0
        require(blend["weight"] == chi, "Incorrect effective equilibrium weight")
        out = blend["values"]
        require(out["traction"] == case["original"]["traction"] and
                out["production"] == case["original"]["production"],
                "EP changed the original momentum/production outputs")
        near(out["epsilon_mean"], (1-chi)*epsilon_mean+chi*p0, "equilibrium_mean")
        near(out["epsilon_P"], (1-chi)*epsilon_a(point_p)+chi*numerator/point_p, "equilibrium_P")
        near(out["epsilon_H"], (1-chi)*epsilon_a(height)+chi*numerator/height, "equilibrium_H")
        mean = integrate(lambda z: (1-chi)*epsilon_a(z)+chi*numerator/max(z, cutoff), extra=(cutoff,))/height
        near(out["epsilon_mean"], mean, "equilibrium_profile_integral")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--mpi-ranks", type=int, default=4)
    args = parser.parse_args()
    require(sys.platform.startswith("linux"), "Authoritative validation must run on Linux")
    allocation = sorted(os.sched_getaffinity(0))[-16:]
    require(1 <= args.mpi_ranks <= min(16, len(allocation)), "Invalid MPI allocation")
    os.sched_setaffinity(0, allocation)
    environment = dict(os.environ)
    threads = THREAD_ENVIRONMENT
    environment.update(threads)
    args.output.mkdir(parents=True, exist_ok=True)
    output = args.output.resolve()
    build = args.build.resolve()
    root = Path(__file__).resolve().parents[4]
    commands = []
    def run(command, label):
        commands.append(command)
        print(f"Running {label}", flush=True)
        result = subprocess.run(command, env=environment, cwd=output, capture_output=True,
                                text=True, timeout=120)
        (output/(label+".log")).write_text(result.stdout+result.stderr)
        require(result.returncode == 0, f"{label} failed; see {label}.log")
        print(result.stdout.strip(), flush=True)
    run([str(build/"RANSKernelProbe"), "--output", str(output/"kernels.json")], "kernels")
    kernels = json.loads((output/"kernels.json").read_text())
    roots = [reference(case) for case in kernels["references"]]
    for case in kernels["bulk"]:
        check_bulk(case)
    for case in kernels["walls"]:
        check_wall(case, *roots[1])
    require(kernels["switch_cases"][0]["weight"] == 1 and
            all(case["weight"] == 0 for case in kernels["switch_cases"][1:]),
            "Switch sensors failed independent departures")
    print("Independent dimensional-index and SciPy quadrature checks passed.", flush=True)
    run([str(build/"RANSWallLayerProbe"), "--output", str(output/"column-serial.json")], "column-serial")
    run(["mpirun", "--bind-to", "none", "--map-by", "slot", "-np", str(args.mpi_ranks),
         str(build/"RANSWallLayerProbe"), "--mpi", "--output", str(output/"column-mpi.json")], "column-mpi")
    serial = json.loads((output/"column-serial.json").read_text())
    mpi = json.loads((output/"column-mpi.json").read_text())
    fine = validate_column(serial)
    validate_column(mpi)
    maximum_difference = 0.
    for a, b in zip(serial["records"], mpi["records"]):
        for key in a:
            if key != "sample_owner":
                maximum_difference = max(maximum_difference, abs(a[key]-b[key]))
                near(a[key], b[key], "column_serial_mpi", floor=1., rtol=3e-10)
    source_files = ["library/SolverUtils/Turbulence/"+name for name in
                    ("RANSModels.h", "RANSModels.cpp", "AnalyticalWallFunction.h", "AnalyticalWallFunction.cpp", "RANSWallStrip.h")]
    source_files += ["library/UnitTests/SolverUtils/"+name for name in
                     ("RANSKernelProbe.cpp", "RANSWallLayerProbe.cpp")]
    source_files += ["library/SolverUtils/CMakeLists.txt", "solvers/IncNavierStokesSolver/CMakeLists.txt"]
    source_files += ["solvers/IncNavierStokesSolver/Tests/RANS/"+name for name in
                     ("CMakeLists.txt", "validate_m1.py", "validate_m2.py")]
    def digest(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()
    summary = {"schema": 1, "status": "passed", "milestone": "neutral_abl_rans_awf_m2",
               "utc": datetime.now(timezone.utc).isoformat(),
               "source_commit": subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip(),
               "source_status": subprocess.check_output(["git", "-C", str(root), "status", "--short"], text=True).splitlines(),
               "source_sha256": {p: digest(root/p) for p in source_files},
               "binary_sha256": {name: digest(build/name) for name in ("RANSKernelProbe", "RANSWallLayerProbe")},
               "cpu_affinity": allocation, "maximum_project_cpus": 16,
               "maximum_test_ranks": args.mpi_ranks, "thread_environment": threads,
               "python": platform.python_version(), "numpy": np.__version__,
               "scipy": __import__("scipy").__version__, "commands": commands,
               "native_summary": kernels["summary"], "oracle_maximum_scaled_errors": errors,
               "references": kernels["references"], "fine_column": fine,
               "column_serial_mpi_maximum_difference": maximum_difference,
               "reports": {name: digest(output/name) for name in
                           ("kernels.json", "column-serial.json", "column-mpi.json")}}
    (output/"validation.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps({"status": "passed", "native": kernels["summary"],
                      "oracle_errors": errors, "column_serial_mpi_maximum_difference": maximum_difference}, indent=2))


if __name__ == "__main__":
    main()
