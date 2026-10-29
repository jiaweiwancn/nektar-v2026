#!/usr/bin/env python3
"""Linear Fourier initial-value/ensemble IMEX2 calibration; NOT a solver result.

Reconstruct the supplied upstream lattice and compare SBDF2/IMEX1-started
advection-diffusion transfers with exact Fourier evolution. No spatial error,
SVV, nonlinearity, pressure, boundary feedback or finite-record centering enters.
"""
import argparse
import bisect
import cmath
import hashlib
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET

TIMESTEPS = (1e-4, 5e-5, 2.5e-5, 1.25e-5)
AGES = (0.28448, 0.65532)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_table(path):
    lines = [line.strip() for line in path.read_text(encoding="ascii").splitlines() if line.strip()]
    require(lines and lines[0] == "# NEKTAR_HIT_SPECTRUM_V1 KIND=Shell3D K_UNITS=rad/m E_UNITS=m3/s2",
            "expected native SI Shell3D spectrum header")
    rows = [tuple(map(float, line.split())) for line in lines[1:] if not line.startswith("#")]
    require(len(rows) >= 2 and all(len(row) == 2 and all(math.isfinite(x) and x > 0 for x in row)
                                 for row in rows), "table requires positive finite k,E pairs")
    ks, es = zip(*rows)
    require(all(a < b for a, b in zip(ks, ks[1:])), "table k must increase strictly")
    logk, loge = list(map(math.log, ks)), list(map(math.log, es))

    def energy(k):
        right = bisect.bisect_left(ks, k)
        if right < len(ks) and ks[right] == k:
            return es[right]
        if right == 0:
            return es[0] * (k / ks[0])**4
        fraction = (math.log(k) - logk[right-1]) / (logk[right] - logk[right-1])
        return math.exp(loge[right-1] + fraction * (loge[right] - loge[right-1]))
    return ks[-1], energy


def calibrate(args):
    source = args.manifest.resolve()
    manifest = json.loads(source.read_text(encoding="utf-8"))
    require(manifest.get("schema") == "NEKTAR_HIT_PILOT_V1", "unsupported pilot manifest schema")
    physics, synthesis = manifest["physics"], manifest["synthesis"]
    require(physics["normalization"] == "Preserve" and physics["interpolation"] == "LogLinear"
            and physics["low_extension"] == "K4", "requires Preserve/LogLinear/K4")
    periods, counts = synthesis["periods_m"], synthesis["counts"]
    cutoff, taper = synthesis["cutoff_rad_m"], synthesis["taper_start_rad_m"]
    speed, nu = physics["mean_velocity_m_s"], physics["kinvis_m2_s"]
    require(len(periods) == len(counts) == 3 and all(math.isfinite(x) and x > 0 for x in periods)
            and all(isinstance(n, int) and n >= 4 and n % 2 == 0 for n in counts), "invalid lattice")
    require(all(math.isfinite(x) for x in (cutoff, taper, speed, nu)) and 0 <= taper < cutoff
            and speed > 0 and nu > 0, "invalid cutoff/taper/mean speed/viscosity")
    case = source.parent / manifest["case"]
    session = ET.parse(case)
    parameters = {p.text.split("=", 1)[0].strip().lower(): p.text.split("=", 1)[1].strip()
                  for p in session.findall("./CONDITIONS/PARAMETERS/P")}
    require(float(parameters["ubar"]) == speed and float(parameters["kinvis"]) == nu,
            "manifest/session mean speed or viscosity mismatch")
    region = session.find("./HITINFLOW/REGION")
    spectrum = region.find("SPECTRUM")
    require(all(spectrum.get(k) == v for k, v in {"TYPE": "Tabulated", "NORMALIZATION": "Preserve",
            "INTERPOLATION": "LogLinear", "LOW": "K4"}.items()), "session spectrum differs from calibration contract")
    session_lattice = region.find("SYNTHESIS")
    require([float(session_lattice.get("L" + d)) for d in "XYZ"] == periods
            and [int(session_lattice.get("N" + d)) for d in "XYZ"] == counts
            and float(region.find("CUTOFF").get("KMAX")) == cutoff
            and float(region.find("CUTOFF").get("TAPERSTART")) == taper, "manifest/session lattice mismatch")
    table = (source.parent / spectrum.get("FILE")).resolve()
    last_k, energy = read_table(table)
    require(cutoff <= last_k, "cutoff above final table row: Gaussian high tail is deliberately unsupported")
    spacing = [2 * math.pi / length for length in periods]
    require(all(cutoff < dk * n / 2 for dk, n in zip(spacing, counts)), "cutoff violates strict Nyquist")
    limits = [min(n // 2 - 1, math.floor(cutoff / dk)) for n, dk in zip(counts, spacing)]
    rows = [{"dt_s": dt, "requested_age_s": age, "steps": round(age / dt),
             "actual_age_s": round(age / dt) * dt, "numerical": [0.0]*3,
             "exact": [0.0]*3, "cross": [0.0]*3} for dt in TIMESTEPS for age in AGES]
    covariance, mode_count = [0.0]*3, 0
    for ix in range(limits[0] + 1):
        for iy in range(0 if ix == 0 else -limits[1], limits[1] + 1):
            for iz in range(1 if ix == iy == 0 else -limits[2], limits[2] + 1):
                wave = (ix*spacing[0], iy*spacing[1], iz*spacing[2])
                k = math.hypot(*wave)
                if k >= cutoff:
                    continue
                e = energy(k)
                if k > taper:
                    e *= math.cos(0.5 * math.pi * (k-taper) / (cutoff-taper))**2
                weight = math.prod(spacing) * e / (4 * math.pi * k*k)
                c = [2 * weight * (1 - (q/k)**2) for q in wave]
                covariance = [a+b for a, b in zip(covariance, c)]
                mode_count += 1
                lam = complex(-nu*k*k, -speed*wave[0])
                for j, dt in enumerate(TIMESTEPS):
                    a, b = nu*k*k*dt, speed*wave[0]*dt
                    aa, bb, cc = 3+2*a, -4+4j*b, 1-2j*b
                    discriminant = cmath.sqrt(bb*bb - 4*aa*cc)
                    r1, r2 = (-bb+discriminant)/(2*aa), (-bb-discriminant)/(2*aa)
                    if abs(r2-cmath.exp(lam*dt)) < abs(r1-cmath.exp(lam*dt)):
                        r1, r2 = r2, r1
                    require(r1 != r2, "repeated characteristic root is outside this calibration")
                    startup = (1-1j*b)/(1+a)
                    amplitude = (startup-r2)/(r1-r2)
                    for row in rows[2*j:2*j+2]:
                        n = row["steps"]
                        numeric = amplitude*r1**n + (1-amplitude)*r2**n
                        exact = cmath.exp(lam*row["actual_age_s"])
                        factors = (abs(numeric)**2, abs(exact)**2, (numeric*exact.conjugate()).real)
                        for key, factor in zip(("numerical", "exact", "cross"), factors):
                            row[key] = [old + component*factor for old, component in zip(row[key], c)]
    kinetic = 0.5 * sum(covariance)
    require(kinetic > 0, "empty or zero-energy calibration lattice")
    checks = {"benchmark_supplied": bool(args.benchmark), "relative_tolerance": 5e-10}
    if args.benchmark:
        benchmark = json.loads(args.benchmark.read_text(encoding="utf-8"))
        expected = benchmark["expected_kinetic_energy"]
        error = abs(kinetic-expected)/abs(expected) if expected else math.inf
        require(math.isfinite(error) and error <= 5e-10, "native expected kinetic energy crosscheck failed")
        checks["kinetic_energy_relative_error"] = error
        if "expected_covariance_diagonal" in benchmark:
            target = benchmark["expected_covariance_diagonal"]
            require(len(target) == 3, "native expected covariance must have three components")
            errors = [abs(a-b)/abs(b) if b else (0.0 if a == 0 else math.inf) for a, b in zip(covariance, target)]
            require(all(math.isfinite(x) and x <= 5e-10 for x in errors), "native covariance crosscheck failed")
            checks["covariance_diagonal_relative_errors"] = errors
        require("modes" not in benchmark or benchmark["modes"] == mode_count, "native mode-count mismatch")
    for row in rows:
        row["relative_component_variance_error"] = [a/b-1 if b else (0.0 if a == 0 else None)
                                                   for a, b in zip(row["numerical"], row["exact"])]
        row["normalized_component_correlation"] = [c/math.sqrt(a*b) if a*b else None for a, b, c in
            zip(row["numerical"], row["exact"], row["cross"])]
        for key in ("numerical", "exact", "cross"):
            row[key + "_covariance_diagonal_m2_s2"] = row.pop(key)
    paths = [source, case.resolve(), table, Path(__file__).resolve()] + ([args.benchmark.resolve()] if args.benchmark else [])
    require(args.output.resolve() not in paths, "output must not overwrite an input or script")
    result = {"schema": "NEKTAR_HIT_LINEAR_IMEX2_CALIBRATION_V1", "scope": __doc__,
              "method": "SBDF2 physical root closest to exp(lambda*dt), exact two-root solution with IMEX1 startup",
              "lambda": "-nu*k^2-i*U*kx", "characteristic": "(3+2*nu*k^2*dt)*r^2+(-4+4i*U*kx*dt)*r+(1-2i*U*kx*dt)=0",
              "normalization": "Ensemble modal covariance 2*DeltaKx*DeltaKy*DeltaKz*E(k)/(4*pi*k^2)*(I-khat*khat); no temporal or plane centering",
              "correlation": "Sum Cii*Re(numeric*conj(exact)) / sqrt(sum Cii*abs(numeric)^2 * sum Cii*abs(exact)^2)",
              "spectrum_interpolation": "Native LogLinear: linear log(E) versus log(k), K4 below first row",
              "independent_modes": mode_count, "initial_expected_kinetic_energy_m2_s2": kinetic,
              "initial_expected_covariance_diagonal_m2_s2": covariance,
              "native_crosscheck": checks, "rows": rows,
              "sha256": {str(path): sha(path) for path in paths}}
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--benchmark", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        calibrate(args)
    except (ValueError, OSError, KeyError, TypeError, AttributeError, OverflowError, ZeroDivisionError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
