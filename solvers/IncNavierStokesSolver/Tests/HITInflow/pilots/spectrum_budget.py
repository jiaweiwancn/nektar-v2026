#!/usr/bin/env python3
"""Continuum SI Shell3D spectrum-policy budgets; requires SciPy explicitly.

These integrals describe the input model, without finite-box sampling, seeded
realization, CFD, downstream measurements, or any energy renormalization.
"""

import argparse
from bisect import bisect_left
import hashlib
import json
import math
from pathlib import Path
import re
import sys


HEADER = "# NEKTAR_HIT_SPECTRUM_V1 KIND=Shell3D K_UNITS=rad/m E_UNITS=m3/s2"
DECIMAL = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z")
MIN_LOG = math.log(float.fromhex("0x0.0000000000001p-1022"))
EPSABS, EPSREL, QUAD_LIMIT = 1e-13, 2e-10, 200
MOMENTS = (("M0", 0, "m2/s2"), ("M2", 2, "1/s2"), ("Mminus1", -1, "m3/s2"))
BUDGET_BREAKS = (80.0, 100.0, 160.0, 200.0)


def log_ratio(a, b):
    """Native LogRatio for a>b>0; avoid forming a/b."""
    difference = a - b
    return math.log1p(difference / b) if difference <= b else math.log(a) - math.log(b)


def log_sum(a, b):
    return math.log(max(a, b)) + math.log1p(min(a, b) / max(a, b))


def bounded_exp(logarithm, upper):
    return 0.0 if logarithm < MIN_LOG else min(math.exp(logarithm), upper)


class Spectrum:
    def __init__(self, raw):
        self.k, self.e = [], []
        header = False
        for line_number, line in enumerate(raw.decode("utf-8").split("\n"), 1):
            line = line.strip(" \t\r\n\f\v")
            if not line:
                continue
            if not header:
                if line != HEADER:
                    raise ValueError("first nonblank line must declare version-1 SI Shell3D exactly")
                header = True
                continue
            if line.startswith("#"):
                continue
            fields = re.split(r"[ \t\r\n\f\v]+", line)
            if len(fields) != 2 or any(not DECIMAL.fullmatch(token) for token in fields):
                raise ValueError(f"line {line_number}: expected exactly two decimal columns")
            k, energy = map(float, fields)
            if not math.isfinite(k) or k <= 0 or (self.k and k <= self.k[-1]):
                raise ValueError(f"line {line_number}: k must be positive, finite, strictly increasing")
            if not math.isfinite(energy) or energy <= 0:
                raise ValueError(f"line {line_number}: selected LogLinear cases require positive finite E")
            self.k.append(k)
            self.e.append(energy)
        if not header or len(self.k) < 2:
            raise ValueError("spectrum requires the native header and at least two rows")
        if not self.e[-2] > self.e[-1]:
            raise ValueError("Gaussian extension requires a decreasing positive final pair")
        # Same stable factorization as native m_gaussianLogBeta.
        self.log_beta = (math.log(log_ratio(self.e[-2], self.e[-1]))
                         - math.log(self.k[-1] - self.k[-2]) - log_sum(self.k[-1], self.k[-2]))
        self.sqrt_beta = math.exp(.5 * self.log_beta)
        if not math.isfinite(self.sqrt_beta) or self.sqrt_beta <= 0:
            raise ValueError("Gaussian tail scale is outside this float64 integration tool's range")

    def evaluate(self, k, policy):
        if not math.isfinite(k) or k < 0:
            raise ValueError("spectrum evaluation k must be finite and nonnegative")
        if k == 0:
            return 0.0
        right = bisect_left(self.k, k)
        if right < len(self.k) and self.k[right] == k:
            return self.e[right]
        if right == 0:
            return (0.0 if policy["low_extension"] == "Zero" else
                    bounded_exp(math.log(self.e[0]) - 4 * log_ratio(self.k[0], k), self.e[0]))
        if right == len(self.k):
            if policy["high_extension"] == "Zero":
                return 0.0
            log_decay = self.log_beta + math.log(k - self.k[-1]) + log_sum(k, self.k[-1])
            if log_decay > 10:  # Native guard: decay is already far below subnormal range.
                return 0.0
            return bounded_exp(math.log(self.e[-1]) - math.exp(log_decay), self.e[-1])
        left = right - 1
        el, er = self.e[left], self.e[right]
        if el == er:
            return el
        if policy["interpolation"] == "Linear":
            fraction = (k - self.k[left]) / (self.k[right] - self.k[left])
            return el + fraction * (er - el) if el <= er else er + (1 - fraction) * (el - er)
        fraction = log_ratio(k, self.k[left]) / log_ratio(self.k[right], self.k[left])
        value = bounded_exp((1 - fraction) * math.log(el) + fraction * math.log(er), max(el, er))
        return max(value, min(el, er))


def taper(k, cutoff, start):
    if cutoff is None or k <= start:
        return 1.0
    if k >= cutoff:
        return 0.0
    return math.cos(.5 * math.pi * (k - start) / (cutoff - start))**2


def integral(spectrum, policy, power, cutoff, start, quad):
    endpoint = (cutoff if cutoff is not None else max(spectrum.k[-1], *BUDGET_BREAKS)
                if policy["high_extension"] == "Gaussian" else spectrum.k[-1])
    points = sorted({0.0, endpoint, *(k for k in spectrum.k if k < endpoint),
                     *(k for k in BUDGET_BREAKS if k < endpoint),
                     *([start] if start is not None and 0 < start < endpoint else [])})
    bounds = list(zip(points[:-1], points[1:]))
    if cutoff is None and policy["high_extension"] == "Gaussian":
        bounds.append((endpoint, None))
    segments = []
    for lower, upper in bounds:
        def integrand(k):
            energy = spectrum.evaluate(k, policy) * taper(k, cutoff, start)
            return energy * k**power if energy and k else 0.0

        if upper is None:
            # Dimensionless Gaussian-tail coordinate avoids an arbitrary finite
            # truncation or a tail width poorly scaled for QUADPACK's infinity map.
            def transformed(t):
                return integrand(spectrum.k[-1] + t / spectrum.sqrt_beta) / spectrum.sqrt_beta
            function, a, b = transformed, (lower - spectrum.k[-1]) * spectrum.sqrt_beta, math.inf
            variable = "t=(k-klast)*sqrt(beta)"
        else:
            function, a, b = integrand, lower, upper
            variable = "k_rad_per_m"
        output = quad(function, a, b, epsabs=EPSABS, epsrel=EPSREL, limit=QUAD_LIMIT, full_output=1)
        if len(output) != 3:
            raise ValueError(f"quadrature did not converge on [{lower},{upper}]: {output[3]}")
        value, error, info = output
        if not all(math.isfinite(v) and v >= 0 for v in (value, error)):
            raise ValueError("quadrature returned invalid nonnegative moment/error")
        region = ("below_first_table_k" if upper is not None and upper <= spectrum.k[0]
                  else "above_last_table_k" if lower >= spectrum.k[-1] else "tabulated_k_range")
        segments.append({"k_limits_rad_m": [lower, upper if upper is not None else "infinity"],
                         "integration_variable": variable, "variable_limits": [a, b if math.isfinite(b) else "infinity"],
                         "value": value, "absolute_error_estimate": error,
                         "function_evaluations": info["neval"], "subintervals": info["last"], "region": region})
    return {"value": math.fsum(s["value"] for s in segments),
            "absolute_error_estimate": math.fsum(s["absolute_error_estimate"] for s in segments),
            "contributions": {region: math.fsum(s["value"] for s in segments if s["region"] == region)
                              for region in ("below_first_table_k", "tabulated_k_range", "above_last_table_k")},
            "segments": segments}


def make_case(spectrum, policy, cutoff, start, viscosity, quad):
    moments = {name: {"power_of_k": power, "units": units,
                      **integral(spectrum, policy, power, cutoff, start, quad)}
               for name, power, units in MOMENTS}
    kinetic = moments["M0"]["value"]
    epsilon = 2 * viscosity * moments["M2"]["value"]
    if kinetic <= 0 or epsilon <= 0:
        raise ValueError("positive K and dissipation are required for the selected length scales")
    return {"policies": policy, "normalization": "Preserve (no scaling)",
            "cutoff_rad_m": cutoff, "taper_start_rad_m": start,
            "requested_k_limits_rad_m": [0, "infinity"],
            "effective_k_limits_rad_m": [0, cutoff if cutoff is not None else
                spectrum.k[-1] if policy["high_extension"] == "Zero" else "infinity"],
            "moments": moments,
            "continuum_input_model": {"K_m2_s2": kinetic, "epsilon_m2_s3": epsilon,
                "epsilon_absolute_quadrature_error_estimate_m2_s3": 2 * viscosity * moments["M2"]["absolute_error_estimate"],
                "isotropic_component_urms_m_s": math.sqrt(2 * kinetic / 3),
                "lambda_m": math.sqrt(10 * viscosity * kinetic / epsilon),
                "Lf_m": 3 * math.pi * moments["Mminus1"]["value"] / (4 * kinetic)},
            "Mminus1_finiteness": "Finite: low-k extension is K4 or Zero; high-k extension is Gaussian or Zero"}


def sum_segments(moment, predicate):
    segments = [s for s in moment["segments"] if predicate(*s["k_limits_rad_m"])]
    return {"value": math.fsum(s["value"] for s in segments),
            "absolute_error_estimate": math.fsum(s["absolute_error_estimate"] for s in segments)}


def excluded_budget(case, reference, reference_name):
    cutoff, start = case["cutoff_rad_m"], case["taper_start_rad_m"]
    result = {"uncut_same_policy_reference": reference_name, "moments": {}}
    for name, _, units in MOMENTS:
        own, full = case["moments"][name], reference["moments"][name]
        above = sum_segments(full, lambda lower, upper: lower >= cutoff)
        in_band = lambda lower, upper: lower >= start and upper != "infinity" and upper <= cutoff
        full_band, own_band = sum_segments(full, in_band), sum_segments(own, in_band)
        taper_removed = {"value": full_band["value"] - own_band["value"],
                         "absolute_error_estimate": full_band["absolute_error_estimate"] + own_band["absolute_error_estimate"]}
        removed = above["value"] + taper_removed["value"]
        result["moments"][name] = {"units": units, "above_cutoff": above,
            "removed_inside_taper_band": taper_removed, "total_removed": removed,
            "total_removed_absolute_error_estimate": above["absolute_error_estimate"] + taper_removed["absolute_error_estimate"],
            "retained_percent": 100 * own["value"] / full["value"],
            "excluded_percent": 100 * removed / full["value"],
            "closure_residual_full_minus_retained_minus_removed": full["value"] - own["value"] - removed}
    result["epsilon_retained_percent"] = result["moments"]["M2"]["retained_percent"]
    result["epsilon_excluded_percent"] = result["moments"]["M2"]["excluded_percent"]
    viscous_factor = case["continuum_input_model"]["epsilon_m2_s3"] / case["moments"]["M2"]["value"]
    m2 = result["moments"]["M2"]
    result["epsilon_removed_m2_s3"] = viscous_factor * m2["total_removed"]
    result["epsilon_removed_absolute_error_estimate_m2_s3"] = viscous_factor * m2["total_removed_absolute_error_estimate"]
    result["epsilon_above_cutoff_m2_s3"] = viscous_factor * m2["above_cutoff"]["value"]
    result["epsilon_removed_inside_taper_band_m2_s3"] = viscous_factor * m2["removed_inside_taper_band"]["value"]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, required=True, help="Native version-1 SI Shell3D table")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--kinvis", type=float, default=1.5e-5, help="Positive molecular viscosity (m2/s)")
    parser.add_argument("--native-source-root", type=Path,
                        help="Optional clean Nektar source root for policy-source hashes")
    args = parser.parse_args()
    try:
        if not math.isfinite(args.kinvis) or args.kinvis <= 0:
            raise ValueError("--kinvis must be finite and positive")
        if args.output.resolve() == args.spectrum.resolve():
            raise ValueError("output must not overwrite the spectrum input")
        import scipy
        from scipy.integrate import quad
        raw = args.spectrum.read_bytes()
        spectrum = Spectrum(raw)
        selected = {"interpolation": "LogLinear", "low_extension": "K4", "high_extension": "Gaussian"}
        policies = {"selected": selected, "lowZero": selected | {"low_extension": "Zero"},
                    "highZero": selected | {"high_extension": "Zero"},
                    "Linear": selected | {"interpolation": "Linear"}}
        cases = {"full-" + name: make_case(spectrum, policy, None, None, args.kinvis, quad)
                 for name, policy in policies.items()}
        for cutoff, start in ((100.0, 80.0), (200.0, 160.0)):
            for name in ("selected", "Linear", "lowZero"):
                case = make_case(spectrum, policies[name], cutoff, start, args.kinvis, quad)
                case["excluded_band_budget"] = excluded_budget(case, cases["full-" + name], "full-" + name)
                cases[f"cutoff{cutoff:g}-taper{start:g}-{name}"] = case
        reference = cases["full-selected"]
        for case in cases.values():
            case["relative_to_full_selected_percent"] = {
                name: 100 * (case["moments"][name]["value"] / reference["moments"][name]["value"] - 1)
                for name, _, _ in MOMENTS}
        source_root = args.native_source_root.resolve() if args.native_source_root else next(
            (p for p in Path(__file__).resolve().parents
             if (p / "library/SolverUtils/Turbulence/HITSpectrum.cpp").is_file()), None)
        native = ([source_root / "library/SolverUtils/Turbulence" / name
                   for name in ("HITSpectrum.cpp", "HITGenerator.cpp")] if source_root else [])
        if args.native_source_root and not all(p.is_file() for p in native):
            raise ValueError("--native-source-root lacks the native spectrum/generator policy sources")
        result = {"schema": "NEKTAR_HIT_CONTINUUM_SPECTRUM_BUDGET_V1",
            "inputs": {"spectrum": str(args.spectrum.resolve()), "spectrum_sha256": hashlib.sha256(raw).hexdigest(),
                "script": str(Path(__file__).resolve()), "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "native_policy_sources_sha256_if_available": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in native if p.is_file()},
                "table_rows": len(spectrum.k), "table_k_limits_rad_m": [spectrum.k[0], spectrum.k[-1]]},
            "implementation": {"python": sys.version, "scipy": scipy.__version__, "arithmetic": "Python float64; native C++ long-double intermediates/fma can round differently",
                "quadrature": "scipy.integrate.quad (QUADPACK), segmented at table knots, 80/100/160/200 budget boundaries, taper and cutoff",
                "epsabs_per_segment_in_each_moment_units": EPSABS, "epsrel_per_segment": EPSREL,
                "maximum_subintervals_per_segment": QUAD_LIMIT,
                "gaussian_infinite_tail_variable": "t=(k-klast)*sqrt(beta), t in [0,infinity); no finite tail truncation",
                "gaussian_log_beta": spectrum.log_beta,
                "error_estimates": "QUADPACK absolute error estimates summed over segments; not rigorous bounds or input-model/experimental uncertainties. Derived length-scale errors are not propagated."},
            "kinvis_m2_s": args.kinvis, "cases": cases,
            "formulas": {"M0": "integral_0^infinity E(k) dk = K", "M2": "integral_0^infinity k^2 E(k) dk",
                "Mminus1": "integral_0^infinity E(k)/k dk", "epsilon": "2*nu*M2", "urms": "sqrt(2*K/3)",
                "lambda": "sqrt(10*nu*K/epsilon)", "Lf": "3*pi*Mminus1/(4*K)",
                "LogLinear": "linear interpolation of log(E) against log(k)", "K4": "E0*(k/k0)^4 below first table k",
                "Gaussian": "Elast*exp(-beta*(k^2-klast^2)); beta=log(Eprev/Elast)/(klast^2-kprev^2)",
                "energy_taper": "1 below taper; cos(pi/2*(k-taper)/(cutoff-taper))^2 in taper band; 0 at or above cutoff"},
            "limitations": ["Continuum input-policy estimates, not seeded finite-box observables or CFD results",
                "Discrete lattice quadrature, finite periods, random realization and solver resolution are not represented",
                "Periodic correlation can retain a nonzero asymptote from kx=0 modes; do not identify its integral with this continuum Lf",
                "K4/Gaussian extensions and interpolation are model choices; their sensitivity is not experimental uncertainty",
                "No normalization or fitting; no downstream experimental data were loaded"]}
        args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    except ImportError as error:
        parser.error(f"explicit SciPy dependency is unavailable: {error}")
    except (ValueError, OSError, OverflowError, TypeError, KeyError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
