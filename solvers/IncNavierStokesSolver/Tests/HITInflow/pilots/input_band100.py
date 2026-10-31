#!/usr/bin/env python3
"""Continuum input budget for the CBC100 pilot; execute on Linux, one CPU.

These quadratures are model predictions, not a finite seeded realization,
boundary projection, solver result, experimental uncertainty, or fidelity pass.
"""
import argparse
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import platform

from scipy.integrate import quad
from spectrum_budget import Spectrum, make_case, taper

HERE = Path(__file__).resolve().parent
POLICY = {"interpolation": "LogLinear", "low_extension": "K4", "high_extension": "Gaussian"}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def projection(spectrum, q, cutoff=None, start=None):
    endpoint = cutoff if cutoff is not None else spectrum.k[-1]
    if cutoff is not None and q >= cutoff:
        return {"value_m3_s2": 0.0, "absolute_quadrature_error_estimate": 0.0}
    points = sorted({q, endpoint, *(k for k in spectrum.k if q < k < endpoint),
                     *([start] if start is not None and q < start < endpoint else [])})

    def integrand(k):
        return spectrum.evaluate(k, POLICY) * taper(k, cutoff, start) * (1-(q/k)**2)/k

    parts = [quad(integrand, a, b, epsabs=1e-13, epsrel=2e-10, limit=200)
             for a, b in zip(points[:-1], points[1:])]
    if cutoff is None:
        parts.append(quad(lambda t: integrand(endpoint+t/spectrum.sqrt_beta)/spectrum.sqrt_beta,
                          0, math.inf, epsabs=1e-13, epsrel=2e-10, limit=200))
    value, error = (math.fsum(p[i] for p in parts) for i in (0, 1))
    if not all(math.isfinite(v) and v >= 0 for v in (value, error)):
        raise ValueError("nonfinite/negative continuum projection")
    return {"value_m3_s2": value, "absolute_quadrature_error_estimate": error}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    source = HERE.parent / "cbc1971-source.json"
    table = HERE.parent / "cbc1971_theta42_spectrum_3d.dat"
    checkout = HERE.parents[4]
    if output.exists() or output == checkout or checkout in output.parents:
        parser.error("output must be fresh and outside the source checkout")
    spectrum = Spectrum(table.read_bytes())
    experimental = json.loads(source.read_text())
    rows = [(float(r[0])*100, float(r[1])*1e-6)
            for r in experimental["table2a"]["rows"] if r[1] is not None and float(r[0])*100 <= 100]
    uncut = {q: projection(spectrum, q) for q, _ in rows}
    cases = []
    for cutoff, start in ((100, 80), (600, 480), (800, 640)):
        budget = make_case(spectrum, POLICY, cutoff, start, 1.5e-5, quad)
        projected = []
        for q, observed in rows:
            own, full = projection(spectrum, q, cutoff, start), uncut[q]
            projected.append({"k1_rad_m": q, "raw_Table2a_F11_m3_s2": observed,
                "uncut_model": full, "tapered_model": own,
                "retained_percent_of_uncut_model": 100*own["value_m3_s2"]/full["value_m3_s2"],
                "signed_percent_difference_from_raw_measurement": 100*(own["value_m3_s2"]/observed-1)})
        modes = 10*.6*.6*cutoff**3/(12*math.pi**2)
        cases.append({"cutoff_rad_m": cutoff, "taper_start_rad_m": start,
            "continuum_budget": budget, "longitudinal_projection": projected,
            "sphere_volume_independent_mode_estimate": modes,
            "mode_payload_estimate_GiB_per_rank": modes*88/2**30,
            "memory_note": "Density estimate for periods(10,.6,.6)m, 88-byte current native mode; not an allocation measurement or native count"})
    result = {"schema": "NEKTAR_HIT_BAND100_INPUT_V1", "status": "ANALYTICAL_MODEL_ONLY",
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "host": platform.node(), "cpu_affinity": sorted(os.sched_getaffinity(0)),
        "source_sha256": {str(path): digest(path) for path in (source, table, Path(__file__), HERE/'spectrum_budget.py')},
        "uncut_continuum_budget": make_case(spectrum, POLICY, None, None, 1.5e-5, quad),
        "cases": cases,
        "limits": ["No downstream data are fitted and no retained-band amplitude is restored.",
            "Shell-table interpolation/extensions and full-band experimental totals are distinct.",
            "Raw Table2a is descriptive; a positive isotropic E implies nonincreasing F11, whereas the lowest published bins are not monotone.",
            "Primary pilot spectral bands start at10rad/m; no uniform experimental uncertainty is assumed.",
            "Analytical truncation error is not solver error or a physical-validation verdict."]}
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x') as stream:
        stream.write(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(output)


if __name__ == "__main__":
    main()
