#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Certify the represented finite increment of the actual Bf initial state."""

import argparse
from pathlib import Path
import time

import numpy as np

import perturbation_certificate as pc
import perturbation_check as checks
import physical_oracle as po
from resources import atomic, identity
from task2r_advancing import read, require
from task2s_local_audit import admissibility


def dependent_trace(base, new, field_noise):
    ref = new.ref
    total_epsilon = new.s["epsilon"]-base.s["epsilon"]
    total_epsilon_norm = po.Norms()
    max_ratio_error = max_actual_ratio = max_target_ratio = cutoff_change = cutoff_noise = 0.
    peak = None
    for e in new.elements:
        if not e["epsilon"]:
            continue
        new.clock()
        xl, xr, rl, rr = e["bounds"]
        a, w = new.gauss
        x, r = np.meshgrid((xl+xr)/2+(xr-xl)*a/2, (rl+rr)/2+(rr-rl)*a/2)
        x, r = x.ravel(), r.ravel()
        weights = np.outer(w,w).ravel()*(xr-xl)*(rr-rl)/4
        delta_e = po.basis(e, "epsilon", x, r)[0] @ total_epsilon[e["epsilon"]]
        total_epsilon_norm.add("outer", "total_epsilon_increment", delta_e*r/ref.A, weights)
        if rl != ref.H:
            continue
        # Resolve both H owners at x=0 and every exterior panel endpoint.
        x = np.r_[xl, (xl+xr)/2+(xr-xl)*a/2, xr]
        eps = po.basis(e, "epsilon", x, np.full(len(x), ref.H))[0]
        e0, e1 = (eps@oracle.s["epsilon"][e["epsilon"]] for oracle in (base, new))
        k0, k1 = (oracle.trace(e, "k", x) for oracle in (base, new))
        r0, r1 = ref.cutoff(k0), ref.cutoff(k1)
        ratio = (k1/k0)**1.5
        max_ratio_error = max(max_ratio_error, float(np.max(abs(e1/e0-ratio))))
        max_actual_ratio = max(max_actual_ratio, float(np.max(abs(e1/e0-1))))
        max_target_ratio = max(max_target_ratio, float(np.max(abs(ratio-1))))
        cutoff_change = max(cutoff_change, float(np.max(abs(r1-r0))))
        cutoff_noise = max(cutoff_noise, float(np.max(abs(
            ref.cutoff(k0+ref.Kin*field_noise)-ref.cutoff(k0-ref.Kin*field_noise))))+1e-12*ref.H)
        i = int(np.argmax((k1-k0)/ref.Kin))
        row = dict(owner=e["geometry"], x=float(x[i]), r=ref.H, K_base=float(k0[i]),
                   K_perturbed=float(k1[i]), delta_K_over_Kin=float((k1[i]-k0[i])/ref.Kin),
                   exact_Ew_ratio=float(ratio[i]), represented_trace_ratio=float(e1[i]/e0[i]),
                   cutoff_base=float(r0[i]), cutoff_perturbed=float(r1[i]))
        if peak is None or row["delta_K_over_Kin"] > peak["delta_K_over_Kin"]:
            peak = row
    return dict(sampled_K_peak=peak, max_target_Ew_relative_change=max_target_ratio,
                max_represented_trace_relative_change=max_actual_ratio,
                max_pointwise_ratio_projection_error=max_ratio_error,
                largest_cutoff_change=cutoff_change, cutoff_noise_bound=cutoff_noise,
                total_epsilon_increment=total_epsilon_norm.report({}))


def audit(base_path, native_path, table, proof_path, deadline):
    base_data, base_run = po.load_native(base_path)
    data, run = po.load_native(native_path)
    proof = read(proof_path)
    ref = po.load_reference(table)
    require(base_run["steady_mode"] and base_run["steady_solver"]["converged"]
            and run["steady_mode"] is False and run["restart"] is False
            and data["step"] == data["time"] == 0 and not data.get("input_history")
            and proof["status"] == "PASS_MATHEMATICAL_SHAPE_ONLY"
            and proof["shape"] == pc.BF_SHAPE and proof["base_state_sha256"] == base_data["state_sha256"]
            and proof["source_identity"] == base_data["source_identity"] == data["source_identity"]
            and proof["reference_identity"] == ref.identity and data["elements"] == base_data["elements"]
            and data["state_sha256"] == data["previous_sha256"]
            and data["advection_sha256"] == data["previous_advection_sha256"],
            "perturbed state, shape, base or fresh history lineage differs")
    base, new = (po.PhysicalOracle(s, ref, 16, deadline) for s in (base_data, data))
    delta = {k: new.s[k]-base.s[k] for k in new.s}
    lift = np.asarray(data["essential_epsilon_lift"])-np.asarray(base_data["essential_epsilon_lift"])
    delta_epsilon_independent = delta["epsilon"]-lift
    normalization = proof["normalization"]
    layout = normalization["layout"]
    packet = dict(elements=data["elements"], sizes=data["sizes"], domain=proof["domain"], H=ref.H,
        signs=[1,1,1], Kin=ref.Kin, Ain=ref.A, velocity=delta["flow"], k_increment=delta["k"],
        epsilon_independent_increment=delta_epsilon_independent,
        normalization=dict(streamfunction_coefficient=.001/normalization["unitPeakEstimate"]))
    norms, diagnostics, norms_ok = checks.represented_norms(packet, proof, deadline, synthetic_base=False)
    divergence = checks.divergence_audit(packet, deadline)
    boundary = dict(base=base.boundary_checks(), perturbed=new.boundary_checks())
    trace = dependent_trace(base, new, diagnostics["global_representation_error_upper"]["k"])
    exterior = dict(U=float(np.max(abs(delta["flow"][new.fixed["u"]]))),
        W=float(np.max(abs(delta["flow"][new.n+new.fixed["w"]]))),
        k=float(np.max(abs(delta["k"][new.fixed["k"]]))/ref.Kin),
        independent_epsilon=float(np.max(abs(delta_epsilon_independent[new.fixed["epsilon"]]))/(ref.A/ref.H)))
    pressure = float(np.max(abs(delta["flow"][2*new.n:])))
    margins = admissibility(new)
    principal = margins["regions"]["whole"]
    noise = diagnostics["global_representation_error_upper"]
    gradient_noise = diagnostics["gradient_representation_error_upper"]["velocity"]
    xlo, xhi = layout["x_knots"][::4]
    rlo, rhi = layout["r_knots"][::4]
    widths = [sum(a >= lo-1e-12 and b <= hi+1e-12 for a,b in zip(edges[:-1],edges[1:]))
              for edges, lo, hi in ((proof["x_edges"],xlo,xhi), (proof["r_edges"],rlo,rhi),
                                    (proof["r_edges"],ref.H,rhi))]
    support = (xhi-xlo)*(rhi-rlo)
    gates = dict(norms=norms_ok, support=bool(support >= 1.5 and min(widths) >= 4 and rlo < ref.H < rhi),
        kinetic_energy=norms["velocity"]["domain_RMS_enclosure"][0]**2/2 >= 5e-11,
        solenoidal=divergence["pass_"], exterior=max(exterior.values()) <= 1e-12 and pressure == 0.,
        dependent_trace=all(row["status"] == "MET" for row in boundary.values())
            and max(trace["max_target_Ew_relative_change"],trace["max_represented_trace_relative_change"]) <= .02
            and abs(trace["sampled_K_peak"]["delta_K_over_Kin"]-.01) <= 1e-12,
        resolved_signal=all(norms[k]["global_peak_enclosure"][0] > 10*noise[k] for k in norms)
            and norms["velocity"]["normalized_field_gradient_RMS_enclosure"][0] > 10*gradient_noise
            and trace["largest_cutoff_change"] > 10*trace["cutoff_noise_bound"],
        admissibility=principal["minimum_covariance_over_k"]["value"] > 1e-10
            and principal["minimum_planar_over_mu"]["value"] > max(1e-8, 10*principal["half_step_change"]))
    return dict(status="MET" if all(gates.values()) else "NOT_MET", gates=gates, layout=layout,
        support_area=support, resolved_intervals=widths, represented_norms=norms,
        polynomial_diagnostics=diagnostics, weak_divergence=divergence, boundary=boundary,
        pressure_increment=pressure, exterior_increment=exterior, current_K_trace_and_cutoff=trace,
        admissibility=margins, state_sha256=data["state_sha256"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("base", "native", "table", "certificate", "out"):
        parser.add_argument("--"+name, type=Path, required=True)
    parser.add_argument("--seconds", type=float, required=True)
    args = parser.parse_args()
    began = time.monotonic()
    result = dict(schema="WJH27R3PhysicalBfPerturbationAuditV1", status="NOT_MET",
                  base=identity(args.base), native=identity(args.native), table=identity(args.table),
                  certificate=identity(args.certificate), analyzer=identity(__file__))
    try:
        result.update(audit(args.base,args.native,args.table,args.certificate,began+args.seconds))
    except Exception as error:
        result.update(reason=str(error), exception=type(error).__name__)
    result["elapsed_seconds"] = time.monotonic()-began
    atomic(args.out,result)
    if result["status"] != "MET":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
