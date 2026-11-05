#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Observe limiter branches and native moving panels in flat-state readbacks.

Panel endpoints are reconstructed from the actual ordered Gauss points, not
from a frozen profile or a differentiated maximum. This observes panel counts,
locations and evaluation elements, not the origin labels of tied breakpoints.
No residual, state, qualification bound or physical history is modified.
"""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import sys

THREADS = {name: "1" for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "MKL_NUM_THREADS", "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
           "NUMEXPR_NUM_THREADS")}
os.environ.update(THREADS)
import numpy as np


def digest(value):
    return hashlib.sha256(json.dumps(value, separators=(",", ":")).encode()).hexdigest()


def inventory(path):
    data = json.loads(path.read_text())
    if data["mapped"] or data["dimension"] != 2 or data["ranks"] != 1:
        raise ValueError("this diagnostic requires the existing serial flat2D scope")
    domain = data["scalar_domain"]
    if not domain["admissible"] or not domain.get("observations"):
        raise ValueError("a complete native scalar evaluation inventory is required")
    order = int(data["quadrature_orders"][1])
    nodes, _ = np.polynomial.legendre.leggauss(order)
    sampled = {tuple(p["coordinate"][:2]): p["value"]
               for p in domain["observations"] if p["inventory"] == "physical_sampled_K"}
    groups = []
    ray_names = {"volume_normalization_ray", "trace_baseline_normalization_ray",
                 "trace_current_normalization_ray"}
    for p in domain["observations"]:
        if p["inventory"] not in ray_names:
            continue
        key = (p["inventory"], *p["coordinate"][:2])
        if not groups or groups[-1][0] != key:
            groups.append((key, []))
        groups[-1][1].append(p)
    pattern = domain["profile_branch_pattern"]
    if len(pattern) != 2*len(groups) or any(c not in "BEM" for c in pattern[:len(groups)]) or \
            any(c not in "CP" for c in pattern[len(groups):]):
        raise ValueError("native profile/cap inventory ordering differs")
    theta = data["robust_parameters"][2]
    rays = []
    for i, (key, points) in enumerate(groups):
        if len(points) % order:
            raise ValueError("native ray does not contain complete Gauss panels")
        sampled_k = sampled[key[1:]]
        if sampled_k <= 0:
            raise ValueError("undefined sampled K in native readback")
        gaps = [p["value"]/sampled_k-theta for p in points]
        limiter = "".join("L" if gap < 0 else "P" for gap in gaps)
        panels = []
        for j in range(0, len(points), order):
            panel = points[j:j+order]
            z = np.array([p["coordinate"][2] for p in panel])
            half = (z[-1]-z[0])/(nodes[-1]-nodes[0])
            mid = .5*(z[-1]+z[0])
            if half <= 0 or np.max(np.abs(z-(mid+half*nodes))) > 2e-13*max(1.,abs(mid)):
                raise ValueError("ordered native points fail the independent Gauss-panel check")
            panels.append(dict(endpoints=[float(mid-half),float(mid+half)],
                               elements=[p["geometry"] for p in panel]))
        rays.append(dict(key=key, profile=pattern[i], cap=pattern[len(groups)+i],
                         minimum_limiter_gap=min(abs(x) for x in gaps),
                         limiter=limiter, panels=panels))
    return data, rays


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("readbacks", nargs="+", type=Path)
    args = parser.parse_args()
    results = []
    base = base_rays = None
    for path in args.readbacks:
        data, rays = inventory(path)
        if base is None:
            base, base_rays = data, rays
        for key in ("history", "alpha", "time", "fields", "parameters", "robust_parameters",
                    "inner_stress_model", "model_variant", "map_parameters", "transition_controls"):
            if data[key] != base[key]:
                raise ValueError("branch comparison changes recorded target data: " + key)
        if [r["key"] for r in rays] != [r["key"] for r in base_rays]:
            raise ValueError("branch comparison changes the independent ray inventory")
        changes = []
        for a, b in zip(base_rays, rays):
            same_size = len(a["panels"]) == len(b["panels"])
            switches = sum(x != y for x,y in zip(a["limiter"],b["limiter"])) if same_size else None
            endpoint_change = max(abs(x-y) for p,q in zip(a["panels"],b["panels"])
                                  for x,y in zip(p["endpoints"],q["endpoints"])) if same_size else None
            elements_changed = not same_size or any(p["elements"] != q["elements"]
                for p,q in zip(a["panels"],b["panels"]))
            changes.append(dict(ray=b["key"], profile_changed=a["profile"]!=b["profile"],
                cap_changed=a["cap"]!=b["cap"], panel_count_change=len(b["panels"])-len(a["panels"]),
                evaluation_elements_changed=elements_changed, limiter_switches=switches,
                maximum_endpoint_displacement=endpoint_change))
        changed = [r for r in changes if r["profile_changed"] or r["cap_changed"] or
                   r["panel_count_change"] or r["evaluation_elements_changed"] or r["limiter_switches"]]
        results.append(dict(path=str(path.resolve()),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            state_sha256=hashlib.sha256(np.array(data["state"],dtype="<f8").tobytes()).hexdigest(),
            history_sha256=hashlib.sha256(np.array(data["history"],dtype="<f8").tobytes()).hexdigest(),
            alpha=data["alpha"],time=data["time"],
            norms=data["norms"]["native_stage"],rays=len(rays),native_panels=sum(len(r["panels"]) for r in rays),
            profile_counts=dict(Counter(r["profile"] for r in rays)),
            cap_counts=dict(Counter(r["cap"] for r in rays)),
            minimum_limiter_gap=min(r["minimum_limiter_gap"] for r in rays),
            limiter_pattern_sha256=digest([r["limiter"] for r in rays]),
            panels_sha256=digest([r["panels"] for r in rays]),
            maximum_endpoint_displacement=max((r["maximum_endpoint_displacement"] or 0 for r in changes)),
            changed_rays=changed))
    args.output.write_text(json.dumps(dict(scope=__doc__,command=sys.argv,thread_limits=THREADS,
        numpy_version=np.__version__,comparisons=results,
        acceptance_inferred=False,stationarity_inferred=False),indent=2,allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
