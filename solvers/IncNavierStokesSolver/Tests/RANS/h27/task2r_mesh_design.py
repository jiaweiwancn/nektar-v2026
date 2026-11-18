#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Allocate the measured steady-solver cell budget to the streamwise resolution."""

import argparse
from functools import lru_cache
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import Legendre

from coupled_fixtures import mesh
from fixtures import modified
from physical_oracle import load_reference
from physical_reference_refinement import axes
from resources import atomic, identity
from task2r_reference import B, C, D0, H, HS, KIN, RT, require


def design(table, base, fine, witness, out, max_cells):
    ref = load_reference(table)
    x, base_z = axes(base)
    _, fine_z = axes(fine)
    radii = sorted(set(round(z-D0, 14) for z in base_z+fine_z))
    # Retain both of the finely resolved one-sided matching-plane intervals
    # and the wall model bands used in the independent comparison.
    required = [0., round(ref.inlet_cutoff, 14), B, C, round(HS, 14), H, RT]
    required.extend(r for r in radii if .128 <= r <= .164)
    required = sorted(set(required))
    radii = sorted(set(radii + required))

    p = Legendre.basis(11)
    nodes = np.r_[-1., p.deriv().roots(), 1.]
    weights = 2/(12*11*p(nodes)**2)
    basis, _ = modified(3, nodes)
    free = basis[:, 2:]
    mass = free.T @ (weights[:, None]*free)
    grid = np.linspace(-1., 1., 101)
    b, dr = modified(3, grid)

    @lru_cache(maxsize=None)
    def errors(lo, hi):
        # Check a projection against the frozen prescribed physical fields
        # using the same 12-point trace rule and an independent radial grid.
        trace = (lo+hi)/2+(hi-lo)*nodes/2
        sample = (lo+hi)/2+(hi-lo)*grid/2
        trace[[0, -1]], sample[[0, -1]] = [lo, hi], [lo, hi]
        values, target = ref.inlet(trace), ref.inlet(sample)
        result = {}
        for key, scale in (("U", 1.), ("k", KIN)):
            y = values[key]
            edge = y[[0, -1]]
            coeffs = np.r_[edge, np.linalg.solve(mass, free.T@(weights*(y-basis[:, :2]@edge)))]
            result[key] = float(max(abs(b@coeffs-target[key]))/scale)
            if key == "U":
                result["G"] = float(max(abs(dr@coeffs*2/(hi-lo)-target["G"])))
        return result

    limits = dict(U=4e-6, G=7e-4, k=5e-5)
    selected = [0.]
    for lo, hi in zip(required[:-1], required[1:]):
        subset = [r for r in radii if lo <= r <= hi]
        best = {lo: (0, 0., [lo])}
        for upper in subset[1:]:
            admissible = []
            for lower, (n, peak, history) in best.items():
                e = errors(lower, upper)
                gate = max(value/limits[field] for field, value in e.items())
                if gate <= 1:
                    admissible.append((n+1, max(peak, gate), history+[upper]))
            if admissible:
                best[upper] = min(admissible)
        require(hi in best, f"no independently admissible radial segmentation in {lo,hi}")
        selected.extend(best[hi][2][1:])

    nr = len(selected)-1
    allowed = max_cells//nr
    require(allowed >= 56, f"only {allowed} x intervals fit the measured resource control")
    # Keep the base x edges and split the widest current interval first.
    # Widths remain comparable while the analysis retains every frozen station.
    while len(x) - 1 < allowed:
        i = max(range(len(x)-1), key=lambda j: x[j+1]-x[j])
        x.insert(i+1, (x[i]+x[i+1])/2)
    out.mkdir()
    template = out/"designed-template.xml"
    mesh(template, x, [r+D0 for r in selected], quadrature=12)
    atomic(out/"forecast.json", dict(scientific_witness=identity(witness), base_mesh=identity(base),
        fine_mesh=identity(fine), designed_mesh=identity(template), table=identity(table), x=x, r=selected,
        cells=[len(x)-1, nr], memory_design_control_cells=max_cells,
        retained_matching_owners=required, independent_projection_limits=limits,
        maximum_predicted_projection={k:max(errors(a,b)[k] for a,b in zip(selected[:-1],selected[1:])) for k in limits},
        reason="Use the resolved one-sided H intervals and projected inlet error reserves while directing the admissible cell count to the unresolved streamwise comparisons. Verify the actual native projection and all field and weak gates."))
    print("allocated native cells", len(x)-1, nr)


def redistribute_outlet(table, designed, coarse_x, localization, out, max_cells):
    x, z = axes(designed)
    original_x, _ = axes(coarse_x)
    witness = json.loads(localization.read_text())
    old_points = {round(v, 12) for v in original_x}

    def bisect(edges, bands):
        return sorted(edges+[(a+b)/2 for a,b in zip(edges[:-1],edges[1:])
                             if any(lo <= a < b <= hi for lo,hi in bands)])

    bands = [(8.75, 11.25), (12.5, 15.)]
    tip = [(14.375, 15.)]
    refined = bisect(bisect(x, bands), tip)
    # The localization comparison independently found k/Kin RMS below
    # 1.1e-4 in these upstream bands. Retain their entire 56-cell parent
    # grid, plus as many additional points as the measured memory permits.
    quiet = [v for v in witness["inner_k_bands"] if v["x"][1] <= 4.75
             and v["RMS_over_Kin"] <= 1.1e-4 and v["x"][0] >= -8.75]
    quiet.sort(key=lambda row: row["RMS_over_Kin"])
    candidates = [v for row in quiet for v in x
                  if row["x"][0] < v < row["x"][1] and round(v,12) not in old_points]
    allowed = max_cells//(len(z)-1)
    remove = max(0, len(refined)-1-allowed)
    require(len(candidates)>=remove, "quiet bands cannot free the measured cells without removing parent edges")
    removed = candidates[:remove]
    redistributed = [v for v in refined if v not in removed]
    require(len(redistributed)-1 <= allowed and set(original_x) <= set(redistributed),
            "streamwise parent mesh or memory design was lost")
    out.mkdir()
    session = out/"designed-template.xml"
    mesh(session,redistributed,z,quadrature=12)
    atomic(out/"forecast.json", dict(scientific_witness=identity(localization), base_mesh=identity(designed),
         coarse_x_mesh=identity(coarse_x), designed_mesh=identity(session), table=identity(table),
         cells=[len(redistributed)-1,len(z)-1], maximum_cells=max_cells, x=redistributed, z=z,
         removed_upstream=removed, observed_quiet_bands=quiet, split_bands=bands, additional_terminal_splits=tip,
         scientific_checks="The design only reallocates cells. Compare every region, owner, and frozen station against the prior uniform-axis witness."))
    print("redistributed outlets", len(redistributed)-1, len(z)-1, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("table", "base", "fine", "witness", "out"):
        parser.add_argument("--"+name, type=Path, required=True)
    parser.add_argument("--outlet-axis", action="store_true")
    parser.add_argument("--max-cells", type=int, required=True)
    a = parser.parse_args()
    operation = redistribute_outlet if a.outlet_axis else design
    operation(a.table, a.base, a.fine, a.witness, a.out, a.max_cells)


if __name__ == "__main__":
    main()
