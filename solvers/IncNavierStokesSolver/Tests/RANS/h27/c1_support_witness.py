#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compare solved WJ fields at the same one-sided support coordinates."""

import argparse
import csv
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq

from c1_margins import response
from c1_oracle import Audit, load_c1_state
from c1_rejected_points import point
from c1_support_controls import coincides, contains, locations
import physical_oracle as po
from resources import identity


FIELD_NAMES = ("margin", "Gxz", "divergence", "pressure", "k_load", "convection")
FIXED_AXIAL = (-.05, -.02, -.005, -.001, 0., .001, .005, .02, .05)
FIXED_RADIAL = (-.0015, -.001, -.00075, -.0005, -.00025, -.0001, -.00005,
                0., .00005, .0001, .00025, .0005, .00075, .001, .0015)
EXTRA_RADIAL = (.002, .0025, .003, .004, .006, .010, .015)
RELATIVE_OFFSETS = (.001, .01, .05, .20)


def owner_side(lo, hi, value, H):
    if coincides(value, lo, H):
        return 1
    if coincides(value, hi, H):
        return -1
    return 0


def section_points(audits, wide):
    first = next(iter(audits.values()))
    L, h, H = first.hill.L, first.ref.h, first.ref.H
    section = {label: a for label, (a, _) in locations(first.geometry).items()
               if label.startswith("support_")}
    left = 0. if wide else .35*L
    if wide:
        section["crown"] = left
    # The field-independent part covers the downstream support segment. Any
    # mesh endpoint in this segment is also evaluated on every other mesh.
    edges = sorted({round(float(v), 13) for audit in audits.values()
                    for lo, hi in audit.columns for v in (lo, hi) if left <= v <= L})
    for v in edges:
        if not any(coincides(v, a, H) for a in section.values()):
            section[f"a={v:.13g}"] = v
    a = sorted(set(np.linspace(left, L, 153 if wide else 129)) |
               {a+off for a in section.values() for off in FIXED_AXIAL if left <= a+off <= L})
    radial = sorted((*FIXED_RADIAL, *EXTRA_RADIAL, *(-v for v in EXTRA_RADIAL))) if wide else FIXED_RADIAL
    r = np.asarray([h+off for off in radial])
    fixed = np.column_stack([v.ravel() for v in np.meshgrid(a, r)])
    return section, fixed, radial


def write_profiles(writer, audit, label, fixed, section):
    h, H = audit.ref.h, audit.ref.H
    records = {}
    relative = []
    for cell in audit.mesh:
        xl, xr, rl, rr = cell["bounds"]
        if cell["epsilon"] or xr < fixed[:, 0].min() or xl > fixed[:, 0].max():
            continue
        audit.clock()
        common = np.array([contains(cell, a, r, H) for a, r in fixed])
        ax, rad = fixed[common].T
        count = len(ax)
        scaled = []
        for a in section.values():
            if not contains(cell, a, h, H):
                continue
            for x_sign, r_sign in ((-1, -1), (-1, 1), (1, -1), (1, 1)):
                span_x = a-xl if x_sign < 0 else xr-a
                span_r = h-rl if r_sign < 0 else rr-h
                if span_x > 1e-13 and span_r > 1e-13:
                    scaled.extend((a+x_sign*t*span_x, h+r_sign*t*span_r, a, x_sign*t, r_sign*t)
                                  for t in RELATIVE_OFFSETS)
        if scaled:
            a_rel, r_rel, *_ = np.asarray(scaled).T
            ax, rad = np.r_[ax, a_rel], np.r_[rad, r_rel]
        if not len(ax):
            continue
        metric = audit.hill.metrics(cell, ax, rad, outer=False)
        _, U, G, _, _, _, scalars, _, _, _ = audit.sample(cell, ax, rad, metric)
        terms = audit.matching_load(cell, ax)
        metric["K"] = terms["K"]
        values, uncertainty = response(audit.ref, G, scalars["k"], scalars["epsilon"], rad,
                                       metric, False, pointwise=True)
        for i, (a, r) in enumerate(zip(ax, rad)):
            x_owner, r_owner = owner_side(xl, xr, a, H), owner_side(rl, rr, r, H)
            values_at_point = dict(zip(FIELD_NAMES, [float(values[i]), float(G[i, 0, 1]),
                float(np.trace(G[i])), float(terms["pressure"][i]), float(terms["k"][i]),
                float(terms["convection"][i])]))
            row = dict(mesh=label, a=float(a), r=float(r), x=float(metric["position"][i, 0]),
                       z=float(metric["position"][i, 1]), owner=cell["geometry"],
                       axial_side=x_owner, radial_side=r_owner, fixed=i < count, r_minus_h=float(r-h),
                       basis=cell["basis"] if "basis" in cell else audit.data["basis"],
                       u=float(U[i, 0]), w=float(U[i, 1]), Gxx=float(G[i, 0, 0]),
                       Gzx=float(G[i, 1, 0]), Gzz=float(G[i, 1, 1]), k=float(scalars["k"][i]),
                       epsilon=float(scalars["epsilon"][i]), Pi_H=float(terms["Pi_H"][i]),
                       K=float(terms["K"][i]), u_t_H=float(terms["u_t"][i]),
                       angular_difference=float(uncertainty[i]), **values_at_point)
            if i < count:
                records.setdefault((float(a), float(r)), []).append(row)
            else:
                _, _, source, fx, fr = scaled[i-count]
                row.update(relative_origin=source, local_axial_fraction=fx, local_radial_fraction=fr)
                relative.append(row)
            writer.writerow(row)
    if len(records) != len(fixed):
        raise ValueError("field does not cover the common physical support grid")
    return records, relative


def limits(audit, section):
    out = {}
    h, H = audit.ref.h, audit.ref.H
    for name, a in section.items():
        rows = []
        for cell in audit.mesh:
            if not contains(cell, a, h, H):
                continue
            row = point(audit, cell, a, h, "support")
            xl, xr, rl, rr = cell["bounds"]
            row["axial_side"] = owner_side(xl, xr, a, H)
            row["radial_side"] = owner_side(rl, rr, h, H)
            row["load_limit_from_below"] = {
                name: float(values[0]) for name, values in audit.matching_load(cell, np.array([a])).items()}
            ax, rad = np.array([a]), np.array([h])
            metric = audit.hill.metrics(cell,ax,rad,outer=False)
            _, _, _, pressure, stress, _, _, diffusivity, gradients, _ = audit.sample(cell,ax,rad,metric)
            sigma = stress[0]-pressure[0]*np.eye(2)
            k_flux = -diffusivity["k"][0]*gradients["k"][0]
            row["raw_tractions"] = {face:(sigma@metric[n][0]).tolist()
                                     for face,n in (("axial","t"),("radial","n"))}
            row["raw_k_diffusive_flux"] = {face:float(k_flux@metric[n][0])
                                           for face,n in (("axial","t"),("radial","n"))}
            rows.append(row)
        out[name] = rows
    return out


def choose(rows, x, r):
    candidates = [p for p in rows if p["axial_side"] in (0, x) and p["radial_side"] in (0, r)]
    if len(candidates) != 1:
        raise ValueError("ambiguous or missing approached one-sided owner")
    return candidates[0]


def raw_jumps(points):
    result = {}
    for name, rows in points.items():
        jumps = []
        for face, pairs in (("axial",[((-1,r),(1,r)) for r in (-1,1)]),
                            ("radial",[((x,-1),(x,1)) for x in (-1,1)])):
            for lo,hi in pairs:
                minus,plus=choose(rows,*lo),choose(rows,*hi)
                if minus["geometry"]==plus["geometry"]:
                    continue
                jumps.append(dict(face=face, owners=[minus["geometry"],plus["geometry"]],
                    traction=(np.asarray(plus["raw_tractions"][face])-
                              minus["raw_tractions"][face]).tolist(),
                    k_diffusive_flux=(plus["raw_k_diffusive_flux"][face]-
                                      minus["raw_k_diffusive_flux"][face])))
        result[name]=jumps
    return result


def extent(rows, relative):
    flat = [r for items in rows.values() for r in items]
    negative = lambda p: p["margin"]+p["angular_difference"] < 0.
    bad = [r for r in flat if negative(r)]
    off_face = [r for r in bad if r["axial_side"]==r["radial_side"]==0]

    def window(points):
        if not points:
            return dict(count=0)
        worst = min(points, key=lambda p: p["margin"])
        return dict(count=len(points), a_range=[min(p["a"] for p in points), max(p["a"] for p in points)],
                    r_minus_h_range=[min(p["r_minus_h"] for p in points), max(p["r_minus_h"] for p in points)],
                    x_range=[min(p["x"] for p in points), max(p["x"] for p in points)],
                    z_range=[min(p["z"] for p in points), max(p["z"] for p in points)],
                    worst={key: worst[key] for key in ("a", "r", "margin", "owner", "Gxz", "divergence",
                                                       "angular_difference", "axial_side", "radial_side")})

    strip = {}
    for v in sorted({p["r_minus_h"] for p in flat}):
        hit = [p for p in bad if p["r_minus_h"] == v]
        if hit:
            strip[f"{v:.12g}"] = window(hit)
    # Report sampled extents, including all one-sided limits. These are finite
    # sample spans and can include gaps between adverse samples.
    return dict(minimum=window(flat)["worst"], adverse=window(bad), strictly_interior=window(off_face),
                adverse_below=window([p for p in off_face if p["r_minus_h"] < 0.]),
                adverse_above=window([p for p in off_face if p["r_minus_h"] > 0.]),
                support_line=window([p for p in bad if p["r_minus_h"] == 0.]),
                relative_offsets=window([p for p in relative if negative(p)]), strip=strip,
                maximum_directional_difference=max(p["angular_difference"] for p in flat))


def directional_effects(section, profiles, h):
    result = {}
    for name, a in section.items():
        rows = []
        for sx, sr in ((-1,-1),(-1,1),(1,-1),(1,1)):
            parent, axial, radial, both = [
                choose(profiles[key][a,h], sx, sr) for key in ("M0", "Mx", "Mr", "Mxr")]
            rows.append(dict(approach=[sx,sr], owners=[p["owner"] for p in (parent, axial, radial, both)],
                values={key: [p[key] for p in (parent, axial, radial, both)] for key in FIELD_NAMES},
                axial={key: axial[key]-parent[key] for key in FIELD_NAMES},
                radial={key: radial[key]-parent[key] for key in FIELD_NAMES},
                interaction={key: both[key]-axial[key]-radial[key]+parent[key] for key in FIELD_NAMES}))
        result[name] = rows
    return result


def variation(profiles, pairs):
    results = {}
    approaches = ((-1,-1),(-1,1),(1,-1),(1,1))
    for a, b in pairs:
        fields = {key: 0. for key in FIELD_NAMES}
        changes = common_adverse = 0
        adverse_difference = 0.
        smallest = math.inf
        for xy, rows in profiles[a].items():
            for sx, sr in approaches:
                A, B = choose(rows, sx, sr), choose(profiles[b][xy], sx, sr)
                for key in fields:
                    fields[key] = max(fields[key], abs(A[key]-B[key]))
                changes += int((A["margin"] < 0.) != (B["margin"] < 0.))
                if max(A["margin"]+A["angular_difference"],
                       B["margin"]+B["angular_difference"]) < 0.:
                    common_adverse += 1
                    adverse_difference = max(adverse_difference, abs(A["margin"]-B["margin"]))
                smallest = min(smallest, abs(A["margin"]), abs(B["margin"]))
        results[f"{a}/{b}"] = dict(maximum_absolute_difference=fields, sampled_sign_changes=changes,
                                   smallest_absolute_margin=smallest,
                                   approach_samples=4*len(profiles[a]),
                                   common_adverse_approaches=common_adverse,
                                   maximum_margin_difference_where_both_adverse=adverse_difference)
    return results


def sampled_action_order(audit):
    """Independent signed row check of the fixed sampled-force support rule."""
    names = ("pressure", "k", "convection")

    def assemble(order):
        result = {k: np.zeros((audit.n,2)) for k in names}
        nodes, weights = np.polynomial.legendre.leggauss(order)
        for cell in audit.mesh:
            audit.clock()
            xl,xr,rl,rr = cell["bounds"]
            upper = min(rr,audit.ref.h)
            if upper <= rl:
                continue
            a = xl+(xr-xl)/2*(nodes+1)
            r = rl+(upper-rl)/2*(nodes+1)
            ax,rad = [v.ravel() for v in np.meshgrid(a,r)]
            metric = audit.hill.metrics(cell,ax,rad,outer=False)
            shape = po.basis(cell,"u",ax,rad)[0]
            terms = audit.matching_load(cell,a)
            measure = np.outer(weights,weights).ravel()*(xr-xl)*(upper-rl)/4*metric["det"]
            for name in names:
                force = measure*np.tile(terms[name],order)
                result[name][cell["u"]] += shape.T @ (force[:,None]*metric["t"])
        return result

    a, b = assemble(16), assemble(28)
    return {name: dict(order16_max=float(abs(a[name]).max()), order28_max=float(abs(b[name]).max()),
                       maximum_difference=float(abs(a[name]-b[name]).max())) for name in names}


def attribution(audit, first, second):
    ref = audit.ref
    r = first["rho"]
    radius = np.array([r])
    t = np.array(first["normal"])[[2, 0]]*np.array([1., -1.])
    metric = {"n": np.array([first["normal"]])[..., [0, 2]], "t": t[None, :]}
    P = np.ix_([0,2],[0,2])
    A, B = (np.asarray(p["gradient"])[P] for p in (first, second))
    k0, k1, K0, K1 = first["k"], second["k"], first["K"], second["K"]
    E0, E1 = first["epsilon"], second["epsilon"]

    def dissipation(K):
        return K**1.5/ref.cl/max(r, float(ref.cutoff(np.array([K]))[0]))

    def evaluate(g, k, K, E):
        metric["K"] = np.array([K])
        values, uncertainty = response(ref,g[None, :, :],np.array([k]),np.array([E]),
                                       radius,metric,False,pointwise=True)
        return dict(response=float(values[0]), difference=float(uncertainty[0]),
                    k_over_epsilon=k/E)

    trials = dict(parent=evaluate(A,k0,K0,E0), comparison=evaluate(B,k1,K1,E1),
                  gradient_only=evaluate(B,k0,K0,E0),
                  scalars_with_rebuilt_dissipation=evaluate(A,k1,K1,dissipation(K1)),
                  isolated_epsilon=evaluate(A,k0,K0,E1),
                  compatible_matching_K=evaluate(A,k0,K1,dissipation(K1)))
    solenoidal = B-np.eye(2)*(np.trace(B)/2)
    trials["local_planar_zero_divergence"]=evaluate(solenoidal,k1,K1,E1)

    paths = {}
    for mode in ("full_gradient_fixed_parent_scalars", "full_gradient_and_compatible_scalars"):
        def at(s):
            g = A+s*(B-A)
            k, K, E = k0, K0, E0
            if mode.endswith("compatible_scalars"):
                k, K = k0+s*(k1-k0), K0+s*(K1-K0)
                E = dissipation(K)
            return evaluate(g,k,K,E)["response"]
        sample = [(float(s),at(float(s))) for s in np.linspace(0,1,17)]
        crossing = next((brentq(at,s,t,xtol=2e-10) for (s,v),(t,w) in zip(sample,sample[1:])
                         if v*w < 0.), None)
        distance = np.linalg.norm(B-A)
        paths[mode] = dict(samples=sample, crossing_fraction=crossing,
                           observed_gradient_change=float(distance),
                           gradient_distance_along_path=None if crossing is None else float(crossing*distance))
    return dict(parent_owner=first["geometry"], comparison_owner=second["geometry"],
                location=[first["a"],r], substitutions=trials, paths=paths)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--field", nargs=2, action="append", required=True, metavar=("LABEL","PACKET"))
    parser.add_argument("--compare", nargs=2, action="append", default=[], metavar=("A","B"))
    parser.add_argument("--attribute", nargs=2, action="append", default=[], metavar=("A","B"))
    parser.add_argument("--table", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--profiles", type=Path, required=True)
    parser.add_argument("--seconds", type=float, default=170)
    parser.add_argument("--check-sampled-order", action="store_true")
    parser.add_argument("--wide", action="store_true",
                        help="include the full downstream half and wider radial offsets")
    args = parser.parse_args()
    end = time.monotonic()+args.seconds
    ref = po.load_reference(args.table)
    audits, identity_rows = {}, {}
    for label, filename in args.field:
        p = Path(filename)
        if label in audits:
            raise ValueError("duplicate field label")
        data, metadata = load_c1_state(p, True)
        if not coincides(data["geometry"]["h"],ref.h,ref.H) or metadata["regularization"]["eta"] != 0.:
            raise ValueError("support observations require one pinned original reference and literal eta=0")
        audits[label] = Audit(data,ref,28,end)
        identity_rows[label] = dict(packet=identity(p), coefficient_sha256=data["state_sha256"],
                                    native_original_residuals=metadata["original_residuals"])
    geometry_keys = ("h","H","d0","hill_height","hill_half_length","normal_map_signature")
    if len({tuple(a.geometry[k] for k in geometry_keys) for a in audits.values()}) != 1:
        raise ValueError("these comparisons require the identical physical hill chart")
    section, fixed, radial = section_points(audits,args.wide)
    profiles, summaries, controls = {}, {}, {}
    with args.profiles.open("x") as f:
        columns = ["mesh","a","r","x","z","owner","axial_side","radial_side","fixed","r_minus_h","basis",
                   "u","w","Gxx","Gzx","Gzz","k","epsilon","Pi_H","K","u_t_H","angular_difference",
                   *FIELD_NAMES,"relative_origin","local_axial_fraction","local_radial_fraction"]
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        for label, audit in audits.items():
            profiles[label], relative = write_profiles(writer,audit,label,fixed,section)
            summaries[label] = extent(profiles[label],relative)
            controls[label] = limits(audit,section)
    result = dict(table=identity(args.table), audit_source=identity(__file__), fields=identity_rows,
                  chart_support_h=ref.h, fixed_axial_offsets=FIXED_AXIAL, fixed_radial_offsets=radial,
                  axial_grid_range=[float(fixed[:,0].min()),float(fixed[:,0].max())],
                  cell_relative_fractions=RELATIVE_OFFSETS, fixed_grid_points=len(fixed), controls=controls,
                  signed_extents=summaries, matching_point_variation=variation(profiles,args.compare),
                  support_face_jumps={label:raw_jumps(points) for label,points in controls.items()},
                  profile_file=identity(args.profiles))
    if args.check_sampled_order:
        result["independent_sampled_action_16_vs_28"]={label:sampled_action_order(audit)
                                                        for label,audit in audits.items()}
    if {"M0","Mx","Mr","Mxr"} <= profiles.keys():
        result["directional_effects"]=directional_effects(section,profiles,ref.h)
    result["attribution"]={}
    for a,b in args.attribute:
        rows=[]
        for label in section:
            for x, r in ((-1,-1),(-1,1),(1,-1),(1,1)):
                first, second = (choose(controls[k][label],x,r) for k in (a,b))
                # Store each owner pair once, including all distinct corner limits.
                key=(label,first["geometry"],second["geometry"])
                if any(row["key"]==key for row in rows):
                    continue
                rows.append(dict(key=key, **attribution(audits[a],first,second)))
        result["attribution"][f"{a}/{b}"]=rows
    with args.out.open("x") as f:
        json.dump(result,f,indent=2,allow_nan=False)
        f.write("\n")


if __name__ == "__main__":
    main()
