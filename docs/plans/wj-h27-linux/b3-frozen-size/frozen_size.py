#!/usr/bin/env python3
"""Propose radius-constrained corrections and summarize fresh native checks."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import time
import warnings

import numpy as np
from scipy.optimize import OptimizeWarning, linprog

from c1_action_gap import read_audit

RADII = (0.00020533, .005, .025, .1, 1, 10, 100, 500, 1000, 2500, 5000, 10000, 17055, 30000)
BLOCKS = ("momentum", "continuity", "gauge", "k", "outer_epsilon", "trace")
UNKNOWNS = ("u", "w", "Pi", "gauge", "k", "outer_epsilon")
BINDINGS = ("base_coefficient_sha256", "reference_identity", "free_id_sha256", "free_counts",
            "quadrature", "correction_mass", "fixed_k_tangent")


def pin(path):
    return {"path": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def write_json(path, data):
    with path.open("x") as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write("\n")


def array_spec(metadata, name):
    return next(item for item in metadata["arrays"] if item["name"] == name)


def vector_norms(parts, coefficients):
    return np.array([np.linalg.norm(part @ coefficients) for part in parts])


class JointModel:
    """Cutting-plane solve of the joint model over one fixed scaled subspace."""
    def __init__(self, q, parts):
        self.q = q
        self.raw_parts = parts
        stacked = np.vstack([p[:, 1:] for p in parts])
        _, singular, vt = np.linalg.svd(stacked, full_matrices=False)
        assert singular.size == q.shape[1] and singular[0] > 0
        self.ranks = singular
        self.scale = np.maximum(singular, singular[0]*1.e-13)
        self.transform = vt.T / self.scale
        self.parts = [(p[:, 0], p[:, 1:] @ self.transform) for p in parts]
        self.rows = q @ self.transform
        self.c_norm = np.linalg.norm(np.concatenate([b for b, _ in self.parts]))
        self.c0 = max(np.linalg.norm(b) for b, _ in self.parts)
        self.norm_cuts = []
        self.add_norm_cuts(np.zeros(q.shape[1]))
        self.initial_rows = set(np.argmax(np.abs(self.rows), axis=0).tolist())

    def add_norm_cuts(self, v):
        for b, a in self.parts:
            r = b + a @ v
            norm = np.linalg.norm(r)
            if norm > 1.e-15:
                normal = r/norm
                self.norm_cuts.append((normal @ a, -float(normal @ b)))

    def solve(self, radius, previous=None, max_iterations=100, *, gap_tolerance=2.e-4,
              relative_gap=None, seconds=None, incumbents=()):
        n = self.q.shape[1]
        tcap = self.c0
        zero = np.zeros(n)
        best_v = zero.copy()
        best = tcap
        best_size = 0.
        lower = None
        best_origin = "zero"
        active = self.initial_rows.copy()
        radius_margin = max(1.e-12*radius, 1.e-14)
        started = time.monotonic()
        model_floor = 8.e-9
        target = gap_tolerance
        stop = "ITERATION_CAP"
        lp_status = None
        stages = []

        def consider(proposal):
            v = np.asarray(proposal)
            size = float(np.max(np.abs(self.rows @ v)))
            if size > radius:
                v = v * ((radius-radius_margin)/size)
                size = float(np.max(np.abs(self.rows @ v)))
            scores = [np.linalg.norm(b+a @ v) for b, a in self.parts]
            quality = max(scores)
            if quality < best:
                return v.copy(), quality, size
            return None

        candidates = ([("previous", previous)] if previous is not None else []) + list(incumbents)
        for name, vector in candidates:
            feasible = consider(vector)
            if feasible:
                best_v, best, best_size = feasible
                best_origin = name
                self.add_norm_cuts(best_v)
        incumbent_score = best
        incumbent_origin = best_origin

        # Both bounds follow from orthonormal Q and ||A y|| <= ||c||+sqrt(7)*tcap.
        by_radius = math.sqrt(self.q.shape[0])*radius*self.scale
        by_merit = np.divide((self.c_norm+math.sqrt(len(self.parts))*tcap)*self.scale,
                             self.ranks, out=np.full(n, np.inf), where=self.ranks > 0)
        caps = np.minimum(by_radius, by_merit)*(1.+1.e-12)
        bounds = [(-b, b) for b in caps]+[(0., tcap*(1.+1.e-12))]
        iterations = 0
        for iteration in range(1, max_iterations+1):
            limits = {"primal_feasibility_tolerance":1.e-9,
                      "dual_feasibility_tolerance":1.e-9, "threads":1}
            if seconds is not None:
                remaining = seconds-(time.monotonic()-started)
                if remaining <= 0.:
                    stop = "TIME_CAP"
                    break
                limits["time_limit"] = remaining
            indices = list(sorted(active))
            rows = self.rows[indices]
            norm = np.maximum(radius, np.max(np.abs(rows), axis=1))
            cuts = np.vstack([row for row, _ in self.norm_cuts])
            Aub = np.vstack((np.column_stack((cuts, -np.ones(len(cuts)))),
                             np.column_stack((rows/norm[:, None], np.zeros(len(rows)))),
                             np.column_stack((-rows/norm[:, None], np.zeros(len(rows))))))
            bub = np.concatenate(([bound for _, bound in self.norm_cuts], radius/norm, radius/norm))
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=OptimizeWarning, message="(Unrecognized options detected|Unknown solver options)")
                result = linprog(np.r_[zero, 1.], A_ub=Aub, b_ub=bub, bounds=bounds,
                                 method="highs-ds", options=limits)
            iterations = iteration
            lp_status = int(result.status)
            if not result.success:
                if result.status == 1 and seconds is not None:
                    stop = "TIME_CAP"
                    break
                if result.status != 4:
                    raise ValueError(f"radius {radius:g}: {result.message}")
                stop = "LP_LIMIT"
                break
            v = result.x[:-1]
            lower = max(0. if lower is None else lower, float(result.x[-1]))
            feasible = consider(v)
            if feasible:
                best_v, best, best_size = feasible
                best_origin = f"LP iteration {iteration}"
            if relative_gap is not None:
                target = min(gap_tolerance, relative_gap*max(tcap-best, model_floor))
            if iteration in (100, 300, max_iterations):
                stages.append({"iteration":iteration, "upper":float(best), "lower":lower, "target":target})
            if best-lower <= target:
                stop = "GAP_TARGET"
                break
            scores = np.abs(self.rows @ v)
            violating = np.flatnonzero(scores > radius+radius_margin)
            if violating.size:
                top = np.argsort(scores[violating], kind="stable")[-12:]
                active.update(violating[top].tolist())
            self.add_norm_cuts(v)
        y = self.transform @ best_v
        scores = vector_norms(self.raw_parts, np.r_[1., y])
        error = max(scores)
        return best_v, y, {"radius":radius, "size":best_size, "model_score":float(error),
                           "model_lower":lower, "model_gap":None if lower is None else float(error-lower),
                           "model_iterations":iterations, "model_size_rows":len(active),
                           "model_ratios":scores.tolist(), "model_gap_target":target,
                           "model_floor":model_floor, "model_stages":stages, "status":stop,
                           "model_lp_status":lp_status, "model_seconds":time.monotonic()-started,
                           "initial_incumbent":incumbent_origin, "initial_incumbent_score":float(incumbent_score),
                           "selected_incumbent":best_origin}


def load_model(path):
    meta, arrays = read_audit(path)
    assert meta["stage"] == "frozen-size-basis" and meta["factor_serial"] == 1
    targets = np.array([meta["measurements"][name+"_target"] for name in ("physical", "scaled")])
    q, response = (arrays[k] for k in ("scaled_correction_basis", "raw_residual_and_actions"))
    sx = arrays["unknown_scale"][:, 0]
    sr = arrays["row_scale"][:, 0]
    assert np.linalg.norm(q.T @ q-np.eye(q.shape[1])) < 1.e-10
    base = array_spec(meta, "raw_base_residual")
    forcing = meta["measurements"]["forcing"]
    assert forcing == .3 and np.array_equal(response[:, 0], arrays["raw_base_residual"][:, 0])
    np.testing.assert_allclose(targets, forcing*np.array([max(base["physical_blocks"]), base["scaled_l2"]]),
                               rtol=1.e-13, atol=0)
    assert meta["correction_mass"] == 0 and meta["fixed_k_tangent"]
    m = arrays["response_column_normalization"][:, 0]
    parts = []
    eigenchecks = []
    for i in range(len(BLOCKS)):
        gram = arrays[f"normalized_physical_gram_{i}"]
        eigenvalues, modes = np.linalg.eigh(gram)
        peak = max(float(eigenvalues[-1]), 1.e-300)
        if eigenvalues[0] < -1.e-10*peak:
            raise ValueError(f"{BLOCKS[i]} violates Hilbert polarization")
        eigenchecks.append([float(eigenvalues[0]), float(eigenvalues[-1])])
        part = (np.sqrt(np.maximum(eigenvalues, 0.))[:, None]*modes.T)*m[None, :]/targets[0]
        parts.append(part)
    parts.append(np.linalg.qr(sr[:, None]*response/targets[1], mode="r"))
    norms = vector_norms(parts, np.r_[1., np.zeros(q.shape[1])])
    np.testing.assert_allclose(norms, np.r_[np.array(base["physical_blocks"])/targets[0], 1./forcing],
                               rtol=2.e-10, atol=1.e-13)
    return meta, arrays, sx, q, targets, JointModel(q, parts), eigenchecks


def propose(args):
    meta, arrays, sx, q, targets, model, eigenchecks = load_model(args.basis)
    old, oa = read_audit(args.archived)
    assert all(old[key] == meta[key] for key in BINDINGS)
    assert all(np.array_equal(oa[key], arrays[key]) for key in ("base", "unknown_scale", "row_scale"))
    directions, candidates = [], []
    previous = None
    radii = args.radii or RADII
    for radius in radii:
        previous, y, item = model.solve(radius, previous)
        direction = sx*(q @ y)
        np.testing.assert_allclose(item["size"], np.max(np.abs(direction/sx)), rtol=5.e-12, atol=1.e-15)
        item.update(label="joint_model", proposed_scaled_l2=float(np.linalg.norm(y)))
        directions.append(direction)
        candidates.append(item)
    matrix = np.column_stack(directions).astype("<f8", copy=False)
    datafile = args.output/"candidates.bin"
    data = matrix.tobytes(order="F")
    with datafile.open("xb") as stream:
        stream.write(data)
    for index, item in enumerate(candidates):
        item.update(index=index, direction_sha256=hashlib.sha256(matrix[:, index].tobytes()).hexdigest())
    write_json(args.output/"candidates.json", {
        "basis": pin(args.basis), "archived_direction": pin(args.archived), "input": pin(datafile),
        "response_singular_values": model.ranks.tolist(), "physical_gram_eigenvalue_extrema": eigenchecks,
        "forcing": meta["measurements"]["forcing"], "targets": targets.tolist(),
        "basis_receipt_identity": {key:meta[key] for key in (*BINDINGS, "signature", "equation_identity",
                                  "source_identity", "factor", "session_sha256", "current_fd")},
        "basis_rank": q.shape[1], "frozen_base": array_spec(meta, "raw_base_residual"), "candidates": candidates,
    })


def summarize(args):
    manifest = json.loads((args.root/"candidates.json").read_text())
    verified = sorted(args.case.glob("state.action-frozen-size-candidate-*.json"),
                      key=lambda path: json.loads(path.read_text())["proposal"])
    targets = np.array(manifest["targets"])
    assert len(verified) == len(manifest["candidates"])
    reports, repeat_shas = [], []
    for path, item in zip(verified, manifest["candidates"]):
        meta, arrays = read_audit(path)
        assert meta["proposal"] == item["index"] and meta["factor_serial"] == 0
        common = manifest["basis_receipt_identity"]
        for key in ("signature", "equation_identity"):
            if key not in common and not manifest["basis_rank"]:
                common[key] = meta[key]
        assert all(meta[key] == common[key] for key in (*BINDINGS, "signature", "equation_identity",
                                                        "source_identity", "current_fd"))
        np.testing.assert_array_equal(targets, [meta["measurements"][name+"_target"] for name in ("physical", "scaled")])
        assert meta["probe_input_sha256"] == manifest["input"]["sha256"]
        for key in ("unknown_scale", "row_scale"):
            assert np.all(arrays[key] > 0)
        # The probe converts d to scaled coordinates and back when publishing.
        d = arrays["direction"][:, 0]
        with Path(manifest["input"]["path"]).open("rb") as source:
            source.seek(8*d.size*item["index"])
            raw = source.read(8*d.size)
        assert hashlib.sha256(raw).hexdigest() == item["direction_sha256"]
        np.testing.assert_allclose(d, np.frombuffer(raw, dtype="<f8"), rtol=5.e-16, atol=0)
        assert array_spec(meta, "raw_base_residual")["sha256"] == manifest["frozen_base"]["sha256"]
        def norms(name):
            spec = array_spec(meta, name)
            return np.array([max(spec["physical_blocks"]), spec["scaled_l2"]])
        pair = [norms(f"raw_model_residual_{i}").tolist()
                if f"raw_model_residual_{i}" in arrays else None for i in range(3)]
        spacing = [row for row in meta["stencils"] if row["method"]]
        ok = len(spacing) >= 2 and all(pair[i] is not None for i in (0,1)) and (
            spacing[1]["parameter"] < .8*spacing[0]["parameter"])
        allowance = 2.*norms("raw_spacing_change_1") if ok else None
        passed = pair[0] is not None and bool(np.all(np.array(pair[0]) <= targets))
        resolved = allowance is not None and bool(np.all(allowance < .5*targets))
        scaled = d/arrays["unknown_scale"][:, 0]
        segments = np.split(scaled, np.cumsum(meta["free_counts"])[:-1])
        direction_blocks = {key: float(np.max(np.abs(block))) if block.size else 0.
                            for key, block in zip(UNKNOWNS, segments)}
        report = item | {"direction_blocks_scaled_infinity": direction_blocks, "receipt": pin(path), "payload": meta["raw_payload"],
                         "native_scaled_infinity": meta["measurements"]["scaled_size"],
                         "residuals": pair, "physical_blocks": [
                             array_spec(meta, f"raw_model_residual_{i}")["physical_blocks"]
                             if pair[i] is not None else None for i in range(3)],
                         "action_allowances": None if allowance is None else allowance.tolist(),
                         "both_targets_pass": passed, "action_resolved": resolved,
                         "pass_and_resolved": passed and resolved,
                         "both_targets_with_allowance": bool(resolved and np.all(np.array(pair[0])+allowance<=targets)),
                         "stencils": meta["stencils"]}
        if item["label"] in ("archived_R9_direction", "repeat_archived_R9_direction"):
            repeat_shas.append({s["name"]:s["sha256"] for s in meta["arrays"] if s["name"].startswith("raw_")})
        reports.append(report)
    sizes = [r for r in reports if r["pass_and_resolved"]]
    result = {"schema":"WJH27FrozenCorrectionSizeV1", "manifest": pin(args.root/"candidates.json"),
              "targets": targets.tolist(), "numerical_identity": manifest["basis_receipt_identity"],
              "basis_rank": manifest["basis_rank"], "physical_time_advanced":0.,
              "repeat_archived_raw_equal": len(repeat_shas)==2 and repeat_shas[0]==repeat_shas[1],
              "smallest_verified_size": min(sizes, key=lambda row:row["native_scaled_infinity"]) if sizes else None,
              "candidates": reports}
    write_json(args.root/"frozen-size-results.json", result)
    with (args.root/"frozen-size-results.csv").open("x") as stream:
        writer = csv.writer(stream)
        writer.writerow(["candidate", "radius", "scaled_Linf", "scaled_L2", "model_score", "model_lower",
                         "physical_residual", "scaled_residual", "physical_allowance", "scaled_allowance",
                         "nominal_both", "resolved", "both_with_allowance", *UNKNOWNS])
        for row in reports:
            writer.writerow([row["label"], row["radius"], row["native_scaled_infinity"],
                             row["proposed_scaled_l2"], row.get("model_score"), row.get("model_lower"),
                             *(row["residuals"][0] or (None, None)), *(row["action_allowances"] or (None, None)),
                             row["both_targets_pass"], row["action_resolved"], row["both_targets_with_allowance"],
                             *(row["direction_blocks_scaled_infinity"][key] for key in UNKNOWNS)])


def self_test():
    q = np.eye(2)
    parts = [np.array([[2., 1., 0.]]), np.array([[1., 0., 2.]])/2.,
             np.array([[2., 1., 0.], [1., 0., 2.]])/2.]
    model = JointModel(q, parts)
    _, direction, info = model.solve(.5)
    assert abs(np.max(np.abs(direction))-.5) <= 1.e-7
    assert abs(info["model_score"]-1.5)<1.e-6 and info["model_gap"]<1.e-6
    _, direction, info = model.solve(3.)
    assert np.linalg.norm(direction-np.array([-2., -.5]), ord=np.inf)<3.e-4
    assert info["model_score"]<2.e-4


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    make = sub.add_parser("propose")
    for arg in ("basis", "archived", "output"):
        make.add_argument(arg, type=Path)
    make.add_argument("--radii", type=float, nargs="+")
    report = sub.add_parser("summarize")
    for arg in ("root", "case"):
        report.add_argument(arg, type=Path)
    sub.add_parser("self-test")
    args = parser.parse_args()
    if args.command == "propose":
        propose(args)
    elif args.command == "summarize":
        summarize(args)
    else:
        self_test()


if __name__ == "__main__":
    main()
