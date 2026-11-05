#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Audit the original reduced merit gradient in recorded response coordinates.

Near-active windows use the saved cut normalization, not a physical positivity
theorem. This measured linearization is not a nonlinear stationarity certificate
and cannot establish global infeasibility outside its tested response space.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import struct
import sys

THREADS = {name: "1" for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "MKL_NUM_THREADS", "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS")}
os.environ.update(THREADS)
import numpy as np
import scipy
from scipy.optimize import nnls


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(root, specification, windows):
    name, iteration = specification.rsplit(":", 1)
    iteration = int(iteration)
    directory = root / name
    snapshot = (directory / "target/solver.log").read_bytes()
    packet = selected = None
    for line in snapshot.splitlines(keepends=True):
        if not line.startswith(b"EP-AWF Reduced ") or not line.endswith(b"\n"):
            continue
        record = json.loads(line[15:])
        if record["event"] == "model_packet":
            packet = record
        elif record["event"] == "direction" and record["iteration"] == iteration:
            selected = (packet, record)
    if selected is None:
        raise ValueError("requested native model has not been published")
    packet, direction = selected
    if direction["restoration"] or direction["residual_weight"] != 1.:
        raise ValueError("use an original-merit model without a restoration slack")
    path = directory / "target" / packet["path"]
    with path.open("rb") as stream:
        n, m = struct.unpack("<QQ", stream.read(16))
        values = np.fromfile(stream, dtype="<f8")
    if n != direction["basis_dimension"] or values.size != n*n+2*n+m*(n+1):
        raise ValueError("native model payload dimensions differ")
    hessian = values[:n*n].reshape(n, n)
    gradient = values[n*n:n*n+n]
    initial = values[n*n+n:n*n+2*n]
    cuts = values[n*n+2*n:].reshape(m, n+1)
    if not np.array_equal(gradient, direction["original_merit_gradient"]) or np.any(initial):
        raise ValueError("payload does not contain the original gradient at the base point")
    length = np.linalg.norm(gradient)
    if not length:
        raise ValueError("zero-gradient models need a separate residual/qualification audit")
    eigen = np.linalg.eigvalsh(.5*(hessian+hessian.T))
    uncertainty = np.linalg.norm(direction["two_scale_gradient_error_estimate"])
    results = []
    for window in windows:
        mask = cuts[:, -1] <= window
        rows = cuts[mask, :-1]
        ids = np.flatnonzero(mask)
        if rows.size:
            weights, _ = nnls(rows.T, -gradient/length, maxiter=10000)
            residual = gradient/length + rows.T @ weights
            dual = rows @ residual
            if np.min(dual) < -1e-8 or np.max(np.abs(weights*dual)) > 1e-8:
                raise ValueError("nonnegative least-squares audit misses its own first-order conditions")
            support = [{"cut": int(ids[i]), "multiplier": float(weights[i]*length),
                        "normalized_base_slack": float(cuts[ids[i], -1])}
                       for i in np.flatnonzero(weights)]
        else:
            residual = gradient/length
            support = []
        results.append(dict(coordinate_window=window, near_active_rows=int(mask.sum()),
            gradient_cone_defect=float(np.linalg.norm(residual)*length),
            empirical_two_size_gradient_error=float(uncertainty), dual_support=support))
    manifest = directory / "manifest.json"
    launch = json.loads(manifest.read_text()) if manifest.exists() else {}
    return dict(case=name, iteration=iteration, operator_alpha=direction["alpha"],
        time=direction["time"], model=str(path), model_sha256=sha(path),
        log_snapshot_sha256=hashlib.sha256(snapshot).hexdigest(), native_source=launch.get("source"),
        solver=launch.get("IncNavierStokesSolver", {}), tested_dimension=int(n),
        measured_response_rank=direction["response_rank"],
        maximum_column_correlation=direction.get("maximum_column_correlation"),
        original_merit_gradient=gradient.tolist(), gradient_norm=float(length),
        empirical_gradient_error=direction["two_scale_gradient_error_estimate"],
        saved_model_symmetric_eigen_range=[float(eigen[0]), float(eigen[-1])],
        saved_model_condition=float(eigen[-1]/eigen[0]) if eigen[0] > 0 else None,
        saved_model_radius_multiplier=direction.get("radius_multiplier"), windows=results)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--case", action="append", required=True,
                        help="Native case directory and numerical iteration, NAME:ITERATION")
    parser.add_argument("--windows", default="1e-10,1e-8,1e-6")
    args = parser.parse_args()
    windows = [float(x) for x in args.windows.split(",")]
    if any(not np.isfinite(x) or x <= 0 for x in windows):
        raise ValueError("invalid near-active coordinate window")
    result = dict(scope=__doc__, status="measured_first_order_audit_not_nonlinear_certificate",
        command=sys.argv, script_sha256=sha(Path(__file__)), thread_limits=THREADS,
        numpy_version=np.__version__, scipy_version=scipy.__version__,
        cases=[audit(args.run_root, case, windows) for case in args.case])
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
