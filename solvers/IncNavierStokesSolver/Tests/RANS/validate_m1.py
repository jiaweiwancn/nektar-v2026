#!/usr/bin/env python3
"""Validate the native Milestone 1 probe and record reproducible Linux evidence.

Uses only Python's standard library. All child processes share an affinity of
at most 16 logical CPUs, with one OpenMP/BLAS thread per process. MPI runs are
sequential with the serial run. This script does not build or install Nektar++.
"""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
from datetime import datetime, timezone


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate(report):
    require(report["strip_row_tests_passed"], "Strip assembly checks failed")
    require(abs(report["corrected_source_mean"]) < 1e-11,
            "The corrected strip mean is not balanced")
    require(abs(report["corrected_source_linear_moment"]) > 1e-3,
            "The naive CG wall-source counterexample was not detected")
    require(abs(report["resolved_UP"] - report["auxiliary_UP"]) > 1e-3,
            "Resolved and auxiliary velocities were conflated")
    records = report["records"]
    fine = next(r for r in records if r["degree"] == 14)
    coarse = records[0]
    for key in ("momentum", "k", "epsilon", "normal_momentum"):
        require(fine[key] < 2e-6, f"Fine residual {key} failed")
        require(fine[key] < coarse[key] * 1e-4,
                f"p convergence of {key} failed")
    h_coarse = next(r for r in records if r["degree"] == 4 and r["subdivisions"] == 1)
    h_fine = next(r for r in records if r["subdivisions"] == 8)
    for key in ("momentum", "k", "epsilon", "normal_momentum"):
        require(h_fine[key] < h_coarse[key] * 0.03,
                f"h convergence of {key} failed")
    for record in records:
        for key, value in record.items():
            require(math.isfinite(value), f"Nonfinite result in {key}")
        require(record["matching_constraint"] < 1e-10,
                "Matching constraints were not satisfied")
        require(record["continuity"] < 1e-13,
                "Flat-reference weak continuity was not satisfied")
        require(record["linear_be_conservation_error"] < 1e-9,
                "Coupled backward-Euler inventory balance failed")
        require(record["linear_be_stationary_defect"] < 1e-9,
                "Frozen linear numerical equilibrium was not preserved")
    for key, reference in (("control_missing_yap", "epsilon"),
                           ("control_zero_top_epsilon_flux", "epsilon"),
                           ("control_wrong_wall_traction", "momentum")):
        require(fine[key] > 100 * fine[reference],
                f"Negative control {key} was not detected")
    require(fine["control_pressure_gradient"] > 1e-3 and
            fine["control_crossflow_divergence"] > 1e-3,
            "Pressure or continuity controls failed")
    quadrature = [r for r in records if r["degree"] == 10]
    for key in ("momentum", "k", "epsilon", "normal_momentum"):
        require(abs(quadrature[0][key] - quadrature[1][key]) < 1e-8,
                f"Additional quadrature changed {key}")
    return fine


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--mpi-ranks", type=int, default=4)
    parser.add_argument("--mpiexec", default="mpirun")
    args = parser.parse_args()
    require(sys.platform.startswith("linux"), "Authoritative validation must run on Linux")
    require(1 <= args.mpi_ranks <= 16, "MPI allocation must be between 1 and 16")
    available = sorted(os.sched_getaffinity(0))
    allocation = available[-min(16, len(available)):]
    require(args.mpi_ranks <= len(allocation), "MPI ranks exceed the available CPU allocation")
    os.sched_setaffinity(0, allocation)
    environment = dict(os.environ)
    threads = {name: "1" for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
               "MKL_NUM_THREADS", "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")}
    environment.update(threads)
    args.output.mkdir(parents=True, exist_ok=True)
    root = Path(__file__).resolve().parents[4]
    probe = args.probe.resolve()
    source_files = (
        "library/SolverUtils/Turbulence/RANSWallStrip.h",
        "library/UnitTests/SolverUtils/RANSWallLayerProbe.cpp",
        "library/SolverUtils/CMakeLists.txt",
        "solvers/IncNavierStokesSolver/CMakeLists.txt",
        "solvers/IncNavierStokesSolver/Tests/RANS/CMakeLists.txt",
        "solvers/IncNavierStokesSolver/Tests/RANS/validate_m1.py",
    )
    commands = []

    def run(command, label):
        commands.append(command)
        print(f"Running {label}: {' '.join(command)}", flush=True)
        result = subprocess.run(command, env=environment, cwd=args.output,
                                text=True, capture_output=True, timeout=120)
        (args.output / f"{label}.log").write_text(result.stdout + result.stderr)
        require(result.returncode == 0, f"{label} failed; see {label}.log")
        print(result.stdout.strip(), flush=True)

    serial_path = args.output / "serial.json"
    run([str(probe), "--output", str(serial_path)], "serial")
    serial = json.loads(serial_path.read_text())
    fine = validate(serial)
    mpi_path = args.output / "mpi.json"
    run([args.mpiexec, "--bind-to", "none", "--map-by", "slot", "-np",
         str(args.mpi_ranks), str(probe), "--mpi", "--output", str(mpi_path)], "mpi")
    mpi = json.loads(mpi_path.read_text())
    validate(mpi)
    require(mpi["ranks"] == args.mpi_ranks, "MPI communicator size mismatch")
    if args.mpi_ranks > 1:
        require(mpi["records"][0]["sample_owner"] != 0,
                "N sample was not placed off the wall-owning rank")
    difference = 0.0
    for a, b in zip(serial["records"], mpi["records"]):
        for key in a:
            if key != "sample_owner":
                difference = max(difference, abs(a[key] - b[key]))
                require(abs(a[key] - b[key]) < 3e-10 * max(1.0, abs(a[key])),
                        f"Serial/MPI discrepancy in {key}")
    summary = {
        "schema": 1,
        "status": "passed",
        "milestone": "neutral_abl_rans_awf_m1",
        "scope": serial["scope"],
        "utc": datetime.now(timezone.utc).isoformat(),
        "source_commit": subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip(),
        "source_status": subprocess.check_output(
            ["git", "-C", str(root), "status", "--short"], text=True).splitlines(),
        "source_sha256": {p: sha256(root / p) for p in source_files},
        "probe_sha256": sha256(probe),
        "python": platform.python_version(),
        "system": platform.platform(),
        "cpu_affinity": allocation,
        "maximum_project_cpus": 16,
        "maximum_test_ranks": args.mpi_ranks,
        "thread_environment": threads,
        "commands": commands,
        "serial_mpi_maximum_difference": difference,
        "wall_source_mean": serial["corrected_source_mean"],
        "wall_source_linear_moment": serial["corrected_source_linear_moment"],
        "fine_column": fine,
        "reports": {p.name: sha256(p) for p in (serial_path, mpi_path)},
    }
    (args.output / "validation.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"status": "passed", "fine_column": fine,
                      "serial_mpi_maximum_difference": difference}, indent=2))


if __name__ == "__main__":
    main()
