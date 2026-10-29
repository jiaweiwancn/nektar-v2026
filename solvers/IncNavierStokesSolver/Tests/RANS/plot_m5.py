#!/usr/bin/env python3
"""Plot recorded corrected-AWF refinement, source balance and switch receipts."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS"):
    os.environ[key] = "1"
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("validation", type=Path); parser.add_argument("output", type=Path)
    args = parser.parse_args(); report = json.loads(args.validation.read_text())
    args.output.mkdir(parents=True, exist_ok=True)

    def save(fig, name):
        for extension in ("png", "pdf"):
            fig.savefig(args.output/(name+"."+extension), dpi=180)
        plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), constrained_layout=True)
    degrees = [4, 6, 8]
    for field in ("u", "v", "k", "eps_tilde", "p"):
        values = [report["cases"][f"corrected_p{p}"]["field_audit"]["fields"][field]["l2_error"] for p in degrees]
        axes[0].semilogy(degrees, values, "o-", label=field)
    residuals = report["checks"]["initial_assembled_stage_residuals"]
    for field in ("u", "v", "k", "eps_tilde"):
        axes[1].semilogy(degrees, [residuals[f"corrected_p{p}"][field] for p in degrees], "o-", label=field)
    axes[0].set(ylabel="Physical L2 error at 0.1 s", title="Corrected reference: field refinement")
    axes[1].set(ylabel="Free global coefficient residual norm", title="Initial assembled spatial defect")
    for axis in axes:
        axis.set(xlabel="Velocity/scalar polynomial degree", xticks=degrees)
        axis.grid(True, which="both", alpha=.2); axis.legend(fontsize=8)
    save(fig, "corrected-refinement")

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), constrained_layout=True)
    for name, label, marker in (("original", "original AWF", "o"), ("corrected", "fixed chi=1", "s"), ("switched_reference", "switched reference", "x")):
        rows = report["cases"][name]["wall_stages"]; times = [row["time"] for row in rows]
        axes[0].plot(times, [row["mean_kP"] for row in rows], marker+"-", label=label)
        axes[1].plot(times, [row["mean_production"]-row["mean_dissipation"] for row in rows], marker+"-", label=label)
    axes[0].axhline(3.190675081624288*.4**2, color="black", linestyle="--", linewidth=1, label="analytical reference")
    axes[0].set(ylabel="Strip representative k (m²/s²)", title="Wall-strip evolution")
    axes[1].axhline(0, color="black", linewidth=1)
    axes[1].set(ylabel="Mean production minus dissipation (m²/s³)", title="The correction removes the source mismatch")
    for axis in axes:
        axis.set(xlabel="Time (s)"); axis.grid(True, alpha=.2); axis.legend(fontsize=8)
    save(fig, "corrected-source-balance")

    # Exact raw receipts remain in the Linux case directory named by the
    # validation file. Its manifest records their SHA256s.
    directory = args.validation.parent/"forcing_cancellation"
    rows = [row for path in directory.glob("forcing_cancellation.awf-points-r*.csv") for row in csv.DictReader(path.open())]
    rows = sorted([row for row in rows if int(row["phase"]) == 0 and float(row["time"]) == .001], key=lambda row: int(row["point"]))
    if not rows:
        raise RuntimeError("The independent-forcing figure requires the hashed raw point receipts")
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), constrained_layout=True)
    point = [int(row["point"]) for row in rows]
    for key, label in (("r_Pressure", "pressure minus body"), ("r_Convection", "convection")):
        axes[0].semilogy(point, [max(1e-12, abs(float(row[key]))) for row in rows], "o-", label=label)
    summed = [2*math.hypot(float(row["forcing_x"]), float(row["forcing_y"]))/(float(row["kP"])/3.190675081624288+(1.5e-5/2)**2) for row in rows]
    axes[0].semilogy(point, [max(1e-12, value) for value in summed], "s-", label="summed forcing indicator")
    axes[0].set(ylabel="Dimensionless forcing indicator", title="Cancelling inputs remain independently sensed")
    axes[1].plot(point, [float(row["chi"]) for row in rows], "o-")
    axes[1].set(ylabel="Effective equilibrium weight", ylim=(-.05, 1.05), title="Finite-departure original-form recovery")
    for axis in axes:
        axis.set(xlabel="Wall quadrature column"); axis.grid(True, which="both", alpha=.2)
    axes[0].legend(fontsize=8)
    save(fig, "independent-forcing-switch")


if __name__ == "__main__":
    main()
