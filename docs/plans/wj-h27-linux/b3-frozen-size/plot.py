#!/usr/bin/env python3
"""Make a compact plot of the two freshly evaluated forcing norms."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plot(path, output, ray):
    study = json.loads(path.read_text())
    target = np.array(study["targets"])
    def series(data, category):
        rows = sorted((row for row in data["candidates"] if row["label"] == category
                       and row["native_scaled_infinity"] > 0 and row["residuals"][0] is not None),
                      key=lambda row: row["native_scaled_infinity"])
        x = np.array([row["native_scaled_infinity"] for row in rows])
        y = np.array([row["residuals"][0] for row in rows])/target
        allowance = np.array([row["action_allowances"] if row["action_allowances"] is not None else
                              (np.nan, np.nan) for row in rows])/target
        return x, y, allowance

    constrained = series(study, "joint_model")
    archived = json.loads(ray.read_text())
    np.testing.assert_array_equal(target, archived["targets"])
    scaled = series(archived, "scaled_archived_R9_direction")
    full = series(archived, "archived_R9_direction")
    fig, axes = plt.subplots(1, 2, sharey=True, figsize=(10.5, 4.2), constrained_layout=True)
    for j, (ax, label) in enumerate(zip(axes, ("Physical block maximum", "Row-scaled L2"))):
        for (x, values, allowance), caption, style, color in (
                (constrained, "Size-constrained candidates", "o-", "#176a7c"),
                (scaled, "Scaled archived direction", "s:", "#8a4b36"),
                (full, "Full archived direction", "D", "#5e6270")):
            y, error = values[:, j], allowance[:, j]
            ax.plot(x, y, style, markersize=4, linewidth=1.25, color=color, label=caption)
            ax.vlines(x, np.maximum(.005, y-error), y+error, color=color, alpha=.4, linewidth=1.)
        ax.axhline(1., color="#313338", ls="--", lw=1.05, label="Original forcing target")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(label, loc="left", fontsize=11, fontweight="semibold")
        ax.set_xlabel(r"Correction size $\|D_x^{-1}d\|_\infty$")
        ax.grid(which="major", alpha=.18)
    axes[0].set_ylabel("Fresh residual / original target")
    axes[0].legend(frameon=False, fontsize=8, loc="best")
    fig.suptitle("Frozen B3 linear residual versus correction size", x=.018, ha="left", fontweight="semibold")
    fig.savefig(output, dpi=180)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--ray", type=Path, required=True)
    args = parser.parse_args()
    plot(args.results, args.out, args.ray)
