#!/usr/bin/env python3
"""Plot completed long-case statistics and native late-checkpoint diagnostics.

Requires matplotlib. Reads immutable JSON/CSV evidence only; no experimental
downstream measurements or solver fields are loaded by this plotting tool.
"""

import argparse
import csv
import datetime
import hashlib
import io
import json
import math
from pathlib import Path

import plot_m6


def profiles(inputs, root, name, case):
    receipt = inputs.read(root / name / "receipt.json")
    if receipt["status"] != "completed" or Path(receipt["case"]).resolve() != (
            Path(case["stats"]["inputs"]["manifest"]).parent.resolve()):
        raise ValueError("profile receipt is incomplete or refers to another case")
    snapshots = []
    for profile in receipt["profiles"]:
        candidates = [Path(p) for p in profile["output_files"] if p.endswith(".slabs.csv")]
        if len(candidates) != 1:
            raise ValueError("expected one native slab CSV per checkpoint")
        path = candidates[0].resolve()
        raw = path.read_bytes()
        record = {"sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw)}
        if profile["output_files"][str(path)] != record:
            raise ValueError("native slab CSV/receipt digest mismatch")
        inputs.files[str(path)] = record
        rows = [{key: float(value) for key, value in row.items()}
                for row in csv.DictReader(io.StringIO(raw.decode("utf-8")))]
        if len(rows) != case["manifest"]["mesh"]["nx"] or any(
                not math.isfinite(value) for row in rows for value in row.values()):
            raise ValueError("incomplete/nonfinite native profiles")
        volume = math.fsum(row["volume"] for row in rows)
        flow = profile["flow"]
        if not math.isclose(volume, flow["volume"], rel_tol=1e-10):
            raise ValueError("slab/global diagnostic volumes differ")
        for column, global_key in (("molecular_strain_dissipation", "instantaneous_resolved_strain_dissipation"),
                                   ("svv_xy", "instantaneous_svv_work_xy"),
                                   ("svv_z", "instantaneous_svv_work_z")):
            average = math.fsum(row["volume"] * row[column] for row in rows) / volume
            if not math.isclose(average, flow[global_key], rel_tol=1e-9, abs_tol=1e-14):
                raise ValueError("slab/global diagnostic work differs")
        snapshots.append({"time_s": profile["time_s"], "rows": rows,
                          "field_build_revision": profile["field_build_revision"], "flow": flow})
    snapshots.sort(key=lambda item: item["time_s"])
    if len(snapshots) != 3 or any(not math.isclose(s["time_s"], t, rel_tol=2e-12, abs_tol=2e-14)
                                   for s, t in zip(snapshots, (.85, .90, .95))):
        raise ValueError("expected completed late checkpoints .85, .90 and .95 s")
    x = [row["xmid"] for row in snapshots[0]["rows"]]
    if any([row["xmid"] for row in snapshot["rows"]] != x for snapshot in snapshots):
        raise ValueError("checkpoint slab grids differ")
    result = {"x_m": x, "times_s": [s["time_s"] for s in snapshots], "columns": {},
              "field_build_revisions": sorted({s["field_build_revision"] for s in snapshots})}
    for key in ("molecular_strain_dissipation", "svv_xy", "svv_z", "divergence_rms"):
        columns = [[row[key] for row in snapshot["rows"]] for snapshot in snapshots]
        result["columns"][key] = {"mean": [math.fsum(values) / 3 for values in zip(*columns)],
                                  "min": [min(values) for values in zip(*columns)],
                                  "max": [max(values) for values in zip(*columns)]}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--long-root", type=Path, required=True)
    parser.add_argument("--statistics-root", type=Path, required=True)
    parser.add_argument("--profiles-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.output.exists():
            raise ValueError("output must be a fresh directory")
        inputs = plot_m6.Inputs()
        cases = {name: plot_m6.load_case(inputs, args.long_root, name)
                 for name in ("long-p9", "long-p10")}
        plot_m6.require_aligned(cases["long-p9"], cases["long-p10"])
        stats, native = {}, {}
        for name, case in cases.items():
            stats[name] = inputs.read(args.statistics_root / (name + ".json"))
            if (any(stats[name]["inputs"][key] != case["stats"]["inputs"][key]
                    for key in ("history_sha256", "manifest_sha256")) or
                    not plot_m6.same_window(stats[name]["window"], case["stats"]["window"])):
                raise ValueError("spectral statistics and completed case inputs/window differ")
            native[name] = profiles(inputs, args.profiles_root, name, case)
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.rcParams.update({"font.size": 10, "svg.fonttype": "none"})
        figure, axes = plt.subplots(2, 2, figsize=(13.2, 8.8))
        energy, spectrum, work, divergence = axes.flat
        for name, color in (("long-p9", "#0072B2"), ("long-p10", "#D55E00")):
            planes = stats[name]["planes"]
            energy.plot([p["station_x_m"] for p in planes], [p["K_m2_s2"] for p in planes],
                        "o-", color=color, label=name)
            profile = native[name]
            values = profile["columns"]["divergence_rms"]
            scale = .024 / .222
            divergence.plot(profile["x_m"], [scale * v for v in values["mean"]], color=color, label=name)
            divergence.fill_between(profile["x_m"], [scale * v for v in values["min"]],
                                    [scale * v for v in values["max"]], color=color, alpha=.12)
        energy.set(title="Temporal-central TKE, t = 0.8–0.95 s", ylabel=r"$K$ ($\mathrm{m^2\,s^{-2}}$)")
        energy.legend(fontsize=8)
        spec = stats["long-p10"]["z_spectra"]
        active = cases["long-p10"]["manifest"]["mesh"]["hom_modes_z"] // 2
        spectra_data = []
        for x, color in ((0., "#0072B2"), (3.612896, "#009E73"), (8.322564, "#D55E00")):
            plane = next(p for p in spec["planes"] if p["station_x_m"] == x)
            values = [math.fsum(v) / 2 for v in zip(*plane["component_variance_per_mode_m2_s2"])]
            spectrum.semilogy(spec["kz_rad_m"][:active], values[:active], "o-", markersize=3,
                              color=color, label=f"x={x:g} m")
            spectra_data.append({"x_m": x, "K_per_mode_m2_s2": values})
        spectrum.axvline(100, color="#777777", linestyle="--", label="Input radial cutoff")
        spectrum.set(title="Degree-10 one-sided z mode energy", xlabel=r"$k_z$ (rad/m)",
                     ylabel=r"$\frac{1}{2}\sum_i P_{ii}$ per mode ($\mathrm{m^2\,s^{-2}}$)", ylim=(1e-9, 1e-2))
        spectrum.legend(fontsize=8)
        for key, label, color in (("molecular_strain_dissipation", "Molecular", "#0072B2"),
                                  ("svv_xy", "SVV xy", "#D55E00"), ("svv_z", "SVV z", "#009E73")):
            for name, style in (("long-p9", "--"), ("long-p10", "-")):
                profile = native[name]
                values = profile["columns"][key]
                work.plot(profile["x_m"], values["mean"], style, color=color,
                          label=f"{name}: {label}", linewidth=1.2)
                if name == "long-p10":
                    work.fill_between(profile["x_m"], values["min"], values["max"], color=color, alpha=.10)
        work.set(title="Native instantaneous work: 3-snapshot mean", ylabel=r"Work ($\mathrm{m^2\,s^{-3}}$)")
        work.legend(fontsize=7, ncol=2)
        divergence.set(title="Native divergence: 3-snapshot mean and range", ylabel=r"$0.024\,D_{\rm rms}/0.222$")
        divergence.legend(fontsize=8)
        for ax in (energy, work, divergence):
            ax.set_xlabel("x (m)")
            ax.axvline(8.322564, color="#777777", linestyle=":", linewidth=.8)
        for ax in axes.flat:
            ax.grid(alpha=.18)
        figure.text(.015, .037, "Native profiles: t=.85, .90, .95 s; shading is snapshot range, not a confidence interval. Work excludes transport/pressure/time terms.", fontsize=9)
        figure.text(.015, .014, "z spectra retain mode0; null Fourier Nyquist bin omitted from display. Values below 1e-9 are clipped by the display axis. No experimental data loaded.", fontsize=9)
        figure.tight_layout(rect=(0, .065, 1, 1))
        args.output.mkdir(parents=True, exist_ok=False)
        plot_m6.save(figure, args.output, "m6-long-diagnostics")
        provenance = {"schema": "NEKTAR_HIT_M6_LONG_DIAGNOSTIC_FIGURE_V1",
                      "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                      "scripts": {str(p.resolve()): hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in (Path(__file__), Path(plot_m6.__file__))},
                      "matplotlib_version": matplotlib.__version__, "backend": "Agg", "input_files": inputs.files,
                      "native_profiles_plotted": native, "degree10_z_mode_energy": spectra_data,
                      "definitions": "Temporal-central K is separate from instantaneous native work. Shaded ranges of three snapshots are descriptive, not CIs or closed energy balances. z mode energies are not shell spectra. Null Nyquist omitted only from display."}
        (args.output / "provenance.json").write_text(json.dumps(provenance, indent=2, allow_nan=False) + "\n")
        print(f"Wrote PNG/SVG and provenance to {args.output}")
    except (OSError, ValueError, KeyError, TypeError, ImportError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
