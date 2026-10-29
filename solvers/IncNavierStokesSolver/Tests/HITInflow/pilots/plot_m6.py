#!/usr/bin/env python3
"""Plot completed M6 pilot summaries; never infer milestone acceptance.

Requires matplotlib on the plotting host. All inputs are existing JSON outputs;
histories, solver fields and downstream experimental measurements are not loaded.
"""

import argparse
import datetime
import hashlib
import json
import math
from pathlib import Path


QUANTITIES = ("u", "v", "w", "K")
LABELS = (r"$\mathrm{var}(u)$", r"$\mathrm{var}(v)$", r"$\mathrm{var}(w)$", r"$K$")
COLORS = ("#0072B2", "#D55E00", "#009E73", "#CC79A7")


class Inputs:
    def __init__(self):
        self.files = {}

    def read(self, path):
        path = path.resolve()
        raw = path.read_bytes()
        receipt = {"sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw)}
        if str(path) in self.files and self.files[str(path)] != receipt:
            raise ValueError(f"input changed while plotting: {path}")
        self.files[str(path)] = receipt
        value = json.loads(raw)
        finite_tree(value)
        return value


def finite_tree(value):
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError("input JSON contains a nonfinite number")
    if isinstance(value, dict):
        for item in value.values():
            finite_tree(item)
    elif isinstance(value, list):
        for item in value:
            finite_tree(item)


def quantities(plane):
    values = plane["variance_m2_s2"] + [plane["K_m2_s2"]]
    if len(values) != 4 or any(not math.isfinite(v) or v < 0 for v in values):
        raise ValueError("expected three finite nonnegative variances and K")
    return values


def planes(stats):
    if stats["schema"] != "NEKTAR_HIT_PILOT_STATISTICS_V1":
        raise ValueError("unsupported statistics schema")
    result = {p["station_x_m"]: quantities(p) for p in stats["planes"]}
    if not result or len(result) != len(stats["planes"]):
        raise ValueError("missing or duplicate statistics stations")
    return result


def canonical(manifest):
    synthesis = manifest["synthesis"]
    sources = sorted(value for key, value in manifest["sha256"].items()
                     if Path(key).suffix == ".dat" and "spectrum" in Path(key).name.lower())
    if not sources:
        raise ValueError("canonical inlet spectrum hash is missing")
    return {key: synthesis[key] for key in (
        "periods_m", "counts", "seed_uint64", "cutoff_rad_m", "taper_start_rad_m")} | {
        "source_spectrum_sha256": sources,
        "spectrum_definition": {key: manifest["physics"][key] for key in (
            "mean_velocity_m_s", "normalization", "interpolation", "low_extension", "high_extension")}}


def same_window(left, right):
    keys = ("requested_start_s", "requested_end_s", "first_sample_s", "last_sample_s", "sample_dt_s")
    return left["samples"] == right["samples"] and all(
        math.isclose(left[key], right[key], rel_tol=2e-12, abs_tol=2e-14) for key in keys)


def load_case(inputs, root, name, flow=False):
    directory = root / name
    manifest_path = directory / "manifest.json"
    manifest = inputs.read(manifest_path)
    stats = inputs.read(directory / "statistics.json")
    run = inputs.read(directory / "run.receipt.json")
    if (run["solver"]["exit_code"] != 0 or run["solver"].get("native_iteration_limit_exceeded", False)
            or not run.get("finished_utc")):
        raise ValueError(f"{directory}: solver receipt is not a completed successful run")
    if stats["inputs"]["manifest_sha256"] != inputs.files[str(manifest_path.resolve())]["sha256"]:
        raise ValueError(f"{directory}: statistics/manifest hash mismatch")
    window = stats["window"]
    if not (math.isclose(window["requested_start_s"], manifest["time"]["record_start_s"], abs_tol=2e-14)
            and math.isclose(window["requested_end_s"], manifest["time"]["end_s"], abs_tol=2e-14)):
        raise ValueError(f"{directory}: statistics do not cover the declared full recording window")
    if not (math.isclose(window["first_sample_s"], window["requested_start_s"], abs_tol=2e-14)
            and math.isclose(window["last_sample_s"], window["requested_end_s"], abs_tol=2e-14)):
        raise ValueError(f"{directory}: sampled recording window is incomplete")
    result = {"name": name, "manifest": manifest, "stats": stats, "planes": planes(stats),
              "canonical": canonical(manifest), "revision": run["revision"]}
    if flow:
        result["flow"] = inputs.read(directory / (Path(manifest["case"]).stem + ".flow.json"))
    return result


def require_aligned(case, reference):
    if case["canonical"] != reference["canonical"]:
        raise ValueError(f"{case['name']} / {reference['name']}: canonical inflow differs")
    if not same_window(case["stats"]["window"], reference["stats"]["window"]):
        raise ValueError(f"{case['name']} / {reference['name']}: sampled recording windows differ")


def relative(case, reference):
    require_aligned(case, reference)
    shared = sorted(case["planes"].keys() & reference["planes"].keys())
    if not shared:
        raise ValueError(f"{case['name']} / {reference['name']}: no shared stations")
    rows = []
    for x in shared:
        if any(v <= 0 for v in reference["planes"][x]):
            raise ValueError(f"{reference['name']}: nonpositive reference variance/K at x={x}")
        rows.append([100 * (a - b) / b for a, b in zip(case["planes"][x], reference["planes"][x])])
    return shared, rows


def block_metrics(stats, x, count):
    if stats["blocks"]["count"] != count:
        raise ValueError(f"expected {count} blocks")
    matches = [p for p in stats["blocks"]["planes"] if p["station_x_m"] == x]
    if len(matches) != 1:
        raise ValueError(f"expected one block-energy station at x={x}")
    blocks = matches[0]["component_variance_m2_s2"] + [matches[0]["K_m2_s2"]]
    if len(blocks) != 4:
        raise ValueError("expected three component block energies and K")
    for block in blocks:
        low, high = block["conditional_95_percent_t_interval"]
        if (len(block["block_values"]) != count or not low <= block["mean_of_blocks"] <= high
                or not math.isfinite(low) or not math.isfinite(high)):
            raise ValueError("invalid block-energy mean/conditional interval")
    return blocks


def sampling_data(inputs, root):
    records = [inputs.read(root / f"record{label}.json") for label in ("01", "03", "07")]
    variants = [inputs.read(root / "record07-blocks4.json"), records[-1],
                inputs.read(root / "record07-blocks16.json")]
    reference, x = records[-1], 0.45
    ref = planes(reference)[x]
    if any(v <= 0 for v in ref):
        raise ValueError("record07 reference variance/K must be positive")
    source = {key: reference["inputs"][key] for key in ("history_sha256", "manifest_sha256")}
    changes, intervals, widths = [], [], []
    for stats, expected_width in zip(records, (0.1, 0.3, 0.7)):
        if any(stats["inputs"][key] != value for key, value in source.items()):
            raise ValueError("recording-width analyses must use the same completed history/manifest")
        window = stats["window"]
        width = window["last_sample_s"] - window["first_sample_s"]
        if (not math.isclose(width, expected_width, rel_tol=2e-12, abs_tol=2e-14)
                or not math.isclose(window["first_sample_s"], 0.1, abs_tol=2e-14)
                or not math.isclose(window["sample_dt_s"], reference["window"]["sample_dt_s"], abs_tol=2e-14)):
            raise ValueError("expected recording windows .1-.2, .1-.4 and .1-.8 s")
        changes.append([100 * (a - b) / b for a, b in zip(planes(stats)[x], ref)])
        intervals.append([{key: [100 * (v - base) / base for v in block[key]]
                            if key == "conditional_95_percent_t_interval" else 100 * (block[key] - base) / base
                           for key in ("mean_of_blocks", "conditional_95_percent_t_interval")}
                          for block, base in zip(block_metrics(stats, x, 8), ref)])
        widths.append(width)
    halfwidths = []
    for stats, count in zip(variants, (4, 8, 16)):
        if (any(stats["inputs"][key] != value for key, value in source.items())
                or not same_window(stats["window"], reference["window"])
                or planes(stats)[x] != ref):
            raise ValueError("block-count analyses must share the fixed .7 s history/window/statistics")
        halfwidths.append([50 * (b["conditional_95_percent_t_interval"][1] -
                                 b["conditional_95_percent_t_interval"][0]) / base
                           for b, base in zip(block_metrics(stats, x, count), ref)])
    return {"station_x_m": x, "widths_s": widths, "full_window_change_percent": changes,
            "conditional_8_block_energy_relative_to_record07_percent": intervals,
            "block_counts": [4, 8, 16], "conditional_ci_halfwidth_percent": halfwidths}


def save(figure, output, name):
    for extension in ("png", "svg"):
        figure.savefig(output / f"{name}.{extension}", dpi=180, bbox_inches="tight")


def plot_short(plt, output, cases, rows):
    figure, (ax, work) = plt.subplots(1, 2, figsize=(13.2, 5.8), gridspec_kw={"width_ratios": [1.7, 1]})
    for component, (label, color) in enumerate(zip(LABELS, COLORS)):
        ax.scatter([row["maximum_absolute_change_percent"][component] for row in rows],
                   [i + (component - 1.5) * .12 for i in range(len(rows))], label=label, color=color, s=30)
    ax.axvline(2, color="#777777", linestyle="--", label="2% sensitivity target")
    ax.set_yticks(range(len(rows)), [row["label"] for row in rows])
    ax.invert_yaxis()
    ax.set_xscale("symlog", linthresh=2)
    ax.set_xlabel(r"Max over shared $x$ of $|$case $-$ reference$|$/reference (%)")
    ax.set_title("Short-domain full-window sensitivity")
    ax.grid(axis="x", alpha=.2)
    ax.legend(fontsize=8, loc="lower right")
    molecular, svv = [], []
    for name in ("base", "p9", "h05"):
        flow = cases[name]["flow"]
        molecular.append(flow["instantaneous_resolved_strain_dissipation"])
        svv.append(flow["instantaneous_svv_work_xy"] + flow["instantaneous_svv_work_z"])
    work.bar([i - .19 for i in range(3)], molecular, .38, label=r"Molecular $2\nu S:S$", color=COLORS[0])
    work.bar([i + .19 for i in range(3)], svv, .38, label="SVV work: xy + z", color=COLORS[1])
    work.set_xticks(range(3), ["base", "p9", "h05"])
    work.set_ylabel(r"Volume-mean instantaneous work ($\mathrm{m^2\,s^{-3}}$)")
    work.set_title("Final-field instantaneous work")
    work.grid(axis="y", alpha=.2)
    work.legend(fontsize=8)
    figure.text(.02, .015, "Descriptive sensitivity; the 2% guide does not close a milestone. DG h/p refinement couples resolution and SVV model effects.", fontsize=9)
    figure.tight_layout(rect=(0, .055, 1, 1))
    save(figure, output, "m6-short-sensitivity")
    plt.close(figure)
    return {"rows": rows, "work_case_order": ["base", "p9", "h05"],
            "molecular_m2_s3": molecular, "svv_xy_plus_z_m2_s3": svv}


def plot_sampling(plt, output, data):
    figure, (ax, blocks) = plt.subplots(1, 2, figsize=(13.2, 5.6))
    for component, (label, color) in enumerate(zip(LABELS, COLORS)):
        centers = [width + (component - 1.5) * .008 for width in data["widths_s"]]
        ax.plot(centers, [row[component] for row in data["full_window_change_percent"]],
                marker="o", color=color, linewidth=1, label=label)
        means = [row[component]["mean_of_blocks"] for row in data["conditional_8_block_energy_relative_to_record07_percent"]]
        ci = [row[component]["conditional_95_percent_t_interval"] for row in data["conditional_8_block_energy_relative_to_record07_percent"]]
        ax.errorbar([x + .006 for x in centers], means,
                    yerr=[[mean - bounds[0] for mean, bounds in zip(means, ci)],
                          [bounds[1] - mean for mean, bounds in zip(means, ci)]],
                    fmt="D", markerfacecolor="white", markersize=4, color=color, capsize=3, linewidth=.9)
        blocks.plot(data["block_counts"], [row[component] for row in data["conditional_ci_halfwidth_percent"]],
                    marker="o", color=color, label=label)
    ax.axhline(0, color="#777777", linewidth=.8)
    ax.set_xticks([.1, .3, .7])
    ax.set_xlabel("Recording width (s), from t = 0.1 s")
    ax.set_ylabel("Change / fixed 0.7 s full-window statistic (%)")
    ax.set_title(r"$x=0.45$ m: full-window values and 8-block energy CIs")
    ax.legend(fontsize=9)
    blocks.set_xticks([4, 8, 16])
    blocks.set_xlabel("Number of contiguous blocks (fixed 0.7 s record)")
    blocks.set_ylabel("Conditional 95% CI halfwidth / full-window statistic (%)")
    blocks.set_title("Block-energy mean CI: block-count sensitivity")
    blocks.legend(fontsize=9)
    for axis in (ax, blocks):
        axis.grid(alpha=.2)
    figure.text(.02, .04, "Filled circles: full-window variance/K. Open diamonds: mean block energy and conditional 95% Student-t CI; the estimators differ.", fontsize=9)
    figure.text(.02, .01, "CIs assume approximate stationarity and independent blocks; nonoverlap is not proof. Reference uncertainty is not propagated.", fontsize=9)
    figure.tight_layout(rect=(0, .095, 1, 1))
    save(figure, output, "m6-recording-sensitivity")
    plt.close(figure)


def plot_long(plt, output, stations, changes):
    figure, ax = plt.subplots(figsize=(9.4, 4.8))
    for component, (label, color) in enumerate(zip(LABELS, COLORS)):
        ax.plot(stations, [row[component] for row in changes], marker="o", color=color, label=label)
    ax.axhline(0, color="#777777", linewidth=.8)
    ax.axhline(2, color="#777777", linestyle="--", label=r"$\pm 2\%$ sensitivity guides")
    ax.axhline(-2, color="#777777", linestyle="--")
    ax.set_xlabel("Downstream station x (m)")
    ax.set_ylabel("100 × (p9 − p10) / p10 (%)")
    ax.set_title("Completed long-domain records: full-window variance/K sensitivity")
    ax.grid(alpha=.2)
    ax.legend(ncol=5, fontsize=9)
    figure.text(.02, .01, "Matched canonical inlet and sampled windows. DG p refinement couples resolution and SVV effects. No milestone decision is inferred.", fontsize=9)
    figure.tight_layout(rect=(0, .065, 1, 1))
    save(figure, output, "m6-long-sensitivity")
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--short-root", required=True, type=Path)
    parser.add_argument("--sampling-root", required=True, type=Path)
    parser.add_argument("--long-root", type=Path)
    parser.add_argument("--skip-long", action="store_true", help="Explicitly omit unfinished long records")
    parser.add_argument("--output", required=True, type=Path, help="Fresh artifact directory")
    args = parser.parse_args()
    try:
        output = args.output.resolve()
        if output.exists():
            raise ValueError("--output must name a fresh directory")
        if not args.skip_long and args.long_root is None:
            raise ValueError("supply --long-root with completed records, or explicitly use --skip-long")
        inputs = Inputs()
        pairs = [("base", "h05"), ("p9", "h05")] + [(name, "base") for name in (
            "dt-half", "z64", "buffer09", "buffer12", "svv-half", "svv-double")]
        cases = {name: load_case(inputs, args.short_root, name, flow=name in ("base", "p9", "h05"))
                 for name in sorted({name for pair in pairs for name in pair})}
        short_rows = []
        for name, reference in pairs:
            stations, changes = relative(cases[name], cases[reference])
            short_rows.append({"label": f"{name}  vs  {reference}", "case": name, "reference": reference,
                               "shared_stations_x_m": stations,
                               "maximum_absolute_change_percent": [max(abs(row[c]) for row in changes) for c in range(4)]})
        sampling = sampling_data(inputs, args.sampling_root)
        long_cases, long_data = {}, None
        if not args.skip_long:
            long_cases = {name: load_case(inputs, args.long_root, name) for name in ("long-p9", "long-p10")}
            if long_cases["long-p9"]["revision"] != long_cases["long-p10"]["revision"]:
                raise ValueError("long-p9 and long-p10 solver source revisions differ")
            stations, changes = relative(long_cases["long-p9"], long_cases["long-p10"])
            long_data = {"shared_stations_x_m": stations, "p9_minus_p10_percent": changes}
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        output.mkdir(parents=True, exist_ok=False)
        plt.rcParams.update({"font.size": 10, "svg.fonttype": "none"})
        short_data = plot_short(plt, output, cases, short_rows)
        plot_sampling(plt, output, sampling)
        if long_data is not None:
            plot_long(plt, output, long_data["shared_stations_x_m"], long_data["p9_minus_p10_percent"])
        names = ["m6-short-sensitivity", "m6-recording-sensitivity"] + (["m6-long-sensitivity"] if long_data is not None else [])
        provenance = {
            "schema": "NEKTAR_HIT_M6_FIGURES_V1", "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "plotting_script": {"path": str(Path(__file__).resolve()), "sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
            "matplotlib_version": matplotlib.__version__, "backend": "Agg", "quantities_in_column_order": list(QUANTITIES),
            "arguments": {key: str(value.resolve()) if isinstance(value, Path) else value for key, value in vars(args).items()},
            "input_files": inputs.files,
            "case_source_revisions": {name: case["revision"] for name, case in (cases | long_cases).items()},
            "figures": {name: [f"{name}.png", f"{name}.svg"] for name in names},
            "plotted_data": {"short": short_data, "sampling": sampling, "long": long_data},
            "long_omission": "Explicit --skip-long; no long results were read or fabricated" if args.skip_long else None,
            "definitions": {"short_relative_change": "100*(case-reference)/reference; absolute maximum over shared stations, separately for var(u),var(v),var(w),K",
                "sampling_reference": "record07 full-window variance/K; conditional intervals describe mean block energy, not the full-window variance estimator",
                "work": "Final-field instantaneous volume means, m^2/s^3; molecular=2*nu*S:S, total SVV work=xy+z",
                "interpretation": "Descriptive sensitivity only. The 2% guides do not imply acceptance. Conditional Student-t intervals assume approximately stationary independent blocks; nonoverlap is not proof of independence. No downstream experimental data were loaded."}}
        (output / "provenance.json").write_text(json.dumps(provenance, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        print(f"Wrote {len(names)} PNG/SVG figure pairs and provenance.json to {output}")
    except (OSError, ValueError, KeyError, TypeError, ImportError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
