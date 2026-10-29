#!/usr/bin/env python3
"""Render scope-qualified M7 comparisons from four completed analysis JSON files.

Run on Linux with Matplotlib. This tool reads immutable analysis/source evidence,
never histories or solver fields. It neither fits parameters nor infers physical
validation from implementation-test success. Each family has two common-random-
input streams; the four curves are not four independent ensemble replicates.
"""

import argparse
import csv
import datetime
import hashlib
import io
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
FAMILIES = ("cbc1971", "von-karman")
FAMILY_LABELS = {"cbc1971": "CBC input", "von-karman": "von Karman input"}
PALETTES = {"cbc1971": ("#0072B2", "#56B4E9"),
            "von-karman": ("#D55E00", "#E69F00")}
QUANTITIES = (r"$\mathrm{var}(u)$", r"$\mathrm{var}(v)$",
              r"$\mathrm{var}(w)$", r"$K$")
SCOPE = "Limited inlet band with SVV; scope-qualified descriptive comparison. No full-band fidelity or 2% accuracy verdict."


def require(condition, message):
    if not condition:
        raise ValueError(message)


def finite_tree(value):
    if isinstance(value, float):
        require(math.isfinite(value), "nonfinite number in JSON evidence")
    elif isinstance(value, dict):
        for item in value.values():
            finite_tree(item)
    elif isinstance(value, list):
        for item in value:
            finite_tree(item)


def reject_constant(value):
    raise ValueError(f"nonstandard JSON constant: {value}")


class Inputs:
    def __init__(self):
        self.files = {}

    def raw(self, path, expected=None):
        path = Path(path).resolve()
        content = path.read_bytes()
        record = {"sha256": hashlib.sha256(content).hexdigest(), "bytes": len(content)}
        if expected is not None:
            expected_hash = expected["sha256"] if isinstance(expected, dict) else expected
            require(record["sha256"] == expected_hash, f"evidence hash mismatch: {path}")
        require(str(path) not in self.files or self.files[str(path)] == record,
                f"input changed while plotting: {path}")
        self.files[str(path)] = record
        return content

    def read(self, path, expected=None):
        value = json.loads(self.raw(path, expected), parse_constant=reject_constant)
        finite_tree(value)
        return value

    def verify_unchanged(self):
        for path, expected in list(self.files.items()):
            self.raw(path, expected)


def number(value):
    result = float(value)
    require(math.isfinite(result), "expected finite numeric data")
    return result


def index_stations(rows, description):
    result = {}
    for row in rows:
        coordinate = number(row["station_x_m"])
        require(coordinate not in result, f"duplicate {description} station {coordinate}")
        result[coordinate] = row
    require(result, f"missing {description} stations")
    return result


def station_at(rows, coordinate):
    matches = [row for x, row in rows.items()
               if math.isclose(x, coordinate, rel_tol=1e-10, abs_tol=1e-8)]
    require(len(matches) == 1, f"expected one station at x={coordinate:g} m")
    return matches[0]


def table_targets(source):
    experiment = source["experiment"]
    speed = number(experiment["measurement_speed_m_s"])
    factor = speed * number(experiment["grid_mesh_m"]) / number(experiment["approach_speed_m_s"])
    theta0 = number(experiment["reference_theta"])
    result = []
    for values in source["table4"]["rows"]:
        row = dict(zip(source["table4"]["columns"], values))
        theta = number(row["theta"])
        rms = number(row["u_rms_cm_s"]) * .01
        result.append({"theta": theta, "x_offset_m": factor * (theta - theta0),
                       "u_rms_m_s": rms, "K_isotropic_inference_m2_s2": 1.5 * rms * rms,
                       "epsilon_m2_s3": number(row["epsilon_cm2_s3"]) * 1e-4})
    require([row["theta"] for row in result] == [42., 98., 171.],
            "expected the three principal CBC Table-4 stations")
    return result


def load_case(inputs, path, source_path, correlation_path, source, table5):
    analysis = inputs.read(path)
    require(analysis["schema"] == "NEKTAR_HIT_M7_ANALYSIS_V1", "unsupported M7 analysis schema")
    family, name = analysis["case_family"], analysis["case_id"]
    require(family in FAMILIES and isinstance(name, str) and name, "invalid case identity")
    manifest_path = Path(analysis["manifest_path"]).resolve()
    manifest = inputs.read(manifest_path, analysis["sha256"][str(manifest_path)])
    require(manifest["schema"] == "NEKTAR_HIT_PILOT_V1", "unsupported case manifest")
    require(Path(analysis["case_directory"]).resolve() == manifest_path.parent,
            f"{name}: case-directory/manifest mismatch")
    require(analysis["physics"] == manifest["physics"] and analysis["synthesis"] == manifest["synthesis"],
            f"{name}: analysis physics/synthesis differs from manifest")
    require(str(analysis["manifest_seed"]) == str(manifest["synthesis"]["seed_uint64"]),
            f"{name}: analysis seed mismatch")
    for source_file in (source_path, correlation_path):
        candidates = [digest for key, digest in analysis["sha256"].items()
                      if Path(key).name == source_file.name]
        require(len(candidates) == 1 and candidates[0] == inputs.files[str(source_file)]["sha256"],
                f"{name}: plotted experimental source differs from analyzed source")
    require(analysis["experimental_source"] == source["source"], f"{name}: experimental metadata mismatch")
    selected = next(grid for grid in table5["grids"] if grid["id"] == "grid_5_08_cm")
    require(analysis["published_Table5_selected_grid"] == selected, f"{name}: wrong Table-5 grid")
    window = analysis["window"]
    require(window["requested_end_s"] > window["requested_start_s"] and window["samples"] >= 32,
            f"{name}: incomplete recording window")
    require(all(math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-10) for a, b in zip(
        window["first_last_samples_s"], (window["requested_start_s"], window["requested_end_s"]))),
        f"{name}: recording does not include the complete requested window")
    require(math.isclose(window["prescribed_convection_speed_m_s"],
                         number(source["experiment"]["measurement_speed_m_s"]), abs_tol=1e-12),
            f"{name}: convection speed differs from experimental convention")
    planes = index_stations(analysis["statistics"]["planes"], "moments")
    spectra = index_stations(analysis["longitudinal_spectra"], "spectral")
    screen = index_stations(analysis["statistics"]["qualification"], "qualification")
    require(planes.keys() == spectra.keys() == screen.keys(), f"{name}: inconsistent station sets")
    for plane in planes.values():
        variances = plane["variance_m2_s2"]
        require(len(variances) == 3 and all(number(v) >= 0 for v in variances), "invalid component variances")
        require(math.isclose(.5 * math.fsum(variances), plane["K_m2_s2"], rel_tol=1e-11, abs_tol=1e-14),
                f"{name}: K does not equal one half of total temporal variance")
    references = {number(pair["upstream_x_m"]) for pair in analysis["mean_convection_delayed_pairs"]}
    require(len(references) == 1, f"{name}: delayed correlations need one reference plane")
    reference = references.pop()
    require(station_at(planes, reference)["K_m2_s2"] > 0,
            f"{name}: normalized decay requires positive reference energy")
    return {"id": name, "family": family, "seed": str(analysis["manifest_seed"]),
            "analysis_path": str(Path(path).resolve()), "analysis": analysis, "manifest": manifest,
            "manifest_path": manifest_path, "directory": manifest_path.parent,
            "planes": planes, "spectra": spectra, "qualification": screen, "reference_x_m": reference}


def validate_suite(cases):
    require(len(cases) == 4 and len({case["id"] for case in cases}) == 4,
            "supply four distinct cases: CBC/von-Karman times two streams")
    seeds = {}
    for family in FAMILIES:
        group = [case for case in cases if case["family"] == family]
        require(len(group) == 2 and len({case["seed"] for case in group}) == 2,
                f"{family}: expected two distinct streams")
        seeds[family] = {case["seed"] for case in group}
    require(seeds[FAMILIES[0]] == seeds[FAMILIES[1]], "CBC/VK must share the same two declared stream identities")
    reference = cases[0]
    for case in cases[1:]:
        require(case["planes"].keys() == reference["planes"].keys(), "four-case station grids differ")
        for key in ("requested_start_s", "requested_end_s", "sample_dt_s", "samples"):
            require(case["analysis"]["window"][key] == reference["analysis"]["window"][key],
                    f"{case['id']}: recording window/sampling differs")
        for key in ("cutoff_rad_m", "taper_start_rad_m", "periods_m", "counts"):
            require(case["analysis"]["synthesis"][key] == reference["analysis"]["synthesis"][key],
                    f"{case['id']}: synthesis geometry/band differs")
        for key in ("mesh", "svv"):
            require(case["manifest"][key] == reference["manifest"][key],
                    f"{case['id']}: numerical mesh/SVV configuration differs")
        require(case["manifest"]["time"]["dt_s"] == reference["manifest"]["time"]["dt_s"],
                f"{case['id']}: integration timestep differs")
        require(case["reference_x_m"] == reference["reference_x_m"], "reference planes differ")
    ordered = sorted(cases, key=lambda case: (FAMILIES.index(case["family"]), int(case["seed"])))
    for case in ordered:
        index = sorted(seeds[case["family"]], key=int).index(case["seed"])
        case["color"] = PALETTES[case["family"]][index]
        case["style"] = "-" if index == 0 else "--"
        case["marker"] = "o" if case["family"] == "cbc1971" else "s"
        case["label"] = f"{FAMILY_LABELS[case['family']]}, seed {case['seed']}"
    return ordered


def line_style(case):
    return {"color": case["color"], "linestyle": case["style"], "linewidth": 1.35,
            "marker": case["marker"], "markersize": 3.5, "label": case["label"]}


def finish(plt, figure, axes, output, name, footer):
    for axis in axes:
        axis.grid(alpha=.18)
    figure.text(.014, .034, SCOPE, fontsize=9)
    figure.text(.014, .011, footer, fontsize=8.5)
    figure.tight_layout(rect=(0, .075, 1, .955))
    for suffix in ("png", "svg"):
        figure.savefig(output / f"{name}.{suffix}", dpi=180, bbox_inches="tight")
    plt.close(figure)


def plot_decay(plt, output, cases, targets):
    figure, axes = plt.subplots(2, 2, figsize=(13.2, 8.5))
    rms, energy, rms_norm, energy_norm = axes.flat
    plotted = {}
    for case in cases:
        x = sorted(case["planes"])
        u = [math.sqrt(case["planes"][station]["variance_m2_s2"][0]) for station in x]
        k = [case["planes"][station]["K_m2_s2"] for station in x]
        ref = station_at(case["planes"], case["reference_x_m"])
        u0, k0 = math.sqrt(ref["variance_m2_s2"][0]), ref["K_m2_s2"]
        require(u0 > 0 and k0 > 0, "positive reference values required for normalized decay")
        normalized_u, normalized_k = [value / u0 for value in u], [value / k0 for value in k]
        for axis, values in ((rms, u), (energy, k), (rms_norm, normalized_u), (energy_norm, normalized_k)):
            axis.plot(x, values, **line_style(case))
        plotted[case["id"]] = {"x_m": x, "u_rms_m_s": u, "K_m2_s2": k,
                               "u_rms_over_own_reference": normalized_u, "K_over_own_reference": normalized_k}
    reference_x = cases[0]["reference_x_m"]
    xp = [reference_x + row["x_offset_m"] for row in targets]
    up, kp = [row["u_rms_m_s"] for row in targets], [row["K_isotropic_inference_m2_s2"] for row in targets]
    for axis, values, label in (
            (rms, up, "CBC Table 4: streamwise rms"),
            (energy, kp, "CBC Table 4: K inferred by isotropy"),
            (rms_norm, [v / up[0] for v in up], "CBC: rms / theta-42 rms"),
            (energy_norm, [v / kp[0] for v in kp], "CBC: K / theta-42 K")):
        axis.plot(xp, values, "D:", color="#222222", markersize=5, linewidth=1, label=label)
        axis.set_xlabel("Downstream x (m)")
        axis.legend(fontsize=7.5, loc="best")
    rms.set(title="Absolute streamwise temporal rms", ylabel=r"$u'_{1,\mathrm{rms}}$ ($\mathrm{m\,s^{-1}}$)")
    energy.set(title="Absolute temporal-central TKE", ylabel=r"$K$ ($\mathrm{m^2\,s^{-2}}$)")
    rms_norm.set(title="Rms decay about each curve's own reference", ylabel=r"$u'_{1,\mathrm{rms}}(x)/u'_{1,\mathrm{rms}}(x_0)$")
    energy_norm.set(title="TKE decay about each curve's own reference", ylabel=r"$K(x)/K(x_0)$")
    figure.suptitle("M7: absolute and normalized decay; input amplitudes are preserved", fontsize=13)
    finish(plt, figure, list(axes.flat), output, "m7-decay",
           "Normalization is for display only. Table-4 K uses 3 u_rms^2 / 2; simulated K uses all three temporal variances. No experimental error bar is assumed.")
    return {"cases": plotted, "experimental_targets": targets,
            "normalization": "Each curve divided by its own reference-plane value; no amplitude changes in absolute data."}


def positive_bins(k, e):
    require(len(k) == len(e) and len(k) >= 2 and all(b > a for a, b in zip(k, k[1:])),
            "invalid increasing spectrum coordinates")
    kept = [(number(a), number(b)) for a, b in zip(k, e) if a > 0 and b > 0]
    require(kept and all(number(value) >= 0 for value in e), "spectrum needs nonnegative density and positive bins")
    return [p[0] for p in kept], [p[1] for p in kept], len(k) - len(kept)


def plot_spectra(plt, output, cases, targets, source):
    figure, axes = plt.subplots(1, 3, figsize=(15.0, 5.7), sharey=True)
    plotted = {}
    for column, (axis, target) in enumerate(zip(axes, targets), 1):
        station = cases[0]["reference_x_m"] + target["x_offset_m"]
        published = [(number(row[0]) * 100, number(row[column]) * 1e-6)
                     for row in source["table2a"]["rows"] if row[column] is not None]
        axis.loglog([p[0] for p in published], [p[1] for p in published], "D:",
                    color="#222222", markersize=4, linewidth=.9, label="CBC Table 2(a), full band")
        pane = {"theta": target["theta"], "station_x_m": station, "Table2a_SI": published, "cases": {}}
        nyquists = set()
        for case in cases:
            estimate = station_at(case["spectra"], station)["welch_primary"]
            if not estimate.get("available", True):
                pane["cases"][case["id"]] = {"available": False, "reason": estimate["reason"]}
                continue
            k, e, omitted = positive_bins(estimate["k1_rad_m"], estimate["E11_m3_s2"])
            axis.loglog(k, e, color=case["color"], linestyle=case["style"], linewidth=1.3, label=case["label"])
            nyquist = 2 * math.pi * case["analysis"]["window"]["sample_Nyquist_Hz"] / case["analysis"]["window"]["prescribed_convection_speed_m_s"]
            nyquists.add(nyquist)
            pane["cases"][case["id"]] = {"k1_rad_m": k, "E11_m3_s2": e,
                "omitted_nonpositive_or_zero_frequency_bins": omitted, "sample_Nyquist_rad_m": nyquist,
                "estimator": {key: estimate.get(key) for key in (
                    "window", "segment_samples", "overlap_samples", "segments", "frequency_resolution_Hz",
                    "equivalent_noise_bandwidth_Hz", "detrending", "parseval_relative_residual", "uncertainty")}}
        require(pane["cases"] and any(v.get("k1_rad_m") for v in pane["cases"].values()),
                "no available Welch estimate at a comparison station")
        cutoff = number(cases[0]["analysis"]["synthesis"]["cutoff_rad_m"])
        taper = number(cases[0]["analysis"]["synthesis"]["taper_start_rad_m"])
        axis.axvspan(taper, cutoff, color="#999999", alpha=.12)
        axis.axvline(cutoff, color="#777777", linestyle="--", linewidth=1, label="Input radial cutoff (not k1 resolution)")
        largest = max([p[0] for p in published] + list(nyquists))
        for i, value in enumerate(sorted(nyquists)):
            axis.axvline(value, color="#666666", linestyle=":", linewidth=1.2,
                         label="Temporal sample Nyquist" if i == 0 else None)
        if len(nyquists) == 1 and next(iter(nyquists)) < largest:
            axis.axvspan(next(iter(nyquists)), largest * 1.08, color="#CCCCCC", alpha=.16)
        axis.set(xlabel=r"Longitudinal $k_1=2\pi f/\overline{U}$ (rad/m)",
                 title=f"theta = {target['theta']:g}; x = {station:.6g} m", xlim=(1, largest * 1.08))
        axis.legend(fontsize=7, loc="lower left")
        plotted[str(int(target["theta"]))] = pane
    axes[0].set_ylabel(r"One-sided $F_{11}(k_1)$ ($\mathrm{m^3\,s^{-2}}$)")
    figure.suptitle("M7: streamwise temporal Welch spectra versus full-band longitudinal measurements", fontsize=13)
    finish(plt, figure, list(axes), output, "m7-longitudinal-spectra",
           "Periodic Hann Welch; full-record own-probe centering; F11 = Ubar Sf/(2 pi). Radial inlet taper is shaded. Low k1 does not imply a shell-band-matched target.")
    return plotted


def available_pairs(case, key):
    result = []
    for pair in case["analysis"]["mean_convection_delayed_pairs"]:
        estimate = pair[key]
        if estimate.get("available", False):
            value = estimate["full_window_normalized_covariance_components"][0]
            if value is not None:
                result.append((number(pair["xi_delta_x_over_M"]), number(value), pair, estimate))
    return sorted(result, key=lambda row: row[0])


def plot_correlations(plt, output, cases, table5):
    figure, axes = plt.subplots(2, 2, figsize=(13.2, 8.7))
    primary, secondary, overlap, normalization = axes.flat
    grid = next(row for row in table5["grids"] if row["id"] == "grid_5_08_cm")
    published = [(number(row[0]), number(row[1])) for row in grid["rows"] if row[0] is not None]
    for axis in (primary, secondary):
        axis.plot([p[0] for p in published], [p[1] for p in published], "D:",
                  color="#222222", markersize=4, linewidth=1, label="CBC Table 5, full band")
    plotted = {}
    for case in cases:
        primary_rows = available_pairs(case, "primary_joint_poststart")
        extended_rows = available_pairs(case, "secondary_stationary_forced_inlet_extension")
        require(primary_rows, f"{case['id']}: no available joint-poststart delayed correlations")
        for axis, rows in ((primary, primary_rows), (secondary, extended_rows)):
            if rows:
                axis.plot([row[0] for row in rows], [row[1] for row in rows], **line_style(case))
        overlap.plot([row[0] for row in primary_rows],
                     [row[3]["actual_paired_sample_span_s"] for row in primary_rows], **line_style(case))
        if extended_rows:
            overlap.plot([row[0] for row in extended_rows],
                         [row[3]["actual_paired_sample_span_s"] for row in extended_rows],
                         color=case["color"], linestyle=":", marker="^", markersize=3, linewidth=1)
        delta = [(row[0], row[3]["overlap_Pearson_components"][0] - row[1])
                 for row in primary_rows if row[3]["overlap_Pearson_components"][0] is not None]
        if delta:
            normalization.plot([row[0] for row in delta], [row[1] for row in delta], **line_style(case))
        plotted[case["id"]] = case["analysis"]["mean_convection_delayed_pairs"]
    for axis in axes.flat:
        axis.set_xscale("symlog", linthresh=1)
        axis.set_xticks((0, 1, 10, 100), ("0", "1", "10", "100"))
        axis.set_xlabel(r"$\xi=\Delta x/M=\overline{U}\Delta t/M$ (not $\Delta\theta$)")
        axis.legend(fontsize=7, loc="best")
    primary.set(title="Primary: both samples after recording start", ylabel="R11: full-window covariance normalization")
    secondary.set(title="Secondary: earlier forced-inlet reference", ylabel="R11: full-window covariance normalization")
    overlap.set(title="Actual paired spans; triangles: earlier inlet", ylabel="Paired sample span (s)")
    normalization.set(title="Primary normalization sensitivity", ylabel="Overlap Pearson minus full-window-normalized R11")
    normalization.axhline(0, color="#777777", linewidth=.8)
    figure.suptitle("M7: fixed-delay correlation versus full-band Table 5; no lag fitting", fontsize=13)
    finish(plt, figure, list(axes.flat), output, "m7-delayed-correlations",
           "Linear temporal interpolation at exact Delta x/Ubar. Earlier-inlet extension is secondary. Full-window covariance ratios need not stay in [-1,1]; no values clipped.")
    return {"published_Table5": grid, "cases": plotted,
            "normalization": "Primary and secondary full-window-normalized covariance; overlap Pearson is a separate sensitivity estimate."}


def plot_sampling(plt, output, cases):
    figure, axes = plt.subplots(2, 4, figsize=(16.0, 7.5), sharex=True)
    plotted, flags = {}, []
    for case in cases:
        x = sorted(case["qualification"])
        rows = [case["qualification"][station] for station in x]
        flags.append(f"{FAMILY_LABELS[case['family']]} {case['seed']}: "
                     f"{sum(row['engineering_screen'] != 'screens_passed_only' for row in rows)}/{len(rows)} flagged")
        for component in range(4):
            axes[0, component].plot(x, [row["full_to_final_half_change_percent"][component] for row in rows],
                                    **line_style(case))
            widths = [row["conditional_interval_halfwidth_percent_of_full"] for row in rows]
            middle = [row["8"][component] for row in widths]
            low, high = [], []
            for row in widths:
                values = [row[count][component] for count in ("4", "8", "16")]
                require(all(v is None or number(v) >= 0 for v in values), "invalid conditional interval halfwidth")
                low.append(min(values) if all(v is not None for v in values) else None)
                high.append(max(values) if all(v is not None for v in values) else None)
            axis = axes[1, component]
            axis.plot(x, middle, **line_style(case))
            kept = [i for i in range(len(x)) if low[i] is not None and high[i] is not None]
            if kept:
                # Fill only contiguous finite runs; never bridge undefined ratios.
                runs = []
                for i in kept:
                    if not runs or i != runs[-1][-1] + 1:
                        runs.append([])
                    runs[-1].append(i)
                for run in runs:
                    axis.fill_between([x[i] for i in run], [low[i] for i in run], [high[i] for i in run],
                                      color=case["color"], alpha=.10)
            undefined = [x[i] for i in range(len(x)) if low[i] is None or middle[i] is None]
            for station in undefined:
                axis.annotate("undefined", xy=(station, .98), xycoords=("data", "axes fraction"),
                              fontsize=7, color=case["color"], rotation=90, va="top")
        plotted[case["id"]] = rows
    for component, title in enumerate(QUANTITIES):
        for row in range(2):
            axes[row, component].axhline(2, color="#777777", linestyle=":", linewidth=1)
            axes[row, component].set_xlabel("x (m)")
        axes[0, component].set_title(title + ": full versus final half")
        axes[1, component].set_title(title + ": conditional halfwidth")
    axes[0, 0].set_ylabel("100*abs(full-half)/max(full,half) (%)")
    axes[1, 0].set_ylabel("CI halfwidth / full-window value (%)")
    axes[0, 0].legend(fontsize=7)
    figure.suptitle("M7: predeclared sampling screens; dotted guides = 2%\n" + "; ".join(flags), fontsize=10.5)
    finish(plt, figure, list(axes.flat), output, "m7-sampling-qualification",
           "Lower lines: 8-block conditional CI halfwidths; shading: 4/8/16-block halfwidth range. Block-energy and full-window estimators differ; passing is no independence proof.")
    return {"cases": plotted, "conditional_intervals": "8-block line; 4/8/16-block range, not an experimental or total numerical uncertainty band."}


def load_native(inputs, path, case):
    receipt = inputs.read(path)
    require(receipt["schema"] == "NEKTAR_HIT_CHECKPOINT_PROFILES_V1" and receipt["status"] == "completed"
            and receipt.get("finished_utc"), "native profile receipt is incomplete")
    require(Path(receipt["case"]).resolve() == case["directory"], "native receipt belongs to another case")
    manifest_hashes = [record["sha256"] for key, record in receipt["source_inputs"].items()
                       if Path(key).name == "manifest.json"]
    require(manifest_hashes == [inputs.files[str(case["manifest_path"])]["sha256"]],
            "native profile/analysis manifest hashes differ")
    snapshots = []
    window = case["analysis"]["window"]
    for profile in receipt["profiles"]:
        time = number(profile["time_s"])
        require(window["requested_start_s"] - 1e-10 <= time <= window["requested_end_s"] + 1e-10,
                "native profile is outside the qualified recording window")
        candidates = [(Path(key), digest) for key, digest in profile["output_files"].items()
                      if key.endswith(".slabs.csv")]
        require(len(candidates) == 1, "expected exactly one slab CSV per native snapshot")
        raw = inputs.raw(*candidates[0])
        rows = sorted(({key: number(value) for key, value in row.items()}
                       for row in csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))), key=lambda row: row["xmid"])
        require(len(rows) == profile["slabs"] == case["manifest"]["mesh"]["nx"]
                and len({row["xmid"] for row in rows}) == len(rows), "incomplete/duplicate native slabs")
        volume = math.fsum(row["volume"] for row in rows)
        require(all(row["volume"] > 0 for row in rows) and
                math.isclose(volume, profile["flow"]["volume"], rel_tol=1e-9), "native slab volume mismatch")
        for column, global_key in (
                ("molecular_strain_dissipation", "instantaneous_resolved_strain_dissipation"),
                ("svv_xy", "instantaneous_svv_work_xy"), ("svv_z", "instantaneous_svv_work_z")):
            average = math.fsum(row["volume"] * row[column] for row in rows) / volume
            require(math.isclose(average, profile["flow"][global_key], rel_tol=1e-8, abs_tol=1e-13),
                    "native slab/global work mismatch")
        snapshots.append({"time_s": time, "rows": rows, "field_build_revision": profile["field_build_revision"]})
    require(snapshots, "empty native profile receipt")
    x = [row["xmid"] for row in snapshots[0]["rows"]]
    require(all([row["xmid"] for row in snapshot["rows"]] == x for snapshot in snapshots),
            "native snapshot grids differ")
    result = {"x_m": x, "snapshot_times_s": sorted(s["time_s"] for s in snapshots),
              "field_build_revisions": sorted({s["field_build_revision"] for s in snapshots}), "columns": {}}
    for key in ("molecular_strain_dissipation", "svv_xy", "svv_z"):
        columns = [[row[key] for row in snapshot["rows"]] for snapshot in snapshots]
        result["columns"][key] = {"mean": [math.fsum(v) / len(v) for v in zip(*columns)],
                                  "min": [min(v) for v in zip(*columns)], "max": [max(v) for v in zip(*columns)]}
    return result


def plot_native(plt, output, cases, profiles, targets):
    figure, axes = plt.subplots(1, 3, figsize=(15.0, 5.5))
    columns = (("molecular_strain_dissipation", r"Molecular $2\nu S:S$"),
               ("svv_xy", "Native x-y SVV work"), ("svv_z", "Native z SVV work"))
    for axis, (key, title) in zip(axes, columns):
        positive = True
        for case in cases:
            if case["id"] not in profiles:
                continue
            profile = profiles[case["id"]]
            values = profile["columns"][key]
            positive = positive and all(value > 0 for value in values["min"])
            axis.plot(profile["x_m"], values["mean"], color=case["color"], linestyle=case["style"],
                      linewidth=1.2, label=case["label"])
            axis.fill_between(profile["x_m"], values["min"], values["max"], color=case["color"], alpha=.12)
        if key == "molecular_strain_dissipation":
            axis.plot([cases[0]["reference_x_m"] + row["x_offset_m"] for row in targets],
                      [row["epsilon_m2_s3"] for row in targets], "D:", color="#222222", markersize=4,
                      label="CBC Table-4 full-band decay inference")
        if positive:
            axis.set_yscale("log")
        axis.set(title=title, xlabel="x slab midpoint (m)", ylabel=r"Instantaneous work ($\mathrm{m^2\,s^{-3}}$)")
        axis.legend(fontsize=7, loc="best")
    figure.suptitle("M7: native late-checkpoint work; means and snapshot ranges", fontsize=13)
    finish(plt, figure, list(axes), output, "m7-native-work",
           "Shading is snapshot min/max, not a confidence interval. Work uses full instantaneous fields; transport, pressure and time-integration terms remain open.")


def inside(path, directory):
    return path == directory or directory in path.parents


def output_preflight(output, cases, inputs, explicit_roots):
    require(not output.exists(), "--output must be a fresh directory")
    roots = {case["directory"] for case in cases} | {HERE}
    roots.update(Path(path).resolve() for path in explicit_roots)
    for anchor in tuple(roots):
        for parent in (anchor, *anchor.parents):
            if (parent / "CMakeLists.txt").is_file() and (parent / "solvers").is_dir():
                roots.add(parent)
                break
    require(all(not inside(output, root) for root in roots), "output must be outside case/source directories")
    require(all(not inside(Path(path), output) for path in inputs.files), "output would contain an input")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis", type=Path, nargs="+", required=True, help="Four analyze_m7.py JSON files")
    parser.add_argument("--experimental-source", type=Path, default=HERE.parent / "cbc1971-source.json")
    parser.add_argument("--correlation-source", type=Path, default=HERE.parent / "cbc1971_full_band_correlation.json")
    parser.add_argument("--native", action="append", default=[], metavar="CASE_ID=RECEIPT_JSON")
    parser.add_argument("--source-root", type=Path, action="append", default=[], help="Additional protected source checkout root")
    parser.add_argument("--output", type=Path, required=True, help="Fresh directory outside case/source trees")
    args = parser.parse_args()
    try:
        inputs = Inputs()
        source_path, correlation_path = args.experimental_source.resolve(), args.correlation_source.resolve()
        source, table5 = inputs.read(source_path), inputs.read(correlation_path)
        require(source["schema"] == "nektar.hit.cbc1971.source.v1" and
                table5["schema"] == "nektar.hit.cbc1971.full_band_correlation.v1", "unsupported experimental source schema")
        require(source["source"]["pdf_sha256"] == table5["source"]["pdf_sha256"], "experimental PDFs differ")
        cases = validate_suite([load_case(inputs, path, source_path, correlation_path, source, table5)
                                for path in args.analysis])
        targets = table_targets(source)
        profiles = {}
        for assignment in args.native:
            name, separator, path = assignment.partition("=")
            require(separator and name not in profiles, "native inputs must be unique CASE_ID=RECEIPT_JSON assignments")
            matches = [case for case in cases if case["id"] == name]
            require(len(matches) == 1, f"native receipt references unknown case: {name}")
            profiles[name] = load_native(inputs, Path(path), matches[0])
        output = args.output.resolve()
        output_preflight(output, cases, inputs, args.source_root + [source_path.parent, correlation_path.parent])
        for target in targets:
            for case in cases:
                station_at(case["planes"], case["reference_x_m"] + target["x_offset_m"])
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.rcParams.update({"font.size": 10, "svg.fonttype": "none", "axes.titlesize": 10.5,
                             "axes.labelsize": 10, "legend.framealpha": .85})
        output.mkdir(parents=True, exist_ok=False)
        plotted = {"decay": plot_decay(plt, output, cases, targets),
                   "longitudinal_spectra": plot_spectra(plt, output, cases, targets, source),
                   "delayed_correlations": plot_correlations(plt, output, cases, table5),
                   "sampling_qualification": plot_sampling(plt, output, cases)}
        names = ["m7-decay", "m7-longitudinal-spectra", "m7-delayed-correlations", "m7-sampling-qualification"]
        if profiles:
            plot_native(plt, output, cases, profiles, targets)
            plotted["native_work"] = profiles
            names.append("m7-native-work")
        inputs.verify_unchanged()
        figure_files = {str(path.name): {"sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "bytes": path.stat().st_size}
                        for name in names for suffix in ("png", "svg") for path in (output / f"{name}.{suffix}",)}
        provenance = {"schema": "NEKTAR_HIT_M7_FIGURES_V1",
            "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "plotting_script": {"path": str(Path(__file__).resolve()), "sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
            "arguments": {key: [str(v) for v in value] if isinstance(value, list) else str(value) if isinstance(value, Path) else value
                          for key, value in vars(args).items()},
            "software": {"matplotlib": matplotlib.__version__, "backend": "Agg"},
            "verified_input_files": inputs.files, "source_pdf_sha256_declared_by_transcriptions": source["source"]["pdf_sha256"],
            "analysis_declared_upstream_hashes_not_rehashed_by_plotter": {case["id"]: case["analysis"]["sha256"] for case in cases},
            "case_order": [{"id": case["id"], "family": case["family"], "seed": case["seed"],
                            "label": case["label"], "color": case["color"], "line_style": case["style"]} for case in cases],
            "case_configurations": {case["id"]: {key: case["manifest"][key]
                                    for key in ("physics", "synthesis", "mesh", "time", "svv")} for case in cases},
            "figures": figure_files, "plotted_data": plotted,
            "native_profiles_omitted_for_cases": [case["id"] for case in cases if case["id"] not in profiles],
            "claim_limits": [SCOPE,
                "Four curves use two stream identities common across spectrum families; no four-replicate ensemble claim.",
                "No amplitude rescaling, downstream fitting, fitted correlation delay or universal experimental error bar.",
                "Normalized decay uses each curve's own reference-plane value and does not repair missing inlet scales.",
                "Welch confidence intervals are not inferred from correlated probes or overlapping segments.",
                "Input radial cutoff differs from temporal Nyquist, longitudinal spectral support and downstream generated scales.",
                "Full-window covariance normalization and paired-overlap Pearson are separate estimators.",
                "Conditional block-energy intervals do not quantify omitted scales, all numerical error or established coverage.",
                "Native work omits open-boundary transport, pressure and time-integration budget terms."]}
        finite_tree(provenance)
        (output / "provenance.json").write_text(json.dumps(provenance, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        print(f"Wrote {len(names)} PNG/SVG figure pairs and provenance to {output}")
    except (OSError, ValueError, KeyError, TypeError, StopIteration, ImportError, OverflowError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
