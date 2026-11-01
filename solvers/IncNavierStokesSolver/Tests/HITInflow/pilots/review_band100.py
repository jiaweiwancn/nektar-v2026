#!/usr/bin/env python3
"""Review immutable CBC100 records without launching or repeating CFD.

Run on Linux with one CPU/thread. Recompute recorded moments and spectra from
the hashed native-history caches, preserve all original outputs, and report
window/probe sensitivity separately from matched numerical refinement.
"""
import argparse
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys

sys.dont_write_bytecode = True


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4*1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    os.replace(tmp, path)


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def same_numbers(left, right):
    """Exact deterministic reanalysis comparison; report any floating difference."""
    report = {"numbers_compared": 0, "different_numbers": 0, "maximum_absolute_difference": 0.0}

    def compare(a, b):
        if isinstance(a, dict):
            require(isinstance(b, dict) and a.keys() == b.keys(), "reanalysis object keys differ")
            for key in a:
                compare(a[key], b[key])
        elif isinstance(a, (list, tuple)):
            require(isinstance(b, (list, tuple)) and len(a) == len(b), "reanalysis array lengths differ")
            for x, y in zip(a, b):
                compare(x, y)
        elif isinstance(a, (int, float)) and not isinstance(a, bool):
            require(isinstance(b, (int, float)) and math.isfinite(a) and math.isfinite(b),
                    "invalid reanalysis number")
            difference = float(abs(a-b))
            report["numbers_compared"] += 1
            report["different_numbers"] += int(difference != 0)
            report["maximum_absolute_difference"] = max(report["maximum_absolute_difference"], difference)
        else:
            require(a == b, "reanalysis nonnumeric field differs")
    compare(left, right)
    report["exact_numeric_match"] = report["different_numbers"] == 0
    return report


def summarize_spectra(spectra, stats, band100):
    import numpy as np
    rows = []
    for station, plane in zip(spectra, stats["planes"]):
        components = {}
        for c, name in enumerate(band100.COMPONENTS):
            estimates = station["components"][name]
            box = estimates["boxcar"]
            k, f = np.asarray(box["k1_rad_m"]), np.asarray(box["F_component_m3_s2"])
            parts = [band100.integrate_density(k, f, low, high)
                     for low, high in ((0, 10), (10, 100), (100, float(k[-1])))]
            variance = plane["variance_m2_s2"][c]
            closure = 100*(sum(parts)/variance-1) if variance else 0.0
            components[name] = {
                "variance_m2_s2": variance,
                "welch_band_powers_m2_s2": [p["power_m2_s2"] for p in estimates["welch"]["band_powers"]],
                "welch_F100_m3_s2": estimates["welch"]["F_component_100_m3_s2"],
                "boxcar_partition_0_10_100_max_m2_s2": parts,
                "partition_trapezoidal_vs_variance_percent": closure,
                "boxcar_parseval_residual": box["parseval_relative_residual"],
            }
        rows.append({"station_x_m": station["station_x_m"], "K_m2_s2": plane["K_m2_s2"],
                     "components": components})
    return rows


def derive(times, values, stations, dt, speed, band100, m7):
    stats, _ = m7.station_statistics(times, values, stations, 256)
    spectra = [{"station_x_m": x, "components": band100.component_spectra(values[:, s], dt, speed)}
               for s, x in enumerate(stations)]
    for s, station in enumerate(spectra):
        for c, component in enumerate(band100.COMPONENTS):
            box = station['components'][component]['boxcar']
            variance = stats['planes'][s]['variance_m2_s2'][c]
            integral = box['F_component_integral_m2_s2']
            residual = abs(integral-variance)/variance if variance else abs(integral)
            require(residual <= 5e-10, 'recomputed boxcar PSD does not recover variance')
            box['relative_difference_from_full_window_variance'] = residual
    return stats, spectra


def comparison_shell(original, stats, spectra, times):
    own = dict(original)
    window = dict(original["window"], samples=len(times),
                  requested_start_s=float(times[0]), requested_end_s=float(times[-1]),
                  first_last_samples_s=[float(times[0]), float(times[-1])])
    own.update(window=window, statistics=stats, spectra=spectra,
               control_identity=dict(original["control_identity"], window=window))
    return own


def sampling_changes(current, reference, band100, scope):
    """Descriptive changes across deliberately different sampling definitions."""
    rows, moments, powers = [], [], []
    for a, b, sa, sb in zip(current['statistics']['planes'], reference['statistics']['planes'],
                             current['spectra'], reference['spectra']):
        require(a['station_x_m'] == b['station_x_m'] == sa['station_x_m'] == sb['station_x_m'],
                'sampling comparison stations differ')
        quantities = {key: band100.difference(x, y) for key, x, y in zip(
            ('u_variance', 'v_variance', 'w_variance', 'K'),
            a['variance_m2_s2'] + [a['K_m2_s2']], b['variance_m2_s2'] + [b['K_m2_s2']])}
        moments.extend(v['signed_relative_change_percent'] for v in quantities.values())
        bands, endpoints = {}, {}
        for c in band100.COMPONENTS:
            own, other = sa['components'][c]['welch'], sb['components'][c]['welch']
            bands[c] = []
            for x, y in zip(own['band_powers'], other['band_powers']):
                require(x['k1_limits_rad_m'] == y['k1_limits_rad_m'], 'sampling physical bands differ')
                change = band100.difference(x['power_m2_s2'], y['power_m2_s2'])
                bands[c].append({'k1_limits_rad_m': x['k1_limits_rad_m'], **change})
                powers.append(change['signed_relative_change_percent'])
            endpoints[c] = band100.difference(own['F_component_100_m3_s2'], other['F_component_100_m3_s2'])
        rows.append({'station_x_m': a['station_x_m'], 'variances_and_K': quantities,
                     'welch_band_power_changes': bands, 'welch_endpoint_100_changes': endpoints})
    def maximum(values):
        return None if any(x is None for x in values) else max(abs(x) for x in values)
    return {'kind': 'sampling_sensitivity', 'scope': scope, 'stations': rows,
            'eligible_for_numerical_refinement_screen': False,
            'maximum_absolute_variance_K_change_percent': maximum(moments),
            'maximum_absolute_band_power_change_percent': maximum(powers)}


def density_changes(current, reference):
    import numpy as np
    rows = []
    for own, other in zip(current["spectra"], reference["spectra"]):
        values = {}
        for c in ("u", "v", "w"):
            a, b = own["components"][c]["welch"], other["components"][c]["welch"]
            k, f, g = map(np.asarray, (a["k1_rad_m"], a["F_component_m3_s2"], b["F_component_m3_s2"]))
            selected = (k >= 10) & (k <= 100) & (g > 0)
            q = 100*(f[selected]/g[selected]-1)
            at = int(np.argmax(np.abs(q)))
            values[c] = {"maximum_absolute_change_percent_at_recorded_bins": float(abs(q[at])),
                         "signed_change_at_max_percent": float(q[at]), "k1_at_max_rad_m": float(k[selected][at])}
        rows.append({"station_x_m": own["station_x_m"], "components": values})
    return rows


def plot(result, full, output):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors = {"base": "#222222", "dt": "#1684a7", "x": "#d46a14"}
    labels = {"base": "Base: hx=0.1m, dt=25us", "dt": "Half dt:12.5us", "x": "Finer x: hx=0.05m"}
    paths = []
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.2), sharex=True, sharey=True)
    for s, axis in enumerate(axes.flat):
        for name in ("base", "dt", "x"):
            spec = full[name]["spectra"][s]["components"]["u"]["welch"]
            k, f = np.asarray(spec["k1_rad_m"]), np.asarray(spec["F_component_m3_s2"])
            keep = (k >= 5) & (k <= 105)
            axis.loglog(k[keep], f[keep], color=colors[name], ls="--" if name == "dt" else "-", label=labels[name])
        if s == 0:
            measured = full["base"]["upstream_experimental_overlay"]["Table2a_F11_SI"]
            axis.loglog([p[0] for p in measured], [p[1] for p in measured], "D:", color="#8a589c", ms=3,
                        label="CBC theta42; raw full-band data")
        axis.axvline(100, color="gray", ls=":", lw=.8)
        axis.set(title=f"x={full['base']['spectra'][s]['station_x_m']:g}m", xlim=(5, 105))
        axis.grid(True, which="both", alpha=.25)
    axes[0, 0].legend(fontsize=7)
    for axis in axes[1]: axis.set_xlabel("Longitudinal k1 (rad/m)")
    for axis in axes[:, 0]: axis.set_ylabel("One-sided F11 (m³/s²)")
    fig.suptitle("CBC100: completed temporal Welch spectra")
    fig.text(.02, .012, "Sampled probes; one realization. Interior pilot stations have no assigned experimental reference.", fontsize=8)
    fig.tight_layout(rect=(0, .04, 1, .96))
    paths.append(("review-spectra", fig))

    fig, axes = plt.subplots(2, 2, figsize=(10, 6.5), sharex=True)
    for s, axis in enumerate(axes.flat):
        reference = full["base"]["spectra"][s]["components"]["u"]["welch"]
        k, g = np.asarray(reference["k1_rad_m"]), np.asarray(reference["F_component_m3_s2"])
        for name in ("dt", "x"):
            f = np.asarray(full[name]["spectra"][s]["components"]["u"]["welch"]["F_component_m3_s2"])
            keep = (k >= 10) & (k <= 105) & (g > 0)
            axis.plot(k[keep], 100*(f[keep]/g[keep]-1), color=colors[name], label=labels[name])
        axis.axhspan(-5, 5, color="gray", alpha=.12)
        axis.axhline(0, color="gray", lw=.8)
        axis.axhline(5, color="gray", lw=.7, ls=":")
        axis.axhline(-5, color="gray", lw=.7, ls=":")
        axis.axvline(100, color="gray", ls=":", lw=.8)
        axis.set(title=f"x={full['base']['spectra'][s]['station_x_m']:g}m", xlim=(10, 105))
        axis.grid(alpha=.25)
    axes[0, 0].legend(fontsize=7)
    for axis in axes[1]: axis.set_xlabel("Longitudinal k1 (rad/m)")
    for axis in axes[:, 0]: axis.set_ylabel("F11 change from base (%)")
    fig.suptitle("Matched refinements: pointwise spectral sensitivity")
    fig.text(.02, .012, "Shading is a5% guide; the registered5% screen uses integrated physical bands, not every spectral bin.", fontsize=8)
    fig.tight_layout(rect=(0, .04, 1, .96))
    paths.append(("review-refinement", fig))

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    xx = [p["station_x_m"] for p in full["base"]["statistics"]["planes"]]
    for name in ("base", "dt", "x"):
        data = result["cases"][name]
        axes[0].plot(xx, [p["K_m2_s2"] for p in full[name]["statistics"]["planes"]],
                     "o-", color=colors[name], label=labels[name])
        axes[1].plot(xx, [q["full_to_final_half_change_percent"][0] for q in data["sampling_qualification"]],
                     "o-", color=colors[name])
        axes[2].plot(xx, [max(v for row in q["conditional_interval_halfwidth_percent_of_full"].values() for v in row)
                          for q in data["sampling_qualification"]], "o-", color=colors[name])
    axes[0].set(ylabel="Sampled K (m²/s²)", title="Full-band kinetic energy")
    axes[1].set(ylabel="Change (%)", title="u variance: full vs final half")
    axes[2].set(ylabel="Conditional halfwidth (%)", title="Largest4/8/16-block halfwidth")
    for axis in axes[1:]: axis.axhline(2, color="gray", ls=":", lw=.8)
    for axis in axes: axis.set_xlabel("Distance from inlet (m)"); axis.grid(alpha=.25)
    axes[0].legend(fontsize=7)
    fig.suptitle("Small low-band differences coexist with energy and sampling failures")
    fig.text(.01, .015, "Intervals are conditional diagnostics; block/probe independence and confidence coverage are not established.", fontsize=8)
    fig.tight_layout(rect=(0, .045, 1, .93))
    paths.append(("review-energy-sampling", fig))
    figures = []
    for stem, fig in paths:
        for extension in ("png", "svg"):
            p = output/(stem + "." + extension)
            fig.savefig(p, dpi=165)
            figures.append({"file": p.name, "sha256": sha(p)})
        plt.close(fig)
    return figures


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--cpu", type=int, required=True)
    args = parser.parse_args()
    require(sys.platform.startswith("linux"), "review calculations run on authoritative Linux")
    require(args.cpu in os.sched_getaffinity(0), "selected review CPU is unavailable")
    os.sched_setaffinity(0, {args.cpu})
    study, source, output = (p.resolve() for p in (args.study, args.source_root, args.output))
    require(not output.exists() and not output.is_relative_to(study) and not output.is_relative_to(source),
            "output must be fresh and outside the source and study")
    require(output.is_relative_to(Path('/media/jay/新加卷')), "review outputs belong on the Linux data disk")
    sys.path.insert(0, str(source/'solvers/IncNavierStokesSolver/Tests/HITInflow/pilots'))
    import numpy as np
    import analyze_band100 as band100
    import analyze_m7 as m7
    output.mkdir(parents=True)
    os.environ["MPLCONFIGDIR"] = str(output/'mplconfig')
    receipt = {"status": "running", "started_utc": utc(), "command": sys.argv,
               "cpu_affinity": sorted(os.sched_getaffinity(0)), "driver_sha256": sha(__file__),
               "thread_environment": {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}}
    write(output/'receipt.json', receipt)
    protected = {}

    def bind(path, expected=None):
        path = Path(path)
        actual = sha(path)
        require(expected is None or actual == expected, "input hash mismatch: " + str(path))
        require(str(path) not in protected or protected[str(path)] == actual,
                'input changed between binding checks: ' + str(path))
        protected[str(path)] = actual
        return actual

    try:
        parent = read(study/'receipt.json')
        require(parent["status"] == "COMPLETED" and set(parent["cases"]) == {"base", "dt", "x"},
                "all three completed pilot cases are required")
        bind(study/'receipt.json')
        for item in parent["commands"]:
            p = Path(item["receipt"])
            bind(p, item["sha256"])
            cmd = read(p)
            require(item["status"] == cmd["status"] == "COMPLETED" and cmd["exit_code"] == 0, "failed command")
            bind(p.with_name(p.name.replace('.receipt.json', '.log')), cmd["log_sha256"])
            for filename, expected in cmd["input_files_sha256"].items(): bind(filename, expected)
        for name, expected in parent['bindings']['source']['files_sha256'].items(): bind(source/name, expected)
        for p in (Path(band100.__file__), Path(m7.__file__), Path(m7.ah.__file__), Path(__file__)): bind(p)
        result = {"schema": "NEKTAR_HIT_BAND100_REVIEW_V1", "created_utc": utc(), "study": str(study),
                  "source_revision": parent['bindings']['source']['revision'], "cases": {}, "comparisons": {},
                  "scope": "Numerical postprocessing reproduced; CFD itself was not rerun. Engineering screens are not experimental fidelity or formal confidence guarantees."}
        full, windows, grid = {}, {}, {}
        for name in ("base", "dt", "x"):
            case = parent["cases"][name]
            directory = Path(case["directory"])
            manifest_path = directory/'manifest.json'
            bind(manifest_path, case["manifest_sha256"])
            manifest = read(manifest_path)
            analysis_path = Path(case["analysis"])
            bind(analysis_path, case["analysis_sha256"])
            analysis = band100.read_json(analysis_path)
            ar = read(analysis_path.parent/'receipt.json')
            require(ar['status'] == 'completed' and ar['analysis_sha256'] == sha(analysis_path), "analysis receipt mismatch")
            bind(analysis_path.parent/'receipt.json')
            bind(analysis_path.parent/'history.npz', ar['cache_sha256'])
            bind(analysis_path.parent/'history.npz.json')
            for filename, expected in analysis['inputs']['sha256'].items(): bind(filename, expected)
            times, velocity, _ = m7.load_cache(analysis_path.parent/'history.npz', manifest_path, manifest)
            require(tuple(velocity.shape[2:4]) == (32, 32), 'review expects the recorded32x32 probe grid')
            clock = analysis["window"]
            ids = np.flatnonzero((times >= clock['requested_start_s']-2e-12) & (times <= clock['requested_end_s']+2e-12))
            times, velocity = times[ids[0]:ids[-1]+1], velocity[ids[0]:ids[-1]+1]
            require(len(times) == clock['samples'] == 1601, "record window mismatch")
            stations = manifest['sampling']['stations_x_m']
            dt, speed = manifest['time']['sample_dt_s'], manifest['physics']['mean_velocity_m_s']
            stats, spectra = derive(times, velocity, stations, dt, speed, band100, m7)
            reproduced = same_numbers({'statistics': stats, 'spectra': spectra},
                                      {'statistics': analysis['statistics'], 'spectra': analysis['spectra']})
            require(reproduced['exact_numeric_match'], "primary deterministic postprocessing differs")
            solver = read(directory/'solver.receipt.json')
            log = (directory/'solver.log').read_text(errors='replace')
            require(solver['final_progress']['step'] == manifest['time']['steps'] and
                    abs(solver['final_progress']['simulated_time_s']-.6) < 1e-12, "solver endpoint mismatch")
            full[name] = comparison_shell(analysis, stats, spectra, times)
            windows[name], grid[name] = {}, {}
            for label, first in (("final_half", len(times)//2), ("last_75_percent", len(times)//4)):
                ss, pp = derive(times[first:], velocity[first:], stations, dt, speed, band100, m7)
                windows[name][label] = comparison_shell(analysis, ss, pp, times[first:])
            for count in (8, 16):
                stride = 32//count
                ss, pp = derive(times, velocity[:, :, ::stride, ::stride], stations, dt, speed, band100, m7)
                grid[name][str(count)] = {'statistics': ss, 'spectra': pp}
            result['cases'][name] = {
                'analysis_path': str(analysis_path), 'reanalysis': reproduced,
                'mesh': manifest['mesh'], 'dt_s': manifest['time']['dt_s'], 'window': clock,
                'solver_wall_hours': solver['wall_seconds']/3600,
                'sampled_peak_rss_gib': solver['peak_sampled_tree_rss_kib']/1024**2,
                'final_progress': solver['final_progress'], 'native_affinity': solver['affinity_audit'],
                'sampling_qualification': stats['qualification'],
                'primary': summarize_spectra(spectra, stats, band100),
                'window_sensitivity': {label: sampling_changes(shell, full[name], band100, 'Different recording windows at identical physical probes.')
                                       for label, shell in windows[name].items()},
                'probe_grid_sensitivity': {count: sampling_changes(shell, full[name], band100, 'Nested probe subsets at the same recording window.')
                                           for count, shell in grid[name].items()},
                'probe_grid_scope': 'Nested8/16/32 grids from the same32-grid record; agreement cannot establish convergence to a denser64/128-grid recording.',
                'welch_metadata': {k: spectra[0]['components']['u']['welch'][k] for k in ('segments','segment_samples','overlap_samples','excluded_tail_samples','frequency_resolution_Hz')},
                'warnings': sorted(set(line.strip() for line in log.splitlines() if 'deprecated' in line or 'failed to converge' in line or 'Exceeded maxIt' in line)),
            }
            print('reproduced', name, reproduced, flush=True)
            del velocity
        for name in ('dt', 'x'):
            primary = band100.compare_analyses(full[name], full['base'])
            saved = read(parent['cases'][name]['analysis'])['reference_comparison']
            require(primary['matched_controls_eligible'], 'unmatched controls')
            for key in ('maximum_absolute_variance_K_change_percent', 'maximum_absolute_band_power_change_percent'):
                require(primary[key] == saved[key], 'saved comparison differs')
            result['comparisons'][name] = {
                'primary': primary, 'pointwise_density': density_changes(full[name], full['base']),
                'window_robustness': {label: band100.compare_analyses(windows[name][label], windows['base'][label])
                                      for label in windows[name]},
            }
        result['figures'] = plot(result, full, output)
        for filename, expected in protected.items():
            require(sha(filename) == expected, 'input changed during review: ' + filename)
        result['inputs_sha256'] = protected
        result['software'] = {'python': sys.version, 'numpy': np.__version__}
        write(output/'review.json', result)
        receipt.update(status='completed', finished_utc=utc(), review_sha256=sha(output/'review.json'))
        write(output/'receipt.json', receipt)
        print('completed', output, flush=True)
    except BaseException as error:
        receipt.update(status='failed', finished_utc=utc(), error=f'{type(error).__name__}: {error}')
        write(output/'receipt.json', receipt)
        raise


if __name__ == '__main__':
    main()
