#!/usr/bin/env python3
"""Prepare an affine, periodic-y/Homogeneous1D-z CBC upstream-spectrum pilot.

This standard-library tool creates inputs only. Numerical and performance
validation belongs on the Linux host specified by AGENTS.md.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import shutil
import sys
import xml.etree.ElementTree as ET


HERE = Path(__file__).resolve().parent
SPECTRUM = HERE.parent / "cbc1971_theta42_spectrum_3d.dat"
U, NU = 12.7, 1.5e-5
VARIABLES = ("u", "v", "w", "p")


def number(value):
    return format(value, ".17g")


def sub(parent, tag, text=None, **attributes):
    child = ET.SubElement(parent, tag, attributes)
    if text is not None:
        child.text = str(text)
    return child


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def finite_positive(value, label, zero=False):
    if not math.isfinite(value) or (value < 0 if zero else value <= 0):
        raise ValueError(f"{label} must be finite and {'nonnegative' if zero else 'positive'}")


def integer_steps(time, dt, label):
    steps = round(time / dt)
    if not math.isclose(steps * dt, time, rel_tol=1e-11, abs_tol=1e-14 * dt):
        raise ValueError(f"{label} must be an integer multiple of dt")
    return steps


def tuple_ints(value, count, label):
    try:
        result = tuple(int(part.strip()) for part in value.split(","))
    except ValueError as error:
        raise ValueError(f"{label} must contain {count} comma-separated integers") from error
    if len(result) != count:
        raise ValueError(f"{label} must contain {count} comma-separated integers")
    return result


def strict_count(kmax, period):
    count = max(4, 2 * (math.floor(kmax * period / (2 * math.pi)) + 1))
    while not kmax < math.pi * count / period:
        count += 2
    return count


def geometry(root, length, width, nx, ny):
    """Use the same edge orientation as generate_m4_cases.py/ChanFlow_3DH1D."""
    geo = sub(root, "GEOMETRY", DIM="2", SPACE="2")
    vertices = sub(geo, "VERTEX")
    for iy in range(ny + 1):
        for ix in range(nx + 1):
            sub(vertices, "V", f"{ix * length / nx:.17e} {iy * width / ny:.17e} {0:.17e}",
                ID=str(iy * (nx + 1) + ix))
    edges = sub(geo, "EDGE")
    horizontal, vertical = {}, {}
    index = 0
    for iy in range(ny + 1):
        for ix in range(nx):
            first = iy * (nx + 1) + ix
            horizontal[ix, iy] = index
            sub(edges, "E", f"{first} {first + 1}", ID=str(index))
            index += 1
    for iy in range(ny):
        for ix in range(nx + 1):
            first = iy * (nx + 1) + ix
            vertical[ix, iy] = index
            sub(edges, "E", f"{first} {first + nx + 1}", ID=str(index))
            index += 1
    elements = sub(geo, "ELEMENT")
    for iy in range(ny):
        for ix in range(nx):
            sides = (horizontal[ix, iy], vertical[ix + 1, iy],
                     horizontal[ix, iy + 1], vertical[ix, iy])
            sub(elements, "Q", " ".join(map(str, sides)), ID=str(iy * nx + ix))
    composites = sub(geo, "COMPOSITE")
    sub(composites, "C", f"Q[0-{nx * ny - 1}]" if nx * ny > 1 else "Q[0]", ID="0")
    borders = ([vertical[0, j] for j in range(ny)],
               [vertical[nx, j] for j in range(ny)],
               [horizontal[i, 0] for i in range(nx)],
               [horizontal[i, ny] for i in range(nx)])
    for index, ids in enumerate(borders, 1):
        sub(composites, "C", "E[" + ",".join(map(str, ids)) + "]", ID=str(index))
    sub(geo, "DOMAIN", "C[0]")


def add_function(conditions, name, values):
    function = sub(conditions, "FUNCTION", NAME=name)
    for variable, value in zip(VARIABLES, values):
        sub(function, "E", VAR=variable, VALUE=value)
    return function


def prepare(args):
    for label in ("length", "width", "h", "dt", "duration", "cutoff", "synthesis_length"):
        finite_positive(getattr(args, label), label)
    finite_positive(args.record_start, "record-start", zero=True)
    finite_positive(args.buffer, "buffer", zero=True)
    if args.modes < 4 or args.hom_modes < 4 or args.hom_modes % 2:
        raise ValueError("modes must be >=4; hom-modes must be even and >=4")
    if args.sample_every < 1 or not 0 <= args.seed < 2**64:
        raise ValueError("sample-every must be positive; seed must fit uint64")
    if args.checkpoint_every < 0:
        raise ValueError("checkpoint-every must be nonnegative")
    steps = integer_steps(args.duration, args.dt, "duration")
    start_step = integer_steps(args.record_start, args.dt, "record-start")
    if steps < 2 or not 0 <= start_step < steps:
        raise ValueError("duration needs >=2 steps and record-start must precede duration")
    if start_step % args.sample_every or steps % args.sample_every:
        raise ValueError("record-start and duration steps must be multiples of sample-every")
    # The solver accumulates floating-point timesteps. Keep its final callback
    # inside a declared record without weakening the core's strict time guard.
    synthetic_duration = args.duration + 2 * args.dt
    if synthetic_duration >= args.synthesis_length / U:
        raise ValueError("REPEAT=Error requires duration + 2*dt < synthesis-length / 12.7")
    taper = 0.8 * args.cutoff if args.taper is None else args.taper
    if not math.isfinite(taper) or not 0 <= taper < args.cutoff:
        raise ValueError("taper is an absolute wavenumber and must satisfy 0<=taper<cutoff")
    periods = (args.synthesis_length, args.width, args.width)
    counts = (tuple(strict_count(args.cutoff, length) for length in periods)
              if args.synthesis_counts.lower() == "auto" else
              tuple_ints(args.synthesis_counts, 3, "synthesis-counts"))
    if any(n < 4 or n % 2 for n in counts):
        raise ValueError("all synthesis counts must be even and >=4")
    if any(not args.cutoff < math.pi * n / length for n, length in zip(counts, periods)):
        raise ValueError("cutoff must be strictly below every synthesis Nyquist wavenumber")
    if not args.cutoff < math.pi * args.hom_modes / args.width:
        raise ValueError("cutoff must be strictly below the solver z Nyquist wavenumber")
    try:
        stations = [float(value) for value in args.stations.split(",")]
    except ValueError as error:
        raise ValueError("stations must be comma-separated x coordinates in metres") from error
    if not stations or len(set(stations)) != len(stations):
        raise ValueError("stations must be nonempty and unique")
    stations.sort()
    station_tolerance = 64 * sys.float_info.epsilon * max(args.length, args.width)
    if args.buffer >= args.length or any(not math.isfinite(x) or x < 0 or
                                         x > args.length - args.buffer + station_tolerance for x in stations):
        raise ValueError("all stations must lie in [0,length-buffer], with buffer<length")
    probe_ny, probe_nz = tuple_ints(args.probe_grid, 2, "probe-grid")
    if min(probe_ny, probe_nz) < 2:
        raise ValueError("each half-open probe-grid count must be >=2")
    for label in ("svv_coeff", "svv_z_coeff"):
        finite_positive(getattr(args, label), label, zero=True)
    for label in ("svv_cutoff", "svv_z_cutoff"):
        if not math.isfinite(getattr(args, label)) or not 0 <= getattr(args, label) < 1:
            raise ValueError(f"{label} must lie in [0,1)")
    # The native cumulative filter supports a fixed output frequency, not an
    # arbitrary list of output times. Bound that I/O explicitly.
    stress_every = start_step if start_step else steps + 1
    periodic_stress_steps = list(range(stress_every, steps + 1, stress_every))
    if len(periodic_stress_steps) > 16:
        raise ValueError("native stress snapshot schedule exceeds 16 writes; increase record-start")
    initial = None
    initial_hashes = {}
    if args.initial_field:
        initial = args.initial_field.expanduser().resolve()
        if initial.is_file():
            initial_hashes[initial.name] = digest(initial)
        elif initial.is_dir() and (initial / "Info.xml").is_file():
            files = sorted((path for path in initial.rglob("*") if path.is_file()),
                           key=lambda path: path.relative_to(initial).as_posix())
            initial_hashes = {path.relative_to(initial).as_posix(): digest(path) for path in files}
        else:
            raise ValueError("initial-field must be an existing native field file or MPI field directory containing Info.xml")
    output = args.output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    targets = ("pilot.xml", "manifest.json", SPECTRUM.name)
    if any((output / name).exists() for name in targets):
        raise ValueError("output already contains pilot inputs; use a new directory")
    nx, ny = max(1, math.ceil(args.length / args.h)), max(1, math.ceil(args.width / args.h))
    root = ET.Element("NEKTAR")
    root.append(ET.Comment(" Generated limited-band physical pilot; not a validated DNS benchmark. "))
    geometry(root, args.length, args.width, nx, ny)
    expansions = sub(root, "EXPANSIONS")
    sub(expansions, "E", COMPOSITE="C[0]", NUMMODES=str(args.modes),
        FIELDS=",".join(VARIABLES), TYPE="MODIFIED")
    conditions = sub(root, "CONDITIONS")
    info = sub(conditions, "SOLVERINFO")
    properties = [("SolverType", "VelocityCorrectionScheme"),
                  ("EQTYPE", "UnsteadyNavierStokes"), ("EvolutionOperator", "Nonlinear"),
                  ("Projection", "Continuous"), ("AdvectionForm", "Convective"),
                  ("Extrapolation", "Standard"), ("HOMOGENEOUS", "1D"),
                  ("GlobalSysSoln", "IterativeStaticCond")]
    if args.backend == "FFTW":
        properties.append(("USEFFT", "FFTW"))
    if args.dealias:
        properties += [("DEALIASING", "True"), ("SPECTRALHPDEALIASING", "True")]
    if args.svv != "none":
        properties += [("SpectralVanishingViscositySpectralHP",
                        "DGKernel" if args.svv == "dg" else "ExpKernel"),
                       ("SpectralVanishingViscosityHomo1D", "ExpKernel")]
    for key, value in properties:
        sub(info, "I", PROPERTY=key, VALUE=value)
    scheme = sub(conditions, "TIMEINTEGRATIONSCHEME")
    sub(scheme, "METHOD", "IMEX")
    sub(scheme, "ORDER", "2")
    parameters = sub(conditions, "PARAMETERS")
    values = [("Ubar", U), ("Kinvis", NU), ("TimeStep", args.dt), ("NumSteps", steps),
              ("IO_InfoSteps", max(1, steps // 20)), ("IO_CheckSteps", args.checkpoint_every),
              ("IO_CheckTime", 0), ("FinTime", 0), ("CFL", 0), ("Flowrate", 0),
              ("HomModesZ", args.hom_modes), ("LZ", args.width),
              ("IterativeSolverTolerance", 1e-12)]
    if args.svv != "none":
        values += [("SVVDiffCoeff", args.svv_coeff), ("SVVCutoffRatio", args.svv_cutoff),
                   ("SVVDiffCoeffHomo1D", args.svv_z_coeff),
                   ("SVVCutoffRatioHomo1D", args.svv_z_cutoff)]
    for key, value in values:
        sub(parameters, "P", f"{key} = {number(value)}")
    if args.pressure_solver != "IterativeStaticCond":
        global_info = sub(conditions, "GLOBALSYSSOLNINFO")
        pressure = sub(global_info, "V", VAR="p")
        sub(pressure, "I", PROPERTY="GlobalSysSoln", VALUE=args.pressure_solver)
    variables = sub(conditions, "VARIABLES")
    for index, variable in enumerate(VARIABLES):
        sub(variables, "V", variable, ID=str(index))
    regions, boundaries = (sub(conditions, key) for key in
                           ("BOUNDARYREGIONS", "BOUNDARYCONDITIONS"))
    for region_id, composite_id in ((7, 1), (11, 2), (21, 3), (22, 4)):
        sub(regions, "B", f"C[{composite_id}]", ID=str(region_id))
    inlet = sub(boundaries, "REGION", REF="7")
    for variable, value in zip(VARIABLES[:3], ("Ubar", "0", "0")):
        sub(inlet, "D", VAR=variable, VALUE=value, USERDEFINEDTYPE="HITInflow")
    sub(inlet, "N", VAR="p", VALUE="0", USERDEFINEDTYPE="H")
    outlet = sub(boundaries, "REGION", REF="11")
    for variable in VARIABLES[:3]:
        sub(outlet, "N", VAR=variable, VALUE="0", USERDEFINEDTYPE="HOutflow")
    sub(outlet, "D", VAR="p", VALUE="0", USERDEFINEDTYPE="HOutflow")
    for region_id, partner in ((21, 22), (22, 21)):
        region = sub(boundaries, "REGION", REF=str(region_id))
        for variable in VARIABLES:
            sub(region, "P", VAR=variable, VALUE=f"[{partner}]")
    if initial:
        function = sub(conditions, "FUNCTION", NAME="InitialConditions")
        sub(function, "F", VAR=",".join(VARIABLES), FILE=str(initial))
    else:
        add_function(conditions, "InitialConditions", ["Ubar", "0", "0", "0"])
    if args.svv == "dg":
        add_function(conditions, "SVVVelocityMagnitude", ["Ubar", "0", "0"])
    hit = sub(root, "HITINFLOW", VERSION="1", UNITS="SI")
    region = sub(hit, "REGION", REF="7")
    sub(region, "ORIGIN", X="0", Y="0", Z="0")
    sub(region, "SYNTHESIS", LX=number(periods[0]), LY=number(periods[1]), LZ=number(periods[2]),
        NX=str(counts[0]), NY=str(counts[1]), NZ=str(counts[2]), SEED=str(args.seed),
        TIMEORIGIN="0", RECORDDURATION=number(synthetic_duration), REPEAT="Error")
    sub(region, "SPECTRUM", TYPE="Tabulated", FILE=SPECTRUM.name,
        INTERPOLATION="LogLinear", LOW="K4", HIGH="Gaussian", NORMALIZATION="Preserve")
    sub(region, "CUTOFF", KMAX=number(args.cutoff), TAPERSTART=number(taper))
    probes = []
    for station_index, x in enumerate(stations):
        for iy in range(probe_ny):
            for iz in range(probe_nz):
                probes.append({"id": len(probes), "station_index": station_index,
                               "iy": iy, "iz": iz,
                               "xyz": [x, iy * args.width / probe_ny, iz * args.width / probe_nz]})
    filters = sub(root, "FILTERS")
    history = sub(filters, "FILTER", TYPE="HistoryPoints")
    for key, value in (("OutputFile", "native"), ("OutputFrequency", args.sample_every),
                       ("OutputOneFile", "True"), ("WaveSpace", "False")):
        sub(history, "PARAM", value, NAME=key)
    sub(history, "PARAM", "\n" + "\n".join(" ".join(map(number, p["xyz"])) for p in probes) + "\n",
        NAME="Points")
    stresses = sub(filters, "FILTER", TYPE="ReynoldsStresses")
    for key, value in (("OutputFile", "native"), ("OutputFrequency", stress_every),
                       ("SampleFrequency", args.sample_every),
                       ("OutputStartTime", number(args.record_start)), ("MovingAverage", "False")):
        sub(stresses, "PARAM", value, NAME=key)
    ET.indent(root, space="    ")
    case = output / "pilot.xml"
    case.write_text(ET.tostring(root, encoding="unicode", xml_declaration=True) + "\n", encoding="utf-8")
    shutil.copyfile(SPECTRUM, output / SPECTRUM.name)
    knobs = vars(args).copy()
    knobs["output"] = str(output)
    knobs["initial_field"] = str(initial) if initial else None
    manifest = {
        "schema": "NEKTAR_HIT_PILOT_V1", "purpose": "Limited-band pilot input; no simulation result",
        "command": [sys.executable, str(Path(__file__).resolve())] + sys.argv[1:],
        "knobs": knobs, "case": "pilot.xml", "history_file": "native.his",
        "physics": {"mean_velocity_m_s": U, "kinvis_m2_s": NU,
                    "spectrum": "CBC1971 Table 3 theta=42 shell-3D; upstream only",
                    "normalization": "Preserve", "interpolation": "LogLinear",
                    "low_extension": "K4", "high_extension": "Gaussian"},
        "mesh": {"nx": nx, "ny": ny, "quadrilaterals": nx * ny,
                 "hx_m": args.length / nx, "hy_m": args.width / ny,
                 "expansion_modes": args.modes, "polynomial_degree": args.modes - 1,
                 "hom_modes_z": args.hom_modes, "z_grid_spacing_m": args.width / args.hom_modes,
                 "period_y_m": args.width, "period_z_m": args.width},
        "synthesis": {"periods_m": periods, "counts": counts, "seed_uint64": str(args.seed),
                      "cutoff_rad_m": args.cutoff, "taper_start_rad_m": taper,
                      "repeat_period_s": args.synthesis_length / U,
                      "strict_nyquist_rad_m": [math.pi * n / length for n, length in zip(counts, periods)],
                      "continuum_independent_mode_estimate": math.prod(periods) * args.cutoff**3 / (12 * math.pi**2),
                      "mode_estimate_note": "Sphere-volume estimate only; not the retained discrete mode count"},
        "time": {"dt_s": args.dt, "steps": steps, "initial_time_s": 0,
                 "end_s": args.duration, "record_start_s": args.record_start,
                 "declared_synthetic_record_duration_s": synthetic_duration,
                 "synthetic_record_end_margin_s": 2 * args.dt,
                 "checkpoint_every_steps": args.checkpoint_every,
                 "sample_every_steps": args.sample_every, "sample_dt_s": args.dt * args.sample_every,
                 "history_includes_initial_state": True,
                 "stress_first_sample_s": args.record_start if start_step else args.dt * args.sample_every,
                 "stress_periodic_snapshot_times_s": [step * args.dt for step in periodic_stress_steps],
                 "stress_final_snapshot_s": args.duration,
                 "stress_snapshot_note": "Cumulative post-start central moments. At positive record-start, first snapshot has one sample/zero variance. Native fixed-frequency output may add intermediate snapshots and a duplicate final snapshot."},
        "sampling": {"stations_x_m": stations, "probe_shape_yz": [probe_ny, probe_nz],
                     "probe_order": "station, then y, then z; uniform half-open periodic grids",
                     "probes": probes, "requested_outlet_buffer_m": args.buffer,
                     "station_bound_roundoff_allowance_m": station_tolerance,
                     "actual_last_station_outlet_distance_m": args.length - max(stations),
                     "estimator": "Temporal central variance at each fixed point, then arithmetic spatial average"},
        "svv": {"model_xy": args.svv, "coefficient_xy": args.svv_coeff,
                "cutoff_ratio_xy": args.svv_cutoff, "model_z": "ExpKernel" if args.svv != "none" else "none",
                "coefficient_z": args.svv_z_coeff, "cutoff_ratio_z": args.svv_z_cutoff,
                "dg_velocity_m_s": [U, 0, 0] if args.svv == "dg" else None,
                "note": "Separate xy/z solver properties; DG xy coefficient scales with native h/p and prescribed constant speed. ExpKernel xy coefficient is direct viscosity, not DG scaling."},
        "initial_conditions": {"kind": "native_field" if initial else "mean_only",
                               "path": str(initial) if initial else None,
                               "storage": ("mpi_directory" if initial.is_dir() else "single_file") if initial else None,
                               "files_sha256": initial_hashes,
                               "files_receipt_sha256": hashlib.sha256(json.dumps(initial_hashes, sort_keys=True,
                                   separators=(",", ":")).encode("utf-8")).hexdigest() if initial else None,
                               "time_note": "Designed for a time-zero initial field. Native nonzero-Time restart fields require a separately reviewed absolute-time schedule."},
        "limitations": ["No downstream experimental input or fitting", "Finite box, cutoff, taper, SVV and startup require explicit convergence/error-budget checks",
                        "8x8 default probe grid resolves statistics at sampled points; it need not resolve all transverse solver modes",
                        "Time-based checkpoints disabled; step checkpoints are optional. Native final solver field and cumulative stress fields remain enabled"],
        "sha256": {"prepare.py": digest(Path(__file__)), "pilot.xml": digest(case),
                   SPECTRUM.name: digest(output / SPECTRUM.name)},
    }
    if initial:
        for relative, checksum in initial_hashes.items():
            manifest["sha256"][str(initial / relative) if initial.is_dir() else str(initial)] = checksum
    analyzer = HERE / "analyze_history.py"
    if analyzer.exists():
        manifest["sha256"]["analyze_history.py"] = digest(analyzer)
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return case


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    for flag in ("length", "width", "h", "dt", "duration", "cutoff"):
        parser.add_argument("--" + flag, type=float, required=True)
    parser.add_argument("--modes", type=int, required=True, help="Spectral/hp expansion modes; polynomial degree=modes-1")
    parser.add_argument("--hom-modes", type=int, required=True, help="Even full Fourier packing size")
    parser.add_argument("--record-start", type=float, default=0.0)
    parser.add_argument("--sample-every", type=int, default=1)
    parser.add_argument("--checkpoint-every", type=int, default=0)
    parser.add_argument("--seed", type=int, default=20261028)
    parser.add_argument("--taper", type=float, help="Absolute taper-start wavenumber, rad/m (default 0.8*cutoff)")
    parser.add_argument("--synthesis-length", type=float, default=32.0)
    parser.add_argument("--synthesis-counts", default="auto", help="auto or NX,NY,NZ; all strict Nyquist")
    parser.add_argument("--buffer", type=float, default=0.0, help="Required distance after last station; metadata, no sponge")
    parser.add_argument("--backend", choices=("FFTW", "Matrix"), default="FFTW")
    parser.add_argument("--pressure-solver", choices=("IterativeStaticCond", "XxtStaticCond", "XxtMultiLevelStaticCond"),
                        default="IterativeStaticCond")
    parser.add_argument("--dealias", action="store_true", help="Enable both Fourier and spectral/hp dealiasing")
    parser.add_argument("--stations", required=True, help="Comma-separated probe x coordinates, m")
    parser.add_argument("--probe-grid", default="8,8", help="Half-open y,z probe counts (default 8,8)")
    parser.add_argument("--initial-field", type=Path)
    parser.add_argument("--svv", choices=("none", "dg", "exp"), default="none")
    parser.add_argument("--svv-coeff", type=float, default=1.0)
    parser.add_argument("--svv-cutoff", type=float, default=0.75)
    parser.add_argument("--svv-z-coeff", type=float, default=0.001)
    parser.add_argument("--svv-z-cutoff", type=float, default=0.75)
    args = parser.parse_args()
    try:
        case = prepare(args)
    except (ValueError, OSError) as error:
        parser.error(str(error))
    print(case)


if __name__ == "__main__":
    main()
