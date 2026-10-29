#!/usr/bin/env python3
"""Linux acceptance of the equilibrium AWF through the production CV-to-CG adapter."""
import argparse
import csv
from datetime import datetime, timezone
from functools import lru_cache
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
THREADS = {key: "1" for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")}
os.environ.update(THREADS)
if sys.platform.startswith("linux"):
    os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[-16:])
from generate_m4 import info, parameter
from generate_m5 import generate_equilibrium, reference


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path):
    return [{key: float(value) for key, value in row.items()} for row in csv.DictReader(path.open())]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--pilot", action="store_true", help="Run only refinement pilots; cannot close the milestone")
    args = parser.parse_args()
    require(sys.platform.startswith("linux"), "Authoritative validation requires Linux")
    source = Path(__file__).resolve().parents[4]
    build = args.build.resolve(); output = args.output.resolve(); output.mkdir(parents=True, exist_ok=True)
    solver = build/"IncNavierStokesSolver"; audit = build/"RANSFieldAudit"
    binaries = {p.name: digest(p) for p in (solver, audit)}
    commands = []; cases = {}; oracle_errors = {}

    def execute(command, directory, name, rejection=None, timeout=1500):
        command = list(map(str, command)); print("Running " + name, flush=True)
        commands.append(dict(command=command, cwd=str(directory), timeout=timeout))
        result = subprocess.run(command, cwd=directory, env=os.environ.copy(), capture_output=True, text=True, timeout=timeout)
        log = result.stdout + result.stderr; (directory/(name+".log")).write_text(log)
        require((result.returncode != 0 and rejection in log) if rejection else result.returncode == 0,
                f"{name} failed: inspect {directory/(name+'.log')}")
        return log

    def run(name, ranks=1, modify=None, rejection=None, **options):
        directory = output/name; directory.mkdir(exist_ok=True)
        session = directory/(name+".xml"); generate_equilibrium(session, **options)
        tree = ET.parse(session)
        if ranks > 1: info(tree.getroot(), "GlobalSysSoln", "IterativeStaticCond")
        if modify: modify(tree.getroot())
        tree.write(session, encoding="utf-8", xml_declaration=True)
        fingerprint = dict(fixture_sha256=digest(session), binary_sha256=binaries, ranks=ranks, rejection=rejection)
        cache = directory/"case.json"
        if cache.exists() and json.loads(cache.read_text()).get("fingerprint") == fingerprint:
            result = json.loads(cache.read_text()); print("Reusing " + name, flush=True)
        else:
            command = (["mpirun", "--bind-to", "none", "--map-by", "slot", "-np", str(ranks)] if ranks > 1 else []) + [solver, session.name]
            log = execute(command, directory, name, rejection)
            result = dict(fingerprint=fingerprint, options=options, solver_command=list(map(str, command)), ranks=ranks)
            if rejection:
                result["expected_rejection"] = rejection
            else:
                require("Total Computation Time" in log, "Solver did not finish")
                raw = execute([audit, session.name, name+".fld", name+".fld"], directory, "audit")
                result["field_audit"] = json.loads(raw[raw.index('{"time"'):])
                result["stages"] = rows(directory/(name+".rans.csv"))
                result["wall_stages"] = rows(directory/(name+".awf.csv"))
                residual = directory/(name+".rans-residual.csv")
                if residual.exists():
                    result["residual_audit"] = [dict(row, **{key: float(row[key]) for key in ("time", "phase", "free_coefficient_l2", "free_coefficient_linf")}) for row in csv.DictReader(residual.open())]
                result["point_diagnostics"] = [row for path in sorted(directory.glob(name+".awf-points-r*.csv")) for row in rows(path)]
                result["artifacts_sha256"] = {p.name: digest(p) for p in directory.iterdir() if p.is_file() and p.suffix in (".xml", ".csv", ".his", ".log")}
            cache.write_text(json.dumps(result, indent=2)+"\n")
        cases[name] = result
        (output/"progress.json").write_text(json.dumps(cases, indent=2)+"\n")
        if "expected_rejection" in result: return result
        stages = result["stages"]; walls = result["wall_stages"]
        require(len(stages) == options.get("steps", 4) and len(walls) == len(stages), "Incomplete stage diagnostics")
        require(all(math.isfinite(v) for row in stages+walls+result["point_diagnostics"] for v in row.values()), "Nonfinite diagnostic")
        require(all(row["min_k"] > 0 and row["min_eps"] > 0 and row["guard_activations"] == 0 for row in stages), "Outer positivity failure")
        require(max(row["live_global_systems"] for row in stages) <= stages[0]["live_global_systems"]+2, "Growing matrix cache")
        residuals = result.get("residual_audit", [])
        require(len(residuals) == 2*len(stages)*(len(result["field_audit"]["fields"])), "Incomplete full-stage residual audit")
        require(all(math.isfinite(row[key]) for row in residuals for key in ("time", "phase", "free_coefficient_l2", "free_coefficient_linf")), "Nonfinite assembled residual")
        first_time = min(row["time"] for row in residuals)
        initial_max = max(row["free_coefficient_l2"] for row in residuals if row["phase"] == 0 and row["time"] == first_time)
        final_max = max(row["free_coefficient_l2"] for row in residuals if row["phase"] == 1)
        require(final_max < max(1e-6, .01*initial_max), "The converged full-stage residual did not decrease")
        for row in walls:
            require(row["min_kP"] > 0 and row["min_epsP"] > 0 and row["guard_activations"] == 0, "Strip positivity failure")
            require(row["epsilonH_error"] < 1e-9 and row["mass_error"] < 1e-7, "Algebraic EH or continuity failed")
            require(row["budget_available"] == 1, "Missing constant-test budget")
            require(max(abs(row["momentum_budget_x"]), abs(row["momentum_budget_y"]))/max(1., abs(row["traction_x"]), abs(row["traction_y"])) < 2e-7, "Momentum inventory/flux budget failed")
            require(abs(row["k_budget"])/max(1., 2*row["area"]*(row["mean_production"]+row["mean_dissipation"])) < 2e-7, "k inventory/source budget failed")
        check_points(result, tree.getroot())
        print(name + ": " + str({key: value["l2_error"] for key, value in result["field_audit"]["fields"].items()}), flush=True)
        return result

    def close(actual, expected, label, tolerance=2e-9):
        error = abs(actual-expected)/max(1., abs(actual), abs(expected))
        oracle_errors[label] = max(oracle_errors.get(label, 0.), error)
        require(error < tolerance, f"Independent {label} mismatch: {actual} vs {expected}")

    def check_points(case, root):
        parameters = {node.text.split("=", 1)[0].strip(): node.text.split("=", 1)[1].strip() for node in root.findall("CONDITIONS/PARAMETERS/P")}
        get = lambda key, default: float(parameters.get(key, default))
        geometry = (get("RANSStripHeight", 2), get("RANSStripPointP", 1), get("RANSRoughnessHeight", .1), get("Kinvis", 1.5e-5), get("Qstar", 3.190675081624288), get("Kappa", .440553905186659))
        high_re_limit = get("RANSAWFMinimumRe", 100)
        weights = [get("RANSAWFWeight"+suffix, 1) for suffix in ("U", "DeltaU", "K", "Epsilon", "Pressure", "Convection", "Normal")]
        on = get("RANSAWFOn", .02); off = get("RANSAWFOff", .1)
        direction = [get("RANSAWFDirectionX", 1), get("RANSAWFDirectionY", 0)]; norm = math.hypot(*direction)
        direction = [v/norm for v in direction]
        settings = {node.attrib["PROPERTY"]: node.attrib["VALUE"] for node in root.findall("CONDITIONS/SOLVERINFO/I")}
        original = settings["RANSWallModel"] == "AWFOriginal"
        @lru_cache(maxsize=None)
        def local_reference(k): return reference(k, *geometry)
        for row in case["point_diagnostics"]:
            ref = local_reference(row["kP"]); height, point_p, _, nu, _, kappa = geometry
            force = [row["forcing_x"], row["forcing_y"]]
            traction = [row[key]/ref["j0"]-ref["j1"]/ref["j0"]*f for key, f in zip(("uH", "vH"), force)]
            production = (sum(t*t for t in traction)*ref["production_weight"]+2*sum(t*f for t, f in zip(traction, force))*ref["production_weight1"]+sum(f*f for f in force)*ref["production_weight2"])/height
            for key, value in zip(("traction_x", "traction_y"), traction): close(row[key], value, key)
            close(row["production"], production, "production")
            if original:
                expected_weight = 0.
            else:
                for key, value in (("reference_u", ref["matching_velocity"]), ("reference_A", ref["numerator"]), ("reference_P0", ref["p0"]), ("lambda", ref["lam"]), ("minimum_lambda", ref["minimum_lambda"]), ("cutoff", ref["cutoff"]), ("log_z0", ref["log_z0"])):
                    close(row[key], value, key)
                us = math.hypot(ref["ustar"], nu/height); point_n = get("RANSStripPointN", 4)
                sensors = [math.hypot(*[(row[key]-ref["matching_velocity"]*d)/us for key, d in zip(("uH", "vH"), direction)]),
                           math.hypot(*[(row[n]-row[p]-ref["ustar"]/kappa*math.log(point_n/point_p)*d)/us for n, p, d in zip(("uN", "vN"), ("uP", "vP"), direction)]),
                           (row["kN"]-row["kP"])/row["kP"], (row["epsN"]-ref["numerator"]/point_n)/(ref["numerator"]/point_n),
                           height*math.hypot(row["pressure_x"], row["pressure_y"])/us**2,
                           height*math.hypot(row["convection_x"], row["convection_y"])/us**2, row["wP"]/us]
                departure = math.sqrt(sum(w*s*s for w, s in zip(weights, sensors)))
                if settings.get("RANSAWFWeight") == "Switched":
                    for key, value in zip(("r_U", "r_DeltaU", "r_K", "r_Epsilon", "r_Pressure", "r_Convection", "r_Normal"), sensors): close(row[key], value, key)
                    close(row["departure"], departure, "departure")
                    transition = (departure**2-on**2)/(off**2-on**2)
                    requested = 1. if departure <= on else (0. if departure >= off else (1-transition)**3*(1+3*transition+6*transition**2))
                else: requested = get("RANSAWFFixedWeight", 1)
                close(row["requested_chi"], requested, "requested_chi")
                wall_re = math.sqrt(row["kP"])*height/nu
                admissible = 0 < ref["cutoff"] <= point_p and math.isfinite(ref["numerator"]/ref["cutoff"])
                close(row["wall_Re"], wall_re, "wall_Re")
                close(row["admissible"], int(admissible), "admissibility")
                expected_weight = requested if admissible and wall_re >= high_re_limit else 0.
                reason = 1 if requested == 0 else (3 if wall_re < high_re_limit else (0 if admissible else 2))
                close(row["reason"], reason, "fallback_reason")
            close(row["chi"], expected_weight, "effective_chi")
            close(row["dissipation"], (1-expected_weight)*ref["epsilon_mean"]+expected_weight*ref["p0"], "mean_blend")
            close(row["epsP"], (1-expected_weight)*ref["epsilon_p"]+expected_weight*ref["numerator"]/point_p, "point_P_blend")
            close(row["epsH"], (1-expected_weight)*ref["epsilon_h"]+expected_weight*ref["numerator"]/height, "point_H_blend")

    for degree in (4, 6, 8):
        run(f"corrected_p{degree}", dt=.025, steps=4, modes=degree+1, points=degree+6)
    comparisons = {}
    if not args.pilot:
        run("corrected", dt=.025, steps=8, modes=7, points=12)
        run("switched_reference", mode="Switched", dt=.025, steps=8, modes=7, points=12)
        run("corrected_mpi", ranks=4, dt=.025, steps=8, modes=7, points=12)
        run("corrected_n5_5", ranks=4, point_n=5.5, dt=.025, steps=4, modes=7, points=12)
        run("corrected_quadrature", dt=.025, steps=4, modes=7, points=16)
        run("corrected_be_dt01", dt=.01, steps=1, order=1, modes=9, points=14)
        run("corrected_be_dt005", dt=.005, steps=1, order=1, modes=9, points=14)
        run("original", original=True, dt=.025, steps=8, modes=7, points=12)
        run("zero_weight", weight=0., dt=.025, steps=8, modes=7, points=12)
        run("partial_weight", weight=.4, dt=.025, steps=2, modes=7, points=12)
        run("corrected_molecular", diffusion="Molecular", dt=.025, steps=4, modes=7, points=12)
        run("corrected_3d", ranks=4, dimension=3, span_angle=.4, mode="Switched", dt=.001, steps=2, modes=5, points=9,
            modify=lambda r: parameter(r, "RANSPressureTolerance", 1e-7))
        run("fixed_3d", ranks=4, dimension=3, span_angle=.4, dt=.001, steps=2, modes=5, points=9,
            modify=lambda r: parameter(r, "RANSPressureTolerance", 1e-7))
        run("switched_3d_p6", ranks=4, dimension=3, span_angle=.4, mode="Switched", dt=.001, steps=2, modes=7, points=12,
            modify=lambda r: parameter(r, "RANSPressureTolerance", 1e-7))

        def epsilon_departure(root, factor):
            node = root.find("CONDITIONS/FUNCTION[@NAME='InitialConditions']/E[@VAR='eps_tilde']")
            node.attrib["VALUE"] = f"({factor})*("+node.attrib["VALUE"]+")"
            for suffix in ("U", "DeltaU", "K", "Pressure", "Convection", "Normal"): parameter(root, "RANSAWFWeight"+suffix, 0)
        run("switched_off", mode="Switched", dt=.001, steps=2, modes=7, points=12, modify=lambda r: epsilon_departure(r, 1.5))
        run("original_departure", original=True, dt=.001, steps=2, modes=7, points=12, modify=lambda r: epsilon_departure(r, 1.5))
        run("switched_transition", mode="Switched", dt=.001, steps=2, modes=7, points=12, modify=lambda r: epsilon_departure(r, 1.06))

        def body(root, expressions):
            names = [v.text for v in root.findall("CONDITIONS/VARIABLES/V")][:-1]
            function = ET.SubElement(root.find("CONDITIONS"), "FUNCTION", NAME="WallForce")
            for i, name in enumerate(names): ET.SubElement(function, "E", VAR=name, VALUE=expressions[i] if i < len(expressions) else "0")
            force = ET.SubElement(ET.SubElement(root, "FORCING"), "FORCE", TYPE="Body")
            ET.SubElement(force, "BODYFORCE").text = "WallForce"
        def pressure_sensor(root):
            body(root, ["-0.02"])
            for suffix in ("U", "DeltaU", "K", "Epsilon", "Convection", "Normal"): parameter(root, "RANSAWFWeight"+suffix, 0)
            parameter(root, "RANSAWFOn", .1); parameter(root, "RANSAWFOff", 1)
        run("pressure_sensor", mode="Switched", dt=.001, steps=2, modes=7, points=12, modify=pressure_sensor)
        def cancellation(root):
            body(root, ["(MatchingVelocity+Ustar/Kappa*log(y/2)+0.01*sin(2*PI*x/300))*0.01*2*PI/300*cos(2*PI*x/300)"])
            for suffix in ("U", "DeltaU", "K", "Epsilon", "Normal"): parameter(root, "RANSAWFWeight"+suffix, 0)
            parameter(root, "RANSAWFOn", .001); parameter(root, "RANSAWFOff", .01)
        run("forcing_cancellation", mode="Switched", perturbation=.01, dt=.001, steps=2, modes=7, points=12, modify=cancellation)
        run("switched_nonuniform", mode="Switched", perturbation=.002, dt=.01, steps=4, modes=7, points=12)
        run("cutoff_fallback", point_p=.001, failure="UseOriginal", dt=.001, steps=2, modes=7, points=12)
        run("original_cutoff", original=True, point_p=.001, dt=.001, steps=2, modes=7, points=12)
        run("regime_fallback", physical_nu=5., roughness=0., failure="UseOriginal", dt=.001, steps=2, modes=5, points=10)
        run("original_regime", original=True, physical_nu=5., roughness=0., dt=.001, steps=2, modes=5, points=10)
        def configured(root):
            parameter(root, "RANSYapCoefficient", 0)
            parameter(root, "Kappa", math.sqrt(.48*1.3/3.190675081624288))
        run("configured_reference", dt=.025, steps=4, modes=7, points=12, modify=configured)
        run("reject_cutoff", ranks=4, point_p=.001, dt=.001, steps=1, modes=5, points=10, rejection="fails cutoff/high-Re preflight")
        run("reject_regime", ranks=4, physical_nu=5., roughness=0., dt=.001, steps=1, modes=5, points=10, rejection="fails cutoff/high-Re preflight")
        run("reject_thresholds", ranks=4, mode="Switched", on=.2, off=.1, dt=.001, steps=1, modes=5, points=10, rejection="off > on")
        def no_sensors(root):
            for suffix in ("U", "DeltaU", "K", "Epsilon", "Pressure", "Convection", "Normal"): parameter(root, "RANSAWFWeight"+suffix, 0)
        run("reject_sensors", ranks=4, mode="Switched", dt=.001, steps=1, modes=5, points=10, modify=no_sensors, rejection="active sensors and a fixed incoming direction")
        run("reject_z0", ranks=4, dt=.001, steps=1, modes=5, points=10, modify=lambda r: parameter(r, "RANSZ0", .01), rejection="do not prescribe RANSZ0")
        run("reject_weight", ranks=4, weight=1.1, dt=.001, steps=1, modes=5, points=10, rejection="Invalid AWF weight or high-Re threshold")

        for first, second in (("corrected", "corrected_mpi"), ("corrected", "switched_reference"), ("original", "zero_weight"), ("switched_off", "original_departure"), ("cutoff_fallback", "original_cutoff"), ("regime_fallback", "original_regime"), ("corrected_p6", "corrected_quadrature")):
            raw = execute([audit, output/first/(first+".xml"), output/first/(first+".fld"), output/second/(second+".fld")], output/second, "comparison")
            comparison = json.loads(raw[raw.index('{"time"'):]); comparisons[first+"__"+second] = comparison
            require(max(v["difference_linf"] for v in comparison["fields"].values()) < 2e-7, "Complete-field recovery/MPI/quadrature comparison failed")

    checks = {}
    for field in ("u", "v", "k", "eps_tilde", "p"):
        errors = [cases[f"corrected_p{p}"]["field_audit"]["fields"][field]["l2_error"] for p in (4, 6, 8)]
        require(errors[-1] < errors[0], f"Corrected reference does not refine: {field}")
        checks[field+"_refinement"] = errors
    initial_residuals = {name: {row["field"]: row["free_coefficient_l2"] for row in case.get("residual_audit", []) if row["phase"] == 0 and row["time"] == .025} for name, case in cases.items() if name.startswith("corrected_p")}
    for field in ("u", "k", "eps_tilde"):
        require(initial_residuals["corrected_p8"][field] < initial_residuals["corrected_p4"][field], "Assembled initial residual does not refine: "+field)
    checks["initial_assembled_stage_residuals"] = initial_residuals
    if not args.pilot:
        for name in ("corrected", "switched_reference"):
            require(all(row["chi_min"] == row["chi_max"] == 1 for row in cases[name]["wall_stages"]), "Reference switch left its plateau")
        require(all(row["chi_min"] == row["chi_max"] == 1 for row in cases["fixed_3d"]["wall_stages"]), "Fixed 3D correction was not active")
        for name in ("corrected_3d", "fixed_3d", "switched_3d_p6"):
            require(all(abs(row["traction_y"]/row["traction_x"]-math.tan(.4)) < 1e-5 for row in cases[name]["wall_stages"]), "Rotated 3D traction lost its direction")
        require(all(row["chi_min"] == row["chi_max"] == 1 for row in cases["switched_3d_p6"]["wall_stages"]), "Refined 3D reference did not reach the switch plateau")
        coarse_departure = max(row["departure"] for row in cases["corrected_3d"]["point_diagnostics"] if row["phase"] == 0 and row["time"] == .001)
        fine_departure = max(row["departure"] for row in cases["switched_3d_p6"]["point_diagnostics"] if row["phase"] == 0 and row["time"] == .001)
        require(fine_departure < coarse_departure, "3D reference sensor does not refine")
        checks["initial_3d_departure"] = dict(degree4=coarse_departure, degree6=fine_departure)
        initial_norms = {}
        for name in ("corrected", "corrected_mpi"):
            residuals = cases[name]["residual_audit"]; first_time = min(row["time"] for row in residuals)
            initial_norms[name] = {row["field"]: row["free_coefficient_l2"] for row in residuals if row["phase"] == 0 and row["time"] == first_time}
        for field in initial_norms["corrected"]:
            close(initial_norms["corrected"][field], initial_norms["corrected_mpi"][field], "assembled_serial_mpi_"+field, 2e-8)
        checks["initial_assembled_serial_mpi"] = initial_norms
        for name in ("zero_weight", "switched_off", "cutoff_fallback", "regime_fallback"):
            require(all(row["chi_max"] == 0 for row in cases[name]["wall_stages"]), "Baseline fallback is not exact")
        require(all(row["cutoff_fallback_points"] > 0 for row in cases["cutoff_fallback"]["wall_stages"]), "Cutoff fallback not reported")
        require(all(row["regime_fallback_points"] > 0 for row in cases["regime_fallback"]["wall_stages"]), "Wall regime fallback not reported")
        require(any(0 < row["chi"] < 1 for row in cases["switched_transition"]["point_diagnostics"]), "Switch transition not exercised")
        initial = [row for row in cases["forcing_cancellation"]["point_diagnostics"] if row["phase"] == 0 and row["time"] == .001]
        require(any(row["r_Pressure"] > .01 and row["r_Convection"] > .01 and math.hypot(row["forcing_x"], row["forcing_y"]) < .02*math.hypot(row["convection_x"], row["convection_y"]) and row["chi"] == 0 for row in initial), "Cancelling forcing was not sensed independently")
        require(max(row["normal_velocity_max"] for row in cases["switched_nonuniform"]["wall_stages"]) > 1e-6, "Nonuniform normal coupling not exercised")
        require(cases["corrected_n5_5"]["wall_stages"][-1]["remote_N_without_wall_samples"] > 0, "Physical nonvertex N did not exercise wall-less rank ownership")
        for field in ("k", "eps_tilde"):
            ratio = cases["original"]["field_audit"]["fields"][field]["l2_error"]/cases["corrected"]["field_audit"]["fields"][field]["l2_error"]
            checks[field+"_original_to_corrected_error_ratio"] = ratio
            require(ratio > 100, "Correction did not remove the dominant original mismatch")

    source_files = ["solvers/IncNavierStokesSolver/EquationSystems/"+name for name in ("RANSVelocityCorrectionScheme.cpp", "RANSVelocityCorrectionScheme.h", "RANSAnalyticalWallStrip.cpp", "RANSAnalyticalWallStrip.h")]
    source_files += ["solvers/IncNavierStokesSolver/Tests/RANS/"+name for name in ("generate_m3.py", "generate_m4.py", "generate_m5.py", "validate_m5.py", "plot_m5.py", "CMakeLists.txt", "AWFEquilibrium.xml")]
    source_files += ["solvers/IncNavierStokesSolver/CMakeLists.txt", "library/UnitTests/SolverUtils/RANSFieldAudit.cpp", "solvers/IncNavierStokesSolver/Tests/TGV3d_2_2.xml"]
    source_files += ["library/SolverUtils/Turbulence/"+name for name in ("RANSModels.cpp", "RANSModels.h", "AnalyticalWallFunction.cpp", "AnalyticalWallFunction.h", "RANSWallStrip.h")]
    # Keep full per-point receipts and native output under the Linux run root.
    # Commit compact summaries with exact hashes of every raw point CSV.
    report_cases = {}
    for name, case in cases.items():
        saved = {key: value for key, value in case.items() if key != "point_diagnostics"}
        points = case.get("point_diagnostics", [])
        if points:
            selected = [points[0], max(points, key=lambda row: row["departure"])]
            transition = next((row for row in points if 0 < row["chi"] < 1), None)
            if transition: selected.append(transition)
            saved["point_summary"] = dict(rows=len(points), ranges={key: [min(row[key] for row in points), max(row[key] for row in points)] for key in points[0]}, examples=selected)
        report_cases[name] = saved
    report = dict(schema=1, milestone="neutral_abl_rans_awf_m5", status="pilot" if args.pilot else "passed", utc=datetime.now(timezone.utc).isoformat(),
                  source_parent=subprocess.check_output(["git", "-C", source, "rev-parse", "HEAD"], text=True).strip(),
                  source_sha256={name: digest(source/name) for name in source_files}, binary_sha256=binaries, python=platform.python_version(),
                  maximum_project_workers=16, cpu_affinity=sorted(os.sched_getaffinity(0)), maximum_mpi_ranks=4, concurrent_cases=1, threads=THREADS,
                  checks=checks, oracle_maximum_scaled_errors=oracle_errors, full_field_comparisons=comparisons, cases=report_cases, commands=commands)
    (output/"validation.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(dict(status=report["status"], case_count=len(cases), checks=checks, oracle_maximum_scaled_errors=oracle_errors), indent=2), flush=True)


if __name__ == "__main__":
    main()
