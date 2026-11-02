#!/usr/bin/env python3
"""Audit the frozen M0 contract and existing Linux baseline without a new solve.

This is an executable configuration/provenance check, not a wall-patch solver.
Future components and tests in the equation matrix remain explicitly planned.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
THREADS = {name: "1" for name in
           ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS")}
os.environ.update(THREADS)
BASE = "6415398d856323c3116e646eaf7c882c66617602"
NOTE_SHA = "d35e437c2f2f8b25a46f7a3f285c6f1ad8b48d7b554aa1d161315f59eeb6de0f"
checks = 0


def require(condition, message):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(command, commands, env=None):
    start = time.monotonic()
    result = subprocess.run(command, capture_output=True, text=True, timeout=30, env=env)
    commands.append({"command": command, "exit_code": result.returncode,
                     "seconds": time.monotonic()-start})
    require(result.returncode == 0, "Command failed: "+repr(command)+"\n"+result.stderr)
    return result.stdout


def contract_check(c, root):
    require(c["schema"] == 1 and c["contract_id"] == "cls-epawf-revised-highre-original-v1",
            "Unexpected contract identity")
    require(c["base_commit"] == BASE, "Changed base revision")
    require(c["note"]["sha256"] == NOTE_SHA and sha(root/c["note"]["path"]) == NOTE_SHA,
            "Changed governing note")
    require(c["comparison_note"]["authority"] == "comparison_only_mac_source_not_part_of_linux_solver",
            "Comparison supplement became governing authority")
    model = c["model"]
    require(model["name"] == "CLSHighRe_original" and model["cubic_coefficients"] ==
            [-.1, .1, .26, "-10*C(q)^2", 0, "-5*C(q)^2", "5*C(q)^2"],
            "Original cubic coefficients changed")
    require(model["transport"] == {"sigma_k": 1, "sigma_epsilon": 1.3,
            "C_epsilon1": 1.44, "C_epsilon2": 1.92, "yap_coefficient": .83,
            "yap_length": 2.5}, "Changed original transport constants")
    for flag in ("outer_molecular_diffusion", "outer_low_re_damping", "outer_low_re_E",
                 "outer_dissipation_split", "constitutive_retuning",
                 "production_clipping", "turbulence_floors"):
        require(model[flag] is False, "Unselected model modification: "+flag)
    require(model["awf"] == {"Cmu_w": .09, "C_length_w": 2.55, "smooth_intercept": 10.7,
            "alpha": .2295, "epsilon_cutoff_star": 5.1,
            "physical_molecular_viscosity_retained": True,
            "negative_virtual_intercept_retained": True}, "Changed original AWF convention")
    fields = c["fields"]
    require(fields["outer"] == ["velocity", "pressure", "epsilon"] and
            fields["full_domain"] == ["k"] and fields["inner_stage"] == ["velocity", "pressure"],
            "Wrong multiregion field layout")
    require(fields["inner_epsilon"] == "dependent E_EP only; no independent PDE row" and
            fields["pressure_continuous_across_model_interface"] is False,
            "Dependent epsilon or discontinuous model pressure lost")
    require(c["pressure"]["variable"] == "pi=p_mechanical/rho+2*k/3" and
            c["pressure"]["flat_reference_mechanical_jump"] == "(p_w-p_o)/rho=-11*K0/75",
            "Changed pressure convention or normal-traction offset")
    fixture = c["reference_fixture"]
    require(fixture["nu"] > 0 and fixture["H"] > 0 and
            0 <= fixture["h"] < fixture["H"], "Invalid physical AWF inputs")
    for key in ("r_K", "r_m"):
        require(0 < fixture[key] <= fixture["H"], "Invalid physical "+key)
    require(fixture["H"] < fixture["r_a"] < fixture["r_b"] < fixture["Z"],
            "Independent sensors must be in the outer region")
    require(fixture["k_sampler"] == "physical_point_evaluation",
            "Changed k sampler")
    d = c["discretization"]
    require((d["velocity_degree"], d["k_degree"], d["outer_epsilon_degree"],
             d["pressure_degree"]) == (4, 4, 4, 2), "Changed starting degree pairing")
    require(d["basis"] == "Modified_A" and d["velocity_scalar_NUMMODES"] == 5 and
            d["pressure_NUMMODES"] == 3 and d["standard_points_per_direction"] == 8,
            "Degree/mode/quadrature convention changed")
    require(d["filters"] is False and d["SVV"] is False and
            d["gradient_trace_projection"] is False, "Undeclared numerical intervention")
    require(c["wall_operator"]["assembly_form"] == "simplified_net_source_once",
            "Wall source assembly may double count redistribution")
    transition = c["transition"]
    require(transition["constant_along_normal_ray"] is True and
            set(transition["state_sensor_weights"]) ==
            {"U", "DeltaU", "Normal", "K", "Epsilon", "Pressure", "Convection"} and
            all(v > 0 for v in transition["state_sensor_weights"].values()),
            "Incomplete independent sensors or normally varying chi")
    t = transition["illustrative_thresholds"]
    for lo, hi in (("d_on", "d_off"), ("m_off", "m_on"), ("R_off", "R_on")):
        require(0 < t[lo] < t[hi], "Invalid ordered gate thresholds")
    require(c["admissibility"]["cutoff_clipping"] is False, "Cutoff clipping selected")
    a = c["anchor"]
    require((a["default_eta"], a["terrain_eta"], a["flat_verification_eta"]) == (0, 0, 1) and
            a["modifies_algebraic_constraints"] is False and a["spatial_eta_allowed"] is False
            and a["reset_solution"] is False, "Invalid discrete anchor scope")
    require(c["time"]["orders"] == [1, 2] and c["time"]["reference_dt"] > 0 and
            (root/c["time"]["tableau_source"]).is_file(), "Unspecified time scheme")
    require(c["resources"]["linux_worker_cap"] == 72 and
            c["resources"]["recorded_compiler_workers"] <= 72 and
            all(c["resources"][k] == 1 for k in THREADS), "Invalid worker/thread contract")
    text = (root/c["note"]["path"]).read_text()
    labels = re.findall(r"\\label\{(eq:[^}]+)\}", text)
    mapped = [label for row in c["equation_matrix"] for label in row["equations"]]
    require(len(labels) == len(set(labels)) == 178 and
            len(mapped) == len(set(mapped)) and set(mapped) == set(labels),
            "Equation mapping has missing, duplicated or unknown labels")
    require(len({row["id"] for row in c["equation_matrix"]}) == len(c["equation_matrix"]),
            "Duplicated equation group")
    for component in c["components"].values():
        require(component["availability"] in ("existing", "planned"), "Unknown component state")
        if component["availability"] == "existing":
            require(all((root/p).is_file() for p in component["paths"]), "Missing existing source")
    for row in c["equation_matrix"]:
        require(row["component"] in c["components"] and
                all(test in c["verification"] for test in row["verification_ids"]),
                "Unmapped component or test")
    for test in c["verification"].values():
        if test["availability"] != "planned":
            require((root/test["entry_point"]).is_file(), "Missing executable verification entry")
    return {"equations": len(labels), "groups": len(c["equation_matrix"]),
            "components": len(c["components"]), "verification_ids": len(c["verification"]),
            "scope": "source-bound contract coverage, not verification of planned equations"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[4]
    contract = args.contract or root/"docs/plans/cls-epawf-revised-contract-v1.json"
    require(sys.platform.startswith("linux"), "Authoritative provenance audit requires Linux")
    args.output.mkdir(parents=True, exist_ok=False)
    report = {"schema": 1, "milestone": "cls_epawf_revised_m0", "status": "running",
              "utc": datetime.now(timezone.utc).isoformat(), "host": platform.node(),
              "cpu_affinity": sorted(os.sched_getaffinity(0)),
              "thread_environment": THREADS, "commands": [],
              "new_numerical_solve_executed": False}
    commands = report["commands"]
    try:
        c = json.loads(contract.read_text())
        report["contract_coverage"] = contract_check(c, root)
        report["contract_sha256"] = sha(contract)
        report["checker_sha256"] = sha(Path(__file__))
        report["source_commit"] = run(["git", "-C", str(root), "rev-parse", "HEAD"], commands).strip()
        status = run(["git", "-C", str(root), "status", "--porcelain"], commands)
        require(not status, "Audit producer checkout is dirty")
        previous = None
        for ref in c["reference_branches"]:
            run(["git", "-C", str(root), "cat-file", "-e", ref["commit"]+"^{commit}"], commands)
            if previous:
                run(["git", "-C", str(root), "merge-base", "--is-ancestor",
                     previous, ref["commit"]], commands)
            previous = ref["commit"]
        run(["git", "-C", str(root), "merge-base", "--is-ancestor",
             BASE, report["source_commit"]], commands)
        receipt_path = root/c["baseline"]["receipt"]
        receipt = json.loads(receipt_path.read_text())
        require(receipt["status"] == "passed" and receipt["source_commit"] == BASE and
                not receipt["source_status"] and all(x["exit_code"] == 0 for x in receipt["commands"]),
                "Baseline receipt is not an accepted exact-base run")
        require(receipt["compiler_workers"] == 16 and receipt["test_workers"] == 1 and
                receipt["mpi_ranks"] == 1 and receipt["thread_environment"] == THREADS,
                "Unexpected recorded baseline allocation")
        source, build, logs = [Path(c["baseline"][key]) for key in
                               ("linux_source", "linux_build", "linux_logs")]
        require(run(["git", "-C", str(source), "rev-parse", "HEAD"], commands).strip() == BASE and
                not run(["git", "-C", str(source), "status", "--porcelain"], commands),
                "Baseline worktree changed")
        log_hashes = {}
        for label, command in zip(("configure", "build", "ctest"), receipt["commands"]):
            p = logs/(label+".log")
            require(sha(p) == command["log_sha256"], "Baseline "+label+" log changed")
            log_hashes[label] = sha(p)
        require("100% tests passed, 0 tests failed out of 16" in (logs/"ctest.log").read_text(),
                "Missing complete baseline pass summary")
        inventory_text = run(["ctest", "--test-dir", str(build), "--show-only=json-v1"], commands)
        inventory = json.loads(inventory_text[inventory_text.index("{"):])
        require(len(inventory["tests"]) == c["baseline"]["expected_tests"] == 16,
                "Baseline test inventory changed")
        native = {}
        libs = {}
        env = dict(os.environ)
        env["LD_LIBRARY_PATH"] = str(build/"geometry-lib")+":"+env.get("LD_LIBRARY_PATH", "")
        for name in ("IncNavierStokesSolver", "RANSWallLayerProbe", "RANSKernelProbe",
                     "RANSWallGeometryProbe", "RANSObservationKernelProbe",
                     "RANSStressTangentAudit", "RANSFieldAudit", "RANSWallFieldAudit",
                     "RANSObservationAudit"):
            binary = build/name
            require(binary.is_file(), "Missing baseline binary "+name)
            native[name] = sha(binary)
            listing = run(["ldd", str(binary)], commands, env)
            require("not found" not in listing, "Unresolved baseline runtime library")
            for path in re.findall(r"(/[^\s]+)", listing):
                p = Path(path)
                if p.is_file():
                    libs[str(p.resolve())] = sha(p)
        m1_path = root/"docs/plans/cls-epawf-revised-evidence/m1/verification.json"
        m1 = json.loads(m1_path.read_text())
        require(m1["verification_status"] == "passed", "M1 kernel record is not accepted")
        for path, expected in m1["source_sha256"].items():
            require(sha(root/path) == expected, "A tested M1 source changed: "+path)
        integration_path = root/"docs/plans/cls-epawf-revised-evidence/m1/integration.json"
        integration = json.loads(integration_path.read_text())
        # Confirm shared installed/system dependencies still match M1.
        common = set(libs) & set(integration["dynamic_libraries_sha256"])
        require(bool(common) and all(libs[p] == integration["dynamic_libraries_sha256"][p]
                                     for p in common), "Shared dependencies changed since M1")
        cache = Path("/home/jay/build/nektar-v2026-mpi/CMakeCache.txt").read_text().splitlines()
        config = {}
        for key in ("NEKTAR_USE_MPI", "CMAKE_BUILD_TYPE", "CMAKE_CXX_COMPILER",
                    "CMAKE_CXX_FLAGS", "CMAKE_CXX_FLAGS_RELEASE"):
            config[key] = next((line.split("=", 1)[1] for line in cache
                                if line.startswith(key+":")), None)
        require(config["NEKTAR_USE_MPI"] == "ON" and config["CMAKE_BUILD_TYPE"] == "Release",
                "Installation is not the declared Release MPI environment")
        version = (source/"VERSION").read_text().strip()
        require(version == c["framework"]["actual_version"] == "5.9.0", "Framework version drift")
        fixture_file = source/"solvers/IncNavierStokesSolver/Tests/RANS/ABLPreservation.xml"
        xml = ET.parse(fixture_file).getroot()
        reference = m1["awf"]["note_reference"]
        f = c["reference_fixture"]
        require((f["H"], f["h"], f["r_K"], f["r_m"], f["u_star"], f["nu"]) ==
                (2, .1, 1, 1, .4, 1.5e-5), "Illustrative fixture changed")
        margin = reference["lambda"]/(1+math.log(f["H"]/f["r_m"]))-1
        reynolds = f["H"]*math.sqrt(m1["reference"]["q"]*f["u_star"]**2)/f["nu"]
        gates = c["transition"]["illustrative_thresholds"]
        require(margin >= gates["m_on"] and reynolds >= gates["R_on"],
                "Starter reference is outside correction admissibility plateaus")
        primary = Path("/media/jay/新加卷/nektar-v2026")
        primary_status = run(["git", "-C", str(primary), "status", "--porcelain"], commands)
        primary_head = run(["git", "-C", str(primary), "rev-parse", "HEAD"], commands).strip()
        report["baseline"] = {"source_commit": BASE,
            "source_tree": run(["git", "-C", str(source), "rev-parse", "HEAD^{tree}"], commands).strip(),
            "source_status": [], "receipt_sha256": sha(receipt_path),
            "recorded_started_utc": receipt["started_utc"], "test_count": 16,
            "test_names": [x["name"] for x in inventory["tests"]],
            "log_sha256": log_hashes, "binary_sha256": native, "dynamic_libraries_sha256": libs,
            "fingerprint_capture": "M0 audit of retained baseline artifacts; numerical execution remains the earlier exact-base run",
            "reference_xml_sha256": sha(fixture_file),
            "reference_expansions": [dict(x.attrib) for x in xml.findall("./EXPANSIONS/E")],
            "reference_time_scheme": {x.tag: x.text for x in xml.findall("./CONDITIONS/TIMEINTEGRATIONSCHEME/*")},
            "baseline_kernel_compile_flags": (build/"CMakeFiles/RANSKernelProbe.dir/flags.make").read_text()}
        report["framework"] = {"version": version, "configuration": config,
            "compiler": run(["c++", "--version"], commands).splitlines()[0],
            "cmake": run(["cmake", "--version"], commands).splitlines()[0],
            "mpi": run(["mpirun", "--version"], commands).splitlines()[0]}
        report["existing_m1_sources_still_match"] = True
        report["shared_dependencies_matching_m1"] = len(common)
        report["reference_gate_inputs"] = {"margin": margin, "Re_K": reynolds,
            "scope": "illustrative fixture admission only, not physical calibration or native patch validation"}
        require(run(["git", "-C", str(primary), "rev-parse", "HEAD"], commands).strip() == primary_head
                and run(["git", "-C", str(primary), "status", "--porcelain"], commands) == primary_status,
                "Unrelated primary checkout changed during the read-only audit")
        report["primary_checkout"] = {"revision": primary_head, "status": primary_status.splitlines(),
                                      "unchanged_during_audit": True}
        report["checks"] = checks
        report["status"] = "passed"
    except Exception as error:
        report["status"] = "failed"
        report["error"] = str(error)
        raise
    finally:
        (args.output/"audit.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({"status": report["status"], "coverage": report["contract_coverage"],
                      "baseline_tests_reused": 16, "new_numerical_solve_executed": False,
                      "report": str(args.output/"audit.json")}, indent=2))


if __name__ == "__main__":
    main()
