#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Prepare and run the unforced outer6 M8 stage via the native IMEX driver."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import struct
import subprocess
import time
import xml.etree.ElementTree as ET

from cls_epawf_m3_cases import set_info, set_parameter, write
from cls_epawf_m8_r2_cases import session

ROOT = Path(__file__).resolve().parents[4]
ORIGINAL = ROOT / "docs/plans/cls-epawf-revised-evidence/m8-r1/archive-session.xml"
H0 = "b8ba7917424057a9f6f16a026d82d53ede43985871537facea170ef6e1841b71"
THREADS = dict.fromkeys(("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                        "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"), "1")
PREFIX = "EP-AWF Picard "


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load(path):
    # Preserve the native sign bit of integer JSON tokens such as -0.
    return json.loads(path.read_text(), parse_int=float)


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def state_hash(state):
    return hashlib.sha256(struct.pack(f"<{len(state)}d", *state)).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def prepare(args):
    cases = {"startup": (0, 0), "rejected": (1, 0), "target": (args.steps, args.sweeps)}
    for name, (steps, sweeps) in cases.items():
        directory = args.output / name
        directory.mkdir()
        xml = session(ORIGINAL, directory / "case.xml", "outer6")
        root = ET.parse(xml).getroot()
        for key, value in {
            "RANSWallNonlinearPolicy": "SegregatedPicard",
            "RANSWallPhysicalRangeGuard": "R2LocalV1",
            "RANSWallSolveDiagnostics": "Off",
            "RANSWallPicardVerify": "On" if name == "target" and args.verify else "Off",
            "RANSWallPicardGroundProjection": "On" if name == "target" and args.ground_projection else "Off",
        }.items():
            set_info(root, key, value)
        if name == "target":
            if getattr(args, "reduced", False):
                set_info(root, "RANSWallPicardReduced", "On")
                set_info(root, "RANSWallReducedVerify", "On" if args.reduced_verify else "Off")
                for key, value in {"Basis": args.reduced_basis, "Radius": args.reduced_radius,
                                   "HTolerance": args.reduced_h_tolerance,
                                   "HUpdates": args.reduced_h_updates,
                                   "Step": args.reduced_step}.items():
                    set_parameter(root, "RANSWallReduced"+key, value)
            if getattr(args, "k_newton", False):
                set_info(root, "RANSWallPicardKNewton", "On")
            if getattr(args, "equilibrate_after_k", False):
                set_info(root, "RANSWallPicardEquilibrateK", "On")
            if getattr(args, "k_radius", None) is not None:
                set_parameter(root, "RANSWallPicardKRadius", args.k_radius)
            if getattr(args, "guess_file", None):
                set_info(root, "RANSWallPicardGuessFile", str(args.guess_file.resolve()))
            if getattr(args, "k_diagnostic", False):
                set_info(root, "RANSWallPicardKDiagnostic", "k-diagnostic.json")
        for key, value in {
            "NumSteps": steps,
            "RANSWallMaxIterations": sweeps,
            "RANSWallPicardFlowRelaxation": args.flow,
            "RANSWallPicardKRelaxation": args.k,
            "RANSWallPicardEpsilonRelaxation": args.epsilon,
            "RANSWallPicardPseudoAlpha": args.pseudo_alpha,
            "RANSWallPicardGuessPerturbation": args.perturbation if name == "target" else 0,
        }.items():
            set_parameter(root, key, value)
        write(root, xml)
    return {name: sha(args.output / name / "case.xml") for name in cases}


def run(command, directory, log, timeout):
    started = time.monotonic()
    with (directory / log).open("w") as stream:
        try:
            exit_code = subprocess.run(command, cwd=directory, env=dict(os.environ, **THREADS),
                                       stdout=stream, stderr=subprocess.STDOUT,
                                       timeout=timeout, check=False).returncode
        except subprocess.TimeoutExpired:
            exit_code = 124
    return dict(command=command, exit_code=exit_code, seconds=time.monotonic() - started,
                log_sha256=sha(directory / log))


def stages(directory):
    records = []
    with (directory / "solver.log").open() as stream:
        for line in stream:
            if line.startswith(PREFIX):
                records.append(json.loads(line[len(PREFIX):]))
    write_json(directory / "sweeps.json", records)
    return records


def readback(args, directory, state, history, time, alpha, manifest):
    packet = directory / "readback-input.json"
    write_json(packet, dict(state=state, history=history, alpha=alpha, time=time,
                           require_history_h0=False, inner_stress_model="SaturatedR2"))
    readback_report = directory / "readback.json"
    manifest["readback"] = run([str(args.resolution_probe), "case.xml", packet.name,
                               readback_report.name], directory, "readback.log", args.timeout)
    require(manifest["readback"]["exit_code"] == 0, "fresh native residual readback failed")
    norms = load(readback_report)["norms"]["native_stage"]
    manifest["readback"]["original_block_norms"] = norms
    return norms


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--solver", type=Path)
    parser.add_argument("--resolution-probe", type=Path)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--verify", action="store_true", help="Run same-state and frozen-saddle comparisons")
    parser.add_argument("--ground-projection", action="store_true", help="Project the k trial against ground cuts")
    parser.add_argument("--k-newton", action="store_true", help="Try a true conditional k correction before the frozen update")
    parser.add_argument("--reduced", action="store_true", help="Use the constrained reduced response search")
    parser.add_argument("--reduced-verify", action="store_true", help="Check repeatability, rollback and two-scale responses")
    parser.add_argument("--reduced-basis", type=int, default=10)
    parser.add_argument("--reduced-radius", type=float, default=.05)
    parser.add_argument("--reduced-h-tolerance", type=float, default=2e-9)
    parser.add_argument("--reduced-h-updates", type=int, default=50)
    parser.add_argument("--reduced-step", type=float, default=2e-5)
    parser.add_argument("--equilibrate-after-k", action="store_true",
                        help="Try fixed-k non-k equilibration when the full sweep would be rejected")
    parser.add_argument("--k-radius", type=float, help="Limit the raw k Newton proposal relative to reference k")
    parser.add_argument("--guess-file", type=Path, help="Use a saved readback packet as the nonlinear guess only")
    parser.add_argument("--k-diagnostic", action="store_true", help="Write the restricted matrices and initial action/projection checks")
    parser.add_argument("--steps", type=int, default=1)
    parser.add_argument("--sweeps", type=int, default=120)
    parser.add_argument("--timeout", type=int, default=1800)
    parser.add_argument("--flow", type=float, default=1.0)
    parser.add_argument("--k", type=float, default=0.5)
    parser.add_argument("--epsilon", type=float, default=0.5)
    parser.add_argument("--pseudo-alpha", type=float, default=0.0)
    parser.add_argument("--perturbation", type=float, default=0.0)
    args = parser.parse_args()
    args.output = args.output.resolve()
    require(args.steps >= 1 and args.sweeps >= 0 and args.timeout >= 1, "invalid run budget")
    require(not args.k_newton or args.ground_projection, "--k-newton requires --ground-projection")
    require(not args.reduced or not (args.k_newton or args.equilibrate_after_k or args.verify),
            "use --reduced-verify for reduced checks and retain R2 selectors as controls")
    require(not args.reduced_verify or args.reduced, "--reduced-verify requires --reduced")
    require(not args.equilibrate_after_k or args.k_newton, "--equilibrate-after-k requires --k-newton")
    require(not args.guess_file or args.steps == 1, "a saved guess binds exactly one stage")
    require(not args.guess_file or args.perturbation == 0, "run saved-guess and perturbation diagnostics separately")
    if not args.prepare_only:
        require(platform.system() == "Linux" and args.solver and args.resolution_probe,
                "native runs require Linux and explicit solver and readback executables")
        require(args.output.is_relative_to(Path("/media/jay/新加卷/nektar-v2026-runs")) and
                not args.output.is_relative_to(ROOT), "use a separate designated Linux run directory")
        args.solver = args.solver.resolve(strict=True)
        args.resolution_probe = args.resolution_probe.resolve(strict=True)
    args.output.mkdir(parents=True, exist_ok=False)
    manifest = dict(source=subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                                                  text=True).strip(),
                    archived_xml_sha256=sha(ORIGINAL), expected_h0_sha256=H0, ranks=1,
                    thread_limits=THREADS, xml_sha256=prepare(args), status="prepared")
    result_file = args.output / "manifest.json"
    write_json(result_file, manifest)
    if args.prepare_only:
        return 0

    try:
        for executable in (args.solver, args.resolution_probe):
            linked = subprocess.check_output(["ldd", str(executable)], text=True)
            require("not found" not in linked, f"unresolved runtime dependency for {executable}")
            libraries = {Path(word).resolve() for word in linked.split() if word.startswith("/")}
            manifest[executable.name] = dict(path=str(executable), sha256=sha(executable),
                                             ldd=linked, libraries={str(p): sha(p) for p in sorted(libraries)})
        for name in ("startup", "rejected", "target"):
            directory = args.output / name
            entry = manifest[name] = run([str(args.solver), "case.xml"], directory, "solver.log", args.timeout)
            records = stages(directory)
            receipt = directory / "case.evolution.json"
            require(receipt.is_file(), f"{name} did not publish a native report")
            report = load(receipt)
            entry.update(steps=report["steps"], time=report["time"], failure=report.get("failure"),
                         native_report_sha256=sha(receipt))
            if name == "startup":
                require(entry["exit_code"] == 0 and report["steps"] == report["time"] == 0,
                        "zero-step native startup failed")
                h0 = report["state"]
                require(len(h0) == 2165 and state_hash(h0) == H0, "historical startup does not match bitwise")
                initial_slots = [state_hash(slot) for slot in report["native_history"]]
                initial_times = report["native_history_times"]
                continue

            inputs = [row for row in records if row["event"] == "stage_input"]
            require(inputs and inputs[0]["alpha"] == 400 and inputs[0]["time"] == 0.0025 and
                    inputs[0]["unknowns"] == 2165 and inputs[0]["tolerance"] == 2e-12,
                    "native first-stage target differs")
            first = inputs[0]
            require(first["incoming_guess_sha256"] == H0, "incoming nonlinear guess differs from h0")
            require(first["history_value_sha256"] == state_hash([x if x else 0.0 for x in h0]),
                    "native callback history differs from original h0 values")
            entry["first_stage_input"] = first
            if name == "rejected":
                require(entry["exit_code"] != 0 and "Picard sweep limit" in entry["failure"] and
                        report["steps"] == report["time"] == 0 and state_hash(report["state"]) == H0 and
                        [state_hash(slot) for slot in report["native_history"]] == initial_slots and
                        report["native_history_times"] == initial_times and not report["stage_audits"],
                        "deliberate rejection changed native time, coupled state or IMEX history")
                continue

            if args.guess_file:
                require(first["guess_sha256"] == state_hash(load(args.guess_file)["state"]),
                        "saved nonlinear guess differs from the native solve input")

            if args.verify:
                checks = [row for row in records if row["event"] == "fixed_point"]
                expected = {f"{base}_{block}" for base in ("initial", "perturbed")
                            for block in ("flow", "k", "epsilon")}
                require(expected <= {row["block"] for row in checks}, "missing same-state checks")
                saddle = next(row for row in records if row["event"] == "frozen_saddle_comparison")
                require(saddle["velocity_correction"] > 0 and
                        max(saddle["original_corrected_continuity_gauge"]) <
                        max(saddle["original_predictor_continuity_gauge"]),
                        "pressure correction did not reduce native continuity and gauge defect")
            if entry["exit_code"] != 0:
                if report["steps"] == 0 and report["rejected_state"]:
                    readback(args, directory, report["rejected_state"], h0, .0025, 400., manifest)
                raise ValueError("unchanged M8 native target did not complete")
            require(report["steps"] == args.steps and
                    abs(report["time"] - args.steps * .0025) < 1e-13 and
                    len(report["stage_audits"]) == args.steps,
                    "native trajectory did not reach the requested horizon")
            require(all(row[4] <= 2e-12 for row in report["stage_audits"]) and
                    report["state_audits"][-1]["r2_qualified"] == [1.0],
                    "original stage tolerance or qualification was missed")
            norms = readback(args, directory, report["state"], report["last_stage_history"],
                             report["last_stage_time"], report["last_stage_alpha"], manifest)
            require(max(norms) <= 2e-12, "saved final state misses the original residual tolerance")
        manifest["status"] = "passed"
        return 0
    except (OSError, ValueError, KeyError, StopIteration, subprocess.CalledProcessError) as error:
        manifest.update(status="failed", failure=str(error))
        print(manifest["failure"])
        return 1
    finally:
        write_json(result_file, manifest)


if __name__ == "__main__":
    raise SystemExit(main())
