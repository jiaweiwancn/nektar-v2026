#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Follow a labelled algebraic mass continuation with the retained Picard control."""
import argparse
import json
from pathlib import Path
import platform
import subprocess
import xml.etree.ElementTree as ET

from cls_epawf_m3_cases import set_info, set_parameter, write
from run_cls_epawf_picard import prepare, run, load, sha, state_hash, write_json, require, H0, ROOT, THREADS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("solver", type=Path)
    parser.add_argument("resolution_probe", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--etas", default="20000,10000,5000,2000,1000,400,0")
    parser.add_argument("--iterations", type=int, default=100)
    parser.add_argument("--timeout", type=int, default=600)
    parser.add_argument("--bracket-iterations", type=int, default=4)
    parser.add_argument("--auxiliary-tolerance", type=float, default=2e-12,
                        help="Declared eta-positive diagnostic tolerance; eta-zero remains 2e-12")
    parser.add_argument("--initial-guess", type=Path,
                        help="Previously saved diagnostic root; physical history remains the original h0")
    args = parser.parse_args()
    require(platform.system() == "Linux", "authoritative continuation runs require Linux")
    args.solver = args.solver.resolve(strict=True)
    args.resolution_probe = args.resolution_probe.resolve(strict=True)
    args.output = args.output.resolve()
    require(args.output.is_relative_to(Path("/media/jay/新加卷/nektar-v2026-runs")),
            "use designated heavy-run storage")
    etas = [float(v) for v in args.etas.split(",")]
    require(0 < args.auxiliary_tolerance <= 1e-8, "invalid auxiliary diagnostic tolerance")
    require(etas and all(v >= 0 for v in etas) and
            all(a > b for a, b in zip(etas, etas[1:])), "etas must descend and stay nonnegative")
    args.output.mkdir()
    summary = dict(scope="Algebraic continuation only: original alpha400, h0 and t.0025 fixed; no physical acceptance.",
                   source=subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
                   solver=str(args.solver), solver_sha256=sha(args.solver),
                   resolution_probe=str(args.resolution_probe), resolution_probe_sha256=sha(args.resolution_probe),
                   ranks=1, thread_limits=THREADS, cases=[], status="running")
    summary.update(original_tolerance=2e-12, eta_positive_auxiliary_tolerance=args.auxiliary_tolerance)
    for executable in [args.solver, args.resolution_probe]:
        linked = subprocess.check_output(["ldd", str(executable)], text=True)
        require("not found" not in linked, "unresolved continuation executable dependency")
        libraries = {Path(word).resolve() for word in linked.split() if word.startswith("/")}
        summary[executable.name] = dict(ldd=linked, libraries={str(p): sha(p) for p in sorted(libraries)})
    history = None
    startup = None

    def save():
        write_json(args.output / "continuation.json", summary)

    def solve(eta, warm, label="ladder"):
        nonlocal history, startup
        case = args.output / f"eta-{eta:.10g}"
        case.mkdir()
        guess = None
        if warm:
            packet = dict(warm, alpha=400.+eta, physical_alpha=400., eta=eta)
            guess = case / "guess.json"
            write_json(guess, packet)
        params = argparse.Namespace(output=case, steps=1, sweeps=args.iterations,
                                    verify=False, ground_projection=True, k_newton=True,
                                    equilibrate_after_k=True, k_radius=.4, flow=1., k=.75,
                                    epsilon=1., pseudo_alpha=0., perturbation=0., guess_file=guess)
        prepare(params)
        if history is None:
            receipt = run([str(args.solver), "case.xml"], case / "startup", "solver.log", args.timeout)
            startup = load(case / "startup/case.evolution.json")
            history = startup["state"]
            require(receipt["exit_code"] == 0 and state_hash(history) == H0, "continuation startup differs from h0")
            summary["startup"] = receipt
        if warm:
            require(warm["history"] == history and warm["time"] == .0025 and
                    len(warm["state"]) == len(history),
                    "continuation guess would replace the original physical history or stage")
        target = case / "target"
        root = ET.parse(target / "case.xml").getroot()
        set_info(root, "RANSWallPicardContinuation", "On")
        set_parameter(root, "RANSWallPicardContinuationEta", eta)
        accuracy = args.auxiliary_tolerance if eta > 0 else 2e-12
        # Pass the request even at eta zero so the native eta-positive guard,
        # rather than only this driver, must retain the original 2e-12 gate.
        set_parameter(root, "RANSWallPicardContinuationTolerance", args.auxiliary_tolerance)
        write(root, target / "case.xml")
        result = dict(eta=eta, operator_alpha=400.+eta, physical_alpha=400., time=.0025,
                      auxiliary_tolerance=accuracy, original_tolerance=2e-12,
                      path=str(case), purpose=label, xml_sha256=sha(target / "case.xml"))
        summary["cases"].append(result)
        save()
        result["native"] = run([str(args.solver), "case.xml"], target, "solver.log", args.timeout)
        report = target / "case.evolution.json"
        if report.exists():
            state = load(report)
            require(state["steps"] == state["time"] == 0 and state_hash(state["state"]) == H0 and
                    state["native_history"] == startup["native_history"] and
                    state["native_history_times"] == startup["native_history_times"],
                    "continuation diagnostic advanced accepted state or history")
            result.update(accepted_steps=0, accepted_time=0., history_preserved=True,
                          native_failure=state["failure"])
        packet_file = target / "picard-continuation-root.json"
        events = [json.loads(line[20:]) for line in (target / "solver.log").read_text().splitlines()
                  if line.startswith("EP-AWF Continuation ")]
        roots = [row for row in events if row["event"] == "root"]
        if not packet_file.exists() or not roots:
            result["status"] = "no_certified_continuation_root"
            save()
            return None
        row, = roots
        packet = load(packet_file)
        require(result.get("history_preserved") and
                "diagnostic completed a root without advancing physical history" in result["native_failure"],
                "continued root lacks its explicit native no-advance rejection receipt")
        require(row["eta"] == eta and row["physical_alpha"] == packet["alpha"] == 400 and
                row["operator_alpha"] == 400.+eta and row["time"] == packet["time"] == .0025 and
                packet["history"] == history and max(row["continued_norms"]) <= accuracy and
                row["original_tolerance"] == 2e-12 and row["auxiliary_tolerance"] == accuracy,
                "continued root has incorrect target/history or incomplete convergence")
        result["root"] = row
        result["readback"] = {}
        for name, alpha in [("original", 400.), ("continued", 400.+eta)]:
            input_file = target / f"{name}-input.json"
            output_file = target / f"{name}-readback.json"
            write_json(input_file, dict(packet, alpha=alpha))
            receipt = run([str(args.resolution_probe), "case.xml", input_file.name, output_file.name],
                          target, name+"-readback.log", args.timeout)
            require(receipt["exit_code"] == 0, "fresh continuation residual reconstruction failed")
            norms = load(output_file)["norms"]["native_stage"]
            result["readback"][name] = dict(receipt, norms=norms)
            if name == "continued":
                require(max(norms) <= accuracy, "fresh continued residual misses the declared auxiliary tolerance")
        result["status"] = "continued_auxiliary_root_only"
        result["state_sha256"] = state_hash(packet["state"])
        print(json.dumps(dict(eta=eta, qualified=row["qualified"], range=row["outer_range"],
                              ground=row["ground"], original=max(row["original_norms"]),
                              continued=max(row["continued_norms"]))), flush=True)
        save()
        return packet

    warm = load(args.initial_guess.resolve(strict=True)) if args.initial_guess else None
    if args.initial_guess:
        summary["initial_guess"] = dict(path=str(args.initial_guess.resolve()),
                                       sha256=sha(args.initial_guess),
                                       state_sha256=state_hash(warm["state"]))
    previous = None
    bracket = None
    try:
        for eta in etas:
            root = solve(eta, warm)
            if root is None:
                summary["status"] = "continued_solve_failed_not_nonexistence"
                break
            current = summary["cases"][-1]
            if previous and previous["root"]["qualified"] and not current["root"]["qualified"] and bracket is None:
                bracket = [current, previous, root, warm]
            previous, warm = current, root
        else:
            summary["status"] = "ladder_complete_diagnostic_only"
        if bracket:
            lo, hi, lo_state, hi_state = bracket
            for _ in range(args.bracket_iterations):
                eta = (lo["eta"]+hi["eta"])/2
                root = solve(eta, hi_state, "qualification_bracket")
                if root is None:
                    break
                row = summary["cases"][-1]
                if row["root"]["qualified"]:
                    hi, hi_state = row, root
                else:
                    lo, lo_state = row, root
            summary["qualification_exit_bracket"] = dict(unqualified_eta=lo["eta"], qualified_eta=hi["eta"],
                                                         scope="Followed nearby root branch only, not global infeasibility.")
    except Exception as error:
        summary.update(status="diagnostic_error", failure=str(error))
        raise
    finally:
        save()


if __name__ == "__main__":
    main()
