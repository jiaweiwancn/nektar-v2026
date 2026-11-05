#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Separate ground-cell refinement diagnostic; never an original M8 success."""
import argparse
from copy import deepcopy
import json
from pathlib import Path
import platform
import subprocess
import xml.etree.ElementTree as ET

from cls_epawf_m3_cases import set_info, set_parameter, write
from cls_epawf_m4_cases import tensor_mesh
from cls_epawf_m8_r1_cases import cost, physics, X
from cls_epawf_m8_r2_cases import RADIAL
from run_cls_epawf_picard import (prepare, run, load, sha, state_hash,
                                  write_json, require, ROOT, THREADS)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("solver", type=Path)
    parser.add_argument("resolution_probe", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--etas", default="600,400,0")
    parser.add_argument("--iterations", type=int, default=120)
    parser.add_argument("--timeout", type=int, default=900)
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    etas = [float(x) for x in args.etas.split(",")]
    require(etas and all(x >= 0 for x in etas) and
            all(a > b for a, b in zip(etas, etas[1:])), "etas must descend")
    args.output = args.output.resolve()
    if not args.prepare_only:
        require(platform.system() == "Linux" and args.output.is_relative_to(
            Path("/media/jay/新加卷/nektar-v2026-runs")), "use designated Linux storage")
        args.solver = args.solver.resolve(strict=True)
        args.resolution_probe = args.resolution_probe.resolve(strict=True)
    args.output.mkdir()
    radial = sorted(RADIAL["outer6"] + [.005])
    specification = cost(radial)
    require(specification["unknowns"] == 2315, "diagnostic tensor count changed")
    summary = dict(scope="Changed discrete problem: split [0,.01] at .005; no original M8 success implied.",
                   source=subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                                                  text=True).strip(),
                   original_unknowns=2165, changed_mesh=specification,
                   ranks=1, thread_limits=THREADS, cases=[])
    if not args.prepare_only:
        for executable in [args.solver, args.resolution_probe]:
            linked = subprocess.check_output(["ldd", str(executable)], text=True)
            require("not found" not in linked, "unresolved diagnostic executable dependency")
            libraries = {Path(s).resolve() for s in linked.split() if s.startswith("/")}
            summary[executable.name] = dict(path=str(executable), sha256=sha(executable),
                ldd=linked, libraries={str(p): sha(p) for p in sorted(libraries)})
    warm = None
    startup = None
    for eta in etas:
        case = args.output / f"eta-{eta:.10g}"
        case.mkdir()
        guess = case / "guess.json" if warm else None
        if warm:
            write_json(guess, dict(warm, alpha=400.+eta, physical_alpha=400., eta=eta))
        params = argparse.Namespace(output=case, steps=1, sweeps=args.iterations,
            verify=False, ground_projection=True, k_newton=True, equilibrate_after_k=True,
            k_radius=.4, flow=1., k=.75, epsilon=1., pseudo_alpha=0., perturbation=0., guess_file=guess)
        prepare(params)
        for name in ["startup", "rejected", "target"]:
            xml = case / name / "case.xml"
            root = ET.parse(xml).getroot()
            protected = deepcopy(physics(root))
            tensor_mesh(root, radial, 3, 2, 6., 1., points=[8, 8], xknots=X)
            require(physics(root) == protected, "mesh diagnostic changed physical inputs")
            if name == "target":
                set_info(root, "RANSWallPicardContinuation", "On")
                set_parameter(root, "RANSWallPicardContinuationEta", eta)
                set_parameter(root, "RANSWallPicardContinuationTolerance", 2e-11 if eta else 0.)
            write(root, xml)
        result = dict(eta=eta, physical_alpha=400., operator_alpha=400.+eta, time=.0025,
                      original_tolerance=2e-12, auxiliary_tolerance=2e-11 if eta else 2e-12,
                      target_xml_sha256=sha(case / "target/case.xml"), path=str(case))
        summary["cases"].append(result)
        if args.prepare_only:
            result["status"] = "prepared_only"
            continue
        if startup is None:
            receipt = run([str(args.solver), "case.xml"], case / "startup", "solver.log", args.timeout)
            startup = load(case / "startup/case.evolution.json")
            require(receipt["exit_code"] == 0 and len(startup["state"]) == 2315 and
                    startup["steps"] == startup["time"] == 0, "changed-mesh startup failed")
            summary["represented_h0_sha256"] = state_hash(startup["state"])
        target = case / "target"
        result["native"] = run([str(args.solver), "case.xml"], target, "solver.log", args.timeout)
        report_path = target / "case.evolution.json"
        if report_path.exists():
            report = load(report_path)
            require(report["steps"] == report["time"] == 0 and report["state"] == startup["state"] and
                    [state_hash(slot) for slot in report["native_history"]] ==
                    [state_hash(slot) for slot in startup["native_history"]] and
                    report["native_history_times"] == startup["native_history_times"],
                    "spatial diagnostic advanced accepted state/history")
            result.update(accepted_steps=0, accepted_time=0., history_preserved=True,
                          failure=report["failure"])
        packet = target / "picard-continuation-root.json"
        result["status"] = "no_certified_root"
        if packet.exists():
            warm = load(packet)
            require(warm["history"] == startup["state"] and warm["time"] == .0025 and
                    result.get("history_preserved"), "spatial diagnostic root changed history")
            events = [json.loads(l[20:]) for l in (target / "solver.log").read_text().splitlines()
                      if l.startswith("EP-AWF Continuation ")]
            result["root"], = [row for row in events if row["event"] == "root"]
            result["readback"] = {}
            for name, alpha in [("original", 400.), ("continued", 400.+eta)]:
                write_json(target / f"{name}-input.json", dict(warm, alpha=alpha))
                receipt = run([str(args.resolution_probe), "case.xml", f"{name}-input.json",
                               f"{name}-readback.json"], target, f"{name}-readback.log", args.timeout)
                require(receipt["exit_code"] == 0, "changed-mesh readback failed")
                norms = load(target / f"{name}-readback.json")["norms"]["native_stage"]
                result["readback"][name] = dict(receipt, norms=norms)
                if name == "continued":
                    require(max(norms) <= result["auxiliary_tolerance"], "changed-mesh root misses accuracy")
            result["status"] = "changed_discretization_root_only"
        write_json(args.output / "ground-mesh.json", summary)
    write_json(args.output / "ground-mesh.json", summary)


if __name__ == "__main__":
    main()
