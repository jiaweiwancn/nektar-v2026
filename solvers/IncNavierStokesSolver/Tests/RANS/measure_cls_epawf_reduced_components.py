#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Separate conditional y responses using original residual readbacks.

Each Y(k) is obtained by the controlled R3 H solve with zero k updates.
Partial-field readbacks are conditional diagnostics, never labelled Y(k).
The tensor mass is assembled independently; no coupled Jacobian is built.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess

from run_cls_epawf_picard import (prepare, run, load, require, sha, state_hash,
                                  write_json, THREADS, ROOT, H0)

os.environ.update(THREADS)
import numpy as np


def basis(x):
    bubble = (1.-x*x)/4.
    return np.array([(1.-x)/2., (1.+x)/2., bubble, 2.*x*bubble])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("solver", type=Path)
    parser.add_argument("resolution_probe", type=Path)
    parser.add_argument("guess", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    require(platform.system() == "Linux" and 1 <= args.workers <= 4, "use an explicit Linux worker allocation")
    args.output = args.output.resolve()
    require(args.output.is_relative_to(Path("/media/jay/新加卷/nektar-v2026-runs")), "use designated storage")
    args.solver = args.solver.resolve(strict=True)
    args.resolution_probe = args.resolution_probe.resolve(strict=True)
    args.output.mkdir()
    packet = load(args.guess)
    require(len(packet["state"]) == 2165 and packet["alpha"] == 400 and packet["time"] == .0025,
            "component diagnostic binds the unchanged original M8 stage")
    summary = dict(scope=__doc__, source=subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
        workers=args.workers, ranks_per_process=1, thread_limits=THREADS,
        guess_sha256=sha(args.guess), evaluations=[], status="running")
    for executable in [args.solver, args.resolution_probe]:
        linked = subprocess.check_output(["ldd", str(executable)], text=True)
        require("not found" not in linked, "unresolved component diagnostic dependency")
        libraries = {Path(x).resolve() for x in linked.split() if x.startswith("/")}
        summary[executable.name] = dict(path=str(executable), sha256=sha(executable),
            ldd=linked, libraries={str(p): sha(p) for p in sorted(libraries)})

    def save():
        write_json(args.output / "components.json", summary)

    def equilibrate(name, state, tolerance):
        case = args.output / name
        case.mkdir()
        guess = case / "guess.json"
        write_json(guess, dict(packet, state=state))
        params = argparse.Namespace(output=case, steps=1, sweeps=0, verify=False,
            reduced=True, reduced_verify=False, reduced_basis=4, reduced_radius=.05,
            reduced_step=2e-5, reduced_h_tolerance=tolerance, reduced_h_updates=60,
            flow=1., k=.5, epsilon=1., pseudo_alpha=0., perturbation=0.,
            ground_projection=False, guess_file=guess)
        prepare(params)
        receipt = run([str(args.solver), "case.xml"], case / "target", "solver.log", 240)
        report = load(case / "target/case.evolution.json")
        events = [json.loads(l[15:]) for l in (case / "target/solver.log").read_text().splitlines()
                  if l.startswith("EP-AWF Reduced ")]
        response, = [e for e in events if e["event"] == "non_k_response" and e["purpose"] == "base"]
        require(response["converged"] and response["attained"] <= tolerance,
                "component measurement lacks a converged H solve")
        require(report["steps"] == report["time"] == 0 and state_hash(report["state"]) == H0 and
                "reduced iteration limit" in report["failure"], "H diagnostic advanced the physical state")
        result = report["rejected_state"]
        require(result[1184:1776] == state[1184:1776], "H response modified an independent k coefficient")
        return result, dict(name=name, receipt=receipt, response=response,
                            physical_steps=0, physical_time=0., fixed_k=True)

    def probe(name, state, xml):
        directory = args.output / name
        directory.mkdir()
        shutil.copyfile(xml, directory / "case.xml")
        write_json(directory / "input.json", dict(packet, state=state, require_history_h0=False))
        receipt = run([str(args.resolution_probe), "case.xml", "input.json", "readback.json"],
                      directory, "readback.log", 180)
        require(receipt["exit_code"] == 0, "component state is outside the native evaluation domain")
        data = load(directory / "readback.json")
        require(data["effective_state"] == state, "readback silently changed an independent coefficient")
        return data, dict(name=name, receipt=receipt, input_sha256=sha(directory / "input.json"),
                           state_sha256=state_hash(state), original_norms=data["norms"]["native_stage"])

    try:
        base, receipt = equilibrate("base-H", packet["state"], 2e-11)
        summary["evaluations"].append(receipt)
        xml = args.output / "base-H/target/case.xml"
        data, receipt = probe("base", base, xml)
        summary["evaluations"].append(receipt)
        offsets = [int(x) for x in data["offsets"]]
        require(offsets[2:4] == [1184,1776], "component diagnostic field layout changed")
        field = data["fields"][2]
        fixed = {int(row) for row, _ in data["constraints"]}
        rows = [i for i in range(offsets[2], offsets[3]) if i not in fixed]
        require(len(rows) == 504, "free-k count changed")
        nodes, weights = np.polynomial.legendre.leggauss(4)
        phi = np.array([basis(x) for x in nodes])
        one = phi.T @ np.diag(weights) @ phi
        mass = np.zeros((int(field["size"]), int(field["size"])))
        volume = 0.
        for element in field["mesh"]:
            require(element["degree"] == 3, "only the original p3 tensor elements are supported")
            ids = np.array(element["ids"], dtype=int)
            signs = np.array(element["sign"])
            width = element["upper"][0]-element["lower"][0]
            height = element["upper"][2]-element["lower"][2]
            volume += width*height
            mass[np.ix_(ids,ids)] += np.kron(one,one)*width*height/4.*signs[:,None]*signs[None,:]
        local = np.array(rows)-offsets[2]
        mass = mass[np.ix_(local,local)]
        np.linalg.cholesky(mass)
        k0 = data["parameters"][5]
        scale = k0*data["ustar"]/data["parameters"][1]*np.sqrt(volume)
        r0 = np.array(data["native"]["stage_constrained"])[rows]
        def norm(v):
            return float(np.sqrt(max(0., v @ np.linalg.solve(mass,v)))/scale)
        require(abs(norm(r0)-data["norms"]["native_stage"][2]) < 1e-9,
                "independent tensor mass disagrees with the native free-k metric")
        sample = min((p for p in data["scalar_domain"]["observations"]
                      if p["inventory"] == "physical_sampled_K"),
                     key=lambda p: abs(p["coordinate"][0]+2.6))
        point = sample["coordinate"]
        element = next(e for e in field["mesh"] if all(
            e["lower"][a]-1e-12 <= point[a] <= e["upper"][a]+1e-12 for a in [0,2]))
        xi = [2.*(point[a]-element["lower"][a])/(element["upper"][a]-element["lower"][a])-1.
              for a in [0,2]]
        values = np.kron(basis(xi[1]),basis(xi[0]))*np.array(element["sign"])
        row = np.zeros(int(field["size"]))
        row[np.array(element["ids"],dtype=int)] += values
        require(abs(row @ np.array(base[offsets[2]:offsets[3]])-sample["value"]) < 1e-12,
                "independent point functional disagrees with the native sampler")
        direction = np.zeros(2165)
        direction[rows] = k0*row[local]/np.linalg.norm(row[local])
        for e in field["mesh"]:
            if e["lower"][2] == 0:
                require(all(direction[offsets[2]+int(i)] == 0. for i in e["ids"][:4]),
                        "sample direction changes the continuous ground polynomial")
        summary.update(point=point, independent_mass_native_norm_difference=norm(r0)-data["norms"]["native_stage"][2],
                       free_k_rows=504, ground_polynomial_unchanged_by_direction=True, responses=[])
        save()
        previous_actions = None
        for h, tolerance in [(1e-5,2e-11),(5e-6,2e-12)]:
            requested = (np.array(base)+h*direction).tolist()
            shifted, receipt = equilibrate(f"H-{h:g}", requested, tolerance)
            summary["evaluations"].append(receipt)
            require(state_hash([shifted[i] for i in rows]) == state_hash([requested[i] for i in rows]),
                    "response changed a free-k coefficient bit")
            states = {"reduced": shifted, "conditional_k": requested}
            for name, start, end in [("velocity",0,offsets[2]),
                                     ("pressure_gauge",offsets[4],2165),
                                     ("epsilon_trace",offsets[3],offsets[4])]:
                state = list(base)
                state[start:end] = shifted[start:end]
                states[name] = state
            def evaluate(job):
                name, state = job
                return name, probe(f"read-{h:g}-{name}",state,xml)
            with ThreadPoolExecutor(max_workers=args.workers) as pool:
                measured = dict(pool.map(evaluate,states.items()))
            actions = {}
            record = dict(h=h, H_tolerance=tolerance, actions={})
            for name, (value, receipt) in measured.items():
                summary["evaluations"].append(receipt)
                action = (np.array(value["native"]["stage_constrained"])[rows]-r0)/h
                actions[name] = action
                components = {key: norm((np.array(v)[rows]-np.array(data["terms"][key])[rows])/h)
                              for key,v in value["terms"].items() if "_k_" in key}
                record["actions"][name] = dict(norm=norm(action),
                    original_merit_slope=float(r0 @ np.linalg.solve(mass,action)/(scale*scale)),
                    k_term_action_norms=components)
            feedback = actions["reduced"]-actions["conditional_k"]
            remainder = feedback-sum(actions[n] for n in ["velocity","pressure_gauge","epsilon_trace"])
            record.update(total_feedback_norm=norm(feedback), component_sum_remainder=norm(remainder))
            if previous_actions is not None:
                summary["two_scale_relative_differences"] = {
                    name: norm(actions[name]-previous_actions[name])/max(1.,norm(actions[name]))
                    for name in actions}
            previous_actions = actions
            summary["responses"].append(record)
            save()
        summary["status"] = "measured_conditional_components_not_stationarity_or_root"
    except Exception as error:
        summary.update(status="failed_diagnostic", failure=str(error))
        raise
    finally:
        save()


if __name__ == "__main__":
    main()
