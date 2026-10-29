#!/usr/bin/env python3
"""Run a bounded, manifest-driven HIT pilot study on the Linux server."""
import argparse
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
import hashlib
import json
import os
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=3)
    parser.add_argument("--ranks", type=int, default=16)
    parser.add_argument("--npz", type=int, default=4)
    args = parser.parse_args()
    if not 1 <= args.jobs <= 3 or not 1 <= args.jobs * args.ranks <= 60:
        parser.error("require 1..3 concurrent jobs and at most 60 MPI ranks in total")
    if args.npz < 1 or args.ranks % args.npz:
        parser.error("npz must divide ranks")
    os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[:64])
    os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                      OMPI_MCA_hwloc_base_binding_policy="none",
                      PYTHONDONTWRITEBYTECODE="1")
    source = Path(__file__).resolve().parents[4]
    tools = Path(__file__).resolve().parent
    study = json.loads(args.study.read_text())
    output = args.output.resolve()
    try:
        output.relative_to(source)
    except ValueError:
        pass
    else:
        parser.error("output must be outside the source checkout")
    output.mkdir(parents=True, exist_ok=True)
    revision = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    if subprocess.check_output(["git", "-C", str(source), "status", "--porcelain"], text=True).strip():
        parser.error("source worktree must be clean")
    receipt = {"source_revision": revision, "study": str(args.study.resolve()),
               "study_sha256": hashlib.sha256(args.study.read_bytes()).hexdigest(),
               "runner_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               "concurrent_jobs": args.jobs, "ranks_per_job": args.ranks,
               "npz": args.npz, "server_cpu_limit": 64,
               "cpu_affinity": sorted(os.sched_getaffinity(0)), "cases": []}

    def execute(command, directory, name):
        with (directory / (name + ".log")).open("w") as log:
            result = subprocess.run(command, cwd=directory, stdin=subprocess.DEVNULL,
                                    stdout=log, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(f"{directory.name}/{name} exited {result.returncode}")

    def case(definition):
        name = definition["id"]
        if not name or any(c not in "abcdefghijklmnopqrstuvwxyz0123456789-_" for c in name):
            raise ValueError("case IDs must be plain lowercase directory names")
        directory = output / name
        directory.mkdir(exist_ok=False)
        knobs = dict(study["baseline"], **definition.get("override", {}))
        command = ["python3", str(tools / "pilots/prepare.py"), "--output", str(directory)]
        for key, value in knobs.items():
            if isinstance(value, bool):
                if value:
                    command.append("--" + key)
            else:
                command += ["--" + key, str(value)]
        execute(command, directory, "prepare")
        command = ["python3", str(tools / "run_m6_case.py"), "pilot.xml",
                   "--ranks", str(args.ranks), "--npz", str(args.npz),
                   "--source-revision", revision]
        if definition.get("mean_initial", False):
            command.append("--mean-initial")
        execute(command, directory, "launch")
        execute(["python3", str(tools / "pilots/analyze_history.py"),
                 "--manifest", "manifest.json", "--output", "statistics.json",
                 "--blocks", "8"], directory, "analysis")
        execute(["HITFlowDiagnostics", "pilot.xml", "pilot.fld"], directory, "diagnostics")
        run = json.loads((directory / "run.receipt.json").read_text())
        statistics = json.loads((directory / "statistics.json").read_text())
        flow = json.loads((directory / "pilot.flow.json").read_text())
        result = {"id": name, "directory": str(directory), "status": "completed",
                  "solver_wall_seconds": run["solver"]["wall_seconds"],
                  "peak_rss_kib": run["solver"]["sampled_process_tree_peak_rss_kib"],
                  "planes": [{"x": p["station_x_m"], "K": p["K_m2_s2"],
                              "variance": p["variance_m2_s2"]} for p in statistics["planes"]],
                  "flow": flow}
        print(json.dumps(result), flush=True)
        return result

    queue = iter(study["cases"])
    failed = False
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        pending = {}
        for _ in range(args.jobs):
            definition = next(queue, None)
            if definition is not None:
                pending[pool.submit(case, definition)] = definition["id"]
        while pending:
            finished, _ = wait(pending, return_when=FIRST_COMPLETED)
            for future in finished:
                name = pending.pop(future)
                try:
                    receipt["cases"].append(future.result())
                except Exception as error:
                    failed = True
                    receipt["cases"].append({"id": name, "status": "failed", "error": str(error)})
                    print(str(error), flush=True)
            if not failed:
                while len(pending) < args.jobs:
                    definition = next(queue, None)
                    if definition is None:
                        break
                    pending[pool.submit(case, definition)] = definition["id"]
            receipt["status"] = "failed" if failed else "running" if pending else "completed"
            (output / "study.receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
