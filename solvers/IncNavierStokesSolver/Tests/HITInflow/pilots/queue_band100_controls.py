#!/usr/bin/env python3
"""Sequential DG-only/Fourier/bandwidth qualification; never launches9m.

Run on Linux in a clean HIT checkout. Each case has a distinct complete native
preflight before its full recording. Failure/rejection stops the queue without
retrying. Existing reference archives are read only, including after remounting.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

sys.dont_write_bytecode = True
import run_band100 as run

SCHEMA = "NEKTAR_HIT_BAND100_CONTROL_QUEUE_V1"
CASES = ("svv", "z", "bandwidth")
PROTOCOL = run.SOURCE / "docs/plans/hit-inflow-band100-qualification-protocol.md"


def read(path):
    return json.loads(Path(path).read_text())


def field_hashes(field):
    paths = sorted(field.rglob("*")) if field.is_dir() else [field]
    return {str(p): run.sha(p) for p in paths if p.is_file()}


def status(root):
    receipt = read(Path(root) / "queue.json")
    print(f"{run.utc()} — CBC100 qualification: {receipt['status']}")
    print("Completed controls: " + (", ".join(receipt.get("completed_cases", [])) or "none"))
    current = receipt.get("current")
    if current:
        print("Current: " + current["stage"])
        progress = Path(current["output"]) / "progress.json"
        if progress.is_file():
            run.status(progress.parent)
    if receipt.get("error"):
        print("Error: " + receipt["error"])
    print("Full-distance production: HELD; qualification review required")


class Queue:
    def __init__(self, args, bindings, allocation):
        self.args, self.bindings, self.allocation = args, bindings, allocation
        self.root = args.root.resolve()
        self.root.mkdir(parents=True)
        self.receipt = {"schema": SCHEMA, "status": "RUNNING", "started_utc": run.utc(),
            "command": sys.argv, "pid": os.getpid(), "bindings": bindings,
            "allocation": allocation, "reference_study": str(args.reference_study),
            "reference_revision": args.reference_revision, "storage_root": str(args.storage_root),
            "queue_sha256": run.sha(__file__), "protocol_sha256": run.sha(PROTOCOL),
            "cases": list(CASES), "completed_cases": [], "stages": [], "current": None,
            "full_distance": "HELD; this queue has no full-distance or automatic spatial/sampling stage"}
        self.save()

    def save(self):
        self.receipt["updated_utc"] = run.utc()
        run.atomic_json(self.root / "queue.json", self.receipt)

    def verify(self):
        run.verify_files(self.bindings)
        run.require(run.sha(__file__) == self.receipt["queue_sha256"]
                    and run.sha(PROTOCOL) == self.receipt["protocol_sha256"], "queue/protocol changed")

    def stage(self, name, output, command, timeout):
        self.verify()
        log = self.root / (name + ".log")
        run.require(not log.exists(), "existing stage log; no retry")
        row = {"stage": name, "status": "RUNNING", "output": str(output), "command": command,
               "log": str(log), "started_utc": run.utc(), "timeout_seconds": timeout}
        self.receipt["current"] = row
        self.save()
        start, process = time.monotonic(), None
        try:
            environment = dict(os.environ, **self.allocation["thread_environment"], PYTHONDONTWRITEBYTECODE="1")
            with log.open("x") as stream:
                process = subprocess.Popen(command, cwd=run.SOURCE, env=environment, stdin=subprocess.DEVNULL,
                                           stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
                row["pid"] = process.pid
                self.save()
                while process.poll() is None:
                    run.require(time.monotonic() - start <= timeout, "queue-stage hard timeout")
                    time.sleep(2)
            row["exit_code"] = process.returncode
            run.require(process.returncode == 0, "stage failed; inspect " + str(log))
            receipt = read(output / "receipt.json")
            run.require(receipt["status"] == "COMPLETED", "stage rejected/incomplete: " + receipt["status"])
            if "budget_gate" in receipt:
                run.require(receipt["budget_gate"]["passed"], "case-specific preflight budget rejected")
                row["budget_gate"] = receipt["budget_gate"]
            row.update(status="COMPLETED", receipt_sha256=run.sha(output / "receipt.json"))
            self.verify()
        except BaseException:
            row["status"] = "FAILED_OR_REJECTED"
            if process is not None and process.poll() is None:
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=20)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
            # The inner controller gives MPI its own process session. Retain
            # its PID/start-time/executable safeguards when an outer deadline
            # interrupts that controller, rather than leaving orphan ranks.
            row["owned_rank_cleanup"] = {}
            if output.exists():
                for directory in output.rglob("*.bindings"):
                    row["owned_rank_cleanup"][str(directory)] = run.cleanup_owned_ranks(directory)
            raise
        finally:
            row.update(finished_utc=run.utc(), wall_seconds=time.monotonic()-start, log_sha256=run.sha(log))
            self.receipt["stages"].append(row)
            self.receipt["current"] = None
            self.save()

    def command(self, case, stage, output, preflight=None, reference=None):
        command = [sys.executable, str(run.HERE / "run_band100.py"), "--recipe", "qualification",
            "--stage", stage, "--run-root", str(output), "--source-revision", self.args.source_revision,
            "--storage-root", str(self.args.storage_root), "--ranks", str(self.args.ranks),
            "--npz", str(self.args.npz), "--cpus", self.args.cpus, "--preflight-case", case]
        if stage == "study":
            command += ["--cases", case, "--preflight", str(preflight), "--reference-study", str(reference)]
            if case != "bandwidth":
                command += ["--compatible-reference-revision", self.args.reference_revision]
        return ["taskset", "-c", self.args.cpus] + command

    def native_work(self, name, case):
        self.verify()
        field = case / "pilot.fld"
        run.require(field.exists(), "completed final field is missing: " + str(field))
        inputs = field_hashes(field)
        run.require(inputs, "field has no partitions")
        directory = self.root / "native-work" / name
        directory.mkdir(parents=True)
        self.receipt["current"] = {"stage": "native-work-" + name, "output": str(directory)}
        self.save()
        for source in [case / "pilot.xml", *sorted(case.glob("*.dat"))]:
            shutil.copyfile(source, directory / source.name)
            inputs[str(source)] = run.sha(source)
        configuration = {"stage": "native-work", "recipe": "qualification", "allocation": self.allocation}
        controller = run.Controller(directory, self.bindings, configuration)
        try:
            controller.execute([self.bindings["executables"]["flow"]["path"], "pilot.xml", str(field)],
                               directory, "diagnostics", 3600)
            for filename, expected in inputs.items():
                run.require(run.sha(filename) == expected, "native diagnostic input changed")
            flow = read(directory / "pilot.flow.json")
            controller.finish("COMPLETED", inputs_sha256=inputs, flow=flow,
                slabs_sha256=run.sha(directory / "pilot.slabs.csv"),
                scope="One final0.6s snapshot; strain/SVV work is not temporal-central dissipation or a precision interval")
        except BaseException as error:
            controller.finish("FAILED", error=str(error))
            raise
        self.receipt.setdefault("native_work", {})[name] = {
            "receipt": str(directory / "receipt.json"), "sha256": run.sha(directory / "receipt.json")}
        self.receipt["current"] = None
        self.save()

    def execute(self):
        # Reference bytes, rather than today's mount name or controller HEAD,
        # establish compatibility with the completed economical study.
        reference = read(self.args.reference_study / "receipt.json")
        run.require(reference.get("status") == "COMPLETED", "reference study incomplete")
        proxy = type("Reference", (), {"config": {"compatible_reference_revision": self.args.reference_revision,
            "recipe": "qualification", "storage_root": str(self.args.storage_root)}, "bindings": self.bindings})()
        # The diagnostic executable is additional, not a changed CFD binary.
        proxy.bindings = dict(self.bindings, executables={k:v for k,v in self.bindings["executables"].items() if k != "flow"})
        self.receipt["reference_bridge"] = run.compatible_reference(proxy, self.args.reference_study, reference)
        self.save()
        for case in CASES:
            if case == "bandwidth":
                z = read(self.root / "z-preflight/receipt.json")
                modes = z["native_benchmark"]["mode_bytes"]
                measured = z["budget_gate"]["sampled_peak_tree_rss_gib"]
                # Planning only: use full per-rank duplication and20% extra
                # total-RSS headroom. Actual800 memory/timing is measured next.
                increment = self.args.ranks * modes * ((800 / 600)**3 - 1) / 1024**3
                estimate = 1.2 * (measured + increment)
                planning = {"schema": SCHEMA, "reference_600_mode_bytes": modes,
                    "reference_600_sampled_tree_rss_gib": measured,
                    "assumed_additional_per_rank_duplicated_mode_storage_gib": increment,
                    "planning_rss_with20percent_headroom_gib": estimate,
                    "maximum_gib": 180, "passed": estimate <= 180,
                    "scope": "Conservative planning estimate, not a rigorous peak-memory bound; fresh800 benchmark/projection/smoke still required"}
                run.atomic_json(self.root / "bandwidth-memory-planning.json", planning)
                self.receipt["bandwidth_memory_planning"] = planning
                self.save()
                run.require(planning["passed"], "bandwidth memory planning rejected; reconsider ranks before800")
            preflight = self.root / (case + "-preflight")
            self.stage(case + "-preflight", preflight, self.command(case, "preflight", preflight), 4*3600)
            study = self.root / (case + "-study")
            ref = self.root / "z-study" if case == "bandwidth" else self.args.reference_study
            self.stage(case + "-study", study, self.command(case, "study", study, preflight, ref), 20*3600)
            self.receipt["completed_cases"].append(case)
            self.save()
        for name, case in [("old-base", self.args.reference_study / "cases/base/run"),
                           ("old-x", self.args.reference_study / "cases/x/run")]:
            self.native_work(name, case)
        for case in CASES:
            self.native_work(case, self.root / (case + "-study") / "cases" / case / "run")
        self.receipt.update(status="COMPLETED_REVIEW_REQUIRED", finished_utc=run.utc())
        self.save()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--status", type=Path)
    parser.add_argument("--root", type=Path)
    parser.add_argument("--source-revision")
    parser.add_argument("--reference-study", type=Path)
    parser.add_argument("--reference-revision", default="f4b70baec10016dd098005bc5a5507b87e4a57bf")
    parser.add_argument("--storage-root", type=Path, default=Path("/media/jay/新加卷1"))
    parser.add_argument("--ranks", type=int, default=32)
    parser.add_argument("--npz", type=int, default=16)
    parser.add_argument("--cpus")
    args = parser.parse_args()
    if args.status:
        status(args.status)
        return
    run.require(sys.platform.startswith("linux"), "execution is Linux-only")
    run.require(all(v is not None for v in (args.root, args.source_revision, args.reference_study, args.cpus)), "missing execution arguments")
    args.root, args.storage_root, args.reference_study = (p.resolve() for p in (args.root, args.storage_root, args.reference_study))
    run.require(args.storage_root.is_mount() and run.within(args.root, args.storage_root)
                and not args.root.exists() and not run.within(args.root, run.SOURCE)
                and not run.within(args.root, args.reference_study), "root must be fresh, external and on the mounted data disk")
    cpus = run.parse_cpus(args.cpus)
    allocation = run.allocation(cpus, args.ranks, args.npz, os.sched_getaffinity(0))
    bindings = {"source": run.source_binding(args.source_revision), "executables": {
        key: run.executable_binding(binary) for key, binary in (
            ("solver", "IncNavierStokesSolver"), ("benchmark", "HITInflowBenchmark"),
            ("probe", "HITBoundaryProbe"), ("flow", "HITFlowDiagnostics"))}}
    os.sched_setaffinity(0, {cpus[-1]})
    queue = Queue(args, bindings, allocation)
    try:
        queue.execute()
    except BaseException as error:
        queue.receipt.update(status="FAILED_OR_REJECTED", error=f"{type(error).__name__}: {error}", finished_utc=run.utc())
        queue.save()
        raise
    finally:
        status(args.root)


if __name__ == "__main__":
    main()
