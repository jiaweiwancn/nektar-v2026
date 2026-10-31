#!/usr/bin/env python3
"""Run fresh, bounded CBC100 pilot jobs on the authoritative Linux host.

The preflight and study are separate invocations. Every MPI rank is pinned to
one explicitly allocated CPU. Commands, running/failed/completed receipts and
native logs are retained outside the checkout; no failed command is retried.
This controller never launches a full-distance production experiment.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import time
import xml.etree.ElementTree as ET


HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[4]
SCHEMA = "NEKTAR_HIT_BAND100_RUN_V1"
CASE_ORDER = ("base", "dt", "z", "bandwidth", "svv", "h")
U, WIDTH = 12.7, .6
SAMPLE_DT, DURATION, RECORD_START = .00025, .6, .2
COUNTS, SEED = "2548,154,154", "20261028"
FLOAT = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
STEP = re.compile(r"\bSteps:\s*(\d+)\s+Time:\s*(" + FLOAT + r").*?CPU Time:\s*(" + FLOAT + r")s")
READY = re.compile(r"HIT_NATIVE_READY\s+rank\s+(\d+)\b")
EVALUATOR = re.compile(r"HIT_EVALUATOR\s+rank\s+(\d+)\s+points\s+(\d+)\s+method\s+(\S+)")
BAD_LOG = re.compile(r"Exceeded maxIt|failed to converge|(?<![A-Za-z])[+-]?nan(?=$|[\s,;\]])|[:=]\s*[+-]?inf(?:inity)?(?=$|[\s,;\]])", re.I)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def within(path, root):
    try:
        Path(path).resolve().relative_to(Path(root).resolve())
        return True
    except ValueError:
        return False


def parse_cpus(text):
    try:
        cpus = [int(v) for v in text.split(",")]
    except ValueError as error:
        raise ValueError("--cpus must be a comma-separated list of CPU integers") from error
    require(cpus and len(cpus) <= 64 and len(cpus) == len(set(cpus)) and all(v >= 0 for v in cpus),
            "--cpus needs 1..64 unique nonnegative CPU IDs")
    return cpus


def allocation(cpus, ranks, npz, allowed):
    require(1 <= ranks <= 48 and ranks <= len(cpus), "require 1..48 MPI ranks and at least one assigned CPU per rank")
    require(npz > 0 and ranks % npz == 0 and 128 % npz == 0 and 192 % npz == 0,
            "npz must divide ranks and both 128/192 Fourier slot counts")
    require(set(cpus) <= set(allowed), "allocated CPUs are outside the calling process affinity")
    return {"cpus": cpus, "ranks": ranks, "npz": npz, "hard_server_worker_limit": 64,
            "thread_environment": {k: "1" for k in (
                "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
                "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")},
            "concurrency": "One native MPI command at a time; controller shares the explicitly assigned CPU pool",
            "external_jobs": "The caller must account for other server jobs within the same64-worker ceiling"}


def case_knobs(name):
    require(name in CASE_ORDER, "unknown case")
    knobs = {"h": .05, "nz": 128, "dt": 12.5e-6, "cutoff": 600, "taper": 480, "svv": 1.0}
    if name == "dt":
        knobs["dt"] /= 2
    elif name == "z":
        knobs["nz"] = 192
    elif name == "bandwidth":
        knobs.update(nz=192, cutoff=800, taper=640)
    elif name == "svv":
        knobs["svv"] = .5
    elif name == "h":
        knobs["h"] = .025
    return knobs


def source_binding(revision):
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=SOURCE, text=True).strip()
    require(re.fullmatch(r"[0-9a-f]{40}", revision) and head == revision,
            "--source-revision must be the exact full40-character HEAD")
    branch = subprocess.check_output(["git", "branch", "--show-current"], cwd=SOURCE, text=True).strip()
    require(branch == "codex/hit-inflow-generator", "the source must remain on codex/hit-inflow-generator")
    status = subprocess.check_output(["git", "status", "--porcelain", "--untracked-files=normal"], cwd=SOURCE, text=True)
    require(not status.strip(), "the Linux checkout must be clean before a stage starts")
    names = [HERE / n for n in ("run_band100.py", "prepare.py", "analyze_band100.py", "analyze_m7.py", "analyze_history.py")]
    names += [HERE.parent / "cbc1971_theta42_spectrum_3d.dat", HERE.parent / "cbc1971-source.json"]
    names += [SOURCE / "docs/plans/hit-inflow-band100-protocol.md"]
    return {"revision": head, "branch": branch, "checkout": str(SOURCE),
            "files_sha256": {str(p.relative_to(SOURCE)): sha(p) for p in names}}


def executable_binding(value):
    found = shutil.which(value)
    require(found is not None, f"native executable not found: {value}; activate the Nektar environment")
    path = Path(found).resolve()
    dependencies = {}
    output = subprocess.check_output(["ldd", str(path)], text=True, timeout=60)
    require("not found" not in output, f"a native dependency is missing for {path}")
    for line in output.splitlines():
        match = re.search(r"(?:=>\s*)?(/\S+)", line)
        if match:
            library = Path(match.group(1)).resolve()
            if library.is_file() and ("nektar" in str(library).lower() or "Nektar" in library.name):
                dependencies[str(library)] = sha(library)
    return {"path": str(path), "sha256": sha(path), "bytes": path.stat().st_size,
            "native_library_sha256": dependencies}


def verify_files(bindings):
    for filename, expected in bindings["source"]["files_sha256"].items():
        require(sha(SOURCE / filename) == expected, f"source changed during execution: {filename}")
    for item in bindings["executables"].values():
        require(sha(item["path"]) == item["sha256"], "native executable changed during execution")
        for filename, expected in item["native_library_sha256"].items():
            require(sha(filename) == expected, "native library changed during execution")


def rank_entry(arguments):
    parser = argparse.ArgumentParser(description="Internal singleton-affinity MPI rank launcher")
    parser.add_argument("--cpus", required=True)
    parser.add_argument("--binding-dir", type=Path, required=True)
    parser.add_argument("--binary", required=True)
    parser.add_argument("arguments", nargs=argparse.REMAINDER)
    args = parser.parse_args(arguments)
    cpus = parse_cpus(args.cpus)
    require("OMPI_COMM_WORLD_LOCAL_RANK" in os.environ, "rank launcher requires OpenMPI local rank metadata")
    rank = int(os.environ["OMPI_COMM_WORLD_LOCAL_RANK"])
    require(0 <= rank < len(cpus), "MPI rank has no CPU allocation")
    os.sched_setaffinity(0, {cpus[rank]})
    selected = sorted(os.sched_getaffinity(0))
    require(selected == [cpus[rank]], "failed to establish singleton rank affinity")
    record = {"rank": rank, "pid": os.getpid(), "cpu_affinity": selected,
              "thread_environment": {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
              "binary": str(Path(args.binary).resolve()), "created_utc": utc()}
    record["process_start_ticks"] = process_start_ticks(os.getpid())
    target = args.binding_dir / f"rank-{rank:03d}.json"
    require(not target.exists(), "rank binding already exists; refusing duplicate rank launch")
    temporary = target.with_suffix(target.suffix + ".tmp")
    with temporary.open("x", encoding="utf-8") as stream:
        json.dump(record, stream, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, target)
    native_args = args.arguments[1:] if args.arguments[:1] == ["--"] else args.arguments
    os.execv(args.binary, [args.binary] + native_args)


def log_progress(text):
    matches = list(STEP.finditer(text))
    if not matches:
        return None
    match = matches[-1]
    return {"step": int(match.group(1)), "simulated_time_s": float(match.group(2)),
            "last_info_chunk_seconds": float(match.group(3))}


def tail(path, count=128 * 1024):
    if not path.exists():
        return ""
    with path.open("rb") as stream:
        stream.seek(max(0, path.stat().st_size - count))
        return stream.read().decode("utf-8", errors="replace")


def resources(pid):
    rows = subprocess.check_output(["ps", "-eo", "pid=,ppid=,rss=,nlwp="], text=True, timeout=10).splitlines()
    table = [tuple(map(int, row.split())) for row in rows]
    descendants = {pid}
    while True:
        added = {child for child, parent, _, _ in table if parent in descendants}
        if added <= descendants:
            break
        descendants |= added
    selected = [row for row in table if row[0] in descendants]
    return sum(row[2] for row in selected), sum(row[3] for row in selected)


def process_start_ticks(pid):
    # After the final ')' fields start at proc-stat field3; starttime is field22.
    return Path(f"/proc/{pid}/stat").read_text().rsplit(")", 1)[1].split()[19]


def task_ids(pid):
    return sorted(int(path.name) for path in Path(f"/proc/{pid}/task").iterdir())


def task_start_ticks(pid, tid):
    return Path(f"/proc/{pid}/task/{tid}/stat").read_text().rsplit(")", 1)[1].split()[19]


class AffinityAudit:
    """Sample native initialization within the pool, then strict rank pinning.

    Some native initialization temporarily changes affinity. This is permitted
    only within the explicitly allocated pool and before the corresponding
    rank's HIT_NATIVE_READY marker. No cause is inferred from that observation.
    """
    def __init__(self, cpus, ranks):
        self.cpus, self.ranks = list(cpus), ranks
        self.ready_ranks = set()
        self.process_checks = self.thread_checks = 0
        self.startup_thread_checks = self.strict_thread_checks = 0
        self.main_thread_checks = self.helper_thread_checks = 0
        self.exited_process_observations = self.disappeared_thread_observations = 0
        self.startup_nonassigned_observations = self.startup_mask_changes = 0
        self.strict_checks_by_rank = [0] * ranks
        self.witness_limit, self.last_mask_limit = 64, 512
        self.witnesses, self.last_masks = [], {}
        self.omitted_witnesses = 0

    def observe_log(self, text, completed=False):
        for match in READY.finditer(text):
            rank = int(match.group(1))
            require(0 <= rank < self.ranks, f"unexpected native readiness rank {rank}")
            self.ready_ranks.add(rank)
        if STEP.search(text) or completed:
            missing = sorted(set(range(self.ranks)) - self.ready_ranks)
            require(not missing, f"native readiness missing for ranks {missing} before " +
                    ("successful completion" if completed else "the first observed solver step"))

    def check_record(self, record):
        rank, pid = record["rank"], record["pid"]
        expected, pool = {self.cpus[rank]}, set(self.cpus)
        try:
            observed_start = process_start_ticks(pid)
            require(observed_start == record["process_start_ticks"],
                    f"MPI PID identity changed for rank {rank}, PID {pid}: "
                    f"observed start {observed_start}, expected {record['process_start_ticks']}")
            threads = task_ids(pid)
        except (FileNotFoundError, ProcessLookupError):
            self.exited_process_observations += 1
            return
        self.process_checks += 1
        for tid in threads:
            try:
                first_start = task_start_ticks(pid, tid)
                mask = set(os.sched_getaffinity(tid))
                # An exiting/reused TID must still belong to this process and
                # have the same start time before its mask is interpreted.
                final_start = task_start_ticks(pid, tid)
            except (FileNotFoundError, ProcessLookupError):
                self.disappeared_thread_observations += 1
                continue
            require(first_start == final_start,
                    f"MPI thread identity changed for rank {rank}, PID {pid}, TID {tid}")
            self.thread_checks += 1
            strict = rank in self.ready_ranks
            if strict:
                self.strict_thread_checks += 1
                self.strict_checks_by_rank[rank] += 1
            else:
                self.startup_thread_checks += 1
                self.startup_nonassigned_observations += mask != expected
            self.main_thread_checks += tid == pid
            self.helper_thread_checks += tid != pid
            phase = "ready" if strict else "native_initialization"
            key, current = (rank, tid, first_start), tuple(sorted(mask))
            previous = self.last_masks.get(key)
            changed = previous is not None and previous != current
            if changed and not strict:
                self.startup_mask_changes += 1
            if changed or (previous is None and mask != expected):
                witness = {"rank": rank, "pid": pid, "tid": tid, "phase": phase,
                           "previous_mask": list(previous) if previous is not None else None,
                           "observed_mask": list(current), "expected_mask": sorted(expected),
                           "observed_utc": utc()}
                if len(self.witnesses) < self.witness_limit:
                    self.witnesses.append(witness)
                else:
                    self.omitted_witnesses += 1
            if key in self.last_masks or len(self.last_masks) < self.last_mask_limit:
                self.last_masks[key] = current
            require(mask and mask <= pool,
                    f"MPI affinity escaped CPU pool for rank {rank}, PID {pid}, TID {tid}, "
                    f"phase {phase}: observed mask {sorted(mask)}, expected mask {sorted(expected)}, "
                    f"allowed pool {self.cpus}")
            if strict:
                require(mask == expected,
                        f"ready MPI affinity changed for rank {rank}, PID {pid}, TID {tid}: "
                        f"observed mask {sorted(mask)}, expected mask {sorted(expected)}")

    def snapshot(self):
        return {"scope": "Sampled main/helper thread affinities; not a continuous scheduler audit",
                "policy": "CPU-pool containment in every phase; exact assigned singleton after each rank's HIT_NATIVE_READY marker",
                "native_initialization_cause": "Not established; no hwloc root-cause claim",
                "ready_ranks": sorted(self.ready_ranks), "expected_rank_count": self.ranks,
                "process_identity_checks": self.process_checks, "thread_affinity_checks": self.thread_checks,
                "main_thread_checks": self.main_thread_checks, "helper_thread_checks": self.helper_thread_checks,
                "startup_thread_checks": self.startup_thread_checks, "strict_thread_checks": self.strict_thread_checks,
                "strict_thread_checks_by_rank": self.strict_checks_by_rank,
                "startup_nonassigned_thread_observations": self.startup_nonassigned_observations,
                "startup_mask_changes": self.startup_mask_changes,
                "exited_process_observations": self.exited_process_observations,
                "disappeared_thread_observations": self.disappeared_thread_observations,
                "affinity_change_witnesses": self.witnesses, "witness_limit": self.witness_limit,
                "omitted_change_witnesses": self.omitted_witnesses}


def owned_native_alive(record):
    try:
        if process_start_ticks(record["pid"]) != record["process_start_ticks"]:
            return False
        actual = os.readlink(f"/proc/{record['pid']}/exe")
        actual = actual.removesuffix(" (deleted)")
        return Path(actual).resolve() == Path(record["binary"]).resolve()
    except (FileNotFoundError, ProcessLookupError):
        return False


def cleanup_owned_ranks(directory):
    """Stop only recorded, identity-checked native ranks after a command failure.

    OpenMPI can put ranks in process groups different from mpirun's. PID start
    time and /proc/exe checks prevent signalling a reused or unrelated PID.
    """
    records = [json.loads(p.read_text()) for p in sorted(directory.glob("rank-*.json"))]
    actions = []
    for record in records:
        if owned_native_alive(record):
            try:
                os.kill(record["pid"], signal.SIGTERM)
                actions.append({"pid": record["pid"], "signal": "SIGTERM"})
            except ProcessLookupError:
                pass
    deadline = time.monotonic() + 10
    while any(owned_native_alive(r) for r in records) and time.monotonic() < deadline:
        time.sleep(.25)
    for record in records:
        if owned_native_alive(record):
            try:
                os.kill(record["pid"], signal.SIGKILL)
                actions.append({"pid": record["pid"], "signal": "SIGKILL"})
            except ProcessLookupError:
                pass
    deadline = time.monotonic() + 5
    while any(owned_native_alive(r) for r in records) and time.monotonic() < deadline:
        time.sleep(.25)
    require(not any(owned_native_alive(r) for r in records), "owned native MPI ranks survived failure cleanup")
    return actions


def validate_bindings(directory, cpus, ranks, live=False, audit=None):
    if live and audit is None:
        # Standalone callers retain strict validation. Only a command's explicit
        # readiness audit enables the bounded native-initialization phase.
        audit = AffinityAudit(cpus, ranks)
        audit.ready_ranks.update(range(ranks))
    records = []
    for path in sorted(directory.glob("rank-*.json")):
        value = json.loads(path.read_text())
        rank = value["rank"]
        require(0 <= rank < ranks and value["cpu_affinity"] == [cpus[rank]], "MPI binding receipt mismatch")
        require(all(v == "1" for v in value["thread_environment"].values()), "MPI threaded library allocation is not1")
        if live:
            audit.check_record(value)
        records.append(value)
    if not live:
        require([v["rank"] for v in records] == list(range(ranks)), "MPI binding receipts are incomplete")
    return records


def evaluator_reports(text, ranks):
    """Validate reports from local inlet owners, independently of MPI readiness.

    Ranks without a local inlet need no evaluator object and emit no report.
    Their native initialization is instead witnessed by HIT_NATIVE_READY.
    """
    reports = [{"rank": int(m.group(1)), "points": int(m.group(2)), "method": m.group(3)}
               for m in EVALUATOR.finditer(text)]
    require(reports, "native accelerated-evaluator reports are missing")
    seen = [r["rank"] for r in reports]
    require(len(seen) == len(set(seen)) and all(0 <= rank < ranks for rank in seen),
            "native evaluator reports have duplicate or out-of-range ranks")
    methods = {"GroupedTransverse", "GroupedSeparable", "EmptyPlane", "MeanOnly"}
    require(all(r["method"] in methods for r in reports),
            "native evaluator report shows unavailable acceleration/direct fallback")
    require(any(r["points"] > 0 and r["method"] in {"GroupedTransverse", "GroupedSeparable"}
                for r in reports), "no nonempty accelerated inlet evaluator was reported")
    return sorted(reports, key=lambda r:r["rank"])


class Controller:
    def __init__(self, root, bindings, config):
        self.root, self.bindings, self.config = Path(root), bindings, config
        self.progress = {"schema": SCHEMA, "stage": config["stage"], "status": "RUNNING",
                         "started_utc": utc(), "source_revision": bindings["source"]["revision"],
                         "current": None, "completed_cases": [], "commands": []}
        self.receipt = {"schema": SCHEMA, "status": "RUNNING", "started_utc": utc(),
                        "bindings": bindings, "configuration": config, "commands": [], "cases": {}}
        atomic_json(self.root / "progress.json", self.progress)
        atomic_json(self.root / "receipt.json", self.receipt)

    def update(self, **values):
        self.progress.update(values, updated_utc=utc())
        atomic_json(self.root / "progress.json", self.progress)

    def execute(self, command, directory, stem, timeout, mpi=False, steps=None):
        directory = Path(directory)
        require(directory.is_dir(), "command working directory does not exist")
        receipt_path, log = directory / f"{stem}.receipt.json", directory / f"{stem}.log"
        require(not receipt_path.exists() and not log.exists(), "command output already exists; no automatic retry")
        verify_files(self.bindings)
        environment = dict(os.environ, **self.config["allocation"]["thread_environment"],
                           PYTHONDONTWRITEBYTECODE="1", TMPDIR=str(directory / "tmp"))
        (directory / "tmp").mkdir(exist_ok=True)
        environment["MPLCONFIGDIR"] = str(directory / "mplconfig")
        cpus = self.config["allocation"]["cpus"]
        command = ["taskset", "-c", ",".join(map(str, cpus if mpi else cpus[:1]))] + list(command)
        result = {"schema": SCHEMA, "status": "RUNNING", "command": command,
                  "working_directory": str(directory), "started_utc": utc(), "timeout_seconds": timeout,
                  "environment": {k: environment[k] for k in self.config["allocation"]["thread_environment"]},
                  "peak_sampled_tree_rss_kib": 0, "peak_sampled_tree_threads": 0,
                  "resource_sample_seconds": 2, "mpi": mpi}
        result["input_files_sha256"] = {str(p): sha(p) for p in sorted(directory.glob("*.xml"))
                                         if p.name == "pilot.xml"}
        result["input_files_sha256"].update({str(p): sha(p) for p in
            (directory / "manifest.json", directory / "cbc1971_theta42_spectrum_3d.dat") if p.is_file()})
        atomic_json(receipt_path, result)
        self.update(current={"directory": str(directory), "command": stem, "log": str(log), "total_steps": steps})
        start, process, failure = time.monotonic(), None, None
        binding_dir = directory / f"{stem}.bindings"
        affinity_audit = AffinityAudit(cpus, self.config["allocation"]["ranks"]) if mpi else None
        if mpi:
            binding_dir.mkdir()
        try:
            with log.open("x", encoding="utf-8") as stream:
                process = subprocess.Popen(command, cwd=directory, env=environment, stdin=subprocess.DEVNULL,
                                           stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
                result["pid"] = process.pid
                atomic_json(receipt_path, result)
                while process.poll() is None:
                    rss, threads = resources(process.pid)
                    result["peak_sampled_tree_rss_kib"] = max(result["peak_sampled_tree_rss_kib"], rss)
                    result["peak_sampled_tree_threads"] = max(result["peak_sampled_tree_threads"], threads)
                    if mpi:
                        affinity_audit.observe_log(tail(log))
                        validate_bindings(binding_dir, cpus, self.config["allocation"]["ranks"], live=True, audit=affinity_audit)
                        result["affinity_audit"] = affinity_audit.snapshot()
                    require(time.monotonic() - start <= timeout, f"hard timeout after {timeout}s")
                    current = dict(self.progress["current"])
                    current.update(log_progress(tail(log)) or {})
                    self.update(current=current)
                    atomic_json(receipt_path, result)
                    time.sleep(2)
            result["exit_code"] = process.returncode
            require(process.returncode == 0, f"{stem} exited {process.returncode}; inspect {log}")
            text = log.read_text(errors="replace")
            problem = BAD_LOG.search(text)
            require(problem is None, f"{stem} has a native failure/nonfinite marker: {problem.group(0) if problem else ''}")
            if mpi:
                affinity_audit.observe_log(text, completed=True)
                result["rank_bindings"] = validate_bindings(binding_dir, cpus, self.config["allocation"]["ranks"])
            if steps is not None:
                last = log_progress(text)
                require(last is not None and last["step"] == steps, "native solver did not reach the declared final step")
                result["final_progress"] = last
            verify_files(self.bindings)
            for filename, expected in result["input_files_sha256"].items():
                require(sha(filename) == expected, "native input changed while command was running")
            result["status"] = "COMPLETED"
        except BaseException as error:
            failure = error
            if process is not None and process.poll() is None:
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
            if mpi:
                try:
                    result["failed_command_rank_cleanup"] = cleanup_owned_ranks(binding_dir)
                except (OSError, ValueError) as cleanup_error:
                    result["rank_cleanup_error"] = str(cleanup_error)
            result.update(status="FAILED", error=f"{type(error).__name__}: {error}",
                          exit_code=process.returncode if process is not None else None)
        finally:
            result.update(finished_utc=utc(), wall_seconds=time.monotonic()-start)
            if mpi:
                result["affinity_audit"] = affinity_audit.snapshot()
            if log.exists():
                result.update(log_sha256=sha(log), log_bytes=log.stat().st_size)
            atomic_json(receipt_path, result)
            self.receipt["commands"].append({"receipt": str(receipt_path), "sha256": sha(receipt_path), "status": result["status"]})
            atomic_json(self.root / "receipt.json", self.receipt)
        if failure is not None:
            raise failure
        return result

    def native_command(self, binary_name, directory, stem, args):
        allocation_ = self.config["allocation"]
        return ["mpirun", "--host", f"localhost:{allocation_['ranks']}", "--bind-to", "none", "-n", str(allocation_["ranks"]),
                sys.executable, str(Path(__file__).resolve()), "rank", "--cpus", ",".join(map(str, allocation_["cpus"])),
                "--binding-dir", str(Path(directory) / f"{stem}.bindings"),
                "--binary", self.bindings["executables"][binary_name]["path"], "--"] + list(args)

    def prepare(self, directory, name, smoke):
        knobs = case_knobs(name)
        duration = 160 * knobs["dt"] if smoke else DURATION
        record = duration / 2 if smoke else RECORD_START
        every = 20 if smoke else round(SAMPLE_DT / knobs["dt"])
        # For the half-dt smoke, retain160steps but use an exact eight-frame clock.
        directory.mkdir(parents=True)
        command = [sys.executable, str(HERE / "prepare.py"), "--output", str(directory), "--length", "1", "--width", ".6",
                   "--h", str(knobs["h"]), "--modes", "11", "--hom-modes", str(knobs["nz"]), "--dt", str(knobs["dt"]),
                   "--duration", str(duration), "--record-start", str(record), "--sample-every", str(every),
                   "--cutoff", str(knobs["cutoff"]), "--taper", str(knobs["taper"]), "--synthesis-length", "10",
                   "--synthesis-counts", COUNTS, "--seed", SEED, "--buffer", ".25", "--stations", "0,.15,.45,.75",
                   "--probe-grid", "2,2" if smoke else "32,32", "--no-stresses", "--dealias", "--svv", "dg",
                   "--svv-coeff", str(knobs["svv"]), "--svv-z-coeff", str(knobs["svv"]*.1*U*WIDTH/knobs["nz"]),
                   "--svv-z-cutoff", ".75", "--pressure-solver", "XxtMultiLevelStaticCond"]
        self.execute(command, directory, "prepare", 120)
        session, manifest_path = directory / "pilot.xml", directory / "manifest.json"
        tree = ET.parse(session)
        parameters = tree.find("./CONDITIONS/PARAMETERS")
        for key, value in (("HITRequireAccelerated", 1), ("HITReportEvaluator", 1), ("HITProbeDirectSamples", 1)):
            ET.SubElement(parameters, "P").text = f"{key} = {value}"
        # Ensure a final status line even when a generic preparer's info interval
        # does not divide the requested step count.
        steps = round(duration / knobs["dt"])
        for item in parameters.findall("P"):
            if item.text.split("=", 1)[0].strip() == "IO_InfoSteps":
                item.text = f"IO_InfoSteps = {8 if smoke else 1000}"
        ET.indent(tree, space="    ")
        tree.write(session, encoding="utf-8", xml_declaration=True)
        manifest = json.loads(manifest_path.read_text())
        manifest["sha256"]["pilot.xml"] = sha(session)
        manifest["band100_controller"] = {"case": name, "smoke": smoke, "mean_only_initialization": True,
                                           "require_accelerated": True, "probe_direct_samples_per_rank_boundary_time": 1,
                                           "controller_sha256": sha(__file__)}
        atomic_json(manifest_path, manifest)
        return manifest, steps

    def solver(self, directory, steps, timeout):
        args = ["--npz", str(self.config["allocation"]["npz"]), "--use-scotch", "pilot.xml"]
        result = self.execute(self.native_command("solver", directory, "solver", args), directory,
                              "solver", timeout, mpi=True, steps=steps)
        reports = evaluator_reports((directory / "solver.log").read_text(errors="replace"),
                                    self.config["allocation"]["ranks"])
        result["evaluator_reports"] = reports
        result["evaluator_methods"] = sorted({r["method"] for r in reports})
        result["nonempty_inlet_evaluator_ranks"] = [r["rank"] for r in reports if r["points"] > 0]
        atomic_json(directory / "solver.receipt.json", result)
        self.receipt["commands"][-1]["sha256"] = sha(directory / "solver.receipt.json")
        atomic_json(self.root / "receipt.json", self.receipt)
        return result

    def finish(self, status, **values):
        self.receipt.update(values, status=status, finished_utc=utc())
        atomic_json(self.root / "receipt.json", self.receipt)
        self.update(status=status, current=None, **values)


def timing_gate(log_text, command_receipt, knobs, max_hours=12, max_memory_gib=180):
    rows = {}
    for match in STEP.finditer(log_text):
        step, seconds = int(match.group(1)), float(match.group(3))
        require(math.isfinite(seconds) and seconds > 0, "invalid native info-chunk timing")
        rows[step] = max(rows.get(step, 0), seconds)
    ordered = sorted(rows)
    require(len(ordered) >= 4, "timing smoke needs at least four native info chunks")
    rates, previous = [], 0
    for step in ordered:
        rate = rows[step] / (step-previous)
        if step > ordered[-1] / 2:
            rates.append(rate)
        previous = step
    require(rates and all(math.isfinite(v) and v > 0 for v in rates), "invalid steady-step estimate")
    rate = max(rates)
    one_off = max(0, command_receipt["wall_seconds"] - sum(rows.values()))
    steps = round(DURATION / knobs["dt"])
    forecast = one_off + rate*steps
    memory = command_receipt["peak_sampled_tree_rss_kib"] / 1024**2
    return {"passed": forecast <= max_hours*3600 and memory <= max_memory_gib,
            "steady_seconds_per_step": rate, "steady_rate_estimator": "maximum native info-chunk seconds/step in the second half of the160-step smoke",
            "steady_chunk_rates_seconds_per_step": rates, "one_off_wall_overhead_seconds": one_off,
            "forecast_duration_s": DURATION, "forecast_steps": steps, "forecast_wall_seconds": forecast,
            "sampled_peak_tree_rss_gib": memory, "maximum_forecast_hours": max_hours,
            "maximum_sampled_memory_gib": max_memory_gib,
            "scope": "Short mean-initial smoke; no guarantee of developed-flow runtime or peak memory"}


def finite_smoke_history(directory, manifest):
    import analyze_history as ah
    path = directory / manifest["history_file"]
    count = 0
    for _, _values in ah.history_frames(path, manifest["sampling"]["probes"], manifest["time"]["sample_dt_s"]):
        count += 1
    expected = manifest["time"]["steps"] // manifest["time"]["sample_every_steps"] + 1
    require(count == expected, "native smoke history is incomplete")
    return {"path": str(path), "sha256": sha(path), "complete_finite_frames": count}


def check_preflight(path, bindings, allocation_):
    path = Path(path).resolve()
    receipt_path = path / "receipt.json"
    receipt = json.loads(receipt_path.read_text())
    require(receipt.get("schema") == SCHEMA and receipt.get("status") == "COMPLETED"
            and receipt["configuration"]["stage"] == "preflight", "preflight is not completed")
    require(receipt["bindings"] == bindings, "preflight source/executable bindings differ from this study")
    require(receipt["configuration"]["allocation"] == allocation_, "preflight allocation differs from this study")
    require(receipt["budget_gate"]["passed"], "preflight did not pass the declared wall/memory budget")
    for record in receipt["commands"]:
        require(sha(record["receipt"]) == record["sha256"] and record["status"] == "COMPLETED", "preflight command receipt changed or failed")
        command = json.loads(Path(record["receipt"]).read_text())
        logfile = Path(record["receipt"]).with_name(Path(record["receipt"]).name.replace(".receipt.json", ".log"))
        require(command["status"] == "COMPLETED" and sha(logfile) == command["log_sha256"],
                "preflight native log changed or failed")
        for filename, expected in command["input_files_sha256"].items():
            require(sha(filename) == expected, "preflight native input changed")
    history = receipt["smoke_history"]
    require(sha(history["path"]) == history["sha256"], "preflight smoke history changed")
    return {"path": str(path), "receipt_sha256": sha(receipt_path), "budget_gate": receipt["budget_gate"]}


def preflight(controller):
    directory = controller.root / "smoke"
    manifest, steps = controller.prepare(directory, "base", smoke=True)
    benchmark = controller.bindings["executables"]["benchmark"]["path"]
    controller.execute([benchmark, "pilot.xml", "2", "2"], directory, "benchmark", 3600)
    native = json.loads((directory / "pilot.benchmark.json").read_text())
    require(native["method"].startswith("Grouped") and native["max_abs_error"] <= 1e-9,
            "native generator benchmark failed acceleration/direct agreement")
    arguments = ["--npz", str(controller.config["allocation"]["npz"]), "--use-scotch", "pilot.xml"]
    controller.execute(controller.native_command("probe", directory, "projection", arguments), directory,
                       "projection", 3600, mpi=True)
    projection_evaluators = evaluator_reports((directory / "projection.log").read_text(errors="replace"),
                                              controller.config["allocation"]["ranks"])
    projection = json.loads((directory / "pilot.probe.json").read_text())
    require(projection["direct_reference_samples"] > 0 and projection["raw_max_error"] <= 1e-9,
            "native boundary direct-reference samples failed")
    require(all(math.isfinite(v) for v in projection["projection_max_error"]), "nonfinite boundary projection error")
    solve = controller.solver(directory, steps, 3600)
    history = finite_smoke_history(directory, manifest)
    gate = timing_gate((directory / "solver.log").read_text(), solve, case_knobs("base"))
    controller.finish("COMPLETED", budget_gate=gate, native_benchmark=native,
                      boundary_projection=projection, boundary_evaluator_reports=projection_evaluators,
                      smoke_history=history,
                      study_eligible=gate["passed"],
                      interpretation="Preflight completed; no downstream physical or precision-qualified spectral result")


def previous_analysis(controller, name, reference_study):
    local = controller.root / "analysis" / name / "analysis.json"
    if local.exists():
        return local
    if reference_study is not None:
        previous = Path(reference_study).resolve()
        receipt = json.loads((previous / "receipt.json").read_text())
        require(receipt["schema"] == SCHEMA and receipt["status"] == "COMPLETED" and
                receipt["bindings"] == controller.bindings, "reference study is incomplete or has different source/executables")
        result = previous / "analysis" / name / "analysis.json"
        require(result.is_file() and name in receipt["cases"], "reference study lacks required matched analysis")
        require(sha(result) == receipt["cases"][name]["analysis_sha256"], "reference analysis hash mismatch")
        return result
    raise ValueError(f"matched reference {name} is missing; include it first or pass --reference-study")


def study(controller, cases, reference_study):
    for name in cases:
        reference = None
        if name != "base":
            reference = previous_analysis(controller, "z" if name == "bandwidth" else "base", reference_study)
        smoke = controller.root / "cases" / name / "smoke"
        manifest, steps = controller.prepare(smoke, name, smoke=True)
        smoke_result = controller.solver(smoke, steps, 3600)
        finite_smoke_history(smoke, manifest)
        gate = timing_gate((smoke / "solver.log").read_text(), smoke_result, case_knobs(name))
        if not gate["passed"]:
            controller.finish("BUDGET_REJECTED", rejected_case=name, budget_gate=gate,
                              interpretation="A fresh smoke rejected this recipe; no full recording was launched or retried")
            return
        directory = controller.root / "cases" / name / "run"
        manifest, steps = controller.prepare(directory, name, smoke=False)
        controller.solver(directory, steps, 16*3600)
        analysis_dir = controller.root / "analysis" / name
        command = [sys.executable, str(HERE / "analyze_band100.py"), "--manifest", str(directory / "manifest.json"),
                   "--output", str(analysis_dir)]
        if reference is not None:
            command += ["--reference-analysis", str(reference)]
        controller.execute(command, controller.root, f"analyze-{name}", 3600)
        analysis_path = analysis_dir / "analysis.json"
        require(analysis_path.is_file(), "analysis did not produce its declared result")
        analysis_receipt = json.loads((analysis_dir / "receipt.json").read_text())
        require(analysis_receipt.get("status") == "completed" and
                analysis_receipt.get("analysis_sha256") == sha(analysis_path),
                "analysis has no matching completed receipt")
        controller.receipt["cases"][name] = {"budget_gate": gate, "directory": str(directory),
                                             "manifest_sha256": sha(directory / "manifest.json"),
                                             "analysis": str(analysis_path), "analysis_sha256": sha(analysis_path)}
        atomic_json(controller.root / "receipt.json", controller.receipt)
        controller.progress["completed_cases"].append(name)
        controller.update(current=None)
    controller.finish("COMPLETED", interpretation="Selected short-pilot cases completed; full-distance production remains unlaunched and unqualified")


def status(root):
    root = Path(root).expanduser().resolve()
    progress = json.loads((root / "progress.json").read_text())
    print(f"{utc()} — CBC100 {progress['stage']}: {progress['status']}")
    print("Completed cases: " + (", ".join(progress["completed_cases"]) or "none"))
    current = progress.get("current")
    if current:
        current = dict(current)
        current.update(log_progress(tail(Path(current["log"]))) or {})
        text = f"Current: {current['command']} in {current['directory']}"
        if "step" in current:
            text += f"; step {current['step']}/{current.get('total_steps', '?')}, t={current['simulated_time_s']:.8g}s"
        print(text)
    if "budget_gate" in progress:
        gate = progress["budget_gate"]
        print(f"Budget gate: {'pass' if gate['passed'] else 'reject'}; forecast {gate['forecast_wall_seconds']/3600:.3f}h, sampled RSS {gate['sampled_peak_tree_rss_gib']:.3f}GiB")
    if "error" in progress:
        print("Error: " + progress["error"])


def main(arguments=None):
    arguments = list(sys.argv[1:] if arguments is None else arguments)
    if arguments[:1] == ["rank"]:
        rank_entry(arguments[1:])
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--status", type=Path, help="Read status only; do not launch any command")
    parser.add_argument("--stage", choices=("preflight", "study"))
    parser.add_argument("--run-root", type=Path)
    parser.add_argument("--source-revision")
    parser.add_argument("--ranks", type=int, default=48, help="MPI ranks, maximum48 (default48)")
    parser.add_argument("--npz", type=int, default=16)
    parser.add_argument("--cpus")
    parser.add_argument("--preflight", type=Path)
    parser.add_argument("--reference-study", type=Path)
    parser.add_argument("--cases", default="base,dt", help="Sequential subset of base,dt,z,bandwidth,svv,h")
    parser.add_argument("--solver", default="IncNavierStokesSolver")
    parser.add_argument("--benchmark", default="HITInflowBenchmark")
    parser.add_argument("--boundary-probe", default="HITBoundaryProbe")
    args = parser.parse_args(arguments)
    if args.status:
        require(not any((args.stage, args.run_root, args.source_revision, args.cpus, args.preflight)),
                "--status cannot be combined with execution options")
        status(args.status)
        return
    require(sys.platform.startswith("linux") and hasattr(os, "sched_setaffinity"), "native pilot execution is Linux-only")
    require(all(value is not None for value in (args.stage, args.run_root, args.source_revision, args.ranks, args.cpus)),
            "execution requires --stage, --run-root, --source-revision, --ranks and --cpus")
    cpus = parse_cpus(args.cpus)
    allocated = allocation(cpus, args.ranks, args.npz, os.sched_getaffinity(0))
    root = args.run_root.expanduser().resolve()
    require(not root.exists() and not within(root, SOURCE), "run-root must be fresh and outside the source checkout")
    require(within(root, Path("/media/jay/新加卷")), "heavy pilot outputs must be on /media/jay/新加卷")
    bindings = {"source": source_binding(args.source_revision), "executables": {
        "solver": executable_binding(args.solver), "benchmark": executable_binding(args.benchmark),
        "probe": executable_binding(args.boundary_probe)}}
    configuration = {"stage": args.stage, "allocation": allocated, "controller_command": [sys.executable, str(Path(__file__).resolve())] + arguments,
                     "host": os.uname().nodename, "python": sys.version, "mean_only_initialization": True}
    if args.stage == "study":
        require(args.preflight is not None, "study requires --preflight from a completed matching preflight")
        configuration["preflight"] = check_preflight(args.preflight, bindings, allocated)
        cases = args.cases.split(",")
        require(cases and len(cases) == len(set(cases)) and all(c in CASE_ORDER for c in cases), "--cases needs unique known case names")
        configuration["cases"] = cases
    else:
        require(args.preflight is None and args.reference_study is None, "preflight does not accept a prior study/preflight")
    root.mkdir(parents=True)
    # The controller itself is bounded and shares one of the rank CPUs. Child
    # taskset widens only to the explicitly allocated pool, before singleton rank pinning.
    os.sched_setaffinity(0, {cpus[-1]})
    controller = Controller(root, bindings, configuration)
    try:
        if args.stage == "preflight":
            preflight(controller)
        else:
            study(controller, cases, args.reference_study)
    except BaseException as error:
        controller.finish("FAILED", error=f"{type(error).__name__}: {error}")
        raise
    finally:
        status(root)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, RuntimeError, OSError, subprocess.SubprocessError) as error:
        print(f"run_band100: {error}", file=sys.stderr)
        raise SystemExit(1)
