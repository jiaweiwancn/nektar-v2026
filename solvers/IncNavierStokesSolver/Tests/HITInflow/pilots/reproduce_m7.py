#!/usr/bin/env python3
"""Sequential Linux reproduction of M7 analysis, native-work joins and figures.

Consumes the original retained archive, normally
/media/jay/新加卷/nektar-v2026-runs/hit-inflow/milestone-7/97014ab.
Its receipts contain absolute paths; relocating the archive is not supported.
No solver, field import or new native profiling is run. NumPy/Matplotlib are
dependencies of the unchanged helper subprocesses, not of this driver.
"""

import argparse
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time


SEEDS = (20261028, 4315228324)
CASES = tuple((f"{prefix}-seed-{seed}", family, seed, f"{prefix}-production")
              for prefix, family in (("cbc", "cbc1971"), ("vk", "von-karman")) for seed in SEEDS)
POST = "postprocessing-da607a9"
RELATIVE_TOOLS = Path("solvers/IncNavierStokesSolver/Tests/HITInflow/pilots")
HELPERS = ("analyze_m7.py", "analyze_history.py", "collect_m7_work.py", "plot_m7.py")
FIGURES = ("m7-decay", "m7-longitudinal-spectra", "m7-delayed-correlations",
           "m7-sampling-qualification", "m7-native-work")
THREAD_ENV = {key: "1" for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS", "NUMEXPR_MAX_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS")}
THREAD_ENV["PYTHONDONTWRITEBYTECODE"] = "1"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def inside(path, root):
    return path == root or root in path.parents


def source_root(explicit):
    candidates = [explicit.expanduser().resolve()] if explicit is not None else Path(__file__).resolve().parents
    for path in candidates:
        if (path/"library/SolverUtils").is_dir() and (path/"solvers/IncNavierStokesSolver").is_dir():
            return path
    raise ValueError("cannot locate source checkout; supply --source-root for an external driver copy")


def output_guard(output, run_root, checkout):
    original = output.expanduser()
    target = original.resolve()
    require(not original.exists() and not original.is_symlink(), "--output must be a fresh directory")
    for root in (run_root, checkout):
        require(not inside(target, root) and not inside(root, target),
                "output must be outside source/case/archive input trees and cannot contain them")
    require(target.parent.is_dir(), "output parent must already exist")
    return target


def finite_tree(value):
    if isinstance(value, float):
        require(math.isfinite(value), "nonfinite JSON receipt number")
    elif isinstance(value, dict):
        for item in value.values():
            finite_tree(item)
    elif isinstance(value, list):
        for item in value:
            finite_tree(item)


def fingerprint(path):
    before = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4*1024*1024), b""):
            digest.update(block)
    after = path.stat()
    require((before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns), f"input changed: {path}")
    return {"sha256": digest.hexdigest(), "bytes": after.st_size}


def identity(path):
    stat = path.stat()
    require(path.is_file(), f"expected a regular input file: {path}")
    return {"bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns, "device": stat.st_dev, "inode": stat.st_ino}


class Inputs:
    def __init__(self):
        self.small, self.large = {}, {}

    def add(self, path, expected=None):
        path = path.resolve()
        record = fingerprint(path)
        require(expected is None or record == expected, f"input hash/size mismatch: {path}")
        require(str(path) not in self.small or self.small[str(path)] == record, f"input changed: {path}")
        self.small[str(path)] = record
        return path

    def read(self, path, expected=None):
        path = self.add(path, expected)
        value = json.loads(path.read_text(encoding="utf-8"))
        finite_tree(value)
        require(fingerprint(path) == self.small[str(path)], f"receipt changed while parsing: {path}")
        return value

    def track_large(self, path, declared_sha256=None):
        path = path.resolve()
        record = identity(path)
        self.large[str(path)] = {"identity": record, "declared_sha256": declared_sha256}
        return record

    def verify(self):
        for path, expected in self.small.items():
            require(fingerprint(Path(path)) == expected, f"immutable input changed: {path}")
        for path, expected in self.large.items():
            require(identity(Path(path)) == expected["identity"], f"immutable large-input identity changed: {path}")


def pin_cpu(cpu):
    require(sys.platform.startswith("linux") and hasattr(os, "sched_getaffinity"), "run this driver on authoritative Linux")
    require(cpu >= 0 and cpu in os.sched_getaffinity(0), "--cpu must be one CPU in the inherited affinity set")
    os.sched_setaffinity(0, {cpu})
    selected = sorted(os.sched_getaffinity(0))
    require(selected == [cpu], "failed to enforce the requested single CPU")
    return selected


def archive_path(path, archive):
    path = path.resolve()
    require(inside(path, archive), f"receipt points outside original retained archive: {path}")
    return path


def preflight(args):
    checkout, archive = source_root(args.source_root), args.run_root.expanduser().resolve()
    require(archive.is_dir(), "--run-root must be the existing original retained archive")
    require(not inside(archive, checkout), "scientific archive must be outside the source checkout")
    output = output_guard(args.output, archive, checkout)
    inputs = Inputs()
    tools = checkout/RELATIVE_TOOLS
    for name in HELPERS:
        path = (tools/name).resolve()
        require(inside(path, checkout), "helper path aliases outside the protected checkout")
        inputs.add(path)
    inputs.add(Path(__file__))
    experimental = tools.parent/"cbc1971-source.json"
    correlation = tools.parent/"cbc1971_full_band_correlation.json"
    require(inputs.read(experimental)["schema"] == "nektar.hit.cbc1971.source.v1", "unsupported experimental source")
    require(inputs.read(correlation)["schema"] == "nektar.hit.cbc1971.full_band_correlation.v1", "unsupported correlation source")
    for production in ("cbc-production", "vk-production"):
        study = inputs.read(archive_path(archive/production/"study.receipt.json", archive))
        expected_ids = {name for name, _, _, group in CASES if group == production}
        rows = study["cases"]
        require(study["status"] == "completed" and len(rows) == 2
                and {row["id"] for row in rows} == expected_ids and all(row["status"] == "completed" for row in rows),
                f"{production}: require a successful completed two-case study receipt")
    cases = []
    for name, family, seed, production in CASES:
        case = archive_path(archive/production/name, archive)
        manifest_path = archive_path(case/"manifest.json", archive)
        require(manifest_path.parent == case, "manifest aliases another case")
        manifest = inputs.read(manifest_path)
        require(manifest["schema"] == "NEKTAR_HIT_PILOT_V1" and str(manifest["synthesis"]["seed_uint64"]) == str(seed),
                f"{name}: invalid manifest/seed")
        require(manifest["time"]["record_start_s"] == .8 and manifest["time"]["end_s"] == 1.6,
                f"{name}: archive is not the frozen M7 recording window")
        session = Path(manifest["case"])
        require(session.name == str(session), "manifest session must be a local filename")
        history_name = Path(manifest["history_file"])
        require(history_name.name == str(history_name), "manifest history must be a local filename")
        history = archive_path(case/history_name, archive)
        session = inputs.add(archive_path(case/session, archive))
        run = inputs.read(archive_path(case/"run.receipt.json", archive))
        require(run.get("finished_utc") and all(run[key]["exit_code"] == 0
                and not run[key].get("native_iteration_limit_exceeded", False) for key in ("initialization", "solver")),
                f"{name}: initialization/solver receipt is not successful")
        require(run["executed_session_sha256"] == inputs.small[str(session)]["sha256"], f"{name}: executed session changed")
        retained = archive_path(archive/POST/name, archive)
        native_path = archive_path(retained/"native-profiles/receipt.json", archive)
        native = inputs.read(native_path)
        require(native["schema"] == "NEKTAR_HIT_CHECKPOINT_PROFILES_V1" and native["status"] == "completed"
                and native.get("finished_utc") and Path(native["case"]).resolve() == case,
                f"{name}: native profiles are incomplete, or archive was relocated")
        require(native["source_inputs"][str(manifest_path)]["sha256"] == inputs.small[str(manifest_path)]["sha256"],
                f"{name}: native manifest binding differs")
        for path, expected in native["source_inputs"].items():
            inputs.add(archive_path(Path(path), archive), expected)
        times = sorted(profile["time_s"] for profile in native["profiles"])
        require(len(times) == 3 and all(math.isclose(a, b, rel_tol=2e-12, abs_tol=1e-14)
                for a, b in zip(times, (1.3, 1.45, 1.6))), f"{name}: require the three late native snapshots")
        for profile in native["profiles"]:
            require(profile["slabs"] == 90, f"{name}: expected 90 native slabs")
            require(sum(path.endswith(".flow.json") for path in profile["output_files"]) == 1
                    and sum(path.endswith(".slabs.csv") for path in profile["output_files"]) == 1,
                    "native receipt needs flow JSON and slab CSV bindings")
            for path, expected in profile["output_files"].items():
                inputs.add(archive_path(Path(path), archive), expected)
            metadata = [path for path in profile["checkpoint_inputs"] if Path(path).name == "Info.xml"]
            require(len(metadata) == 1, "native receipt needs one checkpoint metadata binding")
            inputs.add(archive_path(Path(metadata[0]), archive), profile["checkpoint_inputs"][metadata[0]])
        if args.rebuild_cache:
            inputs.track_large(history)
            cache = output/name/"history.npz"
        else:
            cache = archive_path(retained/"history.npz", archive)
            sidecar = inputs.read(archive_path(cache.with_suffix(".npz.json"), archive))
            require(sidecar["schema"] == "NEKTAR_HIT_M7_HISTORY_CACHE_V1"
                    and sidecar["metadata"]["manifest_sha256"] == inputs.small[str(manifest_path.resolve())]["sha256"]
                    and sidecar["native_history_sha256_before_cache_write"] == sidecar["native_history_sha256_after_cache_write"]
                    == sidecar["metadata"]["native_history_sha256"], f"{name}: cache receipt binding mismatch")
            require(inputs.track_large(cache, sidecar["cache_sha256"])["bytes"] == sidecar["cache_bytes"], f"{name}: cache size changed")
        cases.append({"id": name, "family": family, "seed": seed, "directory": case,
                      "manifest": manifest_path, "history": history, "cache": cache, "native": native_path})
    return {"checkout": checkout, "archive": archive, "output": output, "tools": tools,
            "experimental": experimental, "correlation": correlation, "cases": cases, "inputs": inputs}


def commands(plan, rebuild):
    steps, analyses, natives = [], [], []
    for case in plan["cases"]:
        destination = plan["output"]/case["id"]
        analysis, work = destination/"analysis.json", destination/"work.json"
        command = [sys.executable, str(plan["tools"]/"analyze_m7.py"), "--manifest", str(case["manifest"]),
            "--cache", str(case["cache"]), "--output", str(analysis), "--source-root", str(plan["checkout"]),
            "--case-family", case["family"], "--case-id", case["id"], "--start", "0.8", "--end", "1.6",
            "--probe-batch", "256", "--welch-segment-samples", "512", "--acf-max-lag-s", "0.1",
            "--reference-x", "0", "--forced-inlet-x", "0", "--experimental-source", str(plan["experimental"]),
            "--correlation-source", str(plan["correlation"])]
        if rebuild:
            command += ["--history", str(case["history"]), "--cache-start", "0"]
        steps.append((case["id"]+"-analysis", command, [analysis]))
        steps.append((case["id"]+"-work", [sys.executable, str(plan["tools"]/"collect_m7_work.py"),
            "--analysis", str(analysis), "--profiles-receipts", str(case["native"]),
            "--source-root", str(plan["checkout"]), "--output", str(work)], [work]))
        analyses.append(str(analysis))
        natives += ["--native", case["id"]+"="+str(case["native"])]
    figures = plan["output"]/"figures"
    command = [sys.executable, str(plan["tools"]/"plot_m7.py"), "--analysis", *analyses, *natives,
               "--experimental-source", str(plan["experimental"]), "--correlation-source", str(plan["correlation"]),
               "--source-root", str(plan["checkout"]), "--output", str(figures)]
    products = [figures/f"{name}.{suffix}" for name in FIGURES for suffix in ("png", "svg")]
    steps.append(("figures", command, products+[figures/"provenance.json"]))
    return steps


def atomic_receipt(path, receipt):
    temporary = path.with_name("."+path.name+".tmp")
    try:
        with temporary.open("x", encoding="utf-8") as destination:
            json.dump(receipt, destination, indent=2, allow_nan=False)
            destination.write("\n")
            destination.flush()
            os.fsync(destination.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def reproduce(args, runner=None):
    affinity = pin_cpu(args.cpu)
    plan = preflight(args)
    inputs, output = plan["inputs"], plan["output"]
    steps = commands(plan, args.rebuild_cache)
    for _, _, products in steps:
        require(all(inside(p.resolve(), output) and str(p.resolve()) not in inputs.small
                    and str(p.resolve()) not in inputs.large for p in products), "generated product aliases an input")
    environment = dict(os.environ, **THREAD_ENV, MPLCONFIGDIR=str(output/"matplotlib-config"))
    output.mkdir(exist_ok=False)
    receipt_path = output/"receipt.json"
    receipt = {"schema": "NEKTAR_HIT_M7_REPRODUCTION_V1", "status": "running", "started_utc": now(),
        "command": sys.argv, "python": sys.version, "source_root": str(plan["checkout"]),
        "run_root": str(plan["archive"]), "output": str(output), "cpu_affinity": affinity,
        "thread_environment": THREAD_ENV, "max_concurrent_helpers": 1, "server_worker_limit": 64,
        "cache_mode": "rebuild fresh caches from native histories" if args.rebuild_cache else "reuse immutable retained NPZ caches",
        "verified_small_inputs": inputs.small, "large_inputs": inputs.large,
        "large_input_hash_scope": "Driver checks file identity/size across the run. Unchanged analyzer verifies reused NPZ hashes; rebuilding validates and hashes raw history before/after cache creation. No field partitions are reread.",
        "cases": [{key: str(value) if isinstance(value, Path) else value for key, value in case.items()} for case in plan["cases"]],
        "steps": [], "claim_limits": "Reproduction of retained M7 evidence only; no new CFD, profiling, fitting, ensemble independence or full-band validation claim."}
    atomic_receipt(receipt_path, receipt)
    runner = subprocess.run if runner is None else runner
    try:
        (output/"logs").mkdir()
        for case in plan["cases"]:
            (output/case["id"]).mkdir()
        for index, (name, command, products) in enumerate(steps, 1):
            inputs.verify()
            log = output/"logs"/f"{index:02d}-{name}.log"
            entry = {"name": name, "status": "running", "command": command, "cwd": str(output),
                     "command_sha256": hashlib.sha256(json.dumps(command, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest(),
                     "log": str(log), "started_utc": now()}
            receipt["steps"].append(entry)
            atomic_receipt(receipt_path, receipt)
            start = time.monotonic()
            try:
                with log.open("x") as stream:
                    process = runner(command, cwd=output, env=environment, stdin=subprocess.DEVNULL,
                                     stdout=stream, stderr=subprocess.STDOUT, check=False)
                entry.update(exit_code=process.returncode, wall_seconds=time.monotonic()-start,
                             finished_utc=now(), log_fingerprint=fingerprint(log))
                require(process.returncode == 0, f"{name} exited {process.returncode}; inspect {log}")
                entry["outputs"] = {str(p): fingerprint(p) for p in products}
                product = json.loads(products[-1].read_text(encoding="utf-8"))
                finite_tree(product)
                expected_schema = ("NEKTAR_HIT_M7_ANALYSIS_V1" if name.endswith("-analysis") else
                    "NEKTAR_HIT_M7_NATIVE_WORK_V1" if name.endswith("-work") else "NEKTAR_HIT_M7_FIGURES_V1")
                require(product["schema"] == expected_schema, f"{name}: unexpected output JSON schema")
                inputs.verify()
                entry["status"] = "completed"
            except BaseException as error:
                entry.update(status="failed", error=f"{type(error).__name__}: {error}", finished_utc=now(),
                             wall_seconds=time.monotonic()-start)
                if log.is_file():
                    entry["log_fingerprint"] = fingerprint(log)
                raise
            atomic_receipt(receipt_path, receipt)
        receipt.update(status="completed", finished_utc=now())
        atomic_receipt(receipt_path, receipt)
    except BaseException as error:
        receipt.update(status="failed", finished_utc=now(), error=f"{type(error).__name__}: {error}")
        atomic_receipt(receipt_path, receipt)
        raise
    return receipt_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="Fresh directory outside the archive and source trees")
    parser.add_argument("--cpu", type=int, required=True, help="Explicit single Linux CPU in inherited affinity")
    parser.add_argument("--source-root", type=Path, help="Checkout containing the unchanged M7 helpers; auto-discovered in a checkout")
    parser.add_argument("--rebuild-cache", action="store_true", help="Re-read native histories and write fresh caches under output")
    args = parser.parse_args()
    try:
        print(reproduce(args))
    except (OSError, ValueError, KeyError, TypeError, RuntimeError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
