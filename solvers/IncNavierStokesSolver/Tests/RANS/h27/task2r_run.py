#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Carry forward the H27 ceiling and run serial Task 2R operations on Linux."""

import argparse
import fcntl
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import time

from resources import THREADS, atomic, identity, size, tree, update


UUID = "cb7dc460-f22f-4533-8329-245d10620f47"
BASELINE = "50c57af686e39ab51fd9119eaf126e92264e49dd20f82b90963e2503a9edb016"
GIB = 1024 ** 3
LIMITS = {"analysis": 900, "build": 1800, "science": 7200, "evidence": 900}
WORK = ("simpleFoam", "IncNavierStokes", "RANS", "WJ", "python", "cc1",
        "c++", "g++", "gcc", "clang", "gfortran", "f951", "lto1",
        "orted", "mpirun", "mpiexec", "ninja", "make", "cmake")


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def workers():
    rows = subprocess.check_output(
        ["ps", "-eLo", "pid=,tid=,state=,wchan=,comm="], text=True)
    counted = []
    helpers = 0
    for row in rows.splitlines():
        pid, tid, state, wait, comm = row.split(None, 4)
        if state in ("X", "Z") or comm.startswith("kworker"):
            continue
        if tid != pid and comm.startswith(("simpleFoam", "IncNavierStokes", "WJ")) \
                and state == "S" and wait in ("do_pol", "ep_pol"):
            # Idle OpenMPI event listeners. A runnable helper is counted.
            helpers += 1
            continue
        if state in ("R", "D") or comm.startswith(WORK):
            counted.append(comm)
    return len(counted), helpers


def storage(root):
    return {
        "science_bytes": size(root / "science") + size(root / "tmp"),
        "source_build_runtime_bytes": sum(size(p) for p in root.iterdir()
                                          if p.is_dir() and p.name.startswith(("source", "build", "runtime"))),
        "free_bytes": shutil.disk_usage(root).free,
    }


def check_storage(info):
    require(info["free_bytes"] >= 15 * GIB, "free-space limit")
    require(info["science_bytes"] <= 8 * GIB, "science storage limit")
    require(info["source_build_runtime_bytes"] <= 8 * GIB, "build storage limit")


def balance(book):
    spent = sum(e["seconds"] for e in book["entries"].values() if e["status"] != "RESERVED")
    reserved = sum(e["reserved_seconds"] for e in book["entries"].values() if e["status"] == "RESERVED")
    return dict(
        charged_seconds=book["opening_seconds"] + spent,
        reserved_seconds=reserved,
        nonprotected_remaining=book["ceiling"] - book["protected_seconds"]
        - book["opening_seconds"] - spent - reserved,
        tranche_remaining=book["tranche_seconds"] - spent - reserved,
    )


def initialize(root, book, began):
    baseline = root / "inputs" / "settled-pre-migration.json"
    old = json.loads(baseline.read_text())
    require(identity(baseline)["sha256"] == BASELINE and not old["account_closed"]
            and not old["reservations"], "invalid historical settlement")
    require(not book.exists(), "existing account")
    for directory in ("operations", "science", "tmp", "evidence"):
        (root / directory).mkdir(exist_ok=True)
    check_storage(storage(root))
    host, helpers = workers()
    require(host <= 72, "shared worker limit at migration")
    event = dict(status="PASS", phase="evidence",
                 seconds=time.monotonic() - began + 2.,
                 earlier_linux_inventory_charge_seconds=2.,
                 peak_host_workers=host, idle_mpi_event_helpers=helpers)
    atomic(book, dict(schema="WJH27Task2RContinuationV1", closed=False, ceiling=43200,
                      protected_seconds=old["evidence_protected"],
                      opening_seconds=old["aggregate_seconds"],
                      opening_nonprotected_remaining=old["global_nonprotected_remaining"],
                      initial_tranche_seconds=old["additional_fixture_remaining"],
                      tranche_seconds=old["additional_fixture_remaining"],
                      initial_evidence_seconds=300,
                      prior_account=old["account"], prior_settlement=identity(baseline),
                      entries={"r0-ledger-continuity": event}))
    return event


def terminate(child, groups):
    for group in groups | {child.pid}:
        try:
            os.killpg(group, signal.SIGTERM)
        except ProcessLookupError:
            pass
    try:
        return child.wait(timeout=4)
    except subprocess.TimeoutExpired:
        for group in groups | {child.pid}:
            try:
                os.killpg(group, signal.SIGKILL)
            except ProcessLookupError:
                pass
        return child.wait()


def run(root, book, args, began):
    phase, label = args.phase, args.label
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    require(command and label and label.replace("-", "").isalnum(), "invalid command or label")
    require(0 < args.seconds <= LIMITS[phase], "invalid duration cap")
    ledger = json.loads(book.read_text())
    require(not ledger["closed"] and label not in ledger["entries"], "closed account or repeated label")
    remaining = balance(ledger)
    evidence_spent = sum(e.get("seconds", 0) for e in ledger["entries"].values()
                         if e["phase"] == "evidence")
    hold = 0 if phase == "evidence" else max(0, ledger["initial_evidence_seconds"] - evidence_spent)
    cap = args.seconds + 7  # Allow bounded startup, monitoring and signal grace.
    require(remaining["reserved_seconds"] == 0, "unsettled operation exists")
    require(cap <= remaining["nonprotected_remaining"] and cap + hold <= remaining["tranche_remaining"],
            "operation cap exceeds live tranche or parent allowance")
    operation = root / "operations" / label
    host, helpers = workers()
    limits = {"build": 16, "science": 10, "analysis": 6, "evidence": 5}
    info = storage(root)
    request = dict(label=label, phase=phase, command=command, cap_seconds=args.seconds,
                   reserved_seconds=cap, workers_before=host, idle_mpi_event_helpers=helpers,
                   threads=THREADS, storage_before=info, runner=identity(__file__),
                   input_files=[identity(p) for p in command
                                if p.startswith("/") and "\n" not in p and Path(p).is_file()])
    operation.mkdir()
    atomic(operation / "request.json", request)
    update(book, lambda ledger: ledger["entries"].update(
        {label: {"status": "RESERVED", "phase": phase, "reserved_seconds": cap}}))
    peak, peak_workers, code, reason, groups = 0, host, None, None, set()
    env = dict(os.environ, **THREADS, PYTHONDONTWRITEBYTECODE="1",
               TMPDIR=str(root / "tmp"), H27_ROOT=str(root), H27_OPERATION=str(operation))
    log_path = operation / "command.log"
    with log_path.open("x") as log:
        try:
            check_storage(info)
            require(host + limits[phase] <= 72, "shared worker limit at admission")
            with subprocess.Popen(command, cwd=root, env=env, stdout=log,
                                  stderr=subprocess.STDOUT, start_new_session=True) as child:
                try:
                    while child.poll() is None:
                        rss, child_groups = tree(child.pid)
                        groups |= child_groups
                        peak = max(peak, rss)
                        peak_workers = max(peak_workers, workers()[0])
                        require(peak <= (16 if phase == "build" else 8) * GIB, "RSS_LIMIT")
                        require(peak_workers <= 72, "HOST_WORKER_LIMIT")
                        require(time.monotonic() - began < args.seconds, "TIME_LIMIT")
                        check_storage(storage(root))
                        require(log_path.stat().st_size <= 32 * 1024 ** 2, "LOG_LIMIT")
                        time.sleep(.5)
                    code = child.returncode
                except (RuntimeError, OSError, KeyboardInterrupt) as error:
                    code, reason = terminate(child, groups), str(error)
        except (RuntimeError, OSError, KeyboardInterrupt) as error:
            reason = str(error)
    event = dict(status="PASS" if code == 0 and reason is None else "FAIL", phase=phase,
                 seconds=time.monotonic() - began, exit=code, reason=reason,
                 peak_tree_rss_bytes=peak, peak_host_workers=peak_workers, log=identity(log_path))
    ledger = update(book, lambda ledger: ledger["entries"].update({label: event}))
    atomic(operation / "receipt.json", dict(request=request, result=event,
                                            balance=balance(ledger), storage=storage(root)))
    print(label, event["status"], event["seconds"], reason, flush=True)
    return int(event["status"] != "PASS")


def main():
    began = time.monotonic()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    modes = parser.add_subparsers(dest="action", required=True)
    modes.add_parser("initialize")
    modes.add_parser("balance")
    extend = modes.add_parser("allocate")
    extend.add_argument("seconds", type=float)
    extend.add_argument("forecast", type=Path)
    op = modes.add_parser("run")
    op.add_argument("label")
    op.add_argument("phase", choices=LIMITS)
    op.add_argument("seconds", type=float)
    op.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    mounted = subprocess.check_output(["findmnt", "-rn", "-o", "FSTYPE,UUID", "-T", str(root)],
                                      text=True).split()
    require(root.name.startswith("wj-h27-task2r-") and mounted == ["ext4", UUID], "wrong Linux store")
    book = root / "accounting.json"
    with (root / "execution.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if args.action == "run":
            return run(root, book, args, began)
        if args.action == "initialize":
            print(json.dumps(initialize(root, book, began), indent=2))
        elif args.action == "allocate":
            def allocate(ledger):
                live = balance(ledger)
                extra = args.seconds
                require(not ledger["closed"] and live["reserved_seconds"] == 0
                        and 0 < extra <= live["nonprotected_remaining"] - live["tranche_remaining"],
                        "allocation exceeds available parent balance")
                ledger.setdefault("allocations", []).append(dict(seconds=extra,
                    live_before=live, measured_forecast=identity(args.forecast)))
                ledger["tranche_seconds"] += extra

            update(book, allocate)
        print(json.dumps(balance(json.loads(book.read_text())), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
