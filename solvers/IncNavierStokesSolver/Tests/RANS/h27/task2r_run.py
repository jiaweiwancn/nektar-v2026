#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Carry forward the H27 ceiling and run serial Task 2R operations on Linux."""

import argparse
import fcntl
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import time
import xml.etree.ElementTree as ET

from resources import THREADS, atomic, identity, size, tree, update


UUID = "cb7dc460-f22f-4533-8329-245d10620f47"
BASELINE = "50c57af686e39ab51fd9119eaf126e92264e49dd20f82b90963e2503a9edb016"
GIB = 1024 ** 3
LIMITS = {"analysis": 900, "build": 1800, "science": 7200, "evidence": 900,
          "outlet-science": 7200}
OUTLET_AUTHORIZATION = "9356a31fab90117bf05438281dcf6deef4805b7c18f3323df088e684c94c0037"
S2_SETTLEMENT = "f2ae3f2c814d22dfa4161e976cd1792e2968b7303cbbe597186421383b52b64e"
S2_ACCOUNT = "dec44571419200acdb94c9e63c94c6a9a490f78caff158e0b31f92e5b933a784"
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
        rank_event = comm.startswith(("simpleFoam", "IncNavierStokes", "WJ")) and wait in ("do_pol", "ep_pol")
        launcher_event = comm in ("orted", "mpirun", "mpiexec") and wait in ("do_pol", "ep_pol", "do_sel")
        if tid != pid and state == "S" and (rank_event or launcher_event):
            # Parked OpenMPI event listeners, including launcher transport
            # threads. Main processes and runnable helpers remain counted.
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


def evidence_spend(book):
    return sum(e.get("seconds", 0) for e in book["entries"].values()
               if e["phase"] == "evidence")


def amend_outlet(root, book, authorization, began):
    name = "outlet-cp0-amendment"
    receipt = root / "evidence" / (name + ".json")
    mandate = identity(authorization)
    prior = identity(root / "operations/s2-publish/receipt.json")
    require(mandate["sha256"] == OUTLET_AUTHORIZATION and prior["sha256"] == S2_SETTLEMENT,
            "outlet amendment authority or prior settlement differs")
    ledger = json.loads(book.read_text())
    if name in ledger["entries"]:
        amendment = ledger["entries"][name]
        require(amendment["authorization"] == mandate and amendment["settlement"] == prior,
                "existing amendment has a different binding")
        if not receipt.exists():
            atomic(receipt, amendment)
        require(json.loads(receipt.read_text()) == amendment, "amendment receipt differs")
        return amendment["after"]
    require(identity(book)["sha256"] == S2_ACCOUNT,
            "live book differs from the published S2 settlement")
    before = balance(ledger)
    require(before == json.loads(Path(prior["path"]).read_text())["balance"]
            and not ledger["closed"] and before["reserved_seconds"] == 0
            and ledger["ceiling"] == 43200 and ledger["protected_seconds"] == 3600,
            "invalid starting balance for outlet amendment")
    require(not receipt.exists(), "unexpected outlet amendment receipt")

    def amend(value):
        require(value == ledger, "ledger changed during outlet amendment")
        value["ceiling"] += 10800
        value["outlet_evidence_floor"] = evidence_spend(value) + 300
        event = dict(status="PASS", phase="accounting", seconds=time.monotonic() - began,
                     authorization=mandate, settlement=prior, opening_ledger=identity(book),
                     prior_ceiling=43200, added_seconds=10800, ceiling=value["ceiling"],
                     protected_seconds=value["protected_seconds"], later_publication_hold=300,
                     before=before)
        value["entries"][name] = event
        event["after"] = balance(value)

    ledger = update(book, amend)
    atomic(receipt, ledger["entries"][name])
    return ledger["entries"][name]["after"]


def memory_snapshot(pid="self"):
    mem = {line.split(":")[0]: line.split()[1] for line in
           Path("/proc/meminfo").read_text().splitlines()}
    member, = [line[3:] for line in Path(f"/proc/{pid}/cgroup").read_text().splitlines()
               if line.startswith("0::")]
    mount, = [line.split()[3:5] for line in Path("/proc/self/mountinfo").read_text().splitlines()
              if " - cgroup2 " in line]
    hierarchy, mountpoint = map(Path, mount)
    group = mountpoint / Path(member).relative_to(hierarchy)
    ancestors = []
    while True:
        require(group.is_dir(), f"cannot inspect cgroup {group}")
        usage = group / "memory.current"
        tracked = {"path": str(group)}
        if usage.exists():
            tracked["current"] = int(usage.read_text())
        for limit in ("memory.max", "memory.high"):
            path = group / limit
            value = path.read_text().strip() if path.exists() else "max"
            tracked[limit] = int(value) if value != "max" else None
            require(tracked[limit] is None or "current" in tracked,
                    f"cannot assess {limit} headroom in {group}")
        ancestors.append(tracked)
        if group == mountpoint:
            break
        group = group.parent
    return dict(mem_available_bytes=int(mem["MemAvailable"]) * 1024,
                mem_total_bytes=int(mem["MemTotal"]) * 1024, job_cgroup=member,
                ancestor_limits=ancestors)


def memory_headroom(snapshot):
    return min([snapshot["mem_available_bytes"]] + [
        group[limit] - group["current"] for group in snapshot["ancestor_limits"]
        for limit in ("memory.max", "memory.high") if group[limit] is not None])


def outlet_solve(root, ledger, label, command):
    require("outlet-cp0-amendment" in ledger["entries"], "missing outlet authorization")
    match = re.fullmatch(r"(outlet-(?:b0|b1|a1c|b1c|b2|bx))(?:-continue-[1-9][0-9]*)?-steady", label)
    require(match and len(command) == 7 and command[:3] == ["mpirun", "-np", "1"],
            "RSS amendment is limited to named single-rank outlet solves")
    runtime = root / "runtime-outlet-v1/IncNavierStokesSolver"
    session = Path(command[-1])
    require(Path(command[3]) == runtime and command[4:6] == ["--use-opt-file", str(session) + ".opt"] and
            session == root / "science" / label.removesuffix("-steady") / "run.xml",
            "outlet solver or session outside admitted path")
    conditions = ET.parse(session).getroot().find("CONDITIONS")
    info = {i.attrib["PROPERTY"]: i.attrib["VALUE"] for i in conditions.findall("SOLVERINFO/I")}
    parameters = dict(tuple(v.strip() for v in p.text.split("=", 1))
                      for p in conditions.findall("PARAMETERS/P"))
    policy = "FullTractionV1" if match[1] == "outlet-a1c" else "ZeroTransverseV1"
    require(info["WJNativeMode"] == "Physical" and info["WJNativeSolve"] == "SteadyReference"
            and info["WJReferencePolicy"] == "ExactInnerFiniteNuV1"
            and info["WJOutletPolicy"] == policy and float(parameters["NumSteps"]) == 0
            and parameters["WJSteadySeconds"], "session does not match amended diagnostic")
    return memory_snapshot()


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
    if "outlet_evidence_floor" in ledger:
        final_publication = phase == "evidence" and label.startswith("outlet-cp4-publish")
        hold = 0 if final_publication else 300
    else:
        hold = 0 if phase == "evidence" else max(0, ledger["initial_evidence_seconds"] - evidence_spend(ledger))
    cap = args.seconds + 7  # Allow bounded startup, monitoring and signal grace.
    require(remaining["reserved_seconds"] == 0, "unsettled operation exists")
    require(cap + hold <= min(remaining["nonprotected_remaining"], remaining["tranche_remaining"]),
            "operation cap exceeds live tranche or parent allowance")
    operation = root / "operations" / label
    host, helpers = workers()
    # Evidence commands run serially with a buffer for the process census.
    # The shared 72-worker watchdog also applies after admission.
    limits = {"build": 16, "science": 10, "analysis": 6, "evidence": 3, "outlet-science": 10}
    if phase == "build" and label.startswith("outlet-") and command[:2] == ["cmake", "--build"]:
        require(command[2:5] == [str(root / "build-outlet-v1"), "--parallel", "4"],
                "outlet build requires the four-job isolated target")
        anchor = json.loads((root / "operations/r2-native-build/receipt.json").read_text())
        require(anchor["request"]["command"][3:5] == ["--parallel", "4"]
                and anchor["result"]["status"] == "PASS"
                and anchor["result"]["peak_host_workers"] - anchor["request"]["workers_before"] <= 11,
                "outlet build worker anchor differs")
        limits[phase] = 14  # Previous four-job Make build peaked 11 above admission.
    if phase == "build" and label == "outlet-configure" and command[:2] == ["cmake", "-S"]:
        limits[phase] = 8  # Serial generator and compiler-detection probes.
    memory = outlet_solve(root, ledger, label, command) if phase == "outlet-science" else None
    rss_limit = (16 if phase == "build" else 8) * GIB
    if memory:
        available = memory_headroom(memory)
        require(available >= 10 * GIB, "outlet diagnostic memory pressure at admission")
        rss_limit = (10 if available >= 12 * GIB else 8) * GIB
    require(host + limits[phase] <= 72, "shared worker limit at admission")
    info = storage(root)
    request = dict(label=label, phase=phase, command=command, cap_seconds=args.seconds,
                   reserved_seconds=cap, workers_before=host, idle_mpi_event_helpers=helpers,
                   threads=THREADS, storage_before=info, runner=identity(__file__),
                   rss_limit_bytes=rss_limit, memory_before=memory,
                   input_files=[identity(p) for p in command
                                if p.startswith("/") and "\n" not in p and Path(p).is_file()])
    operation.mkdir()
    atomic(operation / "request.json", request)
    update(book, lambda ledger: ledger["entries"].update(
        {label: {"status": "RESERVED", "phase": phase, "reserved_seconds": cap}}))
    peak, peak_workers, code, reason, groups = 0, host, None, None, set()
    job_memory, lowest_headroom = None, None
    env = dict(os.environ, **THREADS, PYTHONDONTWRITEBYTECODE="1",
               TMPDIR=str(root / "tmp"), H27_ROOT=str(root), H27_OPERATION=str(operation))
    log_path = operation / "command.log"
    with log_path.open("x") as log:
        try:
            check_storage(info)
            with subprocess.Popen(command, cwd=root, env=env, stdout=log,
                                  stderr=subprocess.STDOUT, start_new_session=True) as child:
                try:
                    while child.poll() is None:
                        rss, child_groups = tree(child.pid)
                        groups |= child_groups
                        peak = max(peak, rss)
                        peak_workers = max(peak_workers, workers()[0])
                        require(peak <= rss_limit, "RSS_LIMIT")
                        if memory:
                            current = memory_snapshot()
                            if job_memory is None:
                                job_memory = memory_snapshot(child.pid)
                                require(current["job_cgroup"] == job_memory["job_cgroup"],
                                        "solver moved to a different memory cgroup")
                            available = memory_headroom(current)
                            lowest_headroom = min(available, lowest_headroom if lowest_headroom is not None else available)
                            require(available >= 2 * GIB + max(0, rss_limit - rss),
                                    "outlet diagnostic memory pressure during run")
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
                 peak_tree_rss_bytes=peak, peak_host_workers=peak_workers, log=identity(log_path),
                 rss_limit_bytes=rss_limit, job_memory=job_memory,
                 minimum_memory_headroom_bytes=lowest_headroom)
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
    amendment = modes.add_parser("amend-outlet")
    amendment.add_argument("authorization", type=Path)
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
        elif args.action == "amend-outlet":
            amend_outlet(root, book, args.authorization, began)
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
