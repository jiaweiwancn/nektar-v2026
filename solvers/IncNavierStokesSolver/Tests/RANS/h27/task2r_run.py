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
          "outlet-science": 7200, "r3-advance": 7200, "r3-audit": 3600,
          "r3-stock": 3600, "r3-restart": 7200, "r3-parity": 900}
OUTLET_AUTHORIZATION = "9356a31fab90117bf05438281dcf6deef4805b7c18f3323df088e684c94c0037"
S2_SETTLEMENT = "f2ae3f2c814d22dfa4161e976cd1792e2968b7303cbbe597186421383b52b64e"
S2_ACCOUNT = "dec44571419200acdb94c9e63c94c6a9a490f78caff158e0b31f92e5b933a784"
OUTLET_WORKER_AUTHORIZATION = (
    "You are now allowed to employ 8 workers at most in addition to the current 72 worker limit. Resume"
)
OUTLET_QUAL_AUTHORIZATION = "52982518e25a333b9d4f25d7917be24a0f7963cc3f01285ed7fe47885b60cf01"
OUTLET_CP4_SETTLEMENT = "f1c3c623e6ce17262872b313d132b7305dd585fbdc3580bfe9499ac31b0dcc45"
CAMPAIGN_AUTHORIZATION = "b1595fa9c7cb3ded49accdb6f126b76e75769db52d31b82b9ce718d84786fd19"
CAMPAIGN_SETTLEMENT = "4a40999e5914643c0639b584ad7d432aa8c1c729a61b2dd138842e493c2061da"
CAMPAIGN_AMENDMENT = "outlet-campaign-amendment"
CAMPAIGN_CASE_PATTERN = r"B(?:fl|(?:r|x|h|c|q|t|p|e)[1-9][0-9]*)"
CAMPAIGN_RSS_GIB = 20
EXTENT_RUNTIME = "runtime-outlet-campaign-extent"
EXTENT_BINDING = "outlet-campaign-extent-runtime.json"
QUAL_CASES = {
    "B3": "b3", "B4": "b4", "Bn": "bn", "Bn2": "bn2", "Bf": "bf",
    "Bw": "bw", "Bw2": "bw2", "Btip": "btip",
    "Bq": "bq", "Bt": "bt",
}
WORK = ("simpleFoam", "IncNavierStokes", "RANS", "WJ", "python", "cc1",
        "c++", "g++", "gcc", "clang", "gfortran", "f951", "lto1",
        "orted", "mpirun", "mpiexec", "ninja", "make", "cmake")


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def campaign_case(case):
    return re.fullmatch(CAMPAIGN_CASE_PATTERN, case) is not None


def extent_case(case):
    return campaign_case(case) and case[1] in ("p", "e")


def qualification_slug(case):
    if campaign_case(case):
        return case.lower()
    return QUAL_CASES.get(case)


def worker_census(owned_groups=()):
    rows = subprocess.check_output(
        ["ps", "-eLo", "pid=,tid=,pgid=,state=,wchan=,comm="], text=True)
    counted = helpers = owned = stopped = 0
    for row in rows.splitlines():
        pid, tid, group, state, wait, comm = row.split(None, 5)
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
            counted += 1
            if owned_groups and int(group) in owned_groups:
                owned += 1
                stopped += state in ("T", "t")
    return counted, helpers, owned, stopped


def workers():
    return worker_census()[:2]


def storage(root, advancing=False):
    info = {
        "science_bytes": size(root / "science") + size(root / "tmp"),
        "source_build_runtime_bytes": sum(size(p) for p in root.iterdir()
                                          if p.is_dir() and p.name.startswith(("source", "build", "runtime"))),
        "free_bytes": shutil.disk_usage(root).free,
    }
    if advancing:
        info["r3_science_bytes"] = size(root / "science-r3") + size(root / "tmp-r3")
    return info


def check_storage(info):
    require(info["free_bytes"] >= 15 * GIB, "free-space limit")
    require(info["science_bytes"] <= 8 * GIB, "science storage limit")
    require(info["source_build_runtime_bytes"] <= 8 * GIB, "build storage limit")
    if "r3_science_bytes" in info:
        require(info["r3_science_bytes"] <= 8 * GIB
                and info["r3_science_bytes"] + info["science_bytes"] <= 16 * GIB,
                "R3 science storage limit")


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


def amend_outlet_workers(root, book, authorization, began):
    name = "outlet-worker-amendment"
    receipt = root / "evidence" / (name + ".json")
    require(authorization.read_text().strip() == OUTLET_WORKER_AUTHORIZATION,
            "outlet worker authorization differs")
    mandate = identity(authorization)
    ledger = json.loads(book.read_text())
    require("outlet-cp0-amendment" in ledger["entries"], "missing outlet resource account")
    if name in ledger["entries"]:
        event = ledger["entries"][name]
        require(event["authorization"] == mandate, "existing worker amendment binding differs")
    else:
        before = balance(ledger)
        require(not ledger["closed"] and before["reserved_seconds"] == 0,
                "worker amendment requires a settled outlet account")
        require(not receipt.exists(), "unexpected worker amendment receipt")

        def amend(value):
            require(value == ledger, "ledger changed during outlet worker amendment")
            opening_ledger = identity(book)
            host = workers()[0]
            value["entries"][name] = event = dict(
                status="PASS", phase="accounting", seconds=time.monotonic()-began,
                authorization=mandate, opening_ledger=opening_ledger,
                prior_worker_limit=72, added_workers=8, worker_limit=80,
                ceiling=value["ceiling"], protected_seconds=value["protected_seconds"],
                before=before, workers_at_amendment=host)
            event["after"] = balance(value)

        event = update(book, amend)["entries"][name]
    if not receipt.exists():
        atomic(receipt, event)
    require(json.loads(receipt.read_text()) == event, "worker amendment receipt differs")


def amend_outlet_qualification(root, book, authorization, began):
    name = "outlet-qual-amendment"
    receipt = root / "evidence" / (name + ".json")
    mandate = identity(authorization)
    prior = root / "operations/outlet-cp4-publish/receipt.json"
    settlement = identity(prior)
    require(mandate["sha256"] == OUTLET_QUAL_AUTHORIZATION
            and settlement["sha256"] == OUTLET_CP4_SETTLEMENT,
            "flat qualification authorization or CP4 settlement differs")
    ledger = json.loads(book.read_text())
    old = json.loads(prior.read_text())
    require(ledger["entries"]["outlet-cp4-publish"] == old["result"]
            and ledger["entries"]["outlet-worker-amendment"]["worker_limit"] == 80,
            "CP4 settlement or shared-host allowance differs")
    if name in ledger["entries"]:
        event = ledger["entries"][name]
        require(ledger["ceiling"] == 64800 and event["authorization"] == mandate
                and event["settlement"] == settlement, "existing qualification amendment differs")
    else:
        before = balance(ledger)
        require(not ledger["closed"] and before["reserved_seconds"] == 0
                and ledger["ceiling"] == 54000 and ledger["protected_seconds"] == 3600
                and before["charged_seconds"] >= old["balance"]["charged_seconds"],
                "qualification amendment requires a reconciled CP4 account")
        require(not receipt.exists(), "unexpected qualification amendment receipt")

        def amend(value):
            require(value == ledger, "ledger changed during qualification amendment")
            opening = identity(book)
            value["ceiling"] += 10800
            value["entries"][name] = event = dict(
                status="PASS", phase="accounting", seconds=time.monotonic() - began,
                authorization=mandate, settlement=settlement, opening_ledger=opening,
                prior_ceiling=54000, added_seconds=10800, ceiling=value["ceiling"],
                protected_seconds=value["protected_seconds"], later_publication_hold=300,
                before=before)
            event["after"] = balance(value)

        event = update(book, amend)["entries"][name]
    if not receipt.exists():
        atomic(receipt, event)
    require(json.loads(receipt.read_text()) == event, "qualification amendment receipt differs")


def amend_outlet_campaign(root, book, authorization, began):
    receipt = root / "evidence" / f"{CAMPAIGN_AMENDMENT}.json"
    mandate = identity(authorization)
    settlement = identity(root / "operations/outlet-qual-cp3-publish/receipt.json")
    require(mandate["sha256"] == CAMPAIGN_AUTHORIZATION
            and settlement["sha256"] == CAMPAIGN_SETTLEMENT,
            "campaign authorization or prior settlement differs")
    ledger = json.loads(book.read_text())
    settled = json.loads(Path(settlement["path"]).read_text())
    require(ledger["entries"]["outlet-qual-cp3-publish"] == settled["result"]
            and ledger["entries"]["outlet-worker-amendment"]["worker_limit"] == 80
            and ledger["protected_seconds"] == 3600, "campaign resource anchors differ")
    if CAMPAIGN_AMENDMENT in ledger["entries"]:
        event = ledger["entries"][CAMPAIGN_AMENDMENT]
        require(ledger["ceiling"] == 108000 and event["authorization"] == mandate
                and event["settlement"] == settlement, "existing campaign amendment differs")
    else:
        before = balance(ledger)
        require(not ledger["closed"] and before["reserved_seconds"] == 0
                and ledger["ceiling"] == 64800
                and before["charged_seconds"] >= settled["balance"]["charged_seconds"],
                "campaign amendment requires the reconciled prior account")
        require(not receipt.exists(), "unexpected campaign amendment receipt")

        def amend(value):
            require(value == ledger, "ledger changed during campaign amendment")
            opening = identity(book)
            value["ceiling"] += 43200
            value["entries"][CAMPAIGN_AMENDMENT] = event = dict(
                status="PASS", phase="accounting", seconds=time.monotonic()-began,
                authorization=mandate, settlement=settlement, opening_ledger=opening,
                prior_ceiling=64800, added_seconds=43200, ceiling=value["ceiling"],
                protected_seconds=value["protected_seconds"], later_publication_hold=300,
                worker_limit=80, ordinary_rss_gib=10, named_maximum_rss_gib=CAMPAIGN_RSS_GIB,
                before=before)
            event["after"] = balance(value)

        event = update(book, amend)["entries"][CAMPAIGN_AMENDMENT]
    if not receipt.exists():
        atomic(receipt, event)
    require(json.loads(receipt.read_text()) == event, "campaign amendment receipt differs")


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


def qualification_solve(root, ledger, label, cap, runtime, session, info, parameters):
    require("outlet-qual-amendment" in ledger["entries"], "missing qualification authorization")
    manifest = json.loads((session.parent / "inputs.json").read_text())
    case = manifest["case"]
    slug = qualification_slug(case)
    require(slug is not None, "unknown qualification case")
    base = f"outlet-qual-{slug}"
    require(re.fullmatch(base + r"(?:-continue-[1-9][0-9]*)?-steady", label),
            "qualification label differs from admitted case")
    admission_path = Path(manifest["admission"]["path"])
    require(admission_path.parent == root / "evidence" and
            re.fullmatch(base + r"-admission(?:-[1-9][0-9]*)?\.json", admission_path.name),
            "qualification forecast is outside the admitted evidence path")
    admission = json.loads(admission_path.read_text())
    attempt = admission["attempts"][label.removesuffix("-steady")]
    require(admission["schema"] == "WJH27FlatQualificationAdmissionV1"
            and admission["status"] == "ADMITTED" and admission["case"] == manifest["case"]
            and identity(admission_path) == manifest["admission"]
            and cap == attempt["watch_seconds"], "qualification attempt or cap differs from forecast")
    require(admission["quadrature"] == manifest["quadrature"] ==
                [int(parameters["WJQuadratureX"]), int(parameters["WJQuadratureR"])]
            and admission["tolerances"] == manifest["tolerances"] ==
                [float(parameters["WJLinearTolerance"]), float(parameters["WJCoupledTolerance"])]
            and attempt["native_seconds"] == float(parameters["WJSteadySeconds"])
            and attempt["newton_iterations"] == int(parameters["WJSteadyNewtonIterations"])
            and info["WJConstitutiveBackend"] == "FastFull3DV1"
            and info["WJStressSplit"] == "PositiveShearV1"
            and info["WJSteadyVelocityPC"] == "FullGradientV1"
            and int(parameters["WJSteadyKrylovDimension"]) == 64
            and float(parameters["Time"]) == 0. and float(parameters["Kinvis"]) == 6.423173803526447e-5
            and "WJRestart" not in info, "qualification changed the declared numerical recipe")
    seed_source = manifest.get("seed_source", manifest["source"])
    if extent_case(case):
        proof_path = root/"evidence"/EXTENT_BINDING
        ext = json.loads(proof_path.read_text())
        outlet_x = manifest["x_edges"][-1]
        require(admission["extent_runtime"] == identity(proof_path)
                and manifest["application"] == ext["application"] and manifest["source"] == ext["source"]
                and outlet_x == admission["outlet_x"] == float(parameters["WJDiagnosticOutletX"])
                and seed_source in (admission["seed_source"],manifest["source"]),
                "extent runtime, geometry or typed seed scope differs")
    else:
        require(seed_source == manifest["source"] and float(parameters.get("WJDiagnosticOutletX",15.))==15.,
                "same-domain qualification requires its pinned source and exit")
    require(manifest["session"] == identity(session)
            and manifest["options"] == identity(Path(str(session) + ".opt"))
            and manifest["application"] == admission["application"] == identity(runtime)
            and info["WJInitialFields"] == manifest["initial_guess"]["path"]
            and manifest["initial_guess"] == identity(Path(info["WJInitialFields"]))
            and info["WJInitialSourceIdentity"] == seed_source
            and manifest["source"] == admission["source"]
            and [manifest["x_edges"], manifest["physical_z_edges"]] ==
                [admission["x_edges"], admission["physical_z_edges"]],
            "qualification executable, axes or typed seed binding differs")
    limit = admission["rss_limit_bytes"]
    estimate = admission["forecast_peak_tree_rss_bytes"]
    if campaign_case(case):
        event = ledger["entries"].get(CAMPAIGN_AMENDMENT, {})
        require(event.get("status") == "PASS" and ledger["ceiling"] == 108000
                and admission["campaign_amendment"] == identity(root / "evidence" / f"{CAMPAIGN_AMENDMENT}.json")
                and 10 * GIB <= limit <= event["named_maximum_rss_gib"] * GIB
                and limit % GIB == 0 and limit == max(10 * GIB, ((estimate + GIB-1)//GIB)*GIB)
                and admission["anchors"], "campaign RSS limit requires its measured named admission")
    else:
        require(limit in (10 * GIB, 12 * GIB), "historical qualification RSS ceiling differs")
    require(0 < estimate <= limit and (limit == 10 * GIB or estimate > 10 * GIB),
            "RSS exception requires an individually measured forecast")
    return dict(forecast=manifest["admission"], rss_limit_bytes=limit)


def outlet_solve(root, ledger, label, command, cap):
    require("outlet-cp0-amendment" in ledger["entries"], "missing outlet authorization")
    slugs = "|".join([*QUAL_CASES.values(), CAMPAIGN_CASE_PATTERN.lower()])
    match = re.fullmatch(r"(outlet-(?:b0|b1|a1c|b1c|b2|bx|qual-(?:" + slugs +
                         r")))(?:-continue-[1-9][0-9]*)?-steady", label)
    require(match and len(command) == 7 and command[:3] == ["mpirun", "-np", "1"],
            "RSS amendment is limited to named single-rank outlet solves")
    mode = EXTENT_RUNTIME if re.fullmatch(r"outlet-qual-b[pe][1-9][0-9]*",match[1]) else "runtime-outlet-v1"
    runtime = root / mode / "IncNavierStokesSolver"
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
            and info["WJReferenceTable"] == str(root / "science/r1-qualified/tables-exact-inner-v1.json")
            and info["WJOutputPrefix"] == str(session.parent / "state")
            and info["WJOutletPolicy"] == policy and float(parameters["NumSteps"]) == 0
            and float(parameters["WJSteadySeconds"]) > 0, "session does not match amended diagnostic")
    qualified = (qualification_solve(root, ledger, label, cap, runtime, session, info, parameters)
                 if label.startswith("outlet-qual-") else None)
    return memory_snapshot(), qualified


def choose_rank_cpu():
    allowed = os.sched_getaffinity(0)
    base = min(allowed)
    node, = sorted((Path("/sys/devices/system/cpu") / f"cpu{base}").glob("node[0-9]*"))
    candidates = [cpu for cpu in allowed if
                  (Path("/sys/devices/system/cpu") / f"cpu{cpu}" / node.name).exists()]

    def snapshot():
        ticks = {}
        for line in Path("/proc/stat").read_text().splitlines():
            name, *columns = line.split()
            if name.startswith("cpu") and name[3:].isdigit():
                times = list(map(int, columns[:8]))
                ticks[int(name[3:])] = sum(times), times[3]
        return ticks

    before = snapshot()
    time.sleep(1.)
    after = snapshot()
    idle = {cpu: (after[cpu][1] - before[cpu][1]) / max(1, after[cpu][0] - before[cpu][0])
            for cpu in candidates}
    cpu = max(candidates, key=lambda i: (idle[i], -i))
    require(idle[cpu] >= .9, "no quiet permitted CPU on the launcher's NUMA node")
    return dict(cpu=cpu, reference_cpu=base, numa_node=node.name,
                sampled_idle_fraction=idle[cpu], sample_seconds=1.)


def observe_rank_affinity(groups, executable, cpu, observed):
    candidates = groups if observed is None else (observed["pid"],)
    for pid in candidates:
        try:
            if Path(f"/proc/{pid}/exe").readlink() != executable:
                continue
            cpus = sorted(os.sched_getaffinity(pid))
            if observed is None:
                observed = dict(pid=pid, startup_cpus=cpus, allowed_cpus=cpus)
            if cpus != [cpu]:
                require("rebound_from_cpus" not in observed, "native rank affinity differs after one recovery")
                os.sched_setaffinity(pid, {cpu})
                actual = sorted(os.sched_getaffinity(pid))
                observed.update(rebound_from_cpus=cpus, allowed_cpus=actual)
                require(actual == [cpu], "native rank affinity differs")
            return observed
        except (FileNotFoundError, ProcessLookupError):
            continue
    return observed


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


def signal_groups(groups, signal_number):
    for group in groups:
        try:
            os.killpg(group, signal_number)
        except ProcessLookupError:
            pass


def terminate(child, groups, suspended=False):
    owned = groups | {child.pid}
    signal_groups(owned, signal.SIGTERM)
    if suspended:
        signal_groups(owned, signal.SIGCONT)
    try:
        return child.wait(timeout=4)
    except subprocess.TimeoutExpired:
        signal_groups(owned, signal.SIGKILL)
        return child.wait()


def run(root, book, args, began):
    phase, label = args.phase, args.label
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    require(command and label and label.replace("-", "").isalnum(), "invalid command or label")
    require(0 < args.seconds <= LIMITS[phase], "invalid duration cap")
    ledger = json.loads(book.read_text())
    advancing = label.startswith("r3-")
    if advancing:
        require(ledger["entries"].get("r3-scope", {}).get("status") == "PASS",
                "missing R3 resource admission")
    require(advancing or not phase.startswith("r3-"), "R3 operation needs its own label")
    require(not ledger["closed"] and label not in ledger["entries"], "closed account or repeated label")
    remaining = balance(ledger)
    if "outlet_evidence_floor" in ledger:
        final_publication = phase == "evidence" and label.startswith(
            ("outlet-cp4-publish", "outlet-qual-cp3-publish", "outlet-campaign-final-publish"))
        hold = 0 if final_publication or (advancing and phase == "evidence") else 300
    else:
        hold = 0 if phase == "evidence" else max(0, ledger["initial_evidence_seconds"] - evidence_spend(ledger))
    cap = args.seconds + 7  # Allow bounded startup, monitoring and signal grace.
    require(remaining["reserved_seconds"] == 0, "unsettled operation exists")
    require(cap + hold <= min(remaining["nonprotected_remaining"], remaining["tranche_remaining"]),
            "operation cap exceeds live tranche or parent allowance")
    operation = root / "operations" / label
    host, helpers = workers()
    amendment = (ledger["entries"].get("outlet-worker-amendment")
                 if advancing or label.startswith("outlet-") else None)
    worker_limit = 72 if amendment is None else amendment["worker_limit"]
    require(worker_limit == 72 or (worker_limit == 80 and amendment["status"] == "PASS"),
            "invalid worker ceiling")
    limits = {"build": 16, "science": 10, "analysis": 6, "evidence": 3, "outlet-science": 10}
    if phase == "build" and label.startswith("outlet-") and command[:2] == ["cmake", "--build"]:
        build_path = ("build-outlet-campaign-extent" if label=="outlet-campaign-extent-build"
                      else "build-outlet-v1")
        require(command[2:5] == [str(root / build_path), "--parallel", "4"],
                "outlet build requires the four-job isolated target")
        anchor = json.loads((root / "operations/r2-native-build/receipt.json").read_text())
        require(anchor["request"]["command"][3:5] == ["--parallel", "4"]
                and anchor["result"]["status"] == "PASS"
                and anchor["result"]["peak_host_workers"] - anchor["request"]["workers_before"] <= 11,
                "outlet build worker anchor differs")
        limits[phase] = 14  # Previous four-job Make build peaked 11 above admission.
    if phase == "build" and label in ("outlet-configure","outlet-campaign-extent-configure") and command[:2] == ["cmake", "-S"]:
        limits[phase] = 8  # Serial generator and compiler-detection probes.
    memory, qualification = (outlet_solve(root, ledger, label, command, args.seconds)
                             if phase == "outlet-science" else (None, None))
    rss_limit = (16 if phase == "build" else 8) * GIB
    if memory:
        available = memory_headroom(memory)
        require(available >= 10 * GIB, "outlet diagnostic memory pressure at admission")
        rss_limit = (10 if available >= 12 * GIB else 8) * GIB
    pause_policy = None
    affinity = None
    r3_admission = None
    if advancing:
        from task2r_advancing import admission
        if phase.startswith("r3-"):
            r3_admission = admission(root, label, phase, command, args.seconds)
            memory = memory_snapshot()
            rss_limit = r3_admission["rss_limit_bytes"]
            require(memory_headroom(memory) >= rss_limit + 2 * GIB, "R3 memory headroom")
            if phase != "r3-audit":
                cpu = r3_admission["cpu"]
                require(cpu in os.sched_getaffinity(0), "R3 pinned CPU unavailable")
                affinity = dict(cpu=cpu, admission=r3_admission["admission"])
                pause_policy = dict(trigger_workers=78, resume_workers=76, quiet_seconds=4.,
                                    maximum_paused_seconds=55., count_stopped_ownership=True)
            limits[phase] = 8
        elif phase == "build":
            if command[:2] == ["cmake", "--build"]:
                require(command[2:5] == [str(root / "build-r3"), "--parallel", "4"],
                        "R3 build requires four workers and its isolated build tree")
                limits[phase] = 14
            else:
                require(command[:2] in (["cmake", "-S"], ["cmake", "-E"]),
                        "unexpected R3 build command")
                limits[phase] = 8
        else:
            require(phase in ("analysis", "evidence"), "unexpected R3 phase")
    elif qualification:
        require(worker_limit == 80, "qualification needs the recorded shared-host allowance")
        rss_limit = qualification["rss_limit_bytes"]
        require(memory_headroom(memory) >= rss_limit + 2 * GIB, "qualification memory headroom")
        affinity = choose_rank_cpu()
        host, helpers = workers()
        limits[phase] = 8
        pause_policy = dict(trigger_workers=78, resume_workers=76, quiet_seconds=4.,
                            maximum_paused_seconds=55., admission=qualification["forecast"],
                            count_stopped_ownership=True)
    elif phase == "outlet-science" and worker_limit == 80 and label.startswith("outlet-b2-continue-"):
        admission = root / "evidence/outlet-b2-worker80-admission.json"
        require(label in json.loads(admission.read_text())["steady_labels"],
                "missing scoped worker rescheduling record")
        anchor = root / "operations/outlet-b1c-continue-1-steady/receipt.json"
        prior = json.loads(anchor.read_text())
        require(prior["result"]["status"] == "PASS"
                and prior["request"]["command"][:3] == command[:3]
                and prior["result"]["peak_host_workers"] - prior["request"]["workers_before"] <= 6,
                "paired single-rank worker anchor differs")
        limits[phase] = 8
        pause_policy = dict(trigger_workers=79, resume_workers=77, quiet_seconds=4.,
                            maximum_paused_seconds=55., admission=identity(admission),
                            count_stopped_ownership=True, worker_anchor=identity(anchor))
    elif phase == "outlet-science" and label in ("outlet-b2-continue-2-steady", "outlet-b2-continue-3-steady"):
        admission = root / "evidence/outlet-b2-final-admission.json"
        require(json.loads(admission.read_text())["worker_pause_admitted"] is True,
                "missing scoped worker rescheduling record")
        pause_policy = dict(trigger_workers=68, resume_workers=66, quiet_seconds=5.,
                            maximum_paused_seconds=55., admission=identity(admission))
    require(host + limits[phase] <= worker_limit, "shared worker limit at admission")
    info = storage(root, advancing)
    if r3_admission:
        require(info["r3_science_bytes"] + r3_admission["peak_new_bytes"] <= 8 * GIB,
                "R3 peak overlap exceeds remaining storage")
    peak_r3_bytes = info.get("r3_science_bytes", 0)
    request = dict(label=label, phase=phase, command=command, cap_seconds=args.seconds,
                   reserved_seconds=cap, workers_before=host, idle_mpi_event_helpers=helpers,
                   threads=THREADS, storage_before=info, runner=identity(__file__),
                   rss_limit_bytes=rss_limit, memory_before=memory, worker_pause_policy=pause_policy,
                   worker_limit=worker_limit, worker_admission_margin=limits[phase],
                   input_files=[identity(p) for p in command
                                if p.startswith("/") and "\n" not in p and Path(p).is_file()])
    if r3_admission:
        request["r3_admission"] = r3_admission["admission"]
    launched = command
    if affinity:
        # Bind in the rank wrapper, after OpenMPI has placed the launcher child.
        launched = [*command[:3], "--bind-to", "none", "taskset", "--cpu-list",
                    str(affinity["cpu"]), *command[3:]]
        request["launch"] = dict(command=launched, affinity=affinity)
    operation.mkdir()
    atomic(operation / "request.json", request)
    update(book, lambda ledger: ledger["entries"].update(
        {label: {"status": "RESERVED", "phase": phase, "reserved_seconds": cap}}))
    peak, peak_workers, code, reason, groups = 0, host, None, None, set()
    job_memory, lowest_headroom = None, None
    pauses, suspended_at, quiet_at, suspended_seconds = [], None, None, 0.
    owned_peak = 0
    rank_affinity = None
    env = dict(os.environ, **THREADS, PYTHONDONTWRITEBYTECODE="1",
               TMPDIR=str(root / ("tmp-r3" if advancing else "tmp")),
               H27_ROOT=str(root), H27_OPERATION=str(operation))
    log_path = operation / "command.log"
    with log_path.open("x") as log:
        try:
            check_storage(info)
            with subprocess.Popen(launched, cwd=root, env=env, stdout=log,
                                  stderr=subprocess.STDOUT, start_new_session=True) as child:
                try:
                    while child.poll() is None:
                        rss, child_groups = tree(child.pid)
                        if suspended_at is not None:
                            signal_groups(child_groups-groups, signal.SIGSTOP)
                        groups |= child_groups
                        if affinity:
                            rank_affinity = observe_rank_affinity(
                                child_groups, Path(command[3]), affinity["cpu"], rank_affinity)
                        peak = max(peak, rss)
                        if pause_policy and pause_policy.get("count_stopped_ownership"):
                            observed_workers, _, owned, stopped = worker_census(groups | {child.pid})
                            if suspended_at is None:
                                owned_peak = max(owned_peak, owned)
                            projected_workers = observed_workers-stopped+owned_peak
                        else:
                            observed_workers = projected_workers = workers()[0]
                        peak_workers = max(peak_workers, observed_workers)
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
                        now = time.monotonic()
                        if pause_policy:
                            if suspended_at is None and observed_workers >= pause_policy["trigger_workers"]:
                                signal_groups(groups | {child.pid}, signal.SIGSTOP)
                                suspended_at, quiet_at = now, None
                                pauses.append(dict(start_seconds=now-began, workers_at_stop=observed_workers))
                            elif suspended_at is not None:
                                require(suspended_seconds+now-suspended_at < pause_policy["maximum_paused_seconds"],
                                        "HOST_PRESSURE_REMAINS")
                                if projected_workers <= pause_policy["resume_workers"]:
                                    quiet_at = now if quiet_at is None else quiet_at
                                    if now-quiet_at >= pause_policy["quiet_seconds"]:
                                        signal_groups(groups | {child.pid}, signal.SIGCONT)
                                        suspended_seconds += now-suspended_at
                                        pauses[-1].update(end_seconds=now-began, workers_at_resume=observed_workers,
                                                          projected_workers_at_resume=projected_workers)
                                        suspended_at, quiet_at = None, None
                                else:
                                    quiet_at = None
                        require(observed_workers <= worker_limit
                                or (not advancing and suspended_at is not None), "HOST_WORKER_LIMIT")
                        require(time.monotonic() - began < args.seconds, "TIME_LIMIT")
                        current_storage = storage(root, advancing)
                        check_storage(current_storage)
                        peak_r3_bytes = max(peak_r3_bytes, current_storage.get("r3_science_bytes", 0))
                        if r3_admission:
                            require(peak_r3_bytes - info["r3_science_bytes"] <= r3_admission["peak_new_bytes"],
                                    "R3 operation exceeded its forecast disk allowance")
                        require(log_path.stat().st_size <= 32 * 1024 ** 2, "LOG_LIMIT")
                        time.sleep(.5)
                    code = child.returncode
                    if affinity and code == 0:
                        require(rank_affinity is not None and rank_affinity["allowed_cpus"] == [affinity["cpu"]],
                                "verified native rank affinity is missing")
                except (RuntimeError, OSError, KeyboardInterrupt) as error:
                    code, reason = terminate(child, groups, suspended_at is not None), str(error)
        except (RuntimeError, OSError, KeyboardInterrupt) as error:
            reason = str(error)
    output_files = None
    if phase == "outlet-science" and code == 0 and reason is None:
        case = Path(command[-1]).parent
        try:
            output_files = [identity(case / f"state.{name}") for name in
                            ("result.json", "steady.json", "final.bin")]
        except OSError as error:
            reason = f"cannot seal accepted outlet output: {error}"
    if suspended_at is not None:
        now = time.monotonic()
        suspended_seconds += now-suspended_at
        pauses[-1]["settled_seconds"] = now-began
    event = dict(status="PASS" if code == 0 and reason is None else "FAIL", phase=phase,
                 seconds=time.monotonic() - began, exit=code, reason=reason,
                 peak_tree_rss_bytes=peak, peak_host_workers=peak_workers, log=identity(log_path),
                 rss_limit_bytes=rss_limit, job_memory=job_memory,
                 minimum_memory_headroom_bytes=lowest_headroom,
                 worker_limit=worker_limit, peak_owned_workers=owned_peak,
                 worker_pauses=pauses, worker_paused_seconds=suspended_seconds)
    if output_files is not None:
        event["output_files"] = output_files
    if affinity:
        event["rank_affinity"] = rank_affinity
    settled_storage = storage(root, advancing)
    if advancing:
        event["peak_r3_science_bytes"] = max(peak_r3_bytes, settled_storage["r3_science_bytes"])
        try:
            check_storage(settled_storage)
            if r3_admission:
                require(event["peak_r3_science_bytes"] - info["r3_science_bytes"]
                        <= r3_admission["peak_new_bytes"], "R3 operation exceeded its forecast disk allowance")
        except RuntimeError as error:
            event.update(status="FAIL", reason=str(error))
    ledger = update(book, lambda ledger: ledger["entries"].update({label: event}))
    atomic(operation / "receipt.json", dict(request=request, result=event,
                                            balance=balance(ledger), storage=settled_storage))
    print(label, event["status"], event["seconds"], event["reason"], flush=True)
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
    workers_amendment = modes.add_parser("amend-outlet-workers")
    workers_amendment.add_argument("authorization", type=Path)
    qualification_amendment = modes.add_parser("amend-outlet-qualification")
    qualification_amendment.add_argument("authorization", type=Path)
    campaign_amendment = modes.add_parser("amend-outlet-campaign")
    campaign_amendment.add_argument("authorization", type=Path)
    r3 = modes.add_parser("admit-r3-scope")
    r3.add_argument("authorization", type=Path)
    r3_amendment = modes.add_parser("amend-r3-time")
    r3_amendment.add_argument("forecast", type=Path)
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
        elif args.action == "amend-outlet-workers":
            amend_outlet_workers(root, book, args.authorization, began)
        elif args.action == "amend-outlet-qualification":
            amend_outlet_qualification(root, book, args.authorization, began)
        elif args.action == "amend-outlet-campaign":
            amend_outlet_campaign(root, book, args.authorization, began)
        elif args.action in ("admit-r3-scope", "amend-r3-time"):
            from task2r_advancing import amend_time, admit_scope
            if args.action == "admit-r3-scope":
                admit_scope(root, book, args.authorization, began)
            else:
                amend_time(root, book, args.forecast, began)
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
