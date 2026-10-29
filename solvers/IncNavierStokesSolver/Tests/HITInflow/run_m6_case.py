#!/usr/bin/env python3
"""Execute an already prepared M6 case and retain commands/resource evidence.

Run on the authoritative Linux host with Nektar's environment activated.
All output is placed beside the supplied session, outside the source checkout.
"""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
import xml.etree.ElementTree as ET


def process_tree(root):
    rows = subprocess.check_output(
        ["ps", "-eo", "pid=,ppid=,rss=,nlwp="], text=True).splitlines()
    table = [tuple(map(int, row.split())) for row in rows]
    descendants = {root}
    while True:
        added = {pid for pid, ppid, _, _ in table if ppid in descendants}
        if added <= descendants:
            break
        descendants |= added
    selected = [row for row in table if row[0] in descendants]
    return sum(row[2] for row in selected), sum(row[3] for row in selected)


def run(command, directory, stem, environment):
    start = time.monotonic()
    peak_rss_kib = peak_threads = 0
    with (directory / (stem + ".log")).open("w") as log:
        process = subprocess.Popen(command, cwd=directory, env=environment,
                                   stdin=subprocess.DEVNULL, stdout=log,
                                   stderr=subprocess.STDOUT)
        while process.poll() is None:
            rss, threads = process_tree(process.pid)
            peak_rss_kib = max(peak_rss_kib, rss)
            peak_threads = max(peak_threads, threads)
            time.sleep(2)
    result = {"command": command, "exit_code": process.returncode,
              "wall_seconds": time.monotonic() - start,
              "sampled_process_tree_peak_rss_kib": peak_rss_kib,
              "sampled_process_tree_peak_threads": peak_threads,
              "resource_sample_seconds": 2,
              "log_sha256": hashlib.sha256(
                  (directory / (stem + ".log")).read_bytes()).hexdigest()}
    (directory / (stem + ".receipt.json")).write_text(
        json.dumps(result, indent=2) + "\n")
    if process.returncode:
        raise RuntimeError(f"{stem} exited {process.returncode}; inspect its log")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("session", type=Path)
    parser.add_argument("--ranks", type=int, default=16)
    parser.add_argument("--npz", type=int, default=4)
    parser.add_argument("--mean-initial", action="store_true")
    parser.add_argument("--source-revision", required=True)
    args = parser.parse_args()
    if not 1 <= args.ranks <= 64 or args.npz < 1 or args.ranks % args.npz:
        parser.error("require 1..64 ranks, positive npz dividing ranks")
    session = args.session.resolve()
    directory = session.parent
    source = Path(__file__).resolve().parents[4]
    try:
        directory.relative_to(source)
    except ValueError:
        pass
    else:
        parser.error("run directory must be outside the source checkout")
    environment = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                       TMPDIR=str(directory / "tmp"))
    (directory / "tmp").mkdir(exist_ok=True)
    mpi = ["mpirun", "-n", str(args.ranks)]
    parallel = ["--npz", str(args.npz)] if args.npz > 1 else []
    receipt = {"revision": args.source_revision, "ranks": args.ranks,
               "npz": args.npz, "omp_threads": 1, "openblas_threads": 1,
               "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
               "original_session_sha256": hashlib.sha256(session.read_bytes()).hexdigest()}
    tree = ET.parse(session)
    conditions = tree.getroot().find("CONDITIONS")
    functions = [node for node in conditions.findall("FUNCTION")
                 if node.get("NAME", "").lower() == "initialconditions"]
    for node in functions:
        conditions.remove(node)
    initial_function = ET.SubElement(conditions, "FUNCTION", NAME="InitialConditions")
    for variable, value in [("u", "Ubar"), ("v", "0"), ("w", "0"), ("p", "0")]:
        ET.SubElement(initial_function, "E", VAR=variable, VALUE=value)
    initial_session = directory / "initialization.xml"
    ET.indent(tree, space="  ")
    tree.write(initial_session, encoding="utf-8", xml_declaration=True)
    initial_field = directory / "initial.fld"
    if not args.mean_initial:
        if initial_field.exists():
            parser.error("initial.fld already exists; use a fresh run directory")
        receipt["initialization"] = run(
            mpi + ["HITInitialField"] + parallel + [initial_session.name, initial_field.name],
            directory, "initialization", environment)
        conditions.remove(initial_function)
        initial_function = ET.SubElement(conditions, "FUNCTION", NAME="InitialConditions")
        ET.SubElement(initial_function, "F", VAR="u,v,w,p", FILE=initial_field.name)
    ET.indent(tree, space="  ")
    tree.write(session, encoding="utf-8", xml_declaration=True)
    receipt["executed_session_sha256"] = hashlib.sha256(session.read_bytes()).hexdigest()
    receipt["solver"] = run(mpi + ["IncNavierStokesSolver"] + parallel + [session.name],
                             directory, "solver", environment)
    receipt["finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    (directory / "run.receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2), flush=True)


if __name__ == "__main__":
    main()
