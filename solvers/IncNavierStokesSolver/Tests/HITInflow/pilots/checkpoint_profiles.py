#!/usr/bin/env python3
"""Run native slab diagnostics on immutable checkpoints in separate directories.

Use a completed checkpoint behind the live solver time, or wait for the solver
to finish. The source session and checkpoint are read only. A before/after
digest check rejects changing input, and every output directory must be new.
"""

import argparse
import csv
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time
import xml.etree.ElementTree as ET


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return {"sha256": result.hexdigest(), "bytes": path.stat().st_size}


def checkpoint_inputs(path):
    info = path / "Info.xml"
    document = ET.parse(info).getroot()
    files = [info] + [path / element.attrib["FileName"]
                      for element in document if "FileName" in element.attrib]
    if len(files) < 2 or len(set(files)) != len(files):
        raise ValueError("missing or duplicate checkpoint partitions")
    metadata = document.find("Metadata")
    field_time = float(metadata.findtext("Time"))
    if not math.isfinite(field_time):
        raise ValueError("nonfinite checkpoint time")
    return ({str(p.resolve()): digest(p) for p in files}, field_time,
            metadata.findtext("Provenance/GitSHA1"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("--checkpoints", type=int, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cpu", type=int, default=7)
    args = parser.parse_args()
    if args.cpu not in os.sched_getaffinity(0):
        parser.error("requested CPU is outside the inherited affinity set")
    if len(set(args.checkpoints)) != len(args.checkpoints) or min(args.checkpoints) < 0:
        parser.error("checkpoint indices must be distinct and nonnegative")
    os.sched_setaffinity(0, {args.cpu})
    os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    case, output = args.case.resolve(), args.output.resolve()
    if output == case or case in output.parents:
        parser.error("output must be outside the live case directory")
    manifest = json.loads((case / "manifest.json").read_text())
    if manifest["schema"] != "NEKTAR_HIT_PILOT_V1":
        parser.error("unsupported pilot manifest")
    session = Path(manifest["case"])
    if session.name != str(session):
        parser.error("session must be a local filename")
    executable = shutil.which("HITFlowDiagnostics")
    if executable is None:
        parser.error("activate the Nektar++ installation environment first")
    executable = str(Path(executable).resolve())
    source_files = [case / session, case / "manifest.json", *sorted(case.glob("*.dat"))]
    source_hashes = {str(p): digest(p) for p in source_files}
    output.mkdir(parents=True, exist_ok=False)
    receipt = {"schema": "NEKTAR_HIT_CHECKPOINT_PROFILES_V1", "status": "running",
               "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
               "case": str(case), "script": {"path": str(Path(__file__).resolve()),
                                               **digest(Path(__file__))},
               "executable": {"path": executable, **digest(Path(executable))},
               "cpu_affinity": [args.cpu], "OMP_NUM_THREADS": 1,
               "OPENBLAS_NUM_THREADS": 1, "source_inputs": source_hashes,
               "energy_definition": "Instantaneous energy about prescribed Ubar; not temporal-central TKE",
               "profiles": []}
    for index in args.checkpoints:
        field = case / f"{session.stem}_{index}.chk"
        later = case / f"{session.stem}_{index + 1}.chk" / "Info.xml"
        run_path = case / "run.receipt.json"
        completed = False
        if run_path.is_file():
            run = json.loads(run_path.read_text())
            completed = bool(run.get("finished_utc")) and run["solver"]["exit_code"] == 0
        if not later.is_file() and not completed:
            raise ValueError("require a later checkpoint or a completed successful solver receipt")
        before, field_time, build_revision = checkpoint_inputs(field)
        directory = output / f"checkpoint-{index:03d}"
        directory.mkdir()
        for source in source_files:
            if source.name != "manifest.json":
                shutil.copyfile(source, directory / source.name)
        command = [executable, session.name, str(field)]
        start = time.monotonic()
        with (directory / "diagnostics.log").open("w") as log:
            subprocess.run(command, cwd=directory, stdin=subprocess.DEVNULL,
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        if checkpoint_inputs(field)[0] != before:
            raise ValueError(f"checkpoint changed during diagnostics: {field}")
        flow_path = directory / (session.stem + ".flow.json")
        slabs_path = directory / (session.stem + ".slabs.csv")
        flow = json.loads(flow_path.read_text())
        with slabs_path.open() as stream:
            rows = [{key: float(value) for key, value in row.items()}
                    for row in csv.DictReader(stream)]
        if len(rows) != manifest["mesh"]["nx"] or any(
                not math.isfinite(value) for row in rows for value in row.values()):
            raise ValueError("missing slabs or nonfinite diagnostic output")
        receipt["profiles"].append({"index": index, "time_s": field_time,
                    "field_build_revision": build_revision, "checkpoint_inputs": before,
                    "command": command, "working_directory": str(directory),
                    "wall_seconds": time.monotonic() - start, "flow": flow,
                    "output_files": {str(p): digest(p) for p in (flow_path, slabs_path)},
                    "slabs": len(rows)})
        (output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
        print(f"checkpoint {index}: time={field_time:.9g}, slabs={len(rows)}", flush=True)
    if {str(p): digest(p) for p in source_files} != source_hashes:
        raise ValueError("case input changed during diagnostics")
    receipt["status"] = "completed"
    receipt["finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    (output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    main()
