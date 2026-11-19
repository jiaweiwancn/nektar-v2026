#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Seal compact 2R-S evidence and a complete prospective radial-run admission."""

import argparse
import json
from pathlib import Path
import tarfile

from resources import atomic, identity
from task2s_contract import CAPS


def read(path):
    return json.loads(path.read_text())


def first_radial_forecast(root, out, book):
    floor = 300
    caps = CAPS["S2_first"]
    total_reservation = sum(seconds+7 for seconds in caps.values())
    settled = sum(v["seconds"] for v in book["entries"].values() if v["status"] != "RESERVED")
    account = dict(charged=book["opening_seconds"]+settled,
                   staged=book["tranche_seconds"]-settled,
                   nonprotected=book["ceiling"]-book["protected_seconds"]-book["opening_seconds"]-settled)
    records = root/"science/task2rs-s1"
    replay = read(records/"observation-replay.json")
    wedge = read(records/"wedge-fit.json")
    assert replay["archive_replay"]["status"] == "PASS" and len(replay["rows"]) == 690 and len(replay["failures"]) == 28
    assert wedge["fitted_shape_pencil"]["status"] == "PASS"
    witnesses = [
        root/"operations/r3-terminal-steady/receipt.json",
        root/"operations/r3-terminal-current/receipt.json",
        root/"operations/r3-outlet-refine-steady/receipt.json",
        root/"operations/r3-outlet-lean-steady/receipt.json",
        root/"operations/s1-spatial-replay-corrected/receipt.json",
        root/"operations/s1-boundary-and-admissibility/receipt.json",
        root/"operations/s1-wedge-and-fit-serialization-retry/receipt.json",
    ]
    want = total_reservation+floor
    atomic(out,dict(schema="WJH27Task2RSRadialAdmissionV1",
        live_before_publication_charge=account, capped_operations_seconds=caps,
        total_reservation_seconds=total_reservation, evidence_hold_seconds=floor,
        staged_requirement_including_evidence_hold=want,
        minimum_new_stage_before_publication_charge=max(0.,want-account["staged"]),
        memory=dict(measured_terminal_peak=8005373952,science_rss_cap=8*1024**3,
            estimated_peak_bytes=round(8005373952*62/60+.25*1024**3),
            basis_dimension_previous=40,basis_dimension_next=30,
            largest_unknown="Sparse factor fill and peak coexistence on the changed wall mesh."),
        case=dict(cells=[114,62],polynomial_degrees=[3,2],native_quadrature=[6,8],
                  terminal_x_width=.00244140625,first_radial_width=.000125,
                  internal_solver_cap=980,original_exterior_contract=True,time_code=0.),
        evidence=[identity(p) for p in [*witnesses,records/"observation-replay.json",
                  records/"boundary-audit.json",records/"wedge-fit.json"]],
        interpretation=dict(contracting="Measure all rows against both neighbors and test a bounded sharper x/r level.",
           rss_stop="Use measured quiet strips to fit the local splits, prechecking inlet representation.",
           persistent_strong_divergence="Check near-wall and collar gradients, and require common resolved amplitude before claiming a singular obstruction.",
           worker_preflight="Defer and retry on the same ledger when worker headroom returns.")))


def package(root, stage):
    book = read(root/"accounting.json")
    active = {k for k,v in book["entries"].items() if v["status"] == "RESERVED"}
    if any(not name.startswith(stage+"-publish") for name in active):
        raise ValueError(f"still reserved at publication: {sorted(active)}")
    evidence = root/"evidence"
    if stage == "s1":
        first_radial_forecast(root,evidence/"task2rs-s2-first-forecast.json",book)
    paths = {Path(__file__).resolve(),evidence/"s0-baseline.json"}
    if stage == "s1":
        paths.add(evidence/"task2rs-s2-first-forecast.json")
    paths.update((root/"tools").glob("s-support*/h27/task2s*.py"))
    paths.update(p for p in evidence.glob("task2rs-*-forecast.json"))
    for label,event in book["entries"].items():
        if label.startswith((stage+"-","s0-")) and event["status"]!="RESERVED":
            directory = root/"operations"/label
            paths.update(directory/name for name in ("request.json","receipt.json","command.log"))
    for directory in (root/"science").glob("task2rs-"+stage+"*"):
        for p in directory.iterdir():
            if p.is_file() and p.suffix in (".json",".xml",".opt") and p.name not in (
                    "state.result.json","state.steady-candidate.json","state.failure.json"):
                paths.add(p)
    paths = sorted(paths)
    manifest = evidence/f"task2rs-{stage}-manifest.json"
    atomic(manifest,dict(schema="WJH27Task2RSArchiveV1",stage=stage,
          files={str(p.relative_to(root)):identity(p) for p in paths},
          field_storage="The native states and raw histories stay in the named ext4 case directories."))
    archive = evidence/f"task2rs-{stage}.tar.gz"
    with tarfile.open(archive,"x:gz",compresslevel=5) as stream:
        for p in [*paths,manifest]:
            stream.add(p,arcname=str(p.relative_to(root)),recursive=False)
    atomic(evidence/f"task2rs-{stage}-bundle.json",dict(archive=identity(archive),manifest=identity(manifest)))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("root",type=Path)
    p.add_argument("stage",choices=("s1","s2","s3"))
    args=p.parse_args()
    package(args.root,args.stage)


if __name__=="__main__":
    main()
