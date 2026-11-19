#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bind the reviewed Task 2R inputs and freeze the initial 2R-S operation caps."""

import argparse
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np

from physical_oracle import NUMERICAL, STATIONS, load_native, load_reference, require
from resources import THREADS, atomic, identity


REVIEW = "366627a734b642788e18b9aebe0561074d72e92b"
SOURCE = "061c48fe6805b1b1bb66a8effffc3bfc28c5f86777c5e669f52afd9f8c0d861f"
FILES = {
    "application": ("runtime-r2-v1/IncNavierStokesSolver",
                    "7a94420aaa31f41d9d8b6754d20c16f11d1398e894f243eaa6beaa3af7924819"),
    "table": ("science/r1-qualified/tables-exact-inner-v1.json",
              "3eeb32e6c1f53d6045d9d03707344c65f0f1ec3964dc69116a87e04e104325e1"),
    "pre_radial_archive": ("evidence/task2r-r3-before-radial.tar.gz",
                           "f4b005465b2ec339fa8af810c6c4b438564cfefefd92608db4686023c1ddcffd"),
}
STATE_DIGESTS = {
    "r3-lowq": "2aef745e0f76d61b65682b290c68d4a30a5855bccb8d2cb20bec1f99ce2a08b1",
    "r3-outlet-reallocated": "8d286ddec3197cba12b1800bf76998bac45e27422da12f169484053434ef56ea",
    "r3-terminal-refined": "8f65f3b43e8281a1ed56d29f16997ae257478329da366cc496c748f138bac3d1",
}
CAPS = {
    "S0": dict(baseline_and_publication=75),
    "S1": dict(replay=260, boundary_and_admissibility=190, wedge_and_fit=100, publication=85),
    "S2_first": dict(inputs=45, steady=1080, current=265, terminal_compare=250, parent_compare=250,
                     boundary_fit_and_analysis=150, publication=100),
}


def verified(path, digest):
    item = identity(path)
    require(item["sha256"] == digest, f"reviewed digest mismatch for {path}")
    return item


def verify_state(directory, file_digest, ref):
    state_file = directory / "state.result.json"
    artifact = verified(state_file, file_digest)
    data, meta = load_native(state_file)
    require(meta["steady_mode"] and data["time"] == data["step"] == 0
            and data["source_identity"] == SOURCE and data["reference_identity"] == ref.identity,
            f"invalid reviewed steady state {directory.name}")
    signature = hashlib.sha256()
    for field in ("flow", "k", "epsilon"):
        signature.update(np.asarray(data["coefficients"][field], dtype="<f8").tobytes())
    require(signature.hexdigest() == data["state_sha256"], "state coefficient hash mismatch")
    geometry = data["geometry"]
    expected = dict(left=-10., right=15., d0=ref.d0, H=ref.H, top=ref.top, h=ref.h)
    require(all(abs(geometry[key]-value) <= 1e-12 for key, value in expected.items()),
            "retained geometry mismatch")
    session = directory / "run.xml"
    tree = ET.parse(session)
    vertices = np.asarray([list(map(float, v.text.split())) for v in tree.findall("GEOMETRY/VERTEX/V")])
    require(abs(vertices[:, 1].min() - ref.d0) <= 1e-12
            and abs(vertices[:, 1].max() - 20.) <= 1e-12
            and min(abs(vertices[:, 1] - ref.d0 - ref.H)) <= 1e-12, "XML matching height mismatch")
    elements = data["elements"]
    x = {round(e["bounds"][0], 12) for e in elements}
    r = {round(e["bounds"][2], 12) for e in elements}
    conditions = tree.find("CONDITIONS")
    info = {v.attrib["PROPERTY"]: v.attrib["VALUE"] for v in conditions.findall("SOLVERINFO/I")}
    parameters = dict(tuple(term.strip() for term in v.text.split("=", 1))
                      for v in conditions.findall("PARAMETERS/P"))
    modes = {v.attrib["FIELDS"]: v.attrib["NUMMODES"] for v in tree.findall("EXPANSIONS/E")}
    require(info["WJReferencePolicy"] == "ExactInnerFiniteNuV1"
            and info["WJNativeSolve"] == "SteadyReference" and info["WJNativeMode"] == "Physical"
            and modes == {"u,v,k,epsilon": "4,4", "p": "3,3"}, "reviewed fixture changed")
    require(len(elements) == len(x)*len(r)
            and sum(abs(e["bounds"][2]-ref.H) <= 1e-12 for e in elements) == len(x)
            and all(bool(e["epsilon"]) == (e["bounds"][2] >= ref.H-1e-12) for e in elements),
            "retained matching-plane ownership mismatch")
    return dict(artifact=artifact, session=identity(session),
                coefficient_sha256=data["state_sha256"], cells=[len(x), len(r)],
                polynomial_degrees=dict(velocity_scalars=3, pressure=2),
                quadrature=[int(parameters[k]) for k in ("WJQuadratureX", "WJQuadratureR")],
                actual_simulated_time=0., source=data["source_identity"], reference=data["reference_identity"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    root = args.root
    book = json.loads((root/"accounting.json").read_text())
    require(not book["closed"] and book["ceiling"] == 43200
            and book["protected_seconds"] == 3600, "wrong continuous H27 account")
    active = {name for name, row in book["entries"].items() if row["status"] == "RESERVED"}
    require(active <= {"s0-baseline"}, "prior reservation must be settled")
    spent = sum(row["seconds"] for row in book["entries"].values() if row["status"] != "RESERVED")
    charged = book["opening_seconds"] + spent
    require(charged >= 30739.447072760828, "historical charges were lost")
    manifest = root/"evidence/r2-source.sha256"
    native_files = {}
    for line in manifest.read_text().splitlines():
        digest, location = line.split(maxsplit=1)
        native_files[location] = verified(root/"source-r2-v1"/location, digest)["sha256"]
    compiled = root/"evidence/r2-compiled-source.txt"
    require(SOURCE in compiled.read_text(), "compiled source tag mismatch")
    artifacts = {key:verified(root/name, digest) for key,(name,digest) in FILES.items()}
    ref = load_reference(root/FILES["table"][0])
    states = {case:verify_state(root/"science"/case, digest, ref)
              for case,digest in STATE_DIGESTS.items()}
    atomic(args.out, dict(schema="WJH27Task2RSpatialDiagnosisContractV1", review_git=REVIEW,
        compiled_source=SOURCE, native_files=native_files, source_manifest=identity(manifest),
        compiled_source_record=identity(compiled), artifacts=artifacts, reviewed_states=states,
        geometry=dict(length_unit_m=.04, wall_elevation_above_original_floor_m=.04*ref.d0,
                      matching_height_m=.04*ref.H, matching_height_code=ref.H,
                      matching_plane_elevation_above_original_floor_m=.04*(ref.d0+ref.H),
                      top_elevation_above_original_floor_m=.8,
                      prior_original_floor_elevation_m_key="matching_plane_elevation_above_original_floor_m"),
        station_inventory=STATIONS, numerical_rms_and_max=NUMERICAL, operation_caps_seconds=CAPS,
        limits=dict(extra_reservation_per_operation_seconds=7, new_evidence_hold_seconds=300,
                    one_science_process=True, threads=THREADS, build_workers=4,
                    shared_worker_cap=72, science_rss_bytes=8*1024**3, build_rss_bytes=16*1024**3,
                    per_category_storage_bytes=8*1024**3, free_reserve_bytes=15*1024**3),
        ledger_before_s0_charge=dict(ceiling=book["ceiling"], carried_charge=charged,
            nonprotected_remaining=book["ceiling"]-book["protected_seconds"]-charged,
            staged_remaining=book["tranche_seconds"]-spent, protected=book["protected_seconds"],
            preexisting_reservations=[])))


if __name__ == "__main__":
    main()
