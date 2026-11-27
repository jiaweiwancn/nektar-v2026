#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Replay exact archived physical directions in the nonaccepted finite-trial diagnostic."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET

import numpy as np

import c1_steady_case
from c1_action_gap import read_audit
from resources import identity
from task2r_verify import edit

BASE = "ced15ecef5037f0e638664999a6ec82f46e02585a95895f907aa14ed3704c600"
TABLE = "3eeb32e6c1f53d6045d9d03707344c65f0f1ec3964dc69116a87e04e104325e1"
ARCHIVE = {
    "hybrid-fine/run.xml": "d75919798347090e68e573f59b5c5be99d1446bd392feefc2ca278db1a2506da",
    "hybrid-fine/state.accepted-000001.json": "84dae45ef585f24321137adb4b131ffc3157c15644ed1cd16ac7221b34caabb6",
    "final-rich/state.diagnostic-steady.json": "339a1a3390a12bb1a27d8b28e47b655623fce9ed7e9de7452f365afd747723bb",
    "final-rich/equations.json": "be02ed787ee012b079251d0eb6ee4d594ac0c46f78f10a92bafaad15adfdab66",
}
DIRECTIONS = (
    (.005, ".", 1, "5b669befe0f3f9190ea7b6ec8fac8f9c6dcc5b23620c29fdd980a6cec4e6110d"),
    (.025, ".", 2, "b7839f3d441482a4b43c3fdba8114945d384366537a0b6a500af97ad7cba74d7"),
    (.1, ".", 3, "dadc6f4c3bcdd99c38902cd83f2ae12221ee0de0686aa5dd34e66b059b69f359"),
    (.025, "ray", 1, "fa8bad394ad0e0288d4c1657d27050c8a0b4c662000091aaf9dc98e4b9928994"),
)
SETTINGS = ("geometry", "elements", "model", "basis", "boundary_identity", "outlet_identity",
            "reference_identity", "sizes", "gauge", "free_velocity_ids", "k_free_ids",
            "epsilon_free_ids", "manufactured", "mms_variant", "step_advection")
CASE_GENERATOR = Path(c1_steady_case.__file__).resolve()


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def load(path):
    return json.loads(path.read_text(), parse_int=lambda value: -0. if value == "-0" else int(value))


def save(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def packed(packet, expected_sha=BASE):
    state = packet["state"]
    blocks = [np.array(state["coefficients"][key], dtype="<f8") for key in ("flow", "k", "epsilon")]
    require(state["state_sha256"] == expected_sha == digest(b"".join(block.tobytes() for block in blocks)),
            "canonical coefficient identity differs")
    require(packet["equation_identity"] == digest((state["signature"] + "|steady").encode()) and
            packet["iterate_custody"] == "LAST_ACCEPTED" and not packet["physical_history"] and
            packet["regularization"]["eta"] == 0 and not state["gauge"],
            "expected the original, traction-outlet diagnostic checkpoint")
    ids = [np.asarray(i, dtype="<i4") for i in (*state["free_velocity_ids"],
                                                state["k_free_ids"], state["epsilon_free_ids"])]
    nv, npres, _, _ = state["sizes"]
    free = np.concatenate((blocks[0][ids[0]], blocks[0][nv+ids[1]],
                           blocks[0][2*nv:2*nv+npres], blocks[1][ids[2]], blocks[2][ids[3]]))
    return free, digest(b"".join(i.tobytes() for i in ids))


def physical_xml(path):
    root = ET.parse(path)
    settings = {}
    for block, entries in (("SOLVERINFO", lambda node: (node.attrib["PROPERTY"], node.attrib["VALUE"])),
                           ("PARAMETERS", lambda node: [s.strip() for s in node.text.split("=", 1)])):
        settings[block] = {key: value for key, value in map(entries, root.find("CONDITIONS/"+block))
                           if not key.startswith(("WJDiagnostic", "WJSteady", "WJInitialSource"))
                           and key not in {"WJOutputPrefix", "WJInitialFields", "WJCheckCurrentActions"}}
    settings["GEOMETRY"] = ET.tostring(root.find("GEOMETRY")).split()
    settings["EXPANSIONS"] = ET.tostring(root.find("EXPANSIONS")).split()
    return settings


def prepare(args):
    r9, frozen, case = (p.resolve() for p in (args.r9, args.frozen, args.case))
    checked = {}
    for name, sha in ARCHIVE.items():
        record = identity(r9/name)
        require(record["sha256"] == sha, f"archive bytes differ for {name}")
        checked[name] = record
    seed = r9/"hybrid-fine/state.accepted-000001.json"
    template = r9/"hybrid-fine/run.xml"
    packet, rich = load(seed), load(r9/"final-rich/state.diagnostic-steady.json")
    base, ordering = packed(packet)
    packed(rich)
    state = packet["state"]
    require(all(state[k] == rich["state"][k] for k in SETTINGS) and len(state["elements"]) == 5238 and
            base.size == 175027 and state["reference_identity"] == TABLE and
            (state["quadratureX"], state["quadratureR"]) == (10, 16) and
            (rich["state"]["quadratureX"], rich["state"]["quadratureR"]) == (16, 24),
            "archived native and rich signed equation settings differ")
    old_control = packet["merit_controller"]["state"]
    require(old_control["next_mass"] == 0 and old_control["next_forcing"] == .3 and
            old_control["step"] == 4.728531382401797e-8, "archived merit controller differs")
    report = load(frozen/"combined-results.json")
    require(report["targets"] == [1.9315264416093993e-6, 4.283180208592735e-5],
            "archived original forcing targets differ")
    selected, pieces, common = [], [], None
    for radius, directory, index, sha in DIRECTIONS:
        parent = frozen/directory
        receipt = parent/"verify"/f"state.action-frozen-size-candidate-base-0-proposal-{index}-sequence-{index+1}.json"
        record = identity(receipt)
        published = next(row for row in report["candidates"] if row["receipt"]["path"] == str(receipt))
        require(record == published["receipt"] or record["sha256"] == published["receipt"]["sha256"],
                f"published receipt digest differs for {receipt}")
        metadata, arrays = read_audit(receipt)
        source = load(parent/"candidates.json")
        spec = source["candidates"][index]
        matrix = Path(source["input"]["path"])
        require(identity(matrix)["sha256"] == source["input"]["sha256"] == metadata["probe_input_sha256"],
                "archived matrix identity differs")
        with matrix.open("rb") as stream:
            stream.seek(index*base.size*8)
            column = stream.read(base.size*8)
        require(digest(column) == sha == spec["direction_sha256"] and spec["radius"] == radius,
                "canonical direction column differs")
        verify = {"base":BASE, "ordering":ordering, "quadrature":[10,16], "table":TABLE,
                  "unknown_scale":digest(arrays["unknown_scale"][:,0].tobytes()),
                  "row_scale":digest(arrays["row_scale"][:,0].tobytes()),
                  "raw_base_residual":digest(arrays["raw_base_residual"][:,0].tobytes())}
        require(metadata["base_coefficient_sha256"] == BASE and metadata["free_id_sha256"] == ordering and
                metadata["quadrature"] == verify["quadrature"] and metadata["reference_identity"] == TABLE and
                metadata["correction_mass"] == 0 and metadata["fixed_k_tangent"] and
                metadata["current_fd"] == .000512 and
                np.array_equal(arrays["base"][:,0], base), "archived action field or equation settings differ")
        direction = np.frombuffer(column, dtype="<f8")
        np.testing.assert_allclose(arrays["direction"][:,0], direction, rtol=5.e-16, atol=0.)
        np.testing.assert_allclose(np.max(np.abs(direction/arrays["unknown_scale"][:,0])), radius,
                                   rtol=1.e-11, atol=0.)
        require(common is None or verify == common, "frozen direction receipts do not share one base and scales")
        common = verify
        pieces.append(column)
        selected.append({"case": "bounded" if directory == "." else "old_ray", "radius": radius,
                         "canonical_input_sha256": sha, "receipt": record,
                         "published_direction_sha256": digest(arrays["direction"][:,0].tobytes()),
                         "old_model_physical": published["residuals"][0][0]})
    subprocess.run([sys.executable, str(CASE_GENERATOR), str(template), str(case),
                    "--seed", str(seed), "--eta", "0", "--quadrature", "10", "16",
                    "--evaluate-only", "--seconds", str(args.seconds), "--fd", ".000512",
                    "--steady-forcing", ".3", ".001", ".3"], check=True)
    payload = b"".join(pieces)
    combined = case/"finite-directions.bin"
    combined.write_bytes(payload)
    token = digest("|".join(["WJH27FrozenFiniteBindingV1", BASE, ordering, common["unknown_scale"],
                             common["row_scale"], digest(payload)]).encode())
    session = case/"run.xml"
    edit(session, {"WJDiagnosticActionProbe":combined, "WJDiagnosticActionProbeSHA256":digest(payload),
                   "WJDiagnosticActionProbeBinding":token, "WJDiagnosticFixedKTangent":"True"},
                  {"WJDiagnosticActionProbeCount":len(pieces), "WJDiagnosticFrozenFiniteTrials":4})
    require(physical_xml(session) == physical_xml(template), "prepared finite session changed physical settings")
    inputs = load(case/"inputs.json")
    require(inputs[2]["sha256"] == TABLE, "reference table bytes differ")
    save(case/"frozen-trial-inputs.json", {"archive":checked, "fixed_base":common, "canonical":selected,
                                          "input":identity(combined), "session":identity(session), "binding":token,
                                          "table":inputs[2], "archive_signature":packet["equation_identity"]})


def summarize(case):
    case = case.resolve()
    manifest = load(case/"frozen-trial-inputs.json")
    previous = load(Path(manifest["archive"]["hybrid-fine/state.accepted-000001.json"]["path"]))
    final = load(case/"state.diagnostic-steady.json")
    packed(final)
    require(all(previous["state"][k] == final["state"][k] for k in SETTINGS) and
            final["original_residuals"] == previous["original_residuals"] and
            all(final[k] == previous[k] for k in ("minimum_covariance_over_k", "original_principal_over_mu",
                                                   "physical_whole_cell_scalar_lower_bounds")) and
            final["state"]["interface_K_over_Kin"] == previous["state"]["interface_K_over_Kin"],
            "current source failed frozen original base or scientific-diagnostic replay")
    log = load(case/"state.steady.json")
    require(log["accepted_newton_steps"] == 0 and log["physical_time_advanced"] == 0 and
            "merit_controller" not in final and not list(case.glob("state.accepted-*.json")),
            "audit-only trial advanced state or controller")
    canon = manifest["canonical"]
    base = manifest["fixed_base"]
    trials = []
    for receipt in sorted(case.glob("state.action-frozen-size-finite-trial-*.json")):
        metadata, arrays = read_audit(receipt)
        proposal = metadata["proposal"]
        digests = {s["name"]:s["sha256"] for s in metadata["arrays"]}
        require(metadata["base_coefficient_sha256"] == BASE and
                metadata["probe_input_sha256"] == manifest["input"]["sha256"] and
                metadata["free_id_sha256"] == base["ordering"] and
                all(digests[k] == base[k] for k in ("unknown_scale", "row_scale", "raw_base_residual")) and
                digests["canonical_input_direction"] == canon[proposal]["canonical_input_sha256"],
                "native finite receipt does not bind the canonical input and fixed base")
        np.testing.assert_allclose(arrays["direction"], arrays["canonical_input_direction"],
                                   rtol=5.e-16, atol=0.)
        meas = metadata["measurements"]
        displacement = (arrays["trial_free"]-arrays["base"])/arrays["unknown_scale"]
        require(np.array_equal(arrays["trial_free"], arrays["base"]+
                meas["fraction"]*arrays["canonical_input_direction"]) and
                np.max(np.abs(displacement)) == meas["scaled_size"], "actual trial displacement differs")
        if "raw_remainder" in arrays:
            require(np.array_equal(arrays["raw_remainder"], arrays["raw_actual_residual"]-
                    arrays["raw_model_residual"]), "signed nonlinear remainder differs")
        probes = {s["name"]:s for s in metadata["arrays"] if "physical_blocks" in s}
        trials.append({"proposal":proposal, "candidate":canon[proposal], "receipt":identity(receipt),
                       "payload":metadata["raw_payload"], "trial":metadata["iteration"],
                       "notes":metadata["notes"], "metrics":meas, "normed_arrays":probes,
                       "scaled_displacement_blocks": [float(np.max(np.abs(block))) if block.size else 0.
                           for block in np.split(displacement[:,0], np.cumsum(metadata["free_counts"])[:-1])]})
    trials.sort(key=lambda x:(x["proposal"], x["trial"]))
    require({entry["proposal"] for entry in trials} == set(range(len(canon))), "finite trial panel incomplete")
    save(case/"frozen-trial-results.json", {"base_packet":identity(case/"state.diagnostic-steady.json"),
        "base_residuals":final["original_residuals"], "equation_identity":final["equation_identity"],
        "source_identity":final["state"]["source_identity"], "physical_time_advanced":log["physical_time_advanced"],
        "accepted_updates":log["accepted_newton_steps"], "inputs":identity(case/"frozen-trial-inputs.json"),
        "trials":trials})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prep = commands.add_parser("prepare")
    for name in ("r9", "frozen", "case"):
        prep.add_argument(name, type=Path)
    prep.add_argument("--seconds", type=float, default=750.)
    summary = commands.add_parser("summarize")
    summary.add_argument("case", type=Path)
    args = parser.parse_args()
    if args.command == "prepare":
        prepare(args)
    else:
        summarize(args.case)


if __name__ == "__main__":
    main()
