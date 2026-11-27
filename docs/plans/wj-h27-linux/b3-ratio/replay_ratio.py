#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Re-evaluate signed terminal or reduced-model directions without a factor."""
import argparse
from pathlib import Path
import subprocess
import sys

import numpy as np

from c1_action_gap import read_audit
import c1_bounded_trials as frozen
from c1_recycled_basis import prepare as prepare_frame
from resources import identity
from task2r_verify import edit


SEED = "bounded-active/state.accepted-000003.json"
SESSION = "bounded-active/run.xml"
FRAME = "bounded-active/state.action-bounded-current-basis-base-3-proposal--1-sequence-9.json"
COEFFICIENTS = "bc593d6c54e937bc5ee0986fa9df67d40b39da34cd34c9f1037981406ec4f67b"
PINS = {
    SEED: "6dce7818c2a2db127fe7c0bcf145f597d918e9d7c4ea009d3672a387a1ebbd12",
    SESSION: "8d6bc1bb15a0ce1f0219129237b9dab648da1c163fcdc3bc7e5dc97053a5d801",
    FRAME: "cb81b7ea0d2f93f3841337e6c287e6a1fff40d512451d0e50186e07fb3fc7e5f",
    "endpoint-rich/state.diagnostic-steady.json":
        "9dd6d014900a4fed6451486f3510502bda25dfc9a6828fdc1d496f777c618045",
    "endpoint-rich/equations.json": "e021f18aaa5d291995725f03b0a221d70bcea665a838d9e0fe76d0d68cc478ad",
}
FRAME_PAYLOAD = "4098b36202f86f9a9e38ce455770fec13a6f4c09a71d4c137c95e0345d732045"


def prepare(old, manifest, case, proposals, frame_output, seconds):
    evidence = frozen.load(manifest)
    pins = {name:identity(old/name) for name in PINS}
    assert all(pins[name]["sha256"] == sha for name, sha in PINS.items())
    packet = frozen.load(old/SEED)
    rich = frozen.load(old/"endpoint-rich/state.diagnostic-steady.json")
    base, ordering = frozen.packed(packet, COEFFICIENTS)
    frozen.packed(rich, COEFFICIENTS)
    assert all(packet["state"][key] == rich["state"][key] for key in frozen.SETTINGS)
    state = packet["state"]
    assert len(state["elements"]) == 5238 and base.size == 175027
    assert state["reference_identity"] == frozen.TABLE
    frame, arrays = read_audit(old/FRAME)
    assert frame["raw_payload"]["sha256"] == FRAME_PAYLOAD and frame["factor_serial"] == 0
    assert frame["base_coefficient_sha256"] == COEFFICIENTS and frame["free_id_sha256"] == ordering
    assert frame["equation_identity"] == packet["equation_identity"]
    assert np.array_equal(arrays["base"][:,0], base)
    q = arrays["scaled_correction_basis"]
    assert q.shape == (base.size,18) and np.linalg.norm(q.T@q - np.eye(18)) < 1.e-10
    fixed = {key:frozen.digest(arrays[key].tobytes()) for key in (
        "unknown_scale","row_scale","raw_base_residual")}
    fixed["ordering"] = ordering
    if frame_output is not None:
        prepare_frame(old/SEED, old/FRAME, frame_output, retained_only=True)
    columns, recorded = [], []
    if not proposals:
        old_trials = evidence["lineages"]["bounded-active"]["run"]["trials"][-4:]
        for index, old_trial in enumerate(old_trials):
            name = f"bounded-active/state.action-bounded-direction-base-3-proposal-{index}-sequence-{10+index}.json"
            receipt = old/name
            assert identity(receipt) == old_trial["receipt"]
            trial, data = read_audit(receipt)
            assert trial["raw_payload"] == old_trial["payload"]
            assert all(trial[key] == frame[key] for key in (
                "base_coefficient_sha256", "free_id_sha256", "equation_identity", "reference_identity"))
            assert all(np.array_equal(data[key], arrays[key]) for key in fixed if key != "ordering")
            y = data["optimizer_coordinates"][:,0]
            d = data["direction"][:,0]
            np.testing.assert_allclose(q@y,d/arrays["unknown_scale"][:,0],rtol=1.e-10,atol=5.e-15)
            columns.append(d)
            recorded.append({"receipt":identity(receipt),"native_old":old_trial["measurements"],
                             "historical_direction_sha256":trial["arrays"][4]["sha256"]})
    for proposal in proposals:
        info = frozen.load(Path(str(proposal)+".json"))
        with proposal.open() as stream:
            header = stream.readline().split()
            y = np.asarray([float(x) for x in stream.read().split()])
        assert header[0] in ("WJH27BoundedProposalV1","WJH27BoundedProposalV2")
        assert header[1] == pins[FRAME]["sha256"] and int(header[3]) == len(y) == 18
        assert float(header[2]) == info["radius"] and np.max(np.abs(q@y)) <= info["radius"]*(1+1.e-10)
        columns.append(arrays["unknown_scale"][:,0]*(q@y))
        recorded.append({"optimizer_result":identity(proposal),"optimizer_receipt":identity(Path(str(proposal)+".json")),
                         "optimizer":info})
    subprocess.run([sys.executable,str(frozen.CASE_GENERATOR),str(old/SESSION),str(case),
                    "--seed",str(old/SEED),"--eta","0","--quadrature","10","16","--evaluate-only",
                    "--seconds",str(seconds),"--fd",".000512","--steady-forcing",".3",".001",".3"],
                   check=True)
    payload = b"".join(v.astype("<f8",copy=False).tobytes() for v in columns)
    directions = case/"finite-directions.bin"
    directions.write_bytes(payload)
    for v, record in zip(columns,recorded):
        record["canonical_input_sha256"] = frozen.digest(v.astype("<f8",copy=False).tobytes())
        assert record.get("historical_direction_sha256",record["canonical_input_sha256"]) == record["canonical_input_sha256"]
    token = frozen.digest("|".join(["WJH27FrozenFiniteBindingV1",COEFFICIENTS,ordering,
        fixed["unknown_scale"],fixed["row_scale"],frozen.digest(payload)]).encode())
    session = case/"run.xml"
    edit(session,{"WJDiagnosticActionProbe":directions,"WJDiagnosticActionProbeSHA256":frozen.digest(payload),
                  "WJDiagnosticActionProbeBinding":token,"WJDiagnosticFixedKTangent":"True"},
                 {"WJDiagnosticActionProbeCount":len(columns),"WJDiagnosticFrozenFiniteTrials":1})
    assert frozen.physical_xml(session) == frozen.physical_xml(old/SESSION)
    table = frozen.load(case/"inputs.json")[2]
    assert table["sha256"] == frozen.TABLE
    frozen.save(case/"frozen-trial-inputs.json",{"archive":pins,"fixed_base":fixed,"canonical":recorded,
        "frame_payload":frame["raw_payload"],"input":identity(directions),"session":identity(session),
        "binding":token,"table":table,"archive_signature":packet["equation_identity"]})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case",type=Path)
    parser.add_argument("--old",type=Path)
    parser.add_argument("--manifest",type=Path)
    parser.add_argument("--frame-output",type=Path)
    parser.add_argument("--proposal",action="append",type=Path,default=[])
    parser.add_argument("--seconds",type=float,default=480.)
    args = parser.parse_args()
    if args.old is not None:
        prepare(args.old,args.manifest,args.case,args.proposal,args.frame_output,args.seconds)
    else:
        frozen.summarize(args.case,SEED,COEFFICIENTS)
