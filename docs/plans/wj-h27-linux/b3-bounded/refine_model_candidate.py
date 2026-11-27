#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Recompute one capped model using a verified feasible smaller-radius incumbent."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

import c1_bounded_trials as frozen
import c1_bounded_optimizer as policy
from c1_action_gap import read_audit
from resources import identity
from task2r_verify import edit


def prepare(source, template, case):
    # The zero candidate at .025 and the accepted .0125 trial used this same base.
    seed = source/"state.accepted-000001.json"
    basis, = (p for p in source.glob("state.action-bounded-current-basis-base-1-*.json") if ".json." not in p.name)
    decline, = source.glob("state.action-bounded-direction-base-1-proposal-2-*.json")
    small, = source.glob("state.action-bounded-direction-base-1-proposal-3-*.json")
    archived, a = read_audit(basis)
    rejected, _ = read_audit(decline)
    incumbent, prev = read_audit(small)
    state = frozen.load(seed)
    coef = state["state"]["state_sha256"]
    base, order = frozen.packed(state, coef)
    assert all(r["base_coefficient_sha256"] == coef and r["free_id_sha256"] == order
               for r in (archived, rejected, incumbent))
    assert np.array_equal(base, a["base"][:,0])
    assert all(np.array_equal(a[k], prev[k]) for k in ("row_scale", "unknown_scale", "base", "raw_base_residual"))
    assert rejected["notes"]["decision"]=="NO_NONZERO_BOUNDED_PROPOSAL"
    assert rejected["measurements"]["radius"]==.025
    assert incumbent["notes"]["decision"]=="BOUNDED_ACCEPTED" and incumbent["measurements"]["radius"]==.0125
    assert archived["equation_identity"]==state["equation_identity"]
    assert policy.sha(policy.HELPERS/"frozen_size.py")==policy.HELPER_SHA
    sys.path.insert(0,str(policy.HELPERS))
    from frozen_size import BLOCKS, JointModel, array_spec
    q, f = [a[k] for k in ("scaled_correction_basis","raw_residual_and_actions")]
    sx, sr = a["unknown_scale"][:,0], a["row_scale"][:,0]
    assert np.linalg.norm(q.T@q-np.eye(q.shape[1]))<1.e-10
    base_spec = array_spec(archived,"raw_base_residual")
    denominators = (max(base_spec["physical_blocks"]),base_spec["scaled_l2"])
    normalizer = a["response_column_normalization"][:,0]
    parts = []
    for i in range(len(BLOCKS)):
        gram=a[f"normalized_physical_gram_{i}"]
        values, vectors = np.linalg.eigh(gram)
        assert values[0] >= -1.e-10*max(values[-1],1.e-300)
        parts.append(((np.sqrt(np.maximum(0.,values))[:,None]*vectors.T)*normalizer[None,:])/denominators[0])
    parts.append(np.linalg.qr(sr[:,None]*f/denominators[1],mode="r"))
    model = JointModel(q,parts)
    old_y=prev["optimizer_coordinates"][:,0]
    np.testing.assert_allclose(q@old_y,prev["direction"][:,0]/sx,rtol=3.e-10,atol=2.e-15)
    old_v=np.linalg.solve(model.transform,old_y)
    _, y, info = model.solve(.025,previous=old_v,max_iterations=250)
    physical = sx*(q@y)
    assert 0.<np.max(np.abs(physical/sx))<=.025*(1.+1.e-10)
    opts = [sys.executable,str(frozen.CASE_GENERATOR),str(template),str(case),
            "--seed",str(seed),"--quadrature","10","16","--eta","0","--seconds","340",
            "--evaluate-only","--fd",".000512","--steady-forcing",".3",".001",".3"]
    subprocess.run(opts, check=True)
    directions=(physical.astype("<f8",copy=False),prev["direction"][:,0])
    matrix=case/"trial-directions.bin"
    matrix.write_bytes(b"".join(x.tobytes() for x in directions))
    payload=identity(matrix)
    token=frozen.digest("|".join(["WJH27FrozenFiniteBindingV1",coef,order,
                frozen.digest(sx.tobytes()),frozen.digest(sr.tobytes()),payload["sha256"]]).encode())
    edit(case/"run.xml",{"WJDiagnosticActionProbe":matrix,"WJDiagnosticActionProbeSHA256":payload["sha256"],
                         "WJDiagnosticActionProbeBinding":token,"WJDiagnosticFixedKTangent":"True"},
                        {"WJDiagnosticActionProbeCount":len(directions),"WJDiagnosticFrozenFiniteTrials":1})
    assert frozen.physical_xml(case/"run.xml")==frozen.physical_xml(template)
    frozen.save(case/"refinement-inputs.json",{"current_base":identity(seed),"coefficients":coef,
        "basis":identity(basis),"old_zero_proposal":identity(decline),"feasible_smaller_incumbent":identity(small),
        "model_fit":info,"directions":[frozen.digest(x.tobytes()) for x in directions],"matrix":payload,
        "bound_session":identity(case/"run.xml")})


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("source","template","case"):
        p.add_argument(key,type=Path)
    args=p.parse_args()
    prepare(args.source.resolve(),args.template.resolve(),args.case.resolve())
