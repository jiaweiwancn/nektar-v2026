#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Exercise the native release with scientifically distinct, correctly bound receipts."""

import argparse
import copy
import json
from pathlib import Path
import sys

from c1_steady_continue import call, original_admission, read
from resources import identity


def write(path, report):
    with path.open("x") as f:
        json.dump(report, f, indent=2, allow_nan=False)
        f.write("\n")


def variant(name, response, points, out):
    response, points = copy.deepcopy(response), copy.deepcopy(points)
    if name == "missing_support_region":
        response["regions"].pop("forcing_support")
    elif name == "missing_support_minimum":
        points["local_minima"] = [r for r in points["local_minima"] if r["region"] != "forcing_support"]
    elif name == "missing_reported_support_point":
        support = next(r for r in points["local_minima"] if r["region"] == "forcing_support")
        support["values"] = [v for v in support["values"] if v["side"] != "reported-minimum"]
    elif name == "adverse_support_point":
        support = next(r for r in points["local_minima"] if r["region"] == "forcing_support")
        reported = next(v for v in support["values"] if v["side"] == "reported-minimum")
        reported["minima"][0]["over_mu"] = -.001
    elif name == "missing_fixed_owner":
        index = next(i for i, v in enumerate(points["fixed_controls"]) if v.get("control") == "support_original")
        points["fixed_controls"].pop(index)
    elif name == "adverse_fixed_owner":
        fixed = next(v for v in points["fixed_controls"] if v.get("control") == "support_original")
        fixed["minima"][0]["over_mu"] = -.001
    out.mkdir()
    response_file, points_file = out/"response.json", out/"points.json"
    write(response_file,response)
    for local in points["local_minima"]:
        local["receipt"] = identity(response_file)
    write(points_file,points)
    return response_file, points_file


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("template","packet","equations","response","points","table","solver","out"):
        p.add_argument("--"+name, type=lambda s:Path(s).resolve(), required=True)
    args = p.parse_args()
    args.out.mkdir()
    expected = {
        "valid": ("QUALIFIED", ""),
        "missing_support_region": ("AUDIT_IDENTITY_FAILURE", "missing independent original response for forcing_support"),
        "missing_support_minimum": ("AUDIT_IDENTITY_FAILURE", "missing local original response audit for forcing_support"),
        "missing_reported_support_point": ("AUDIT_IDENTITY_FAILURE", "missing reported one-sided forcing support witness"),
        "adverse_support_point": ("ADVERSE_CONVERGED_ROOT", "independent original point response failed for forcing_support"),
        "missing_fixed_owner": ("AUDIT_IDENTITY_FAILURE", "missing fixed original response owner for support_original"),
        "adverse_fixed_owner": ("ADVERSE_CONVERGED_ROOT", "independent original point response failed for support_original"),
    }
    result = {}
    for name, (admission, message) in expected.items():
        case = args.out/name
        response, points = args.response,args.points
        if name!="valid":
            response,points=variant(name,read(args.response),read(args.points),args.out/(name+"-proofs"))
        decision=original_admission(args.packet,args.table,
                                    dict(equations=args.equations,margins=response,points=points))
        if decision!=admission:
            raise ValueError(f"{name} returned {decision}, expected {admission}")
        prepare=[sys.executable,Path(__file__).resolve().with_name("c1_release_original.py"),
                 args.template,case,"--packet",args.packet,"--equations",args.equations,
                 "--response",response,"--points",points,"--steps",1]
        assert call(prepare,args.out/(name+".setup.log"),20.,args.out)==0
        native=[args.solver,"run.xml"]
        log=case/"native.log"
        status=call(native,log,100.,case)
        artifacts={part.name:identity(part) for part in case.glob("state*.bin")}
        if name=="valid":
            assert status==0 and "state.final.bin" in artifacts
            time=read(case/"state.result.json")["state"]["time"]
            assert time>0.
        else:
            assert status==1 and message in log.read_text() and not artifacts
            assert not (case/"state.result.json").exists() and not (case/"state.steps.jsonl").exists()
            time=None
        result[name]=dict(python_admission=decision,exit_status=status,native_scientific_message=message,
                          physical_time=time,history_files=artifacts,native_log=identity(log),
                          response=identity(response),points=identity(points),
                          executed_receipt=identity(case/"release-receipt.json"),
                          session=identity(case/"run.xml"))
    write(args.out/"receipt-regression.json",result)


if __name__=="__main__":
    main()
