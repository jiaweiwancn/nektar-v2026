#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bind independent original-equation receipts to a fresh or cold hill run."""

import argparse
import json
from pathlib import Path
import xml.etree.ElementTree as ET

from resources import identity
from task2r_verify import edit, read


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("template", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--packet", type=Path)
    parser.add_argument("--equations", type=Path)
    parser.add_argument("--response", type=Path)
    parser.add_argument("--points", type=Path)
    parser.add_argument("--restart", type=Path)
    parser.add_argument("--steps", type=int, required=True)
    parser.add_argument("--recovery", type=int, default=1)
    parser.add_argument("--pulse-cell", type=int)
    parser.add_argument("--k-peak", type=float, default=0.)
    args = parser.parse_args()
    proofs = (args.packet, args.equations, args.response, args.points)
    if args.restart:
        if any(proofs) or args.k_peak or args.pulse_cell is not None:
            parser.error("cold restart has no coefficient packet, proof, or disturbance")
    elif not all(proofs):
        parser.error("fresh release requires the typed packet and the three independent audits")
    if bool(args.k_peak) != (args.pulse_cell is not None):
        parser.error("a signed k peak requires a native pulse cell")
    session = ET.parse(args.template)
    info = session.find("CONDITIONS/SOLVERINFO")
    parameters = session.find("CONDITIONS/PARAMETERS")
    for item in list(info):
        name = item.get("PROPERTY")
        if (name.startswith(("WJDiagnostic", "WJOriginalSteady", "WJGeometrySeed")) or
                name in ("WJInitialFields", "WJInitialSourceIdentity", "WJRestart")):
            info.remove(item)
    for item in list(parameters):
        name = item.text.split("=")[0].strip()
        if name.startswith(("WJDiagnostic", "WJOriginalSteady", "WJGeometrySeed")):
            parameters.remove(item)

    out = args.out.resolve()
    out.mkdir(parents=True)
    settings = dict(WJNativeSolve="Advance", WJOutputPrefix=out/"state",
                    WJCheckCurrentActions="False", WJDiagnosticCandidateOutput="True")
    values = dict(NumSteps=args.steps, IO_InfoSteps=1, IO_CheckSteps=args.recovery,
                  WJScientificSteps=0)
    if args.restart:
        settings["WJRestart"] = args.restart.resolve()
        files = [args.restart]
    else:
        packet = read(args.packet)["state"]
        receipt = out/"release-receipt.json"
        with receipt.open("x") as stream:
            json.dump(dict(schema="WJH27OriginalSteadyReleaseV1",
                           packet=identity(args.packet.resolve()),
                           coefficient_sha256=packet["state_sha256"],
                           equations=identity(args.equations.resolve()),
                           response=identity(args.response.resolve()),
                           points=identity(args.points.resolve())),
                      stream, indent=2, allow_nan=False)
            stream.write("\n")
        settings.update(WJInitialFields=args.packet.resolve(),
                        WJInitialSourceIdentity=packet["source_identity"],
                        WJOriginalSteadyReceipt=receipt,
                        WJOriginalSteadyReceiptSHA256=identity(receipt)["sha256"])
        if args.pulse_cell is not None:
            values.update(WJOriginalSteadyKPulseCell=args.pulse_cell, WJOriginalSteadyKPulse=args.k_peak)
        files = [*proofs, receipt]
    target = out/"run.xml"
    session.write(target, encoding="unicode")
    edit(target, settings, values)
    (out/"inputs.json").write_text(json.dumps([identity(path) for path in
                                               (args.template, target, *files)], indent=2)+"\n")


if __name__ == "__main__":
    main()
