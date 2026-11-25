#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Pin an immutable original-residual or diffusion-homotopy steady case."""

import argparse
import json
from pathlib import Path
import xml.etree.ElementTree as ET

from resources import identity
from task2r_verify import edit


def prepare(args):
    session = ET.parse(args.template)
    info = session.find("CONDITIONS/SOLVERINFO")
    parameters = session.find("CONDITIONS/PARAMETERS")
    for item in list(info):
        if item.get("PROPERTY") in {
            "WJInitialFields", "WJInitialSourceIdentity", "WJRestart", "WJGeometrySeedPacketSHA256",
            "WJDiagnosticSeedSHA256", "WJDiagnosticCandidateOutput", "WJTerminalOutput"
        }:
            info.remove(item)
    for item in list(parameters):
        if item.text.split("=", 1)[0].strip().startswith("WJGeometrySeed"):
            parameters.remove(item)
    table = Path(info.find("I[@PROPERTY='WJReferenceTable']").get("VALUE"))
    output = args.out.resolve()
    fields = dict(WJNativeSolve="DiagnosticSteady", WJOutputPrefix=output/"state",
                  WJDiagnosticIncompleteLinear="True" if args.incomplete_linear else "False",
                  WJSteadyVelocityPC=args.velocity_pc,
                  WJCheckCurrentActions="True" if args.check_actions else "False")
    controls = dict(NumSteps=0, WJSteadySeconds=args.seconds,
                    WJSteadyFDStep=args.fd, WJSteadyNewtonIterations=args.newton,
                    WJSteadyKrylovDimension=args.krylov, WJSteadyLinearCycles=args.cycles,
                    WJSteadyPreconditionerMass=args.inverse_mass, WJDiagnosticPseudoTime=args.pseudo_time,
                    WJDiagnosticCorrectionDiffusion=args.beta, WJDiagnosticEta=args.eta,
                    WJDiagnosticRegularizationViscosity=args.regularization_viscosity)
    if args.height is not None:
        if info.find("I[@PROPERTY='WJGeometry']") is None:
            raise ValueError("height override requires a normal chart")
        controls["WJHillHeight"] = args.height
    if args.quadrature:
        controls["WJQuadratureX"], controls["WJQuadratureR"] = args.quadrature
    if args.seed:
        packet = json.loads(args.seed.read_text())["state"]
        sha = identity(args.seed)["sha256"]
        fields.update(WJInitialFields=args.seed.resolve(),
                      WJInitialSourceIdentity=packet["source_identity"], WJDiagnosticSeedSHA256=sha)
        if args.remap:
            fields["WJGeometrySeedPacketSHA256"] = sha
            controls.update(WJGeometrySeedFromHeight=packet["geometry"]["hill_height"],
                            WJGeometrySeedTaper=args.remap_taper)
    elif args.remap:
        raise ValueError("geometry remap requires a seed packet")
    output.mkdir(parents=True)
    target = output/"run.xml"
    session.write(target, encoding="unicode")
    edit(target, fields, controls)
    files = [args.template, target, table]
    if args.seed:
        files.append(args.seed)
    (output/"inputs.json").write_text(json.dumps([identity(p) for p in files], indent=2)+"\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("template", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--seed", type=Path)
    parser.add_argument("--remap", action="store_true")
    parser.add_argument("--remap-taper", type=float, default=.30)
    parser.add_argument("--height", type=float)
    parser.add_argument("--quadrature", type=int, nargs=2)
    parser.add_argument("--seconds", type=float, required=True)
    parser.add_argument("--fd", type=float, default=4e-6)
    parser.add_argument("--newton", type=int, default=48)
    parser.add_argument("--krylov", type=int, default=40)
    parser.add_argument("--cycles", type=int, default=4)
    parser.add_argument("--inverse-mass", type=float, default=0.)
    parser.add_argument("--pseudo-time", type=float, default=0.)
    parser.add_argument("--beta", type=float, default=0.)
    parser.add_argument("--eta", type=float, default=0.)
    parser.add_argument("--regularization-viscosity", type=float, default=0.)
    parser.add_argument("--velocity-pc", choices=("FullGradientV1", "PositiveFrozenV1"),
                        default="FullGradientV1")
    parser.add_argument("--incomplete-linear", action="store_true")
    parser.add_argument("--check-actions", action="store_true")
    prepare(parser.parse_args())


if __name__ == "__main__":
    main()
