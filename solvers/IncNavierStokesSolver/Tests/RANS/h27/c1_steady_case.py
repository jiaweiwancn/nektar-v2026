#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Pin an immutable original-residual or diffusion-homotopy steady case."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import struct
import xml.etree.ElementTree as ET

from resources import identity
from task2r_verify import edit


def action_probe_difference(base, peer):
    """Map a signed, same-mesh packet difference to independent coefficients."""
    def load(path):
        packet = json.loads(path.read_text(), parse_int=lambda text: -0. if text == "-0" else int(text))
        state = packet["state"]
        if (packet.get("role") != "STEADY_SEARCH" or packet.get("iterate_custody") != "LAST_ACCEPTED"
                or packet.get("physical_history") or packet["regularization"]["eta"] != 0
                or packet["equation_identity"] != hashlib.sha256((state["signature"]+"|steady").encode()).hexdigest()):
            raise ValueError("direction peer requires bound accepted original diagnostic coefficients")
        coeff = state["coefficients"]
        raw = (v for block in ("flow", "k", "epsilon") for v in coeff[block])
        digest = hashlib.sha256(b"".join(struct.pack("<d", v) for v in raw)).hexdigest()
        if digest != state["state_sha256"]:
            raise ValueError("action probe coefficient digest differs")
        return state

    current, previous = load(base), load(peer)
    bindings = ("geometry", "elements", "sizes", "free_velocity_ids", "k_free_ids", "epsilon_free_ids",
                "gauge", "reference_identity", "outlet_identity", "quadratureX", "quadratureR")
    if any(current[key] != previous[key] for key in bindings):
        raise ValueError("action probe packets must use the identical signed free inventory and rule")
    nu, np, _, _ = current["sizes"]
    flow_ids = [c*nu+i for c, ids in enumerate(current["free_velocity_ids"]) for i in ids]
    flow_ids.extend(range(2*nu, 2*nu+np))
    direction = []
    for block, ids in (("flow", flow_ids), ("k", current["k_free_ids"]), ("epsilon", current["epsilon_free_ids"])):
        a, b = current["coefficients"][block], previous["coefficients"][block]
        direction.extend(a[i]-b[i] for i in ids)
        if block == "flow" and current["gauge"]:
            direction.append(0.)
    if not all(math.isfinite(value) for value in direction) or not any(direction):
        raise ValueError("action probe requires a finite nonzero difference")
    return struct.pack(f"<{len(direction)}d", *direction)


def prepare(args):
    session = ET.parse(args.template)
    info = session.find("CONDITIONS/SOLVERINFO")
    parameters = session.find("CONDITIONS/PARAMETERS")
    for item in list(info):
        if item.get("PROPERTY") in {
            "WJInitialFields", "WJInitialSourceIdentity", "WJRestart", "WJGeometrySeedPacketSHA256",
            "WJDiagnosticSeedSHA256", "WJDiagnosticCandidateOutput", "WJTerminalOutput",
            "WJDiagnosticTransferParent", "WJDiagnosticTransferParentSHA256",
            "WJDiagnosticActionProbe", "WJDiagnosticActionProbeSHA256", "WJDiagnosticActionGapAudit"
        }:
            info.remove(item)
    for item in list(parameters):
        if item.text.split("=", 1)[0].strip().startswith("WJGeometrySeed"):
            parameters.remove(item)
    table = Path(info.find("I[@PROPERTY='WJReferenceTable']").get("VALUE"))
    output = args.out.resolve()
    fields = dict(WJNativeSolve="DiagnosticSteady", WJOutputPrefix=output/"state",
                  WJDiagnosticIncompleteLinear="True" if args.incomplete_linear else "False",
                  WJDiagnosticSteadyInverse=args.inverse,
                  WJDiagnosticCorrectionAudit="True" if args.correction_audit else "False",
                  WJDiagnosticEvaluateOnly="True" if args.evaluate_only else "False",
                  WJSteadyVelocityPC=args.velocity_pc,
                  WJCheckCurrentActions="True" if args.check_actions else "False")
    if args.export_rejections:
        fields["WJDiagnosticCandidateOutput"] = "True"
    # Clear inherited numerical controls when reproducing the rebuild route.
    for name in ("WJDiagnosticDirectionModel", "WJDiagnosticGlobalization",
                 "WJDiagnosticUncommittedGuess", "WJDiagnosticAcceptedCheckpoint",
                 "WJDiagnosticMeritPolicy", "WJDiagnosticMeritEntry"):
        for item in info.findall(f"I[@PROPERTY='{name}']"):
            info.remove(item)
    for item in list(parameters):
        if item.text.split("=", 1)[0].strip() in (
                "WJDiagnosticInverseMaximumAge", "WJDiagnosticReuseLinearIterations",
                "WJDiagnosticCheckpointProcessSeconds", "WJDiagnosticCheckpointReserveSeconds",
                "WJDiagnosticCheckpointInverseSeconds", "WJDiagnosticMeritInitialMass",
                "WJDiagnosticMeritReanchorMass", "WJDiagnosticMeritWork",
                "WJSteadyForcingInitial", "WJSteadyForcingMinimum", "WJSteadyForcingMaximum",
                "WJSteadyForcingGamma", "WJSteadyForcingExponent", "WJDiagnosticFactorFDStep"):
            parameters.remove(item)
    controls = dict(NumSteps=0, WJSteadySeconds=args.seconds,
                    WJSteadyFDStep=args.fd, WJSteadyNewtonIterations=args.newton,
                    WJSteadyKrylovDimension=args.krylov, WJSteadyLinearCycles=args.cycles,
                    WJSteadyPreconditionerMass=args.inverse_mass, WJDiagnosticPseudoTime=args.pseudo_time,
                    WJDiagnosticCorrectionDiffusion=args.beta, WJDiagnosticEta=args.eta,
                    WJDiagnosticRegularizationViscosity=args.regularization_viscosity,
                    WJDiagnosticRegularizationStartRadius=args.regularization_radius[0],
                    WJDiagnosticRegularizationFullRadius=args.regularization_radius[1],
                    WJDiagnosticScalarSubdivisionDepth=args.scalar_subdivision)
    if args.maximum_inverse_age:
        controls.update(WJDiagnosticInverseMaximumAge=args.maximum_inverse_age,
                        WJDiagnosticReuseLinearIterations=args.reuse_linear_iterations)
    if args.accepted_checkpoint:
        outer, reserve, estimate = args.accepted_checkpoint
        if args.evaluate_only or args.eta != 0 or not outer > reserve + estimate or reserve <= 0 or estimate < 0:
            raise ValueError("accepted checkpoints need original correction and positive process and output budgets")
        fields["WJDiagnosticAcceptedCheckpoint"] = "True"
        controls.update(WJDiagnosticCheckpointProcessSeconds=outer, WJDiagnosticCheckpointReserveSeconds=reserve,
                        WJDiagnosticCheckpointInverseSeconds=estimate)
    if args.steady_forcing:
        initial, low, high = args.steady_forcing
        if not 0 < low <= initial <= high < 1:
            raise ValueError("steady forcing requires 0 < minimum <= initial <= maximum < 1")
        controls.update(WJSteadyForcingInitial=initial, WJSteadyForcingMinimum=low,
                        WJSteadyForcingMaximum=high)
    if args.merit_mass is not None:
        if (not all(math.isfinite(x) and x >= 0 for x in (args.merit_mass, args.merit_reanchor))
                or args.merit_work <= 0 or args.eta or args.beta or args.pseudo_time
                or args.evaluate_only or args.incomplete_linear or args.remap or args.uncommitted_guess):
            raise ValueError("original merit requires a fixed-height original diagnostic correction")
        if args.merit_entry != "Fresh" and not args.seed:
            raise ValueError("merit resume or re-anchor requires a pinned seed")
        fields.update(WJDiagnosticMeritPolicy="WJH27OriginalMeritV1",
                      WJDiagnosticMeritEntry=args.merit_entry)
        controls.update(WJDiagnosticMeritInitialMass=args.merit_mass,
                        WJDiagnosticMeritReanchorMass=args.merit_reanchor,
                        WJDiagnosticMeritWork=args.merit_work)
    elif args.merit_entry != "Fresh" or args.merit_reanchor or args.merit_work != 96:
        raise ValueError("merit settings require --merit-mass")
    if args.direction_model:
        fields["WJDiagnosticDirectionModel"] = "True"
    if args.globalization != "HalvingV1":
        fields["WJDiagnosticGlobalization"] = args.globalization
    if args.uncommitted_guess:
        if not args.seed or args.remap:
            raise ValueError("an uncommitted guess needs a digest-bound fixed-height seed")
        fields["WJDiagnosticUncommittedGuess"] = "True"
    if args.height is not None:
        if info.find("I[@PROPERTY='WJGeometry']") is None:
            raise ValueError("height override requires a normal chart")
        controls["WJHillHeight"] = args.height
    if args.quadrature:
        controls["WJQuadratureX"], controls["WJQuadratureR"] = args.quadrature
    if args.seed:
        bound = json.loads(args.seed.read_text())
        if bound.get("iterate_custody") == "UNCOMMITTED_TRIAL" and not args.uncommitted_guess:
            raise ValueError("uncommitted diagnostic seed requires --uncommitted-guess")
        if args.merit_mass is not None:
            controller = bound.get("merit_controller")
            if (controller and args.merit_entry=="Fresh") or (not controller and args.merit_entry=="Resume"):
                raise ValueError("merit continuation requires a saved controller and explicit Resume, or deliberate Reanchor")
        packet = bound["state"]
        sha = identity(args.seed)["sha256"]
        fields.update(WJInitialFields=args.seed.resolve(),
                      WJInitialSourceIdentity=packet["source_identity"], WJDiagnosticSeedSHA256=sha)
        if args.remap:
            fields["WJGeometrySeedPacketSHA256"] = sha
            controls.update(WJGeometrySeedFromHeight=packet["geometry"]["hill_height"],
                            WJGeometrySeedTaper=args.remap_taper)
    elif args.remap:
        raise ValueError("geometry remap requires a seed packet")
    if args.transfer_audit_parent:
        if not (args.seed and args.evaluate_only and args.check_actions) or (
                args.remap or args.uncommitted_guess or args.eta or args.beta or args.pseudo_time):
            raise ValueError("transfer audit requires a bound unshifted evaluate-only root with action checks")
        fields.update(WJDiagnosticTransferParent=args.transfer_audit_parent.resolve(),
                      WJDiagnosticTransferParentSHA256=identity(args.transfer_audit_parent)["sha256"])
    if args.action_gap_audit or args.action_probe_peer or args.factor_fd:
        if not args.seed or args.eta or args.beta or args.pseudo_time or args.remap or args.uncommitted_guess:
            raise ValueError("action audit settings require an exact original diagnostic base")
    if args.action_gap_audit:
        if args.evaluate_only or args.inverse != "ColoredSparseV1" or args.merit_mass != 0 or args.merit_reanchor:
            raise ValueError("action gap audit requires unshifted merit policy and ColoredSparseV1")
        fields["WJDiagnosticActionGapAudit"] = "True"
    if args.factor_fd is not None:
        if not math.isfinite(args.factor_fd) or args.factor_fd <= 0 or args.inverse != "ColoredSparseV1" or args.evaluate_only:
            raise ValueError("factor FD requires a finite positive step and a colored solve")
        controls["WJDiagnosticFactorFDStep"] = args.factor_fd
    direction = None
    if args.action_probe_peer:
        if not args.evaluate_only or args.action_gap_audit:
            raise ValueError("the action peer requires a separate evaluate-only probe")
        direction = action_probe_difference(args.seed, args.action_probe_peer)
        fields.update(WJDiagnosticActionProbe=output/"action-probe.bin",
                      WJDiagnosticActionProbeSHA256=hashlib.sha256(direction).hexdigest())
    output.mkdir(parents=True)
    if direction is not None:
        (output/"action-probe.bin").write_bytes(direction)
    target = output/"run.xml"
    session.write(target, encoding="unicode")
    edit(target, fields, controls)
    files = [args.template, target, table]
    if args.seed:
        files.append(args.seed)
    if args.transfer_audit_parent:
        files.append(args.transfer_audit_parent)
    if args.action_probe_peer:
        files.extend((args.action_probe_peer, output/"action-probe.bin"))
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
    parser.add_argument("--factor-fd", type=float, help="independent assembly spacing for the colored factor")
    parser.add_argument("--action-gap-audit", action="store_true",
                        help="capture one material failed Krylov check and audit the retained factor")
    parser.add_argument("--action-probe-peer", type=Path,
                        help="audit base minus peer in the identical free coordinates without correction")
    parser.add_argument("--newton", type=int, default=48)
    parser.add_argument("--krylov", type=int, default=40)
    parser.add_argument("--cycles", type=int, default=4)
    parser.add_argument("--inverse-mass", type=float, default=0.)
    parser.add_argument("--pseudo-time", type=float, default=0.)
    parser.add_argument("--beta", type=float, default=0.)
    parser.add_argument("--eta", type=float, default=0.)
    parser.add_argument("--regularization-viscosity", type=float, default=0.)
    parser.add_argument("--regularization-radius", type=float, nargs=2, default=(0.,0.))
    parser.add_argument("--scalar-subdivision", type=int, choices=range(7), default=0)
    parser.add_argument("--velocity-pc", choices=("FullGradientV1", "PositiveFrozenV1"),
                        default="FullGradientV1")
    parser.add_argument("--incomplete-linear", action="store_true")
    parser.add_argument("--check-actions", action="store_true")
    parser.add_argument("--inverse", choices=("BlockV1", "ColoredSparseV1"),
                        default="BlockV1")
    parser.add_argument("--correction-audit", action="store_true")
    parser.add_argument("--evaluate-only", action="store_true")
    parser.add_argument("--export-rejections", action="store_true")
    parser.add_argument("--maximum-inverse-age", type=int, default=0)
    parser.add_argument("--reuse-linear-iterations", type=int, default=12)
    parser.add_argument("--direction-model", action="store_true")
    parser.add_argument("--steady-forcing", type=float, nargs=3,
                        metavar=("INITIAL", "MINIMUM", "MAXIMUM"))
    parser.add_argument("--merit-mass", type=float,
                        help="enable original-merit recovery with this first correction mass")
    parser.add_argument("--merit-reanchor", type=float, default=0.,
                        help="explicit stronger same-base candidate before bounded changes")
    parser.add_argument("--merit-work", type=int, default=96,
                        help="maximum current Krylov iterations per nonlinear base")
    parser.add_argument("--merit-entry", choices=("Fresh", "Resume", "Reanchor"), default="Fresh")
    parser.add_argument("--globalization", choices=("HalvingV1", "ModelBacktrackV1"),
                        default="HalvingV1")
    parser.add_argument("--uncommitted-guess", action="store_true")
    parser.add_argument("--accepted-checkpoint", type=float, nargs=3,
                        metavar=("OUTER_SECONDS", "OUTPUT_RESERVE", "INVERSE_ESTIMATE"),
                        help="publish accepted corrections before later factor preparations")
    parser.add_argument("--transfer-audit-parent", type=Path,
                        help="same-height original steady parent for a nested-transfer residual/Jv audit")
    prepare(parser.parse_args())


if __name__ == "__main__":
    main()
