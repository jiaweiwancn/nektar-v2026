#!/usr/bin/env python3
"""Verify the R3 bindings and prepare the bounded R4 experiment on Linux."""
import argparse
from decimal import Decimal, localcontext
import importlib.util
import json
from pathlib import Path
import platform
import struct
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("predictor", HERE.parent / "m8-r3/prepare-predictor.py")
predictor = importlib.util.module_from_spec(spec)
spec.loader.exec_module(predictor)
pin, require, load = predictor.identity, predictor.require, predictor.load
WITNESS = Decimal("1.6046496089512394")


def bound(item):
    path = Path(item["path"])
    require(pin(path) == item, "Changed retained payload " + str(path))
    return path


def vector(item):
    return list(struct.unpack("<2165d", bound(item).read_bytes()))


def save(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def wall_rows(raw, state, precision):
    """Compute extrema from the signed native coefficient map in Decimal."""
    with localcontext() as ctx:
        ctx.prec = precision
        offset = raw["offsets"][2]
        results = []
        for e in raw["fields"][2]["mesh"]:
            if e["lower"][2] != 0:
                continue
            require(e["degree"] == 3, "Expected the frozen cubic wall trace")
            lo, hi = map(Decimal, (e["lower"][0], e["upper"][0]))
            l, r, even, odd = (Decimal(state[offset+i])*sign for i, sign in zip(e["ids"][:4], e["sign"][:4]))
            p = ((l+r)/2+even/4, (r-l+odd)/2, -even/4, -odd/2)
            A, B, C = 3*p[3], 2*p[2], p[1]
            stationary = []
            if A and B*B-4*A*C >= 0:
                root = (B*B-4*A*C).sqrt()
                stationary = [(-B-root)/(2*A), (-B+root)/(2*A)]
            elif not A and B:
                stationary = [-C/B]
            points = [("left", Decimal(-1)), ("right", Decimal(1))]
            points += [("stationary", s) for s in sorted(set(stationary)) if -1 < s < 1]
            if lo <= WITNESS <= hi:
                points.append(("original_GLL", 2*(WITNESS-lo)/(hi-lo)-1))
            samples = []
            for name, z in points:
                value = p[0]+z*(p[1]+z*(p[2]+z*p[3]))
                samples.append(dict(kind=name, x=str(lo+(z+1)*(hi-lo)/2), value=str(value)))
            results.append(dict(geometry=e["geometry"], bounds=[str(lo), str(hi)],
                modal=[str(x) for x in (l, r, even, odd)], power=[str(x) for x in p],
                samples=samples, minimum=min(samples, key=lambda point: Decimal(point["value"]))))
        require(results, "Missing ground elements")
        return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("epoch", "r3-ledger", "output"):
        parser.add_argument("--"+name, type=Path, required=True)
    args = parser.parse_args()
    require(platform.system() == "Linux", "Use the designated Linux system")
    epoch, output = args.epoch.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    old = load(args.r3_ledger)
    raw = load(bound(old["target"]["native_report"]))
    collection = load(bound(old["target"]["durable_collection"]))
    previous = load(bound(old["inputs"]["predictor"]))
    history = vector(old["inputs"]["original_h0"])
    snapshots = collection["completed_iterates"]
    states = [vector(row["state"]) for row in snapshots]
    findings = []
    for record, state in zip(snapshots, states):
        bound(record["snapshot"])
        require(record["history"]["sha256"] == predictor.H0 and record["alpha"] == 400 and
                record["time"] == .0025 and not record["history_advanced"], "Changed saved target")
        rows = wall_rows(raw, state, 80)
        check = wall_rows(raw, state, 110)
        differences = [abs(Decimal(x["value"])-Decimal(y["value"]))
            for a, b in zip(rows, check) for x, y in zip(a["samples"], b["samples"], strict=True)]
        require(max(differences) < Decimal("1e-65"), "Wall extrema disagree on precision ladder")
        lowest = min((e["minimum"] for e in check), key=lambda p: Decimal(p["value"]))
        minimum = Decimal(lowest["value"])
        sign = "positive" if minimum > Decimal("1e-60") else "negative" if minimum < Decimal("-1e-60") else "unresolved"
        findings.append(dict(index=record["iteration"], state=record["state"], snapshot=record["snapshot"],
            historical_inventory_admissible=record["scalar_domain"]["admissible"],
            historical_native_inner_minimum=record["scalar_domain"]["minimum_inner_k"],
            wall_minimum=lowest, wall_sign=sign, wall=check, precision_digits=[80, 110]))
    eligible = [r for r in findings[:5] if r["wall_sign"] == "positive" and r["historical_inventory_admissible"]]
    require(eligible, "No sign-resolved saved wall trace in checked recovery sequence")
    selected = max(r["index"] for r in eligible)
    base = snapshots[selected]
    matrix = next(r["event"]["matrix"] for r in collection["jacobians"] if r["iteration"] == selected)
    bound(matrix)
    bound(base["residual"])
    bound(base["constraints"])
    runtime = epoch / "runtime-fingerprint.json"
    exe = pin(epoch / "build/IncNavierStokesSolver")
    current = load(runtime)
    require(current["success"] and current["runtime_artifacts"]["build/IncNavierStokesSolver"] == exe["sha256"],
            "Current runtime/executable binding differs")
    reference_xml = bound(old["inputs"]["target_xml"])
    preflight_xml = bound(old["inputs"]["current_preflight_xml"])
    preflight = bound(old["inputs"]["current_preflight_report"])
    for phase in ("local", "pilot"):
        case = output / phase
        case.mkdir()
        tree = ET.parse(reference_xml)
        info = tree.find(".//SOLVERINFO")
        selectors = {r.get("PROPERTY"): r for r in info}
        selectors["RANSWallStageTest"].set("VALUE", "DirectionExperimentR4" if phase == "local" else "ExactTargetR3")
        selectors["RANSWallR3PredictorFile"].set("VALUE", str(case / "predictor.json"))
        selectors["RANSWallR3ProgressFile"].set("VALUE", str(case / "case.r3-progress.jsonl"))
        ET.SubElement(info, "I", PROPERTY="RANSWallGlobalization", VALUE="CoupledR4")
        xml = case / "case.xml"
        tree.write(xml, encoding="unicode")
        digest = predictor.protected_xml(xml)
        require(digest == predictor.protected_xml(reference_xml) == predictor.protected_xml(preflight_xml),
                "Changed frozen XML settings")
        manifest = dict(previous, role="retained_r3_internal_iterate_only", state=base["state"],
            history=base["history"], current_executable=exe, current_target_xml=pin(xml),
            current_preflight=pin(preflight), current_preflight_session=pin(preflight_xml),
            current_runtime=pin(runtime), producing_source_commit=current["source_commit"],
            protected_physical_xml_sha256=digest, origin=dict(alpha=400, time=.0025, converged=False,
                inner_stress_model="SaturatedR2", iteration=selected, snapshot=base["snapshot"],
                residual=base["residual"], constraints=base["constraints"], matrix=matrix))
        save(case / "predictor.json", manifest)
    probe = ET.parse(reference_xml.with_name("case-m4-native-rank0.xml"))
    ET.SubElement(probe.find(".//SOLVERINFO"), "I", PROPERTY="RANSWallGlobalization", VALUE="CoupledR4")
    probe.write(output / "probe.xml", encoding="unicode")
    packet = dict(state=states[selected], history=history, alpha=400, time=.0025)
    save(output / "saved-states.json", dict(packet,
        audit_states=[dict(index=i, state=s) for i, s in enumerate(states)]))
    save(output / "direction-state.json", dict(packet, only_external_direction=True))
    save(output / "wall-check.json", dict(scope="All five ground cubic traces, including limiting element and its neighbor",
        selected=selected, states=findings, original_gll=str(WITNESS),
        report=old["target"]["native_report"], collection=old["target"]["durable_collection"],
        base=base["snapshot"], matrix=matrix, h0=old["inputs"]["original_h0"],
        script=pin(__file__), r3_terminal_norm=max(snapshots[-1]["scaled_block_norms"])))


if __name__ == "__main__":
    main()
