#!/usr/bin/env python3
"""Check Milestone 1 examples and rejection cases, without running Nektar++.

This fixture checker handles literal numbers and named literal parameters in
the supplied examples. The future C++ parser uses the full session interpreter.
"""

import argparse
import copy
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import prepare_cbc1971 as cbc


HERE = Path(__file__).resolve().parent
U64_MAX = (1 << 64) - 1
INT_MAX = (1 << 31) - 1


def require(condition, message):
    if not condition:
        raise ValueError(message)


def attrs(element, required, optional=()):
    require(element is not None, "Missing element")
    names = set(element.attrib)
    require(set(required) <= names <= set(required) | set(optional),
            "Missing or unknown attributes on " + element.tag)


def decimal_integer(text, maximum=U64_MAX):
    require(re.fullmatch(r"[0-9]+", text) is not None, "Expected unsigned decimal text")
    value = int(text)
    require(value <= maximum, "Unsigned integer overflow")
    return value


def validate(root):
    require(root.tag == "NEKTAR", "Expected NEKTAR root")
    conditions = root.findall("CONDITIONS")
    require(len(conditions) == 1, "Expected one CONDITIONS block")
    condition = conditions[0]
    parameters = {}
    for item in condition.findall("PARAMETERS/P"):
        name, text = item.text.split("=", 1)
        name = name.strip()
        require(name not in parameters, "Duplicate fixture parameter")
        parameters[name] = text.strip()

    def number(text):
        result = Decimal(parameters.get(text, text))
        require(result.is_finite() and math.isfinite(float(result)), "Nonfinite number")
        return result

    props = {}
    for item in condition.findall("SOLVERINFO/I"):
        name = item.get("PROPERTY")
        require(name not in props, "Duplicate solver property")
        props[name] = item.get("VALUE")
    expected = {"SolverType": "VelocityCorrectionScheme", "EQTYPE": "UnsteadyNavierStokes",
                "EvolutionOperator": "Nonlinear", "Projection": "Continuous",
                "Extrapolation": "Standard"}
    require(all(props.get(key) == value for key, value in expected.items()),
            "Unsupported solver combination")
    require(props.get("HOMOGENEOUS") in (None, "1D"), "Unsupported expansion")
    require(props.get("ModeType", "MultipleModes") == "MultipleModes", "Unsupported modes")
    schemes = condition.findall("TIMEINTEGRATIONSCHEME")
    require(len(schemes) == 1, "Expected one integration block")
    scheme = schemes[0]
    require(scheme.findtext("METHOD", "").strip() == "IMEX",
            "Expected IMEX")
    require(scheme.findtext("VARIANT", "").strip() == "" and
            scheme.findtext("ORDER", "").strip() == "2", "Unsupported IMEX variant/order")
    require([item.text.strip() for item in condition.findall("VARIABLES/V")] ==
            ["u", "v", "w", "p"], "Expected u,v,w,p")
    require(number(parameters.get("Flowrate", "0")) == 0, "Flowrate must be disabled")
    for key in ("CFL", "CFLEnd", "IO_CheckTime", "FinTime"):
        require(number(parameters.get(key, "0")) == 0, "Unsupported time-step control " + key)
    for key in ("CFLGrowth", "TimeIncrementFactor"):
        require(number(parameters.get(key, "1")) == 1, "Unsupported time-step growth " + key)
    timestep = number(parameters["TimeStep"])
    require(timestep > 0, "Positive timestep required")
    steps = number(parameters["NumSteps"])
    require(steps == int(steps) and steps > 0, "Positive integer NumSteps required")

    elements = root.findall("HITINFLOW")
    require(len(elements) == 1, "Expected one HITINFLOW block")
    hit = elements[0]
    attrs(hit, ("VERSION", "UNITS"))
    require(hit.attrib == {"VERSION": "1", "UNITS": "SI"}, "Unsupported version or units")
    require([item.tag for item in hit] == ["REGION"], "Expected one region")
    region = hit[0]
    attrs(region, ("REF",))
    ref = decimal_integer(region.get("REF"), INT_MAX)
    require(len(list(region)) == 4 and {item.tag for item in region} ==
            {"ORIGIN", "SYNTHESIS", "SPECTRUM", "CUTOFF"}, "Unexpected region elements")
    global_refs = [decimal_integer(item.get("ID"), INT_MAX) for item in condition.findall("BOUNDARYREGIONS/B")]
    require(ref in global_refs and len(global_refs) == len(set(global_refs)), "Unknown global REF")
    all_boundaries = condition.findall("BOUNDARYCONDITIONS/REGION")
    hit_boundaries = [item for item in all_boundaries
                      if any(bc.get("USERDEFINEDTYPE") == "HITInflow" for bc in item)]
    require(len(hit_boundaries) == 1 and
            decimal_integer(hit_boundaries[0].get("REF"), INT_MAX) == ref,
            "Only the configured single HIT inlet is supported")
    boundaries = [item for item in all_boundaries
                  if decimal_integer(item.get("REF"), INT_MAX) == ref]
    require(len(boundaries) == 1, "Missing or duplicate inlet boundary")
    bcs = list(boundaries[0])
    require(len(bcs) == 4 and [item.get("VAR") for item in bcs] == ["u", "v", "w", "p"],
            "Expected four inlet fields")
    require(all(item.tag == "D" and item.get("USERDEFINEDTYPE") == "HITInflow"
                for item in bcs[:3]), "Expected vector HIT Dirichlet boundary")
    require(bcs[3].tag == "N" and bcs[3].get("USERDEFINEDTYPE") == "H" and
            number(bcs[3].get("VALUE")) == 0, "Expected H pressure boundary")
    speed = number(bcs[0].get("VALUE"))
    require(speed > 0 and all(number(item.get("VALUE")) == 0 for item in bcs[1:3]),
            "Unsupported mean velocity")

    origin = region.find("ORIGIN")
    attrs(origin, ("X", "Y", "Z"))
    for value in origin.attrib.values():
        number(value)
    synthesis = region.find("SYNTHESIS")
    attrs(synthesis, ("LX", "LY", "LZ", "NX", "NY", "NZ", "SEED",
                      "TIMEORIGIN", "RECORDDURATION", "REPEAT"))
    periods = [number(synthesis.get(key)) for key in ("LX", "LY", "LZ")]
    require(all(value > 0 for value in periods), "Box lengths must be positive")
    counts = [number(synthesis.get(key)) for key in ("NX", "NY", "NZ")]
    require(all(value == int(value) and 4 <= value <= (1 << 31) - 1 and int(value) % 2 == 0
                for value in counts), "Grid counts must be even integers at least four")
    count_product = math.prod(int(value) for value in counts)
    require(6 * (count_product - 1) + 5 <= U64_MAX, "Mode counter overflow")
    decimal_integer(synthesis.get("SEED"))
    number(synthesis.get("TIMEORIGIN"))
    duration = number(synthesis.get("RECORDDURATION"))
    require(duration > 0 and duration >= number(parameters["NumSteps"]) * timestep,
            "Record duration does not cover example run")
    repeat = synthesis.get("REPEAT")
    require(repeat in ("Error", "Allow"), "Unknown repeat policy")
    require(repeat == "Allow" or duration < periods[0] / speed, "Record repeats synthetic box")
    if props.get("HOMOGENEOUS") == "1D":
        require(periods[2] == number(parameters["LZ"]) and
                counts[2] == number(parameters["HomModesZ"]), "Inconsistent homogeneous grid")

    cutoff = region.find("CUTOFF")
    attrs(cutoff, ("KMAX", "TAPERSTART"))
    kmax, start = number(cutoff.get("KMAX")), number(cutoff.get("TAPERSTART"))
    require(0 <= start < kmax, "Invalid cutoff/taper")
    require(float(kmax) < min(math.pi * float(n / length) for n, length in zip(counts, periods)),
            "Cutoff reaches an excluded Nyquist plane")
    require(float(speed * kmax * timestep) < math.pi, "Temporal alias bound exceeded")

    spectrum = region.find("SPECTRUM")
    kind, normalization = spectrum.get("TYPE"), spectrum.get("NORMALIZATION")
    require(normalization in ("Preserve", "TargetRMS"), "Unknown normalization")
    optional = ("RMS",) if normalization == "TargetRMS" else ()
    if kind == "Tabulated":
        attrs(spectrum, ("TYPE", "FILE", "INTERPOLATION", "LOW", "HIGH", "NORMALIZATION"), optional)
        require(spectrum.get("INTERPOLATION") in ("Linear", "LogLinear"), "Unknown interpolation")
        require(spectrum.get("LOW") in ("Zero", "K4") and
                spectrum.get("HIGH") in ("Zero", "Gaussian"), "Unknown tail policy")
        filename = spectrum.get("FILE")
        require(Path(filename).name == filename, "Example expects a local fixture filename")
        text = (HERE / filename).read_text()
        require(text.splitlines()[0] ==
                "# NEKTAR_HIT_SPECTRUM_V1 KIND=Shell3D K_UNITS=rad/m E_UNITS=m3/s2",
                "Expected SI shell spectrum, not a longitudinal PSD")
        rows = [[Decimal(token) for token in line.split()]
                for line in text.splitlines() if line.strip() and not line.lstrip().startswith("#")]
        require(len(rows) >= 2 and all(len(row) == 2 for row in rows), "Bad table shape")
        require(all(k.is_finite() and energy.is_finite() and k > 0 and energy >= 0
                    for k, energy in rows), "Bad table values")
        require(all(left[0] < right[0] for left, right in zip(rows, rows[1:])), "Unordered k")
        if spectrum.get("INTERPOLATION") == "LogLinear":
            require(all(row[1] > 0 for row in rows), "LogLinear requires positive energy")
        if spectrum.get("HIGH") == "Gaussian":
            require(rows[-2][1] > rows[-1][1] > 0, "Gaussian tail needs a decreasing positive pair")
    elif kind == "VonKarman":
        required = ("TYPE", "LENGTH", "NORMALIZATION")
        if normalization == "Preserve":
            required += ("AMPLITUDE",)
        attrs(spectrum, required, optional)
        require(number(spectrum.get("LENGTH")) > 0, "Positive length required")
        if normalization == "Preserve":
            require(number(spectrum.get("AMPLITUDE")) >= 0, "Nonnegative amplitude required")
    else:
        raise ValueError("Unknown spectrum")
    if normalization == "TargetRMS":
        require(spectrum.get("RMS") is not None and number(spectrum.get("RMS")) >= 0,
                "Nonnegative target RMS required")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="Emit a machine-readable result")
    args = parser.parse_args()
    source_bytes = cbc.SOURCE.read_bytes()
    for name, expected in cbc.products(json.loads(source_bytes), source_bytes).items():
        require((HERE / name).read_bytes() == expected, "Stale generated fixture: " + name)
    examples = ("hit-inflow-tabulated.xml", "hit-inflow-vonkarman.xml")
    roots = [ET.parse(HERE / name).getroot() for name in examples]
    for root in roots:
        validate(root)
    mean_only = copy.deepcopy(roots[1])
    mean_only.find("HITINFLOW/REGION/SPECTRUM").set("RMS", "0")
    validate(mean_only)
    invalid = [
        (0, "HITINFLOW", "UNITS", "CGS"),
        (0, "HITINFLOW/REGION", "REF", "0"),
        (0, "HITINFLOW/REGION", "REF", str(1 << 31)),
        (0, "HITINFLOW/REGION/SYNTHESIS", "SEED", "-1"),
        (0, "HITINFLOW/REGION/SYNTHESIS", "SEED", str(1 << 64)),
        (0, "HITINFLOW/REGION/SYNTHESIS", "SEED", "1.5"),
        (0, "HITINFLOW/REGION/SYNTHESIS", "NX", "127"),
        (0, "HITINFLOW/REGION/SYNTHESIS", "LX", "0"),
        (0, "HITINFLOW/REGION/SYNTHESIS", "NZ", "16"),
        (0, "HITINFLOW/REGION/SYNTHESIS", "RECORDDURATION", "0.2"),
        (0, "HITINFLOW/REGION/CUTOFF", "KMAX", "400"),
        (0, "HITINFLOW/REGION/CUTOFF", "TAPERSTART", "201"),
        (0, "HITINFLOW/REGION/SPECTRUM", "FILE", "cbc1971_theta42_longitudinal_1d.dat"),
        (0, "HITINFLOW/REGION/SPECTRUM", "LOW", "Guess"),
        (0, "HITINFLOW/REGION/SPECTRUM", "TYPO", "1"),
        (1, "HITINFLOW/REGION/SPECTRUM", "AMPLITUDE", "1"),
        (0, "CONDITIONS/BOUNDARYCONDITIONS/REGION/N", "USERDEFINEDTYPE", "TimeDependent"),
    ]
    for index, path, key, value in invalid:
        root = copy.deepcopy(roots[index])
        root.find(path).set(key, value)
        try:
            validate(root)
        except ValueError:
            continue
        raise ValueError("Failed to reject " + path + " " + key + "=" + value)
    for tag, value in (("ORDER", "3"), ("VARIANT", "DIRK")):
        root = copy.deepcopy(roots[0])
        scheme = root.find("CONDITIONS/TIMEINTEGRATIONSCHEME")
        node = scheme.find(tag)
        if node is None:
            node = ET.SubElement(scheme, tag)
        node.text = value
        try:
            validate(root)
        except ValueError:
            continue
        raise ValueError("Failed to reject unsupported IMEX " + tag)
    additional = []
    root = copy.deepcopy(roots[0])
    root.append(copy.deepcopy(root.find("CONDITIONS")))
    additional.append(root)
    root = copy.deepcopy(roots[0])
    condition = root.find("CONDITIONS")
    condition.append(copy.deepcopy(condition.find("TIMEINTEGRATIONSCHEME")))
    additional.append(root)
    root = copy.deepcopy(roots[0])
    conditions = root.find("CONDITIONS")
    ET.SubElement(conditions.find("BOUNDARYREGIONS"), "B", ID="8").text = "C[2]"
    duplicate = copy.deepcopy(conditions.find("BOUNDARYCONDITIONS/REGION"))
    duplicate.set("REF", "8")
    conditions.find("BOUNDARYCONDITIONS").append(duplicate)
    additional.append(root)
    for parameter, value in (("CFL", "0.5"), ("FinTime", "1"),
                             ("IO_CheckTime", "0.1"), ("TimeIncrementFactor", "1.1")):
        root = copy.deepcopy(roots[0])
        ET.SubElement(root.find("CONDITIONS/PARAMETERS"), "P").text = parameter + "=" + value
        additional.append(root)
    for root in additional:
        try:
            validate(root)
        except ValueError:
            continue
        raise ValueError("Failed to reject extra configuration or variable time step")
    rejected = len(invalid) + 2 + len(additional)
    result = {"status": "PASS", "examples": len(examples), "mean_only_case": "PASS",
              "rejection_cases": rejected,
              "data_manifest_sha256": hashlib.sha256((HERE / "cbc1971-manifest.json").read_bytes()).hexdigest(),
              "example_sha256": {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest() for name in examples}}
    print(json.dumps(result, sort_keys=True, indent=2) if args.json else
          "PASS: two XML contracts, mean-only limit, {} rejection cases, and prepared fixture digests".format(rejected))


if __name__ == "__main__":
    main()
