#!/usr/bin/env python3
"""Prepare an exact single-pair HITInflow manufactured regression.

This constructs inputs and an independent frozen-RNG oracle. It does not run a
solver or establish numerical accuracy. The unequal synthesis periods isolate
one transverse pair deliberately; this fixture is not a physical HIT ensemble.
"""

import argparse
import cmath
import importlib.util
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET

FIXTURES = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("hit_m3_fixtures", FIXTURES / "generate_m3_cases.py")
m3 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m3)

PERIODS = (m3.TWO_PI, 1.0, 1.0)
INDEX = (1, 0, 0)
RMS = 0.05
SEED = m3.SEED
FINAL_TIME = 0.1


def reference():
    mode_id, words, xi = m3.gaussian(INDEX)
    # A single independent member contributes K_expected=2*weight*scale.
    # TargetRMS thus gives weight*scale=3*RMS^2/4. The x-directed wave
    # has xi cross e_x = (0, xi_z, -xi_y).
    factor = math.sqrt(3.0 * RMS * RMS / 4.0)
    return (0j, factor * xi[2], -factor * xi[1]), mode_id, words


def parameter(conditions, name, value):
    entries = conditions.find("PARAMETERS")
    for item in entries:
        if item.text.split("=", 1)[0].strip() == name:
            item.text = f"{name} = {m3.number(value)}"
            return
    m3.sub(entries, "P", f"{name} = {m3.number(value)}")


def wave(component, derivative=False):
    if component == 0:
        return "0"
    real, imag = f"A{component}r", f"A{component}i"
    phase = "x-Ubar*t"
    if derivative:
        return f"(-2)*({real}*sin({phase})+{imag}*cos({phase}))"
    return f"2*({real}*cos({phase})-{imag}*sin({phase}))"


def make_case(dt=0.005, modes=8):
    steps = round(FINAL_TIME / dt)
    if dt <= 0 or not math.isclose(steps * dt, FINAL_TIME, abs_tol=1e-15):
        raise ValueError("The positive timestep must divide 0.1 seconds")
    if modes < 4:
        raise ValueError("At least four polynomial modes are required")
    root = ET.fromstring(m3.make_case("mean-only", dt=0.005, num_modes=modes))
    root.insert(0, ET.Comment(" Generated single Fourier-pair manufactured test; not an isotropic statistical ensemble. "))
    for vertex in root.findall("./GEOMETRY/VERTEX/V"):
        x, y, z = map(float, vertex.text.split())
        vertex.text = " ".join(map(m3.number, (x, y / m3.TWO_PI, z / m3.PERIODS[2])))
    conditions = root.find("CONDITIONS")
    m3.sub(conditions.find("SOLVERINFO"), "I", PROPERTY="GlobalSysSoln", VALUE="DirectStaticCond")
    for name, value in (("TimeStep", dt), ("NumSteps", steps), ("IO_InfoSteps", steps)):
        parameter(conditions, name, value)
    coefficients, _, _ = reference()
    for component, coefficient in enumerate(coefficients):
        for suffix, value in (("r", coefficient.real), ("i", coefficient.imag)):
            parameter(conditions, f"A{component}{suffix}", value)
    for function in list(conditions.findall("FUNCTION")):
        conditions.remove(function)
    velocity = ["Ubar", wave(1), wave(2)]
    m3.function(conditions, "InitialConditions", velocity + ["0"])
    m3.function(conditions, "ExactSolution", velocity + ["0"])
    m3.function(conditions, "BodyForce", ["0", "Kinvis*(" + wave(1) + ")", "Kinvis*(" + wave(2) + ")"])
    forcing = m3.sub(root, "FORCING")
    m3.sub(m3.sub(forcing, "FORCE", TYPE="Body"), "BODYFORCE", "BodyForce")
    outlet = conditions.find("./BOUNDARYCONDITIONS/REGION[@REF='11']")
    for component, condition in enumerate(outlet.findall("N")):
        condition.set("USERDEFINEDTYPE", "TimeDependent")
        condition.set("VALUE", wave(component, derivative=True))
    outlet.find("D").attrib.pop("USERDEFINEDTYPE", None)
    synthesis = root.find("./HITINFLOW/REGION/SYNTHESIS")
    synthesis.set("LY", "1")
    synthesis.set("LZ", "1")
    synthesis.set("RECORDDURATION", m3.number(FINAL_TIME))
    root.find("./HITINFLOW/REGION/SPECTRUM").set("RMS", m3.number(RMS))
    ET.indent(root, space="    ")
    return ET.tostring(root, encoding="unicode", xml_declaration=True) + "\n"


def write_reference(path):
    coefficients, mode_id, words = reference()
    document = {
        "purpose": "Independent input oracle; not a solver validation result",
        "physical_scope": "One transverse pair isolated on an unequal-period lattice; not a physical HIT ensemble",
        "algorithm": "SplitMix64BoxMuller52V1", "seed": str(SEED),
        "counts": (4, 4, 4), "periods": PERIODS,
        "canonical_index": INDEX, "wavevector": (1.0, 0.0, 0.0),
        "mode_id": mode_id, "counter_words_hex": [f"{word:016x}" for word in words],
        "coefficient_real_imag": [[a.real, a.imag] for a in coefficients],
        "phase": "x-Ubar*t", "mean_speed": 1.0,
        "velocity": "(Ubar, 2*Re(a_y*exp(i*phase)), 2*Re(a_z*exp(i*phase)))",
        "pressure": 0, "body_force": "(0, Kinvis*v, Kinvis*w)",
        "expected_kinetic_energy": 1.5 * RMS * RMS,
        "target_ensemble_rms": RMS, "samples": [],
        "acceptance": {
            "active_independent_modes": 1,
            "scaled_native_coefficient_tolerance": 2e-13,
            "normalized_modal_transversality_maximum": 1e-12,
            "normalized_sample_divergence_maximum": 1e-12,
            "normalized_finite_difference_derivative_error_maximum": 1e-8,
            "solver_L2_u_p_maximum": 1e-9,
            "solver_Linf_u_p_maximum": 1e-9,
            "solver_L2_v_w_maximum": 2e-6,
            "solver_Linf_v_w_maximum": 3e-6,
        },
    }
    for x, time in ((0.0, 0.0), (0.0, 0.05), (0.0, 0.1), (0.37, 0.073)):
        exponential = cmath.exp(1j * (x - time))
        values = [1.0, 0.0, 0.0]
        gradient = [[0.0] * 3 for _ in range(3)]
        for component, a in enumerate(coefficients):
            value = a * exponential
            values[component] += 2.0 * value.real
            gradient[component][0] = -2.0 * value.imag
        document["samples"].append({"point": (x, 0.37, 0.2), "time": time,
                                    "velocity": values, "gradient": gradient,
                                    "time_derivative": [-row[0] for row in gradient]})
    path.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=FIXTURES)
    parser.add_argument("--dt", type=float, default=0.005)
    parser.add_argument("--modes", type=int, default=8)
    args = parser.parse_args()
    try:
        xml = make_case(args.dt, args.modes)
    except ValueError as error:
        parser.error(str(error))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "m5-single-mode.xml").write_text(xml, encoding="utf-8")
    write_reference(args.output_dir / "m5-single-mode-reference.json")
    print(f"Prepared one single-mode session and independent oracle in {args.output_dir}")


if __name__ == "__main__":
    main()
