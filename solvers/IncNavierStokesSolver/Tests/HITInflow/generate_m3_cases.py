#!/usr/bin/env python3
"""Generate small, complete full-3D HITInflow integration fixtures.

Only the Python standard library is needed. This generates inputs and an
independent reference; it does not run a solver or establish numerical accuracy.
The existing Milestone 1 experimental fixtures are never read or overwritten.
"""

import argparse
import cmath
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET


SEED = 20261028
RMS = 0.05
FINAL_TIME = 0.1
COUNTS = (4, 4, 4)
TWO_PI = 2.0 * math.acos(-1.0)
PERIODS = (TWO_PI, TWO_PI, TWO_PI / 4.0)
MASK64 = (1 << 64) - 1
VARIABLES = ("u", "v", "w", "p")


def number(value):
    return format(value, ".17g")


def sub(parent, name, text=None, **attributes):
    element = ET.SubElement(parent, name, attributes)
    if text is not None:
        element.text = str(text)
    return element


def mix64(value):
    value = (value + 0x9E3779B97F4A7C15) & MASK64
    value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & MASK64
    return value ^ (value >> 31)


def gaussian(index):
    """Independent implementation of the frozen SplitMix64BoxMuller52V1 rule."""
    mode_id = 0
    for component, count in zip(index, COUNTS):
        mode_id = mode_id * count + component + count // 2
    words = [mix64((SEED + 6 * mode_id + q) & MASK64) for q in range(6)]
    uniforms = [((word >> 12) + 0.5) * 2.0**-52 for word in words]
    xi = []
    for j in range(3):
        radius = math.sqrt(-math.log(uniforms[2 * j]))
        angle = TWO_PI * uniforms[2 * j + 1]
        xi.append(radius * complex(math.cos(angle), math.sin(angle)))
    return mode_id, words, xi


def reference_modes():
    """Only |k|=sqrt(2) survives this table and anisotropic lattice.

    Two independent modes have identical scalar weight. TargetRMS therefore
    gives weight * energyScale = 3*RMS**2/8, independently of table amplitude.
    Python/C++ libm and long-double rounding need not be bitwise identical.
    """
    result = []
    factor = math.sqrt(3.0 * RMS * RMS / 8.0)
    for index in ((1, -1, 0), (1, 1, 0)):
        mode_id, words, xi = gaussian(index)
        k = math.hypot(*index)
        unit = [component / k for component in index]
        coefficients = [
            factor * (xi[1] * unit[2] - xi[2] * unit[1]),
            factor * (xi[2] * unit[0] - xi[0] * unit[2]),
            factor * (xi[0] * unit[1] - xi[1] * unit[0]),
        ]
        result.append((index, coefficients, mode_id, words))
    return result


def geometry(root, z_period):
    """Four affine hexes: two cells in x, two in y, one in z."""
    geo = sub(root, "GEOMETRY", DIM="3", SPACE="3")
    vertices = sub(geo, "VERTEX")
    coordinates = []
    node_ids = {}
    for iz in range(2):
        for iy in range(3):
            for ix in range(3):
                point = (ix / 2.0, iy * math.pi, iz * z_period)
                node_ids[ix, iy, iz] = len(coordinates)
                coordinates.append(point)
                sub(vertices, "V", " ".join(map(number, point)),
                    ID=str(len(coordinates) - 1))
    edge_block, face_block, element_block = (sub(geo, name) for name in
                                            ("EDGE", "FACE", "ELEMENT"))
    edges, faces = {}, {}
    # Hex face numbering follows Hex_channel_m3.xml and HexGeom's convention.
    local_faces = ((0, 1, 2, 3), (0, 1, 5, 4), (1, 2, 6, 5),
                   (2, 3, 7, 6), (3, 0, 4, 7), (4, 5, 6, 7))
    boundary_faces = {name: [] for name in
                      ("inlet", "outlet", "ylow", "yhigh", "zlow", "zhigh")}
    for iy in range(2):
        for ix in range(2):
            nodes = [node_ids[ix + dx, iy + dy, dz]
                     for dx, dy, dz in ((0, 0, 0), (1, 0, 0), (1, 1, 0),
                                        (0, 1, 0), (0, 0, 1), (1, 0, 1),
                                        (1, 1, 1), (0, 1, 1))]
            element_faces = []
            for local in local_faces:
                cycle = [nodes[n] for n in local]
                key = tuple(sorted(cycle))
                if key not in faces:
                    face_edges = []
                    for a, b in zip(cycle, cycle[1:] + cycle[:1]):
                        edge_key = tuple(sorted((a, b)))
                        if edge_key not in edges:
                            edges[edge_key] = len(edges)
                            sub(edge_block, "E", f"{a} {b}",
                                ID=str(edges[edge_key]))
                        face_edges.append(edges[edge_key])
                    face_id = faces[key] = len(faces)
                    sub(face_block, "Q", " ".join(map(str, face_edges)),
                        ID=str(face_id))
                    for axis, lower, upper, width in (
                            (0, "inlet", "outlet", 1.0),
                            (1, "ylow", "yhigh", TWO_PI),
                            (2, "zlow", "zhigh", z_period)):
                        if all(coordinates[n][axis] == 0.0 for n in cycle):
                            boundary_faces[lower].append(face_id)
                        if all(coordinates[n][axis] == width for n in cycle):
                            boundary_faces[upper].append(face_id)
                element_faces.append(faces[key])
            sub(element_block, "H", " ".join(map(str, element_faces)),
                ID=str(iy * 2 + ix))
    composites = sub(geo, "COMPOSITE")
    sub(composites, "C", "H[0-3]", ID="0")
    for composite_id, ids in enumerate(boundary_faces.values(), 1):
        sub(composites, "C", "F[" + ",".join(map(str, ids)) + "]",
            ID=str(composite_id))
    sub(geo, "DOMAIN", "C[0]")


def expression(modes, component, derivative=None, k_squared=False):
    terms = []
    for mode, (_, coefficients, _, _) in enumerate(modes):
        phase = "x-y-Ubar*t" if mode == 0 else "x+y-Ubar*t"
        real, imag = (f"A{mode}{component}{suffix}" for suffix in ("r", "i"))
        if derivative is None:
            term = f"2*({real}*cos({phase})-{imag}*sin({phase}))"
            if k_squared:
                term = "2*(" + term + ")"
        else:
            wave = modes[mode][0][derivative]
            if wave == 0:
                continue
            term = f"(-2*({wave}))*({real}*sin({phase})+{imag}*cos({phase}))"
        terms.append(term)
    return "(" + "+".join(terms) + ")" if terms else "0"


def function(conditions, name, expressions):
    block = sub(conditions, "FUNCTION", NAME=name)
    for var, value in zip(VARIABLES, expressions):
        sub(block, "E", VAR=var, VALUE=value)


def make_case(kind, dt=0.001, num_modes=6):
    manufactured = kind == "manufactured"
    mean_only = kind == "mean-only"
    z_period = PERIODS[2] if manufactured or mean_only else TWO_PI
    end_time = FINAL_TIME if manufactured else 0.02
    modes = reference_modes() if manufactured else []
    root = ET.Element("NEKTAR")
    root.append(ET.Comment(" Generated by generate_m3_cases.py; full 3D, sparse inlet REF=7. "))
    geometry(root, z_period)
    expansions = sub(root, "EXPANSIONS")
    sub(expansions, "E", COMPOSITE="C[0]", NUMMODES=str(num_modes),
        FIELDS=",".join(VARIABLES), TYPE="MODIFIED")
    conditions = sub(root, "CONDITIONS")
    info = sub(conditions, "SOLVERINFO")
    for key, value in (("SolverType", "VelocityCorrectionScheme"),
                       ("EQTYPE", "UnsteadyNavierStokes"),
                       ("EvolutionOperator", "Nonlinear"),
                       ("Projection", "Continuous"),
                       ("AdvectionForm", "Convective"),
                       ("Extrapolation", "Standard")):
        sub(info, "I", PROPERTY=key, VALUE=value)
    scheme = sub(conditions, "TIMEINTEGRATIONSCHEME")
    sub(scheme, "METHOD", "IMEX")
    sub(scheme, "ORDER", "2")
    parameters = sub(conditions, "PARAMETERS")
    count = round(end_time / dt)
    if not math.isclose(count * dt, end_time, abs_tol=1e-15):
        raise ValueError("Timestep must divide the record duration")
    for name, value in (("Ubar", 1.0), ("Kinvis", 0.01), ("TimeStep", dt),
                        ("NumSteps", count), ("IO_InfoSteps", count),
                        ("IO_CheckSteps", 0), ("FinTime", 0), ("CFL", 0),
                        ("IO_CheckTime", 0), ("Flowrate", 0)):
        sub(parameters, "P", f"{name} = {number(value)}")
    for mode, (_, coefficients, _, _) in enumerate(modes):
        for component, value in enumerate(coefficients):
            for suffix, part in (("r", value.real), ("i", value.imag)):
                sub(parameters, "P", f"A{mode}{component}{suffix} = {number(part)}")
    variables = sub(conditions, "VARIABLES")
    for n, name in enumerate(VARIABLES):
        sub(variables, "V", name, ID=str(n))
    regions = sub(conditions, "BOUNDARYREGIONS")
    for region, composite in ((7, 1), (11, 2), (21, 3), (22, 4), (31, 5), (32, 6)):
        sub(regions, "B", f"C[{composite}]", ID=str(region))
    boundaries = sub(conditions, "BOUNDARYCONDITIONS")
    inlet = sub(boundaries, "REGION", REF="7")
    for var, value in zip(VARIABLES[:3], ("Ubar", "0", "0")):
        sub(inlet, "D", VAR=var, USERDEFINEDTYPE="HITInflow", VALUE=value)
    sub(inlet, "N", VAR="p", USERDEFINEDTYPE="H", VALUE="0")
    outlet = sub(boundaries, "REGION", REF="11")
    for component, var in enumerate(VARIABLES[:3]):
        if manufactured:
            sub(outlet, "N", VAR=var, USERDEFINEDTYPE="TimeDependent",
                VALUE=expression(modes, component, derivative=0))
        else:
            sub(outlet, "N", VAR=var, USERDEFINEDTYPE="HOutflow", VALUE="0")
    pressure_options = {} if manufactured else {"USERDEFINEDTYPE": "HOutflow"}
    sub(outlet, "D", VAR="p", VALUE="0", **pressure_options)
    for region, partner in ((21, 22), (22, 21), (31, 32), (32, 31)):
        lateral = sub(boundaries, "REGION", REF=str(region))
        for var in VARIABLES:
            sub(lateral, "P", VAR=var, VALUE=f"[{partner}]")
    velocity = [expression(modes, i) for i in range(3)]
    velocity[0] = "Ubar+" + velocity[0] if manufactured else "Ubar"
    function(conditions, "InitialConditions", velocity + ["0"])
    if manufactured or mean_only:
        function(conditions, "ExactSolution", velocity + ["0"])
    if manufactured:
        force = []
        for component in range(3):
            advection = "+".join(
                expression(modes, j) + "*" + expression(modes, component, derivative=j)
                for j in range(2))
            force.append(advection + "+Kinvis*" +
                         expression(modes, component, k_squared=True))
        function(conditions, "BodyForce", force)
        forcing = sub(root, "FORCING")
        sub(sub(forcing, "FORCE", TYPE="Body"), "BODYFORCE", "BodyForce")
    hit = sub(root, "HITINFLOW", VERSION="1", UNITS="SI")
    region = sub(hit, "REGION", REF="7")
    sub(region, "ORIGIN", X="0", Y="0", Z="0")
    sub(region, "SYNTHESIS", LX=number(TWO_PI), LY=number(TWO_PI),
        LZ=number(z_period), NX="4", NY="4", NZ="4", SEED=str(SEED),
        TIMEORIGIN="0", RECORDDURATION=number(end_time), REPEAT="Error")
    if manufactured:
        sub(region, "SPECTRUM", TYPE="Tabulated", FILE="m3-narrow-band.dat",
            INTERPOLATION="Linear", LOW="Zero", HIGH="Zero",
            NORMALIZATION="TargetRMS", RMS=number(RMS))
    else:
        sub(region, "SPECTRUM", TYPE="VonKarman", LENGTH="0.8",
            NORMALIZATION="TargetRMS", RMS="0" if mean_only else number(RMS))
    sub(region, "CUTOFF", KMAX="1.55", TAPERSTART="1.5")
    ET.indent(root, space="    ")
    return ET.tostring(root, encoding="unicode", xml_declaration=True) + "\n"


def write_reference(output):
    modes = reference_modes()
    covariance = [[sum(2 * (a[i] * a[j].conjugate()).real for _, a, _, _ in modes)
                   for j in range(3)] for i in range(3)]
    document = {
        "purpose": "Independent fixture oracle; not a solver validation result",
        "algorithm": "SplitMix64BoxMuller52V1",
        "seed": str(SEED), "target_ensemble_rms": RMS,
        "periods": PERIODS, "counts": COUNTS,
        "retained_modes": [{"index": k, "mode_id": mode_id,
                            "counter_words_hex": [f"{word:016x}" for word in words],
                            "coefficient_real_imag": [[v.real, v.imag] for v in a]}
                           for k, a, mode_id, words in modes],
        "expected_kinetic_energy": 1.5 * RMS * RMS,
        "realized_covariance": covariance,
        "realized_kinetic_energy": sum(covariance[i][i] for i in range(3)) / 2,
        "samples": [],
        "acceptance": {
            "native_coefficient_scaled_absolute_tolerance": 2e-14,
            "temporal_order_minimum": 1.8,
            "mean_only_l2_linf_maximum": 1e-10,
        },
    }
    for point, time in (((0.0, 0.37, 0.2), 0.0), ((0.0, 1.2, 0.4), 0.073),
                        ((0.41, 2.1, 0.0), FINAL_TIME)):
        values = [1.0, 0.0, 0.0]
        derivatives = [[0.0] * 3 for _ in range(3)]
        for k, coefficients, _, _ in modes:
            exponential = cmath.exp(1j * (sum(q * x for q, x in zip(k, point)) - time))
            for i, coefficient in enumerate(coefficients):
                value = coefficient * exponential
                values[i] += 2.0 * value.real
                for j in range(3):
                    derivatives[i][j] -= 2.0 * k[j] * value.imag
        document["samples"].append({"point": point, "time": time,
                                    "velocity": values, "gradient": derivatives})
    (output / "m3-manufactured-reference.json").write_text(
        json.dumps(document, indent=2) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--mms-modes", type=int, default=12,
                        help="Expansion NUMMODES for all manufactured cases (default: 12)")
    args = parser.parse_args()
    if args.mms_modes < 4:
        parser.error("--mms-modes must be at least four")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    spectrum = ("# NEKTAR_HIT_SPECTRUM_V1 KIND=Shell3D K_UNITS=rad/m E_UNITS=m3/s2\n"
                "# Manufactured integration-test spectrum; no experimental values.\n"
                "# Linear interpolation, zero outside table; only |k|=sqrt(2) survives.\n"
                f"1.1 0\n{number(math.sqrt(2.0))} 1\n1.6 0\n")
    (args.output_dir / "m3-narrow-band.dat").write_text(spectrum, encoding="utf-8")
    for dt in (0.01, 0.005, 0.0025, 0.00125):
        name = f"m3-manufactured-dt{dt:.5f}.xml"
        (args.output_dir / name).write_text(
            make_case("manufactured", dt=dt, num_modes=args.mms_modes), encoding="utf-8")
    for kind in ("mean-only", "stochastic"):
        (args.output_dir / f"m3-{kind}.xml").write_text(
            make_case(kind, num_modes=6 if kind == "mean-only" else 8), encoding="utf-8")
    write_reference(args.output_dir)
    print(f"Generated 6 complete sessions, a spectrum, and an oracle in {args.output_dir}")


if __name__ == "__main__":
    main()
