#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Prepare and independently audit the affected native Task 2R paths."""

import argparse
import json
import math
from pathlib import Path
import struct
import time
import xml.etree.ElementTree as ET

import numpy as np
from numpy.polynomial import chebyshev as cheb
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad

import coupled_fixtures as sessions
import coupled_oracle as core
import coupled_return_run as returns
import physical_oracle as po
from resources import atomic, identity
from task2r_reference import AIN, B, C, CD, D0, H, HS, KAPPA, KIN, NE, NU, POLICY, QE, RT, TAU0, UTAU
from task2r_reference import Column, blend, digest, require, tensor


def read(path):
    return json.loads(Path(path).read_text(), parse_int=lambda s: -0. if s == "-0" else int(s))


def edit(path, info, parameters):
    root = ET.parse(path)
    si = root.find("CONDITIONS/SOLVERINFO")
    for key, value in info.items():
        for item in list(si.findall(f"I[@PROPERTY='{key}']")):
            si.remove(item)
        ET.SubElement(si, "I", PROPERTY=key, VALUE=str(value))
    block = root.find("CONDITIONS/PARAMETERS")
    for item in list(block):
        if item.text.split("=")[0].strip() in parameters:
            block.remove(item)
    for key, value in parameters.items():
        ET.SubElement(block, "P").text = f"{key} = {value}"
    ET.indent(root)
    root.write(path, encoding="unicode")


def mms_mesh(path, face, refined=False):
    xs, rs = [0., .5, 1.], [0., face, H, .325, .5]
    if refined:
        xs = sorted(xs + [(a + b) / 2 for a, b in zip(xs[:-1], xs[1:])])
        rs = sorted(rs + [(a + b) / 2 for a, b in zip(rs[:-1], rs[1:])])
    sessions.mesh(path, xs, rs, quadrature=12)
    root = ET.parse(path)
    for vertex in root.findall("GEOMETRY/VERTEX/V"):
        x, z, y = map(float, vertex.text.split())
        vertex.text = f"{x:.17g} {z + D0:.17g} {y:.17g}"
    root.write(path, encoding="unicode")


def prepare(out, table, historical, geometry):
    out.mkdir()
    ref = po.load_reference(table)
    require(po.load_reference(historical).identity == po.TABLE_SHA, "historical fixture table differs")
    when = .28
    cross_K = KIN * (1 + .012 * math.sin(.6 * when) + .004 * math.cos(.8 * when)) ** 2
    face = ref.cutoff(np.array([cross_K]))[0]
    template = out / "crossing-template.xml"
    refined = out / "crossing-refined-template.xml"
    mms_mesh(template, face)
    mms_mesh(refined, face, refined=True)
    source = dict(exact=table, historical=historical)
    recipe = {}

    def fixture(name, candidate="exact", dt=.02, steps=2, recovery=0,
                variant="Historical", restart="", q=(8, 12), tolerances=(3e-13, 3e-12),
                physical=False, mesh=template):
        xml, prefix = out / (name + ".xml"), out / name
        recipe[name] = sessions.generate(xml, source[candidate], prefix,
            mode="Physical" if physical else "Manufactured", dt=dt, steps=steps, recovery=recovery,
            quadrature_x=q[0], quadrature_r=q[1], restart=restart, template=geometry if physical else mesh,
            backend="FastFull3DV1", mms_body_amplitude=0 if physical else 1, mms_variant=variant,
            linear=tolerances[0], nonlinear=tolerances[1], check_current_actions=not physical)
        options = {"WJImplicitSolver": "CoupledJFNKV1", "WJStepPreconditioner": "CurrentDataDriftV1"}
        if candidate == "exact":
            options["WJReferencePolicy"] = POLICY
        params = {"WJStepNewtonSeconds": 30., "WJStepPCRelativeDrift": .02, "WJStepPCMaximumUses": 8}
        edit(xml, options, params)
        recipe[name]["session"] = identity(xml)

    fixture("initial", steps=0, physical=True)
    fixture("hist", steps=18, recovery=10)
    fixture("old-reference", candidate="historical")
    fixture("old-reference-installed", candidate="historical")
    fixture("old-to-new", steps=1, restart=out / "old-reference.final.bin")
    fixture("new-to-old", candidate="historical", steps=1, restart=out / "hist.final.bin")
    fixture("same-policy-stale-epsilon", steps=1, restart=out / "stale-epsilon.bin")
    fixture("wrong-policy-artifact")
    fixture("wrong-calibration-artifact")
    edit(out / "wrong-policy-artifact.xml", {"WJReferenceTable": historical}, {})
    wrong = read(table)
    wrong["geometry"]["force_support"] *= .999
    wrong["bindings"]["geometry"] = digest(wrong["geometry"])
    atomic(out / "wrong-calibration.json", wrong)
    edit(out / "wrong-calibration-artifact.xml", {"WJReferenceTable": out / "wrong-calibration.json"}, {})
    cases = [("temporal-08", .08), ("temporal-04", .04), ("temporal-02", .02)]
    for name, dt in cases:
        fixture(name, dt=dt, steps=round(.4 / dt), variant="PolynomialScalarsV1")
    fixture("temporal-q", steps=20, variant="PolynomialScalarsV1", q=(12, 18))
    fixture("temporal-tight", steps=20, variant="PolynomialScalarsV1", tolerances=(1e-13, 1e-12))
    fixture("temporal-h", steps=20, variant="PolynomialScalarsV1", mesh=refined)
    controls = [("mesh", "temporal-h"), ("quadrature", "temporal-q"), ("tolerance", "temporal-tight")]
    atomic(out / "temporal-study.json", dict(schema="WJH27MMSStudyManifestV1", terminal_time=.4,
        runs=[dict(label=name, role=role, native=str(out / (name + ".result.json")))
              for name, role in [(name, "temporal") for name, _ in cases] +
              [(name, "control") for _, name in controls]],
        comparisons=[dict(label=kind, kind=kind, coarse="temporal-02", fine=name)
                     for kind, name in controls]))
    atomic(out / "fixtures.json", dict(schema="WJH27Task2RAffectedInputsV1", crossing_face=float(face),
        retained_geometry=identity(geometry), exact_table=identity(table), historical_table=identity(historical),
        fixtures=recipe, observations="All native projection owner endpoints and interior Gauss samples."))


def reference(table, exact, legacy, out):
    ref = po.load_reference(table)
    actual, prior = read(exact), read(legacy)
    require(actual["table_sha256"] == ref.identity and actual["reference_policy"] == POLICY,
            "native exact-inner reference binding differs")
    errors = {}

    def check(label, val, target, scale, limit=1e-10):
        error = float(np.max(np.abs(np.asarray(val) - target) / scale))
        errors[label] = max(error, errors.get(label, 0.))
        require(error <= limit, f"native fixed scale {label} error {error:.6g} exceeds {limit}")

    x0, x1 = ref.lo, ref.hi
    slope = cheb.chebder(ref.coefficients) * 2 / (x1 - x0)
    for row in actual["dynamic"]:
        K, r = row["K"], row["r"]
        c = ref.cutoff(np.array([K]))[0]
        dc = c * cheb.chebval((2 * K - x0 - x1) / (x1 - x0), slope)
        A, Ap = K ** 1.5 / ref.cl, 1.5 * math.sqrt(K) / ref.cl
        check("cutoff", row["c"], c, H)
        check("cutoff_K", row["dc"], dc, c / KIN)
        check("cutoff_inverse", row["K_at_c"], K, KIN, 1e-9)
        check("amplitude", [row["A"], row["dA"], row["ddA"]],
              np.array([A, Ap, .75 / (math.sqrt(K) * ref.cl)]),
              np.array([A, A / KIN, A / KIN ** 2]))

        def aux(k):
            x = HS * np.sqrt(k) / (70 * NU)
            return KAPPA * np.sqrt(k / QE) * (r + 10.7 * NU / np.sqrt(k) *
                       np.expm1((1 - .79 * x ** (-.28)) * np.log(x)))

        a = aux(K)
        check("auxiliary", row["aux"], np.array([a, aux(K + 1e-28j).imag / 1e-28,
                                                 KAPPA * math.sqrt(K / QE)]),
              np.array([a, a / KIN, a / (r + .002)]))
        # An exact represented kink owns the plateau side. Adjacent one-sided
        # radii refer to the emitted native cutoff, whose value is checked above.
        E = A / max(r, c)
        expected = (E, E * (Ap / A - dc / c), 0.) if r <= row["c"] else (E, Ap / r, -E / r)
        check("dependent_epsilon", row["epsilon"], np.array(expected), np.array([E, E / KIN, E / H]))

    radii = np.array([row["r"] for row in actual["inlet"]])
    values = ref.inlet(radii)
    for field, key, scale in (("U", "U", 1), ("k", "k", KIN), ("U_r", "G", 1),
                              ("k_r", "kr", KIN / H), ("Pi", "Pi", TAU0)):
        check(field, [row[field] for row in actual["inlet"]], values[key], scale)
    check("outlet_traction", [row["traction"] for row in actual["inlet"]], ref.outlet(radii), TAU0)
    check("top_tangent", actual["top_tangential_traction"], TAU0 + NU * UTAU / (KAPPA * RT), TAU0)
    check("top_epsilon_flux", actual["top_epsilon_outward_flux"], ref.top_flux, ref.top_flux)
    old_radii = {row["r"]: row for row in prior["inlet"] if row["r"] > H}
    for row in actual["inlet"]:
        r = row["r"]
        if r in old_radii:
            old = old_radii[r]
            for key, scale in (("U", 1.), ("k", KIN), ("U_r", 1.), ("k_r", KIN / H),
                               ("epsilon", AIN / r), ("epsilon_r", AIN / r ** 2)):
                check("outer_parity_" + key, row[key], old[key], scale)
    atomic(out, dict(status="PASS", actual=identity(exact), historical=identity(legacy),
                     fixed_scale_errors=errors, dynamic_cases=len(actual["dynamic"]),
                     inlet_cases=len(actual["inlet"]), rejected=actual["rejections"]))


def projection(table, native, out, *, outlet_x=15.):
    data, metadata = po.load_native(native)
    ref = po.load_reference(table)
    oracle = po.PhysicalOracle(data, ref, 28, time.monotonic() + 60, outlet_x=outlet_x)
    require(data["step"] == 0 and not data["manufactured"], "projection needs the fresh registered state")
    norms = po.Norms()
    samples = []
    gauss, weights = leggauss(34)
    for owner, side in oracle.faces["u"]:
        if side != 0:
            continue
        _, _, rl, rr = owner["bounds"]
        edges = oracle.cuts(rl, rr, ref.inlet_panel_bounds)
        for lo, hi in zip(edges[:-1], edges[1:]):
            r = np.r_[lo, (lo + hi) / 2 + (hi - lo) * gauss / 2,
                      lo + (hi - lo) * np.linspace(.001, .999, 63), hi]
            w = np.r_[0., weights * (hi - lo) / 2, np.zeros(64)]
            x = np.full_like(r, owner["bounds"][0])
            values = ref.inlet(r)
            U, k = oracle.s["flow"][owner["u"]], oracle.s["k"][owner["k"]]
            u, _, ur = po.basis(owner, "u", x, r)
            kbase = po.basis(owner, "k", x, r)[0]
            d = dict(U=u @ U - values["U"], k=(kbase @ k - values["k"]) / KIN,
                     G=ur @ U - values["G"])
            for key, value in d.items():
                norms.add("inlet", key, value, w)
                at = int(np.argmax(abs(value)))
                samples.append(dict(field=key, radius=float(r[at]), error=float(value[at]),
                                    owner=owner["geometry"], bounds=[rl, rr]))
    maxima = {key: max((r for r in samples if r["field"] == key), key=lambda row: abs(row["error"]))
              for key in ("U", "k", "G")}
    limits = dict(U=1e-5, k=1e-4, G=1e-3)
    result = norms.report(limits)
    atomic(out, dict(status="PASS" if all(v["status"] == "MET" for v in result.values()) else "NOT_MET",
                     errors=result, maxima=maxima, limits=limits, initial=identity(native),
                     H=oracle.geometry["H"], d0=oracle.geometry["d0"], physical_H_mm=40*H,
                     interface_owners=sum(abs(e["bounds"][2]-H)<1e-12 for e in oracle.elements),
                     basis="Independent Modified_A", gradient_scale="Uref/L0 in code coordinates"))
    print(json.dumps(result), flush=True)


def stale_histories(folder):
    a, b = folder / "hist.recovery-10.bin", folder / "hist.final.bin"
    a_header, a_parts = returns.history_parts(a)
    _, b_parts = returns.history_parts(b)
    require(a_parts[2] != b_parts[2] and len(a_parts[2]) == len(b_parts[2]), "stale explicit control has no change")

    def write(name, parts):
        import hashlib
        path = folder / name
        content = a_header + b"".join(parts)
        with path.open("xb") as stream:
            stream.write(content + hashlib.sha256(content).hexdigest().encode("ascii"))
        return identity(path)

    stale = write("stale-explicit-A.bin", [a_parts[0], a_parts[1], b_parts[2], a_parts[3]])
    data = read(folder / "hist.result.json")["state"]
    outer = next(e for e in data["elements"] if abs(e["bounds"][2]-H)<1e-12 and abs(e["bounds"][0]-.5)<1e-12)
    dof = outer["epsilon"][0]
    packet = bytearray(a_parts[0])
    offset = 0
    for _ in range(2):
        length, = struct.unpack_from("<Q", packet, offset)
        offset += 8 * (length + 1)
    start = offset + 8 * (1 + dof)
    value, = struct.unpack_from("<d", packet, start)
    struct.pack_into("<d", packet, start, value + .001 * AIN / H)
    write("stale-epsilon.bin", [bytes(packet), *a_parts[1:]])
    atomic(folder / "cache-inputs.json", dict(schema="WJH27CoupledReturnInputsV1",
                                             A=identity(a), B=identity(b), stale_A=stale))


def grade_projection(table, template, out):
    from numpy.polynomial.legendre import Legendre
    from physical_reference_refinement import axes
    from fixtures import modified

    out.mkdir()
    xs, zs = axes(template)
    ref = po.load_reference(table)
    polynomial = Legendre.basis(11)
    nodes = np.r_[-1., polynomial.deriv().roots(), 1.]
    weights = 2 / (12 * 11 * polynomial(nodes) ** 2)
    base = modified(3, nodes)[0]
    fixed_mass = base[:, 2:].T @ (weights[:, None] * base[:, 2:])
    check = np.linspace(-1., 1., 211)
    b, d = modified(3, check)
    results = []

    def refine(lo, hi, depth=0):
        radii = (lo + hi) / 2 + (hi - lo) * nodes / 2 - D0
        radii[[0, -1]] = [lo - D0, hi - D0]
        v = ref.inlet(radii)["U"]
        free = np.linalg.solve(fixed_mass, base[:, 2:].T @ (weights * (v - base[:, :2] @ v[[0, -1]])))
        coeffs = np.r_[v[[0, -1]], free]
        sample = (lo + hi) / 2 + (hi - lo) * check / 2 - D0
        sample[[0, -1]] = radii[[0, -1]]
        target = ref.inlet(sample)
        errors = dict(U=float(np.max(abs(b @ coeffs - target["U"]))),
                      G=float(np.max(abs(d @ coeffs * 2 / (hi - lo) - target["G"]))))
        accepted = errors["U"] <= .8e-5 and errors["G"] <= .8e-3
        results.append(dict(lo=lo-D0, hi=hi-D0, errors=errors, depth=depth, accepted=accepted))
        if accepted:
            return [hi]
        require(depth < 6, "projection grading exceeds bounded refinement")
        mid = (lo + hi) / 2
        return refine(lo, mid, depth + 1) + refine(mid, hi, depth + 1)

    graded = [zs[0]] + [v for lo, hi in zip(zs[:-1], zs[1:]) for v in refine(lo, hi)]
    geometry = out / "graded-template.xml"
    sessions.mesh(geometry, xs, graded, quadrature=12)
    cfg = out / "initial.xml"
    sessions.generate(cfg, table, out / "initial", mode="Physical", dt=.02, steps=0,
                      template=geometry, backend="FastFull3DV1")
    edit(cfg, {"WJReferencePolicy": POLICY}, {})
    atomic(out / "mesh-grading.json", dict(baseline=identity(template), graded=identity(geometry),
        cells=[len(xs)-1, len(graded)-1], x=xs, z=graded, projection=results,
        controls="Independent endpoint-constrained Q3 on GLL12. Selected with an 0.8 reserve and verified natively."))
    print("Candidate mesh", len(xs)-1, len(graded)-1, flush=True)


def continuum_actions(table, out):
    ref = po.load_reference(table)
    column = Column(KIN)
    edges = [p["domain"][0] for p in ref.inlet_panels] + [H]
    scalars = rsz = 0.
    for lo, hi in zip(edges[:-1], edges[1:]):
        z, w = leggauss(48)
        r, weight = (lo + hi) / 2 + (hi - lo) * z / 2, w * (hi - lo) / 2
        # Frozen fields with direct independent full-tensor and diffusion actions.
        v = ref.inlet(r)
        e = AIN / np.maximum(r, ref.inlet_cutoff)
        r2, prod = tensor(v["k"], e, v["G"], column.auxiliary(r), blend(r))
        d = NU + CD * v["k"] ** 2 / e
        scalars += weight @ (d * v["kr"] / H - (prod - e) * r / H)
        rsz += weight @ (2 * v["k"] / 3 - r2[:, 2, 2] - v["Pi"]) / H
    integrate = lambda fn: quad(fn, H, RT, epsabs=1e-17, epsrel=1e-12)[0]
    tx = lambda r: TAU0 + NU * UTAU / (KAPPA * r)
    momentum = column.shear * H / RT + integrate(lambda r: tx(r) / RT) - ref.top_tangent
    m_strong = integrate(lambda r: NU * UTAU / (KAPPA * r * r) * r / RT)
    phi = lambda r: (r - H) / (RT - H)
    eq = lambda r: -NU * AIN / r ** 2 - CD * KIN ** 2 / (1.3 * r)
    eweak = integrate(lambda r: eq(r) / (RT - H) - (1.44 - 1.82) * (AIN / r) ** 2 / KIN * phi(r)) + ref.top_flux
    estrong = integrate(lambda r: -2 * NU * AIN / r ** 3 * phi(r))
    outflow = ref.outlet(np.array([H, np.nextafter(H, RT), RT]))
    defects = dict(inner_k=abs(scalars) / AIN, inner_normal_stress=abs(rsz - KIN / NE) / TAU0,
                   axial_momentum=abs(momentum - m_strong) / TAU0,
                   outer_epsilon=abs(eweak - estrong) / AIN,
                   interface_traction=float(np.max(abs(outflow[0] - outflow[1])) / TAU0),
                   corner_tangent=abs(outflow[-1, 1] - ref.top_tangent) / TAU0)
    require(max(defects.values()) <= 1e-10, f"independent unforced continuum weak actions {defects}")
    atomic(out, dict(status="PASS", fixed_scale_defects=defects, loaded_outlet=outflow.tolist(),
                     weak_actions=dict(momentum=momentum, k=scalars, epsilon=eweak),
                     strong_rhs_integrals=dict(momentum=-integrate(lambda r: NU*UTAU/(KAPPA*r*r)),
                                               k=0., epsilon=integrate(lambda r: 2*NU*AIN/r**3)),
                     convention="Weak steady residual is minus strong RHS, including the exterior weak loads."))


def direct_fixture(out, table):
    out.mkdir()
    path = out / "case.xml"
    sessions.generate(path, table, out / "case", nx=1, nr=1, dt=.005, steps=2,
                      backend="FastFull3DV1", mms_body_amplitude=1., check_current_actions=True,
                      linear=3e-13, nonlinear=3e-12, mms_variant="PolynomialScalarsV1")
    edit(path, {"WJReferencePolicy": POLICY, "WJImplicitSolver": "CoupledJFNKV1",
                "WJStepPreconditioner": "CurrentDataDriftV1"}, {"WJStepNewtonSeconds": 30.})


def temporal_mesh_control(out, table):
    """Repeat the inherited temporal mesh and preserve the separate crossing case."""
    out.mkdir()
    cases = [("dt08", .08, 4, (8, 12), (1e-12, 3e-11), "temporal"),
             ("dt04", .04, 4, (8, 12), (1e-12, 3e-11), "temporal"),
             ("dt02", .02, 4, (8, 12), (1e-12, 3e-11), "temporal"),
             ("h", .02, 8, (8, 12), (1e-12, 3e-11), "control"),
             ("q", .02, 4, (12, 18), (1e-12, 3e-11), "control"),
             ("tight", .02, 4, (8, 12), (3e-13, 3e-12), "control")]
    manifest = []
    for name, dt, n, q, tolerance, role in cases:
        xml = out / (name + ".xml")
        sessions.generate(xml, table, out / name, nx=n, nr=n, dt=dt, steps=round(.4/dt),
                          backend="FastFull3DV1", mms_variant="PolynomialScalarsV1",
                          quadrature_x=q[0], quadrature_r=q[1], linear=tolerance[0], nonlinear=tolerance[1])
        edit(xml, {"WJReferencePolicy": POLICY, "WJImplicitSolver": "CoupledJFNKV1",
                   "WJStepPreconditioner": "CurrentDataDriftV1"}, {"WJStepNewtonSeconds": 45.})
        manifest.append(dict(label=name, role=role, session=identity(xml),
                             native=str(out / (name + ".result.json"))))
    atomic(out / "study.json", dict(schema="WJH27MMSStudyManifestV1", terminal_time=.4, runs=manifest,
        comparisons=[dict(label=kind, kind=kind, coarse="dt02", fine=name)
                     for kind, name in (("mesh", "h"), ("quadrature", "q"), ("tolerance", "tight"))]))


def assess_coupling(table, folder, out):
    from types import SimpleNamespace
    import coupled_return_check as coupled_return
    from implicit_algorithm_check import coefficient_delta

    limit = time.monotonic() + 155
    frozen = read(table)
    ref = po.load_reference(table)
    compare_data = lambda name: read(folder / (name + ".result.json"))["state"]
    old, baseline = compare_data("old-reference"), compare_data("old-reference-installed")
    require(old["reference_identity"] == baseline["reference_identity"] == po.TABLE_SHA,
            "historical parity requires the inherited table")
    def difference(a, b):
        require(a["sizes"] == b["sizes"], "cache or parity coefficient inventory differs")
        sizes = a["sizes"]
        return coefficient_delta(core.state(a["coefficients"], sizes), core.state(b["coefficients"], sizes), sizes, ref)

    parity = difference(old, baseline)
    require(all(value["scaled_MAX"] <= 1e-10 for value in parity.values()), "historical preset parity failed")

    fixture = read(folder / "fixtures.json")
    probe = folder / "cache"
    data = {name: read(probe / (name + ".result.json")) for name in coupled_return.ACCEPTED + ("failed-a",)}
    A, after = data["warm-a-first"]["state"], data["warm-a-return"]["state"]
    same_source = {row["state"]["source_identity"] for row in data.values()}
    require(len(same_source) == 1 and all(row["state"]["reference_identity"] == ref.identity for row in data.values()),
            "cache returns do not share one native reference and source")
    require(all(data[name]["accepted"] for name in coupled_return.ACCEPTED)
            and not data["failed-a"]["accepted"] and data["failed-a"]["history_unchanged"],
            "replay/failed-step disposition differs")
    require((probe / "failed-a.input.bin").read_bytes() == (probe / "failed-a.final.bin").read_bytes(),
            "failed accepted history bytes changed")
    require(A["input_history"] == after["input_history"], "A return history changed")
    require("restored explicit history differs" in read(probe / "stale-rejection.json")["reason"],
            "stale action guard was not exercised")
    comparisons = {}
    for label, base in (("warm-a-return", "warm-a-first"), ("cold-a", "warm-a-first"),
                        ("cold-b", "warm-b"), ("retry-a", "warm-a-first")):
        a, b = data[label]["state"], data[base]["state"]
        comparisons[label] = difference(a, b)
        require(all(v["scaled_MAX"] <= 1e-8 for v in comparisons[label].values()), "failed cache/cold parity: "+label)
    fixed = {}
    for label in ("warm-a-first", "warm-b", "warm-a-return"):
        packet = data[label]
        d = packet["state"]
        h = d["input_history"]
        state = core.state(h["current"], d["sizes"])
        fixed[label] = coupled_return.fixed_time_equations(d, frozen, state, h["time"], fixture["crossing_face"], limit)
        exact = fixed[label]["assemblies"][-1]["residuals"]
        sequence = ("momentum", "continuity", "gauge", "k", "epsilon", "trace")
        differences = {key: abs(actual - exact[key]) for key, actual in
                       zip(sequence, packet["exact_input_current_audit"]["residuals"])}
        require(max(differences.values()) <= 1e-10, f"native versus independent fixed current actions {differences}")
        fixed[label]["native_current_residual_difference"] = differences
    roots = [fixed[key]["assemblies"][-1]["roots"][0]["x"] for key in ("warm-a-first", "warm-b", "warm-a-return")]
    require(roots[0] > .5 > roots[1] and roots[2] == roots[0], "A/B/A did not cross/recover original owner pair")
    current = {}
    for label in ("warm-a-return", "warm-b"):
        d = data[label]["state"]
        measured = {}
        core.audit(d, frozen, SimpleNamespace(nx=20, nr=24, refine_nx=24, refine_nr=32, direct=False),
                   measured, limit)
        gates = measured["negative_controls"]
        require(measured["current_coupled_step_gate_pass"] and measured["nonzero_prescribed_body_coverage_pass"]
                and gates["frozen_prior_scalar_and_source"]["status"] == "DETECTED",
                "current-state, stale cutoff/scalar, or sampled-force guard failed")
        current[label] = {key: measured[key] for key in ("independent_assemblies", "negative_controls",
                             "native_prescribed_body_action_comparison", "current_coupled_step_gate_pass")}
    atomic(out, dict(status="PASS", source_identity=A["source_identity"], table=identity(table),
                     historical_parity=parity, cache_return_comparisons=comparisons,
                     failed_history_byte_preserved=True, roots=roots, fixed=fixed, current=current))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("prepare", "reference", "projection", "stale", "continuum", "grade", "cache", "temporal", "direct"))
    parser.add_argument("paths", type=Path, nargs="+")
    args = parser.parse_args()
    stages = dict(prepare=prepare, reference=reference, projection=projection,
                  stale=stale_histories, continuum=continuum_actions, grade=grade_projection, cache=assess_coupling, temporal=temporal_mesh_control, direct=direct_fixture)
    stages[args.stage](*args.paths)


if __name__ == "__main__":
    main()
