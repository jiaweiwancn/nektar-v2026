#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Prepare and independently compare the fixed-height C1 resolution witnesses."""

import argparse
import json
from pathlib import Path
import time
import xml.etree.ElementTree as ET

import numpy as np

from c1_oracle import Audit, load_c1_state
from c1_rejected_points import neighborhood, point
from coupled_fixtures import generate, mesh
import physical_oracle as po
from physical_transient_oracle import state_hash
from resources import identity
from task2_outlet_baseline import mesh_axes
from task2r_verify import read


def write(path, payload):
    with path.open("x") as stream:
        json.dump(payload, stream, indent=2, allow_nan=False)
        stream.write("\n")


def prepare(args):
    session = ET.parse(args.template)
    xs, zs, _ = mesh_axes(args.template)
    info = session.find("CONDITIONS/SOLVERINFO")
    parameters = session.find("CONDITIONS/PARAMETERS")
    required = {item.get("PROPERTY"): item.get("VALUE") for item in info}
    values = dict(item.text.replace(" ", "").split("=", 1) for item in parameters)
    if required["WJNativeMode"] != "Physical" or required["WJGeometry"] != "NormalCosineHillV1":
        raise ValueError("nested witness requires an original physical normal-hill template")
    d0 = float(values["WJPhysicalDisplacement"])

    def insert(edges, points):
        if any(not edges[0] < v < edges[-1] or min(abs(v-e) for e in edges) < 1e-12
               for v in points):
            raise ValueError("new cuts must be strictly inside existing chart cells")
        return sorted(set(edges + points))

    args.out.mkdir(parents=True)
    geometry = args.out / "geometry.xml"
    boundary = args.out / "boundary-template.xml"
    table = Path(required["WJReferenceTable"])
    mesh(geometry, insert(xs, args.streamwise), insert(zs, [r+d0 for r in args.radial]),
         quadrature=12)
    generate(boundary, table, args.out / "unused", template=geometry,
             mode="Physical", outlet_policy=required["WJOutletPolicy"])
    rebuilt = ET.parse(boundary)
    for key in ("GEOMETRY", "EXPANSIONS", "CONDITIONS/BOUNDARYREGIONS",
                "CONDITIONS/BOUNDARYCONDITIONS"):
        parent = session.find("CONDITIONS") if "/" in key else session.getroot()
        current = session.find(key)
        index = list(parent).index(current)
        parent.remove(current)
        parent.insert(index, rebuilt.find(key))
    target = args.out / "run.xml"
    session.write(target, encoding="unicode")
    write(args.out / "mesh-inputs.json",
          dict(template=identity(args.template), session=identity(target),
               geometry=identity(geometry), table=identity(table),
               inserted_streamwise=args.streamwise, inserted_radius=args.radial))


def inspect(args):
    deadline = time.monotonic()+args.seconds
    ref = po.load_reference(args.table)
    data, _ = load_c1_state(args.native, True)
    audit = Audit(data, ref, args.order, deadline)
    controls = []
    for a, r in args.point:
        for cell in audit.mesh:
            xl, xr, rl, rr = cell["bounds"]
            if xl <= a <= xr and rl <= r <= rr:
                controls.append(point(audit, cell, a, r, "fixed-comparison"))
    local = []
    for label, filename in args.minima:
        report = read(filename)
        if report["native"]["sha256"] != identity(args.native)["sha256"]:
            raise ValueError("principal witness does not bind the inspected coefficients")
        witness = report["regions"][label]["witness"]
        owner, = (c for c in audit.mesh if c["geometry"] == witness["geometry"])
        a, r = witness["a"], witness["rho"]
        sides, dense, count = neighborhood(audit, owner, a, r, args.samples)
        _, denser, _ = neighborhood(audit, owner, a, r, 2*args.samples-1)
        checks = [(owner, a, r, "reported-minimum"), *sides, dense, denser]
        local.append(dict(region=label, receipt=identity(filename), nearby_owners=count,
                          grid_orders=[args.samples, 2*args.samples-1],
                          values=[point(audit, *sample) for sample in checks]))
    force = np.zeros(2)
    tangential = 0.
    for cell in audit.mesh:
        if cell["bounds"][2] != 0.:
            continue
        audit.clock()
        a, r, weight = audit.face(cell, 2)
        length, normal = audit.hill.face(cell, 2, a, r)
        metric = audit.hill.metrics(cell, a, r)
        _, _, _, pressure, stress, *_ = audit.sample(cell, a, r, metric)
        traction = np.einsum("qij,qj->qi", stress-pressure[:, None, None]*np.eye(2), normal)
        weight *= length
        force += weight @ traction
        tangential += float(weight @ np.einsum("qi,qi->q", traction, metric["t"]))
    write(args.out, dict(native=identity(args.native), table=identity(args.table),
                         coefficient_sha256=data["state_sha256"], order=args.order,
                         wall_force_on_fluid=force.tolist(), tangential_wall_force_on_fluid=tangential,
                         fixed_controls=controls, local_minima=local))


def rejected_input(args):
    # Deliberately move one free matching vertex outside the certified family.
    # The packet remains a typed diagnostic input and retains its exact bytes.
    packet = read(args.native)
    state = packet["state"]
    ref = po.load_reference(args.table)
    eligible = set(state["k_free_ids"])
    for element in state["elements"]:
        if abs(element["bounds"][3]-ref.H) > 1e-12:
            continue
        vertex = next((j for j in po.EDGES[3][:2] if element["k"][j] in eligible), None)
        if vertex is not None:
            state["coefficients"]["k"][element["k"][vertex]] = (
                element["ks"][vertex]*1.12*ref.Kin)
            break
    else:
        raise ValueError("missing free matching vertex")
    state["state_sha256"] = state_hash(state["coefficients"])
    write(args.out, packet)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest="mode", required=True)
    nested = modes.add_parser("nested")
    nested.add_argument("template", type=Path)
    nested.add_argument("out", type=Path)
    nested.add_argument("--streamwise", nargs="+", type=float, default=[])
    nested.add_argument("--radial", nargs="+", type=float, default=[])
    nested.set_defaults(run=prepare)
    compare = modes.add_parser("inspect")
    rejection = modes.add_parser("rejected-input")
    for mode in (compare, rejection):
        for key in ("native", "table", "out"):
            mode.add_argument("--"+key, type=Path, required=True)
    compare.add_argument("--point", nargs=2, action="append", type=float, default=[])
    compare.add_argument("--minima", nargs=2, action="append", default=[])
    compare.add_argument("--order", type=int, default=32)
    compare.add_argument("--samples", type=int, default=25)
    compare.add_argument("--seconds", type=float, default=180.)
    compare.set_defaults(run=inspect)
    rejection.set_defaults(run=rejected_input)
    args = parser.parse_args()
    args.run(args)


if __name__ == "__main__":
    main()
