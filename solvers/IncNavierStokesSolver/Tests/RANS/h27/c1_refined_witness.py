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
from c1_paired_grid import subdivisions
from c1_rejected_points import neighborhood, point
from c1_support_controls import coincides, contains, locations
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
    xs, zs, parent_geometry = mesh_axes(args.template)
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
    splits, mesh_axis, radius_axis = dict(), None, None
    if args.bisect_radial_at is not None and (args.streamwise or args.radial or args.axial_envelope):
        raise ValueError("one witnessed radial bisection cannot be combined with other subdivisions")
    if args.axial_envelope:
        if args.streamwise or args.radial:
            raise ValueError("a frozen axial envelope already defines all subdivisions")
        endpoints, midpoints = subdivisions(xs, args.axial_envelope)
        mesh_axis = sorted(xs+endpoints+midpoints)
        radius_axis = zs
        splits = dict(axial_envelope=args.axial_envelope, inserted_endpoints=endpoints,
                      bisected_midpoints=midpoints, parent_axial=xs, child_axial=mesh_axis,
                      cells=(len(mesh_axis)-1)*(len(zs)-1))
    else:
        mesh_axis = insert(xs, args.streamwise)
        radius_axis = insert(zs, [r+d0 for r in args.radial])
    if args.bisect_radial_at is not None:
        z = d0 + args.bisect_radial_at
        owners = [(lo,hi) for lo,hi in zip(zs,zs[1:]) if lo+1e-12 < z < hi-1e-12]
        if len(owners) != 1:
            raise ValueError("the radial witness must identify one original interval interior")
        (lo,hi), = owners
        middle = lo + (hi-lo)/2
        radius_axis = insert(zs,[middle])
        splits = dict(radial_witness=args.bisect_radial_at,parent_vertex_z_bounds=[lo,hi],
                      parent_radial_bounds=[lo-d0,hi-d0],inserted_vertex_z=middle,
                      inserted_radial_midpoint=middle-d0,
                      parent_axial=xs,child_axial=mesh_axis,
                      cells=(len(mesh_axis)-1)*(len(radius_axis)-1))
    mesh(geometry, mesh_axis, radius_axis, quadrature=12)
    if args.axial_envelope or args.bisect_radial_at is not None:
        # Preserve original tokens on both axes, including at a new cut.
        vertices = {}
        axial, radial = {}, {}
        for v in parent_geometry.findall("VERTEX/V"):
            a,z,w = v.text.split()
            if float(w)!=0.:
                raise ValueError("expected a planar chart vertex")
            vertices[float(a),float(z)] = v.text
            if axial.setdefault(float(a),a)!=a or radial.setdefault(float(z),z)!=z:
                raise ValueError("chart coordinate uses multiple original tokens")
        gen = ET.parse(geometry)
        for v in gen.findall("GEOMETRY/VERTEX/V"):
            a,z,w = v.text.split()
            x, height = float(a),float(z)
            v.text = vertices.get((x,height),f"{axial.get(x,a)} {radial.get(height,z)} {w}")
        gen.write(geometry,encoding="unicode")
        splits["radial_xml_tokens"] = [radial[z] for z in zs]
        recorded_x,recorded_z,recorded = mesh_axes(geometry)
        if recorded_x!=mesh_axis or recorded_z!=radius_axis or any(
                v.text!=vertices[key] for v in recorded.findall("VERTEX/V")
                if (key:=tuple(map(float,v.text.split()[:2]))) in vertices):
            raise ValueError("mesh output does not preserve the parent's chart bytes and intended cuts")
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
               inserted_streamwise=args.streamwise, inserted_radius=args.radial, **splits))


def inspect(args):
    deadline = time.monotonic()+args.seconds
    ref = po.load_reference(args.table)
    data, _ = load_c1_state(args.native, True)
    audit = Audit(data, ref, args.order, deadline)
    controls = []
    coordinates = [(None, a, r) for a, r in args.point]
    if args.release_controls:
        if not coincides(data["geometry"]["h"], ref.h, ref.H):
            raise ValueError("support coordinate differs from the pinned reference")
        coordinates += [(label, a, r) for label, (a, r) in locations(data["geometry"]).items()]
    for label, a, r in coordinates:
        for cell in audit.mesh:
            if contains(cell, a, r, ref.H):
                row = point(audit, cell, a, r, "fixed-comparison")
                if label is not None:
                    row["control"] = label
                controls.append(row)
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
                         source_identity=data["source_identity"],
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
    nested.add_argument("--bisect-radial-at", type=float,
                        help="bisect the parent radial interval containing this chart radius")
    nested.add_argument("--axial-envelope",nargs=2,type=float,
                        help="bisect each interval in this fixed envelope, preserving radial XML coordinates")
    nested.set_defaults(run=prepare)
    compare = modes.add_parser("inspect")
    rejection = modes.add_parser("rejected-input")
    for mode in (compare, rejection):
        for key in ("native", "table", "out"):
            mode.add_argument("--"+key, type=Path, required=True)
    compare.add_argument("--point", nargs=2, action="append", type=float, default=[])
    compare.add_argument("--release-controls", action="store_true",
                         help="check required chart locations and every containing owner")
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
