"""Reproducible native-IMEX M4 sessions; solver execution belongs on Linux.

The neutral default fixes H=2, Z=30 and the physical horizontal dimensions.
In 2D the native mesh coordinates (x,y) represent physical (x,z). Manufactured
sources are the independent strong conservative expressions from M3.
"""
import sys

sys.dont_write_bytecode = True

import math
from pathlib import Path
import xml.etree.ElementTree as ET

from cls_epawf_m3_cases import (
    evolving_session as _m3_evolving_session,
    set_info,
    set_parameter,
    write,
)


def _knots(values, name):
    values = [float(value) for value in values]
    if (len(values) < 2 or not all(math.isfinite(value) for value in values)
            or any(b <= a for a, b in zip(values, values[1:]))):
        raise ValueError(f"{name} must contain finite, strictly increasing knots")
    return values


def _remove_info(root, name):
    info = root.find("CONDITIONS/SOLVERINFO")
    for node in list(info):
        if node.get("PROPERTY", "").lower() == name.lower():
            info.remove(node)


def tensor_mesh(root, zknots, degree, dimension, length, span, points=None,
                xknots=None, yknots=None):
    """Replace topology and expansions by a conforming affine tensor mesh.

    Shared edges/faces have one native ID. Hex faces use native order bottom,
    front, right, back, left, top; the corresponding local vertices have x as
    the fastest coordinate. Velocity/k/epsilon have degree p and pressure p-2.
    Explicit x/y knots, when supplied, determine the horizontal mesh bounds.
    """
    if dimension not in (2, 3):
        raise ValueError("dimension must be 2 or 3")
    if int(degree) != degree or degree < 3:
        raise ValueError("degree must be an integer >= 3")
    degree = int(degree)
    z = _knots(zknots, "zknots")
    x = _knots([0., length] if xknots is None else xknots, "xknots")
    y = _knots([0., span] if yknots is None else yknots, "yknots")
    if dimension == 2 and yknots is not None:
        raise ValueError("yknots are only meaningful for a 3D mesh")
    if points is None:
        points = [max(8, degree+3)]*dimension
    elif isinstance(points, int):
        points = [points]*dimension
    else:
        points = list(points)
    if (len(points) != dimension
            or any(int(value) != value or value < degree+1 for value in points)):
        raise ValueError("points must provide at least p+1 points per axis")

    conditions = root.find("CONDITIONS")
    for tag in ("BOUNDARYREGIONS", "BOUNDARYCONDITIONS"):
        for node in list(conditions.findall(tag)):
            conditions.remove(node)
    variables = conditions.find("VARIABLES")
    if variables is None:
        variables = ET.SubElement(conditions, "VARIABLES")
    variables.clear()
    names = ["u", "v"] + (["w"] if dimension == 3 else [])
    names += ["k", "eps_tilde", "p"]
    for index, name in enumerate(names):
        ET.SubElement(variables, "V", {"ID": str(index)}).text = name
    expansions = root.find("EXPANSIONS")
    if expansions is None:
        expansions = ET.SubElement(root, "EXPANSIONS")
    expansions.clear()
    for fields, p in ((names[:-1], degree), (["p"], degree-2)):
        ET.SubElement(expansions, "E", {
            "COMPOSITE": "C[0]",
            "NUMMODES": ",".join([str(p+1)]*dimension),
            "BASISTYPE": ",".join(["Modified_A"]*dimension),
            "POINTSTYPE": ",".join(["GaussLobattoLegendre"]*dimension),
            "NUMPOINTS": ",".join(str(int(value)) for value in points),
            "FIELDS": ",".join(fields),
        })
    geometry = root.find("GEOMETRY")
    if geometry is None:
        geometry = ET.SubElement(root, "GEOMETRY")
    geometry.clear()
    geometry.attrib.update({"DIM": str(dimension), "SPACE": str(dimension)})
    vertices = ET.SubElement(geometry, "VERTEX")
    edges = ET.SubElement(geometry, "EDGE")
    faces = ET.SubElement(geometry, "FACE") if dimension == 3 else None
    elements = ET.SubElement(geometry, "ELEMENT")
    edge_ids, face_ids = {}, {}

    def edge(a, b):
        key = tuple(sorted((a, b)))
        if key not in edge_ids:
            edge_ids[key] = len(edge_ids)
            ET.SubElement(edges, "E", {"ID": str(edge_ids[key])}).text = f"{a} {b}"
        return edge_ids[key]

    def face(v):
        key = tuple(sorted(v))
        if key not in face_ids:
            face_ids[key] = len(face_ids)
            ids = [edge(a, b) for a, b in zip(v, v[1:]+v[:1])]
            ET.SubElement(faces, "Q", {"ID": str(face_ids[key])}).text = " ".join(map(str, ids))
        return face_ids[key]

    nx, ny = len(x), len(y)
    if dimension == 2:
        for k, height in enumerate(z):
            for i, tangent in enumerate(x):
                ET.SubElement(vertices, "V", {"ID": str(k*nx+i)}).text = f"{tangent:.17g} {height:.17g} 0"
        for k in range(len(z)-1):
            for i in range(nx-1):
                v = [k*nx+i, k*nx+i+1, (k+1)*nx+i+1, (k+1)*nx+i]
                ids = [edge(a, b) for a, b in zip(v, v[1:]+v[:1])]
                ET.SubElement(elements, "Q", {"ID": str(len(elements))}).text = " ".join(map(str, ids))
    else:
        def vertex(i, j, k):
            return (k*ny+j)*nx+i
        for k, height in enumerate(z):
            for j, cross in enumerate(y):
                for i, tangent in enumerate(x):
                    ET.SubElement(vertices, "V", {"ID": str(vertex(i, j, k))}).text = f"{tangent:.17g} {cross:.17g} {height:.17g}"
        local_faces = ((0, 1, 2, 3), (0, 1, 5, 4), (1, 2, 6, 5),
                       (2, 3, 7, 6), (3, 0, 4, 7), (4, 5, 6, 7))
        for k in range(len(z)-1):
            for j in range(ny-1):
                for i in range(nx-1):
                    v = [vertex(i, j, k), vertex(i+1, j, k),
                         vertex(i+1, j+1, k), vertex(i, j+1, k),
                         vertex(i, j, k+1), vertex(i+1, j, k+1),
                         vertex(i+1, j+1, k+1), vertex(i, j+1, k+1)]
                    ids = [face([v[index] for index in f]) for f in local_faces]
                    ET.SubElement(elements, "H", {"ID": str(len(elements))}).text = " ".join(map(str, ids))
    composite = ET.SubElement(geometry, "COMPOSITE")
    shape = "Q" if dimension == 2 else "H"
    interval = "0" if len(elements) == 1 else f"0-{len(elements)-1}"
    ET.SubElement(composite, "C", {"ID": "0"}).text = f"{shape}[{interval}]"
    ET.SubElement(geometry, "DOMAIN").text = "C[0]"


def _evolution(root, steps, order, dt, time, condensed, restart, checkpoint):
    if int(steps) != steps or steps < 0:
        raise ValueError("steps must be a nonnegative integer")
    if order not in (1, 2):
        raise ValueError("native IMEX order must be 1 or 2")
    if not math.isfinite(dt) or dt <= 0 or not math.isfinite(time):
        raise ValueError("dt must be positive and dt/time must be finite")
    set_info(root, "SolverType", "RANSEquilibriumWallEvolution")
    set_info(root, "RANSWallCondense", "True" if condensed else "False")
    set_info(root, "RANSWallRestart", "None" if restart is None else str(restart))
    set_info(root, "RANSWallCheckpoint", "None" if checkpoint is None else str(checkpoint))
    _remove_info(root, "RANSWallStageMode")
    set_parameter(root, "NumSteps", int(steps))
    set_parameter(root, "TimeStep", repr(float(dt)))
    set_parameter(root, "Time", repr(float(time)))
    conditions = root.find("CONDITIONS")
    for node in list(conditions.findall("TIMEINTEGRATIONSCHEME")):
        conditions.remove(node)
    scheme = ET.SubElement(conditions, "TIMEINTEGRATIONSCHEME")
    ET.SubElement(scheme, "METHOD").text = "IMEX"
    ET.SubElement(scheme, "ORDER").text = str(order)


def neutral_session(template, destination, degree=3, inner_knots=None,
                    outer_knots=None, dimension=2, steps=0, order=2, dt=.05,
                    eta=0, chi=1, direction=(1., 0.), condensed=False,
                    length=300., span=20., tangential_order=4,
                    normal_order=16, restart=None, checkpoint=None,
                    xknots=None, yknots=None):
    """Generate the isolated flat neutral case without analytic FUNCTIONs.

    Knot overrides permit h studies while preserving their common matching
    height. Restart sessions retain Time=0: native checkpoint time/history are
    restored by the solver, and the physical signature remains unchanged.
    """
    inner = _knots([0., .1, 2.] if inner_knots is None else inner_knots, "inner_knots")
    outer = _knots([2., 4., 8., 16., 30.] if outer_knots is None else outer_knots, "outer_knots")
    if inner[0] != 0. or inner[-1] != outer[0]:
        raise ValueError("inner knots must start at zero and meet the outer knots")
    if len(direction) != 2 or not all(math.isfinite(value) for value in direction):
        raise ValueError("direction must contain two finite horizontal components")
    if math.hypot(*direction) == 0.:
        raise ValueError("direction must be nonzero")
    if dimension == 2 and direction[1] != 0.:
        raise ValueError("a 2D mesh has only the physical x horizontal direction")
    if eta not in (0, 1) or chi not in (0, 1):
        raise ValueError("eta and chi must be 0 or 1")
    root = ET.parse(template).getroot()
    tensor_mesh(root, inner+outer[1:], degree, dimension, length, span,
                xknots=xknots, yknots=yknots)
    conditions = root.find("CONDITIONS")
    for node in list(conditions.findall("FUNCTION")):
        conditions.remove(node)
    _evolution(root, steps, order, dt, 0., condensed, restart, checkpoint)
    set_info(root, "RANSWallReference", "Neutral")
    for name, value in (("RANSWallHeight", inner[-1]), ("RANSWallEta", eta),
                        ("RANSWallChi", chi), ("RANSWallDirectionX", direction[0]),
                        ("RANSWallDirectionY", direction[1]),
                        ("RANSWallTangentialOrder", tangential_order),
                        ("RANSWallNormalOrder", normal_order)):
        set_parameter(root, name, value)
    return write(root, destination)


def evolving_session(template, destination, degree=6, time=.3, steps=4,
                     order=2, dt=.05, condensed=False, restart=None,
                     checkpoint=None, dimension=2):
    """Generate the unperturbed time-dependent 2D strong manufactured case.

    Four radial elements admit four MPI ranks. M3's native StageHistory
    functions are removed because native IMEX supplies the complete history.
    Keep the same initial ``time`` when generating a restart session.
    """
    if dimension != 2:
        raise ValueError("the strong manufactured evolution case is defined only in 2D")
    destination = Path(destination)
    _m3_evolving_session(template, destination, degree=degree, time=time,
                         condensed=condensed, perturb=False)
    root = ET.parse(destination).getroot()
    tensor_mesh(root, [0., .5, 2., 2.02, 4.], degree, 2, 1., 1.)
    conditions = root.find("CONDITIONS")
    for node in list(conditions.findall("FUNCTION")):
        if node.get("NAME", "").startswith("StageHistory"):
            conditions.remove(node)
    _evolution(root, steps, order, dt, time, condensed, restart, checkpoint)
    set_info(root, "RANSWallReference", "Functions")
    set_parameter(root, "RANSWallHeight", 2)
    set_parameter(root, "RANSWallEta", 0)
    set_parameter(root, "RANSWallChi", 1)
    set_parameter(root, "RANSWallNormalOrder", 16)
    return write(root, destination)
