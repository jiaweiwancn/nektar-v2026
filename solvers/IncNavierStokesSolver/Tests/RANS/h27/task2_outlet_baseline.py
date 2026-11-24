#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Verify and freeze the two archived steady states and their native inputs."""

import hashlib
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

from resources import atomic, identity
from task2r_run import memory_snapshot, require, workers


def read(path):
    return json.loads(Path(path).read_text())


def verified(item):
    require(identity(item["path"]) == item, f"archived content differs: {item['path']}")
    return Path(item["path"])


def source_identity(root):
    eq = root / "solvers/IncNavierStokesSolver/EquationSystems"
    turb = root / "library/SolverUtils/Turbulence"
    names = ["WJEARSM", "RANSModels", "AnalyticalWallFunction"]
    recipe = eq.parent / "Tests/RANS/h27/source_identity.cmake"
    if recipe.is_file() and "RANSWallNormalMap.cpp" in recipe.read_text():
        names.append("RANSWallNormalMap")
    files = sorted([*eq.glob("WJNative*.h"), *eq.glob("WJNative*.cpp"),
                    *(turb / f"{name}.{ext}" for name in names for ext in ("h", "cpp"))])
    inventory = [identity(path) for path in files]
    material = "WJH27-v2|C++17|no-fast-math|fp-contract=off|" + "|".join(
        entry["sha256"] for entry in inventory)
    return hashlib.sha256(material.encode()).hexdigest(), inventory


def mesh_axes(path):
    geometry = ET.parse(path).getroot().find("GEOMETRY")
    vertices = [tuple(map(float, v.text.split())) for v in geometry.find("VERTEX")]
    xs, zs = (sorted({v[i] for v in vertices}) for i in (0, 1))
    require(len(vertices) == len(xs) * len(zs) and
            len(geometry.find("ELEMENT")) == (len(xs)-1) * (len(zs)-1),
            "expected a complete rectangular tensor mesh")
    return xs, zs, geometry


def main(root, summary_path, bundle_path, out):
    summary = read(summary_path)
    bundle = read(bundle_path)
    for item in bundle.values():
        verified(item)
    executable = verified(summary["application"])
    verified(summary["reference"])
    verified(summary["original_690_witness"])
    source, inventory = source_identity(root / "source-r2-v1")
    require(source == summary["source"] and source.encode() in executable.read_bytes(),
            "retained source does not bind the archived executable")
    result = dict(schema="WJH27OutletDiagnosticsBaselineV1",
                  git_baseline="cabca8e91a76ad3fa365692612dc16c9312fcaf4",
                  s2_summary=identity(summary_path), s2_bundle=bundle,
                  source_identity=source, source_files=inventory,
                  application=summary["application"], table=summary["reference"], cases={})
    axes = []
    for name, path in (("A0", "task2rs-s2-tip-resume"), ("A1", "task2rs-s2-nested-finish")):
        case = summary["cases"][path]
        inputs = read(verified(case["inputs"]))
        geometry, session = (verified(inputs[key]) for key in ("geometry", "session"))
        xs, zs, geom = mesh_axes(geometry)
        _, _, session_geom = mesh_axes(session)
        for tag in ("VERTEX", "EDGE", "ELEMENT"):
            require([(e.attrib, e.text.split()) for e in geom.find(tag)] ==
                    [(e.attrib, e.text.split()) for e in session_geom.find(tag)],
                    "native session differs from archive mesh")
        state = read(verified(case["native_state"]))["state"]
        verified(case["native_steady"])
        require(state["source_identity"] == source and
                state["reference_identity"] == summary["reference"]["sha256"] and
                state["step"] == state["time"] == 0 and not state["gauge"] and
                state["inner_k_transport"] == "ScalarCdV2" and
                [len(xs)-1, len(zs)-1] == case["cells"] and
                [state["quadratureX"], state["quadratureR"]] == [6, 8] and
                [state["linear_tolerance"], state["nonlinear_tolerance"]] == [1e-12, 1e-10] and
                max(state["residuals"]) <= 1e-10, "archive state misses frozen scientific binding")
        result["cases"][name] = dict(
            cells=case["cells"], input_manifest=case["inputs"], geometry=inputs["geometry"],
            session=inputs["session"], session_opt=identity(Path(str(session)+".opt")),
            state=case["native_state"], steady=case["native_steady"],
            coefficient_sha256=state["state_sha256"], signature=state["signature"],
            boundary_identity=state["boundary_identity"], residuals=state["residuals"],
            x_edges=xs, physical_z_edges=zs)
        axes.append((xs, zs))
    require(all(set(coarse) < set(fine) for coarse, fine in zip(*axes)),
            "archived M1 must be nested in both coordinates")
    xs, zs = axes[1]
    collar = [i for i, z in enumerate(zs) if abs(z-.1875) <= 1e-13]
    require(len(collar) == 1 and 0 < collar[0] < len(zs)-1, "missing r=.10 collar face")
    j = collar[0]
    midpoints = [(a+b)/2 for a, b in zip(zs[j-1:j+1], zs[j:j+2])]
    result["M1c"] = dict(x_edges=xs, physical_z_edges=sorted(zs+midpoints),
                         original_collar_owner_r=[z-.0875 for z in zs[j-1:j+2]],
                         inserted_radial_edges=[z-.0875 for z in midpoints])
    result["measurements"] = read(verified(summary["cases"]["task2rs-s2-nested-finish"]
                                           ["inputs"]))["observations"]
    result["host"] = dict(memory=memory_snapshot(), workers=workers())
    atomic(out, result)


if __name__ == "__main__":
    main(*map(Path, sys.argv[1:]))
