#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Prepare and time matched conventional VCSWeakPressure controls on Bf."""

import argparse
import math
from pathlib import Path
import statistics
import xml.etree.ElementTree as ET

import numpy as np
from numpy.polynomial.chebyshev import cheb2poly

import physical_oracle as po
from physical_cost_screen import STEP
from resources import atomic, identity
from task2r_advancing import read, require, verified
from task2r_run import GIB
from task2r_verify import edit


def inlet_expression(ref):
    radius = f"(y-{ref.d0:.17g})"
    pieces = []
    for segment in ref.inlet_panels:
        lo, hi = segment["domain"]
        t = f"max(-1,min(1,((2*{radius}-{lo:.17g}-{hi:.17g})/{(hi-lo):.17g})))"
        polynomial = "0"
        for coefficient in reversed(cheb2poly(segment["U"])):
            polynomial = f"(({coefficient:.17g})+({t})*({polynomial}))"
        weight = f"(sign({radius}-({lo:.17g}))-sign({radius}-({hi:.17g})))*0.5"
        pieces.append(f"({weight})*({polynomial})")
    outer = f"{ref.utau/ref.kappa:.17g}*log(max({radius},{ref.H:.17g})/{ref.z0:.17g})"
    pieces.append(f"(1+sign({radius}-{ref.H:.17g}))*0.5*({outer})")
    return "+".join(pieces)


def prepare(root, label, coupled_path, cap, rss_gib, svv):
    coupled = read(coupled_path)
    require(coupled["role"] == "startup" and coupled["phase"] == "r3-advance",
            "stock control must match a retained-inlet coupled run")
    bf_path = verified(coupled["bf_binding"])
    bound = read(bf_path)
    fixture = verified(coupled["fixture"])
    option = verified(coupled["options"])
    app = verified(coupled["application"])
    table = verified(coupled["table"])
    geometry = verified(coupled["geometry"])
    ref = po.load_reference(table)
    out = root / "science-r3" / label
    require(label.startswith("r3-stock-") and svv in ("Off", "ExpKernel"), "invalid stock label or SVV control")
    out.mkdir()
    session = out / "run.xml"
    doc = ET.parse(fixture)
    conditions = doc.find("CONDITIONS")
    si = conditions.find("SOLVERINFO")
    for e in list(si):
        if e.get("PROPERTY").startswith("WJ"):
            si.remove(e)
    params = conditions.find("PARAMETERS")
    for e in list(params):
        key = e.text.split("=",1)[0].strip()
        if key.startswith("WJ") or key == "IO_CheckSteps":
            params.remove(e)
    expansions = doc.find("EXPANSIONS")
    for e in expansions:
        if e.get("FIELDS") != "p":
            e.set("FIELDS", "u,v")
        e.set("NUMPOINTS", "6,8")
        e.set("POINTSTYPE", "GaussGaussLegendre,GaussGaussLegendre")
    for group in ("VARIABLES", "FUNCTION"):
        block = conditions.find(f"{group}" if group == "VARIABLES" else "FUNCTION[@NAME='InitialConditions']")
        for e in list(block):
            if e.text in ("k","epsilon") or e.get("VAR") in ("k","epsilon"):
                block.remove(e)
    for i, e in enumerate(conditions.find("VARIABLES")):
        e.set("ID", str(i))
    startup = inlet_expression(ref)
    for region in conditions.findall("BOUNDARYCONDITIONS/REGION"):
        for e in list(region):
            if e.get("VAR") in ("k","epsilon"):
                region.remove(e)
        side = int(region.get("REF"))
        pressure = region.find("N[@VAR='p']")
        if side == 1:
            pressure.tag = "D"
        else:
            pressure.set("USERDEFINEDTYPE","H")
        if side == 0:
            region.find("D[@VAR='u']").set("VALUE",startup)
    conditions.find("FUNCTION[@NAME='InitialConditions']/E[@VAR='u']").set("VALUE",startup)
    frequency = coupled["recovery_stride"]
    if frequency:
        filters = doc.find("FILTERS")
        if filters is None:
            filters = ET.SubElement(doc.getroot(),"FILTERS")
        check = ET.SubElement(filters,"FILTER",TYPE="Checkpoint")
        ET.SubElement(check,"PARAM",NAME="OutputFrequency").text = str(frequency)
    doc.write(session, encoding="unicode")
    numerical = dict(CheckAbortSteps=1)
    if svv == "ExpKernel":
        numerical.update(SVVDiffCoeff=.01, SVVCutoffRatio=.8)
    stock_switch = "False" if svv == "Off" else svv
    edit(session,dict(SolverType="VCSWeakPressure",SpectralVanishingViscosity=stock_switch),numerical)
    opt = out / "run.xml.opt"
    opt.write_bytes(option.read_bytes())
    command = ["mpirun", "-np", "1", str(app), "--use-opt-file", str(opt), str(session)]
    atomic(out / "admission.json", dict(schema="WJH27R3AdmissionV1", label=label, phase="r3-stock",
        role="conventional_retained_inlet", seconds=cap, peak_new_bytes=128*1024**2, rss_limit_bytes=rss_gib*GIB,
        cpu=bound["affinity"]["cpu"], bf_binding=identity(bf_path), command=command, geometry=identity(geometry),
        fixture=identity(session), options=identity(opt), application=identity(app),
        timestep=coupled["timestep"], requested_steps=coupled["requested_steps"], stock_svv=svv,
        recovery_stride=frequency, inputs=[identity(p) for p in (coupled_path,geometry,table)],
        physics=dict(viscosity=ref.nu, order=[3,2], gauss=[6,8],
            inlet="Exact-inner polynomial panels kinematically continued to the retained outer log profile",
            transverse="Zero at inlet, ground, top and exit",
            wall_loads="No-slip ground, zero diffusive axial trace at top and outlet",
            pressure="Zero outlet Dirichlet and weak high-order Neumann data elsewhere",
            excluded="RANS transport, current-K wall cutoff and stress/production/AWF panels",
            output="Native stock velocity/pressure fields and on-cadence checkpoint filter")))


def report(root, label, out):
    case = root / "science-r3" / label
    admission = read(case / "admission.json")
    operation = root / "operations" / label / "receipt.json"
    data = read(operation)
    elapsed = data["result"]
    log = verified(elapsed["log"]).read_text()
    steps = [(int(n),float(t),float(s)) for n,t,s in STEP.findall(log)]
    count, dt = admission["requested_steps"], admission["timestep"]
    files = sorted(case.glob("*.fld")) + sorted(case.glob("*.chk"))
    stable = (elapsed["status"] == "PASS" and len(steps) == count
              and all(n == i and abs(t-i*dt) < 1e-6 and math.isfinite(s) and s > 0
                      for i,(n,t,s) in enumerate(steps,1))
              and any(p.suffix == ".fld" for p in files))
    atomic(out, dict(schema="WJH27R3ConventionalCostV1", status="MET" if stable else "NOT_MET",
        admission=identity(case / "admission.json"), operation=identity(operation), solver_status=elapsed["status"],
        process_seconds=elapsed.get("child_elapsed_seconds"), step_report=steps,
        warmed_mean_seconds=statistics.mean(s for _,_,s in steps[1:]) if stable and count>1 else None,
        outputs=[identity(p) for p in files if p.is_file()]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root",type=Path)
    sub=parser.add_subparsers(dest="action",required=True)
    prepare_parser=sub.add_parser("prepare")
    prepare_parser.add_argument("label")
    prepare_parser.add_argument("coupled_path",type=Path)
    prepare_parser.add_argument("--cap",type=int,required=True)
    prepare_parser.add_argument("--rss-gib",type=int,default=8)
    prepare_parser.add_argument("--svv",choices=("Off","ExpKernel"),default="Off")
    report_parser=sub.add_parser("report")
    report_parser.add_argument("label")
    report_parser.add_argument("out",type=Path)
    args=vars(parser.parse_args())
    action=args.pop("action")
    (prepare if action=="prepare" else report)(**args)


if __name__=="__main__":
    main()
