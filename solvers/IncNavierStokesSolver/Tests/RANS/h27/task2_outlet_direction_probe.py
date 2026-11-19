#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check overlap of the direct radial and streamwise responses at signed extrema."""

import argparse
from pathlib import Path
import time

import physical_oracle as po
from resources import atomic, identity
from task2_outlet_baseline import read
from task2_outlet_case import get_context
from task2_outlet_owner_probe import is_owner
from task2_outlet_qualification import load_case, pass_audit, raw_owner
from task2r_run import require


def common_owner(e, base, x, r):
    return (is_owner(e,x,r) and all(
        min(e["bounds"][d+1],base["bounds"][d+1])-max(e["bounds"][d],base["bounds"][d])>1e-12
        for d in (0,2)))


def quantity(raw, field, component):
    if field=="R":
        i,j = ("xyr".index(c) for c in component[1:])
        return raw["R"][i][j]
    return raw["Vxr" if field=="traction" else "Pi"]


def check(root, out, seconds):
    radial, streamwise = ("Br3","Br4"), ("Bx2","Bx3")
    context = get_context(root)
    ref = po.load_reference(context[2])
    deadline = time.monotonic()+seconds
    bases, bindings, comparisons = {}, {}, {}
    for name in ("Bf",*radial,*streamwise):
        audit,_ = pass_audit(root,name)
        bindings[name] = audit["provenance"]
    for name in (*radial,*streamwise):
        path = root/"evidence"/f"outlet-qual-{name}-Bf.json"
        data = read(path)["comparison"]
        require(data["first"]==bindings["Bf"] and data["second"]==bindings[name] and not data["failures"],
                "directional comparison differs from the passing signed result")
        comparisons[name] = (data,identity(path))
    canonical,binding = load_case(root,"Bf",context,ref,20,deadline)
    require(binding==bindings["Bf"],"canonical directional state differs from its independent audit")
    samples = []
    for origin,field in (("Br3","traction"),("Br4","R"),("Br4","Pi"),
                         ("Bx2","traction"),("Bx3","R"),("Bx3","Pi")):
        result,comparison = comparisons[origin]
        row = result["rows"][f"whole/{field}"]
        location = row["component_locations"][row["component_MAX"].index(row["MAX"])]
        x,r,component = (location[k] for k in ("x","r","component"))
        base, = (owner for owner in canonical.elements if owner["geometry"]==location["coarse_owner"])
        value = quantity(raw_owner(canonical,base,x,r),field,component)
        bases[origin,field] = (base,value,location)
        samples.append(dict(origin=origin,comparison=comparison,field=field,component=component,x=x,r=r,
            base_owner=base["geometry"],base_bounds=base["bounds"],effects={},
            inherited_limit=po.NUMERICAL[field]))
    del canonical
    for name in (*radial,*streamwise):
        state,binding = load_case(root,name,context,ref,20,deadline)
        require(binding==bindings[name],"directional state differs from the passing independent audit")
        for sample in samples:
            origin,field,component,x,r = (sample[key] for key in ("origin","field","component","x","r"))
            base,value,location = bases[origin,field]
            matches = [e for e in state.elements if common_owner(e,base,x,r)]
            require(matches,"signed point has no owner on the retained side")
            effects = [dict(owner=e["geometry"],bounds=e["bounds"],
                signed=(quantity(raw_owner(state,e,x,r),field,component)-value)/ref.tau) for e in matches]
            if name==origin:
                require(any(v["owner"]==location["fine_owner"] and
                            abs(v["signed"]-location["signed_difference"])<=2e-12 for v in effects),
                        "directional replay differs from the signed component extremum")
            sample["effects"][name] = effects
        del state
    for sample in samples:
        magnitudes = [max(abs(v["signed"]) for name in group for v in sample["effects"][name])
                      for group in (radial,streamwise)]
        sample["sum_of_directional_magnitudes"] = sum(magnitudes)
    atomic(out,dict(schema="WJH27OutletDirectionalOwnerProbeV1",bindings=bindings,
        sampling="Original signed stress, pressure, and shear extrema on their retained one-sided physical owners.",
        samples=samples,coupled_follow_up_indicated=any(v["sum_of_directional_magnitudes"]>v["inherited_limit"]
                                                       for v in samples)))
    print([(v["origin"],v["component"],v["sum_of_directional_magnitudes"]) for v in samples],flush=True)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root",type=Path)
    parser.add_argument("out",type=Path)
    parser.add_argument("seconds",type=int)
    args = parser.parse_args()
    check(args.root,args.out,args.seconds)
