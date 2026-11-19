#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check raw outlet shear and physical derivatives with a second Q3 evaluator."""

import argparse
from decimal import Decimal, localcontext
from pathlib import Path
import time

import numpy as np

import physical_oracle as po
from resources import atomic
from task2_outlet_campaign import SHEAR_TRACK
from task2_outlet_case import get_context
from task2_outlet_qualification import load_case, raw_owner
from task2r_run import require


def modes(t):
    return (1-t,t,t-t*t,-t+3*t*t-2*t**3), (-1,1,1-2*t,-1+6*t-6*t*t)


def gradient(state, owner, x, r):
    with localcontext() as work:
        work.prec = 45
        as_decimal = lambda x: Decimal.from_float(float(x))
        xl,xr,rl,rr = map(as_decimal,owner["bounds"])
        width_x,width_r = xr-xl,rr-rl
        bx,dx = modes((as_decimal(x)-xl)/width_x)
        br,dr = modes((as_decimal(r)-rl)/width_r)
        result = {}
        for field,offset in (("U",0),("W",state.n)):
            coeff = [as_decimal(state.s["flow"][offset+index])*int(sign)
                     for index,sign in zip(owner["u"],owner["us"])]
            result[field+"x"] = float(sum(coeff[4*j+i]*br[j]*dx[i]
                                         for j in range(4) for i in range(4))/width_x)
            result[field+"r"] = float(sum(coeff[4*j+i]*dr[j]*bx[i]
                                         for j in range(4) for i in range(4))/width_r)
        return result


def probe(state, owner, x, r, label):
    raw = raw_owner(state,owner,x,r)
    back = gradient(state,owner,x,r)
    defect = max(abs(raw[k]-back[k]) for k in back)
    require(defect<=1e-10, "signed independent Q3 derivative check failed")
    nu_A = float(state.ref.aux(np.array([raw["K"]]),np.array([r]))[0])
    return dict(x=x,r=r,label=label,raw=raw,reconstructed_gradient=back,
                max_derivative_error=defect,nu_A=nu_A)


def is_owner(owner, x, r):
    xl,xr,rl,rr = owner["bounds"]
    return xl-1e-12<=x<=xr+1e-12 and rl-1e-12<=r<=rr+1e-12


def samples(state):
    rows, tracked = [], []
    for x,r in SHEAR_TRACK:
        matches = [probe(state,e,x,r,"tracked")
                   for e in state.elements if is_owner(e,x,r)]
        require(matches,"tracked observation has no physical owner")
        tracked.extend(matches)
    columns = sorted({tuple(e["bounds"][:2]) for e in state.elements})
    relevant = set(columns[:2]+columns[-2:])
    for e in state.elements:
        a,b,c,d = e["bounds"]
        if (a,b) not in relevant or c>=.071875+1e-12 or d<=.05-1e-12:
            continue
        for t in (0.,.5,1.):
            for s in (0.,.5,1.):
                if .05-1e-12<=c+s*(d-c)<=.071875+1e-12:
                    rows.append(probe(state,e,a+t*(b-a),c+s*(d-c),"local"))
    require(rows,"adjacent x or radial owner probes are absent")
    return dict(tracked=tracked,local=rows,max_derivative_error=max(
        row["max_derivative_error"] for row in [*rows,*tracked]))


def shear(a, b, reference, radial_side):
    selected = []
    for v in (a,b):
        index = 3 if radial_side == "lower" else 2
        row, = [item for item in v["tracked"] if item["x"]==15. and item["r"]==.06
                and abs(item["raw"]["bounds"][index]-.06)<1e-12]
        require(row["raw"]["chi"]==0. and not row["raw"]["E_type"]=="transported",
                "wall shear identity is scoped to chi=0 in the inner band")
        q = row["raw"]
        effective = reference.nu+row["nu_A"]
        value = effective*(q["Ur"]+q["Wx"])
        require(abs(value-q["Vxr"])<1e-12*reference.tau, "raw shear identity differs")
        selected.append(dict(owner=q["owner"],bounds=q["bounds"],Ur=q["Ur"],Wx=q["Wx"],
                             viscosity=effective,nu_A=row["nu_A"],Vxr=q["Vxr"],K=q["K"]))
    p,q = selected
    average_mu = (p["viscosity"]+q["viscosity"])/2
    average_shear = (p["Ur"]+p["Wx"]+q["Ur"]+q["Wx"])/2
    changes = dict(Ur=average_mu*(q["Ur"]-p["Ur"])/reference.tau,
        Wx=average_mu*(q["Wx"]-p["Wx"])/reference.tau,
        viscosity=(q["viscosity"]-p["viscosity"])*average_shear/reference.tau)
    observed = (q["Vxr"]-p["Vxr"])/reference.tau
    require(abs(sum(changes.values())-observed)<1e-12, "signed centred shear product difference differs")
    return dict(first=p,second=q,terms_over_tau0=changes,signed_Vxr_over_tau0=observed)


def inspect(root, base, case, out, seconds):
    context = get_context(root)
    ref = po.load_reference(context[2])
    deadline = time.monotonic()+seconds
    a,pa = load_case(root,base,context,ref,20,deadline)
    first = samples(a)
    b,pb = load_case(root,case,context,ref,20,deadline)
    second = samples(b)
    results = {side:shear(first,second,ref,side) for side in ("lower","upper")}
    atomic(out,dict(schema="WJH27OutletSignedQ3OwnerProbeV1",cases=[base,case],
        first=pa,second=pb,samples={base:first,case:second},chi_zero_shear_x15_r006=results,
        method="Signed tensor coefficients evaluated in unit affine coordinates at 45 decimal digits.",
        derivative_verification_absolute_limit=1e-10,
        acceptance_gradients="Raw physical-coordinate gradients",endpoint_in_MAX=True,
        shear_normalization=ref.tau))
    print(base,case,{side:v["signed_Vxr_over_tau0"] for side,v in results.items()},flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root",type=Path)
    parser.add_argument("base")
    parser.add_argument("case")
    parser.add_argument("output",type=Path)
    parser.add_argument("seconds",type=float)
    args=parser.parse_args()
    inspect(args.root,args.base,args.case,args.output,args.seconds)
