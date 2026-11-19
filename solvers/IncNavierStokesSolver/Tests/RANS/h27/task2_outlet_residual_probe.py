#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Localize the independent momentum defect on a fixed typed flat state."""

import argparse
import math
from pathlib import Path
import time

import numpy as np
from numpy.polynomial.legendre import leggauss

import physical_oracle as po
from resources import atomic, identity
from task2_outlet_baseline import read, verified
from task2_outlet_case import get_context
from task2_outlet_qualification import GRACE, PUBLICATION_SECONDS, load_case
from task2r_run import balance, require


def freeze(root, case, out):
    path = root/"evidence"/f"outlet-{case}-audit.json"
    audit = read(path)
    state = audit["provenance"]["state"]
    verified(state)
    require(audit["status"]["independent_current"]=="NOT_MET"
            and audit["core"]["residuals"]["momentum"]>po.RESIDUAL_LIMITS["momentum"],
            "the fixed-state diagnostic requires the measured momentum miss")
    receipt = root/"operations"/f"outlet-qual-{case.lower()}-audit"/"receipt.json"
    cost = read(receipt)["result"]["seconds"]
    orders = [30,8]
    cap = math.ceil(1.24*cost*sum((q/24)**2 for q in orders)+25.)
    required = cap+GRACE+57+PUBLICATION_SECONDS
    live = balance(read(root/"accounting.json"))
    require(required<=min(live["nonprotected_remaining"],live["tranche_remaining"]),
            "full residual diagnostic and publication require more time")
    atomic(out,dict(schema="WJH27FixedStateMomentumProbeForecastV1",case=case,
        state=state,audit=identity(path),cost_anchor=identity(receipt),orders=orders,
        prior_independent_order=24,once_built_Riesz_order=orders[0],
        prediction="If independent quadrature convergence is already reached, higher order will retain the missed norm. The dual energy location will identify the owners used in the next numerical remedy.",
        analysis_cap_seconds=cap,required_seconds=required,live=live))
    print(case,"fixed-state momentum probe",cap,required,flush=True)


def local_energy(state, residual):
    rhs = state.total["flow"]
    representer = np.zeros((2,state.n))
    for component, ids in enumerate(state.ufree):
        representer[component,ids] = state.factors["u"][component].solve(rhs[component*state.n+ids])
    a,w = leggauss(4)
    weights = np.outer(w,w).ravel()
    rows, groups = [], dict(inlet=0.,terminal=0.,remainder=0.)
    for owner in state.elements:
        state.clock()
        xl,xr,rl,rr = owner["bounds"]
        x = np.repeat((xl+xr)/2+(xr-xl)*a/2,len(a))
        r = np.tile((rl+rr)/2+(rr-rl)*a/2,len(a))
        v = po.basis(owner,"u",x,r)[0]@representer[:,owner["u"]].T
        energy = float(weights@np.sum(v*v,axis=1)*(xr-xl)*(rr-rl)/4)
        rows.append(dict(owner=owner["geometry"],bounds=owner["bounds"],energy=energy))
        group = "inlet" if xr<=-9.84375+1e-12 else "terminal" if xl>=14.9951171875-1e-12 else "remainder"
        groups[group] += energy
    exact = residual**2*state.volume
    observed = sum(groups.values())
    require(abs(exact-observed)<=1e-8*max(exact,1e-30),"partitioned momentum dual energy differs")
    top = sorted(rows,key=lambda row:row["energy"],reverse=True)[:20]
    work = {}
    for term,part in state.parts.items():
        flow = part["flow"][:2*state.n].reshape(2,state.n)
        work[term] = float(np.sum(flow*representer))
    return dict(total_energy=observed,relative_energy_defect=abs(observed-exact)/exact,
        fractions={key:value/observed for key,value in groups.items()},top_owners=top,
        signed_weak_term_work=work)


def probe(root, design, out, seconds):
    plan=read(design)
    context=get_context(root)
    ref=po.load_reference(context[2])
    deadline=time.monotonic()+seconds
    orders=plan["orders"]
    state,binding=load_case(root,plan["case"],context,ref,orders[0],deadline)
    require(binding["state"]==plan["state"],"independent quadrature probe changed the state")
    peak={}
    for q in orders:
        state.order=q
        state.gauss=leggauss(q)
        result=state.assemble()
        residuals=result["residuals"]
        peak[str(q)]=dict(residuals=residuals,accounting=result["accounting_identity_status"],
                         momentum_dual=local_energy(state,residuals["momentum"]))
    audit=read(verified(plan["audit"]))
    atomic(out,dict(schema="WJH27FixedStateMomentumProbeV1",state=binding["state"],case=plan["case"],
        forecast=identity(design),prior_order=plan["prior_independent_order"],
        prior_residuals=audit["core"]["residuals"],measurements=peak,once_built_Riesz_order=orders[0],
        native_quadrature=binding["quadrature"],native_tolerances=binding["tolerances"],
        scientific_scope="Independent diagnostic on the archived field. The declared acceptance audit is retained."))
    print(plan["case"],{q:v["residuals"]["momentum"] for q,v in peak.items()},flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root",type=Path)
    modes=parser.add_subparsers(dest="mode",required=True)
    plan=modes.add_parser("freeze")
    plan.add_argument("case")
    plan.add_argument("out",type=Path)
    check=modes.add_parser("probe")
    check.add_argument("design",type=Path)
    check.add_argument("out",type=Path)
    check.add_argument("seconds",type=float)
    args=parser.parse_args()
    if args.mode=="freeze":
        freeze(args.root,args.case,args.out)
    else:
        probe(args.root,args.design,args.out,args.seconds)
