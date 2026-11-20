#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Reconstruct the two polynomial scalar caches in exact reference geometry."""

import argparse
from pathlib import Path
import time

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse.linalg import splu

import physical_oracle as po
from physical_transient_oracle import state_hash
from resources import atomic, identity


def tensor(x, r):
    bx, dx = po.modified(3, x)
    br, dr = po.modified(3, r)
    combine = lambda a,b: (a[:,:,None]*b[:,None,:]).reshape(len(x),16)
    return combine(br,bx), combine(br,dx), combine(dr,bx)


def polynomial_transport(model, order, deadline):
    q, w = (np.asarray(a,dtype=np.longdouble) for a in leggauss(order))
    x, r = np.meshgrid(q,q)
    val, dx, dr = tensor(x.ravel(),r.ravel())
    ww = np.outer(w,w).ravel()
    edge_values = tuple(tensor(a,b)[0] for a,b in (
        (np.full_like(q,-1),q), (np.ones_like(q),q),
        (q,np.full_like(q,-1)), (q,np.ones_like(q))))
    result = {k:np.zeros_like(model.s[k],dtype=np.longdouble) for k in ("k","epsilon")}

    def local(e,key,shift=0):
        field = "flow" if key == "u" else key
        return np.asarray(e[po.SIGNS[key]],dtype=np.longdouble)*model.s[field][shift+np.asarray(e[key])]

    for e in model.elements:
        po.require(time.monotonic() < deadline,"scalar cache control exhausted its time allocation")
        xl,xr,rl,rr = e["bounds"]
        u = val@local(e,"u")
        v = val@local(e,"u",model.n)
        for key in result:
            if not e[key]:
                continue
            scalar = val@local(e,key)
            action = -(rr-rl)/2*(dx.T@(ww*u*scalar))-(xr-xl)/2*(dr.T@(ww*v*scalar))
            result[key][e[key]] += np.asarray(e[po.SIGNS[key]])*action
    for key in result:
        for e, side in model.faces[key]:
            xl,xr,rl,rr = e["bounds"]
            b = edge_values[side]
            shift = model.n if side >= 2 else 0
            sign = -1 if side in (0,2) else 1
            length = rr-rl if side < 2 else xr-xl
            flux = sign*(b@local(e,"u",shift))*(b@local(e,key))
            result[key][e[key]] += np.asarray(e[po.SIGNS[key]])*(length/2)*(b.T@(w*flux))
    return result


def measure(native, table, out, seconds):
    began = time.monotonic()
    deadline = began+seconds
    ref = po.load_reference(table)
    state, _ = po.load_native(native)
    saved = state["current_advection"]
    po.require(state_hash(saved) == state["advection_sha256"],"invalid native cache checksum")
    model = po.PhysicalOracle(state,ref,20,deadline)
    old, fine = (polynomial_transport(model,q,deadline) for q in (8,12))
    model.build_riesz()
    diagnostics = {}
    for key in fine:
        po.require(time.monotonic() < deadline,"scalar cache control exhausted its time allocation")
        factor = splu(model.masses[key])
        scale = model.volume if key == "k" else model.evolume
        delta = np.asarray(np.asarray(saved[key],dtype=np.longdouble)-fine[key],dtype=float)
        quadrature = np.asarray(fine[key]-old[key],dtype=float)
        norm = lambda d:np.sqrt(max(0.,d@factor.solve(d))/scale)
        peaks = np.argsort(abs(delta))[-6:][::-1]
        rows = []
        for i in peaks:
            owners = [dict(geometry=e["geometry"],bounds=e["bounds"]) for e in model.elements if i in e[key]]
            rows.append(dict(dof=int(i),difference=float(delta[i]),owners=owners))
        diagnostics[key] = dict(native_minus_reference_dual=float(norm(delta)),
            reference_order_difference_dual=float(norm(quadrature)),peak_coefficients=rows)
    atomic(out,dict(schema="WJH27R3ReferenceScalarCacheControlV1",native=identity(native),table=identity(table),
        rule_orders=[8,12],limit=1e-10,elapsed_seconds=time.monotonic()-began,cache_differences=diagnostics,
        note="Degree at most nine on complete affine cells and exterior faces. Local reference coordinates "
             "give direct tensor integrals with signed global assembly. Scalar transported fields contain no "
             "cutoff-dependent source or wall trace evaluation inside the transport integrand."))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("native","table","out"):
        p.add_argument("--"+key,type=Path,required=True)
    p.add_argument("--seconds",type=float,required=True)
    measure(**vars(p.parse_args()))
