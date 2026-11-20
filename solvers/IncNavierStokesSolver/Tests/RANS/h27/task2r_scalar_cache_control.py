#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Reconstruct the two polynomial scalar caches in exact reference geometry."""

import argparse
from pathlib import Path
import time

import numpy as np
from scipy.sparse.linalg import splu

import physical_oracle as po
from physical_transient_oracle import state_hash
from resources import atomic, identity
from scalar_reference_transport import polynomial_transport


def measure(native, table, out, seconds):
    began = time.monotonic()
    deadline = began+seconds
    ref = po.load_reference(table)
    state, _ = po.load_native(native)
    saved = state["current_advection"]
    po.require(state_hash(saved) == state["advection_sha256"],"invalid native cache checksum")
    model = po.PhysicalOracle(state,ref,20,deadline)
    old, fine = (polynomial_transport(model,q) for q in (8,12))
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
