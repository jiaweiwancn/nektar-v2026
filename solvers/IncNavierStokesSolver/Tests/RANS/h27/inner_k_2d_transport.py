#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Resolve the remaining inner-k production terms and spatial sampling error."""

import argparse
from bisect import bisect_right
from collections import defaultdict
import math
from pathlib import Path
import time

import physical_oracle as po
from inner_k_2d import checked_state
from inner_k_fixture_oracle import FixtureReference, atomic_out
from inner_k_mechanism import read
import numpy as np
from numpy.polynomial.legendre import leggauss


TERMS = ('P_minus_E','parallel_aux_minus_E','aux_shear_shift','isotropic_work','weighted_full_WJ_shift')
REGIONS = ('below_chi','blended_collar','full_WJ_below_H')


def signed_terms(data, reference, order, deadline):
    oracle = po.PhysicalOracle(data,reference,order,deadline)
    columns = defaultdict(lambda: np.zeros((len(REGIONS),len(TERMS))))
    for e in oracle.elements:
        if e['epsilon']:
            continue
        x,r,w = oracle.rule(e)
        a = oracle.sample(e,x,r)
        aux = reference.aux(a['K'],r)
        parallel = aux*((a['K']/reference.qe)/(reference.nu+aux))**2
        shear = aux*(a['G'][:,0,2]+a['G'][:,2,0])**2
        iso = -2*a['k']*a['div']/3
        actual = a['P']-a['E']
        values = np.array([actual,parallel-a['E'],shear-parallel,iso,a['P']-shear-iso])
        po.require(np.max(abs(actual-values[1:].sum(axis=0)))<2e-14,'signed split failed to reconstruct actual production')
        xspan = tuple(round(p,12) for p in e['bounds'][:2])
        for i, mask in enumerate((r<.4*reference.H,(r>=.4*reference.H)&(r<.75*reference.H),r>=.75*reference.H)):
            columns[xspan][i] += values[:,mask]@w[mask]
    total = sum(columns.values())
    scale = 25*reference.A
    rows = [dict(x=x, terms=v.tolist(), normalized=(v/((x[1]-x[0])*reference.A)).tolist())
            for x,v in sorted(columns.items())]
    po.require(abs(total[:,1].sum())/scale<=1e-10,'parallel auxiliary reference no longer balances dependent E')
    return dict(order=order,regions=REGIONS,terms=TERMS,columns=rows,total=total.tolist(),
                integrated_total_over_25A=(total/scale).tolist())


def cell_k(oracle, cell, x, r):
    return (po.basis(cell,'k',x,r)[0]@oracle.s['k'][cell['k']])/oracle.ref.Kin


def polynomial_control(coarse, sharp, deadline):
    bounds = {tuple(round(p,12) for p in e['bounds']):e for e in coarse.elements if not e['epsilon']}
    xs = sorted(set(v for b in bounds for v in b[:2]))
    rule, w = leggauss(8)
    probe_orders = (20,40)
    probes = [np.r_[-1.,leggauss(n)[0],1.] for n in probe_orders]
    square, area = 0., 0.
    peaks = [dict(signed=0.),dict(signed=0.)]
    for e in sharp.elements:
        if e['epsilon']:
            continue
        po.require(time.monotonic()<deadline,'k-only spatial analysis exceeded its cap')
        xl,xr,rl,rr = e['bounds']
        j = bisect_right(xs,(xl+xr)/2)-1
        original = bounds.get((xs[j],xs[j+1],round(rl,12),round(rr,12)))
        po.require(original and original['bounds'][0]<=xl+1e-12 and original['bounds'][1]>=xr-1e-12,
                   'directional check needs shared radial cells and nested x faces')
        def difference(points):
            x,r = np.meshgrid((xl+xr)/2+(xr-xl)*points/2,(rl+rr)/2+(rr-rl)*points/2)
            x,r = x.ravel(),r.ravel()
            return x,r,cell_k(sharp,e,x,r)-cell_k(coarse,original,x,r)
        da = (xr-xl)*(rr-rl)
        d = difference(rule)[2]
        square += float(np.outer(w,w).ravel()@(d*d))*da/4
        area += da
        for i, points in enumerate(probes):
            x,r,d = difference(points)
            q = np.argmax(abs(d))
            if abs(d[q])>abs(peaks[i]['signed']):
                peaks[i] = dict(x=float(x[q]),r=float(r[q]),signed=float(d[q]),
                                coarse_owner=original['geometry'],finer_owner=e['geometry'])
    po.require(abs(area-25*coarse.ref.H)<=1e-10,'inner common-owner measure differs')
    return dict(exact_polynomial_RMS_over_Kin=math.sqrt(square/area),sampled_peaks=peaks,
                probe_orders=probe_orders,growth=abs(peaks[1]['signed'])-abs(peaks[0]['signed']),
                sampling='Gauss points and exact endpoints on every common rectangle, including all original owners')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('table','gate','base','directional','out'):
        p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--seconds',type=float,required=True)
    args = p.parse_args()
    deadline = time.monotonic()+args.seconds
    po.require(args.seconds>0 and po.identity(args.table)['sha256']==po.TABLE_SHA,'invalid transport analysis inputs')
    ref,gate = FixtureReference(read(args.table)),read(args.gate)
    breakdown = {}
    for name, path in (('base',args.base),('directional',args.directional)):
        data,_ = checked_state(path,gate,ref)
        q20,q32 = [signed_terms(data,ref,q,deadline) for q in (20,32)]
        error = float(np.max(np.abs(np.asarray(q20['integrated_total_over_25A'])-q32['integrated_total_over_25A'])))
        po.require(error<=1e-10,'signed production quadrature differs')
        breakdown[name] = dict(native=po.identity(path),lower_order=q20,full_order=q32,fixed_scale_difference=error)
    a,b = [po.PhysicalOracle(checked_state(p,gate,ref)[0],ref,20,deadline) for p in (args.base,args.directional)]
    body = dict(schema='WJInnerKActualTransportDecompositionV1',breakdown=breakdown,
                common_owner_polynomial_control=polynomial_control(a,b,deadline), source=po.identity(__file__))
    atomic_out(args.out,body)
    print('INNER_K_2D_TRANSPORT_DECOMPOSED',flush=True)


if __name__=='__main__':
    main()
