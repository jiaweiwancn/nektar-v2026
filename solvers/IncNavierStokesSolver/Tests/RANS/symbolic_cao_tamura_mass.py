#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded Q3 signed-CG support graphs; no >3000-slot physical solver case."""
import argparse,time,json
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix,block_diag,diags
from scipy.sparse.linalg import splu
from continue_cls_epawf_pimple import atomic

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();out=[]
    for label,nz,no in [('current',13,7),('retained-localized',21,10),('retained-full-h',25,14)]:
        nx=6;px=3*nx+1;pz=3*nz+1;n=px*pz;rows=[];cols=[];started=time.perf_counter()
        for ex in range(nx):
            for ez in range(nz):
                ids=np.array([(3*ex+i)*pz+3*ez+j for i in range(4) for j in range(4)]);rows.extend(np.repeat(ids,16));cols.extend(np.tile(ids,16))
        pattern=coo_matrix((np.ones(len(rows)),(rows,cols)),shape=(n,n)).tocsr();pattern.sum_duplicates();pattern.data[:]=1
        ufree=np.array([i*pz+j for i in range(px) for j in range(pz) if i and j]);wfree=np.array([i*pz+j for i in range(px) for j in range(pz) if i and j and j<pz-1]);graph=block_diag([pattern[ufree][:,ufree],pattern[wfree][:,wfree]],format='csc');seconds=time.perf_counter()-started
        # Fill-only SPD surrogate has this graph and a strictly dominant diagonal.
        # Its values and numerical cost are NOT a target physical-mass result.
        surrogate=graph+diags(np.asarray(graph.sum(axis=1)).ravel()+1);t=time.perf_counter();factor=splu(surrogate,diag_pivot_thresh=0);fill_seconds=time.perf_counter()-t
        ni=nz-no;slots=3*n+19*(3*no+1)+7*(ni+1)+7*(no+1)+1
        out.append(dict(label=label,radial_panels=nz,outer_panels=no,estimated_coupled_slots=slots,velocity_field_size=n,free_velocity_size=len(ufree)+len(wfree),support_graph_nnz=int(graph.nnz),CSR_upper_bytes=int(graph.data.nbytes+graph.indices.nbytes+graph.indptr.nbytes),graph_seconds=seconds,fill_only_surrogate_L_U_nnz=int(factor.L.nnz+factor.U.nnz),fill_only_surrogate_seconds=fill_seconds,dense_free_mass_bytes=8*(len(ufree)+len(wfree))**2,dense_coupled_array_bytes=8*slots*slots))
    atomic(a.output,dict(version='H2fSymbolicMassGraphV1',production_guard=3000,physical_target_factorizations=0,cases=out,scope='support graphs of retained tensor Q3/P1 allocations; modified-basis signs do not change adjacency; curved/affine numerical zeros can reduce actual nnz; surrogate fill is an estimate, not physical target factor timing'))
if __name__=='__main__':main()
