#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent reader for the frozen active-v2 inner profile."""

import physical_oracle as po
import numpy as np
from numpy.polynomial.chebyshev import chebder, chebval


SCHEMA = 'WJH27ActiveV2InnerInletV1'
MODEL = 'FiniteNuIsotropicWallFullWJCollarV2'
TRANSPORT = 'ScalarCdV2'


class Profile:
    def __init__(self, doc, reference):
        po.require(doc['schema']==SCHEMA and doc['reference_table_sha256']==po.TABLE_SHA and
                   doc['model']==MODEL and doc['transport']==TRANSPORT and doc['H']==reference.H and
                   doc['Kin']==reference.Kin and doc['nu']==reference.nu and
                   abs(doc['cutoff']-reference.cutoff(reference.Kin))<1e-14,
                   'inlet profile belongs to a different model/reference')
        self.H = reference.H
        self.panels = [(p['lo'],p['hi'],np.asarray(p['U']),np.asarray(p['k'])) for p in doc['panels']]
        self.edges = [p[0] for p in self.panels]+[self.panels[-1][1]]
        po.require(self.edges[0]==0. and self.edges[-1]==self.H and
                   all(a<b and b==c for (a,b,_,_),c in zip(self.panels,self.edges[1:])),
                   'invalid physical inner interpolation panels')

    def evaluate(self, r):
        r = np.asarray(r,dtype=float)
        po.require(np.all(np.isfinite(r)&(r>=0)&(r<=self.H)), 'inner profile radius outside [0,H]')
        U,k,G,kr = (np.empty_like(r) for _ in range(4))
        for i,(lo,hi,u,q) in enumerate(self.panels):
            mask = (r>=lo)&((r<hi) if i+1<len(self.panels) else (r<=hi))
            t = (2*r[mask]-lo-hi)/(hi-lo)
            U[mask],k[mask] = chebval(t,u),chebval(t,q)
            G[mask],kr[mask] = (2/(hi-lo))*chebval(t,chebder(u)),(2/(hi-lo))*chebval(t,chebder(q))
        return dict(U=U,k=k,G=G,kr=kr)
