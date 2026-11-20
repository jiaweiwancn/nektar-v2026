# SPDX-License-Identifier: MIT
"""Independent conservative Q3 scalar transport in affine reference coordinates."""

import numpy as np
from numpy.polynomial.legendre import leggauss


SIGNS = {"u": "us", "k": "ks", "epsilon": "epsilon_signs"}


def tensor(x, r):
    def modes(t):
        bubble = (1-t*t)/4
        return (np.column_stack(((1-t)/2,(1+t)/2,bubble,2*t*bubble)),
                np.column_stack((np.full_like(t,-.5),np.full_like(t,.5),-t/2,(1-3*t*t)/2)))

    bx, dx = modes(x)
    br, dr = modes(r)

    def combine(a, b):
        return (a[:,:,None]*b[:,None,:]).reshape(len(x),16)

    return combine(br,bx), combine(br,dx), combine(dr,bx)


def polynomial_transport(model, order):
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
        return np.asarray(e[SIGNS[key]],dtype=np.longdouble)*model.s[field][shift+np.asarray(e[key])]

    for e in model.elements:
        model.clock()
        xl,xr,rl,rr = e["bounds"]
        u = val@local(e,"u")
        v = val@local(e,"u",model.n)
        for key in result:
            if not e[key]:
                continue
            scalar = val@local(e,key)
            action = -(rr-rl)/2*(dx.T@(ww*u*scalar))-(xr-xl)/2*(dr.T@(ww*v*scalar))
            result[key][e[key]] += np.asarray(e[SIGNS[key]])*action
    for key in result:
        for e, side in model.faces[key]:
            xl,xr,rl,rr = e["bounds"]
            b = edge_values[side]
            shift = model.n if side >= 2 else 0
            sign = -1 if side in (0,2) else 1
            length = rr-rl if side < 2 else xr-xl
            flux = sign*(b@local(e,"u",shift))*(b@local(e,key))
            result[key][e[key]] += np.asarray(e[SIGNS[key]])*(length/2)*(b.T@(w*flux))
    return result
