#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Seed the signed current model with descent from its active block gradients."""
import argparse
import hashlib
from pathlib import Path

import numpy as np
from scipy.optimize import minimize, minimize_scalar

import c1_bounded_optimizer as bounded

MODEL_SHA256 = "d6ed3495eec4109d9ddbeace67125646e8aa3e79dc954a719197397d7cd0c811"


def active_seed(model, radius):
    norms = [np.linalg.norm(b) for b, _ in model.parts]
    top = max(norms)
    gradients = np.vstack([a.T @ b / n for (b, a), n in zip(model.parts, norms)
                           if n >= top*(1.-64.*np.finfo(float).eps)])
    gram = gradients @ gradients.T
    uniform = np.full(len(gram), 1./len(gram))
    weights = minimize(lambda x: .5*x @ gram @ x, uniform, jac=lambda x: gram @ x,
                       method="SLSQP", bounds=[(0.,1.)]*len(gram),
                       constraints={"type":"eq","fun":lambda x:np.sum(x)-1.,"jac":lambda x:np.ones_like(x)},
                       options={"ftol":1.e-14,"maxiter":50})
    axis = -gradients.T @ (weights.x if weights.success else uniform)
    size = np.max(np.abs(model.rows @ axis))
    if not np.isfinite(size) or size == 0.:
        return None
    axis *= radius*(1.-1.e-11)/size

    def score(fraction):
        return max(np.linalg.norm(b+fraction*(a @ axis)) for b,a in model.parts)

    line = minimize_scalar(score, method="bounded", bounds=(0.,1.), options={"xatol":1.e-10})
    choices = [0., 1., .25, .0625, .00390625]
    if line.success:
        choices.append(float(line.x))
    fraction = min(choices,key=score)
    if fraction <= 0. or score(fraction)>=top:
        return None
    return fraction*axis


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("basis",type=Path)
    p.add_argument("radius",type=float)
    p.add_argument("output",type=Path)
    args=p.parse_args()
    if hashlib.sha256(Path(bounded.__file__).read_bytes()).hexdigest()!=MODEL_SHA256:
        raise ValueError("pinned bounded model implementation differs")
    bounded.choose(args.basis,args.radius,args.output,seed_model=active_seed)
