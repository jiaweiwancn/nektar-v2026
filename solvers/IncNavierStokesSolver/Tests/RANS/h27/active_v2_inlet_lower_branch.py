#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Track the other end of the canonical fixed-K inlet column branch."""

import argparse
import math
from pathlib import Path
import sys
import time

import physical_oracle as po
from active_v2_inlet import Column
from active_v2_inlet_obstruction import lower_match
from inner_k_columns import planar_response
from inner_k_fixture_oracle import atomic_out
from inner_k_mechanism import read
import numpy as np
from scipy.optimize import brentq


def follow(column, t, guess):
    a = column.ref
    stop = .4*a.H
    T = t*a.tau
    target,_,_ = lower_match(column,T)
    upper = target+math.log(a.H/stop)
    lower = upper-(a.H-stop)*T*T/(4*a.nu*a.A)
    attempts = []

    def test(q):
        try:
            end,_ = column.integrate([t,q],tight=True,lower=stop)
        except ValueError as error:
            attempts.append(dict(q=q,rejected=str(error)))
            return None
        f = float(end[1]-target)
        attempts.append(dict(q=q,defect=f,slope=float(end[6])))
        return f

    start = max(lower,guess-.04)
    f0 = test(start)
    while (f0 is None or f0>0) and start>lower:
        start = max(lower,start-.05)
        f0 = test(start)
    result = dict(T_over_tau=t,trials=attempts)
    if f0 is None or f0>0:
        return {**result,'status':'NO_VALID_LOWER_BRACKET'}
    last = start
    bound = None
    for q in np.linspace(start,upper,81)[1:]:
        f = test(float(q))
        if f is None:
            bound = float(q)
            for _ in range(26):
                q = (last+bound)/2
                f = test(q)
                if f is None:
                    bound = q
                elif f>=0:
                    break
                else:
                    last = q
            if f is None or f<0:
                return {**result,'status':'NO_ZERO_BEFORE_CONNECTED_SHOT_LIMIT',
                        'connected_q_limit_bracket':[last,bound]}
        if f>=0:
            root = brentq(test,last,q,xtol=2e-14)
            break
        last = float(q)
    else:
        return {**result,'status':'NO_ZERO_IN_ANALYTIC_BOUND'}

    try:
        state,sample = column.integrate([t,root],tight=True)
    except ValueError as error:
        return {**result,'status':'FIXED_SHEAR_ROOT_FAILS_FULL_COLUMN',
                'Q_H_over_A':root,'reason':str(error)}
    r = np.unique(np.r_[np.linspace(0.,a.H,257),column.panels])
    k = a.Kin*np.array([sample(x)[0] for x in r])
    v = np.array([column.shear(x,y/a.Kin,T) for x,y in zip(r,k)])
    guard,_ = planar_response(k,v[:,1],v[:,0],v[:,5],v[:,6],a.nu)
    po.require(abs(state[1])<=1e-10,'lower branch failed ground flux identity')
    return {**result,'status':'CHECKED_ROOT','Q_H_over_A':root,'U_H':float(-state[2]),
            'min_k':float(np.min(k)),'Q_0_over_A':float(state[1]),'guard':guard}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('table','branch','out'):
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--seconds',type=float,required=True)
    args = p.parse_args()
    table = po.identity(args.table)
    prior = read(args.branch)
    po.require(sys.platform.startswith('linux') and args.seconds>0 and table['sha256']==po.TABLE_SHA
               and prior['table']==table and prior['status']=='BLOCKED','requires reviewed fixed-K branch')
    ref = po.Reference(read(args.table))
    column = Column(ref,time.monotonic()+args.seconds)
    args.out.mkdir()
    guess = .19094601735826994
    results = []
    try:
        for t in np.r_[np.linspace(.94,.87,15),.855,.84]:
            row = follow(column,float(t),guess)
            results.append(row)
            if 'Q_H_over_A' in row:
                guess = row['Q_H_over_A']
        roots = [r for r in results if r['status']=='CHECKED_ROOT']
        po.require(roots and max(r['U_H'] for r in roots)<prior['target_U_H'],
                   'lower-shear continuation supplies an alternate candidate test')
        atomic_out(args.out/'lower-branch-assessment.json',dict(schema='WJH27MatchedInletLowerBranchV1',
            table=table,prior=po.identity(args.branch),target=prior['target_U_H'],
            maximum_observed_U_H=max(r['U_H'] for r in roots),
            observations=[{k:v for k,v in r.items() if k!='trials'} for r in results],
            scope='Finite follow of the canonical shooting interval toward lower shear.'))
    finally:
        atomic_out(args.out/'lower-branch-trials.json',results)
    print('ACTIVE_V2_INLET_LOWER_BRANCH_CHECKED',flush=True)


if __name__=='__main__':
    main()
