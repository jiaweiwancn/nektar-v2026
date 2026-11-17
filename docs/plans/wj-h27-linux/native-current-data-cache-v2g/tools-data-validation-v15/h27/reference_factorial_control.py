#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Inactive Linux-only transport/cutoff column controls; never a runtime table."""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
import argparse
from dataclasses import replace
import hashlib
import json
import math
from pathlib import Path
import sys
import time
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
import fast_reference as fast
import reference_oracle as reference
from physical_oracle import Reference, TABLE_SHA, identity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--table', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--seconds', type=float, required=True)
    args = parser.parse_args()
    if not sys.platform.startswith('linux') or not 0 < args.seconds <= 120:
        raise ValueError('Linux and explicit bounded analysis allowance required')
    if args.out.exists() or identity(args.table)['sha256'] != TABLE_SHA:
        raise ValueError('create-only output and pinned v2 table required')
    began = time.monotonic()
    deadline = began + args.seconds
    generator = fast.Generator(args.seconds)
    old = Reference(json.loads(args.table.read_text()))
    historical_logc = math.log(old.old_c / old.H)
    trials, rows = [], []
    report = dict(schema='WJH27InactiveReferenceFactorialV1', status='NOT_MET',
                  proposal_only=True, active_model_changed=False,
                  runtime_table_produced=False, physical_reference_qualified=False,
                  full_H27_qualification=False, TASK3_READY=False,
                  scope='Independent Kin columns only; no new actual-flow PDE or history',
                  sources=dict(control=identity(__file__), fast=identity(fast.__file__),
                               oracle=identity(reference.__file__)),
                  existing_table=identity(args.table), historical_logc=historical_logc,
                  trials=trials, rows=rows)

    def clock():
        if time.monotonic() >= deadline:
            raise TimeoutError('inactive reference control deadline')

    def trial(law, logc, qh):
        clock()
        record = dict(diffusion=law, logc=logc, reference_QH_over_A=qh)
        trials.append(record)
        try:
            p = replace(generator.p, inner_diffusion=law, endpoint_flux_over_A=qh)
            column = reference.Column(p)
            if not column.zv < 0:
                raise ValueError('historical negative roughness branch lost')
            answer = column.solve_trial(logc)
            record.update(status='PASS_IVP_ONLY', ground_flux_over_A=answer['ground_flux'])
            return column, answer
        except Exception as error:
            record.update(status='REJECTED', reason=str(error))
            raise

    def assess(law, cutoff_policy, logc, qh, column, answer):
        breaks = sorted(set(answer['breaks'] + [old.old_c / old.H]))
        integrated = sum(quad(lambda y: column.state(y, answer['sample'](y)[0], logc)['production'],
                              a, b, epsabs=column.A*2e-12, epsrel=2e-11)[0]
                         for a, b in zip(breaks[:-1], breaks[1:])) * old.H
        budget = 1-logc-integrated/column.A-qh
        gauss, weights = np.polynomial.legendre.leggauss(64)
        radius = np.concatenate([(a+b)/2+(b-a)*gauss/2 for a, b in zip(breaks[:-1], breaks[1:])])
        weight = np.concatenate([weights*(b-a)/2 for a, b in zip(breaks[:-1], breaks[1:])])
        endpoint = np.unique(np.r_[np.linspace(0, 1, 257), breaks])
        radius = np.r_[radius, endpoint]
        weight = np.r_[weight, np.zeros(len(endpoint))]
        state = np.array([answer['sample'](float(y)) for y in radius])
        historical = old.inlet(radius*old.H)
        delta_k = state[:, 0]-historical['k']/old.Kin
        U = answer['anchor']+column.us*state[:, 2]
        delta_u = U-historical['U']
        physical = [column.state(float(y), float(k), logc) for y, k in zip(radius, state[:, 0])]
        minimum_k = float(np.min(state[:, 0]))
        minimum_tangent = min(v['shear_tangent'] for v in physical)
        max_traction = max(abs(v['traction_error'])/column.t for v in physical)
        row = dict(diffusion=law, cutoff_policy=cutoff_policy, logc=logc,
                   c_over_H=math.exp(logc), reference_QH_over_A=qh,
                   ground_flux_over_A=answer['ground_flux'], budget_over_A=budget,
                   minimum_k_over_K=minimum_k, minimum_shear_tangent=minimum_tangent,
                   traction_error_over_tau=max_traction,
                   ground_k_over_K=float(answer['sample'](0.)[0]),
                   ground_k_departure_over_K=float(answer['sample'](0.)[0])-float(old.inlet(np.array([0.]))['k'][0])/old.Kin,
                   inner_k_departure_RMS_over_K=math.sqrt(float(weight@(delta_k*delta_k))),
                   inner_k_departure_sampled_MAX_over_K=float(np.max(np.abs(delta_k))),
                   inner_U_departure_RMS=math.sqrt(float(weight@(delta_u*delta_u))),
                   inner_U_departure_sampled_MAX=float(np.max(np.abs(delta_u))),
                   velocity_anchor=answer['anchor'],
                   interpretation='Counterfactual column versus unchanged inlet; not a 2D application result')
        row['column_controls_pass'] = (abs(answer['ground_flux']) <= 1e-10 and abs(budget) <= 1e-10
                                       and minimum_k > 0 and minimum_tangent > 0 and max_traction <= 1e-10)
        rows.append(row)
        if not row['column_controls_pass']:
            raise ValueError('inactive column controls failed')

    try:
        for law in ('scalar', 'auxiliary'):
            generator.p = replace(generator.p, inner_diffusion=law, endpoint_flux_over_A=0.)
            generator.cache.clear()
            calibration = generator.fast(1.)
            logc = calibration['logc']
            column, answer = trial(law, logc, 0.)
            assess(law, 'full-WJ zero-endpoint-flux calibration', logc, 0., column, answer)
            # Separate inactive control: hold the historical cutoff and solve
            # its counterfactual endpoint flux; this never constrains actual QH.
            values = []
            for qh in (0., .025, .05, .075, .10, .125, .15, .20, .3, .4):
                try:
                    col, result = trial(law, historical_logc, qh)
                except ValueError:
                    continue
                values.append((qh, result['ground_flux']))
                if len(values) > 1 and values[-1][1]*values[-2][1] <= 0:
                    break
            bracket = next(((a[0], b[0]) for a, b in zip(values[:-1], values[1:]) if a[1]*b[1] <= 0), None)
            if bracket is None:
                raise ValueError('no bracket in explicitly admitted reference QH/A interval [0,.4]')
            qh = brentq(lambda q: trial(law, historical_logc, q)[1]['ground_flux'], *bracket,
                        xtol=2e-13, rtol=2e-13)
            column, answer = trial(law, historical_logc, qh)
            assess(law, 'historical cutoff, independently solved reference QH', historical_logc, qh, column, answer)
        report['status'] = 'PASS_INACTIVE_COLUMN_CONTROLS_ONLY'
    except Exception as error:
        report.update(exception=type(error).__name__, reason=str(error))
    report['elapsed_seconds'] = time.monotonic()-began
    with args.out.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({key: report[key] for key in ('status', 'elapsed_seconds')}))
    if 'exception' in report:
        raise SystemExit(2)


if __name__ == '__main__':
    main()
