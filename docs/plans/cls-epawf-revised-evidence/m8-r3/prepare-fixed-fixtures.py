#!/usr/bin/env python3
"""Construct fixed operator inputs on Linux; never a trajectory or field repair."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

sys.dont_write_bytecode = True
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'


def pin(path):
    path = Path(path).resolve()
    return dict(path=str(path), bytes=path.stat().st_size,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('source', 'reference', 'session', 'history', 'output'):
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args()
    if platform.system() != 'Linux' or a.output.exists():
        raise ValueError('Linux and a fresh external output directory required')
    source = a.source.resolve()
    if subprocess.check_output(['git', '-C', str(source), 'status', '--porcelain'], text=True).strip():
        raise ValueError('Clean analysis source required')
    sys.path.insert(0, str(source/'solvers/IncNavierStokesSolver/Tests/RANS'))
    import numpy as np
    from cls_epawf_m8_r2_stage_analysis import Analyzer
    raw = json.loads(a.reference.read_text(), parse_int=lambda x: -0.0 if x == '-0' else int(x))
    h0 = np.fromfile(a.history, dtype='<f8')
    if len(h0) != 2165 or h0.tobytes() != np.asarray(raw['initial'], dtype='<f8').tobytes():
        raise ValueError('Original outer6 h0 binding differs')
    a.output.mkdir(parents=True)
    ledger = dict(schema=1, scope='Fixed-state operator verification only; no physical qualification',
        source_commit=subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip(),
        inputs=[pin(a.reference), pin(a.session), pin(a.history)], outputs=[])

    def save(name, state, **extra):
        record = dict(schema=1, state=np.asarray(state, float).tolist(), history=h0.tolist(),
            alpha=400., time=.0025, inner_stress_model='SaturatedR2',
            constructed_input=True, physical_time_sample=False, **extra)
        path = a.output/(name+'.json')
        path.write_text(json.dumps(record, separators=(',', ':'), allow_nan=False)+'\n')
        ledger['outputs'].append(pin(path))

    save('h0', h0, scope='Exact original represented initial state; no state changes')
    save('constructed-contact', h0, construct_contact=True, compute_scalar_rate=True,
         scope='One preregistered contact polynomial is constructed by the native probe; not an accepted state')

    # A declared positive polynomial test INPUT is not a model field floor or
    # a modification of a computed update. It deliberately lies outside the
    # restricted outer box to exercise cap/source/trace derivatives.
    K0, Q, H, nu = (raw['parameters'][i] for i in (5, 6, 1, 0))
    offsets, fields = raw['offsets'], raw['fields']
    state = h0.copy(); one = np.asarray([1., 1., 0., 0.])
    values = {}
    def fx(x):
        return x+4 if x < -3 else 2-x if x > 1 else 1.
    for e in fields[2]['mesh']:
        if e['degree'] != 3:
            raise ValueError('Frozen fixture requires p3 k')
        xl, xu, rl, ru = (e[edge][axis] for edge, axis in
            [('lower', 0), ('upper', 0), ('lower', 2), ('upper', 2)])
        f = np.asarray([fx(xl), fx(xu), 0., 0.])
        if ru <= .5+1e-14:
            g = .01*one
        elif abs(rl-.5)<1e-14 and abs(ru-1)<1e-14:
            g = np.asarray([.01, 1., -1.5*.99, -.25*.99])
        elif abs(rl-1)<1e-14 and abs(ru-2)<1e-14:
            g = np.asarray([1., .01, -1.5*.99, .25*.99])
        elif abs(rl-2)<1e-14:
            if abs(ru-7/3)>1e-12:
                raise ValueError('Unexpected first outer interval')
            g = np.asarray([.01, 1., 0., 0.])
        elif rl >= 7/3-1e-12:
            g = one.copy()
        else:
            raise ValueError('Unexpected radial fixture interval')
        for j, (index, sign) in enumerate(zip(e['ids'], e['sign'])):
            ix, ir = j % 4, j // 4
            value = K0*(one[ix]*one[ir]-f[ix]*(one[ir]-g[ir]))*sign
            index = int(index)
            if index in values and abs(values[index]-value)>1e-15:
                raise ValueError('Nonconforming constructed polynomial')
            values[index] = value
    for index, value in values.items():
        state[offsets[2]+index] = value

    scale2 = K0/Q+(nu/H)**2
    pressure_values = {}
    for e in fields[5]['mesh']:
        if e['degree'] != 1:
            raise ValueError('Frozen fixture requires p1 reduced pressure')
        for j, (index, sign) in enumerate(zip(e['ids'], e['sign'])):
            x = e['lower' if j % 2 == 0 else 'upper'][0]
            kH = K0*(1-.99*fx(x))
            increment = (2/3)*(kH-K0)+.06*scale2/H*(x+1)
            value = sign*increment
            index = int(index)
            if index in pressure_values and abs(pressure_values[index]-value)>1e-15:
                raise ValueError('Nonconforming constructed pressure')
            pressure_values[index] = value
    for index, increment in pressure_values.items():
        state[offsets[5]+index] += increment
    analyzer = Analyzer(raw, a.session, source)
    constraints = analyzer.constraints(state)
    # Only the constructed input's dependent matching trace is prescribed;
    # neither saved h0 nor a Newton update is changed.
    for row in analyzer.trace_rows:
        state[row] = constraints[row]
    save('active-sampledK-source-trace', state,
         scope='Positive constructed operator fixture; outside the fixed physical applicability box',
         polynomial='k=K0*(1-f(x)*(1-g(r))); f piecewise linear, g declared piecewise cubic/linear',
         prescribed_polynomial_minimum_ratio=.01, pressure_sensor_dose=.06,
         trace_construction='Independent actual projected dependent epsilon; only matching coefficients set',
         native_scalar_domain_gate_required=True, active_cap_and_sensor_measured_not_assumed=True)
    (a.output/'receipt.json').write_text(json.dumps(ledger, indent=2)+'\n')
    print(json.dumps(dict(success=True, packets=len(ledger['outputs']))))


if __name__ == '__main__':
    main()
