#!/usr/bin/env python3
"""Check native transverse tangents against an independent parallel-shear formula.

This is a constitutive principal-part diagnostic, not a global stability proof
or a calibrated flow. Both physical sides of max(S,Omega) are retained.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import subprocess
import tempfile


def coefficient(q):
    if q > 8.:
        f, df = 1., 0.
    else:
        e = math.exp(.75*q)
        f, df = -math.expm1(-.36*e), .27*e*math.exp(-.36*e)
    d = 1.+.35*q**1.5
    return .3*f/d, .3*(df*d-f*.525*math.sqrt(q))/d**2


def transverse(q, theta, side, k, epsilon, nu):
    """a_i n_j dT_ij/dh for g_xz=G, dg=h*a*n, a perpendicular to n.

    For delta=k*h/epsilon, S'=cos(2theta), Omega'=-1 at q>0.
    Positive delta selects S; negative delta selects Omega. In the x-z
    trace-free block S^2 is a scalar matrix, so its commutator vanishes.
    Isotropic terms vanish in a_i n_j contraction. The remaining quadratic
    term is 0.2*k*Cmu*q*(q-delta)*sin(2theta); the cubic derivative is
    20*k*Cmu^3*q^2*cos(2theta)*(1+cos(2theta)).
    """
    c, dc = coefficient(q)
    x, y = math.cos(2*theta), math.sin(2*theta)
    dc_side = dc*x if side > 0 else -dc
    return nu+k*k/epsilon*(c*(1-.2*q*y)+dc_side*(q*x+.2*q*q*y)
                            +20*c**3*q*q*x*(1+x))


def check(audit, output):
    output.mkdir(parents=True, exist_ok=False)
    k, epsilon, nu = .02, .1, 3.8e-6
    qs = (0., .5, 1., 2., 2.8, 2.9, 3., 3.190675081624288, 3.2, 4., 5., 8., 9., 10.)
    source = output/'input.csv'
    with source.open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(('id', 'eta', 'k', 'epsilon', 'shear'))
        for i, q in enumerate(qs):
            writer.writerow((f'q-{i}', .05, k, epsilon, q*epsilon/k))
    results = {}
    for model in ('standard', 'cls'):
        destination = output/(model+'.csv')
        command = [str(audit.resolve()), model, str(nu), str(source), str(destination)]
        subprocess.run(command, check=True, capture_output=True, text=True, timeout=60)
        rows = list(csv.DictReader(destination.open()))
        assert len(rows) == len(qs)
        maximum = 0.
        for q, row in zip(qs, rows):
            exact = nu+.09*k*k/epsilon if model == 'standard' else min(
                transverse(q, math.pi*i/360, side, k, epsilon, nu)
                for i in range(361) for side in (-1, 1))
            coarse = float(row['min_transverse_tangent_h1e5'])
            fine = float(row['min_transverse_tangent_h1e6'])
            scale = nu+float(row['eddy_viscosity'])
            error = abs(fine-exact)/scale
            assert math.isfinite(error) and error < 2e-5, (model, q, exact, fine)
            if model == 'cls':
                assert abs(fine-exact) <= .3*abs(coarse-exact)+5e-12, (q, exact, coarse, fine)
            maximum = max(maximum, error)
            row.update(analytic_minimum=exact, error_scaled_by_nu_plus_nut=error)
        results[model] = dict(command=command, maximum_scaled_error=maximum, rows=rows)
    # The stated retained-nu, fixed-k/epsilon probe brackets a sign change.
    cls = results['cls']['rows']
    assert float(cls[5]['analytic_minimum']) > 0 > float(cls[6]['analytic_minimum'])
    result = dict(accepted=True, k=k, epsilon=epsilon, resolved_nu=nu, results=results,
        interpretation='local transverse principal-part diagnostic at fixed positive k/epsilon; not a global flow-instability attribution')
    (output/'results.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--audit', type=Path, required=True)
    p.add_argument('--output', type=Path)
    args = p.parse_args()
    if args.output:
        result = check(args.audit, args.output.resolve())
    else:
        with tempfile.TemporaryDirectory() as temporary:
            result = check(args.audit, Path(temporary)/'audit')
    print('Transverse stress tangent oracle PASS:',
          {m: v['maximum_scaled_error'] for m, v in result['results'].items()})


if __name__ == '__main__':
    main()
