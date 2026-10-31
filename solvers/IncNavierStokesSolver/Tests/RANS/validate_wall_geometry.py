#!/usr/bin/env python3
"""Independent Linux quadrature of normal-metric AWF sources (no solver reset)."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not sys.platform.startswith('linux'):
        raise RuntimeError('Use the authoritative Linux environment')
    from scipy.integrate import quad
    args.output.mkdir(parents=True, exist_ok=True)
    receipt = args.output/'geometry-kernel.json'
    run = subprocess.run([str(args.probe.resolve()), str(receipt)], capture_output=True, text=True, timeout=120)
    (args.output/'probe.log').write_text(run.stdout+run.stderr)
    if run.returncode:
        raise RuntimeError('Geometry probe failed; inspect '+str(args.output/'probe.log'))
    result = json.loads(receipt.read_text())
    maximum = 0
    for c in result['metric_cases']:
        H, nu, h, k, curvature, chi = .008, 1.5e-5, c['roughness'], c['k'], c['curvature'], c['chi']
        root = math.sqrt(k)
        r = h*root/(70*nu)
        intercept = 10.7 if h == 0 else 10.7*(1-r**max(.5-.4*r**.7, 1-.79*r**(-.28)))
        zv, slope, ze = nu*intercept/root, .09*2.55*root, 2*2.55*nu/root
        breaks = {0., H}
        for p in (h, zv, ze, c['cutoff']):
            if 0 < p < H:
                breaks.add(p)
        p = nu/slope
        while p < H:
            breaks.add(p)
            p *= 2
        def integrate(f):
            b = sorted(breaks)
            return sum(quad(f, a, b, epsabs=1e-12, epsrel=2e-12)[0] for a, b in zip(b, b[1:]))
        nut = lambda z: slope*max(0., z-zv)
        total = lambda z: nu+nut(z)
        j0 = integrate(lambda z: 1/total(z))
        j1 = integrate(lambda z: max(0., z-h)/total(z))
        velocity, forcing = (4., -.2), (.03, -.01)
        tau = [(u-f*j1)/j0 for u, f in zip(velocity, forcing)]
        volume = integrate(lambda z: 1-curvature*z)
        production = integrate(lambda z: (1-curvature*z)*nut(z)*sum(
            ((t+f*max(0., z-h))/total(z))**2 for t, f in zip(tau, forcing)))/volume
        dissipation = integrate(lambda z: (1-curvature*z)*((1-chi)*k**1.5/(2.55*max(z, ze))
            +chi*c['A']/max(z, c['cutoff'])))/volume
        shear = integrate(lambda z: tau[0]+forcing[0]*max(0., z-h))
        for key, expected in [('P', production), ('D', dissipation), ('shear_x', shear)]:
            error = abs(c[key]-expected)/max(abs(expected), 1e-30)
            maximum = max(maximum, error)
            if error > 5e-10:
                raise RuntimeError(f'Independent metric {key}: {error}: {c}')
    evidence = dict(status='PASS', cases=len(result['metric_cases']),
                    maximum_relative_quadrature_error=maximum,
                    kernel_checks=result['checks'], invalid_rejections=result['rejections'],
                    command=[str(args.probe.resolve()), str(receipt)],
                    binary_sha256=hashlib.sha256(args.probe.read_bytes()).hexdigest())
    (args.output/'evidence.json').write_text(json.dumps(evidence, indent=2)+'\n')
    print(json.dumps(evidence, indent=2))


if __name__ == '__main__':
    main()
