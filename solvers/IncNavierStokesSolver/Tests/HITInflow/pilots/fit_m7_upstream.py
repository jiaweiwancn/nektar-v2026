#!/usr/bin/env python3
"""Fit the uncut von Karman model using CBC station-42 RMS and Lf only.

Run numerical quadrature on the authoritative Linux host. SciPy is required.
The prescribed radial energy taper is applied after fitting, without rescaling.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path


def beta(a, b):
    return math.exp(math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, default=Path(__file__).resolve().parents[1] / 'cbc1971-source.json')
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    from scipy.integrate import quad
    raw = args.source.read_bytes()
    source = json.loads(raw)
    table = source['table4']
    index = {name: n for n, name in enumerate(table['columns'])}
    rows = [r for r in table['rows'] if r[index['theta']] == '42' and r[index['M_cm']] == '5.08']
    if len(rows) != 1:
        raise ValueError('require exactly one principal upstream Table-4 row')
    row = rows[0]
    urms = float(row[index['u_rms_cm_s']]) / 100
    lf = float(row[index['Lf_cm']]) / 100
    kinetic = 1.5 * urms**2
    i4, i3 = beta(2.5, 1/3)/2, beta(2, 5/6)/2
    length = lf / ((3*math.pi/4)*i3/i4)
    amplitude = kinetic / (length**(2/3)*i4)
    def energy(k):
        x = length*k
        return amplitude*length**(5/3)*x**4/(1+x*x)**(17/6)
    def taper(k):
        return 1 if k <= 80 else math.cos(math.pi*(k-80)/40)**2 if k < 100 else 0
    moments = {}
    for power in [-1, 0, 2]:
        parts = [quad(lambda k: energy(k)*taper(k)*k**power, a, b, epsabs=1e-13, epsrel=2e-10)
                 for a, b in [(0,80),(80,100)]]
        moments[str(power)] = {'value':math.fsum(v for v,e in parts),'absolute_quadrature_error':math.fsum(e for v,e in parts)}
    retained = moments['0']['value']; viscosity=1.5e-5
    result = {'schema':'NEKTAR_HIT_M7_UPSTREAM_FIT_V1','source':str(args.source.resolve()),
              'source_sha256':hashlib.sha256(raw).hexdigest(),
              'used_data':{'grid_M_m':0.0508,'theta':42,'u_rms_m_s':urms,'Lf_m':lf},
              'downstream_data_used':False,'upstream_dissipation_used_in_fit':False,
              'model':'E(k)=A*L^(5/3)*(L*k)^4/(1+(L*k)^2)^(17/6)',
              'I4':i4,'I3':i3,'length_m':length,'amplitude_m4_3_s2':amplitude,
              'uncut_K_m2_s2':kinetic,'uncut_Lf_m':lf,'uncut_dissipation':'divergent inertial tail',
              'cutoff_rad_m':100,'energy_taper_start_rad_m':80,'normalization':'Preserve',
              'represented_continuum':{'moments':moments,'K_m2_s2':retained,
                  'epsilon_m2_s3':2*viscosity*moments['2']['value'],
                  'Lf_m':3*math.pi*moments['-1']['value']/(4*retained),
                  'retained_fraction_of_upstream_K':retained/kinetic},
              'limitations':['Continuum input model, not finite-lattice or CFD measurement',
                  'Lf fitting uses isotropy and the unbounded uncut model; filtered scales differ',
                  'No restoration of energy lost to cutoff/taper; no downstream fitting']}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as out: json.dump(result,out,indent=2,allow_nan=False); out.write('\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
