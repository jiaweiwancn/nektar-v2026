#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent Linux-only checker for native H27 reference_probe output."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from numpy.polynomial.chebyshev import chebval, chebder
from scipy.integrate import quad

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from wj_h2_profiles import wall


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--table',type=Path,required=True)
    parser.add_argument('--probe',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():parser.error('immutable checker output exists')
    table=json.loads(args.table.read_text()); probe=json.loads(args.probe.read_text())
    assert table['status']=='PASS' and probe['schema']=='WJH27NativeReferenceProbeV1'
    assert probe['table_sha256']==hashlib.sha256(args.table.read_bytes()).hexdigest()
    p=table['parameters']; const=table['constants']; Kin=p['K']; H=p['H']; nu=p['nu']; h=p['h']
    qe=const['q_e']; kappa=const['kappa']; cl=const['c_length']; A=Kin**1.5/cl
    us=table['inlet_binding']['u_tau']; top=table['inlet_binding']['top_model_radius']
    tau=Kin/qe; old=wall(Kin,H,h,nu=nu); coeff=table['cutoff_table']['coefficients']
    lo,hi=table['cutoff_table']['domain']; z0=.005676270321216798
    errors={}
    def check(name,actual,expected,scale=1.,tolerance=2e-11):
        error=abs(float(actual)-float(expected))/scale
        errors[name]=max(errors.get(name,0.),error)
        assert math.isfinite(error) and error<=tolerance,(name,actual,expected,error,tolerance)
    def aux(K,r):
        x=h*np.sqrt(K)/(70*nu)
        # Select exponent by the real base state, then complex-step its exact
        # smooth branch independently of the native analytic derivative.
        m1=.5-.4*x**.7; m2=1-.79*x**(-.28)
        exponent=m1 if m1.real>m2.real else m2
        b=-10.7*np.expm1(exponent*np.log(x)); zv=nu*b/np.sqrt(K)
        return kappa*np.sqrt(K/qe)*(r-zv)
    for row in probe['dynamic']:
        K,r=row['K'],row['r']; x=(2*K-lo-hi)/(hi-lo)
        logc=float(chebval(x,coeff)); c=H*math.exp(logc)
        dc=c*float(chebval(x,chebder(coeff)))*2/(hi-lo)
        a=K**1.5/cl; ap=1.5*math.sqrt(K)/cl
        for key,expected in [('A',a),('dA',ap),('ddA',.75/(math.sqrt(K)*cl)),('c',c),('dc',dc)]:
            check(key,row[key],expected,max(abs(expected),1e-4))
        actual_aux=row['aux']
        check('aux',actual_aux[0],aux(K,r),float(aux(K,r)))
        check('aux_dK',actual_aux[1],aux(K+1e-28j,r).imag/1e-28,max(abs(actual_aux[1]),1e-4))
        check('aux_dr',actual_aux[2],kappa*math.sqrt(K/qe),kappa*math.sqrt(K/qe))
        expected=(a/c,ap/c-a*dc/c**2,0.) if r<=c else (a/r,ap/r,-a/r**2)
        for j,name in enumerate(('Ew','Ew_dK','Ew_dr')):
            check(name,row['epsilon'][j],expected[j],max(abs(expected[j]),1e-4))
    segments=table['outlet_traction']['segments']
    for row in probe['inlet']:
        r=row['r']; D=nu+old['a']*max(r-old['zv'],0.)
        U=us**2*float(old['J'](r)) if r<=H else us/kappa*math.log(r/z0)
        k=old['kref'](r) if r<H else Kin
        check('historical_U',row['U'],U,us,tolerance=5e-12)
        check('historical_k',row['k'],k,Kin,tolerance=5e-12)
        check('historical_U_r',row['U_r'],us**2/D if r<=H else us/(kappa*r),max(abs(row['U_r']),1e-4))
        check('historical_k_r',row['k_r'],old['Q'](r)/D if r<H else 0.,Kin)
        if r>=H:
            check('historical_epsilon',row['epsilon'],A/r,A/r)
            check('historical_epsilon_r',row['epsilon_r'],-A/r**2,A/r**2)
        else:
            assert row['epsilon'] is None and row['epsilon_r'] is None
        expected=Kin/const['N_e']
        if r>.4*H:
            segment=segments[0] if r<=.75*H else segments[1] if r<=H else segments[2]
            argument=math.log(r/H) if r>H else r
            left,right=segment['domain']; expected=float(chebval((2*argument-left-right)/(right-left),segment['tx_coefficients']))
        check('outlet_tx',row['traction'][0],expected,tau,tolerance=2e-13)
        check('outlet_tz',row['traction'][1],tau,tau,tolerance=2e-13)
    check('historical_cutoff',probe['historical_inlet_cutoff'],old['rc'],H)
    flux=nu*A/top**2+p['cd']*Kin**2/(table['guard']['sigma_e']*top)
    check('top_epsilon_flux',probe['top_epsilon_outward_flux'],flux,flux)
    assert len(probe['rejections'])==8
    result=dict(schema='WJH27NativeReferenceCheckV1',status='PASS',maximum_scaled_errors=errors,
        dynamic_rows=len(probe['dynamic']),inlet_rows=len(probe['inlet']),rejection_count=len(probe['rejections']),
        table_sha256=probe['table_sha256'],scope='Native lookup/analytic derivatives/historical inlet only; no coupled qualification')
    with args.output.open('x') as output:json.dump(result,output,indent=2);output.write('\n')
    print(json.dumps(result))


if __name__=='__main__':main()
