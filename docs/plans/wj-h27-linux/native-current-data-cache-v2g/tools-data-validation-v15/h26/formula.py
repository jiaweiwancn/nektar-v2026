#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independently differentiate the physical PDE and compare compiled samples."""
import argparse,json
from pathlib import Path
import sympy as s
import numpy as np

def verify(source,out):
    x,z,t=s.symbols('x z t');nu=s.Float('6.423173803526447e-5',24);mu=nu+(1+s.sin(s.pi*t)/10)*(1+x+z)/100
    pressure=s.exp(-2*t)*(x*x+z*z)/100
    T=s.exp(-t)*s.Matrix([[x*z,(x+z)/2],[(x+z)/2,-x*z]])/1000
    compiled=json.loads(source.read_text());functions={};maxima={}
    for shape in ('polynomial','trigonometric'):
        f=lambda a:a*a*(1-a)if shape=='polynomial'else s.sin(s.pi*a)**2
        psi=s.exp(-t)*f(x)*f(z);u=s.Matrix([s.diff(psi,z),-s.diff(psi,x)]);G=u.jacobian([x,z]);stress=mu*(G+G.T)+T
        force=s.diff(u,t)+G*u+s.Matrix([s.diff(pressure,x),s.diff(pressure,z)])-s.Matrix([sum(s.diff(stress[i,j],(x,z)[j])for j in range(2))for i in range(2)])
        assert s.simplify(G.trace())==0
        functions[shape]=s.lambdify((x,z,t),force,'numpy');maxima[shape]=0.
    for record in compiled['samples']:
        ref=np.array(functions[record['shape']](record['x'],record['z'],record['t']),float).ravel();actual=np.array(record['force']);error=float(np.linalg.norm(ref-actual)/(np.linalg.norm(ref)))
        maxima[record['shape']]=max(maxima[record['shape']],error)
        assert error<=1e-11,(record,error,ref)
    out.write_text(json.dumps(dict(status='PASS',maximum_relative_errors=maxima,samples=len(compiled['samples']),derivation='independent SymPy continuous physical derivatives; compiled Force does not read generated expressions'),indent=2)+'\n');print(maxima)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();verify(a.source,a.out)
