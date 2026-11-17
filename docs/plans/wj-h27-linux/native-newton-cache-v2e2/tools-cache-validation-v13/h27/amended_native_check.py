#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent amended point values/principal checks, not coupled qualification."""
import argparse,hashlib,json,math,sys
from pathlib import Path
import numpy as np

def main():
    p=argparse.ArgumentParser();p.add_argument('--native',type=Path,required=True)
    p.add_argument('--source-dir',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    sys.path.insert(0,str(a.source_dir));from proposal_principal import chi,proposal_tensor,full_wj
    d=json.loads(a.native.read_text());assert d['amended']
    rows=[]
    for v in d['states']:
        k,e,g=v['k'],v['epsilon'],v['shear'];nu,aux=v['nu'],v['nuA'];blend=chi(v['r']/v['H'])
        gradient=np.zeros((3,3));gradient[0,2]=g;r=proposal_tensor(k,e,gradient,aux,blend)
        mu=nu+(1-blend)*aux+blend*(-full_wj(k,e,gradient)[0,2]/g)
        errors=dict(tensor_over_k=float(np.max(abs(np.asarray(v['R'])-r))/k),
            mu_relative=abs(v['mu']-mu)/mu,
            production_over_E=abs(v['production']+float(np.sum(r*gradient)))/e)
        theta=v['theta'];n=np.array([math.cos(theta),0.,math.sin(theta)]);l=np.array([-math.sin(theta),0.,math.cos(theta)]);h=np.outer(l,n)
        steps=[]
        for relative in(2e-4,1e-4,5e-5):
            delta=relative*e/k
            dr=(proposal_tensor(k,e,gradient+delta*h,aux,blend)-proposal_tensor(k,e,gradient-delta*h,aux,blend))/(2*delta)
            value=float(np.sum(h*(nu*(h+h.T)-dr)));steps.append(dict(relative=relative,value=value,error_over_mu=(value-v['planar_min'])/mu))
        assert max(errors.values())<=1e-10,(v['label'],errors)
        assert abs(steps[-1]['error_over_mu'])<=2e-8,(v['label'],steps)
        assert v['planar_min_over_mu']>0,v['label']
        spectrum=np.linalg.eigvalsh(r)/k;assert spectrum[0]>0
        exact=None
        if blend==0:
            exact=(nu+aux*math.cos(2*theta)**2)/mu
            assert abs(exact-v['planar_min_over_mu'])<=1e-10
        rows.append(dict(label=v['label'],errors=errors,independent_FD=steps,
            native_principal_over_mu=v['planar_min_over_mu'],chi_zero_exact=exact,
            minimum_covariance_over_k=float(spectrum[0])))
    result=dict(status='PASS_NATIVE_AMENDED_POINT_VALUES_AND_PLANAR_DIRECTIONS',
        full_coupled_gate_pass=False,source=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        native_input=hashlib.sha256(a.native.read_bytes()).hexdigest(),rows=rows,
        scope='23 prescribed point states; no evolving scalar, collar continuum or PDE qualification')
    with a.out.open('x')as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(result['status'],len(rows))

if __name__=='__main__':main()
