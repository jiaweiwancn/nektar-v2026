#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Linux independent finite-difference check of the native planar continuum minimum."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from proposal_principal import chi,proposal_tensor


def wj_mu(k,e,g):
    s=(k/e)*((g+g.T)/2-np.trace(g)*np.eye(3)/3);o=(k/e)*(g-g.T)/2
    si=float(np.trace(s@s));wi=float(np.trace(o@o))
    roots=np.roots([1.,-1.8,-2*wi-2.7*si,3.6*wi])
    n=max(v.real for v in roots if abs(v.imag)<1e-10 and v.real>=1.8*(1-1e-12))
    return 3*n/(5*(n*n-2*wi))*k*k/e


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--native',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():p.error('immutable output exists')
    native=json.loads(a.native.read_text());assert native['schema']=='WJH27PlanarAdmissibilityProbeV1'
    basis=np.zeros((3,3,3));basis[0,0,2]=-.5;basis[0,2,0]=.5
    basis[1,0,2]=basis[1,2,0]=.5;basis[2,0,0]=-.5;basis[2,2,2]=.5
    rows=[]
    for row in native['states']:
        k,e,nu,aux,r,H=(row[s] for s in ('k','epsilon','nu','auxiliary','r','H'))
        g=np.asarray(row['gradient']);blend=chi(r/H);mu=nu+(1-blend)*aux+blend*wj_mu(k,e,g)
        tensor=proposal_tensor(k,e,g,aux,blend)
        errors=dict(covariance_tensor=float(np.max(abs(np.asarray(row['R'])-tensor))/k),
            covariance_eigenvalue=abs(row['covariance_over_k']-np.linalg.eigvalsh(tensor)[0]/k),
            mu=abs(row['mu']-mu)/mu,
            minimum_native_action=abs(row['native_at_minimum']-row['principal_raw_over_mu']))
        assert max(errors.values())<2e-11,(row['name'],errors)
        assert row['principal_over_mu']<=min(row['grid'])+2e-12
        assert abs(row['principal_raw_over_mu']-row['principal_reserve_over_mu']-row['principal_over_mu'])<2e-14
        assert len(row['grid'])>=96 and 0<=row['theta']<math.pi
        steps=[]
        for relative in (2e-5,1e-5):
            delta=relative*e/k; actions=[]
            for direction in basis:
                derivative=(proposal_tensor(k,e,g+delta*direction,aux,blend)
                           -proposal_tensor(k,e,g-delta*direction,aux,blend))/(2*delta)
                actions.append(nu*(direction+direction.T)-derivative)
            q=np.einsum('iab,jab->ij',basis,np.asarray(actions))/mu
            def value(theta):
                v=np.array([1.,math.cos(2*theta),math.sin(2*theta)])
                return float(v@q@v)
            angles=np.arange(512)*math.pi/512;grid=np.asarray([value(t)for t in angles])
            candidate=[]
            for i in range(len(angles)):
                if grid[i]<=grid[i-1] and grid[i]<=grid[(i+1)%len(grid)]:
                    result=minimize_scalar(value,bounds=(angles[i]-math.pi/512,angles[i]+math.pi/512),
                        method='bounded',options={'xatol':1e-13})
                    candidate.append(float(result.fun))
            independent_minimum=min(candidate+[float(grid.min())])
            native_angles=np.arange(len(row['grid']))*math.pi/len(row['grid'])
            grid_error=max(abs(value(t)-v)for t,v in zip(native_angles,row['grid']))
            minimum_error=abs(independent_minimum-row['principal_raw_over_mu'])
            assert max(grid_error,minimum_error)<5e-8,(row['name'],relative,grid_error,minimum_error)
            steps.append(dict(relative_step=relative,grid_error_over_mu=grid_error,
                              independent_minimum=independent_minimum,minimum_error_over_mu=minimum_error))
        if blend==0:
            assert abs(row['principal_raw_over_mu']-nu/mu)<2e-14
            assert row['principal_reserve_over_mu']==0
        # Joint k/E/gradient/aux/r direction: independent full tensor values,
        # including the moving blend, never the native tangent implementation.
        direction=row['direction'];dg=np.asarray(direction['dg']);step=1e-5
        def shifted(sign):
            kk=k+sign*step*direction['dk'];ee=e+sign*step*direction['de']
            gg=g+sign*step*dg;aa=aux+sign*step*direction['da'];rr=r+sign*step*direction['dr']
            ch=chi(rr/H);stress=proposal_tensor(kk,ee,gg,aa,ch)
            vv=nu*(gg+gg.T)-stress+2*kk*np.eye(3)/3
            mm=nu+(1-ch)*aa+ch*wj_mu(kk,ee,gg)
            return stress,vv,mm,-float(np.sum(stress*gg))
        plus,minus=shifted(1),shifted(-1)
        direction_errors={}
        for index,(name,scale) in enumerate((('R',k),('V',k),('mu',mu),('production',e))):
            expected=(plus[index]-minus[index])/(2*step)
            error=float(np.max(abs(expected-np.asarray(direction[name]))))/scale
            direction_errors[name]=error
        assert max(direction_errors.values())<5e-8,(row['name'],direction_errors)
        rows.append(dict(name=row['name'],value_errors=errors,principal_FD=steps,
                         joint_direction_errors=direction_errors,reserve_over_mu=row['principal_reserve_over_mu'],
                         stationary_residual=row['stationary_residual']))
    assert native['negative_covariance_count']>0 and native['negative_principal_count']>0
    assert native['unsupported_spanwise_rejected']
    result=dict(schema='WJH27PlanarAdmissibilityCheckV1',status='PASS',states=rows,
        native_sha256=hashlib.sha256(a.native.read_bytes()).hexdigest(),
        scope='Finite independent checks of continuum planar minimizer, negative witnesses, and joint tangent; no global 3D claim')
    with a.output.open('x')as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print('PASS',len(rows),'planar states')


if __name__=='__main__':main()
