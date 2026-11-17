#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Linux independent old-native parity and full-3D finite-difference checks."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from proposal_principal import full_wj


def independent(k,e,g):
    R=full_wj(k,e,g)
    S=(k/e)*((g+g.T)/2-np.trace(g)*np.eye(3)/3)
    O=(k/e)*(g-g.T)/2
    s=float(np.trace(S@S));w=float(np.trace(O@O))
    roots=np.roots([1.,-1.8,-2*w-2.7*s,3.6*w])
    n=max(v.real for v in roots if abs(v.imag)<1e-10*max(1.,abs(v.real))and v.real>=1.8*(1-1e-12))
    mu=3*n/(5*(n*n-2*w))*k*k/e
    return dict(R=R,eddyViscosity=mu,production=-float(np.sum(R*g)),root=n)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--native',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():p.error('immutable output exists')
    native=json.loads(a.native.read_text());assert native['schema']=='WJH27FastConstitutiveProbeV1'
    assert native['invalid_inputs_rejected']==5
    rows=[];worst_value=0.;worst_direction=0.;worst_fd=0.
    for row in native['states']:
        k,e=row['k'],row['epsilon'];g=np.asarray(row['G']);dg=np.asarray(row['dG']);old=row['old']
        production_scale=max(e,abs(old['production']),k*float(np.max(abs(g))),1e-300)
        eps_scale=max(e*e/k,abs(old['epsilonProduction']),abs(old['epsilonDestruction']),1e-300)
        scalar_scale=max(old['scalarViscosity'],1e-300);mu_scale=max(old['eddyViscosity'],1e-300)
        root=old['root'];s=old['s'];w=old['w']
        root_residual_scale=abs(root**3)+1.8*root**2+2*abs(w)*(root+1.8)+2.7*abs(s)*root
        def scale(name):
            if name in ('R','deviatoricStress','quadraticStress','cubicStress'):return k
            if name=='stressRemainder':return k+mu_scale*float(np.max(abs(g)))
            if name in ('eddyViscosity','momentumDiffusion'):return mu_scale
            if name in ('scalarViscosity','kDiffusion','epsilonDiffusion'):return scalar_scale
            if name in ('production','productionDev','kSource'):return production_scale
            if name in ('epsilonProduction','epsilonDestruction','epsilonSource'):return eps_scale
            if name=='rootResidual':return root_residual_scale
            return max(1.,float(np.max(abs(np.asarray(old[name])))))
        errors={}
        for label,left,right in [('value','fast','old'),('direction','fast_direction','old_direction'),
                                 ('paired_value','paired_value','fast'),('paired_direction','paired_direction','fast_direction')]:
            errors[label]={}
            for name in old:
                difference=float(np.max(abs(np.asarray(row[left][name])-np.asarray(row[right][name]))))/scale(name)
                assert math.isfinite(difference) and difference<=2e-11,(row['name'],label,name,difference)
                errors[label][name]=difference
        worst_value=max(worst_value,max(errors['value'].values()))
        worst_direction=max(worst_direction,max(errors['direction'].values()))
        for status in ('value_status','direction_status','paired_status'):
            state=row[status]
            if state['fallback']=='none':
                assert state['normalized_root_residual']<=64*np.finfo(float).eps
                assert state['relative_root_condition']>128*np.finfo(float).eps
        assert row['fast']['root']>=1.8*(1-1e-14)
        assert abs(np.trace(np.asarray(row['fast']['R']))-2*k)/k<=2e-11
        assert abs(row['fast']['production']+np.sum(np.asarray(row['fast']['R'])*g))/production_scale<=2e-12
        finite_difference=[];independent_value={}
        if not row['name'].startswith('fallback_'):
            exact=independent(k,e,g)
            for name,value in exact.items():
                error=float(np.max(abs(np.asarray(value)-np.asarray(row['fast'][name]))))/scale(name)
                independent_value[name]=error;assert error<2e-11,(row['name'],name,error)
            for step in (3e-5,1e-5):
                plus=independent(k+step*row['dk'],e+step*row['de'],g+step*dg)
                minus=independent(k-step*row['dk'],e-step*row['de'],g-step*dg)
                result={}
                for name in exact:
                    derivative=(plus[name]-minus[name])/(2*step)
                    error=float(np.max(abs(derivative-np.asarray(row['fast_direction'][name]))))/scale(name)
                    result[name]=error;assert error<=2e-8,(row['name'],step,name,error)
                    worst_fd=max(worst_fd,error)
                finite_difference.append(dict(step=step,scaled_errors=result))
        else:
            assert row['value_status']['fallback']!='none' and row['direction_status']['fallback']!='none'
        rows.append(dict(name=row['name'],native_parity_max={key:max(value.values())for key,value in errors.items()},
                         independent_value=independent_value,joint_FD=finite_difference,
                         value_status=row['value_status'],direction_status=row['direction_status']))
    assert native['value_fallbacks']>=3 and native['direction_fallbacks']>=3
    benchmark=native['benchmark']
    result=dict(schema='WJH27FastConstitutiveCheckV1',status='PASS',states=rows,
        maximum_value_parity_error=worst_value,maximum_direction_parity_error=worst_direction,
        maximum_joint_FD_error=worst_fd,value_fallbacks=native['value_fallbacks'],direction_fallbacks=native['direction_fallbacks'],
        kernel_only_speedup=dict(value=benchmark['old_value']['seconds']/benchmark['fast_value']['seconds'],
                                direction=benchmark['old_direction']['seconds']/benchmark['fast_direction']['seconds']),
        scope='Kernel-only equivalence/speed evidence; coupled recipe requires its own qualification',
        native_sha256=hashlib.sha256(a.native.read_bytes()).hexdigest())
    with a.output.open('x')as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({key:value for key,value in result.items()if key!='states'}))


if __name__=='__main__':main()
