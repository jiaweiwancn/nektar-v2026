#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent split-quadrature tests of the guarded CT6 Accepted kernel."""
import argparse, json, math, os
from pathlib import Path
from rans_qualification import THREADS, write_json
os.environ.update(THREADS)
import numpy as np
from scipy.optimize import brentq
from ct_accepted_feasibility import Column, parameters, require
from ct_accepted_open import run
from rans_request_manifest import RequestManifest, runtime_identity, sha


def native(runtime, output, label, inputs):
    root=output/label;root.mkdir()
    write_json(root/'input.json',inputs)
    execution=run([runtime/'RANSCurvedAcceptedKernelProbe',root/'input.json',root/'native.json'],root,'native.log')
    return json.loads((root/'native.json').read_text()),execution


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runtime',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    manifest=RequestManifest(a.output,dict(schema='CT6CurvedKernelQualificationV1',
        runtime=runtime_identity(a.runtime/'RANSCurvedAcceptedKernelProbe'),
        source_sha256=sha(Path(__file__)),models=['A93','A95','E94E95'],gamma=[-.25,0.,.25],
        inputs=dict(H=.5,hr=.05,nu=1.5e-5,K='q*.16'),orders=[64,96],
        tolerance=2e-10,curvature='B_g eigenvalue gamma; J=1+gamma*r; oracle curvature=-gamma'))
    rows=[];negative=[]
    for model in ('A93','A95','E94E95'):
        K=parameters(model)['q']*.16
        for gamma in (-.25,0.,.25):
            cols=[Column(model,K,.5,.05,1.5e-5,curvature=-gamma,order=o) for o in (64,96)]
            c=cols[-1];require(c.cutoff is not None,'predeclared_rough_kernel_case_inadmissible')
            locations=sorted(set([0.,c.cutoff*.5,c.cutoff,c.cutoff*1.5,.05,.125,.25,.5]))
            inputs=dict(model=model,K=K,H=.5,hr=.05,nu=1.5e-5,gamma=gamma,order=64,r=locations)
            got,execution=native(a.runtime,a.output,model+'-'+str(gamma),inputs)
            require(execution['exit_code']==0,'unexpected_kernel_rejection')
            expected=[[r,float(c.J(r)),float(c.production(r)),c.A/max(r,c.cutoff),c.k(r),
                       float(c.Q(r)/(c.J(r)*c.D(r))),float(c.Q(r)),float(c.D(r))] for r in locations]
            scales=np.array([.5,1.,c.A/.5,c.A/c.cutoff,K,K/.5,c.A,.4*.5])
            defects=np.max(np.abs(np.array(got['samples'])-np.array(expected))/scales,axis=0)
            delta=max(abs(got['cutoff']-c.cutoff)/.5,abs(got['budget']-c.PJ)/c.A,float(max(defects)))
            uncertainty=max(abs(cols[0].PJ-c.PJ)/c.A,abs(cols[0].cutoff-c.cutoff)/.5,
                max(abs(cols[0].k(r)-c.k(r))/K for r in locations))
            result=c.result()
            rows.append(dict(model=model,gamma=gamma,status='pass' if max(delta,uncertainty)<=2e-10 else 'fail',
                scaled_defect=delta,independent_quadrature_uncertainty=uncertainty,
                sample_scaled_defects=defects.tolist(),reference=result,execution=execution,
                covariance_scope='curved column positivity is not a curved-flow covariance certificate'))
            write_json(a.output/'progress.json',dict(cases=rows,negative=negative))
        # Independently locate the budget boundary, with a fixed bracket, not
        # a search for a passing physical input. Its two sides are controls.
        def margin(nu):
            col=Column(model,K,.5,.05,nu,curvature=-.25,order=96)
            return col.PJ-col.minimum_budget
        boundary=brentq(margin,1e-9,10.,xtol=1e-15)
        for label,nu,accepted in [('boundary-admissible',boundary*(1-1e-10),True),
                                  ('boundary-inadmissible',boundary*(1+1e-8),False)]:
            inputs=dict(model=model,K=K,H=.5,hr=.05,nu=nu,gamma=.25,order=96,r=[0.,.5])
            got,execution=native(a.runtime,a.output,model+'-'+label,inputs)
            reason=got.get('reason','')
            valid=(execution['exit_code']==0 and abs(got['cutoff']/.5-1)<2e-8) if accepted else (
                execution['exit_code']!=0 and 'budget' in reason.lower())
            negative.append(dict(model=model,kind=label,status='pass' if valid else 'fail',
                exact_independent_boundary_nu=boundary,margin=margin(nu),native=got,execution=execution))
        controls=[('unapproved-curvature',dict(gamma=.25,verification=False),'curved'),
                  ('endpoint-mismatch',dict(sample=.45),'endpoint'),
                  ('invalid-metric',dict(gamma=-3.),'metric')]
        for label,changes,reasonword in controls:
            inputs=dict(model=model,K=K,H=.5,hr=.05,nu=1.5e-5,gamma=.25,order=64,r=[0.,.5]);inputs.update(changes)
            got,execution=native(a.runtime,a.output,model+'-'+label,inputs)
            ok=execution['exit_code']!=0 and reasonword in got.get('reason','').lower()
            negative.append(dict(model=model,kind=label,status='pass' if ok else 'fail',native=got,execution=execution))
    out=dict(schema='CT6CurvedKernelQualificationV1',status='pass' if all(r['status']=='pass' for r in rows+negative) else 'fail',
        cases=rows,negative_controls=negative,scope='production profile; independent physical-breakpoint quadrature and moments; no native curved evolution claim')
    write_json(a.output/'results.json',out);manifest.seal(a.output,['results.json'])
    return 0 if out['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
