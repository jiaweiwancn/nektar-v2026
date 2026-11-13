#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent original-law checks of static variable-H compatible collars."""
import argparse,json,math,os
from pathlib import Path
from rans_qualification import THREADS,write_json
os.environ.update(THREADS)
import numpy as np
from ct_accepted_feasibility import Column,parameters
from ct_accepted6_operator import geometry
from ct_accepted_open import require
from validate_a93_a95_e94e95 import oracle
from rans_request_manifest import RequestManifest,sha


def check(packet,model):
    maximum=0.;mincov=1e99;minrad=1.;flux_signal=0.;rows=0
    for row in packet['rows']:
        q=row['chart'][0];H=row['H'];j=row['compatible_collar'];K=j['K']
        a,J,gamma,t,n=geometry(q,H);Hp=.03*math.pi/2*math.cos(math.pi*q/2);qh=Hp/(a*J)
        N=(n-qh*t)/math.sqrt(1+qh*qh);area=a*J*(n-qh*t)
        c=Column(model,K,H,.05,1.5e-5,curvature=-gamma,order=96)
        require(c.cutoff is not None,'static independent metric cutoff inadmissible')
        Go,Gi=np.asarray(j['outer_G']),np.asarray(j['inner_G']);law=oracle(model,K,c.A/H,Go)
        nt=float(c.nt(H));DI=1.5e-5+nt;Ri=2*K/3*np.eye(3)-nt*(Gi+Gi.T)
        Ti=DI*(Gi+Gi.T)-j['inner_Pi']*np.eye(3)
        To=-law['R']+(2*K/3-j['outer_Pi'])*np.eye(3)
        U,ki,ko=map(np.asarray,(j['U'],j['inner_k_gradient'],j['outer_k_gradient']))
        fi=(float(U@area)*K-DI*float(ki@area))/a
        fo=(float(U@area)*K-law['nut']*float(ko@area))/a
        defects=[np.max(abs(N-row['normal'])),np.max(abs(area-row['cofactor'])),
            np.max(abs(Ri-j['inner_R'])),np.max(abs(law['R']-j['outer_R'])),np.max(abs(law['Q']-j['outer_Q'])),
            np.max(abs(Ti@N-j['inner_traction'])),np.max(abs(To@N-j['outer_traction'])),
            np.max(abs((Ti-To)@N)),abs(fi-fo),abs(fi-j['inner_conservative_k_flux']),abs(fo-j['outer_conservative_k_flux']),
            abs(j['inner_signed_production']+float(np.sum(Ri*Gi))),abs(j['outer_signed_production']-law['production']),
            abs(j['outer_epsilon_trace']-c.A/H),abs(j['dependent_inner_epsilon_trace']-c.A/H),
            abs(j['inner_mechanical_pressure']-(j['inner_Pi']-2*K/3)),abs(j['outer_mechanical_pressure']-(j['outer_Pi']-2*K/3))]
        maximum=max(maximum,*map(float,defects));flux_signal=max(flux_signal,abs(fi));rows+=1
        mincov=min(mincov,float(np.linalg.eigvalsh(Ri)[0]),float(np.linalg.eigvalsh(law['R'])[0]))
        p=parameters(model);S=(Go+Go.T)/2-np.trace(Go)/3*np.eye(3)
        if model=='A95':minrad=min(minrad,1-9*law['cmu']**2*(K/(c.A/H))**2*float(np.sum(S*S)))
    return dict(status='pass' if maximum<=2e-10 and mincov>=0 and minrad>=0 and flux_signal>1e-8 else 'fail',
        maximum_primitive_defect=maximum,minimum_covariance=mincov,minimum_A95_radicand=minrad if model=='A95' else None,
        nonzero_conservative_flux=flux_signal,rows=rows,
        scope='independent original-law and metric-column calculations on native common traces and constructed Hermite collar jets; full variable-height evolution remains guarded')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--native',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    path=a.native/'results.json';data=json.loads(path.read_text());m=RequestManifest(a.output,dict(schema='CT6IndependentStaticCollarsV1',native_receipt_sha256=sha(path),source_sha256=sha(Path(__file__))))
    checks={}
    for model in ('A95','A93','E94E95'):
        for degree in (3,4):
            label=model+'-static-p'+str(degree);r=data['checks'][label]
            require(r['execution']['exit_code']==0 and r['packet']['status']=='pass','native static fixture failed')
            checks[label]=check(r['packet'],model)
    out=dict(schema='CT6IndependentStaticCollarsV1',status='pass' if all(r['status']=='pass' for r in checks.values()) else 'fail',checks=checks)
    write_json(a.output/'results.json',out);m.seal(a.output,['results.json']);return 0 if out['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
